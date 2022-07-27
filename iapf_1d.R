library(mvtnorm)
library(MASS)
library(profvis)
library(FKF)

####settings####
set.seed(456)
alpha = A = 0.42
d <- 1
k <- 5
kappa <- 0.5
tau <- 0.5
m <- 0
N <- vector()
N[1] <- 1000
Time <- 100
cov = B = C = D = 1
X <- matrix(0, nrow = Time, ncol = 30000 )#particles t horizontal; N vertical, d
X_true <- vector()
obs <- vector()
w <- matrix(NA, nrow = Time, ncol = 30000 ) #weight
Z <- vector()  #approximation

# FKF

dt <- ct <- 0
Tt <- as.matrix(A)
P0 <- Zt <- Ht <- Gt <- as.matrix(1)
a0 <- 0

f <- function(x){
  return (rnorm(d) + as.vector(A*x))   #trans prob
}

g <- function(y, x){  
  return (det(diag(2*pi, nrow = d, ncol = d))^(-1/2)*exp((-1/2)*(y-x)^2)) #obs prob  C%*%x = x
  #dnorm(X[t,1:N[l]],obs[t])
}

psi <- matrix(NA, nrow = Time, ncol = 30000) #iterated psi to decide psi_t for each l

mu_aux <- function(psi_pa, l){
  return(rnorm(N[l], mean =  psi_pa[1,1]/(1+psi_pa[1,2]), #change mean and var!
                 sd = psi_pa[1,2]/(1+psi_pa[1,2])))
}

g_aux <- function(y, x, t, psi_pa){  
  if(t == 1){
    return(dnorm(x, y)*psi_tilda(x, psi_pa, 1)*(2*pi*(psi_pa[1, d+1]^2+1))^
                                              (-1/2)*exp((-1/2)*(-psi_pa[1, 1])^2/(psi_pa[1,  d+1]^2+1))/psi_t(x, psi_pa, 1))  #g_1
  }else{
    return(dnorm(x, y)*psi_tilda(x, psi_pa, t)/psi_t(x, psi_pa, t))  #g_2:T 
  }
}

f_aux <- function(x, psi_pa, t){    #改！
  return(rnorm(1, (psi_pa[t,2]^2*A*x +psi_pa[t,1])/(psi_pa[t,2]^2+1), 
                 sqrt(psi_pa[t,2]^2 / (psi_pa[t,2]^2+1))))  #f_2:T 
}#rnorm(N,(psi[t,2]^2 * a*X[t-1,anc] + sigma^2*psi[t,1])/(psi[t,2]^2+sigma^2), sd=sqrt(psi[t,2]^2 * sigma^2/ (psi[t,2]^2+sigma^2)))

ESS <- function(t,l,w){
  return(sum(w[t-1,1:N[l]])^2/sum(w[t-1,1:N[l]]^2))
}

psi_pa <- matrix(NA, nrow = Time, ncol = 2)

Num <- function(Z, l, k){
  return(sd(exp(Z[max(l-k,1):l]-max(Z)))/mean(exp(Z[max(l-k,1):l]-max(Z))))
}

Obs <- function(){
  X_true[1] <- rnorm(d) + m    #real writing!
  for(t in 2:Time){  #observations
    X_true[t] <- f(X_true[t-1])   #t(rnorm(d) + A%*%x)
  }
  return(rnorm(Time) + X_true)
}


#1. Initialize
#l=0  large outer iteration
#psi_t: 1-T+1  psi_T+1 =1
#psi_tilda: 0-T  index +1

#use psi_t to calulcate psi_tilda[t]

psi_tilda <- function(x, psi_pa, t){  #from 0 to T. 0,T = 1  ?????
  if(t == Time){
    psi_tilda <- 1
  }else{   #psi_pa_t = psi_t
    psi_tilda <- (2*pi*(psi_pa[t+1, d+1]^2+1))^   #dnorm(A*x, psi_pa[t+1, 1], sqrt(psi_pa[t+1,2]^2+sigma^2)))
      (-1/2)*exp((-1/2)*(A*x-psi_pa[t+1, 1])^2/(psi_pa[t+1, d+1]^2+1))  #f(xt, Ït+1 )   #var should be psi^2+1?
  }
  return(psi_tilda)
}

psi_t <- function(x, psi_pa, t){ #from 1 to T+1. 1, T+1 = 1  
  if(t == Time + 1){
    psi_t <- 1
  }else{
    psi_t <- (2*pi*psi_pa[t, d+1]^2)^(-1/2)*						  #psi_pa is sd not var!!!
      exp((-1/2)*(x-psi_pa[t, 1])^2/(psi_pa[t, d+1]^2))  #dnorm(x, psi_pa[t, 1], psi_pa[t, d+1])
  }
  
  return(psi_t)
}


####APF function####
APF <- function(psi_pa, l){
  #l >= 2
  Z[l] <- 0
  
  X[1,1:N[l]] <- mu_aux(psi_pa, l)  #particles
  
  for(i in 1:N[l]){
    w[1,i] <- g_aux(obs[1], X[1,i],1, psi_pa) #weights g(obs[1,], X[1,i,])*psi_tilda(X[1,i,], psi_pa, 2)  
  }
  re=0
  #t=2:T
  #2. conditional sample
  for(t in 2:Time){
    
    print(t)
    
    #a)
    
    if(ESS(t,l,w) <= kappa*N[l]){
      re = re+1
      w_ <- w[t-1,1:N[l]]/sum(w[t-1,1:N[l]])   #each t
      mix <- sample(1:N[l],N[l], replace = TRUE, prob = w_)
      
      for(i in 1:N[l]){
        X[t,i] <- f_aux(X[t-1, mix[i]], psi_pa, t)
        w[t,i] <- g_aux(obs[t], X[t,i], t, psi_pa)  
      }
    }else{
      
      #b)
      for(i in 1:N[l]){
        
        X[t,i] <- f_aux(X[t-1,i],psi_pa, t) 
        w[t,i] <- w[t-1,i]*g_aux(obs[t], X[t,i],t, psi_pa)  
      }
    }
    
  }
  
  for(t in 1:Time){
    sum_g = 0
    for(i in 1:N[l]){
      sum_g = sum_g + g_aux(obs[t],X[t,i],t,psi_pa)
    }
    Z[l] = Z[l] + log(sum_g/N[l])
  }
  print(paste0('re=',re))
  
  return(list(obs, X, w, Z))
}

####psi function####
Psi <- function(l, obs, X){
  
  
  for(t in Time:1){
    
    print(t)
    
    if(t == Time){
      psi[t,1:N[l]] <- dnorm(X[t,1:N[l]], obs[t])  
      
    }else{
      
      for(i in 1:N[l]){
        psi[t,i] <- dnorm(X[t,i], obs[t])*dnorm(A*X[t,i],psi_pa[t+1, 1],sqrt(psi_pa[t+1,2]^2+1))
        #*(2*pi*(psi_pa[t+1, d+1]+1))^(-1/2)*
          #exp((-1/2)*(A*X[t,i]-psi_pa[t+1, 1])^2/(psi_pa[t+1, d+1]+1))   #?
      
      }
    }
    
    
    #2. calculate psi_t
    #calculate min
    fn <- function(x, X, psi){
      lambda <-  2*sum(dnorm(X[t,1:N[l]],mean=x[1],sd=x[2]) * psi[t,1:N[l]]) / sum(psi[t,1:N[l]]^2)
      return(sum((psi[t,1:N[l]] - (1/lambda)*dnorm(X[t,1:N[l]],mean=x[1],sd=x[2]))^2))
      #exp(-1/2*((X[t,1:N[l]]-x[1])/x[2])^2)/(x[2]*sqrt(2*pi))
      
      #sum_arg = 0						
      #for(i in 1:N[l]){						
        #sum_arg = sum_arg + (x[2*d+1]*dnorm(X[t,i],mean = x[1],sd=x[2])-	psi[t,i])^2						
                               
                               					
      #}						
      #return(sum_arg)
    }
    
    #get the distribution of psi_t
    if(t == Time){
      psi_pa[t,] <- optim(par = c(mean(X[t,1:N[l]]),1),
                          fn = fn, X = X, psi = psi, method='L-BFGS-B',lower=c(-Inf,0,0),upper=c(Inf,Inf,Inf))$par
    }else{
      psi_pa[t,] <- optim(par = c(X[t,which.max(psi[t,1:N[l]])],1), 
                          fn = fn, X = X, psi = psi, method='L-BFGS-B',lower=c(-Inf,0,0),upper=c(Inf,Inf,Inf))$par
    }
    
    print(psi_pa[t,])
    print(obs[t])
  }
  
  return(psi_pa)
  
  }


####2. repeat####
obs <- Obs()
fkf.obj <- fkf(a0, P0, dt, ct, Tt, Zt, Ht, Gt, yt = t(obs))$logLik

index = 1
l = 1  #actually 0

#initialization l=1, run a psi^0 apf with N[1]
#psi_1...psi_T+1 = 1
#psi.tilda_0...psi.tilda_T = 1

#t = 1
#generate particles and weights

X[1,1:N[l]] <- rnorm(N[l]*d)  #particles
w[1,1:N[l]] <- dnorm(X[1,1:N[l]], obs[1])   #problem!!!!
#for(i in 1:N[l]){
  #w[1,i] <- g(obs[1], X[1,i])  #weights
#}
re = 0
#t=2:T
#2. conditional sample

for(t in 2:Time){
  print(t)
  
  #a)
  
  if(ESS(t,l,w) <= kappa*N[l]){
    
    re = re + 1
    
    w_ <- w[t-1,1:N[l]]/sum(w[t-1,1:N[l]])
    mix <- sample(1:N[l], N[l], replace = TRUE, prob = w_)
    X[t,1:N[l]] <- rnorm(Time) + A*X[t-1, mix]   #rnorm should be different!
    w[t,1:N[l]] <- dnorm(X[t,1:N[l]], obs[t]) 
    
  #$for(i in 1:N[l]){
      
      #w[t,i] <- g(obs[t], X[t,i])  
    #}
    
  }else{
    
    #b)
    
    for(i in 1:N[l]){
      X[t,i] <- f(X[t-1,i])   
      w[t,i] <- w[t-1,i]*dnorm(X[t,i], obs[t])  
    }
  }
  
}

Z[1] = 0
for(t in 1:Time){
  Z[l] = Z[l] + log(mean(apply(as.matrix(X[t,1:N[l]]),1,function(x) dnorm(x, obs[t]))))
}

while(index){
  
  print(l)
  #a)
  output <- list()
  
  if(l != 1){
    output <-  APF(psi_pa, l)
    obs <-output[[1]]
    X <- output[[2]]
    w <- output[[3]]
    Z <- output[[4]]
  }
  
  #b)
  
  if(l <= k | (Num(Z, l, k) >= tau)){
    #psi^{l+1}
    
    psi_pa <- Psi(l, obs, X) 
    
    
    if(l > k & N[max(l-k,1)] == N[l] & is.unsorted(Z[max(l-k,1):l])){
      N[l+1] <- 2*N[l]
      
    }else{
      N[l+1] <- N[l]
    }
    
    print(paste0('Z[l]=',Z[l]))
    
    print(paste0('Z=',fkf.obj))
    
    l <- l+1
  }else break
}

#3.
output <-  APF(psi_pa, l)
Z <- output[[4]]
Z_appro <- Z[l]



