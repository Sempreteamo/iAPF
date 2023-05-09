library(mvtnorm)
library(MASS)
library(profvis)
library(FKF)

####settings####
set.seed(672)
alpha <- 0.42
d <- 5
k <- 5
kappa <- 0.5
tau <- 0.5
m <- rep(0, d)
N <- vector()
N[1] <- 100
Time <- 100
cov = B = C = D = diag(1, nrow = d, ncol = d)
A <- matrix(nrow = d, ncol = d)
for (i in 1:d){
  for (j in 1:d){
    A[i,j] = alpha^(abs(i-j) + 1)
  }
}
X <- array(NA, dim = c(Time, N, d))
#X <- array(0, dim = c(Time, 30000, d)) #particles t horizontal; N vertical, d
X_true <- matrix(0, nrow = Time, ncol = d )
obs <- matrix(0, nrow = Time, ncol = d )
w <- matrix(NA, Time, N)
#w <- matrix(NA, nrow = Time, ncol = 30000 ) #weight
Z <- vector()  #approximation

# FKF

dt <- ct <- matrix(0,d,1)
Tt <- A
P0 <- Zt <- Ht <- Gt <- diag(1,d,d)
a0 <- rep(0,d)

f <- function(x){
  return (rnorm(d) + as.vector(A%*%x))   #trans prob
}

g <- function(y, x){  
  return (as.numeric(det(diag(2*pi, nrow = d, ncol = d))^(-1/2)*exp((-1/2)*t(y-x)%*%(y-x)))) #obs prob  C%*%x = x
  #det(diag(2*pi, nrow = d, ncol = d))^(-1/2) = 0.01010533
  
}

#psi <- matrix(NA, nrow = Time, ncol = 30000) #iterated psi to decide psi_t for each l

mu_aux <- function(psi_pa, l){   #???
  return(mvrnorm(N[l], mu =  as.vector(diag(((psi_pa[1, (d+1):(d+d)])^(-1)+1)^(-1), nrow=d,ncol=d)%*%
                                         (diag((psi_pa[1, (d+1):(d+d)])^(-1), nrow=d,ncol=d)%*%psi_pa[1,1:d])), 
                 Sigma = diag(((psi_pa[1, (d+1):(d+d)])^(-1)+1)^(-1), nrow=d,ncol=d)))
}

g_aux <- function(y, x, t, psi_pa){   
  if(t == 1){
    return(g(y, x)*psi_tilda(x, psi_pa, 1)*((2*pi)^(-d/2)*det(diag((psi_pa[1, (d+1):(d+d)]+1), nrow=d,ncol=d))^
                                              (-1/2)*exp((-1/2)*t(-psi_pa[1, 1:d])%*%
                                                           diag((psi_pa[1, (d+1):(d+d)]+1)^(-1), nrow=d,ncol=d)%*%
                                                           (-psi_pa[1, 1:d])))/psi_t(x, psi_pa, 1))  #g_1
  }else{
    return(g(y, x)*psi_tilda(x, psi_pa, t)/psi_t(x, psi_pa, t))  #g_2:T 
  }
}

f_aux <- function(x, psi_pa, t){  #?? 
  return(mvrnorm(1, mu = as.vector(diag(((psi_pa[t, (d+1):(d+d)])^(-1)+1)^(-1), nrow=d,ncol=d)%*%
                                     (A%*%x + diag(psi_pa[t, (d+1):(d+d)]^(-1), nrow=d,ncol=d)%*%psi_pa[t,1:d])), 
                 Sigma = diag(((psi_pa[t, (d+1):(d+d)])^(-1)+1)^(-1), nrow=d,ncol=d)))
}                   
#diag((psi_pa[t, (d+1):(d+d)]+1)^(-1), nrow = d, ncol = d)%*%
#diag(psi_pa[t, (d+1):(d+d)], nrow = d, ncol = d)))  #f_2:T 

#as.vector(diag((psi_pa[t, (d+1):(d+d)]+1)^(-1), nrow=d,ncol=d)%*%
#(diag(psi_pa[t, (d+1):(d+d)], nrow=d,ncol=d)%*%A%*%x + psi_pa[t,1:d]))


ESS <- function(t,l,w, is.log=FALSE){
  if(is.log) {
    mx <- max(w[t-1,1:N[l]])
    s <- sum(exp(w[t-1,1:N[l]]-mx))
    ess <- 1/sum((exp(w[t-1,1:N[l]]-mx)/s)^2)
  }else{
    s <- sum(w[t-1,1:N[l]])
    ess <- 1/sum((w[t-1,1:N[l]]/s)^2) 
  }
  return(ess)  
}
#sum(w[t-1,1:N[l]])^2/sum(w[t-1,1:N[l]]^2)

#psi_pa <- matrix(NA, nrow = Time, ncol = 2*d)

Num <- function(Z, l, k){
  return(sd(exp(Z[max(l-k,1):l]-max(Z[max(l-k,1):l])))/mean(exp(Z[max(l-k,1):l]-max(Z[max(l-k,1):l]))))
}

Obs <- function(){
  X_true[1,] <- rnorm(d) + m    #real writing!
  for(t in 2:Time){  #observations
    X_true[t,] <- f(X_true[t-1,])   #t(rnorm(d) + A%*%x)
  }
  return(matrix(rnorm(Time*d),Time,d) + X_true)
}


#1. Initialize
#l=0  large outer iteration
#psi_t: 1-T+1  psi_T+1 =1

#psi_tilda: 0-T  index +1

#use psi_t to calulcate psi_tilda[t]

psi_tilda <- function(x, psi_pa, t){  #from 0 to T. 0,T = 1
  if(t == Time){
    psi_tilda <- 1
  }else{   #psi_pa_t = psi_t
    psi_tilda <- (2*pi)^(-d/2)*det(diag(psi_pa[t+1, (d+1):(d+d)]+1, nrow=d, ncol=d))^(-1/2)*
      exp((-1/2)*t(as.vector(A%*%x) - psi_pa[t+1, 1:d])%*%diag((psi_pa[t+1, (d+1):(d+d)]+1)^(-1), nrow=d,ncol=d)%*%
            (as.vector(A%*%x)-psi_pa[t+1, 1:d]))
  }
  return(psi_tilda)
}

psi_t <- function(x, psi_pa, t){ #from 1 to T+1. 1, T+1 = 1  
  if(t == Time + 1){
    psi_t <- 1
  }else{
    psi_t <- (2*pi)^(-d/2)*det(diag(psi_pa[t, (d+1):(d+d)], nrow=d,ncol=d))^(-1/2)*						
      exp((-1/2)*t(x-psi_pa[t, 1:d])%*%diag((psi_pa[t, (d+1):(d+d)])^(-1), nrow=d,ncol=d)%*%						
            (x-psi_pa[t, 1:d]))
  }
  
  return(psi_t)
}


####APF function####
APF <- function(psi_pa, l, N){
  #l >= 2
  X <- array(NA, dim = c(Time, N[l], d))
  w <- matrix(NA, Time, N[l])
  Z[l] <- 0
  
  X[1,1:N[l],] <- mu_aux(psi_pa, l)  #particles
  
  for(i in 1:N[l]){
    w[1,i] <- log(g_aux(obs[1,], X[1,i,],1, psi_pa)) #weights g(obs[1,], X[1,i,])*psi_tilda(X[1,i,], psi_pa, 2)  
  }
  re=0
  
  #t=2:T
  #2. conditional sample
  for(t in 2:Time){
    
    #print(t)
    
    #a)
    
    if(ESS(t,l,w, is.log=TRUE) <= kappa*N[l]){
      mx <- max(w[t-1, 1:N[l]])
      Z[l] <- Z[l] + log(mean(exp(w[t-1, 1:N[l]]-mx))) + mx	
      re = re+1
      w_ <- exp(w[t-1,1:N[l]]-mx)/sum(exp(w[t-1,1:N[l]]-mx))   #each t
      mix <- sample(1:N[l],N[l], replace = TRUE, prob = w_)
      
      for(i in 1:N[l]){
        X[t,i,] <- f_aux(X[t-1, mix[i],], psi_pa, t)
        w[t,i] <- log(g_aux(obs[t,], X[t,i,], t, psi_pa))  
      }
    }else{
      
      #b)
      for(i in 1:N[l]){
        
        X[t,i,] <- f_aux(X[t-1,i,],psi_pa, t) 
        w[t,i] <- w[t-1,i] + log(g_aux(obs[t,], X[t,i,],t, psi_pa))  
      }
    }
    
  }
  mx <- max(w[t,1:N])
  Z[l] <- Z[l] + log(mean(exp(w[t,1:N]-mx))) + mx
  
  #for(t in 1:Time){
  #sum_g = 0
  #for(i in 1:N[l]){
  #sum_g = sum_g + g_aux(obs[t,],X[t,i,],t,psi_pa)
  #}
  #Z[l] = Z[l] + log(sum_g/N[l])
  #}
  #print(paste0('re=',re))
  
  return(list(obs, X, w, Z))
}

####psi function####
Psi <- function(l, obs, X){
  psi <- matrix(NA, nrow = Time, ncol = N[l])
  psi_pa <- matrix(NA, nrow = Time, ncol = 2*d)
  
  for(t in Time:1){
    
    #print(t)
    
    if(t == Time){
      psi[t,1:N[l]] <- dmvnorm(X[t,1:N[l],], obs[t,])
      
    }else{
      
      for(i in 1:N[l]){
        psi[t,i] <- g(obs[t,],X[t,i,])*dmvnorm(as.vector(A%*%X[t,i,]), psi_pa[t+1, 1:d], diag(psi_pa[t+1, (d+1):(d+d)]+1, nrow=d,ncol=d))
        
      }
      #det(2*pi*(diag(psi_pa[t+1, (d+1):(d+d)]+1, nrow=d,ncol=d)))^
      #(-1/2)*exp((-1/2)*t(A%*%X[t,i,]-psi_pa[t+1, 1:d])%*%
      # diag((psi_pa[t+1, (d+1):(d+d)]+1)^(-1), nrow=d,ncol=d)%*%
      #(A%*%X[t,i,]-psi_pa[t+1, 1:d])) 
    }
    
    
    #2. calculate psi_t
    #calculate min
    fn <- function(x, X, psi){
      lambda <-  sum(dmvnorm(X[t,1:N[l],],x[1:d],diag(x[(d+1):(d+d)], nrow=d,ncol=d))%*%psi[t,1:N[l]])/sum(psi[t,1:N[l]]^2)
      return(sum((psi[t,1:N[l]] - (1/lambda)*dmvnorm(X[t,1:N[l],],x[1:d],diag(x[(d+1):(d+d)], nrow=d,ncol=d)))^2))
      #sum_arg = 0						
      #for(i in 1:N[l]){						
      #sum_arg = sum_arg + (det(diag(2*pi*x[(d+1):(d+d)], nrow=d,ncol=d))^(-1/2)*						
      #exp((-1/2)*t(X[t,i,]-x[1:d])%*%diag((2*pi*x[(d+1):(d+d)])^(-1), nrow=d,ncol=d)%*%						
      #(X[t,i,]-x[1:d]))-x[2*d+1]*psi[t,i])^2						
      #}						
      #return(sum_arg)
    }
    
    #get the distribution of psi_t
    if(t == Time){
      psi_pa[t,] <- optim(par = c(colMeans(X[t,1:N[l],]), rep(1, d)),
                          fn = fn, X = X, psi = psi, method='L-BFGS-B',lower=c(rep(-Inf, d),rep(0, d)),upper=rep(Inf, 2*d))$par
    }else{
      psi_pa[t,] <- optim(par = c(X[t,which.max(psi[t,1:N[l]]),], rep(1, d)), 
                          fn = fn, X = X, psi = psi, method = 'L-BFGS-B',lower=c(rep(-Inf, d),rep(0, d)),upper=rep(Inf, 2*d))$par
    }#X[t,which.max(psi[t,1:N[l]]),]
    
    print(psi_pa[t, 1:d])
    print(obs[t,])
  }
  
  return(psi_pa)
  
}

####bpf####
obs <- Obs()
fkf.obj <- fkf(a0, P0, dt, ct, Tt, Zt, Ht, Gt, yt = t(obs))$logLik

index = 1
l = 1  #actually 0

#initialization l=1, run a psi^0 apf with N[1]
#psi_1...psi_T+1 = 1
#psi.tilda_0...psi.tilda_T = 1

#t = 1
#generate particles and weights

X[1,1:N[l],] <- rnorm(N[l]*d)  #particles
for(i in 1:N[l]){
  w[1,i] <- log(g(obs[1,], X[1,i,]))  #weights
}
Z[l] <- 0

#t=2:T
#2. conditional sample
re = 0
for(t in 2:Time){
  #print(t)
  
  #a)
  
  if(ESS(t,l,w, is.log = TRUE) <= kappa*N[l]){
    
    re = re + 1
    mx <- max(w[t-1, 1:N[l]])
    Z[l] <- Z[l] + log(mean(exp(w[t-1, 1:N[l]]-mx))) + mx
    w_ <- exp(w[t-1,1:N[l]]-mx)/sum(exp(w[t-1,1:N[l]] - mx))
    mix <- sample(1:N[l], N[l], replace = TRUE, prob = w_)
    X[t,1:N[l],] <- apply(X[t-1, mix,], 1, f )
    
    for(i in 1:N[l]){
      
      w[t,i] <- log(g(obs[t,], X[t,i,]))  
    }
    
  }else{
    
    #b)
    
    for(i in 1:N[l]){
      X[t,i,] <- f(X[t-1,i,])   
      w[t,i] <- w[t-1,i] + log(g(obs[t,], X[t,i,]))  
    }
  }
  
}

mx <- max(w[Time,1:N[l]])
Z[l] <- Z[l] + log(mean(exp(w[t,1:N[l]]-mx))) + mx
#print(paste0('re=',re))

####iapf####
while(index){
  
  #print(l)
  #a)
  output <- list()
  
  if(l != 1){
    output <-  APF(psi_pa, l, N)
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
    
    #print(paste0('Z[l]=',Z[l]))
    
    #print(paste0('Z=',fkf.obj))
    
    l <- l+1
  }else break
  
}

#3.
output <-  APF(psi_pa, l, N)
Z <- output[[4]]
Z_appro <- Z[l]
exp(Z_appro-fkf.obj)


