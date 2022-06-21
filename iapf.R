library(mvtnorm)
library(MASS)
library(profvis)

####settings####
set.seed(123)
alpha <- 0.42
d <- 5
k <- 5
kappa <- 0.5
tau <- 0.5
m <- rep(0, d)
N <- vector()
N[1] <- 1000
Time <- 100
cov = B = C = D = diag(rep(1,5), nrow = d, ncol = d)
A <- matrix(nrow = d, ncol = d)
for (i in 1:d){
  for (j in 1:d){
    A[i,j] = alpha^(abs(i-j) + 1)
  }
}
X <- array(0, dim = c(Time, 20000, d)) #particles t horizontal; N vertical, d
X_true <- matrix(0, nrow = Time, ncol = 5 )
obs <- matrix(0, nrow = Time, ncol = 5 )
w <- matrix(NA, nrow = Time, ncol = 20000 ) #weight
Z <- vector()  #approximation

# FKF

dt <- ct <- matrix(0,5,1)
Tt <- A
P0 <- Zt <- Ht <- Gt <- diag(1,5,5)
a0 <- rep(0,d)

f <- function(x){
  return (rnorm(d) + as.vector(A%*%x))   #trans prob
}

g <- function(x, state){  
  return (0.01010533*exp((-1/2)*t(x-state)%*%(x-state))) #obs prob  C%*%state = state
  #det(diag(2*pi, nrow = d, ncol = d))^(-1/2) = 0.01010533
}

psi <- matrix(NA, nrow = Time, ncol = 20000) #iterated psi to decide psi_t for each l

g_aux <- function(x, state, t, best, l){   
  if(t == 1){
    return(g(x, state)*psi_tilda(state, best, 1, l)
           )  #g_1
  }else{
    return(g(x, state)*psi_tilda(state, best, t, l)/psi_t(state, best, t, l))  #g_2:T 
  }
}


# this is Gaussian above how to simplify? can z and z diminish?
f_aux <- function(x, best, t){   
  return(mvrnorm(1, mu = diag(((best[t, (d+1):(d+5)])^(-1)+1)^(-1), nrow=d,ncol=d)%*%
                   (diag((best[t, (d+1):(d+5)])^(-1), nrow=d,ncol=d)%*%best[t,1:d]+A%*%x), 
                 Sigma = diag(((best[t, (d+1):(d+5)])^(-1)+1)^(-1), nrow=d,ncol=d)))  #f_2:T 
}

ESS <- function(t,l,w){
  return(sum(w[t-1,1:N[l]])^2/sum(w[t-1,1:N[l]]^2))
}

best <- matrix(NA, nrow = Time, ncol = 2*d+1)

Num <- function(Z, l, k){
  return(sd(exp(Z[max(l-k,1):l]-max(Z)))/mean(exp(Z[max(l-k,1):l]-max(Z))))
  }

#1. Initialize
#l=0  large outer iteration
#psi_t: 1-T+1  psi_T+1 =1
#psi_tilda: 0-T  index +1

#use psi_t to calulcate psi_tilda[t]

psi_tilda <- function(x, best, t, l){  #from 0 to T. 0,T = 1
  if(t == Time){
    psi_tilda <- 1
  }else{   #best_t = psi_t
    psi_tilda <- det(2*pi*(diag(best[t, (d+1):(d+5)]+1, nrow=d,ncol=d)))^
      (-1/2)*exp((-1/2)*t(A%*%x-best[t, 1:d])%*%
        diag((best[t, (d+1):(d+5)]+1)^(-1), nrow=d,ncol=d)%*%
          (A%*%x-best[t, 1:d])) + 1/N[l]  #f(xt, Ït+1 )   #need +1??
  }
  
  return(psi_tilda)
}

psi_t <- function(x, best, t, l){ #from 1 to T+1. 1, T+1 = 1  
  if(t == (Time + 1) | 1){
    psi_t <- 1
  }else{
    psi_t <- det(diag(2*pi*best[t, (d+1):(d+5)], nrow=d,ncol=d))^(-1/2)*						
      exp((-1/2)*t(x-best[t, 1:d])%*%diag((2*pi*best[t, (d+1):(d+5)])^(-1), nrow=d,ncol=d)%*%						
            (x-best[t, 1:d])) + 1/N[l]
  }
  
  
  return(psi_t)
}

####APF function####
APF <- function(best, l){
  #l >= 2
  
  Z[l] <- 0
  X_true[1,] <- rnorm(d) 
  for(t in 2:Time){  #observations
    X_true[t,] <- f(X_true[t-1,])   #t(rnorm(d) + A%*%x)
  }
  obs <- rnorm(d) + X_true
  
  X[1,1:N[l],] <- rnorm(N[l]*d)   #particles
  
  for(i in 1:N[l]){
    w[1,i] <- g_aux(obs[1,], X[1,i,],1, best, l) #weights g(obs[1,], X[1,i,])*psi_tilda(X[1,i,], best, 2)  
  }
  
  #t=2:T
  #2. conditional sample
  for(t in 2:Time){
    
    print(t)
    
    #a)
    
    if(ESS(t,l,w) <= kappa*N[l]){
      
      w_ <- w[t-1,1:N[l]]/sum(w[t-1,1:N[l]])   #each t
      mix <- sample(1:N[l],N[l], replace = TRUE, prob = w_)
      X[t,1:N[l],] <- apply(X[t-1, mix,], 1, function(x) rnorm(d)%*%
                              diag(((best[t, (d+1):(d+5)])^(-1)+1)^(-1), nrow=d,ncol=d) + 
                              as.vector(diag(((best[t, (d+1):(d+5)])^(-1)+1)^(-1), nrow=d, ncol=d)%*%
                                          (diag((best[t, (d+1):(d+5)])^(-1), nrow=d,ncol=d)%*%best[t,1:d]+A%*%x)))
      
      for(i in 1:N[l]){
        
        w[t,i] <- g_aux(obs[t,], X[t,i,], t, best, l)  
      }
      
    }else{
      
      #b)
      for(i in 1:N[l]){
        
        X[t,i,] <- f_aux(X[t-1,i,],best, t) 
        w[t,i] <- w[t-1,i]*g_aux(obs[t,], X[t,i,],t, best, l)  
      }
    }
    
  }
  Z[l] = 0
  for(t in 1:Time){
    Z[l] = Z[l] + log(mean(apply(X[t,1:N[l],],1,function(state) 
      0.01010533*exp((-1/2)*t(obs[t,]-state)%*%(obs[t,]-state)))))
  }
  
  fkf.obj <- -fkf(a0, P0, dt, ct, Tt, Zt, Ht,
                  Gt, yt = t(obs))$logLik
  
  return(list(obs, X, w, Z, fkf.obj))
}

####psi function####
Psi <- function(l, obs, state, X){
  
  
  for(t in Time:1){
    
    print(t)
    
    if(t == Time){
      psi[t,1:N[l]] <- apply(X[t,1:N[l],],1,function(state) 
        0.01010533*exp((-1/2)*t(obs[t,]-state)%*%(obs[t,]-state)))
    }else{
      
      for(i in 1:N[l]){
        psi[t,i] <- g(obs[t,],X[t,i,])*
          (det(2*pi*(B+diag(best[t+1, (d+1):(d+5)], nrow=d,ncol=d)))^
             (-1/2)*exp((-1/2)*t(A%*%X[t,i,]-best[t+1, 1:d])%*%solve(
               B+diag(best[t+1, (d+1):(d+5)], nrow=d,ncol=d))%*%
                 (A%*%X[t,i,]-best[t+1, 1:d]))+1/N[l])
      }
    }
    
    
    #2. calculate psi_t
    #calculate min
    
    
    fn <- function(x, X, psi){
      sum_arg = 0						
      for(i in 1:N[l]){						
        sum_arg = sum_arg + (det(diag(2*pi*x[(d+1):(d+5)], nrow=d,ncol=d))^(-1/2)*						
                               exp((-1/2)*t(X[t,i,]-x[1:d])%*%diag((2*pi*x[(d+1):(d+5)])^(-1), nrow=d,ncol=d)%*%						
                                     (X[t,i,]-x[1:d]))-x[2*d+1]*psi[t,i])^2						
      }						
      return(sum_arg)
    }
    #get the distribution of psi_t
    if(t == Time){
      best[t,] <- optim(par = c(rep(0, d), rep(10, d), 3),
                        fn = fn, X = X, psi = psi, method = "BFGS")$par
    }else{
      best[t,] <- optim(par = c(best[t+1,1:d], best[t+1,(d+1):(d+5)], best[t+1,ncol(best)]),
                        fn = fn, X = X, psi = psi, method = "BFGS")$par
    }
     
    
  }
 
  return(best)
  
}

####2. repeat####
index = 1
Sum <- matrix(NA, nrow = 20000, ncol = 5)
l = 1  #actually 0
Z[l] = 0

#initialization l=1, run a psi^0 apf with N[1]
#psi_1...psi_T+1 = 1
#psi.tilda_0...psi.tilda_T = 1

X_true[1,] <- rnorm(d) + m    #real writing!
for(t in 2:Time){  #observations
  X_true[t,] <- f(X_true[t-1,])   #t(rnorm(d) + A%*%x)
}
obs <- rnorm(d) + X_true


#t = 1
#generate particles and weights

X[1,1:N[l],] <- rnorm(N[l]*d)  #particles
for(i in 1:N[l]){
  w[1,i] <- g(obs[1,], X[1,i,])  #weights
}

a=0
#t=2:T
#2. conditional sample

for(t in 2:Time){
  print(t)
  
  #a)
  
  if(ESS(t,l,w) <= kappa*N[l]){
    
    w_ <- w[t-1,1:N[l]]/sum(w[t-1,1:N[l]])
    mix <- sample(1:N[l],N[l], replace = TRUE, prob = w_)
    X[t,1:N[l],] <- apply(X[t-1, mix,], 1, f )
    
    w[t, 1:N[l]] <- apply(X[t, 1:N[l],], 1, function(state) 0.01010533*exp((-1/2)*t(obs[t,]-state)%*%(obs[t,]-state)))
      
    a=a+1
    
  }else{
    
    #b)
   
    for(i in 1:N[l]){
      X[t,i,] <- f(X[t-1,i,])   
      w[t,i] <- w[t-1,i]*g(obs[t,], X[t,i,])  
    }
  }
  
}
for(t in 1:Time){
  Z[l] = Z[l] + log(mean(apply(X[t,1:N[l],],1,function(state) 0.01010533*exp((-1/2)*t(obs[t,]-state)%*%(obs[t,]-state)))))
}
print(Z[l])

fkf.obj <- -fkf(a0, P0, dt, ct, Tt, Zt, Ht, Gt, yt = t(obs))$logLik
print(exp(Z[l] + fkf.obj))

while(index){
  
  print(l)
  #a)
  output <- list()
  
  if(l != 1){
    output <-  APF(best, l)
    obs <-output[[1]]
    X <- output[[2]]
    w <- output[[3]]
    Z <- output[[4]]
    fkf.obj <- output[[5]]
  }
  
  
  
  #b)

  if(l <= k | (Num(Z, l, k) >= tau)){
    #psi^{l+1}
    best <- Psi(l, obs, state, X) 
    
    if(l > k & N[max(l-k,1)] == N[l] & is.unsorted(Z[max(l-k,1):l])){
      N[l+1] <- 2*N[l]
      
    }else{
      N[l+1] <- N[l]
    }
    
    print(Z[l])   #Z[l]
    
    print(exp(Z[l] + fkf.obj))
    
    l <- l+1
  }else break
}

#3.
Z_appro <- Z[l]
