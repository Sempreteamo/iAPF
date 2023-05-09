
seed.label<-as.numeric(Sys.getenv("SGE_TASK_ID"))
set.seed(123)
library(mvtnorm)
library(MASS)
library(FKF)
####settings####
final <- vector()
N <- vector()
N[1] <- 100
Time = 100
alpha = 0.42
d = 3
k <- 3
tau <- 0.5
kappa = 0.5

A <- matrix(nrow = d, ncol = d)
for (i in 1:d){
  for (j in 1:d){
    A[i,j] = alpha^(abs(i-j) + 1)
  }
}
B = C = D = diag(1, nrow = d, ncol = d)
X_true <-  matrix(0, nrow = Time, ncol = d )
obs <-  matrix(0, nrow = Time, ncol = d )

dt <- ct <- matrix(0,d,1)
Tt <- A
P0 <- Zt <- Ht <- Gt <- diag(1,d,d)
a0 <- rep(0,d)

index <- 1

Obs <- function(){
  X_true[1,] <- rnorm(d)     
  for(t in 2:Time){  #observations
    X_true[t,] <- rnorm(d) + A%*%X_true[t-1,]  #t(rnorm(d) + A%*%x)
  }
  return(matrix(rnorm(Time*d, X_true, 1), ncol = d))
}
obs <- Obs()

for(qq in 1:3){
set.seed(qq)

X <- array(NA, dim = c(Time, N[1], d))
X_ <- array(NA, dim = c(Time, N[1], d))
w <- matrix(NA, Time, N[1])
Z <- vector()

fkf.obj <- fkf(a0, P0, dt, ct, Tt, Zt, Ht, Gt, yt = t(obs))$logLik

f <- function(x){
  return (rnorm(d) + as.vector(A%*%x))   #trans prob
}

g <- function(y, x){  
  return ((2*pi)^(-d/2)*exp((-1/2)*t(y-x)%*%(y-x))) #obs prob  C%*%x = x
}

mu_aux <- function(psi_pa, l, N, t){  
  return(mvrnorm(N[l], mu =  as.vector(diag(((psi_pa[t, (d+1):(d+d)])^(-1)+1)^(-1), nrow=d,ncol=d)%*%
                                         (diag((psi_pa[t, (d+1):(d+d)])^(-1), nrow=d,ncol=d)%*%psi_pa[t,1:d])), 
                 Sigma = diag(((psi_pa[t, (d+1):(d+d)])^(-1)+1)^(-1), nrow=d,ncol=d)))
}

g_aux <- function(y, x, t, psi_pa){  
  if(t == 1){
    return(g(y, x)*psi_tilda(x, psi_pa, 1)*((2*pi)^(-d/2)*det(diag(psi_pa[t, (d+1):(d+d)]+1, nrow=d,ncol=d))^
                                              (-1/2)*exp((-1/2)*t(-psi_pa[t, 1:d])%*%
                                                           diag((psi_pa[t, (d+1):(d+d)]+1)^(-1), nrow=d,ncol=d)%*%
                                                           (-psi_pa[t, 1:d])))/psi_t(x, psi_pa, 1))  #g_1
    }else{
    return(g(y,x)*psi_tilda(x, psi_pa, t)/psi_t(x, psi_pa, t))  #g_2:T 
  }
}

f_aux <- function(x, psi_pa, t){
  return(mvrnorm(1, mu = as.vector(diag(((psi_pa[t, (d+1):(d+d)])^(-1)+1)^(-1), nrow=d,ncol=d)%*%
                                     (A%*%x + diag(psi_pa[t, (d+1):(d+d)]^(-1), nrow=d,ncol=d)%*%psi_pa[t,1:d])), 
                 Sigma = diag(((psi_pa[t, (d+1):(d+d)])^(-1)+1)^(-1), nrow=d,ncol=d)))  #f_2:T 
}#rnorm(N,(psi[t,2]^2 * a*X[t-1,anc] + sigma^2*psi[t,1])/(psi[t,2]^2+sigma^2), sd=sqrt(psi[t,2]^2 * sigma^2/ (psi[t,2]^2+sigma^2)))

ESS <- function(t,w, is.log=FALSE){
  if(is.log) {
    mx <- max(w[t-1,])
    s <- sum(exp(w[t-1,]-mx))
    ess <- 1/sum((exp(w[t-1,]-mx)/s)^2)
  }else{
    s <- sum(w[t-1,])
    ess <- 1/sum((w[t-1,]/s)^2) 
  }
  return(ess)  
}

Num <- function(Z, l, k){
  return(sd(exp(Z[max(l-k,1):l]-max(Z[max(l-k,1):l])))/mean(exp(Z[max(l-k,1):l]-max(Z[max(l-k,1):l]))))
}

psi_tilda <- function(x, psi_pa, t){  #from 0 to T. 0,T = 1  ?????
  if(t == Time){
    psi_tilda <- 1
  }else{   #psi_pa_t = psi_t #dnorm(A*x, psi_pa[t+1, 1], sqrt(psi_pa[t+1,2]^2+sigma^2)))
    psi_tilda <- (2*pi)^(-d/2)*det(diag(psi_pa[t+1, (d+1):(d+d)]+1, nrow=d, ncol=d))^(-1/2)*
      exp((-1/2)*t(as.vector(A%*%x) - psi_pa[t+1, 1:d])%*%diag((psi_pa[t+1, (d+1):(d+d)]+1)^(-1), nrow=d,ncol=d)%*%
            (as.vector(A%*%x)-psi_pa[t+1, 1:d]))  #f(xt, ÃÂt+1 )   #var should be psi^2+1?
  }
  return(psi_tilda)
}

psi_t <- function(x, psi_pa, t){ #from 1 to T+1. 1, T+1 = 1  
  if(t == Time + 1){
    psi_t <- 1
  }else{
    psi_t <-  (2*pi)^(-d/2)*det(diag(psi_pa[t, (d+1):(d+d)], nrow=d,ncol=d))^(-1/2)*						
      exp((-1/2)*t(x-psi_pa[t, 1:d])%*%diag((psi_pa[t, (d+1):(d+d)])^(-1), nrow=d,ncol=d)%*%						
            (x-psi_pa[t, 1:d]))
    #dnorm(x, psi_pa[t, 1], psi_pa[t, d+1])
  }
  
  return(psi_t)
}


####APF function####
APF <- function(psi_pa, l, Z, N){
  #l >= 2
  X <- array(NA, dim = c(Time, N[l], d))
  w <- matrix(NA, Time, N[l])
  Z[l] <- 0
  
  X[1,1:N[l],] <- mu_aux(psi_pa, l, N, t)  #particles
  
  for(i in 1:N[l]){
    w[1,i] <- log(g_aux(obs[1,], X[1,i,], 1, psi_pa)) #weights g(obs[1,], X[1,i,])*psi_tilda(X[1,i,], psi_pa, 2)  
  }
  
  for(t in 2:Time){
    
    if(ESS(t,w, is.log = TRUE) <= kappa*N[l]){
      mx <- max(w[t-1,])
      w_ <- exp(w[t-1,1:N[l]]-mx)/sum(exp(w[t-1, 1:N[l]] - mx))
      Z[l] = Z[l] + log(mean(exp(w[t-1,]-mx))) + mx
      mix <- sample(1:N[l],N[l], replace = TRUE, prob = w_)
      
      for(i in 1:N[l]){
        X[t,i,] <- f_aux(X[t-1, mix[i],], psi_pa, t)
        w[t,i] <- log(g_aux(obs[t,], X[t,i,], t, psi_pa))  
      }
    }else{
      
      for(i in 1:N[l]){
        
        X[t,i,] <- f_aux(X[t-1,i,],psi_pa, t) 
        w[t,i] <- w[t-1,i] + log(g_aux(obs[t,], X[t,i,],t, psi_pa))  
      }
    }
    
  }
  
  mx <- max(w[Time, 1:N[l]])
  Z[l] <- Z[l] + log(mean(exp(w[Time, 1:N[l]]-mx))) + mx
 
  return(list(obs, X, w, Z))
}

####psi function####
Psi <- function(l, X, N){
  psi <- matrix(NA, nrow = Time, ncol = N[l])
  psi_pa <- matrix(NA, nrow = Time, ncol = 2*d)
  
  for(t in Time:1){
    
    #print(t)
    
    if(t == Time){
      psi[t,1:N[l]] <- dmvnorm(X[t,1:N[l],], obs[t,])  
      
    }else{
      
      for(i in 1:N[l]){
        psi[t,i] <- g(obs[t,],X[t,i,])%*%dmvnorm(as.vector(A%*%X[t,i,]), 
                                                 psi_pa[t+1, 1:d], diag(psi_pa[t+1, (d+1):(d+d)]+1, nrow=d,ncol=d))
      }
    }
    
    fn <- function(x, X, psi){
      
      lambda <-  sum(dmvnorm(X[t,1:N[l],],x[1:d],
                             diag(x[(d+1):(d+d)], nrow=d,ncol=d))%*%psi[t,1:N[l]])/sum(psi[t,1:N[l]]^2)
      return( sum((psi[t,1:N[l]] - (1/lambda)*dmvnorm(X[t,1:N[l],],
                                                      x[1:d],diag(x[(d+1):(d+d)], nrow=d,ncol=d)))^2))
     
    }
    
    if(t == Time){
      psi_pa[t,] <- optim(par = c(colMeans(X[t,1:N[l],]), rep(1, d)),
                          fn = fn, X = X, psi = psi, method = 'L-BFGS-B',lower = c(rep(-Inf, d), rep(0.1, d)),upper = rep(Inf, 2*d))$par
    }else{
      psi_pa[t,] <- optim(par = c(X[t,which.max(psi[t,1:N[l]]),], rep(1, d)), 
                          fn = fn, X = X, psi = psi, method='L-BFGS-B',lower=c(rep(-Inf, d),rep(0.1, d)), upper=rep(Inf, 2*d))$par
    }
    
    print(psi_pa[t, 1:d])
    print(obs[t,])
  }
  
  return(psi_pa)
  
}


####Init####

l = 1  #actually 0


X[1,1:N[l],] <- rnorm(N[l]*d)  #particles  
for(i in 1:N[l]){
  w[1,i] <- log(g(obs[1,], X[1,i,]))  
}
Z[1] <- 0

for(t in 2:Time){

  if(ESS(t,w, is.log=TRUE) <= kappa*N[l]){
    
    mx <- max(w[t-1,])
    w_ <- exp(w[t-1,]-mx)/sum(exp(w[t-1,] - mx))
    Z[1] = Z[1] + log(mean(exp(w[t-1,] - mx))) + mx
    mix <- sample(1:N[l],N[l], replace = TRUE, prob = w_)
    
    for(i in 1:N[l]){
      X[t,i,] <- f(X[t-1,mix[i],]) 
      w[t,i] <- log(g(obs[t,], X[t,i,]))  
    }
    
  }else{
    
    for(i in 1:N[l]){
      X[t,i,] <- f(X[t-1,i,])   
      w[t,i] <- w[t-1,i] + log(g(obs[t,], X[t,i,]))  
    }
  }
  
}

mx <- max(w[Time, 1:N[l]])
Z[1] <- Z[1] + log(mean(exp(w[Time, 1:N[l]]-mx))) + mx

####iAPF####
while(TRUE){

  output <- list()
  
  if(l != 1){
    output <-  APF(psi_pa, l, Z, N)
    X <- output[[2]]
    w <- output[[3]]
    Z <- output[[4]]
  }
  
  #b)
  
  if(l <= k | (Num(Z, l, k) >= tau)){
     
    psi_pa <- Psi(l, X, N) 
    
    
    if(l > k & N[max(l-k,1)] == N[l] & is.unsorted(Z[max(l-k,1):l])){
      N[l+1] <- 2*N[l]
      
    }else{
      N[l+1] <- N[l]
    }
  
    l <- l+1
  }else break
}

output <-  APF(psi_pa, l, Z, N)
Z <- output[[4]]
Z_appro <- Z[l]
print(exp(Z_appro-fkf.obj))
final[index] <- exp(Z_appro-fkf.obj)
index = index + 1
}
