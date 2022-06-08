
library(mvtnorm)
library(MASS)

#settings
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
cov = B = C = D = diag(1, nrow = d, ncol = d)
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

f <- function(x){
  return (mvrnorm(n = 1, mu = A%*%x, Sigma = B)) #trans prob
}

g <- function(x, state){
  return (mvtnorm::dmvnorm(x, mean = C%*%state, sigma = D)) #obs prob
}

psi <- matrix(NA, nrow = Time, ncol = 20000) #iterated psi to decide psi_t for each l

g_aux <- function(x, state, t, best, l){   
  if(t == 1){
    return(g(x, state)*psi_tilda(state, best, 2))  #g_1
  }else{
    return(g(x, state)*psi_tilda(state, best, t+1)/psi_t(state, best, t, l))  #g_2:T 
  }
}

f_aux <- function(x, best, t, l){   
  return(f(x)*psi_t(x, best, t, l)/psi_tilda(x, best, t))  #f_2:T 
}

ESS <- function(t,l){
  return(sum(w[t-1,1:N[l]])^2/sum(w[t-1,1:N[l]]^2))
}

best <- matrix(NA, nrow = Time, ncol = d+d*d+1)

#1. Initialize
#l=0  large outer iteration
#psi_t: 1-T+1  psi_T+1 =1
#psi_tilda: 0-T  index +1

#use psi_t to calulcate psi_tilda[t]

psi_tilda <- function(x, best, t){  #from 0 to T. 0,T = 1
  if(t == (Time + 1) | 1){
    psi_tilda <- 1
  }else{   #best_t = psi_t
    psi_tilda <- det(2*pi*(B+matrix(best[t, (d+1):(d+d^2)], nrow=d,ncol=d)))^
      (-1/2)*exp((-1/2)*t(A%*%x-best[t, 1:d])%*%solve(
        B+matrix(best[t, (d+1):(d+d^2)], nrow=d,ncol=d))%*%
          (A%*%x-best[t, 1:d]))  #f(xt, ψt+1 )   #need +1??
  }
  
  return(psi_tilda)
}

psi_t <- function(x, best, t, l){ #from 1 to T+1. 1, T+1 = 1  
  if(t == (Time + 1) | 1){
    psi_t <- 1
  }else{
    psi_t <- mvtnorm::dmvnorm(x, mean = best[t, 1:d], sigma = matrix(best[t, (d+1):(d+d^2)], nrow=d,ncol=d)) + 1/N[l]
  }
  
  return(psi_t)
}

#2. repeat
index = 1
Sum <- matrix(NA, nrow = 20000, ncol = 5)
Z_true <- vector()
l = 1  #actually 0

#initialization l=1, run a psi^0 apf with N[1]
#psi_1...psi_T+1 = 1
#psi.tilda_0...psi.tilda_T = 1

X_true[1,] <- mvrnorm(n = 1, mu = m, Sigma = cov)
obs[1,] <- mvrnorm(n = 1, mu = C%*%X_true[1,], Sigma = D)

for(t in 2:Time){  #observations
  X_true[t,] <- f(X_true[t-1,])
  obs[t,] <- mvrnorm(n = 1, mu = C%*%X_true[t,], Sigma = D)
}

#t = 1
#generate particles and weights

Fun <- function(state){
  return(g(obs[1,], state))
}
X[1,,] <- mvrnorm(n = N[l], mu = m, Sigma = cov)  #particles
w[1,] <- apply(X[1,,],1,Fun)  #weights

#approximation
Z[l] <-sum(w[1,1:N[l]])/N[l]  

#t=2:T
#2. conditional sample

for(t in 2:Time){
  
  #a)
  
  if(ESS(t,l) <= kappa*N[l]){
    
    for(i in 1:N[l]){
      
      for(j in 1:N[l]){
        Sum[j,] <- w[t-1,j]*f(X[t-1,j,])
      }
      X[t,i,] <- colSums(Sum[1:N[l],])/sum(w[t-1,1:N[l]])  
    }
    w[t,] <- apply(X[t,,],1,Fun)
    
    #approximation
    Z[l] <- Z[l]*sum(w[t,1:N[l]])/N[l] 
    
  }else{
    
    #b)
    X[t,,] <- apply(X[t-1,,],1,f)
    for(i in 1:N[l]){
      w[t,i] <- w[t-1,i]*g(obs[t,], X[t,i,])  
    }
  }
  
}
Z[l] <- Z[l]*sum(w[Time,1:N[l]])/N[l]

while(index){
  #a)
  if(l != 1){
    source("APF.R")
  }
  
  print(l)
  
  #b)

  if(l <= k | (sd(Z[max(l-k,1):l])/mean(Z[max(l-k,1):l]) >= tau)){
    #psi^{l+1}
    source("psi.R") 
    
    if(l > k & N[max(l-k,1)] == N[l] & is.unsorted(Z[max(l-k,1):l])){
      N[l+1] <- 2*N[l]
    }else{
      N[l+1] <- N[l]
    }
    
    print(Z[l])
    print(Z[l]/mean(Z_true))
          
    l <- l+1
  }else break
}

#3.
Z_appro <- Z[l]
