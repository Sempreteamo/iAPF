
#l >= 2
X_true[1,] <- mvrnorm(n = 1, mu = m, Sigma = cov)
obs[1,] <- mvrnorm(n = 1, mu = C%*%X_true[1,], Sigma = D)

for(t in 2:Time){  #observations
  X_true[t,] <- f(X_true[t-1,])
  obs[t,] <- mvrnorm(n = 1, mu = C%*%X_true[t,], Sigma = D)
}

X[1,,] <- mvrnorm(n = N[l], mu = m, Sigma = cov)   #particles
for(i in 1:N[l]){
  w[1,i] <- g_aux(obs[1,], X[1,i,],1, best, l) #weights   
}

Z[l] <-sum(w[1,1:N[l]])/N[l]

#t=2:T
#2. conditional sample
for(t in 2:Time){

  #a)
  
  if(ESS(t,l) <= kappa*N[l]){
    
    for(i in 1:N[l]){
      
      for(j in 1:N[l]){
        Sum[j,] <- w[t-1,j]*f_aux(X[t-1,j,],best, t, l)/sum(w[t-1,1:N[l]])  
      }
      X[t,i,] <- colSums(Sum[1:N[l],])
      w[t,i] <- g_aux(obs[t,], X[t,i,], t, best, l)  
    }
    
    Z[l] <- Z[l]*sum(w[t,1:N[l]])/N[l] 
    
  }else{
    
    #b)
    for(i in 1:N[l]){
      X[t,i,] <- f_aux(X[t-1,i,],best, t, l) 
      w[t,i] <- w[t-1,i]*g_aux(obs[t,], X[t,i,],t, best, l)  
    }
  }
  
}
Z[l] <- Z[l]*sum(w[Time,1:N[l]])/N[l]