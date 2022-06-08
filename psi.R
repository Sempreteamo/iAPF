
for(t in Time:1){
  
  for(i in 1:N[l]){
    
   #1. set iterated psi_tilda^i  
    
    if(t == Time){
      psi[t,i] <- g(obs[t,],X[t,i,])
    }else{
      psi[t,i] <- g(obs[t,],X[t,i,])*
        (det(2*pi*(B+matrix(best[t+1, (d+1):(d+d^2)], nrow=d,ncol=d)))^
        (-1/2)*exp((-1/2)*t(A%*%X[t,i,]-best[t+1, 1:d])%*%solve(
          B+matrix(best[t+1, (d+1):(d+d^2)], nrow=d,ncol=d))%*%
            (A%*%X[t,i,]-best[t+1, 1:d]))+1/N[l])
    }
 
  }
  #2. calculate psi_t
  #calculate min
  
  
  fn <- function(x){
    sum_arg = 0
    for(i in 1:N[l]){
      sum_arg = sum_arg + (mvtnorm::dmvnorm(
        X[t,i,], mean = x[1:d], sigma = 
          diag(x[c(d+1,2*d+2,3*d+3,4*d+4,5*d+5)], nrow=d,ncol=d)) - 
          x[d+d^2+1]*psi[t,i])^2
    }
    return(sum_arg)
  }
  
  #get the distribution of psi_t
  best[t,] <- optim(par = c(rep(0, d), c(diag(1, nrow = d, ncol = d)), 1),
                    fn = fn, method = "L-BFGS-B")$par   
  

}
