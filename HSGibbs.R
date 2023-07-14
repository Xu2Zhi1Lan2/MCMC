HSGibbs <- function(x,lambdasq0,eta0,tausq0,xi0,sigmasq0,burn_in,N){
  
  p = length(x) # dimension of parameters (theta)
  theta0 = numeric(p)
  theta.sample = matrix(0,p,N) # chain of theta
  lambdasq.sample = matrix(0,p,N) # chain of lambda square
  tausq.sample = numeric(N) # chain of tau square
  sigmasq.sample = numeric(N) # chain of sigma square
  
  ## burn-in
  
  for (b in 1:burn_in) {
    
    ## sample theta
    for (i in 1:p) {
      V = 1/(1/sigmasq0 + 1/(lambdasq0[i]^2*tausq0))
      E = V*x[i]/sigmasq0
      theta0[i] = rnorm(1,E,sqrt(V))
    }
    
    ## sample lambda and eta
    for (i in 1:p) {
      R.lambda = theta0[i]^2/(2*tausq0^2) + 1/eta0[i]
      lambdasq0[i] = 1/rgamma(1,1,R.lambda)
      R.eta = 1 + 1/lambdasq0[i]^2
      eta0[i] = 1/rgamma(1,1,R.eta)
    }
    
    ## sample tau and xi
    R.tau = 1/xi0 + 0.5*sum(theta0^2/lambdasq0^2)
    tausq0 = 1/rgamma(1,(1+0.5*p),R.tau)
    R.xi = 1 + 1/tausq0
    xi0 = 1/rgamma(1,1,R.xi)
    
    ## sample sigmasq
    R.sigmasq = 1 + 0.5*sum((x-theta0)^2)
    sigmasq0 = 1/rgamma(1,(1+0.5*p),R.sigmasq)
    
    if(b%%1000==0){print(b)}
  }
  
  ## MCMC
  
  for (n in 1:N) {
    ## sample theta
    for (i in 1:p) {
      V = 1/(1/sigmasq0 + 1/(lambdasq0[i]^2*tausq0))
      E = V*x[i]/sigmasq0
      theta0[i] = rnorm(1,E,sqrt(V))
      theta.sample[i,n] = theta0[i]
    }
    
    ## sample lambda and eta
    for (i in 1:p) {
      R.lambda = theta0[i]^2/(2*tausq0^2) + 1/eta0[i]
      lambdasq0[i] = 1/rgamma(1,1,R.lambda)
      lambdasq.sample[i,n] = lambdasq0[i]
      R.eta = 1 + 1/lambdasq0[i]^2
      eta0[i] = 1/rgamma(1,1,R.eta)
    }
    
    ## sample tau and xi
    R.tau = 1/xi0 + 0.5*sum(theta0^2/lambdasq0^2)
    tausq0 = 1/rgamma(1,(1+0.5*p),R.tau)
    tausq.sample[n] = tausq0
    R.xi = 1 + 1/tausq0
    xi0 = 1/rgamma(1,1,R.xi)
    
    ## sample sigmasq
    R.sigmasq = 1 + 0.5*sum((x-theta0)^2)
    sigmasq0 = 1/rgamma(1,(1+0.5*p),R.sigmasq)
    sigmasq.sample[n] = sigmasq0
  }
  
  MCMCsample = list(theta = theta.sample,
                    lambdasq = lambdasq.sample,
                    tausq = tausq.sample,
                    sigmasq = sigmasq.sample)
  return(MCMCsample)
}


