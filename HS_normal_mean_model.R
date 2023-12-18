################################################################
### functions to implement the sampling algorithm for HS
################################################################

## Adaptive Metropolis-Hastings for lambda
MH.lambda <- function(theta0,lambda0,tau0,ss_lambda){
  lambda.star = lambda0+rnorm(1,0,ss_lambda)
  T0 = log(1+exp(-1000*lambda0))- log(1+exp(-1000*lambda.star)) # sigmoid approximation to force lambda to be positive
  T1 = log(1 + lambda0^2) - log(1 + lambda.star^2) # prior ratio
  T2 = log(abs(lambda0)) - log(abs(lambda.star)) + 0.5*(theta0^2/(tau0^2))*(1/lambda0^2 - 1/lambda.star^2) # likelihood ratio
  MH = min(1,exp(T0+T1+T2))
  U = runif(1,0,1)
  if(U<=MH){return(list(lambda = lambda.star,acc = 1))}
  else{return(list(lambda = lambda0, acc = 0))}
}

## Adaptive Metropolis-Hastings for tau
MH.tau <- function(theta0,lambda0,tau0,ss_tau){
  n = length(theta0)
  tau.star = tau0+rnorm(1,0,ss_tau)
  T0 = log(1+exp(-1000*tau0)) - log(1+exp(-1000*tau.star)) # sigmoid approximation to force tau to be positive
  T1 = log(1/n^2 + tau0^2) - log(1/n^2 + tau.star^2) # prior ratio
  T2 = n*(log(abs(tau0))-log(abs(tau.star))) + 0.5*(1/tau0^2 - 1/tau.star^2)*sum(theta0^2/(lambda0^2)) # likelihood ratio
  MH = min(1,exp(T0+T1+T2))
  U = runif(1,0,1)
  if(U<=MH){return(list(tau = tau.star,acc = 1))}
  else{return(list(tau = tau0, acc = 0))}
}

## Main function
HS <- function(y,theta0,lambda0,tau0,sigmasq0,batch_size,max_batch,burn_in,nmc,n_thin){
  ## outputs
  n = length(y)
  J = rep(1,n)
  theta.nmc = matrix(0,n,nmc/n_thin)
  lambda.nmc = matrix(0,n,nmc/n_thin)
  kappa.nmc = matrix(0,n,nmc/n_thin)
  tau.nmc = rep(0,nmc/n_thin)
  sigmasq.nmc = rep(0,nmc/n_thin)
  
  ## tuning procedure
  ss_tau = 1
  ss_lambda = rep(1,n)
  
  if (max_batch > 0) {
    cat("tuning begins...\n")
    for (index_batch in 1:max_batch) {
      
      ar_tau = 0
      ar_lambda = rep(0,n)
      
      for (index_iter in 1:batch_size) {
        
        A = lambda0^2*tau0^2
        theta0 = rnorm(n,y*A/(1+A),sqrt(A*sigmasq0/(1+A)))
        
        for (j in 1:n) {
          lambda_res <- MH.lambda(theta0[j],lambda0[j],tau0,ss_lambda[j])
          lambda0[j] <- lambda_res$lambda
          ar_lambda[j] <- ar_lambda[j] + lambda_res$acc
        }
        
        tau_res <- MH.tau(theta0,lambda0,tau0,ss_tau)
        tau0 <- tau_res$tau
        ar_tau <- ar_tau + tau_res$acc
      }
      
      ar_lambda <- ar_lambda / batch_size
      for (j in 1:n) {
        if (ar_lambda[j] < 0.35) ss_lambda[j] <- ss_lambda[j] * 0.9
        else if (ar_lambda[j] > 0.45) ss_lambda[j] <- ss_lambda[j] * 1.1
      }
      
      ar_tau <- ar_tau / batch_size
      if (ar_tau < 0.35) ss_tau <- ss_tau * 0.9
      else if (ar_tau > 0.45) ss_tau <- ss_tau * 1.1
      
      sigmasq0 = 1/rgamma(1,0.5*n,0.5*sum((y-theta0)^2))
      
      cat("batch", index_batch,"ar_lambda1:", ar_lambda[1], "ss_lambda1:", ss_lambda[1], "ss_tau:", ss_tau, "ar_tau:", ar_tau, "\n")
      
    }
    cat("tuning ends\n")
  }
  
  ## burn-in procedure
  for(i in 1:burn_in){
    A = lambda0^2*tau0^2
    theta0 = rnorm(n,y*A/(1+A),sqrt(A*sigmasq0/(1+A)))
    for (j in 1:n) {
      lambda_res <- MH.lambda(theta0[j],lambda0[j],tau0,ss_lambda[j])
      lambda0[j] <- lambda_res$lambda
    }
    
    tau_res <- MH.tau(theta0,lambda0,tau0,ss_tau)
    tau0 <- tau_res$tau
    
    sigmasq0 = 1/rgamma(1,0.5*n,0.5*sum((y-theta0)^2))
    
    if(i%%1000==0){print(i)}
  }
  
  ## sampling procedure
  for (i in 1:nmc) {
    A = lambda0^2*tau0^2
    theta0 = rnorm(n,y*A/(1+A),sqrt(A*sigmasq0/(1+A)))
    
    for (j in 1:n) {
      lambda_res <- MH.lambda(theta0[j],lambda0[j],tau0,ss_lambda[j])
      lambda0[j] <- lambda_res$lambda
    }
    
    tau_res <- MH.tau(theta0,lambda0,tau0,ss_tau)
    tau0 <- tau_res$tau
    
    sigmasq0 = 1/rgamma(1,0.5*n,0.5*sum((y-theta0)^2))
    
    if(i%%1000==0){print(i)}
    
    if(i%%n_thin==0){
      theta.nmc[,i/n_thin] = theta0
      lambda.nmc[,i/n_thin] = lambda0
      kappa.nmc[,i/n_thin] = 1/(1+lambda0^2*tau0^2)
      tau.nmc[i/n_thin] = tau0
      sigmasq.nmc[i/n_thin] = sigmasq0
    }
  }
  
  ## return
  samples = list(theta = theta.nmc, lambda = lambda.nmc,
                 kappa = kappa.nmc, tau = tau.nmc, sigmasq = sigmasq.nmc)
}


set.seed(2023)
y = rnorm(200,0,3.26)*rbinom(200,1,0.2) + rnorm(200,0,1)
batch_size = 200
max_batch = 20
burn_in = 6000
nmc = 20000
n_thin = 10
theta0 = rep(0.2,200)
lambda0 = rep(1,200)
tau0 = 0.1
sigmasq0 = 1
HS.MCMC = HS(y,theta0,lambda0,tau0,sigmasq0,batch_size,max_batch,burn_in,nmc,n_thin)

plot(HS.MCMC$sigmasq,type = "l")
plot(rowMeans(HS.MCMC$theta)~y)


## repeated measure
set.seed(2023)
M = 20
Type = rbinom(200,1,0.2)
Signals = rnorm(200,0,3.26)*Type
y_bar = Signals + rnorm(200,0,1/sqrt(M))
batch_size = 200
max_batch = 20
burn_in = 6000
nmc = 20000
n_thin = 10
theta0 = rep(0.2,200)
lambda0 = rep(1,200)
tau0 = 0.1
sigmasq0 = 1/M
HS.MCMC = HS(y,theta0,lambda0,tau0,sigmasq0,batch_size,max_batch,burn_in,nmc,n_thin)

plot(HS.MCMC$sigmasq/M,type = "l")
plot(rowMeans(HS.MCMC$theta)~y)
