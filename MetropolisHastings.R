### model
### Likelihood: x \sim N(theta,1)
### Prior: theta \sim Laplace(0,1)

### Independence Chain Metropolis-Hastings
IndMH <- function(theta0,x,burnin,N){
  n = length(x)
  theta.indchain = numeric(N)
  
  ## burn-in
  for (i in 1:burnin) {
    theta.star = rnorm(1,0,2) # proposal
    MH1 = -0.5*sum((x-theta.star)^2)-0.2*abs(theta.star)-theta0^2/8
    MH2 = -0.5*sum((x-theta0)^2)-0.2*abs(theta0)-theta.star^2/8
    alpha = min(1,exp(MH1-MH2)) # acceptance rate
    u = runif(1,0,1)
    if(u<=alpha){theta0 = theta.star}
  }
  
  ## sampling
  for (i in 1:N) {
    theta.star = rnorm(1,0,2) # proposal
    MH1 = -0.5*sum((x-theta.star)^2)-0.2*abs(theta.star)-theta0^2/8
    MH2 = -0.5*sum((x-theta0)^2)-0.2*abs(theta0)-theta.star^2/8
    alpha = min(1,exp(MH1-MH2)) # acceptance rate
    u = runif(1,0,1)
    if(u<=alpha){theta.indchain[i] = theta.star}
    else{theta.indchain[i] = theta0}
    theta0 = theta.indchain[i]
  }
  
  return(theta.indchain)
}

### Random Walk Metropolis-Hastings
RWMH <- function(theta0,x,burnin,N){
  n = length(x)
  theta.rwchain = numeric(N)
  
  ## burn-in
  for (i in 1:burnin) {
    theta.star = theta0 + rnorm(1,0,1) # proposal
    MH1 = -0.5*sum((x-theta.star)^2)-0.2*abs(theta.star)
    MH2 = -0.5*sum((x-theta0)^2)-0.2*abs(theta0)
    alpha = min(1,exp(MH1-MH2)) # acceptance rate
    u = runif(1,0,1)
    if(u<=alpha){theta0 = theta.star}
  }
  
  ## sampling
  for (i in 1:N) {
    theta.star = theta0 + rnorm(1,0,1) # proposal
    MH1 = -0.5*sum((x-theta.star)^2)-0.2*abs(theta.star)
    MH2 = -0.5*sum((x-theta0)^2)-0.2*abs(theta0)
    alpha = min(1,exp(MH1-MH2)) # acceptance rate
    u = runif(1,0,1)
    if(u<=alpha){theta.rwchain[i] = theta.star}
    else{theta.rwchain[i] = theta0}
    theta0 = theta.rwchain[i]
  }
  
  return(theta.indchain)
}

###### Simple Comparison
x = rbinom(10,1,0.2)*rnorm(10,0,4) + rnorm(10,0,1)
theta.indchain = IndMH(1,x,5000,10000)
theta.rwchain = RWMH(1,x,5000,10000)
plot(theta.indchain, type = "l")
plot(theta.rwchain, type = "l")


