## Model: x_i \sim N(\theta_i,1)
## prior: \theta_i \sim N(0,\sigma_i^2), \sigma_i^2 \sim Gamma(0.1,0.01)
## Try Auto-tuning HMC to estimate this model

################################### Functions

### Data Generating
DGP <- function(n,p,sigma){
  Signal = rbinom(n,1,p)
  psi = sqrt(2*log(n))
  theta = Signal*rnorm(n,0,psi)
  x = theta + rnorm(n,0,sigma)
  simdata = list(x = x, theta = theta, Signal = Signal)
  return(simdata)
}

### direct sampling for theta
sampling.theta <- function(x,sigmasq0){
  n = length(sigmasq0)
  E = sigmasq0*x/(1+sigmasq0)
  S = sqrt(sigmasq0/(1+sigmasq0))
  theta = numeric(n)
  for (i in 1:n) {
    theta[i] = rnorm(1,E[i],S[i])
  }
  return(theta)
}

### Auto-tuning HMC for sigmasq
U <- function(theta0,sigmasq0){
  return(-1.4*log(abs(sigmasq0))-0.01*sigmasq0-0.5*theta0^2/sigmasq0)
}

Hamiltonian <- function(theta0,sigmasq0,p,m){
  return(U(theta0,sigmasq0)+0.5*m*p^2) # m is the inverse of mass matrix
}

grad_U <- function(theta0,sigmasq0){
  return(-1.4/sigmasq0-0.01+0.5*theta0^2/sigmasq0^2)
}

Hessian_U <- function(theta0,sigmasq0){
  return(1.4/sigmasq0^2-theta0^2/sigmasq0^3)
}

sampling.sigmasq <- function(theta0,sigmasq0,L = 51){
  n = length(theta0)
  sigmasq = numeric(n)
  for (i in 1:n) {
    m0 = 1/Hessian_U(theta0[i],sigmasq0[i])
    epsilon0 = 1
    q = numeric(L);q[1] = sigmasq0[i]
    p0 = rnorm(1,0,1)
    p = numeric(L);p[1] = p0 - 0.5*epsilon0*grad_U(theta0[i],q[1])
    m = numeric(L);m[1] = m0
    for (j in 2:L) {
      q[j] = q[j-1] + epsilon0*m[j-1]*p[j-1]
      if(q[j]<=0){p[j] = -p[j-1]}
      else{p[j] = p[j-1] - epsilon0*grad_U(theta0[i],q[j])}
      if(abs(q[j]-q[j-1])<=0.1){m[j] = m[j-1]}
      else{m[j] = m[j-1] + (grad_U(theta0[i],q[j]) - grad_U(theta0[i],q[j-1]) - m[j-1]*(q[j]-q[j-1]))/(q[j]-q[j-1])}
      if(j==L){
        q[j] = q[j-1] + epsilon0*m[j-1]*p[j-1]
        p[j] = p[j-1] - 0.5*epsilon0*grad_U(theta0[i],q[j])
        if(abs(q[j]-q[j-1])<=0.1){m[j] = m[j-1]}
        else{m[j] = m[j-1] + (grad_U(theta0[i],q[j]) - grad_U(theta0[i],q[j-1]) - m[j-1]*(q[j]-q[j-1]))/(q[j]-q[j-1])} }
    }
    sigmasq.star = q[L]
    p.star = -p[L]
    m1 = m[L]
    H0 = Hamiltonian(theta0[i],sigmasq0[i],p0,m0)
    Hstar = Hamiltonian(theta0[i],sigmasq.star,p.star,m1)
    accept = exp(H0-Hstar)
    U = runif(1,0,1)
    if(U<=accept){sigmasq[i] = sigmasq.star}
    else{sigmasq[i] = sigmasq0[i]}
  }
  return(sigmasq)
}

### main function
ATHMC <- function(x,burn_in,nmc,theta0,sigmasq0){
  n = length(x)
  
  theta.nmc = matrix(0,n,nmc)
  sigmasq.nmc = matrix(0,n,nmc)
  
  # burn-in procedure
  for (k in 1:burn_in) {
    theta0 = sampling.theta(x,sigmasq0)
    sigmasq0 = sampling.sigmasq(theta0,sigmasq0)
    if(k%%1000==0){print(k)}
  }
  cat("burn-in completed","\n")
  
  ## sampling procedure
  for (k in 1:nmc) {
    theta0 = sampling.theta(x,sigmasq0)
    theta.nmc[,k] = theta0
    sigmasq0 = sampling.sigmasq(theta0,sigmasq0)
    sigmasq.nmc[,k] = sigmasq0
    if(k%%1000==0){print(k)}
  }
  
  samples = list(theta = theta.nmc,sigmasq = sigmasq.nmc)
  return(samples)
}

################################### Example
set.seed(2023)
n = 1000
p = 0.2
sigma = 1
simdata = DGP(n,p,sigma)

theta0 = rep(1,n)
sigmasq0 = rep(1,n)
burn_in = 10000
nmc = 10000

ATHMC.sample = ATHMC(simdata$x,burn_in,nmc,theta0,sigmasq0)
