### Example 1: value the hyperparameter in Horseshoe model
library(horseshoe)

## Data Generation
set.seed(2023)
m = 200 # dimension
p = 0.2 # degree of sparsity
psi = sqrt(2*log(m)) # signal strength
sigma0 = 1 # noise strength

Z = rbinom(m,0,p)
theta = Z*rnorm(m,0,psi)
x = theta + rnorm(m,0,sigma)

## candidate hyperparameter
tau = seq(0.01,1,0.05); L = length(tau)

## function to evaluate DIC
loglikelihood <- function(theta,x,sigma=sigma0){sum(dnorm(x, theta, sigma, log=T))}

DIC.HS <- function(x,tau){
  HSmodel <- HS.normal.means(x,
                             method.tau = "fixed",
                             tau = tau,
                             method.sigma = "fixed",
                             burn = 1000,
                             nmc = 5000)
  
  theta.Bayes = HSmodel$BetaHat
  pD = 2*loglikelihood(theta.Bayes,x)
  for (i in 1:5000) {
    pD = pD - 2*loglikelihood(HSmodel$BetaSamples[i,],x)/5000
  }
  DIC = -2*loglikelihood(theta.Bayes,x) + 2*pD
  return(DIC)
}

DIC.tau = rep(0,L)
for (k in 1:L) {
  DIC.tau[k] = DIC.HS(x,tau[k])
}

#########################################################################
### Example 2: choose a MCMC algorithm
plot.post <- function(x,alpha=0.5,beta=2){
  theta = c(seq(-10,-0.01,0.01),seq(0.01,10,0.01))
  logP = (alpha-1)*log(abs(theta)) -beta*abs(theta) - 0.5*theta^2 + x*theta
  Post.kernel = exp(logP)
  plot(Post.kernel~theta,type = "l")
}

### Method 1: adaptive random-walk Metropolis-Hastings
arwMH <- function(theta0,x,burn_in,nmc,alpha=0.5,beta=5,Tuning=100){
  ## set output
  theta.nmc = rep(0,nmc)
  
  ## burn-in procedure
  step.theta = 1 # set initial step size
  theta.burn = rep(0,burn_in) # record burn-in samples
  step.theta.path = rep(0,burn_in/Tuning) # record the tuning path of step size
  accept.theta.path = rep(0,burn_in/Tuning) # record the tuning path of acceptance rate
  
  for (i in 1:burn_in) {
    ## Metropolis-Hastings
    theta.star = theta0 + rnorm(1,0,step.theta) # proposal
    MH = min(exp((alpha-1)*(log(abs(theta.star))-log(abs(theta0)))-
                   beta*(abs(theta.star)-abs(theta0))+x*(theta.star-theta0)+
                   0.5*(theta0^2-theta.star^2)),1) # acceptance probability
    U = runif(1,0,1)
    if(U<=MH){theta0 = theta.star}
    theta.burn[i] = theta0
    
    ## adaptive tuning
    if(i%%Tuning==0){
      # get the tuning window
      I0 = i/Tuning
      Index.range0 = (Tuning*(I0-1)+1):(Tuning*I0-1)
      Index.range1 = (Tuning*(I0-1)+2):(Tuning*I0)
      # compute the acceptance rate
      A.rate = 1-mean(theta.burn[Index.range0] == theta.burn[Index.range1])
      # targeting the acceptance rate between 0.3-0.5
      if(A.rate>0.5){step.theta=exp(log(step.theta)+0.1)}
      if(A.rate<0.3){step.theta=exp(log(step.theta)-0.1)}
      # record the path of tuning
      step.theta.path[I0] = step.theta
      accept.theta.path[I0] = A.rate
    }
  }
  
  ## Sampling procedure
  for (i in 1:burn_in) {
    ## Metropolis-Hastings
    theta.star = theta0 + rnorm(1,0,step.theta) # proposal
    MH = min(exp((alpha-1)*(log(abs(theta.star))-log(abs(theta0)))-
                   beta*(abs(theta.star)-abs(theta0))+x*(theta.star-theta0)+
                   0.5*(theta0^2-theta.star^2)),1) # acceptance probability
    U = runif(1,0,1)
    if(U<=MH){theta0 = theta.star}
    theta.nmc[i] = theta0
  }
  
  ## return
  results = list(theta = theta.nmc,theta.burn = theta.burn,
                 step.theta.path = step.theta.path,
                 accept.theta.path = accept.theta.path)
  return(results)
}

### Method 2: Hamiltonian Monte Carlo
U <- function(theta,x,alpha=0.5,beta=2){
  return(-(alpha-1)*log(abs(theta)) + beta*abs(theta) + 0.5*theta^2 - x*theta)
}

grad_U <- function(theta,x,alpha=0.5,beta=2){
  return(-(alpha-1)/theta + beta*(2*(theta>0)-1) + theta - x)
}

HMC <- function(theta0,x,burn_in,nmc,alpha=0.5,beta=5,Tuning=100,L=50){
  ## set output
  theta.nmc = rep(0,nmc)
  
  ## burn-in procedure
  step.theta = 1 # set initial leapfrog step size
  theta.burn = rep(0,burn_in) # record burn-in samples
  step.theta.path = rep(0,burn_in/Tuning) # record the tuning path of step size
  accept.theta.path = rep(0,burn_in/Tuning) # record the tuning path of acceptance rate
  
  for (i in 1:burn_in) {
    theta.traj = rep(0,L);theta.traj[1] = theta0
    p0 = rnorm(1,0,1)
    p.star = p0 - 0.5*step.theta*grad_U(theta.traj[1],x)*theta.traj[1] # half-step jump in p-direction
    for (k in 1:(L-1)) {
      theta.traj[k+1] = theta.traj[k] + step.theta*p.star
      if(k<L){
        p.star = p.star - step.theta*grad_U(theta.traj[k+1],x)*theta.traj[k+1] # one-step jump in p-direction
      }
      else{
        p.star = p.star - 0.5*step0*grad_U(theta.traj[k+1],x)*theta.traj[k+1] # half-step jump in p-direction
      }
    }
    p.star = -p.star # negative momentum to make trajectory reversible (symmetric)
    theta.star = theta.traj[L]
    
    H.star = U(theta.star,x) + 0.5*sum(p.star^2) # proposal Hamiltonian
    H0 = U(theta0,x) + 0.5*sum(p0^2) # current Hamiltonian
    
    alpha = min(1,exp(H0-H.star))
    U = runif(1,0,1)
    if(U<alpha){theta0 = theta.star}
    theta.burn[i] = theta0
    
    ## adaptive tuning
    if(i%%Tuning==0){
      # get the tuning window
      I0 = i/Tuning
      Index.range0 = (Tuning*(I0-1)+1):(Tuning*I0-1)
      Index.range1 = (Tuning*(I0-1)+2):(Tuning*I0)
      # compute the acceptance rate
      A.rate = 1-mean(theta.burn[Index.range0] == theta.burn[Index.range1])
      # targeting the acceptance rate between 0.3-0.7
      if(A.rate>0.7){step.theta=exp(log(step.theta)+0.1)}
      if(A.rate<0.3){step.theta=exp(log(step.theta)-0.1)}
      # record the path of tuning
      step.theta.path[I0] = step.theta
      accept.theta.path[I0] = A.rate
    }
  }
  
  ## Sampling procedure
  for (i in 1:burn_in) {
    theta.traj = rep(0,L);theta.traj[1] = theta0
    p0 = rnorm(1,0,1)
    p.star = p0 - 0.5*step.theta*grad_U(theta.traj[1],x)*theta.traj[1] # half-step jump in p-direction
    for (k in 1:(L-1)) {
      theta.traj[k+1] = theta.traj[k] + step.theta*p.star
      if(k<L){
        p.star = p.star - step.theta*grad_U(theta.traj[k+1],x)*theta.traj[k+1] # one-step jump in p-direction
      }
      else{
        p.star = p.star - 0.5*step0*grad_U(theta.traj[k+1],x)*theta.traj[k+1] # half-step jump in p-direction
      }
    }
    p.star = -p.star # negative momentum to make trajectory reversible (symmetric)
    theta.star = theta.traj[L]
    
    H.star = U(theta.star,x) + 0.5*sum(p.star^2) # proposal Hamiltonian
    H0 = U(theta0,x) + 0.5*sum(p0^2) # current Hamiltonian
    
    alpha = min(1,exp(H0-H.star))
    U = runif(1,0,1)
    if(U<alpha){theta0 = theta.star}
    theta.nmc[i] = theta0
  }
  
  ## return
  results = list(theta = theta.nmc,theta.burn = theta.burn,
                 step.theta.path = step.theta.path,
                 accept.theta.path = accept.theta.path)
  return(results)
}

### Method 3: adaptive Langevin Metropolis-Hastings
aLMH <- function(theta0,x,burn_in,nmc,alpha=0.5,beta=5,Tuning=100){
  ## set output
  theta.nmc = rep(0,nmc)
  
  ## burn-in procedure
  step.theta = 1 # set initial step size
  theta.burn = rep(0,burn_in) # record burn-in samples
  step.theta.path = rep(0,burn_in/Tuning) # record the tuning path of step size
  accept.theta.path = rep(0,burn_in/Tuning) # record the tuning path of acceptance rate
  
  for (i in 1:burn_in) {
    ## Langevin Metropolis-Hastings
    theta.star = theta0 - 0.5*step.theta^2*grad_U(theta0,x) + step.theta*rnorm(1,0,1) # proposal
    post.ratio = (alpha-1)*(log(abs(theta.star))-log(abs(theta0)))-
      beta*(abs(theta.star)-abs(theta0))+x*(theta.star-theta0)+
      0.5*(theta0^2-theta.star^2)
    trans.ratio = (0.5/step.theta^2)*(theta.star-theta0-0.5*step.theta^2*grad_U(theta0,x))^2-
      (0.5/step.theta^2)*(theta0-theta.star-0.5*step.theta^2*grad_U(theta.star,x))^2
    MH = min(exp(trans.ratio + post.ratio),1) # acceptance probability
    U = runif(1,0,1)
    if(U<=MH){theta0 = theta.star}
    theta.burn[i] = theta0
    
    ## adaptive tuning
    if(i%%Tuning==0){
      # get the tuning window
      I0 = i/Tuning
      Index.range0 = (Tuning*(I0-1)+1):(Tuning*I0-1)
      Index.range1 = (Tuning*(I0-1)+2):(Tuning*I0)
      # compute the acceptance rate
      A.rate = 1-mean(theta.burn[Index.range0] == theta.burn[Index.range1])
      # targeting the acceptance rate between 0.3-0.7
      if(A.rate>0.6){step.theta=exp(log(step.theta)+0.1)}
      if(A.rate<0.4){step.theta=exp(log(step.theta)-0.1)}
      # record the path of tuning
      step.theta.path[I0] = step.theta
      accept.theta.path[I0] = A.rate
    }
  }
  
  ## Sampling procedure
  for (i in 1:nmc) {
    ## Langevin Metropolis-Hastings
    theta.star = theta0 - 0.5*step.theta^2*grad_U(theta0,x) + step.theta*rnorm(1,0,1) # proposal
    post.ratio = (alpha-1)*(log(abs(theta.star))-log(abs(theta0)))-
      beta*(abs(theta.star)-abs(theta0))+x*(theta.star-theta0)+
      0.5*(theta0^2-theta.star^2)
    trans.ratio = (0.5/step.theta^2)*(theta.star-theta0-0.5*step.theta^2*grad_U(theta0,x))^2-
      (0.5/step.theta^2)*(theta0-theta.star-0.5*step.theta^2*grad_U(theta.star,x))^2
    MH = min(exp(trans.ratio + post.ratio),1) # acceptance probability
    U = runif(1,0,1)
    if(U<=MH){theta0 = theta.star}
    theta.nmc[i] = theta0
  }
  
  ## return
  results = list(theta = theta.nmc,theta.burn = theta.burn,
                 step.theta.path = step.theta.path,
                 accept.theta.path = accept.theta.path)
  return(results)
}

### use x = 4 as example
x3 = 4
theta0 = 1
burn_in = 20000
nmc = 10000
plot.post(x3)

set.seed(2023)
theta.arwMH <- arwMH(theta0,x3,burn_in,nmc)
theta.HMC <- HMC(theta0,x3,burn_in,nmc)
theta.aLMH <- aLMH(theta0,x3,burn_in,nmc)

DIC.MCMC <- function(theta,theta.Bayes,x = x3, n = nmc){
  pD = 2*loglikelihood(theta.Bayes,x)
  for (i in 1:n) {
    pD = pD - 2*loglikelihood(theta[i],x)/n
  }
  DIC = -2*loglikelihood(theta.Bayes,x) + 2*pD
  return(DIC)
}

### arwMH
theta1 = theta.arwMH$theta
theta1.Bayes = mean(theta1)
(DIC1 = DIC.MCMC(theta1,theta1.Bayes))

### HMC
theta2 = theta.HMC$theta
theta2.Bayes = mean(theta2)
(DIC2 = DIC.MCMC(theta2,theta2.Bayes))

### aLMH
theta3 = theta.aLMH$theta
theta3.Bayes = mean(theta3)
(DIC3 = DIC.MCMC(theta3,theta3.Bayes))
