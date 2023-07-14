## To get this script
## https://pan.baidu.com/s/1fRlCC2OtBcoi9w7ick2TGA
## qg9k

## Model: x \sim N(theta,1)
## Prior: theta \sim Gamma(alpha,beta), alpha = 0.5,beta = 2

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

### example 1: x = 0.2 (small observation, prior dominates)
x1 = 0.2
theta0 = 1
burn_in = 20000
nmc = 10000
plot.post(x1)

# adaptive random-walk Metropolis-Hastings
set.seed(2023)
theta.arwMH <- arwMH(theta0,x1,burn_in,nmc)
plot(theta.arwMH$theta,type = "l") # traceplot
plot(density(theta.arwMH$theta),xlim = c(-10,10)) # posterior density
acf(theta.arwMH$theta)
cat("acceptance rate is:",1-mean(theta.arwMH$theta[1:(nmc-1)]==theta.arwMH$theta[2:nmc]),"\n")

# Hamiltonian Monte Carlo
set.seed(2023)
theta.HMC <- HMC(theta0,x1,burn_in,nmc)
plot(theta.HMC$theta,type = "l") # traceplot
plot(density(theta.HMC$theta),xlim = c(-10,10)) # posterior density
acf(theta.HMC$theta)
cat("acceptance rate is:",1-mean(theta.HMC$theta[1:(nmc-1)]==theta.HMC$theta[2:nmc]),"\n")

# adaptive Langevin Metropolis-Hastings
set.seed(2023)
theta.aLMH <- aLMH(theta0,x1,burn_in,nmc)
plot(theta.aLMH$theta,type = "l") # traceplot
plot(density(theta.aLMH$theta),xlim = c(-10,10)) # posterior density
acf(theta.aLMH$theta)
cat("acceptance rate is:",1-mean(theta.aLMH$theta[1:(nmc-1)]==theta.aLMH$theta[2:nmc]),"\n")


### example 2: x = 2.4 (most difficult problem)
x2 = 2.4
theta0 = 1
burn_in = 20000
nmc = 10000
plot.post(x2)

# adaptive random-walk Metropolis-Hastings
set.seed(2023)
theta.arwMH <- arwMH(theta0,x2,burn_in,nmc)
plot(theta.arwMH$theta,type = "l") # traceplot
plot(density(theta.arwMH$theta),xlim = c(-10,10)) # posterior density
acf(theta.arwMH$theta)
cat("acceptance rate is:",1-mean(theta.arwMH$theta[1:(nmc-1)]==theta.arwMH$theta[2:nmc]),"\n")

# Hamiltonian Monte Carlo
set.seed(2023)
theta.HMC <- HMC(theta0,x2,burn_in,nmc)
plot(theta.HMC$theta,type = "l") # traceplot
plot(density(theta.HMC$theta),xlim = c(-10,10)) # posterior density
acf(theta.HMC$theta)
cat("acceptance rate is:",1-mean(theta.HMC$theta[1:(nmc-1)]==theta.HMC$theta[2:nmc]),"\n")

# adaptive Langevin Metropolis-Hastings
set.seed(2023)
theta.aLMH <- aLMH(theta0,x2,burn_in,nmc)
plot(theta.aLMH$theta,type = "l") # traceplot
plot(density(theta.aLMH$theta),xlim = c(-10,10)) # posterior density
acf(theta.aLMH$theta)
cat("acceptance rate is:",1-mean(theta.aLMH$theta[1:(nmc-1)]==theta.aLMH$theta[2:nmc]),"\n")


### example 3: x = 4 (large observation, bi-modal)
x3 = 4
theta0 = 1
burn_in = 20000
nmc = 10000
plot.post(x3)

# adaptive random-walk Metropolis-Hastings
set.seed(2023)
theta.arwMH <- arwMH(theta0,x3,burn_in,nmc)
plot(theta.arwMH$theta,type = "l") # traceplot
plot(density(theta.arwMH$theta),xlim = c(-10,10)) # posterior density
acf(theta.arwMH$theta)
cat("acceptance rate is:",1-mean(theta.arwMH$theta[1:(nmc-1)]==theta.arwMH$theta[2:nmc]),"\n")

# Hamiltonian Monte Carlo
set.seed(2023)
theta.HMC <- HMC(theta0,x3,burn_in,nmc)
plot(theta.HMC$theta,type = "l") # traceplot
plot(density(theta.HMC$theta),xlim = c(-10,10)) # posterior density
acf(theta.HMC$theta)
cat("acceptance rate is:",1-mean(theta.HMC$theta[1:(nmc-1)]==theta.HMC$theta[2:nmc]),"\n")

# adaptive Langevin Metropolis-Hastings
set.seed(2023)
theta.aLMH <- aLMH(theta0,x3,burn_in,nmc)
plot(theta.aLMH$theta,type = "l") # traceplot
plot(density(theta.aLMH$theta),xlim = c(-10,10)) # posterior density
acf(theta.aLMH$theta)
cat("acceptance rate is:",1-mean(theta.aLMH$theta[1:(nmc-1)]==theta.aLMH$theta[2:nmc]),"\n")


### example 4: x = 6 (large observation, likelihood dominates)
x4 = 6
theta0 = 1
burn_in = 20000
nmc = 10000
plot.post(x4)

# adaptive random-walk Metropolis-Hastings
set.seed(2023)
theta.arwMH <- arwMH(theta0,x4,burn_in,nmc)
plot(theta.arwMH$theta,type = "l") # traceplot
plot(density(theta.arwMH$theta),xlim = c(-10,10)) # posterior density
acf(theta.arwMH$theta)
cat("acceptance rate is:",1-mean(theta.arwMH$theta[1:(nmc-1)]==theta.arwMH$theta[2:nmc]),"\n")

# Hamiltonian Monte Carlo
set.seed(2023)
theta.HMC <- HMC(theta0,x4,burn_in,nmc)
plot(theta.HMC$theta,type = "l") # traceplot
plot(density(theta.HMC$theta),xlim = c(-10,10)) # posterior density
acf(theta.HMC$theta)
cat("acceptance rate is:",1-mean(theta.HMC$theta[1:(nmc-1)]==theta.HMC$theta[2:nmc]),"\n")

# adaptive Langevin Metropolis-Hastings
set.seed(2023)
theta.aLMH <- aLMH(theta0,x4,burn_in,nmc)
plot(theta.aLMH$theta,type = "l") # traceplot
plot(density(theta.aLMH$theta),xlim = c(-10,10)) # posterior density
acf(theta.aLMH$theta)
cat("acceptance rate is:",1-mean(theta.aLMH$theta[1:(nmc-1)]==theta.aLMH$theta[2:nmc]),"\n")





