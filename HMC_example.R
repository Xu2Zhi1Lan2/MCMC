## Assume we are sampling from a bivariate normal distribution with mean vector mu = c(0,0)
## and covariance matrix Sigma = matrix(c(1,0.25,0.25,1),2,2). We will try

## sampling with HMC, use H(x,p) = x'Sigma.inv x + p'p/2
## U(x) = - log bi-normal = x'Sigma.inv x, grad_U = 2Sigma.inv x
bn_HMC <- function(x0,burn_in,nmc){
  
  Sigma = matrix(c(1,0.25,0.25,1),2,2) # set the covariance matrix
  Sigma.inv = solve(Sigma)
  
  ## set outputs
  x.sample = matrix(0,2,nmc)
  
  ## burn_in procedure
  step0 = 1 # initial time step size for leapfrog
  
  for (i in 1:burn_in) {
    p0 = rnorm(2,0,1)
    p.star = p0 - step0*Sigma.inv%*%x0 # half-step jump in p-direction
    x.star = x0 + step0*p.star # one-step jump in x-direction
    p.star = p.star - step0*Sigma.inv%*%x.star # half-step jump in p-direction
    p.star = -p.star # negative momentum to make trajectory reversible (symmetric)
    
    H.star = t(x.star)%*%Sigma.inv%*%x.star + 0.5*sum(p.star^2) # proposal Hamiltonian
    H0 = t(x0)%*%Sigma.inv%*%x0 + 0.5*sum(p0^2) # current Hamiltonian
    
    alpha = min(1,exp(H0-H.star))
    U = runif(1,0,1)
    if(U<alpha){x0 = x.star}
    if(i%%1000==0){print(i)}
  }
  
  cat("Burn-in completed","\n")
  
  ## sampling procedure
  for (i in 1:nmc) {
    p0 = rnorm(2,0,1)
    p.star = p0 - step0*Sigma.inv%*%x0 # half-step jump in p-direction
    x.star = x0 + step0*p.star # one-step jump in x-direction
    p.star = p.star - step0*Sigma.inv%*%x.star # half-step jump in p-direction
    p.star = -p.star # negative momentum to make trajectory reversible (symmetric)
    
    H.star = t(x.star)%*%Sigma.inv%*%x.star + 0.5*sum(p.star^2) # proposal Hamiltonian
    H0 = t(x0)%*%Sigma.inv%*%x0 + 0.5*sum(p0^2) # current Hamiltonian
    
    alpha = min(1,exp(H0-H.star))
    U = runif(1,0,1)
    if(U<alpha){x0 = x.star;p0 = p.star}
    x.sample[,i]=x0
    
    if(i%%1000==0){print(i)}
  }
  
  return(x.sample)
}

### Simple test
set.seed(2023)
x0 = c(0,0)
burn_in = 10000
nmc = 10000
sample.HMC = bn_HMC(x0,burn_in,nmc)

## HMC
plot(sample.HMC[1,1:nmc],type = "l")
1 - mean(sample.HMC[1,1:(nmc-1)] == sample.HMC[1,2:nmc]) 
acf(sample.HMC[1,1:nmc])
plot(density(sample.HMC[1,1:nmc]),xlim = c(-4,4))
ess(sample.HMC[1,1:nmc])

plot(sample.HMC[2,1:nmc],type = "l")
1 - mean(sample.HMC[2,1:(nmc-1)] == sample.HMC[2,2:nmc]) 
acf(sample.HMC[2,1:nmc])
plot(density(sample.HMC[2,1:nmc]),xlim = c(-4,4))
ess(sample.HMC[2,1:nmc])

cor(sample.HMC[1,1:nmc],sample.HMC[2,1:nmc])

### illustration of trajectory 
L = 100
step0 = 1
Sigma = matrix(c(1,-0.4,-0.4,1),2,2)
Sigma.inv = solve(Sigma)
x.traj = matrix(2,2,L)
p0 = rnorm(2,0,1)
p.star = p0 - step0*Sigma.inv%*%x.traj[,1] # half-step jump in p-direction
for (i in 1:(L-1)) {
  x.traj[,(i+1)] = x.traj[,i] + step0*p.star
  if(i<L){
    p.star = p.star - 2*step0*Sigma.inv%*%x.traj[,(i+1)] # one-step jump in p-direction
  }
  else{
    p.star = p.star - step0*Sigma.inv%*%x.traj[,(i+1)] # half-step jump in p-direction
  }
}
p.star = -p.star # negative momentum to make trajectory reversible (symmetric)

Tj = c(seq(1,100,10),100)
X1 = seq(-4,4,0.01)
X2 = seq(-4,4,0.01)
Z = matrix(0,length(X1),length(X2))
for (i in 1:length(X1)) {
  for (j in 1:length(X2)) {
    mu = c(X1[i],X2[j])
    Z[i,j] = t(mu)%*%Sigma.inv%*%mu
  }
}
contour(X1,X2,Z)
lines(x = x.traj[1,Tj],y = x.traj[2,Tj],col = "red")
points(x = x.traj[1,Tj],y = x.traj[2,Tj],col = "red")



































