Horselike_norm.mean_rwMH <- function(x,n,theta0,Burn_in,nmc,stepsize=1){
  ### random-walk M-H for normal mean model with Horseshoe-like prior
  ### Input
  ### x: data
  ### n: length of data
  ### theta0: initial value
  ### Burn_in: length of burn-in
  ### nmc: length of Markov chain
  
  ### set return
  theta.sample = matrix(0,n,nmc)
  
  ### Burn-in
  for (i in 1:Burn_in) {
    for (j in 1:n) {
      theta.star = theta0[j] + rnorm(1,0,stepsize)
      A = ((x[j]-theta0[j])^2 - (x[j]-theta.star)^2)/2
      A = A + log(log(1+1/theta.star^2)) - log(log(1+1/(theta0[j])^2))
      A = min(1,exp(A))
      U = runif(1,0,1)
      if(U<=A){theta0[j]=theta.star}
    }
  }
  cat("Burn-in completed","\n")
  
  ### sample from stationary distribution 
  for (i in 1:nmc) {
    for (j in 1:n) {
      theta.star = theta0[j] + rnorm(1,0,stepsize)
      A = ((x[j]-theta0[j])^2 - (x[j]-theta.star)^2)/2
      A = A + log(log(1+1/theta.star^2)) - log(log(1+1/(theta0[j])^2))
      A = min(1,exp(A))
      U = runif(1,0,1)
      if(U<=A){theta0[j]=theta.star;theta.sample[j,i]=theta.star}
      else{theta.sample[j,i]=theta0[j]}
    }
    if(i%%1000==0){print(i)}
  }
  
  ### return
  return(theta.sample)
}

##################################################################

### data generating
n = 200
theta.true = rbinom(n,1,0.2)*rnorm(n,0,sqrt(2*log(n)))
x = theta.true + rnorm(n,0,1)

### set initial value
theta0 = rep(1,n)

### set sample size
Burn_in = 10000
nmc = 10000

### run rwMH
theta.sample = Horselike_norm.mean_rwMH(x,n,theta0,Burn_in,nmc)

#####################################################################

### traceplot and acf of theta200
plot(theta.sample[200,],type = "l")
plot(density(theta.sample[200,]))
acf(theta.sample[200,])

### traceplot and acf of theta198
plot(theta.sample[198,],type = "l")
plot(density(theta.sample[198,]))
acf(theta.sample[198,])

### acceptance rate
acc.rate = rep(0,n)
for (j in 1:200) {
  acc.rate[j] = 1-mean(theta.sample[j,2:Burn_in]==theta.sample[j,1:(Burn_in-1)])
}
summary(acc.rate)

### theta estimation
theta.mean = rowMeans(theta.sample)
plot(theta.mean~x);abline(a=0,b=1)

#########################################################################
## adptive MCMC
#########################################################################

Horselike_norm.mean_rwaMH <- function(x,n,theta0,Burn_in,nmc,Tuning=100,stepsize=1){
  ### random-walk M-H for normal mean model with Horseshoe-like prior
  ### Input
  ### x: data
  ### n: length of data
  ### theta0: initial value
  ### Burn_in: length of burn-in
  ### nmc: length of Markov chain
  
  ### set return
  theta.burnin = matrix(0,n,Burn_in)
  theta.sample = matrix(0,n,nmc)
  stepsize = rep(stepsize,n)
  
  ### Burn-in & tuning stepsize
  for (i in 1:Burn_in) {
    
    ## M-H step
    for (j in 1:n) {
      theta.star = theta0[j] + rnorm(1,0,stepsize[j])
      A = ((x[j]-theta0[j])^2 - (x[j]-theta.star)^2)/2
      A = A + log(log(1+1/theta.star^2)) - log(log(1+1/(theta0[j])^2))
      A = min(1,exp(A))
      U = runif(1,0,1)
      if(U<=A){theta0[j]=theta.star;theta.burnin[j,i]=theta.star}
      else{theta.burnin[j,i]=theta0[j]}
    }
    
    ## tuning stepsize
    if(i%%Tuning==0){
      I0 = i/Tuning
      Index.range0 = (Tuning*(I0-1)+1):(Tuning*I0-1)
      Index.range1 = (Tuning*(I0-1)+2):(Tuning*I0)
      for (j in 1:n) {
        A.rate = 1-mean(theta.burnin[j,Index.range0] == theta.burnin[j,Index.range1])
        if(A.rate>0.5){stepsize[j]=exp(log(stepsize[j])+0.01)}
        if(A.rate<0.3){stepsize[j]=exp(log(stepsize[j])-0.01)}
      }
    }
  }
  cat("Burn-in completed","\n")
  
  ### sample from stationary distribution 
  for (i in 1:nmc) {
    for (j in 1:n) {
      theta.star = theta0[j] + rnorm(1,0,stepsize[j])
      A = ((x[j]-theta0[j])^2 - (x[j]-theta.star)^2)/2
      A = A + log(log(1+1/theta.star^2)) - log(log(1+1/(theta0[j])^2))
      A = min(1,exp(A))
      U = runif(1,0,1)
      if(U<=A){theta0[j]=theta.star;theta.sample[j,i]=theta.star}
      else{theta.sample[j,i]=theta0[j]}
    }
    if(i%%1000==0){print(i)}
  }
  
  ### return
  return(theta.sample)
}

####################################################################

theta.sample.aMH = Horselike_norm.mean_rwaMH(x,n,theta0,Burn_in,nmc)
  
### acceptance rate
acc.rate = rep(0,n)
for (j in 1:200) {
  acc.rate[j] = 1-mean(theta.sample.aMH[j,2:Burn_in]==theta.sample.aMH[j,1:(Burn_in-1)])
}
summary(acc.rate)

### traceplot and acf of theta200
plot(theta.sample.aMH[200,],type = "l")
plot(density(theta.sample.aMH[200,]))
acf(theta.sample.aMH[200,])

### traceplot and acf of theta198
plot(theta.sample.aMH[198,],type = "l")
plot(density(theta.sample.aMH[198,]))
acf(theta.sample.aMH[198,])

### theta estimation
theta.mean.aMH = rowMeans(theta.sample.aMH)
plot(theta.mean.aMH~x);abline(a=0,b=1)
