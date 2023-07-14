## X \sim Gamma(1/2,1)
plot(seq(0.001,5,0.001),dgamma(seq(0.001,5,0.001),0.5,1),type = "l",
     xlab = "x", ylab = "density")

## Concentration probability P(X<\epsilon)
epsilon = c(seq(0.00001,0.001,0.00001),seq(0.001,0.05,0.0001))
ConP = pgamma(epsilon,0.5,1)
plot(ConP~epsilon,type = "l",xlab = "epsilon",ylab = "Concentration Probability")

## adaptive random-walk Metropolis-Hastings with sigmoid approximation
N = 100000
burn = 20000
thin = 8
max_batch = 100
batch_size = 200

arwMH.gamma <- function(x0,alpha,beta,max_batch,batch_size,N,burn,thin){
  x.mcmc = rep(0,(N-burn)/thin)
  
  ## tuning procedure
  ss = 0.5
  for (j in 1:max_batch) {
    ar = 0
    for (i in 1:batch_size){
      x.star = x0 + rnorm(1,0,ss)
      Sigmoid = (1+exp(-1000*x0))/(1+exp(-1000*x.star))
      MH = min(1, Sigmoid*exp(alpha*(log(abs(x0))-log(abs(x.star))) + beta*(x0-x.star)))
      U = runif(1,0,1)
      if(U<=MH){x0 = x.star; ar = ar + 1}
    }
    ar = ar/batch_size
    if(ar<0.3){ss = ss*0.9}
    if(ar>0.5){ss = ss*1.1}
  }
  
  ## sampling procedure
  ar = 0
  for (i in 1:N) {
    x.star = x0 + rnorm(1,0,ss)
    Sigmoid = (1+exp(-1000*x0))/(1+exp(-1000*x.star))
    MH = min(1, Sigmoid*exp(alpha*(log(abs(x0))-log(abs(x.star))) + beta*(x0-x.star)))
    U = runif(1,0,1)
    if(U<=MH){x0 = x.star; ar = ar + 1}
    
    if(i>burn & i%%thin==0){
      x.mcmc[(i-burn)/thin] = x0
    }
  }
  ar = ar/N
  return(list(gamma = x.mcmc, accept = ar, step_size = ss))
}

set.seed(2023)
gamma.sample = arwMH.gamma(0.5,0.5,1,max_batch,batch_size,N,burn,thin)

plot(gamma.sample$gamma,type = "l") # traceplot
acf(gamma.sample$gamma) #acf
mcmcse::ess(gamma.sample$gamma) # effective sample size
plot(density(gamma.sample$gamma),xlim = c(0,5)) # posterior density

## Empirical Concentration probability P(X<\epsilon)
ConP.empirical = rep(0,length(epsilon))
for (i in 1:length(epsilon)) {
  ConP.empirical[i] = length(which(abs(gamma.sample$gamma)<epsilon[i]))/((N-burn)/thin)
}
plot(ConP.empirical~epsilon,type = "l",
     xlab = "epsilon",ylab = "Concentration Probability",col = "red",
     xlim = c(0,0.05),ylim = c(0,0.3))
par(new = T)
plot(ConP~epsilon,type = "l",xlab = "",ylab = "",xlim = c(0,0.05),ylim = c(0,0.3))

#### other applicable methods
#### method 1: transformation
#### If Y \sim beta(alpha,1-alpha) and Z \sim Exp(1), X = YZ \sim gamma(alpha,1)
set.seed(2023)
Y = rbeta((N-burn)/thin,1/2,1/2)
Z = rexp((N-burn)/thin,1)
X = Y*Z
plot(density(X),xlim = c(0,5)) # posterior density

ConP.empirical2 = rep(0,length(epsilon))
for (i in 1:length(epsilon)) {
  ConP.empirical2[i] = length(which(abs(X)<epsilon[i]))/((N-burn)/thin)
}
plot(ConP.empirical~epsilon,type = "l",
     xlab = "epsilon",ylab = "Concentration Probability",col = "red",
     xlim = c(0,0.05),ylim = c(0,0.3))
par(new = T)
plot(ConP.empirical2~epsilon,type = "l",
     xlab = "epsilon",ylab = "Concentration Probability",col = "blue",
     xlim = c(0,0.05),ylim = c(0,0.3))
par(new = T)
plot(ConP~epsilon,type = "l",xlab = "",ylab = "",xlim = c(0,0.05),ylim = c(0,0.3))

# beta(1/2,1/2) also troublesome
# U, V iid unif(0,1) and U^alpha/(U^alpha+V^\beta)|(U^alpha+V^\beta)<1 \sim beta(\alpha,\beta)
# Johnk (1964) algorithm (not good if alpha or beta large)

#### method 2: R
set.seed(2023)
X = rgamma((N-burn)/thin,0.5,1) 

#### How the estimation affected by algorithm
# mean = 0.5
bias.MH = mean(gamma.sample$gamma) - 0.5

# consider gamma(alpha,1), mean = alpha so sample mean is the estimator
alpha = seq(0.1,0.9,0.05)
alpha.esti = rep(0,length(alpha))
set.seed(2023)
for (i in 1:length(alpha)) {
  gamma.sample = arwMH.gamma(0.05,alpha[i],1,max_batch,batch_size,N,burn,thin)
  alpha.esti[i] = mean(gamma.sample$gamma)
}
plot(alpha.esti~alpha)



