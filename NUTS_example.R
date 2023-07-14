source("C:/MyPhD/Year 4/Rref/HMC/NUTS.R")

## consider standard bivariate normal
f <- function(x){
  return(-0.5*sum(t(x)*x))
}

grad_f <- function(x){
  return(-x)
}

theta.nuts <- NUTS(c(1,1), f, grad_f, n_iter = 10000, M_diag = NULL, M_adapt = 5000, delta = 0.5,
     max_treedepth = 50, eps = 1, verbose = TRUE)

plot(theta.nuts[,1],type = "l");acf(theta.nuts[,1])
plot(theta.nuts[,2],type = "l");acf(theta.nuts[,2])

