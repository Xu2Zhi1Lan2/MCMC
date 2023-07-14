library(horseshoe)

set.seed(2023)
truetheta = rbinom(200,1,0.2)*rnorm(200,0,3.26)
x = truetheta + rnorm(200,0,1)

start.time = Sys.time()
HS = HS.normal.means(x,
                     method.tau = "halfCauchy",
                     method.sigma = "fixed",
                     burn = 5000,
                     nmc = 10000)
end.time = Sys.time()
time.used = end.time - start.time

theta.sample = t(HS$BetaSamples)
tau.sample = HS$TauSamples

## Estimation
theta.mean = rowMeans(theta.sample)
plot(theta.mean~x, xlim = c(-10,10), ylim = c(-10,10),
     xlab = "Observation",
     ylab = "Posterior Mean")
abline(a = 0, b = 1);abline(a = 0, b = 1/2, col = "red")

## traceplot
plot(tau.sample, type = "l")

## ACF (autocorrelation function)
ACF = acf(tau.sample, lag.max = 100)

## Effective Sample Size
ESS = 10000/(1+2*sum(ACF$acf[ACF$acf>=0.1]))
ESS/2