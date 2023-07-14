set.seed(2022)

# create grid
gd = seq(0, 1, .001)

# define colors for plotting
blue.3 = scales::alpha('blue', .3)
red.3 = scales::alpha('red', .3)
purple.5 = scales::alpha('purple', .5)

########################################################
## visualize slice sampler
par(mfrow = c(1, 2), mar = rep(1, 4))

### Unimodal Beta(2,5) density ###

plot(gd, dbeta(gd, 2, 5), type = 'l', ylim = c(0, 3),
     axes = F, frame.plot = F, col = 'grey50')
abline(h = 0)

# starting value
x0 = 0.15
abline(v = x0, col = blue.3)

# sample y uniformly from [0, f(x0)]
y1 = runif(1, 0, dbeta(x0, 2, 5))
abline(h = y1, col = red.3)

# find end-point of slice
m.roots = sort(
  rootSolve::uniroot.all(
    function(x) { dbeta(x, 2, 5) - y1 }, 
    lower = 0, 
    upper = 1
  )
)

# plot slice
segments(m.roots[1], y1, m.roots[2], y1, col = purple.5, lwd = 5)

### Bimodal mixture of Beta densities ###

# define density function
mix.beta = function(x) {
  .4 * dbeta(x, 2, 10) + .4 * dbeta(x, 10, 2) + .2 * dbeta(x, 5, 5)
}

# plot density
plot(gd, mix.beta(gd), type = 'l', ylim = c(0, 2), col = 'grey50',
     axes = F, frame.plot = F)
abline(h = 0)

# starting value
x0 = 0.15
abline(v = x0, col = blue.3)

# sample y uniformly from [0, f(x0)]
y2 = runif(1, 0, mix.beta(x0))
abline(h = y2, col = red.3)

# find end-points of slices 
m.roots = sort(
  rootSolve::uniroot.all(
    function(x) { mix.beta(x) - y2 }, 
    lower = 0, 
    upper = 1
  )
)

# draw slice
for (rr in 1:(length(m.roots) / 2)) {
  segments(m.roots[2 * rr - 1], y2, m.roots[2 * rr], y2, 
           col = purple.5, lwd = 5)
}

############################################################################
##### diagnostics of slice sampler

SliceSampler = function(x, pstar, w, n.sims) {
  ## x is the initial value
  ## pstar is a function proportional to target distribution
  ## w is the window size for stepping out
  ## n.sims is the number of desired random samples
  
  # define function to check whether proposal is in slice    
  in.slice = function(x, y) { y <= pstar(x) }
  
  # container for samples
  sims = matrix(NA, n.sims, 2)
  
  for (s in 1:n.sims) {
    ### sample y ### 
    y = runif(1, 0, pstar(x)) 
    
    ### sample x ###
    
    # initial window
    l = x - w * runif(1)
    u = l + w
    
    # expand lower-limit if necessary
    l.in = in.slice(l, y)
    if (l.in) {
      while (l.in) {
        l = l - w
        # check whether lower bound is hit
        if (l < 0) { 
          l = 0
          break
        }
        l.in = in.slice(l, y)
      }
    }
    
    # expand upper-limit if necessary
    u.in = in.slice(u, y)
    if (u.in) {
      while (u.in) {
        u = u + w
        # check whether upper bound is hit
        if (u > 1) { 
          u = 1
          break
        }
        u.in = in.slice(u, y)
      }
    }
    
    # sample x from y-slice and shrink
    x.old = x 
    x.in = FALSE
    while (!x.in) {
      # sample x    
      x = runif(1, l, u)
      # check whether x is in slice
      x.in = in.slice(x, y)
      # shrink interval    
      if (x > x.old) {
        u = x
      } else {
        l = x
      }
    }
    
    # store samples   
    sims[s, 1] = x
    sims[s, 2] = y
  }
  return(sims)
}

##### set parameters
# number of simulation draws to take
n.sims = 30000

# initial value of x
x.init = .15

# beta(2, 5) density multiplied by constant
beta.prop = function(x, k = 1) { k * dbeta(x, 2, 5) }

# new mixture distribution, with smaller density in the middle
mix.beta = function(x) {
  .45 * dbeta(x, 2, 10) + .45 * dbeta(x, 10, 2) + .1 * dbeta(x, 3, 3)
}

# window size
w = .2

######## Time cost
start.time = Sys.time()
slice.samps.beta = SliceSampler(
  x.init,
  beta.prop,
  w, 
  n.sims)
print(Sys.time() - start.time)


start.time = Sys.time()
slice.samps.mix = SliceSampler(
  x.init,
  mix.beta,
  w,
  n.sims)
print(Sys.time() - start.time)


######## path of sampler
gen.path = function(sims, n.draw) {
  if (nrow(sims) > n.draw) {
    z = sims[1:(n.draw + 1), ]
  } else {
    z = rbind(sims, sims[nrow(sims), ])
  }
  
  path.mat = lapply(1:n.draw, function(w) {
    z1 = z2 = z[w, ] 
    z2[2] = z[w + 1, 2]
    return(rbind(z1, z2))
  }) 
  
  do.call('rbind', path.mat)
}

# subsample first 100 draws and generate path
n.draw = 100
path.mat.beta = gen.path(slice.samps.beta, n.draw)
path.mat.mix = gen.path(slice.samps.mix, n.draw)

par(mfrow = c(2, 2), mar = c(1, 2, 1, 1))
plot(gd, beta.prop(gd), type = 'l',
     axes = F, frame.plot = F, col = 'grey50',
     ylab = '')
abline(h = 0)
lines(path.mat.beta, col = purple.5)
points(slice.samps.beta[1:n.draw, ], col = purple.5, 
       pch = 19, cex = .75)
mtext(paste0('First ', n.draw, ' samples'), outer = F, side = 2)

plot(gd, mix.beta(gd), type = 'l',
     axes = F, frame.plot = F, col = 'grey50')
abline(h = 0)
lines(path.mat.mix, type = 'l', col = purple.5)
points(slice.samps.mix[1:n.draw, ], col = purple.5, 
       pch = 19, cex = .75)

plot(gd, dbeta(gd, 2, 5), type = 'l',
     axes = F, frame.plot = F, col = 'grey50',
     ylab = '')
abline(h = 0)
hist(slice.samps.beta[ , 1], border = purple.5, 
     fill ='white', freq = F, add = T, breaks = 100)
mtext(paste0('All samples'), outer = F, side = 2)

plot(gd, mix.beta(gd), type = 'l',
     axes = F, frame.plot = F, col = 'grey50')
abline(h = 0)
hist(slice.samps.mix[ , 1], border = purple.5, 
     fill ='white', freq = F, add = T, breaks = 100)

####################### Autocorrelation and Traceplot
par(mfrow = c(3, 2), mar=c(2, 4, 2, 1))
blue.7 = scales::alpha('blue', .7)

coda::autocorr.plot(slice.samps.beta[ , 1], 
                    auto.layout = F, col = 'blue',
                    main = 'Beta(2,5) Density')
coda::autocorr.plot(slice.samps.mix[ , 1], 
                    auto.layout = F, col = 'blue',
                    main = 'Mixture of Beta Densities')
plot(slice.samps.beta[1:n.draw, 1], ylim = c(0, 1), 
     type = 'l', col = 'blue',
     ylab = 'First 100 Iterations')
plot(slice.samps.mix[1:n.draw, 1], ylim = c(0, 1), 
     type = 'l', col = 'blue',
     ylab = 'First 100 Iterations')
plot(slice.samps.beta[ , 1], 
     ylim = c(0, 1), 
     col = blue.7,
     pch = 16, cex = .1,
     ylab = 'All Iterations')
plot(slice.samps.mix[ , 1], ylim = c(0, 1), 
     col = blue.7, pch = 16, cex = .1,
     ylab = 'All Iterations')


###################### Effective Sample Size
esize = c(round(coda::effectiveSize(slice.samps.beta[ , 1]), 2),
          round(coda::effectiveSize(slice.samps.mix[ , 1]), 2))
names(esize) = c('Unimodal', 'Bimodal')
print(esize)

