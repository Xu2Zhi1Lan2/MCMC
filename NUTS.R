leapfrog_step = function(theta, r, eps, grad_f, M_diag){
  r_tilde <- r + 0.5 * eps * grad_f(theta)
  theta_tilde <- theta + eps * r_tilde / M_diag
  r_tilde <- r_tilde + 0.5 * eps * grad_f(theta_tilde)
  list(theta = theta_tilde, r = r_tilde)
}

joint_log_density = function(theta, r, f, M_diag){
  f(theta) - 0.5*sum(r**2 / M_diag)
}

find_reasonable_epsilon = function(theta, f, grad_f, M_diag, eps = 1, verbose = TRUE){
  r <- rnorm(length(theta), 0, sqrt(M_diag))
  proposed <- leapfrog_step(theta, r, eps, grad_f, M_diag)
  log_ratio <- joint_log_density(proposed$theta, proposed$r, f, M_diag) - joint_log_density(theta, r, f, M_diag)
  alpha <- ifelse(exp(log_ratio) > 0.5, 1, -1)
  if(is.nan(alpha)) alpha <- -1
  count <- 1
  while(is.nan(log_ratio) || alpha * log_ratio > (-alpha)*log(2)){
    eps <- 2**alpha * eps
    proposed <- leapfrog_step(theta, r, eps, grad_f, M_diag)
    log_ratio <- joint_log_density(proposed$theta, proposed$r, f, M_diag) - joint_log_density(theta, r, f, M_diag)
    count <- count + 1
    if(count > 100) {
      stop("Could not find reasonable epsilon in 100 iterations!")
    }
  }
  if(verbose) message("Reasonable epsilon = ", eps, " found after ", count, " steps")
  eps
}

check_NUTS = function(s, theta_plus, theta_minus, r_plus, r_minus){
  if(is.na(s)) return(0)
  condition1 <- crossprod(theta_plus - theta_minus, r_minus) >= 0
  condition2 <- crossprod(theta_plus - theta_minus, r_plus) >= 0
  s && condition1 && condition2
}

build_tree = function(theta, r, u, v, j, eps, theta0, r0, f, grad_f, M_diag, Delta_max = 1000){
  if(j == 0){
    proposed <- leapfrog_step(theta, r, v*eps, grad_f, M_diag)
    theta <- proposed$theta
    r <- proposed$r
    log_prob <- joint_log_density(theta, r, f, M_diag)
    log_prob0 <- joint_log_density(theta0, r0, f, M_diag)
    n <- (log(u) <= log_prob)
    s <- (log(u) < Delta_max + log_prob)
    alpha <- min(1, exp(log_prob - log_prob0))
    if(is.nan(alpha)) stop()
    if(is.na(s) || is.nan(s)){
      s <- 0
    }
    if(is.na(n) || is.nan(n)){
      n <- 0
    }
    return(list(theta_minus=theta, theta_plus=theta, theta=theta, r_minus=r,
                r_plus=r, s=s, n=n, alpha=alpha, n_alpha=1))
  } else{
    obj0 <- build_tree(theta, r, u, v, j-1, eps, theta0, r0, f, grad_f, M_diag)
    theta_minus <- obj0$theta_minus
    r_minus <- obj0$r_minus
    theta_plus <- obj0$theta_plus
    r_plus <- obj0$r_plus
    theta <- obj0$theta
    if(obj0$s == 1){
      if(v == -1){
        obj1 <- build_tree(obj0$theta_minus, obj0$r_minus, u, v, j-1, eps, theta0, r0, f, grad_f, M_diag)
        theta_minus <- obj1$theta_minus
        r_minus <- obj1$r_minus
      } else{
        obj1 <- build_tree(obj0$theta_plus, obj0$r_plus, u, v, j-1, eps, theta0, r0, f, grad_f, M_diag)
        theta_plus <- obj1$theta_plus
        r_plus <- obj1$r_plus
      }
      n <- obj0$n + obj1$n
      if(n != 0){
        prob <- obj1$n / n
        if(runif(1) < prob){
          theta <- obj1$theta
        }
      }
      s <- check_NUTS(obj1$s, theta_plus, theta_minus, r_plus, r_minus)
      alpha <- obj0$alpha + obj1$alpha
      n_alpha <- obj0$n_alpha + obj1$n_alpha
      
    } else{
      n <- obj0$n
      s <- obj0$s
      alpha <- obj0$alpha
      n_alpha <- obj0$n_alpha
    }
    if(is.na(s) || is.nan(s)){
      s <- 0
    }
    if(is.na(n) || is.nan(n)){
      n <- 0
    }
    return(list(theta_minus=theta_minus, theta_plus=theta_plus, theta=theta,
                r_minus=r_minus, r_plus=r_plus, s=s, n=n, alpha=alpha, n_alpha=n_alpha))
  }
}

NUTS <- function(theta, f, grad_f, n_iter, M_diag = NULL, M_adapt = 50, delta = 0.5, max_treedepth = 10, eps = 1, verbose = TRUE){
  theta_trace <- matrix(0, n_iter, length(theta))
  par_list <- list(M_adapt = M_adapt)
  for(iter in 1:n_iter){
    nuts <- NUTS_one_step(theta, iter, f, grad_f, par_list, delta = delta, max_treedepth = max_treedepth, eps = eps, verbose = verbose)
    theta <- nuts$theta
    par_list <- nuts$pars
    theta_trace[iter, ] <- theta
  }
  theta_trace
}


NUTS_one_step <- function(theta, iter, f, grad_f, par_list, delta = 0.5, max_treedepth = 10, eps = 1, verbose = TRUE){
  kappa <- 0.75
  t0 <- 10
  gamma <- 0.05
  M_adapt <- par_list$M_adapt
  if(is.null(par_list$M_diag)){
    M_diag <- rep(1, length(theta))
  } else{
    M_diag <- par_list$M_diag
  }
  
  if(iter == 1){
    eps <- find_reasonable_epsilon(theta, f, grad_f, M_diag, eps = eps, verbose = verbose)
    mu <- log(10*eps)
    H <- 0
    eps_bar <- 1
  } else{
    eps <- par_list$eps
    eps_bar <- par_list$eps_bar
    H <- par_list$H
    mu <- par_list$mu
  }
  
  r0 <- rnorm(length(theta), 0, sqrt(M_diag))
  u <- runif(1, 0, exp(f(theta) - 0.5 * sum(r0**2 / M_diag)))
  if(is.nan(u)){
    warning("NUTS: sampled slice u is NaN")
    u <- runif(1, 0, 1e5)
  }
  theta_minus <- theta
  theta_plus <- theta
  r_minus <- r0
  r_plus <- r0
  j=0
  n=1
  s=1
  if(iter > M_adapt){
    eps <- runif(1, 0.9*eps_bar, 1.1*eps_bar)
  }
  while(s == 1){
    # choose direction {-1, 1}
    direction <- sample(c(-1, 1), 1)
    if(direction == -1){
      temp <- build_tree(theta_minus, r_minus, u, direction, j, eps, theta, r0, f, grad_f, M_diag)
      theta_minus <- temp$theta_minus
      r_minus <- temp$r_minus
    } else{
      temp <- build_tree(theta_plus, r_plus, u, direction, j, eps, theta, r0, f, grad_f, M_diag)
      theta_plus <- temp$theta_plus
      r_plus <- temp$r_plus
    }
    if(is.nan(temp$s)) temp$s <- 0
    if(temp$s == 1){
      if(runif(1) < temp$n / n){
        theta <- temp$theta
      }
    }
    n <- n + temp$n
    s <- check_NUTS(temp$s, theta_plus, theta_minus, r_plus, r_minus)
    j <- j + 1
    if(j > max_treedepth){
      warning("NUTS: Reached max tree depth")
      break
    }
  }
  if(iter <= M_adapt){
    H <- (1 - 1/(iter + t0))*H + 1/(iter + t0) * (delta - temp$alpha / temp$n_alpha)
    log_eps <- mu - sqrt(iter)/gamma * H
    eps_bar <- exp(iter**(-kappa) * log_eps + (1 - iter**(-kappa)) * log(eps_bar))
    eps <- exp(log_eps)
  } else{
    eps <- eps_bar
  }
  
  return(list(theta = theta,
              pars = list(eps = eps, eps_bar = eps_bar, H = H, mu = mu, M_adapt = M_adapt, M_diag = M_diag)))
}



