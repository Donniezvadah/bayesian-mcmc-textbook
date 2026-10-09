# Slice sampling (Chapter 25) --------------------------------------------------------

#' One univariate slice-sampling update (Neal, 2003): stepping-out followed by
#' shrinkage. Leaves the target invariant for any initial width `w` > 0.
#'
#' 1. Draw the auxiliary level log y = log pi(x0) - Exp(1)   (i.e. y ~ U(0, pi(x0))).
#' 2. Place an interval of width w randomly around x0 and step out until both
#'    ends lie outside the slice {x : log pi(x) > log y} (at most `max_steps`).
#' 3. Sample uniformly from the interval, shrinking it towards x0 on rejection.
slice_update <- function(log_target, x0, w = 1, max_steps = 100L, logp0 = log_target(x0)) {
  log_y <- logp0 - stats::rexp(1)
  left <- x0 - w * stats::runif(1); right <- left + w
  j <- floor(max_steps * stats::runif(1)); k <- (max_steps - 1L) - j
  while (j > 0 && eval_log_density(log_target, left) > log_y) { left <- left - w; j <- j - 1 }
  while (k > 0 && eval_log_density(log_target, right) > log_y) { right <- right + w; k <- k - 1 }
  repeat {
    x1 <- stats::runif(1, left, right)
    logp1 <- eval_log_density(log_target, x1)
    if (logp1 > log_y) return(list(x = x1, logp = logp1))
    if (x1 < x0) left <- x1 else right <- x1
  }
}

#' Univariate slice sampler.
slice_sampler <- function(log_target, initial, n_iter, w = 1, max_steps = 100L) {
  check_function(log_target, "log_target"); check_scalar(initial, "initial")
  n_iter <- check_count(n_iter, "n_iter"); check_scalar(w, "w", positive = TRUE)
  x <- initial; logp <- eval_log_density(log_target, x)
  if (!is.finite(logp)) course_stop("Initial value must have finite log target density.")
  draws <- numeric(n_iter)
  for (t in seq_len(n_iter)) {
    s <- slice_update(log_target, x, w, max_steps, logp)
    x <- s$x; logp <- s$logp; draws[t] <- x
  }
  draws
}

#' Elliptical slice sampling (Murray, Adams and MacKay, 2010) for posteriors of
#' the form  pi(f) propto L(f) N(f | 0, Sigma).
#' Each update draws nu ~ N(0, Sigma), defines the ellipse
#'   f'(phi) = f cos(phi) + nu sin(phi),
#' and slice-samples the angle phi with a shrinking bracket. There are no tuning
#' parameters, and every update produces a move (no rejections).
elliptical_slice <- function(log_lik, initial, n_iter, Sigma) {
  check_function(log_lik, "log_lik"); check_finite_vector(initial, "initial")
  n_iter <- check_count(n_iter, "n_iter")
  d <- length(initial); check_spd_matrix(Sigma, "Sigma", d = d)
  U <- chol(Sigma)
  f <- initial; ll <- eval_log_density(log_lik, f)
  draws <- matrix(NA_real_, n_iter, d)
  for (t in seq_len(n_iter)) {
    nu <- as.vector(crossprod(U, stats::rnorm(d)))
    log_y <- ll - stats::rexp(1)
    phi <- stats::runif(1, 0, 2 * pi); lo <- phi - 2 * pi; hi <- phi
    repeat {
      f_new <- f * cos(phi) + nu * sin(phi)
      ll_new <- eval_log_density(log_lik, f_new)
      if (ll_new > log_y) break
      if (phi < 0) lo <- phi else hi <- phi
      phi <- stats::runif(1, lo, hi)
    }
    f <- f_new; ll <- ll_new; draws[t, ] <- f
  }
  draws
}

#' Preconditioned Crank-Nicolson (pCN) sampler (Cotter et al., 2013) for
#'   pi(f) propto L(f) N(f | 0, Sigma).
#' Proposal f' = sqrt(1 - beta^2) f + beta xi, xi ~ N(0, Sigma), which is reversible
#' with respect to the prior N(0, Sigma); hence the Metropolis-Hastings acceptance
#' probability involves only the likelihood ratio, min{1, L(f') / L(f)}, and the
#' acceptance rate does not degenerate as the dimension of f grows.
pcn_sampler <- function(log_lik, initial, n_iter, Sigma, beta = 0.3) {
  check_function(log_lik, "log_lik"); check_finite_vector(initial, "initial")
  n_iter <- check_count(n_iter, "n_iter")
  if (!(beta > 0 && beta <= 1)) course_stop("`beta` must lie in (0, 1].")
  d <- length(initial); check_spd_matrix(Sigma, "Sigma", d = d); U <- chol(Sigma)
  f <- initial; ll <- eval_log_density(log_lik, f)
  draws <- matrix(NA_real_, n_iter, d); acc <- 0L
  for (t in seq_len(n_iter)) {
    prop <- sqrt(1 - beta^2) * f + beta * as.vector(crossprod(U, stats::rnorm(d)))
    ll_p <- eval_log_density(log_lik, prop)
    if (log(stats::runif(1)) < ll_p - ll) { f <- prop; ll <- ll_p; acc <- acc + 1L }
    draws[t, ] <- f
  }
  list(draws = draws, acceptance_rate = acc / n_iter)
}
