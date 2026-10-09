# Adaptive and gradient-based Metropolis methods (Chapter 16) ---------------------

#' Adaptive Metropolis (Haario, Saksman and Tamminen, 2001).
#'
#' Random-walk Metropolis whose proposal covariance at iteration t is
#'   s_d * (Cov(x_0, ..., x_{t-1}) + eps I),  s_d = 2.38^2 / d,
#' once t exceeds `adapt_start`. The empirical covariance is updated
#' recursively in O(d^2) per iteration. Adaptation diminishes automatically
#' (the covariance estimate changes by O(1/t)), which together with the eps I
#' regularisation is the setting of the ergodicity result of Haario et al.
#' Setting `adapt_stop` freezes the proposal so that later draws come from a
#' fixed Markov kernel (the safest practice).
adaptive_metropolis <- function(log_target, initial, n_iter, adapt_start = 500,
                                adapt_stop = n_iter, init_sd = 0.1, eps = 1e-6) {
  check_function(log_target, "log_target"); check_finite_vector(initial, "initial")
  n_iter <- check_count(n_iter, "n_iter")
  d <- length(initial); s_d <- 2.38^2 / d
  x <- initial; logp <- eval_log_density(log_target, x)
  if (!is.finite(logp)) course_stop("Initial value must have finite log target density.")
  draws <- matrix(NA_real_, n_iter, d)
  # running mean and (unnormalised) scatter matrix of the chain so far
  m <- x; S <- matrix(0, d, d); prop_cov <- diag(init_sd^2, d)
  accepted <- 0L
  for (t in seq_len(n_iter)) {
    if (t > adapt_start && t <= adapt_stop) {
      prop_cov <- s_d * (S / (t - 1) + eps * diag(d))
    }
    y <- x + as.vector(crossprod(chol(prop_cov), stats::rnorm(d)))
    logp_y <- eval_log_density(log_target, y)
    if (log(stats::runif(1)) < logp_y - logp) { x <- y; logp <- logp_y; accepted <- accepted + 1L }
    draws[t, ] <- x
    # Welford-type recursive update of mean and scatter (numerically stable).
    delta <- x - m
    m <- m + delta / (t + 1)
    S <- S + tcrossprod(delta, x - m)
  }
  list(draws = draws, acceptance_rate = accepted / n_iter, final_proposal_cov = prop_cov)
}

#' Metropolis-adjusted Langevin algorithm (MALA).
#'
#' Proposal: x' = x + (h/2) grad log pi(x) + sqrt(h) Z, a single Euler step of
#' the Langevin diffusion with invariant law pi. The Euler step does not preserve
#' pi exactly, so the asymmetric Hastings correction is required:
#'   log q(x' | x) = -||x' - x - (h/2) grad log pi(x)||^2 / (2h) + const.
mala_sampler <- function(log_target, grad_log_target, initial, n_iter, step) {
  check_function(log_target, "log_target"); check_function(grad_log_target, "grad_log_target")
  check_finite_vector(initial, "initial"); n_iter <- check_count(n_iter, "n_iter")
  check_scalar(step, "step", positive = TRUE)
  log_q <- function(to, from, g_from) -sum((to - from - 0.5 * step * g_from)^2) / (2 * step)
  x <- initial; logp <- eval_log_density(log_target, x); g <- grad_log_target(x)
  if (!is.finite(logp)) course_stop("Initial value must have finite log target density.")
  d <- length(x); draws <- matrix(NA_real_, n_iter, d); accepted <- 0L
  for (t in seq_len(n_iter)) {
    y <- x + 0.5 * step * g + sqrt(step) * stats::rnorm(d)
    logp_y <- eval_log_density(log_target, y)
    if (is.finite(logp_y)) {
      g_y <- grad_log_target(y)
      log_alpha <- logp_y + log_q(x, y, g_y) - logp - log_q(y, x, g)
      if (log(stats::runif(1)) < log_alpha) {
        x <- y; logp <- logp_y; g <- g_y; accepted <- accepted + 1L
      }
    }
    draws[t, ] <- x
  }
  list(draws = draws, acceptance_rate = accepted / n_iter)
}

#' Random-walk Metropolis with Robbins-Monro adaptation of the proposal scale
#' during warm-up, followed by sampling with the scale frozen.
#'
#' During warm-up: log s_{t+1} = log s_t + t^{-kappa} (accept_t - target_accept),
#' a stochastic-approximation recursion (Robbins and Monro, 1951) whose fixed
#' point has expected acceptance probability equal to `target_accept`. The
#' proposal is s * chol(prop_cov). Warm-up draws are discarded; the retained
#' draws come from a fixed, pi-invariant kernel.
rwm_adapt_scale <- function(log_target, initial, n_warmup, n_iter, prop_cov = diag(length(initial)),
                            target_accept = 0.234, kappa = 0.6, init_scale = 1) {
  check_function(log_target, "log_target"); check_finite_vector(initial, "initial")
  n_warmup <- check_count(n_warmup, "n_warmup"); n_iter <- check_count(n_iter, "n_iter")
  d <- length(initial); check_spd_matrix(prop_cov, "prop_cov", d = d)
  U <- chol(prop_cov)
  x <- initial; lp <- eval_log_density(log_target, x)
  if (!is.finite(lp)) course_stop("Initial value must have finite log target density.")
  log_s <- log(init_scale); scale_trace <- numeric(n_warmup)
  for (t in seq_len(n_warmup)) {
    y <- x + exp(log_s) * as.vector(crossprod(U, stats::rnorm(d)))
    lp_y <- eval_log_density(log_target, y)
    a <- if (is.finite(lp_y)) min(1, exp(lp_y - lp)) else 0
    if (stats::runif(1) < a) { x <- y; lp <- lp_y }
    log_s <- log_s + t^(-kappa) * (a - target_accept)
    scale_trace[t] <- exp(log_s)
  }
  out <- rw_metropolis(log_target, x, n_iter, proposal_cov = exp(2 * log_s) * prop_cov)
  c(out, list(scale = exp(log_s), scale_trace = scale_trace))
}
