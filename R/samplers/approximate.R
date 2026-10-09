# Approximate and scalable inference (Chapter 27) --------------------------------------
#
# These methods trade exactness for applicability or speed. Each targets a
# *different* distribution from the posterior, and the chapter quantifies how.

#' Rejection ABC (Pritchard et al., 1999; Beaumont et al., 2002).
#' Simulates theta ~ prior and y* ~ p(. | theta), and keeps the fraction `keep` of
#' parameter draws whose summaries S(y*) are closest to the observed summaries.
#' The accepted draws target pi_eps(theta) propto pi(theta) P(rho(S(Y*), s_obs) <= eps | theta),
#' with eps the `keep`-quantile of the simulated distances.
#'
#' @param r_prior function(n) returning n prior draws (vector, or n x d matrix).
#' @param simulate function(theta) returning one simulated data set.
#' @param summary_fn function(data) returning a numeric vector of summaries.
#' @param s_obs observed summaries, summary_fn(y_obs).
#' @param scale divide each summary's differences by the MAD of its simulated values.
abc_rejection <- function(n_sims, r_prior, simulate, summary_fn, s_obs, keep = 0.01, scale = TRUE) {
  n_sims <- check_count(n_sims, "n_sims"); check_probability(keep, "keep")
  check_function(r_prior, "r_prior"); check_function(simulate, "simulate"); check_function(summary_fn, "summary_fn")
  theta <- r_prior(n_sims)
  theta_mat <- as.matrix(theta)
  if (nrow(theta_mat) != n_sims) course_stop("`r_prior(n)` must return n draws.")
  S <- t(vapply(seq_len(n_sims), function(i) {
    th <- if (is.matrix(theta)) theta[i, ] else theta[i]
    as.numeric(summary_fn(simulate(th)))
  }, numeric(length(s_obs))))
  if (length(s_obs) == 1L) S <- matrix(S, ncol = 1)
  sc <- if (scale) apply(S, 2, stats::mad) else rep(1, ncol(S))
  sc[!is.finite(sc) | sc == 0] <- 1
  dist <- sqrt(rowSums(sweep(sweep(S, 2, s_obs), 2, sc, "/")^2))
  n_keep <- max(1L, floor(keep * n_sims))
  idx <- order(dist)[seq_len(n_keep)]
  list(draws = if (is.matrix(theta)) theta[idx, , drop = FALSE] else theta[idx],
       eps = max(dist[idx]), distances = dist, summaries = S, n_sims = n_sims)
}

#' Mean-field Gaussian automatic differentiation variational inference (ADVI;
#' Kucukelbir et al., 2017) on an unconstrained parameter space.
#' q(x) = prod_i N(x_i | m_i, exp(2 omega_i)); maximises the ELBO
#'   L(m, omega) = E_q[log p(x, y)] + sum_i omega_i + const
#' by stochastic gradient ascent with reparameterisation gradients
#'   dL/dm = E[grad log p(m + s * z)],  dL/domega = E[grad log p(m + s * z) * z * s] + 1,
#' Adam step sizes, and Polyak averaging over the last half of the iterations.
advi_meanfield <- function(log_density, grad, initial, n_iter = 5000L, n_mc = 5L, step = 0.05,
                           init_log_sd = -1) {
  check_function(log_density, "log_density"); check_function(grad, "grad")
  check_finite_vector(initial, "initial"); n_iter <- check_count(n_iter, "n_iter")
  n_mc <- check_count(n_mc, "n_mc"); check_scalar(step, "step", positive = TRUE)
  d <- length(initial)
  lam <- c(initial, rep(init_log_sd, d))
  mom1 <- mom2 <- numeric(2 * d); b1 <- 0.9; b2 <- 0.999
  elbo <- numeric(n_iter); avg <- numeric(2 * d); n_avg <- 0L; avg_from <- floor(n_iter / 2) + 1L
  for (t in seq_len(n_iter)) {
    m <- lam[1:d]; s <- exp(lam[d + 1:d])
    g <- numeric(2 * d); lp <- 0
    for (k in seq_len(n_mc)) {
      z <- stats::rnorm(d); x <- m + s * z
      gx <- grad(x)
      g <- g + c(gx, gx * z * s) / n_mc
      lp <- lp + eval_log_density(log_density, x) / n_mc
    }
    g[d + 1:d] <- g[d + 1:d] + 1
    elbo[t] <- lp + sum(lam[d + 1:d])
    mom1 <- b1 * mom1 + (1 - b1) * g; mom2 <- b2 * mom2 + (1 - b2) * g^2
    lam <- lam + step * (mom1 / (1 - b1^t)) / (sqrt(mom2 / (1 - b2^t)) + 1e-8)
    if (t >= avg_from) { n_avg <- n_avg + 1L; avg <- avg + (lam - avg) / n_avg }
  }
  mean <- avg[1:d]; sd <- exp(avg[d + 1:d])
  list(mean = mean, sd = sd, elbo = elbo,
       draw = function(n) matrix(stats::rnorm(n * d, rep(mean, each = n), rep(sd, each = n)), n, d),
       log_q = function(x) sum(stats::dnorm(x, mean, sd, log = TRUE)))
}

#' Stochastic gradient Langevin dynamics (Welling and Teh, 2011) with a constant step.
#'   theta' = theta + (step / 2) [grad log prior + (n / b) sum_{i in B} grad log p(y_i | theta)] + N(0, step)
#' The minibatch B of size b is drawn with replacement. There is no Metropolis
#' correction, so the chain's stationary distribution differs from the posterior:
#' both the discretisation and the minibatch noise inflate its spread.
#'
#' @param grad_log_lik function(theta, data_rows) returning the *sum* of the
#'   log-likelihood gradients over the rows supplied.
sgld_sampler <- function(grad_log_prior, grad_log_lik, data, initial, n_iter, batch_size, step) {
  check_function(grad_log_prior, "grad_log_prior"); check_function(grad_log_lik, "grad_log_lik")
  check_finite_vector(initial, "initial"); n_iter <- check_count(n_iter, "n_iter")
  check_scalar(step, "step", positive = TRUE)
  data <- as.matrix(data); n <- nrow(data); b <- check_count(batch_size, "batch_size")
  d <- length(initial); theta <- initial
  draws <- matrix(NA_real_, n_iter, d)
  for (t in seq_len(n_iter)) {
    idx <- sample.int(n, b, replace = TRUE)
    g <- grad_log_prior(theta) + (n / b) * grad_log_lik(theta, data[idx, , drop = FALSE])
    theta <- theta + 0.5 * step * g + sqrt(step) * stats::rnorm(d)
    draws[t, ] <- theta
  }
  draws
}

#' Stationary mean and variance of SGLD for the Gaussian-mean model
#' y_i ~ N(theta, 1), theta ~ N(0, s0^2), with minibatches drawn with replacement
#' (@prp-sgld-variance in Chapter 27). Requires step < 4 / P.
sgld_gaussian_stationary <- function(y, s0, batch_size, step) {
  n <- length(y); P <- n + 1 / s0^2; s2y <- mean((y - mean(y))^2)
  if (step >= 4 / P) course_stop("The SGLD recursion is unstable: need step < 4 / P.")
  c(mean = sum(y) / P,
    var = (1 + step * n^2 * s2y / (4 * batch_size)) / (P * (1 - step * P / 4)),
    posterior_var = 1 / P, ar_coefficient = 1 - step * P / 2)
}
