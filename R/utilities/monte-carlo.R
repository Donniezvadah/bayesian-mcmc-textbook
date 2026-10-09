# Independent Monte Carlo estimation (Chapter 9) -------------------------------------

#' Monte Carlo estimate of E[f(X)] from independent evaluations fx = f(X_1..X_N),
#' with its Monte Carlo standard error s / sqrt(N) and a CLT-based interval.
#' Valid only for INDEPENDENT draws; for MCMC output use mcse_mean() (Chapter 17).
mc_estimate <- function(fx, level = 0.95) {
  check_finite_vector(fx, "fx")
  N <- length(fx)
  if (N < 2L) course_stop("Need at least two evaluations.")
  est <- mean(fx); se <- stats::sd(fx) / sqrt(N)
  z <- stats::qnorm(1 - (1 - level) / 2)
  c(estimate = est, mcse = se, lower = est - z * se, upper = est + z * se, n = N)
}

#' Control-variate estimator of E[f(X)] using h(X) with known mean mu_h.
#' beta is estimated by least squares, beta = Cov(f, h) / Var(h); the estimator
#' mean(f - beta (h - mu_h)) has variance (1 - rho^2) Var(f) / N asymptotically,
#' with rho = Cor(f, h). Estimating beta from the same draws adds an O(1/N) bias.
control_variate_estimate <- function(fx, hx, mu_h) {
  check_finite_vector(fx, "fx"); check_finite_vector(hx, "hx", len = length(fx))
  check_scalar(mu_h, "mu_h")
  beta <- stats::cov(fx, hx) / stats::var(hx)
  adj <- fx - beta * (hx - mu_h)
  c(estimate = mean(adj), mcse = stats::sd(adj) / sqrt(length(adj)), beta = beta,
    variance_ratio = stats::var(adj) / stats::var(fx))
}

#' Antithetic estimator of E[f(U)], U ~ U(0,1)^d, using pairs (U, 1 - U).
#' Returns the estimate and its standard error computed from the N/2 pair means.
antithetic_estimate <- function(f, n_pairs, d = 1L) {
  check_function(f, "f"); n_pairs <- check_count(n_pairs, "n_pairs")
  u <- matrix(stats::runif(n_pairs * d), n_pairs, d)
  pair_means <- (apply(u, 1, f) + apply(1 - u, 1, f)) / 2
  c(estimate = mean(pair_means), mcse = stats::sd(pair_means) / sqrt(n_pairs), n_evaluations = 2 * n_pairs)
}
