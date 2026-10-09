# Benchmark target distributions with known answers ---------------------------
#
# A sampler is validated by running it on targets whose moments are known
# exactly (Chapter 19). Each constructor returns a list with
#   log_density(x)  log of an (un)normalised density; -Inf outside the support
#   grad(x)         gradient of log_density (for HMC / MALA)
#   dim             dimension
#   mean, cov       exact first and second moments when available
#   name            label for tables and plots

#' Univariate normal N(mu, sigma^2).
target_normal <- function(mu = 0, sigma = 1) {
  check_scalar(mu, "mu"); check_scalar(sigma, "sigma", positive = TRUE)
  list(
    name = sprintf("N(%g, %g^2)", mu, sigma), dim = 1L,
    log_density = function(x) -0.5 * ((x - mu) / sigma)^2,
    grad = function(x) -(x - mu) / sigma^2,
    mean = mu, cov = matrix(sigma^2)
  )
}

#' Multivariate normal N(mu, Sigma). The log density omits the constant.
target_mvn <- function(mu, Sigma) {
  d <- length(mu)
  check_spd_matrix(Sigma, "Sigma", d = d)
  P <- solve(Sigma)
  list(
    name = sprintf("%d-dimensional Gaussian", d), dim = d,
    log_density = function(x) { r <- x - mu; -0.5 * sum(r * (P %*% r)) },
    grad = function(x) -as.vector(P %*% (x - mu)),
    mean = mu, cov = Sigma, precision = P
  )
}

#' Two-dimensional Gaussian with unit variances and correlation rho.
target_correlated_gaussian <- function(rho = 0.9) {
  if (!is.numeric(rho) || abs(rho) >= 1) course_stop("`rho` must lie in (-1, 1).")
  t <- target_mvn(c(0, 0), matrix(c(1, rho, rho, 1), 2))
  t$name <- sprintf("bivariate Gaussian, rho = %g", rho)
  t
}

#' Gamma(shape, rate) on (0, Inf): a target with a hard support boundary at 0.
target_gamma <- function(shape = 2, rate = 1) {
  check_scalar(shape, "shape", positive = TRUE); check_scalar(rate, "rate", positive = TRUE)
  list(
    name = sprintf("Gamma(%g, %g)", shape, rate), dim = 1L,
    log_density = function(x) if (x <= 0) -Inf else (shape - 1) * log(x) - rate * x,
    grad = function(x) (shape - 1) / x - rate,
    mean = shape / rate, cov = matrix(shape / rate^2)
  )
}

#' Univariate Gaussian mixture sum_k w_k N(mu_k, s_k^2): a multimodal target.
target_mixture <- function(weights = c(0.5, 0.5), means = c(-3, 3), sds = c(1, 1)) {
  if (abs(sum(weights) - 1) > 1e-12 || any(weights <= 0)) course_stop("`weights` must be positive and sum to one.")
  K <- length(weights)
  if (length(means) != K || length(sds) != K || any(sds <= 0)) course_stop("Mismatched mixture components.")
  lw <- log(weights)
  m <- sum(weights * means)
  v <- sum(weights * (sds^2 + means^2)) - m^2
  list(
    name = "Gaussian mixture", dim = 1L,
    log_density = function(x) log_sum_exp(lw + stats::dnorm(x, means, sds, log = TRUE)),
    grad = function(x) {
      lc <- lw + stats::dnorm(x, means, sds, log = TRUE)
      r <- normalise_log_weights(lc)                   # posterior component responsibilities
      sum(r * (-(x - means) / sds^2))
    },
    mean = m, cov = matrix(v),
    density = function(x) sapply(x, function(xi) sum(weights * stats::dnorm(xi, means, sds)))
  )
}

#' Neal's funnel: v ~ N(0, 3^2), x_i | v ~ N(0, exp(v)), i = 1..(d-1).
#' Coordinates are ordered (v, x_1, ..., x_{d-1}). Exact moments: E v = 0,
#' Var v = 9, E x_i = 0, Var x_i = E exp(v) = exp(9/2).
target_funnel <- function(d = 10) {
  check_count(d, "d", min = 2L)
  list(
    name = sprintf("Neal's funnel (d = %d)", d), dim = d,
    log_density = function(z) {
      v <- z[1]; x <- z[-1]
      stats::dnorm(v, 0, 3, log = TRUE) + sum(stats::dnorm(x, 0, exp(v / 2), log = TRUE))
    },
    grad = function(z) {
      v <- z[1]; x <- z[-1]
      c(-v / 9 - (d - 1) / 2 + 0.5 * sum(x^2) * exp(-v), -x * exp(-v))
    },
    mean = rep(0, d),
    cov = diag(c(9, rep(exp(4.5), d - 1)))
  )
}

#' The "banana" target used by Haario et al. (2001): x1 ~ N(0, 100), x2 = y + b (x1^2 - 100),
#' y ~ N(0, 1). Exact mean (0, 0); Var x2 = 1 + 2 b^2 10^4; Cov(x1, x2) = 0.
target_banana <- function(b = 0.03) {
  check_scalar(b, "b")
  list(
    name = sprintf("banana (b = %g)", b), dim = 2L,
    log_density = function(x) -x[1]^2 / 200 - 0.5 * (x[2] - b * (x[1]^2 - 100))^2,
    grad = function(x) {
      r <- x[2] - b * (x[1]^2 - 100)
      c(-x[1] / 100 + 2 * b * x[1] * r, -r)
    },
    mean = c(0, 0),
    cov = diag(c(100, 1 + 2 * b^2 * 1e4))
  )
}
