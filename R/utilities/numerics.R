# Numerically stable building blocks -----------------------------------------
#
# Bayesian computation is done on the log scale: densities of realistic data
# sets underflow double precision (the smallest positive double is about
# 4.9e-324, i.e. a log value of about -744). The functions below combine
# log-scale quantities without leaving the log scale.

#' log(sum(exp(x))) computed without overflow or underflow.
#'
#' Uses the identity log sum exp(x_i) = m + log sum exp(x_i - m), m = max x_i.
#' Every term exp(x_i - m) lies in (0, 1] and at least one equals 1, so the
#' sum lies in [1, n] and its logarithm is well conditioned.
#' Returns -Inf when all inputs are -Inf (an empty sum of zero densities).
log_sum_exp <- function(x) {
  if (!is.numeric(x) || length(x) == 0L) course_stop("`x` must be a non-empty numeric vector.")
  if (anyNA(x)) return(NA_real_)
  m <- max(x)
  if (m == -Inf) return(-Inf)
  if (m == Inf) return(Inf)
  m + log(sum(exp(x - m)))
}

#' log(mean(exp(x))), e.g. the log of a Monte Carlo average of likelihoods.
log_mean_exp <- function(x) log_sum_exp(x) - log(length(x))

#' Convert log weights to normalised weights that sum to one (softmax).
normalise_log_weights <- function(log_w) {
  lse <- log_sum_exp(log_w)
  if (!is.finite(lse)) course_stop("All log weights are -Inf or non-finite; cannot normalise.")
  exp(log_w - lse)
}

#' log(1 - exp(-a)) for a > 0, accurate for both small and large a.
#'
#' Following Maechler (2012): for a <= log 2 use log(-expm1(-a)); otherwise
#' use log1p(-exp(-a)). Each branch avoids catastrophic cancellation in its
#' own range.
log1mexp <- function(a) {
  if (any(a < 0, na.rm = TRUE)) course_stop("log1mexp() requires a >= 0.")
  ifelse(a <= log(2), log(-expm1(-a)), log1p(-exp(-a)))
}

#' log(1 + exp(x)) (the "softplus" function), stable for large |x|.
log1pexp <- function(x) {
  ifelse(x > 0, x + log1p(exp(-x)), log1p(exp(x)))
}

#' Logistic function and its logarithm.
inv_logit <- function(x) stats::plogis(x)
log_inv_logit <- function(x) stats::plogis(x, log.p = TRUE)          # log(1 / (1 + e^{-x}))
log1m_inv_logit <- function(x) stats::plogis(x, lower.tail = FALSE, log.p = TRUE)  # log(1 - inv_logit(x))
logit <- function(p) stats::qlogis(p)

#' Central finite-difference gradient. Used only to *check* analytic gradients
#' (for HMC and MALA); never inside a sampler. Error is O(h^2) plus rounding
#' error of order eps * |f| / h, so h around 1e-5 balances the two terms for
#' well-scaled functions.
numerical_gradient <- function(f, x, h = 1e-5) {
  g <- numeric(length(x))
  for (i in seq_along(x)) {
    e <- numeric(length(x)); e[i] <- h
    g[i] <- (f(x + e) - f(x - e)) / (2 * h)
  }
  g
}

#' Composite Simpson rule on an equally spaced grid with an odd number of points.
#' Exact for cubic polynomials; error O(h^4) for smooth integrands.
simpson <- function(y, x) {
  n <- length(x)
  if (n < 3L || n %% 2L == 0L) course_stop("Simpson's rule needs an odd number (>= 3) of grid points.")
  h <- diff(range(x)) / (n - 1)
  if (max(abs(diff(x) - h)) > 1e-8 * max(1, abs(h))) course_stop("Grid must be equally spaced.")
  w <- rep(2, n)                       # weights 1, 4, 2, 4, ..., 2, 4, 1
  w[seq(2, n - 1, by = 2)] <- 4
  w[c(1, n)] <- 1
  sum(w * y) * h / 3
}

#' Normalise an unnormalised log density evaluated on a 1-D equally spaced grid.
#' Returns the density values and the log normalising constant. This "grid
#' approximation" is exact in the limit of a fine grid covering the support and
#' is the reference against which we validate samplers in one dimension.
grid_posterior <- function(log_density_values, grid) {
  m <- max(log_density_values)
  unnorm <- exp(log_density_values - m)
  if (length(grid) %% 2L == 1L) {
    z <- simpson(unnorm, grid)
  } else {
    z <- sum(diff(grid) * (head(unnorm, -1) + tail(unnorm, -1)) / 2)   # trapezoid
  }
  list(density = unnorm / z, log_normalising_constant = m + log(z))
}

#' Draw from a multivariate normal N(mu, Sigma) using the Cholesky factor.
#' If Sigma = L L^T and z ~ N(0, I), then mu + L z ~ N(mu, Sigma).
#' Cost: O(d^3) for the factorisation once, then O(n d^2).
rmvnorm_chol <- function(n, mu, Sigma) {
  check_count(n, "n")
  d <- length(mu)
  check_spd_matrix(Sigma, "Sigma", d = d)
  U <- chol(Sigma)                                   # upper triangular, Sigma = U^T U
  z <- matrix(stats::rnorm(n * d), nrow = n, ncol = d)
  sweep(z %*% U, 2, mu, "+")
}

#' Log density of N(mu, Sigma) at the rows of x, via the Cholesky factor.
#' log N(x | mu, Sigma) = -d/2 log(2 pi) - sum(log diag(U)) - ||U^{-T}(x - mu)||^2 / 2.
dmvnorm_log <- function(x, mu, Sigma) {
  x <- matrix(x, ncol = length(mu))
  U <- chol(Sigma)
  r <- backsolve(U, t(x) - mu, transpose = TRUE)     # solves U^T r = x - mu
  -0.5 * ncol(x) * log(2 * pi) - sum(log(diag(U))) - 0.5 * colSums(r^2)
}

#' Evaluate `code` with a temporary random seed, restoring the previous RNG
#' state afterwards (so a demonstration does not disturb later random numbers).
with_seed <- function(seed, code) {
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had_seed) old <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit({
    if (had_seed) assign(".Random.seed", old, envir = globalenv())
    else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
      rm(".Random.seed", envir = globalenv())
  })
  set.seed(seed)
  code
}

#' Draws from a multivariate Student t with nu degrees of freedom, location mu
#' and scale matrix Sigma: mu + Z / sqrt(W / nu), Z ~ N(0, Sigma), W ~ chi^2_nu.
#' Heavier tails than the normal make it a safer importance-sampling proposal.
rmvt_chol <- function(n, mu, Sigma, nu) {
  check_scalar(nu, "nu", positive = TRUE)
  z <- rmvnorm_chol(n, rep(0, length(mu)), Sigma)
  w <- stats::rchisq(n, df = nu)
  sweep(z / sqrt(w / nu), 2, mu, "+")
}

#' Log density of the multivariate t at the rows of x.
dmvt_log <- function(x, mu, Sigma, nu) {
  d <- length(mu); x <- matrix(x, ncol = d)
  U <- chol(Sigma)
  r <- backsolve(U, t(x) - mu, transpose = TRUE)
  q <- colSums(r^2)
  lgamma((nu + d) / 2) - lgamma(nu / 2) - 0.5 * d * log(nu * pi) - sum(log(diag(U))) -
    0.5 * (nu + d) * log1p(q / nu)
}

#' Gauss-Hermite quadrature for integrals against the standard normal density:
#'   int f(u) phi(u) du  ~=  sum_k w_k f(x_k).
#' Nodes and weights from the Golub-Welsch algorithm: eigenvalues of the symmetric
#' tridiagonal Jacobi matrix of the (probabilists') Hermite polynomials, with
#' off-diagonal entries sqrt(k), k = 1..n-1. Exact for polynomials of degree <= 2n - 1.
gauss_hermite <- function(n) {
  n <- check_count(n, "n")
  if (n == 1L) return(list(nodes = 0, weights = 1))
  J <- matrix(0, n, n)
  off <- sqrt(seq_len(n - 1))
  J[cbind(1:(n - 1), 2:n)] <- off; J[cbind(2:n, 1:(n - 1))] <- off
  e <- eigen(J, symmetric = TRUE)
  ord <- order(e$values)
  list(nodes = e$values[ord], weights = e$vectors[1, ord]^2)
}
