# Predictive checks and scoring (Chapter 8) ------------------------------------------

#' Posterior predictive p-value: the posterior probability that the replicated
#' test quantity is at least as large as the observed one,
#'   p = P( T(y_rep, theta) >= T(y, theta) | y ).
#' `t_obs` may be a scalar (a test statistic T(y)) or a vector of the same length
#' as `t_rep` (a realised discrepancy T(y, theta_s), evaluated draw by draw).
ppc_pvalue <- function(t_obs, t_rep) {
  if (length(t_obs) != 1L && length(t_obs) != length(t_rep)) {
    course_stop("`t_obs` must have length 1 or the same length as `t_rep`.")
  }
  mean(t_rep >= t_obs)
}

#' Probability integral transform values from predictive draws.
#' @param y     observed values (length n).
#' @param draws S x n matrix of predictive draws, column i for observation i.
#' @param discrete if TRUE, use the randomised PIT  P(Y < y) + U P(Y = y),
#'   which is exactly uniform for a correct discrete predictive distribution.
pit_values <- function(y, draws, discrete = FALSE) {
  draws <- as.matrix(draws)
  if (ncol(draws) != length(y)) course_stop("ncol(draws) must equal length(y).")
  below <- colMeans(sweep(draws, 2, y, "<"))
  if (!discrete) return(colMeans(sweep(draws, 2, y, "<=")))
  equal <- colMeans(sweep(draws, 2, y, "=="))
  below + stats::runif(length(y)) * equal
}

#' Continuous ranked probability score of a predictive distribution represented
#' by draws x, at observation y (lower is better):
#'   CRPS = E|X - y| - (1/2) E|X - X'|,
#' with E|X - X'| computed in O(S log S) from the sorted draws.
crps_sample <- function(y, x) {
  x <- sort(x); S <- length(x)
  e_xy <- mean(abs(x - y))
  e_xx <- 2 * sum((2 * seq_len(S) - S - 1) * x) / S^2
  e_xy - 0.5 * e_xx
}

#' Log pointwise predictive density from a matrix of log-likelihood values
#' (S draws x n observations): lppd_i = log( (1/S) sum_s p(y_i | theta_s) ).
lppd_pointwise <- function(log_lik) apply(as.matrix(log_lik), 2, log_mean_exp)
