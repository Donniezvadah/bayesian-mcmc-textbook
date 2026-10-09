# Predictive model comparison (Chapter 23) --------------------------------------------
#
# All functions take an S x n matrix of pointwise log-likelihood values,
# log_lik[s, i] = log p(y_i | theta^(s)), computed from S posterior draws.

#' Widely applicable information criterion (Watanabe, 2010), on the elpd scale:
#'   lppd   = sum_i log( (1/S) sum_s p(y_i | theta_s) )
#'   p_waic = sum_i Var_s( log p(y_i | theta_s) )
#'   elpd_waic = lppd - p_waic   (WAIC = -2 elpd_waic)
waic <- function(log_lik) {
  log_lik <- as.matrix(log_lik)
  lppd_i <- apply(log_lik, 2, log_mean_exp)
  p_i <- apply(log_lik, 2, stats::var)
  elpd_i <- lppd_i - p_i
  list(elpd_waic = sum(elpd_i), p_waic = sum(p_i), se = sqrt(ncol(log_lik) * stats::var(elpd_i)),
       pointwise = elpd_i)
}

#' Plain importance-sampling leave-one-out (Gelfand-style): with posterior draws as
#' proposals for the leave-one-out posterior, the weights are 1 / p(y_i | theta_s), and
#'   p(y_i | y_{-i}) ~= 1 / ( (1/S) sum_s 1 / p(y_i | theta_s) ),
#' a harmonic mean that can have infinite variance. Shown for teaching; use PSIS
#' (the loo package) in practice. Also returns Kish's ESS of the weights for each i.
loo_is <- function(log_lik) {
  log_lik <- as.matrix(log_lik)
  elpd_i <- -apply(-log_lik, 2, log_mean_exp)
  ess_i <- apply(-log_lik, 2, function(lw) 1 / sum(normalise_log_weights(lw)^2))
  list(elpd_loo = sum(elpd_i), se = sqrt(ncol(log_lik) * stats::var(elpd_i)), pointwise = elpd_i, ess = ess_i)
}

#' Difference in elpd between two models with its standard error, computed from
#' pointwise values (the paired differences across observations).
elpd_difference <- function(pointwise_a, pointwise_b) {
  d <- pointwise_a - pointwise_b
  c(elpd_diff = sum(d), se_diff = sqrt(length(d) * stats::var(d)))
}
