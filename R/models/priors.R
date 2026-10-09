# Prior elicitation and prior predictive tools (Chapter 7) ---------------------------

#' Beta(a, b) parameters with a given mean m and standard deviation s.
#' From E = a / (a + b) and Var = m (1 - m) / (a + b + 1):
#'   a + b = m (1 - m) / s^2 - 1, which must be positive.
beta_from_mean_sd <- function(m, s) {
  check_scalar(m, "m"); check_scalar(s, "s", positive = TRUE)
  if (m <= 0 || m >= 1) course_stop("`m` must lie in (0, 1).")
  n0 <- m * (1 - m) / s^2 - 1
  if (n0 <= 0) course_stop("No Beta distribution has this mean and standard deviation (s too large).")
  c(a = m * n0, b = (1 - m) * n0)
}

#' Beta(a, b) parameters whose quantiles at probabilities p match the values q
#' (two quantiles determine two parameters). Solved numerically on the log scale.
beta_from_quantiles <- function(q, p = c(0.05, 0.95)) {
  if (length(q) != 2L || length(p) != 2L || any(q <= 0 | q >= 1) || q[1] >= q[2]) {
    course_stop("Need two increasing quantiles in (0, 1).")
  }
  obj <- function(lab) {
    ab <- exp(lab)
    sum((stats::qbeta(p, ab[1], ab[2]) - q)^2)
  }
  start <- log(beta_from_mean_sd(mean(q), diff(q) / 3.3))
  fit <- stats::optim(start, obj, control = list(reltol = 1e-14, maxit = 5000))
  ab <- exp(fit$par)
  if (max(abs(stats::qbeta(p, ab[1], ab[2]) - q)) > 1e-4) {
    warning("Quantiles matched only approximately; check the result.")
  }
  c(a = unname(ab[1]), b = unname(ab[2]))
}

#' Simulate from the prior predictive distribution: theta ~ prior, y ~ p(y | theta).
#' @param r_prior function() returning one parameter draw.
#' @param r_data  function(theta) returning one simulated data set (any object).
prior_predictive <- function(n_sims, r_prior, r_data) {
  n_sims <- check_count(n_sims, "n_sims")
  lapply(seq_len(n_sims), function(i) {
    th <- r_prior()
    list(theta = th, y = r_data(th))
  })
}

#' Prior-data conflict check of Evans and Moshonov (2006) for the Beta-Binomial model.
#' The prior predictive of the sufficient statistic y is Beta-Binomial(n, a, b);
#' the check is the prior predictive probability of observing a value whose
#' prior predictive probability is no larger than that of the observed y:
#'   P( m(Y) <= m(y_obs) ),  m = Beta-Binomial pmf.
#' Small values indicate that the observed data are surprising under the prior.
prior_data_conflict_betabinom <- function(y, n, a, b) {
  pmf <- exp(dbetabinom_log(0:n, n, a, b))
  sum(pmf[pmf <= pmf[y + 1] * (1 + 1e-12)])
}
