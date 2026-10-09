# Exact sampling methods (Chapter 3) --------------------------------------------------

#' Inverse-transform sampling: if U ~ U(0, 1) and F^{-1} is the (generalised)
#' quantile function of a distribution F, then F^{-1}(U) ~ F.
inverse_cdf_sample <- function(n, quantile_fun) {
  n <- check_count(n, "n"); check_function(quantile_fun, "quantile_fun")
  quantile_fun(stats::runif(n))
}

#' Rejection sampling from a density known up to a constant.
#'
#' Requires an envelope: tilde_pi(x) <= M g(x) for all x, where g is a density
#' we can sample. Propose x ~ g and accept with probability
#'   tilde_pi(x) / (M g(x)).
#' Accepted draws are exact, independent draws from pi = tilde_pi / Z, and the
#' expected acceptance rate is Z / M.
#'
#' @param log_target   vectorised function: log tilde_pi(x).
#' @param r_proposal   function(k) returning k draws from g.
#' @param log_proposal vectorised function: log g(x).
#' @param log_M        log of the envelope constant M.
#' @return list(draws, acceptance_rate, n_proposed)
#' The function stops if it observes tilde_pi(x) > M g(x), because the output
#' would then not follow pi (the envelope condition is violated).
rejection_sample <- function(n, log_target, r_proposal, log_proposal, log_M) {
  n <- check_count(n, "n")
  check_function(log_target, "log_target"); check_function(r_proposal, "r_proposal")
  check_function(log_proposal, "log_proposal"); check_scalar(log_M, "log_M")
  draws <- numeric(0); n_prop <- 0; rate_guess <- 0.5
  while (length(draws) < n) {
    k <- max(100L, ceiling(1.2 * (n - length(draws)) / rate_guess))
    x <- r_proposal(k)
    log_ratio <- log_target(x) - log_proposal(x) - log_M
    log_ratio[is.na(log_ratio)] <- -Inf
    if (any(log_ratio > 1e-10)) {
      course_stop("Envelope violated: tilde_pi(x) > M g(x) at x = %g. Increase M.",
                  x[which.max(log_ratio)])
    }
    accept <- log(stats::runif(k)) < log_ratio
    draws <- c(draws, x[accept]); n_prop <- n_prop + k
    rate_guess <- max(length(draws) / n_prop, 1e-4)
  }
  list(draws = draws[seq_len(n)], acceptance_rate = length(draws) / n_prop, n_proposed = n_prop)
}
