# Metropolis-Hastings samplers (Chapter 12) -------------------------------------
#
# Educational implementations: every step of the algorithm is visible, inputs
# are validated, and all acceptance decisions are made on the log scale. They
# are designed for clarity and testing, not speed; production work should use
# established software such as Stan (Chapter 16).

#' Univariate random-walk Metropolis sampler.
#'
#' Target: a density pi on the real line known up to a constant through
#' `log_target`. Proposal: theta' = theta + proposal_sd * Z, Z ~ N(0, 1), which
#' is symmetric, so the Hastings correction cancels and
#'   log alpha = log pi(theta') - log pi(theta).
#'
#' @param log_target  function returning log pi(theta) up to an additive
#'                    constant; -Inf outside the support (NaN is treated as -Inf).
#' @param initial     starting value; must have finite log target density.
#' @param n_iter      number of iterations (each one produces one draw).
#' @param proposal_sd standard deviation of the Gaussian increment (> 0).
#' @return list with `draws` (length n_iter), `acceptance_rate`, `log_target`
#'         values along the chain, and the settings used.
#' Complexity: O(n_iter) target evaluations; memory O(n_iter).
metropolis_hastings <- function(log_target, initial, n_iter, proposal_sd) {
  check_function(log_target, "log_target")
  check_scalar(initial, "initial")
  n_iter <- check_count(n_iter, "n_iter")
  check_scalar(proposal_sd, "proposal_sd", positive = TRUE)

  current <- initial
  current_logp <- eval_log_density(log_target, current)
  if (!is.finite(current_logp)) {
    course_stop("Initial value must have finite log target density (got %s).", current_logp)
  }

  draws <- numeric(n_iter)
  logp_trace <- numeric(n_iter)
  accepted <- 0L

  for (i in seq_len(n_iter)) {
    proposal <- current + proposal_sd * stats::rnorm(1)
    proposal_logp <- eval_log_density(log_target, proposal)   # -Inf outside the support

    # Symmetric proposal: q(theta | theta') / q(theta' | theta) = 1.
    log_alpha <- proposal_logp - current_logp
    # log(U) <= 0, so comparing with log_alpha is equivalent to comparing with
    # min(0, log_alpha). If proposal_logp = -Inf then log_alpha = -Inf: rejected.
    if (log(stats::runif(1)) < log_alpha) {
      current <- proposal
      current_logp <- proposal_logp
      accepted <- accepted + 1L
    }
    draws[i] <- current          # on rejection the current state is repeated
    logp_trace[i] <- current_logp
  }

  list(draws = draws,
       acceptance_rate = accepted / n_iter,
       log_target = logp_trace,
       settings = list(initial = initial, n_iter = n_iter, proposal_sd = proposal_sd))
}

#' Multivariate random-walk Metropolis with Gaussian proposal N(theta, Sigma_q).
#'
#' @param proposal_cov d x d proposal covariance (or a positive scalar s, meaning s^2 I).
#' @return list with `draws` (n_iter x d matrix), `acceptance_rate`, `log_target`.
#' Complexity per iteration: one target evaluation plus O(d^2) for the proposal.
rw_metropolis <- function(log_target, initial, n_iter, proposal_cov) {
  check_function(log_target, "log_target")
  check_finite_vector(initial, "initial")
  n_iter <- check_count(n_iter, "n_iter")
  d <- length(initial)
  if (length(proposal_cov) == 1L) {
    check_scalar(proposal_cov, "proposal_cov", positive = TRUE)
    proposal_cov <- diag(proposal_cov^2, d)
  }
  check_spd_matrix(proposal_cov, "proposal_cov", d = d)
  U <- chol(proposal_cov)                         # proposal increment = t(U) %*% z

  current <- initial
  current_logp <- eval_log_density(log_target, current)
  if (!is.finite(current_logp)) course_stop("Initial value must have finite log target density.")

  draws <- matrix(NA_real_, n_iter, d, dimnames = list(NULL, names(initial)))
  logp_trace <- numeric(n_iter)
  accepted <- 0L
  for (i in seq_len(n_iter)) {
    proposal <- current + as.vector(crossprod(U, stats::rnorm(d)))
    proposal_logp <- eval_log_density(log_target, proposal)
    if (log(stats::runif(1)) < proposal_logp - current_logp) {
      current <- proposal; current_logp <- proposal_logp; accepted <- accepted + 1L
    }
    draws[i, ] <- current
    logp_trace[i] <- current_logp
  }
  list(draws = draws, acceptance_rate = accepted / n_iter, log_target = logp_trace)
}

#' General Metropolis-Hastings with an arbitrary (possibly asymmetric) proposal.
#'
#' @param r_proposal function(current) returning a proposed state.
#' @param log_q      function(to, from) returning log q(to | from), up to a
#'                   constant that does not depend on `to` or `from`.
#' The acceptance ratio is
#'   log alpha = log pi(x') + log q(x | x') - log pi(x) - log q(x' | x).
mh_sampler <- function(log_target, initial, n_iter, r_proposal, log_q) {
  check_function(log_target, "log_target")
  check_function(r_proposal, "r_proposal"); check_function(log_q, "log_q")
  check_finite_vector(initial, "initial")
  n_iter <- check_count(n_iter, "n_iter")

  current <- initial
  current_logp <- eval_log_density(log_target, current)
  if (!is.finite(current_logp)) course_stop("Initial value must have finite log target density.")
  d <- length(initial)
  draws <- matrix(NA_real_, n_iter, d)
  accepted <- 0L
  for (i in seq_len(n_iter)) {
    proposal <- r_proposal(current)
    proposal_logp <- eval_log_density(log_target, proposal)
    if (is.finite(proposal_logp)) {
      log_alpha <- proposal_logp + log_q(current, proposal) -
        current_logp - log_q(proposal, current)
      if (log(stats::runif(1)) < log_alpha) {
        current <- proposal; current_logp <- proposal_logp; accepted <- accepted + 1L
      }
    }
    draws[i, ] <- current
  }
  list(draws = if (d == 1L) as.vector(draws) else draws, acceptance_rate = accepted / n_iter)
}

#' Independence Metropolis-Hastings: proposals are drawn from a fixed density g
#' that does not depend on the current state. The acceptance ratio compares
#' importance weights w(x) = pi(x) / g(x):  alpha = min{1, w(x') / w(x)}.
independence_sampler <- function(log_target, initial, n_iter, r_g, log_g) {
  mh_sampler(log_target, initial, n_iter,
             r_proposal = function(current) r_g(1),
             log_q = function(to, from) log_g(to))
}

#' Exact transition matrix of Metropolis-Hastings on a finite state space.
#'
#' For a target vector `pi` (possibly unnormalised) and proposal matrix Q
#' (rows sum to one), returns P with P[i, j] = Q[i, j] alpha(i, j) for i != j
#' and P[i, i] = 1 - sum_{j != i} P[i, j]. Used to *verify* detailed balance
#' and invariance exactly (Chapters 10 and 12).
mh_transition_matrix <- function(pi, Q) {
  K <- length(pi)
  if (!is.matrix(Q) || any(dim(Q) != K)) course_stop("Q must be a K x K matrix.")
  if (any(abs(rowSums(Q) - 1) > 1e-12)) course_stop("Rows of Q must sum to one.")
  P <- matrix(0, K, K)
  for (i in seq_len(K)) for (j in seq_len(K)) {
    if (i != j && Q[i, j] > 0) {
      ratio <- (pi[j] * Q[j, i]) / (pi[i] * Q[i, j])
      P[i, j] <- Q[i, j] * min(1, ratio)
    }
  }
  diag(P) <- 1 - rowSums(P)
  P
}
