# Tempering and sequential Monte Carlo (Chapters 14 and 25) -------------------------

#' Parallel tempering (replica exchange; Swendsen and Wang, 1986; Geyer, 1991).
#'
#' Runs K random-walk Metropolis chains targeting pi_k(x) propto pi(x)^{beta_k}
#' with 1 = beta_1 > beta_2 > ... > beta_K > 0. After each sweep, a swap between
#' a randomly chosen adjacent pair (k, k+1) is proposed and accepted with
#' probability min{1, [pi(x_{k+1}) / pi(x_k)]^(beta_k - beta_{k+1})}, which keeps
#' the product distribution prod_k pi_k invariant. Only the beta = 1 chain
#' targets pi; the hotter chains exist to move between modes.
parallel_tempering <- function(log_target, initial, n_iter, betas, proposal_sd) {
  check_function(log_target, "log_target"); n_iter <- check_count(n_iter, "n_iter")
  K <- length(betas)
  if (K < 2L || betas[1] != 1 || any(diff(betas) >= 0) || any(betas <= 0)) {
    course_stop("`betas` must start at 1 and decrease strictly towards 0.")
  }
  proposal_sd <- rep_len(proposal_sd, K)
  x <- rep_len(initial, K)
  lp <- vapply(x, function(xi) eval_log_density(log_target, xi), numeric(1))
  if (any(!is.finite(lp))) course_stop("Initial values must have finite log target density.")
  cold <- numeric(n_iter); all_states <- matrix(NA_real_, n_iter, K)
  acc_move <- numeric(K); swap_try <- numeric(K - 1L); swap_acc <- numeric(K - 1L)
  for (t in seq_len(n_iter)) {
    for (k in seq_len(K)) {                         # within-temperature Metropolis moves
      y <- x[k] + proposal_sd[k] * stats::rnorm(1)
      lp_y <- eval_log_density(log_target, y)
      if (log(stats::runif(1)) < betas[k] * (lp_y - lp[k])) {
        x[k] <- y; lp[k] <- lp_y; acc_move[k] <- acc_move[k] + 1
      }
    }
    k <- sample.int(K - 1L, 1L)                     # propose swapping states of k and k+1
    swap_try[k] <- swap_try[k] + 1
    log_r <- (betas[k] - betas[k + 1L]) * (lp[k + 1L] - lp[k])
    if (log(stats::runif(1)) < log_r) {
      x[c(k, k + 1L)] <- x[c(k + 1L, k)]; lp[c(k, k + 1L)] <- lp[c(k + 1L, k)]
      swap_acc[k] <- swap_acc[k] + 1
    }
    cold[t] <- x[1]; all_states[t, ] <- x
  }
  list(draws = cold, states = all_states, move_acceptance = acc_move / n_iter,
       swap_acceptance = ifelse(swap_try > 0, swap_acc / swap_try, NA_real_))
}

#' Importance sampling with an un-normalised target.
#'
#' Draws x_1..x_n ~ g and forms log weights log w_i = log pi~(x_i) - log g(x_i).
#' Returns the self-normalised estimate of E_pi[f], the estimate of log Z
#' (log normalising constant of pi~), Kish's effective sample size
#' (sum w)^2 / sum w^2, and a delta-method standard error for the
#' self-normalised estimator.
importance_sampling <- function(f, log_target, r_proposal, log_proposal, n) {
  n <- check_count(n, "n")
  x <- r_proposal(n)
  xs <- if (is.matrix(x)) split(x, row(x)) else as.list(x)
  log_w <- vapply(xs, function(xi) log_target(xi), numeric(1)) - log_proposal(x)
  w <- normalise_log_weights(log_w)
  fx <- vapply(xs, f, numeric(1))
  est <- sum(w * fx)
  list(estimate = est,
       log_Z = log_mean_exp(log_w),
       ess = 1 / sum(w^2),
       se = sqrt(sum(w^2 * (fx - est)^2)),
       weights = w, log_weights = log_w, x = x)
}

#' Kish's effective sample size of a vector of (unnormalised) weights.
kish_ess <- function(w) sum(w)^2 / sum(w^2)

#' Adaptive-tempering sequential Monte Carlo sampler for a posterior
#'   pi(theta) propto L(theta) p0(theta)
#' through the bridge pi_b propto L^b p0, 0 = b_0 < b_1 < ... < b_T = 1.
#' Each new temperature is chosen by bisection so that the conditional ESS of
#' the incremental weights is `ess_frac` times n (an adaptive schedule; see
#' Chopin and Papaspiliopoulos, 2020); particles are resampled (multinomial) and moved with `n_mh`
#' random-walk Metropolis steps that leave pi_b invariant, with proposal
#' scale set from the current particle spread. Also returns an estimate of
#' log of the marginal likelihood  log integral L p0.
smc_tempering <- function(log_lik, log_prior, r_prior, n = 1000, ess_frac = 0.5,
                          n_mh = 5, scale = 2.38) {
  n <- check_count(n, "n")
  th <- r_prior(n); th <- if (is.matrix(th)) th else matrix(th, ncol = 1)
  d <- ncol(th)
  rows <- function(m) split(m, row(m))
  ll <- vapply(rows(th), log_lik, numeric(1))
  b <- 0; log_Z <- 0; schedule <- 0
  while (b < 1) {
    ess_at <- function(b_new) { lw <- (b_new - b) * ll; w <- normalise_log_weights(lw); 1 / sum(w^2) }
    b_new <- if (ess_at(1) >= ess_frac * n) 1 else
      stats::uniroot(function(z) ess_at(z) - ess_frac * n, c(b, 1), tol = 1e-10)$root
    lw <- (b_new - b) * ll
    log_Z <- log_Z + log_mean_exp(lw)              # incremental normalising-constant ratio
    w <- normalise_log_weights(lw)
    idx <- sample.int(n, n, replace = TRUE, prob = w)
    th <- th[idx, , drop = FALSE]; ll <- ll[idx]
    b <- b_new; schedule <- c(schedule, b)
    # MH moves targeting pi_b, proposal covariance from the particle cloud
    Sig <- (scale^2 / d) * (stats::cov(th) + 1e-10 * diag(d))
    U <- chol(Sig)
    lp <- vapply(rows(th), log_prior, numeric(1)) + b * ll
    for (s in seq_len(n_mh)) {
      prop <- th + matrix(stats::rnorm(n * d), n, d) %*% U
      ll_p <- vapply(rows(prop), log_lik, numeric(1))
      lp_p <- vapply(rows(prop), log_prior, numeric(1)) + b * ll_p
      lp_p[is.nan(lp_p)] <- -Inf
      acc <- log(stats::runif(n)) < lp_p - lp
      acc[is.na(acc)] <- FALSE
      th[acc, ] <- prop[acc, ]; ll[acc] <- ll_p[acc]; lp[acc] <- lp_p[acc]
    }
  }
  list(draws = th, log_marginal_likelihood = log_Z, temperatures = schedule)
}
