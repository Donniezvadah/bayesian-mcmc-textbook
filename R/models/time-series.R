# Time-series models (Chapter 22) ------------------------------------------------------

#' Design matrix for an AR(p) model with intercept, conditioning on the first p
#' observations: y_t = c + phi_1 y_{t-1} + ... + phi_p y_{t-p} + e_t, t = p+1..n.
ar_design <- function(y, p) {
  check_finite_vector(y, "y"); p <- check_count(p, "p")
  n <- length(y)
  if (n <= p + 1) course_stop("Series too short for AR(%d).", p)
  X <- cbind(1, sapply(seq_len(p), function(k) y[(p + 1 - k):(n - k)]))
  colnames(X) <- c("c", paste0("phi", seq_len(p)))
  list(X = X, y = y[(p + 1):n])
}

#' Is an AR(p) coefficient vector stationary? (all roots of 1 - phi_1 z - ... - phi_p z^p
#' lie outside the unit circle)
ar_is_stationary <- function(phi) all(Mod(polyroot(c(1, -phi))) > 1)

#' Kalman filter for the local level model
#'   y_t = mu_t + eps_t,   eps_t ~ N(0, s2_eps)
#'   mu_t = mu_{t-1} + eta_t, eta_t ~ N(0, s2_eta),   mu_1 ~ N(m0, C0).
#' Returns filtered means m_t, variances C_t, one-step predictive means f_t and
#' variances Q_t, and the log-likelihood log p(y_1..y_n | s2_eps, s2_eta) by the
#' prediction-error decomposition  sum_t log N(y_t | f_t, Q_t).
kalman_local_level <- function(y, s2_eps, s2_eta, m0 = 0, C0 = 1e7) {
  n <- length(y)
  m <- C <- f <- Q <- numeric(n)
  a <- m0; R <- C0                                  # prior for mu_1
  loglik <- 0
  for (t in seq_len(n)) {
    if (t > 1) { a <- m[t - 1]; R <- C[t - 1] + s2_eta }   # predict state
    f[t] <- a; Q[t] <- R + s2_eps                          # predict observation
    K <- R / Q[t]                                          # Kalman gain
    m[t] <- a + K * (y[t] - f[t]); C[t] <- (1 - K) * R      # update
    loglik <- loglik + stats::dnorm(y[t], f[t], sqrt(Q[t]), log = TRUE)
  }
  list(m = m, C = C, f = f, Q = Q, loglik = loglik)
}

#' Forward-filtering backward-sampling (Carter and Kohn, 1994; Fruhwirth-Schnatter, 1994):
#' one joint draw of mu_1..mu_n from p(mu | y, s2_eps, s2_eta). Backward step:
#'   mu_n ~ N(m_n, C_n);  mu_t | mu_{t+1}, y_{1:t} ~ N(m_t + B_t (mu_{t+1} - m_t), C_t (1 - B_t)),
#'   B_t = C_t / (C_t + s2_eta).
ffbs_local_level <- function(y, s2_eps, s2_eta, m0 = 0, C0 = 1e7, kf = NULL) {
  if (is.null(kf)) kf <- kalman_local_level(y, s2_eps, s2_eta, m0, C0)
  n <- length(y); mu <- numeric(n)
  mu[n] <- stats::rnorm(1, kf$m[n], sqrt(kf$C[n]))
  for (t in (n - 1):1) {
    B <- kf$C[t] / (kf$C[t] + s2_eta)
    mu[t] <- stats::rnorm(1, kf$m[t] + B * (mu[t + 1] - kf$m[t]), sqrt(kf$C[t] * (1 - B)))
  }
  mu
}
