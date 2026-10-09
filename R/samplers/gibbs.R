# Gibbs samplers (Chapter 13) ------------------------------------------------------
#
# A Gibbs sampler updates blocks of parameters by drawing each block from its
# full conditional distribution given the current values of all other blocks.
# Each such update leaves the joint posterior invariant (Chapter 13), so the
# composition of the updates does too.

#' Systematic-scan Gibbs sampler for a bivariate normal with unit variances and
#' correlation rho. Full conditionals: x1 | x2 ~ N(rho x2, 1 - rho^2) and
#' symmetrically. The x1-marginal of the chain is an AR(1) process with
#' coefficient rho^2, which makes the cost of strong correlation explicit.
gibbs_bivariate_normal <- function(rho, n_iter, initial = c(0, 0)) {
  if (!is.numeric(rho) || length(rho) != 1L || abs(rho) >= 1) course_stop("`rho` must lie in (-1, 1).")
  n_iter <- check_count(n_iter, "n_iter")
  check_finite_vector(initial, "initial", len = 2L)
  s <- sqrt(1 - rho^2)
  out <- matrix(NA_real_, n_iter, 2, dimnames = list(NULL, c("x1", "x2")))
  x <- initial
  for (i in seq_len(n_iter)) {
    x[1] <- stats::rnorm(1, rho * x[2], s)
    x[2] <- stats::rnorm(1, rho * x[1], s)
    out[i, ] <- x
  }
  out
}

#' Gibbs sampler for the semi-conjugate normal model
#'   y_i | mu, s2 ~ N(mu, s2),  mu ~ N(m0, v0),  s2 ~ InvGamma(a0, b0)
#' (independent priors, so the joint posterior is NOT of closed form).
#' Full conditionals:
#'   mu | s2, y ~ N( (m0/v0 + n ybar/s2) / (1/v0 + n/s2), 1 / (1/v0 + n/s2) )
#'   s2 | mu, y ~ InvGamma( a0 + n/2, b0 + sum (y_i - mu)^2 / 2 )
gibbs_normal <- function(y, n_iter, m0 = 0, v0 = 100, a0 = 1, b0 = 1,
                         initial = c(mu = mean(y), sigma2 = stats::var(y))) {
  check_finite_vector(y, "y")
  n_iter <- check_count(n_iter, "n_iter")
  check_scalar(v0, "v0", positive = TRUE); check_scalar(a0, "a0", positive = TRUE)
  check_scalar(b0, "b0", positive = TRUE)
  n <- length(y); ybar <- mean(y)
  mu <- initial[[1]]; s2 <- initial[[2]]
  if (!(s2 > 0)) course_stop("Initial variance must be positive.")
  out <- matrix(NA_real_, n_iter, 2, dimnames = list(NULL, c("mu", "sigma2")))
  for (i in seq_len(n_iter)) {
    prec <- 1 / v0 + n / s2
    mu <- stats::rnorm(1, (m0 / v0 + n * ybar / s2) / prec, sqrt(1 / prec))
    s2 <- rinvgamma(1, a0 + n / 2, b0 + 0.5 * sum((y - mu)^2))
    out[i, ] <- c(mu, s2)
  }
  out
}

#' Gibbs sampler for the semi-conjugate linear regression
#'   y | beta, s2 ~ N(X beta, s2 I),  beta ~ N(m0, S0),  s2 ~ InvGamma(a0, b0).
#' Full conditionals (Chapter 20):
#'   beta | s2, y ~ N(Vn (S0^{-1} m0 + X'y / s2), Vn),  Vn = (S0^{-1} + X'X / s2)^{-1}
#'   s2 | beta, y ~ InvGamma(a0 + n/2, b0 + ||y - X beta||^2 / 2)
gibbs_linear_regression <- function(X, y, n_iter, m0 = rep(0, ncol(X)),
                                    S0 = diag(100, ncol(X)), a0 = 1, b0 = 1) {
  X <- as.matrix(X); n <- nrow(X); p <- ncol(X)
  n_iter <- check_count(n_iter, "n_iter")
  check_spd_matrix(S0, "S0", d = p)
  S0_inv <- solve(S0); XtX <- crossprod(X); Xty <- as.vector(crossprod(X, y))
  s2 <- stats::var(y)
  cn <- colnames(X); if (is.null(cn)) cn <- paste0("beta[", seq_len(p) - 1L, "]")
  out <- matrix(NA_real_, n_iter, p + 1L, dimnames = list(NULL, c(cn, "sigma2")))
  for (i in seq_len(n_iter)) {
    Vn_inv <- S0_inv + XtX / s2
    U <- chol(Vn_inv)                                   # Vn_inv = U'U
    mn <- backsolve(U, backsolve(U, S0_inv %*% m0 + Xty / s2, transpose = TRUE))
    beta <- as.vector(mn + backsolve(U, stats::rnorm(p)))   # cov = U^{-1} U^{-T} = Vn
    resid <- y - as.vector(X %*% beta)
    s2 <- rinvgamma(1, a0 + n / 2, b0 + 0.5 * sum(resid^2))
    out[i, ] <- c(beta, s2)
  }
  out
}

#' Centred Gibbs sampler for the normal hierarchical model (eight schools).
#'
#' theta_j and mu have normal full conditionals. With a half-normal prior on
#' tau the full conditional of tau is not a standard distribution:
#'   p(tau | theta, mu) propto tau^{-J} exp(-S / (2 tau^2) - tau^2 / (2 s_tau^2)),
#'   S = sum_j (theta_j - mu)^2,
#' so tau is updated by a random-walk Metropolis step on log(tau)
#' ("Metropolis-within-Gibbs"). Each step leaves the posterior invariant.
gibbs_hier_centered <- function(y, sigma, n_iter, s_mu = 5, s_tau = 5,
                                initial = list(mu = 0, tau = 1), log_tau_step = 0.5) {
  J <- length(y); n_iter <- check_count(n_iter, "n_iter")
  mu <- initial$mu; tau <- initial$tau; theta <- y
  out <- matrix(NA_real_, n_iter, J + 2L,
                dimnames = list(NULL, c("mu", "tau", paste0("theta[", seq_len(J), "]"))))
  acc <- 0L
  log_cond_logtau <- function(u, S) -J * u - S * exp(-2 * u) / 2 - exp(2 * u) / (2 * s_tau^2) + u
  for (i in seq_len(n_iter)) {
    # theta_j | mu, tau, y: precision-weighted combination of y_j and mu
    prec <- 1 / sigma^2 + 1 / tau^2
    theta <- stats::rnorm(J, (y / sigma^2 + mu / tau^2) / prec, sqrt(1 / prec))
    # mu | theta, tau
    prec_mu <- J / tau^2 + 1 / s_mu^2
    mu <- stats::rnorm(1, (sum(theta) / tau^2) / prec_mu, sqrt(1 / prec_mu))
    # log(tau) | theta, mu: one Metropolis step (includes the log-Jacobian + u)
    S <- sum((theta - mu)^2)
    u <- log(tau); u_new <- u + log_tau_step * stats::rnorm(1)
    if (log(stats::runif(1)) < log_cond_logtau(u_new, S) - log_cond_logtau(u, S)) {
      tau <- exp(u_new); acc <- acc + 1L
    }
    out[i, ] <- c(mu, tau, theta)
  }
  attr(out, "tau_acceptance") <- acc / n_iter
  out
}

#' Non-centred Gibbs sampler for the same model, using parameter expansion.
#'
#' Write theta_j = mu + tau * eta_j with eta_j ~ N(0, 1) and let tau range over
#' the whole real line with prior N(0, s_tau^2); then |tau| has exactly the
#' half-normal(0, s_tau) prior and theta_j = mu + tau eta_j has the same joint
#' law as in the centred model (the sign of (tau, eta) is not identified, but
#' theta, mu and |tau| are). Every full conditional is then normal:
#'   eta_j | . ~ N( tau (y_j - mu) / sigma_j^2 / P_j, 1 / P_j ),  P_j = tau^2 / sigma_j^2 + 1
#'   mu    | . ~ N( sum (y_j - tau eta_j) / sigma_j^2 / P_mu, 1 / P_mu ),  P_mu = sum 1/sigma_j^2 + 1/s_mu^2
#'   tau   | . ~ N( sum eta_j (y_j - mu) / sigma_j^2 / P_tau, 1 / P_tau ),  P_tau = sum eta_j^2/sigma_j^2 + 1/s_tau^2
gibbs_hier_noncentered <- function(y, sigma, n_iter, s_mu = 5, s_tau = 5,
                                   initial = list(mu = 0, tau = 1)) {
  J <- length(y); n_iter <- check_count(n_iter, "n_iter")
  mu <- initial$mu; tau <- initial$tau; eta <- rep(0, J)
  out <- matrix(NA_real_, n_iter, J + 2L,
                dimnames = list(NULL, c("mu", "tau", paste0("theta[", seq_len(J), "]"))))
  w <- 1 / sigma^2
  for (i in seq_len(n_iter)) {
    P <- tau^2 * w + 1
    eta <- stats::rnorm(J, tau * (y - mu) * w / P, sqrt(1 / P))
    P_mu <- sum(w) + 1 / s_mu^2
    mu <- stats::rnorm(1, sum((y - tau * eta) * w) / P_mu, sqrt(1 / P_mu))
    P_tau <- sum(eta^2 * w) + 1 / s_tau^2
    tau <- stats::rnorm(1, sum(eta * (y - mu) * w) / P_tau, sqrt(1 / P_tau))
    out[i, ] <- c(mu, abs(tau), mu + tau * eta)
  }
  out
}
