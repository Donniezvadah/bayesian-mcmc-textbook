# Bayesian regression models -----------------------------------------------------
#
# Log posteriors and gradients for generalised linear models with independent
# Gaussian priors on the coefficients, plus the conjugate Normal-Inverse-Gamma
# linear model. Notation: X is the n x p design matrix (first column usually 1),
# beta the coefficient vector, eta = X beta the linear predictor.

#' Log posterior (up to a constant) and gradient for logistic regression:
#'   y_i | beta ~ Bernoulli(inv_logit(eta_i)),  beta_j ~ N(prior_mean_j, prior_sd_j^2).
#' Uses log p(y_i = 1) = log_inv_logit(eta_i) and log p(y_i = 0) = log1m_inv_logit(eta_i),
#' which do not overflow for large |eta|.
logistic_regression_model <- function(X, y, prior_sd = 2.5, prior_mean = 0) {
  X <- as.matrix(X)
  if (!is.null(dim(y))) course_stop("`y` must be a vector, not a matrix.")
  if (!all(y %in% c(0, 1))) course_stop("`y` must be 0/1.")
  if (nrow(X) != length(y)) course_stop("nrow(X) must equal length(y).")
  p <- ncol(X)
  prior_sd <- rep_len(prior_sd, p); prior_mean <- rep_len(prior_mean, p)
  list(
    dim = p,
    log_density = function(beta) {
      eta <- as.vector(X %*% beta)
      sum(y * log_inv_logit(eta) + (1 - y) * log1m_inv_logit(eta)) +
        sum(stats::dnorm(beta, prior_mean, prior_sd, log = TRUE))
    },
    grad = function(beta) {
      eta <- as.vector(X %*% beta)
      as.vector(crossprod(X, y - inv_logit(eta))) - (beta - prior_mean) / prior_sd^2
    },
    log_lik_pointwise = function(beta) {
      eta <- as.vector(X %*% beta)
      y * log_inv_logit(eta) + (1 - y) * log1m_inv_logit(eta)
    }
  )
}

#' Log posterior and gradient for Poisson regression with log link:
#'   y_i | beta ~ Poisson(exp(eta_i + offset_i)),  beta_j ~ N(prior_mean_j, prior_sd_j^2).
poisson_regression_model <- function(X, y, prior_sd = 2.5, prior_mean = 0,
                                     offset = rep(0, length(y))) {
  X <- as.matrix(X)
  if (!is.null(dim(y))) course_stop("`y` must be a vector, not a matrix.")
  if (any(y < 0 | y != round(y))) course_stop("`y` must contain non-negative integers.")
  if (nrow(X) != length(y)) course_stop("nrow(X) must equal length(y).")
  p <- ncol(X)
  prior_sd <- rep_len(prior_sd, p); prior_mean <- rep_len(prior_mean, p)
  list(
    dim = p,
    log_density = function(beta) {
      eta <- as.vector(X %*% beta) + offset
      sum(y * eta - exp(eta) - lgamma(y + 1)) +
        sum(stats::dnorm(beta, prior_mean, prior_sd, log = TRUE))
    },
    grad = function(beta) {
      eta <- as.vector(X %*% beta) + offset
      as.vector(crossprod(X, y - exp(eta))) - (beta - prior_mean) / prior_sd^2
    },
    log_lik_pointwise = function(beta) {
      eta <- as.vector(X %*% beta) + offset
      y * eta - exp(eta) - lgamma(y + 1)
    }
  )
}

#' Conjugate Bayesian linear regression:
#'   y | beta, s2 ~ N(X beta, s2 I),  beta | s2 ~ N(m0, s2 V0),  s2 ~ InvGamma(a0, b0).
#' Posterior: beta | s2, y ~ N(mn, s2 Vn), s2 | y ~ InvGamma(an, bn) with
#'   Vn = (V0^{-1} + X'X)^{-1},  mn = Vn (V0^{-1} m0 + X'y),
#'   an = a0 + n/2,  bn = b0 + (y'y + m0' V0^{-1} m0 - mn' Vn^{-1} mn) / 2.
#' The marginal posterior of beta is multivariate t with 2 an degrees of freedom,
#' location mn and scale matrix (bn / an) Vn.
conjugate_linear_regression <- function(X, y, m0 = rep(0, ncol(X)), V0 = diag(1e2, ncol(X)),
                                        a0 = 1, b0 = 1) {
  X <- as.matrix(X); n <- nrow(X); p <- ncol(X)
  check_spd_matrix(V0, "V0", d = p)
  V0_inv <- solve(V0)
  Vn_inv <- V0_inv + crossprod(X)
  Vn <- solve(Vn_inv)
  mn <- as.vector(Vn %*% (V0_inv %*% m0 + crossprod(X, y)))
  an <- a0 + n / 2
  bn <- b0 + 0.5 * (sum(y^2) + sum(m0 * (V0_inv %*% m0)) - sum(mn * (Vn_inv %*% mn)))
  list(m = mn, V = Vn, a = an, b = bn,
       draw = function(n_draws) {
         s2 <- rinvgamma(n_draws, an, bn)
         U <- chol(Vn)
         z <- matrix(stats::rnorm(n_draws * p), n_draws, p)
         beta <- sweep((z %*% U) * sqrt(s2), 2, mn, "+")
         cn <- colnames(X); if (is.null(cn)) cn <- paste0("beta[", seq_len(p) - 1L, "]")
         colnames(beta) <- cn
         cbind(beta, sigma2 = s2)
       })
}
