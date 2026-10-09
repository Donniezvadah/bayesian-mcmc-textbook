# Normal hierarchical model ("eight schools") ----------------------------------
#
# Data: estimated treatment effects y_j with known standard errors sigma_j for
# J = 8 schools (Rubin, 1981; Gelman et al., 2013, Section 5.5).
# Model (used in Chapters 13, 16, 21 and Project 4):
#   y_j | theta_j       ~ N(theta_j, sigma_j^2)       (sigma_j known)
#   theta_j | mu, tau   ~ N(mu, tau^2)
#   mu                  ~ N(0, s_mu^2)
#   tau                 ~ half-normal(0, s_tau)       (density 2 N(tau | 0, s_tau^2), tau > 0)

eight_schools <- data.frame(
  school = LETTERS[1:8],
  y      = c(28, 8, -3, 7, -1, 1, 18, 12),
  sigma  = c(15, 10, 16, 11, 9, 11, 10, 18)
)

#' Centred parameterisation on an unconstrained space.
#' Parameter vector: (mu, log_tau, theta_1, ..., theta_J). The log density
#' includes the log-Jacobian log(tau) of the transformation tau = exp(log_tau).
eight_schools_centered <- function(y = eight_schools$y, sigma = eight_schools$sigma,
                                   s_mu = 5, s_tau = 5) {
  J <- length(y)
  list(
    dim = J + 2L,
    names = c("mu", "log_tau", paste0("theta[", seq_len(J), "]")),
    log_density = function(z) {
      mu <- z[1]; lt <- z[2]; tau <- exp(lt); th <- z[-(1:2)]
      sum(stats::dnorm(y, th, sigma, log = TRUE)) +
        sum(stats::dnorm(th, mu, tau, log = TRUE)) +
        stats::dnorm(mu, 0, s_mu, log = TRUE) +
        log(2) + stats::dnorm(tau, 0, s_tau, log = TRUE) + lt
    },
    grad = function(z) {
      mu <- z[1]; lt <- z[2]; tau <- exp(lt); th <- z[-(1:2)]
      r <- th - mu
      g_mu <- sum(r) / tau^2 - mu / s_mu^2
      # d/d(log tau) of [-J log tau - sum r^2 / (2 tau^2) - tau^2 / (2 s_tau^2) + log tau]
      g_lt <- -J + sum(r^2) / tau^2 - tau^2 / s_tau^2 + 1
      g_th <- (y - th) / sigma^2 - r / tau^2
      c(g_mu, g_lt, g_th)
    },
    to_natural = function(z) c(mu = z[1], tau = exp(z[2]), stats::setNames(z[-(1:2)], paste0("theta[", seq_len(J), "]")))
  )
}

#' Non-centred parameterisation: theta_j = mu + tau * eta_j with eta_j ~ N(0, 1).
#' Parameter vector: (mu, log_tau, eta_1, ..., eta_J).
eight_schools_noncentered <- function(y = eight_schools$y, sigma = eight_schools$sigma,
                                      s_mu = 5, s_tau = 5) {
  J <- length(y)
  list(
    dim = J + 2L,
    names = c("mu", "log_tau", paste0("eta[", seq_len(J), "]")),
    log_density = function(z) {
      mu <- z[1]; lt <- z[2]; tau <- exp(lt); eta <- z[-(1:2)]
      th <- mu + tau * eta
      sum(stats::dnorm(y, th, sigma, log = TRUE)) +
        sum(stats::dnorm(eta, 0, 1, log = TRUE)) +
        stats::dnorm(mu, 0, s_mu, log = TRUE) +
        log(2) + stats::dnorm(tau, 0, s_tau, log = TRUE) + lt
    },
    grad = function(z) {
      mu <- z[1]; lt <- z[2]; tau <- exp(lt); eta <- z[-(1:2)]
      e <- (y - mu - tau * eta) / sigma^2          # d loglik / d theta_j
      g_mu <- sum(e) - mu / s_mu^2
      g_lt <- tau * sum(e * eta) - tau^2 / s_tau^2 + 1
      g_eta <- tau * e - eta
      c(g_mu, g_lt, g_eta)
    },
    to_natural = function(z) {
      tau <- exp(z[2])
      c(mu = z[1], tau = tau, stats::setNames(z[1] + tau * z[-(1:2)], paste0("theta[", seq_len(J), "]")))
    }
  )
}

#' Exact independent draws from the eight-schools posterior (normal prior on mu,
#' half-normal prior on tau), by the factorisation of Gelman et al. (2013, Section 5.4):
#'   p(tau | y) on a fine grid  ->  mu | tau, y ~ Normal  ->  theta | mu, tau, y ~ Normal.
#' With mu ~ N(0, s_mu^2) and v_j = sigma_j^2 + tau^2, integrating mu out gives
#'   p(y | tau) propto V^{1/2} prod_j v_j^{-1/2} exp{ -1/2 sum_j y_j^2 / v_j + mu_hat^2 / (2 V) },
#'   V = 1 / (sum_j 1/v_j + 1/s_mu^2),  mu_hat = V sum_j y_j / v_j.
#' The only approximation is the discretisation of tau, which is negligible on a fine grid.
eight_schools_exact <- function(n_draws, y = eight_schools$y, sigma = eight_schools$sigma,
                                s_mu = 5, s_tau = 5, tau_max = 10 * s_tau, n_grid = 20000) {
  n_draws <- check_count(n_draws, "n_draws")
  J <- length(y)
  h <- tau_max / n_grid
  tau_grid <- seq(h / 2, tau_max - h / 2, length.out = n_grid)
  log_post_tau <- vapply(tau_grid, function(tau) {
    v <- sigma^2 + tau^2; V <- 1 / (sum(1 / v) + 1 / s_mu^2); mu_hat <- V * sum(y / v)
    0.5 * log(V) - 0.5 * sum(log(v)) - 0.5 * sum(y^2 / v) + 0.5 * mu_hat^2 / V +
      stats::dnorm(tau, 0, s_tau, log = TRUE)
  }, numeric(1))
  w <- normalise_log_weights(log_post_tau)
  tau <- tau_grid[sample.int(n_grid, n_draws, replace = TRUE, prob = w)] + stats::runif(n_draws, -h / 2, h / 2)
  tau <- abs(tau)
  out <- matrix(NA_real_, n_draws, J + 2L, dimnames = list(NULL, c("mu", "tau", paste0("theta[", seq_len(J), "]"))))
  for (s in seq_len(n_draws)) {
    v <- sigma^2 + tau[s]^2; V <- 1 / (sum(1 / v) + 1 / s_mu^2)
    mu <- stats::rnorm(1, V * sum(y / v), sqrt(V))
    prec <- 1 / sigma^2 + 1 / tau[s]^2
    out[s, ] <- c(mu, tau[s], stats::rnorm(J, (y / sigma^2 + mu / tau[s]^2) / prec, sqrt(1 / prec)))
  }
  out
}

#' Poisson regression with an observation-level random effect (a Poisson-lognormal
#' model for overdispersed counts), in the non-centred parameterisation:
#'   y_i ~ Poisson(exp(x_i' beta + sigma_u z_i)),  z_i ~ N(0, 1),
#'   beta_j ~ N(0, prior_sd^2),  sigma_u ~ half-normal(0, s_u).
#' Parameter vector: (beta, log sigma_u, z); log density includes the Jacobian.
poisson_lognormal_model <- function(X, y, prior_sd = 2.5, s_u = 1) {
  X <- as.matrix(X); n <- nrow(X); p <- ncol(X)
  prior_sd <- rep_len(prior_sd, p)
  idx_b <- seq_len(p); idx_ls <- p + 1L; idx_z <- p + 1L + seq_len(n)
  list(
    dim = p + 1L + n,
    log_density = function(par) {
      b <- par[idx_b]; su <- exp(par[idx_ls]); z <- par[idx_z]
      eta <- as.vector(X %*% b) + su * z
      sum(y * eta - exp(eta) - lgamma(y + 1)) + sum(stats::dnorm(b, 0, prior_sd, log = TRUE)) +
        log(2) + stats::dnorm(su, 0, s_u, log = TRUE) + par[idx_ls] + sum(stats::dnorm(z, log = TRUE))
    },
    grad = function(par) {
      b <- par[idx_b]; su <- exp(par[idx_ls]); z <- par[idx_z]
      eta <- as.vector(X %*% b) + su * z; r <- y - exp(eta)
      c(as.vector(crossprod(X, r)) - b / prior_sd^2,
        su * sum(z * r) - su^2 / s_u^2 + 1,
        su * r - z)
    }
  )
}
