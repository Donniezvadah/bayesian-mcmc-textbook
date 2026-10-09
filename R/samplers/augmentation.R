# Latent-variable augmentation (Chapter 25) --------------------------------------------

#' Draws from N(mean, 1) truncated to (0, Inf) if positive = TRUE, else to (-Inf, 0).
#' Uses the inverse CDF on the log scale, which is accurate even far in the tails:
#' for W = Z - mean truncated to (-mean, Inf),  Phi(-W) = V Phi(mean), V ~ U(0,1).
rtruncnorm_sign <- function(mean, positive) {
  n <- length(mean)
  s <- ifelse(positive, 1, -1)                 # reflect negative truncation to positive
  m <- s * mean
  log_p <- log(stats::runif(n)) + stats::pnorm(m, log.p = TRUE)
  w <- -stats::qnorm(log_p, log.p = TRUE)       # W > -m
  s * (m + w)
}

#' Albert and Chib (1993) data-augmentation Gibbs sampler for probit regression:
#'   y_i = 1{z_i > 0},  z_i ~ N(x_i' beta, 1),  beta ~ N(0, prior_sd^2 I).
#' Full conditionals:
#'   z_i | beta, y ~ N(x_i' beta, 1) truncated to (0, Inf) if y_i = 1, (-Inf, 0) if y_i = 0;
#'   beta | z ~ N(V X'z, V),  V = (X'X + I / prior_sd^2)^{-1}.
albert_chib_probit <- function(X, y, n_iter, prior_sd = 10, initial = rep(0, ncol(X))) {
  X <- as.matrix(X); n_iter <- check_count(n_iter, "n_iter")
  if (!all(y %in% c(0, 1))) course_stop("`y` must be 0/1.")
  p <- ncol(X)
  U <- chol(crossprod(X) + diag(1 / prior_sd^2, p))          # V^{-1} = U'U
  beta <- initial
  out <- matrix(NA_real_, n_iter, p, dimnames = list(NULL, colnames(X)))
  for (t in seq_len(n_iter)) {
    z <- rtruncnorm_sign(as.vector(X %*% beta), positive = y == 1)
    mean_beta <- backsolve(U, backsolve(U, crossprod(X, z), transpose = TRUE))
    beta <- as.vector(mean_beta + backsolve(U, stats::rnorm(p)))
    out[t, ] <- beta
  }
  out
}
