# Finite-state Markov chains (Chapters 10-11) --------------------------------------

#' Check that P is a transition (row-stochastic) matrix.
check_transition_matrix <- function(P, name = "P") {
  if (!is.matrix(P) || nrow(P) != ncol(P) || any(!is.finite(P)) || any(P < 0) ||
      any(abs(rowSums(P) - 1) > 1e-10)) {
    course_stop("`%s` must be a square matrix with non-negative entries and rows summing to 1.", name)
  }
  invisible(P)
}

#' Stationary distribution of a finite transition matrix by solving
#' pi (P - I) = 0 together with sum(pi) = 1 (a least-squares solve of the
#' stacked system, which is exact when the solution is unique).
stationary_distribution <- function(P) {
  check_transition_matrix(P)
  K <- nrow(P)
  A <- rbind(t(P) - diag(K), rep(1, K))
  b <- c(rep(0, K), 1)
  pi <- as.vector(qr.solve(A, b))
  pi[abs(pi) < 1e-15] <- 0
  pi
}

#' n-step transition matrix P^n by repeated squaring (O(K^3 log n)).
matrix_power <- function(P, n) {
  n <- check_count(n, "n", min = 0L)
  result <- diag(nrow(P)); base <- P
  while (n > 0) {
    if (n %% 2 == 1) result <- result %*% base
    base <- base %*% base; n <- n %/% 2
  }
  result
}

#' Total variation distance between two probability vectors: (1/2) sum |p - q|.
tv_distance <- function(p, q) 0.5 * sum(abs(p - q))

#' Detailed balance check: pi_i P_ij == pi_j P_ji for all i, j.
is_reversible <- function(P, pi, tol = 1e-10) {
  flux <- pi * P
  max(abs(flux - t(flux))) < tol
}

#' Simulate a path X_0, ..., X_n of a finite Markov chain with states 1..K.
simulate_chain <- function(P, n, x0) {
  check_transition_matrix(P); n <- check_count(n, "n")
  K <- nrow(P); x <- integer(n + 1L); x[1] <- x0
  # inverse-CDF sampling from each row, using precomputed cumulative rows
  C <- t(apply(P, 1, cumsum))
  u <- stats::runif(n)
  for (t in seq_len(n)) x[t + 1L] <- findInterval(u[t], C[x[t], ], left.open = TRUE) + 1L
  pmin(x, K)
}

#' Simulate a stationary Gaussian AR(1) chain X_t = phi X_{t-1} + sigma e_t.
#' Its stationary law is N(0, sigma^2 / (1 - phi^2)) for |phi| < 1.
simulate_ar1 <- function(n, phi, sigma = 1, x0 = NULL) {
  n <- check_count(n, "n")
  if (abs(phi) >= 1) course_stop("|phi| must be < 1 for a stationary AR(1) chain.")
  x <- numeric(n)
  x[1] <- if (is.null(x0)) stats::rnorm(1, 0, sigma / sqrt(1 - phi^2)) else x0
  if (n > 1) for (t in 2:n) x[t] <- phi * x[t - 1] + sigma * stats::rnorm(1)
  x
}
