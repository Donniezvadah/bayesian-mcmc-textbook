test_that("Kalman log-likelihood equals the dense multivariate normal log density", {
  set.seed(140)
  y <- cumsum(rnorm(30, 0, 0.5)) + rnorm(30, 0, 1)
  s2e <- 1.2; s2n <- 0.3; m0 <- 0.5; C0 <- 4
  n <- length(y); Tm <- outer(1:n, 1:n, pmin) - 1                  # Cov of random-walk increments
  Sigma <- C0 + s2n * Tm + diag(s2e, n)
  kf <- kalman_local_level(y, s2e, s2n, m0, C0)
  expect_equal(kf$loglik, dmvnorm_log(y, rep(m0, n), Sigma), tolerance = 1e-10)
})

test_that("FFBS draws have the exact smoothing mean and variance", {
  set.seed(141)
  y <- cumsum(rnorm(15, 0, 0.5)) + rnorm(15, 0, 1)
  s2e <- 1; s2n <- 0.25; m0 <- 0; C0 <- 10; n <- 15
  Tm <- outer(1:n, 1:n, pmin) - 1
  Smu <- C0 + s2n * Tm                                              # prior covariance of states
  # exact posterior of the states (Gaussian conditioning)
  P <- solve(solve(Smu) + diag(n) / s2e)
  mpost <- as.vector(P %*% (solve(Smu, rep(m0, n)) + y / s2e))
  kf <- kalman_local_level(y, s2e, s2n, m0, C0)
  draws <- t(replicate(20000, ffbs_local_level(y, s2e, s2n, m0, C0, kf = kf)))
  z <- (colMeans(draws) - mpost) / (sqrt(diag(P)) / sqrt(nrow(draws)))   # in Monte Carlo standard errors
  expect_lt(max(abs(z)), 4.5)
  expect_equal(apply(draws, 2, var), diag(P), tolerance = 0.05)
})

test_that("AR design and stationarity check", {
  d <- ar_design(1:10, 2)
  expect_equal(d$X[1, ], c(c = 1, phi1 = 2, phi2 = 1))
  expect_equal(d$y[1], 3)
  expect_true(ar_is_stationary(c(0.5, 0.3)))
  expect_false(ar_is_stationary(c(0.5, 0.6)))
  expect_false(ar_is_stationary(1))
})
