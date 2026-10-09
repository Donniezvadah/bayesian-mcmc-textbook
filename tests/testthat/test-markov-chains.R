test_that("stationary distribution, matrix powers and reversibility of finite chains", {
  P <- matrix(c(0.5, 0.5, 0, 0.25, 0.5, 0.25, 0, 0.5, 0.5), 3, byrow = TRUE)   # birth-death chain
  p <- stationary_distribution(P)
  expect_equal(p, c(0.25, 0.5, 0.25))
  expect_equal(as.vector(p %*% P), p)
  expect_true(is_reversible(P, p))
  expect_equal(matrix_power(P, 7), P %*% P %*% P %*% P %*% P %*% P %*% P)
  expect_equal(matrix_power(P, 0), diag(3))
  # a non-reversible chain: deterministic-ish cycle with a lazy step
  Q <- matrix(c(0.1, 0.9, 0, 0, 0.1, 0.9, 0.9, 0, 0.1), 3, byrow = TRUE)
  q <- stationary_distribution(Q)
  expect_equal(q, rep(1 / 3, 3))
  expect_false(is_reversible(Q, q))
  expect_error(stationary_distribution(matrix(c(0.5, 0.6, 0.5, 0.4), 2, byrow = TRUE)), "rows summing")
})

test_that("two-state chain: geometric convergence and ergodic averages", {
  a <- 0.3; b <- 0.1
  P <- matrix(c(1 - a, a, b, 1 - b), 2, byrow = TRUE)
  p <- c(b, a) / (a + b)
  lambda2 <- 1 - a - b
  for (n in c(1, 5, 20)) {
    expect_equal(tv_distance(matrix_power(P, n)[1, ], p), abs(lambda2)^n * a / (a + b), tolerance = 1e-8)
  }
  x <- with_seed(120, simulate_chain(P, 2e5, 1))
  f <- (x == 2) * 1
  expect_equal(mean(f), p[2], tolerance = 0.01)
  # IAT of the indicator equals (1 + lambda2) / (1 - lambda2)
  expect_equal(iat(f), (1 + lambda2) / (1 - lambda2), tolerance = 0.1)
})

test_that("AR(1) simulation has the stationary variance", {
  x <- with_seed(121, simulate_ar1(1e5, 0.8))
  expect_equal(var(x), 1 / (1 - 0.8^2), tolerance = 0.05)
  expect_error(simulate_ar1(10, 1), "phi")
})
