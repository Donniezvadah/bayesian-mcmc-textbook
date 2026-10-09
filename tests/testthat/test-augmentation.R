test_that("truncated normal draws are correct in the bulk and far in the tails", {
  set.seed(160)
  x <- rtruncnorm_sign(rep(0.5, 1e5), positive = TRUE)
  expect_true(all(x > 0))
  # mean of N(m,1) truncated to (0, Inf) is m + phi(m) / Phi(m)
  expect_equal(mean(x), 0.5 + dnorm(0.5) / pnorm(0.5), tolerance = 0.01)
  far <- rtruncnorm_sign(rep(-30, 1000), positive = TRUE)      # Phi(-30) underflows a naive method
  expect_true(all(is.finite(far) & far > 0))
  expect_equal(mean(far), 1 / 30, tolerance = 0.05)             # approx. Exp(30) excess
  neg <- rtruncnorm_sign(rep(1, 1e4), positive = FALSE)
  expect_true(all(neg < 0))
})

test_that("Albert-Chib probit sampler matches a grid posterior", {
  set.seed(161)
  X <- cbind(1, rnorm(40)); y <- rbinom(40, 1, pnorm(X %*% c(0.3, 1)))
  d <- albert_chib_probit(X, y, 30000, prior_sd = 5)[-(1:1000), ]
  lp <- function(b) sum(pnorm(X %*% b, log.p = TRUE) * y + pnorm(-(X %*% b), log.p = TRUE) * (1 - y)) +
    sum(dnorm(b, 0, 5, log = TRUE))
  g1 <- seq(-1.5, 2, length.out = 251); g2 <- seq(-0.5, 3.5, length.out = 251)
  L <- outer(g1, g2, Vectorize(function(a, b) lp(c(a, b)))); W <- exp(L - max(L)); W <- W / sum(W)
  expect_within_mcse(mean(d[, 1]), sum(W * g1), mcse_mean(d[, 1]) + 0.002)
  expect_within_mcse(mean(d[, 2]), sum(t(W) * g2), mcse_mean(d[, 2]) + 0.002)
})
