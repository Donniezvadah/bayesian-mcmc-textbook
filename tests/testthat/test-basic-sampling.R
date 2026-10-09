test_that("inverse-CDF sampling reproduces an exponential distribution", {
  x <- with_seed(70, inverse_cdf_sample(1e5, function(u) -log(1 - u) / 2))
  expect_within_mcse(mean(x), 0.5, sd(x) / sqrt(length(x)))
  expect_gt(ks.test(x, "pexp", rate = 2)$p.value, 0.001)
})

test_that("rejection sampling gives exact Beta draws at rate Z / M", {
  a <- 2.7; b <- 6.3
  mode <- (a - 1) / (a + b - 2)
  log_M <- dbeta(mode, a, b, log = TRUE)            # envelope M * U(0,1) density
  out <- with_seed(71, rejection_sample(5e4, function(x) dbeta(x, a, b, log = TRUE),
                                        runif, function(x) dunif(x, log = TRUE), log_M))
  expect_length(out$draws, 5e4)
  expect_within_mcse(mean(out$draws), a / (a + b), sd(out$draws) / sqrt(5e4))
  expect_equal(out$acceptance_rate, exp(-log_M), tolerance = 0.02)   # Z = 1 here
  expect_gt(ks.test(out$draws, "pbeta", a, b)$p.value, 0.001)
})

test_that("rejection sampling detects a violated envelope", {
  expect_error(with_seed(72, rejection_sample(100, function(x) dnorm(x, log = TRUE), rnorm,
                                              function(x) dnorm(x, sd = 0.5, log = TRUE), log_M = 0)),
               "Envelope violated")
})
