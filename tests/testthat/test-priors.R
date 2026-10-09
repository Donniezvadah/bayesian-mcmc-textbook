test_that("Beta elicitation from mean/sd and from quantiles", {
  ab <- beta_from_mean_sd(0.1, 0.05)
  expect_equal(unname(ab), c(3.5, 31.5))
  expect_error(beta_from_mean_sd(0.5, 0.6), "No Beta")
  ab2 <- beta_from_quantiles(c(0.05, 0.30), c(0.05, 0.95))
  expect_equal(qbeta(c(0.05, 0.95), ab2[1], ab2[2]), c(0.05, 0.30), tolerance = 1e-5)
})

test_that("prior predictive simulation returns parameter-data pairs", {
  sims <- with_seed(90, prior_predictive(500, function() rbeta(1, 2, 2), function(t) rbinom(1, 10, t)))
  expect_length(sims, 500)
  y <- vapply(sims, `[[`, numeric(1), "y")
  expect_within_mcse(mean(y), 5, sd(y) / sqrt(500))
})

test_that("prior-data conflict p-value is uniform-ish under the prior and small under conflict", {
  # no conflict: y drawn from the prior predictive itself
  p <- with_seed(91, replicate(2000, prior_data_conflict_betabinom(rbetabinom(1, 30, 2, 8), 30, 2, 8)))
  expect_gt(mean(p <= 0.05), 0.01); expect_lt(mean(p <= 0.05), 0.08)   # discrete: at most nominal
  expect_lt(prior_data_conflict_betabinom(28, 30, 2, 20), 0.001)
})
