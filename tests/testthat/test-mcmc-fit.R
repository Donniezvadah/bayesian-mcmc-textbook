test_that("run_chains is reproducible and independent of parallelism", {
  tg <- target_normal(1, 2)
  smp <- function(m) { d <- metropolis_hastings(tg$log_density, 0, 500, 4)$draws; cbind(theta = d) }
  a <- run_chains(smp, 3, seed = 7)
  b <- run_chains(smp, 3, seed = 7)
  expect_identical(a$draws, b$draws)
  expect_false(identical(a$draws[, 1, 1], a$draws[, 2, 1]))       # chains use different streams
  skip_on_os("windows")
  c2 <- run_chains(smp, 3, seed = 7, cores = 2)
  expect_identical(a$draws, c2$draws)
  expect_s3_class(a, "mcmc_fit")
  expect_equal(nrow(summary(a)), 1L)
  expect_identical(RNGkind()[1], "Mersenne-Twister")               # RNG kind restored
})
