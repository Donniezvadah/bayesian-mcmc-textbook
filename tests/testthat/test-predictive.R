test_that("CRPS from draws matches the closed form for a normal predictive", {
  crps_normal <- function(y, m, s) { z <- (y - m) / s; s * (z * (2 * pnorm(z) - 1) + 2 * dnorm(z) - 1 / sqrt(pi)) }
  x <- with_seed(100, rnorm(1e5, 1, 2))
  for (y in c(-2, 1, 4)) expect_equal(crps_sample(y, x), crps_normal(y, 1, 2), tolerance = 0.01)
})

test_that("PIT values are uniform for a correct predictive distribution (continuous and discrete)", {
  set.seed(101)
  y <- rnorm(2000); draws <- matrix(rnorm(500 * 2000), 500)
  unif_p <- function(u) chisq.test(table(cut(u, seq(0, 1, 0.1), include.lowest = TRUE)))$p.value
  expect_gt(unif_p(pit_values(y, draws)), 0.001)   # binned test: PIT from draws has ties
  yc <- rpois(2000, 3); dc <- matrix(rpois(2000 * 2000, 3), 2000)
  expect_gt(unif_p(pit_values(yc, dc, discrete = TRUE)), 0.001)
})

test_that("posterior predictive p-values and lppd behave as defined", {
  expect_equal(ppc_pvalue(2, c(1, 2, 3, 4)), 0.75)
  expect_equal(ppc_pvalue(c(5, 0), c(1, 1)), 0.5)
  ll <- matrix(log(c(0.2, 0.4, 0.1, 0.3)), 2)
  expect_equal(lppd_pointwise(ll), log(c(0.3, 0.2)))
})
