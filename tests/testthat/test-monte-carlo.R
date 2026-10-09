test_that("Monte Carlo intervals cover the truth at the nominal rate", {
  cover <- with_seed(110, replicate(2000, { e <- mc_estimate(rexp(200)); e["lower"] <= 1 && 1 <= e["upper"] }))
  expect_within_mcse(mean(cover), 0.95, sqrt(0.95 * 0.05 / 2000))
  expect_error(mc_estimate(1), "two")
})

test_that("antithetic variates reduce variance for a monotone integrand", {
  a <- with_seed(111, antithetic_estimate(exp, 5000))
  expect_within_mcse(a[["estimate"]], exp(1) - 1, a[["mcse"]])
  # exact variance of a pair mean: (Var + Cov) / 2 = 0.003912
  v_pair <- ((exp(2) - 1) / 2 - (exp(1) - 1)^2 + exp(1) - (exp(1) - 1)^2) / 2
  expect_equal(a[["mcse"]]^2 * 5000, v_pair, tolerance = 0.1)
})

test_that("control variates are unbiased and reduce variance by 1 - rho^2", {
  x <- with_seed(112, rnorm(1e4))
  cv <- control_variate_estimate(exp(x), x, 0)            # E[e^X] = e^{1/2}
  expect_within_mcse(cv[["estimate"]], exp(0.5), cv[["mcse"]])
  rho2 <- cor(exp(x), x)^2
  expect_equal(cv[["variance_ratio"]], 1 - rho2, tolerance = 1e-6)
})
