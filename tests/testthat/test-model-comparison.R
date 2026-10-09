test_that("IS-LOO and WAIC approximate exact leave-one-out for a conjugate normal model", {
  set.seed(150)
  y <- rnorm(30, 1, 1); n <- length(y)
  # model: y_i ~ N(mu, 1), mu ~ N(0, 10^2); exact LOO predictive is normal
  exact <- sapply(seq_len(n), function(i) {
    p <- normal_known_var_update(y[-i], 1, 0, 100)
    dnorm(y[i], p$mean, sqrt(1 + p$var), log = TRUE)
  })
  post <- normal_known_var_update(y, 1, 0, 100)
  mu <- rnorm(20000, post$mean, sqrt(post$var))
  ll <- sapply(y, function(yi) dnorm(yi, mu, 1, log = TRUE))
  expect_equal(loo_is(ll)$elpd_loo, sum(exact), tolerance = 0.01)
  expect_equal(waic(ll)$elpd_waic, sum(exact), tolerance = 0.01)
  skip_if_not_installed("loo")
  ref <- suppressWarnings(loo::loo(ll, r_eff = rep(1, n)))
  expect_equal(waic(ll)$elpd_waic, unname(suppressWarnings(loo::waic(ll))$estimates["elpd_waic", "Estimate"]), tolerance = 1e-8)
  expect_equal(ref$estimates["elpd_loo", "Estimate"], sum(exact), tolerance = 0.01)
})
