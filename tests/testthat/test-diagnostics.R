ar1 <- function(n, phi, sd = 1) {
  x <- numeric(n); x[1] <- rnorm(1, 0, sd / sqrt(1 - phi^2))
  for (t in 2:n) x[t] <- phi * x[t - 1] + rnorm(1, 0, sd)
  x
}

test_that("FFT autocovariance equals the direct definition", {
  x <- with_seed(40, rnorm(200))
  direct <- sapply(0:20, function(k) sum((x[1:(200 - k)] - mean(x)) * (x[(1 + k):200] - mean(x))) / 200)
  expect_equal(autocovariance(x)[1:21], direct)
  expect_equal(autocorrelation(x)[1:21], as.vector(acf(x, lag.max = 20, plot = FALSE)$acf))
})

test_that("IAT of an AR(1) chain matches (1 + phi) / (1 - phi)", {
  for (phi in c(0.5, 0.9)) {
    x <- with_seed(41, ar1(2e5, phi))
    expect_equal(iat(x), (1 + phi) / (1 - phi), tolerance = 0.08)
  }
  # antithetic chain: ESS larger than the number of draws
  x <- with_seed(42, ar1(1e5, -0.5))
  expect_gt(ess_chains(x), 1e5)
})

test_that("ESS, R-hat and MCSE agree with the posterior package", {
  skip_if_not_installed("posterior")
  sims <- with_seed(43, cbind(ar1(1000, 0.7), ar1(1000, 0.7), ar1(1000, 0.7) + 0.1, ar1(1000, 0.8)))
  expect_equal(rhat_basic(sims), posterior::rhat_basic(sims), tolerance = 1e-8)
  expect_equal(rhat_rank(sims), posterior::rhat(sims), tolerance = 1e-8)
  expect_equal(ess_bulk(sims), posterior::ess_bulk(sims), tolerance = 0.03)
  expect_equal(ess_tail(sims), posterior::ess_tail(sims), tolerance = 0.03)
  expect_equal(ess_mean(sims), posterior::ess_mean(sims), tolerance = 0.03)
  expect_equal(mcse_mean(sims), posterior::mcse_mean(sims), tolerance = 0.03)
  expect_equal(mcse_quantile(sims, 0.9), posterior::mcse_quantile(sims, 0.9), tolerance = 0.05,
               ignore_attr = TRUE)
})

test_that("R-hat detects chains stuck in different places", {
  good <- with_seed(44, matrix(rnorm(4000), ncol = 4))
  bad <- good; bad[, 4] <- bad[, 4] + 2
  expect_lt(rhat_rank(good), 1.01)
  expect_gt(rhat_rank(bad), 1.1)
  # different scale but same location: caught by the folded (tail) R-hat
  scale_bad <- good; scale_bad[, 4] <- 4 * scale_bad[, 4]
  expect_lt(rhat_basic(z_scale(split_chains(scale_bad)), split = FALSE), 1.02)
  expect_gt(rhat_rank(scale_bad), 1.05)
})

test_that("batch-means and spectral MCSE estimate the asymptotic variance of AR(1)", {
  phi <- 0.8
  x <- with_seed(45, ar1(1e5, phi))
  sigma2_true <- 1 / (1 - phi)^2            # asymptotic variance of the mean: sigma_e^2/(1-phi)^2
  expect_equal(mcse_batch_means(x)^2 * length(x), sigma2_true, tolerance = 0.15)
  expect_equal(spectrum0_ar(x), sigma2_true, tolerance = 0.1)
  skip_if_not_installed("coda")
  expect_equal(spectrum0_ar(x), coda::spectrum0.ar(x)$spec, tolerance = 1e-8)
  expect_equal(geweke_z(x), unname(coda::geweke.diag(coda::mcmc(x))$z), tolerance = 1e-8)
})

test_that("simulation-based calibration ranks are uniform for an exact posterior", {
  ranks <- with_seed(46, sbc_ranks(
    n_sims = 1000, L = 99,
    r_prior = function() rbeta(1, 2, 2),
    r_data = function(theta) rbinom(1, 20, theta),
    posterior_draws = function(y, L) rbeta(L, 2 + y, 2 + 20 - y)))
  expect_true(all(ranks >= 0 & ranks <= 99))
  expect_gt(sbc_uniformity_test(ranks, 99)$p_value, 0.001)
  # A deliberately wrong, overconfident posterior (data counted twice) is detected
  wrong <- with_seed(47, sbc_ranks(
    n_sims = 1000, L = 99,
    r_prior = function() rbeta(1, 2, 2),
    r_data = function(theta) rbinom(1, 20, theta),
    posterior_draws = function(y, L) rbeta(L, 2 + 2 * y, 2 + 2 * (20 - y))))
  expect_lt(sbc_uniformity_test(wrong, 99)$p_value, 0.001)
})
