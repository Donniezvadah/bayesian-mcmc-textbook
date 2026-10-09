test_that("bivariate-normal Gibbs has the right moments and AR(1) mixing", {
  rho <- 0.9
  d <- with_seed(20, gibbs_bivariate_normal(rho, 40000))[-(1:500), ]
  expect_within_mcse(mean(d[, 1]), 0, mcse_mean(d[, 1]))
  expect_within_mcse(mean(d[, 1] * d[, 2]), rho, mcse_mean(d[, 1] * d[, 2]))
  # x1-chain is AR(1) with coefficient rho^2: lag-1 autocorrelation and IAT
  expect_equal(autocorrelation(d[, 1])[2], rho^2, tolerance = 0.02)
  expect_equal(iat(d[, 1]), (1 + rho^2) / (1 - rho^2), tolerance = 0.15)
  expect_error(gibbs_bivariate_normal(1, 10), "\\(-1, 1\\)")
})

test_that("semi-conjugate normal Gibbs agrees with a 2-D grid posterior", {
  y <- with_seed(21, rnorm(15, 4, 2))
  d <- with_seed(22, gibbs_normal(y, 30000, m0 = 0, v0 = 100, a0 = 2, b0 = 2))[-(1:500), ]
  mu_g <- seq(0, 8, length.out = 401); s2_g <- seq(0.3, 20, length.out = 801)
  lj <- outer(mu_g, s2_g, Vectorize(function(m, s2)
    sum(dnorm(y, m, sqrt(s2), log = TRUE)) + dnorm(m, 0, 10, log = TRUE) + dinvgamma_log(s2, 2, 2)))
  J <- exp(lj - max(lj)); J <- J / sum(J)
  e_mu <- sum(J * mu_g); e_s2 <- sum(t(J) * s2_g)
  expect_within_mcse(mean(d[, "mu"]), e_mu, mcse_mean(d[, "mu"]))
  expect_within_mcse(mean(d[, "sigma2"]), e_s2, mcse_mean(d[, "sigma2"]))
})

test_that("Gibbs for linear regression agrees with the conjugate posterior in the limit of a vague prior", {
  X <- cbind(1, with_seed(23, rnorm(40)))
  y <- as.vector(X %*% c(1, 2)) + with_seed(24, rnorm(40, 0, 0.5))
  d <- with_seed(25, gibbs_linear_regression(X, y, 20000, S0 = diag(1e6, 2), a0 = 0.001, b0 = 0.001))[-(1:500), ]
  ols <- lm.fit(X, y)$coefficients
  expect_within_mcse(mean(d[, 1]), ols[1], mcse_mean(d[, 1]) + 1e-3)
  expect_within_mcse(mean(d[, 2]), ols[2], mcse_mean(d[, 2]) + 1e-3)
})

test_that("centred and non-centred hierarchical Gibbs target the same posterior", {
  y <- eight_schools$y; s <- eight_schools$sigma
  ch_c <- with_seed(26, replicate(4, gibbs_hier_centered(y, s, 15000)[-(1:2000), c("mu", "theta[1]")], simplify = FALSE))
  ch_n <- with_seed(27, replicate(4, gibbs_hier_noncentered(y, s, 15000)[-(1:2000), c("mu", "theta[1]")], simplify = FALSE))
  for (v in c("mu", "theta[1]")) {
    a <- sapply(ch_c, function(m) m[, v]); b <- sapply(ch_n, function(m) m[, v])
    se <- sqrt(mcse_mean(a)^2 + mcse_mean(b)^2)
    expect_within_mcse(mean(a), mean(b), se)
  }
})

test_that("exact eight-schools draws agree with the non-centred Gibbs sampler", {
  ex <- with_seed(28, eight_schools_exact(40000))
  gi <- with_seed(29, gibbs_hier_noncentered(eight_schools$y, eight_schools$sigma, 60000)[-(1:2000), ])
  for (v in c("mu", "tau", "theta[1]")) {
    se <- sqrt(sd(ex[, v])^2 / 40000 + mcse_mean(gi[, v])^2)
    expect_within_mcse(mean(gi[, v]), mean(ex[, v]), se)
  }
})

test_that("Poisson-lognormal model gradient is correct", {
  X <- cbind(1, c(0, 1, 0, 1, 1)); y <- c(2, 5, 0, 9, 4)
  m <- poisson_lognormal_model(X, y)
  par <- with_seed(30, rnorm(m$dim, 0, 0.3))
  expect_equal(m$grad(par), numerical_gradient(m$log_density, par), tolerance = 1e-6)
})
