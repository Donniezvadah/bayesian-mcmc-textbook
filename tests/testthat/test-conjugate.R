# Every conjugate update is checked against brute-force numerical integration of
# prior x likelihood, so the algebra in Chapter 5 is verified independently.

test_that("Beta-Binomial posterior matches numerical normalisation", {
  post <- beta_binomial_update(y = 7, n = 20, a = 2, b = 3)
  expect_equal(post, list(a = 9, b = 16))
  g <- seq(0, 1, length.out = 4001)
  lp <- dbinom(7, 20, g, log = TRUE) + dbeta(g, 2, 3, log = TRUE)
  gp <- grid_posterior(lp, g)
  expect_equal(gp$density[-c(1, 4001)], dbeta(g, 9, 16)[-c(1, 4001)], tolerance = 1e-6)
  expect_error(beta_binomial_update(y = 21, n = 20))
  expect_error(beta_binomial_update(y = 2, n = 5, a = 0))
})

test_that("Beta-Binomial predictive pmf sums to one and matches simulation", {
  pmf <- exp(dbetabinom_log(0:10, 10, 2.5, 4))
  expect_equal(sum(pmf), 1)
  sims <- with_seed(2, rbetabinom(2e5, 10, 2.5, 4))
  expect_equal(as.vector(table(factor(sims, levels = 0:10))) / 2e5, pmf, tolerance = 0.02)
  # mean of the Beta-Binomial is m a / (a + b)
  expect_equal(sum(0:10 * pmf), 10 * 2.5 / 6.5)
})

test_that("Gamma-Poisson posterior and Negative-Binomial predictive are correct", {
  y <- c(3, 0, 2, 5, 1); t <- c(1, 2, 1, 1.5, 1)
  post <- gamma_poisson_update(y, a = 2, b = 1, exposure = t)
  expect_equal(post, list(a = 2 + 11, b = 1 + 6.5))
  g <- seq(1e-6, 8, length.out = 4001)
  lp <- sapply(g, function(l) sum(dpois(y, l * t, log = TRUE))) + dgamma(g, 2, 1, log = TRUE)
  expect_equal(grid_posterior(lp, g)$density, dgamma(g, post$a, post$b), tolerance = 1e-5)
  # predictive: integrate Poisson(k | lambda t*) against the posterior
  for (k in 0:4) {
    f <- function(l) dpois(k, l * 2) * dgamma(l, post$a, post$b)
    expect_equal(exp(dgamma_poisson_pred_log(k, post$a, post$b, t = 2)),
                 integrate(f, 0, Inf)$value, tolerance = 1e-6)
  }
})

test_that("Gamma-Exponential posterior and Lomax predictive are correct", {
  y <- c(0.4, 1.2, 0.1, 2.3)
  post <- gamma_exponential_update(y, a = 1.5, b = 2)
  expect_equal(post, list(a = 5.5, b = 2 + sum(y)))
  for (yn in c(0, 0.5, 3)) {
    f <- function(l) dexp(yn, l) * dgamma(l, post$a, post$b)
    expect_equal(exp(dgamma_exponential_pred_log(yn, post$a, post$b)),
                 integrate(f, 0, Inf)$value, tolerance = 1e-6)
  }
  expect_equal(integrate(function(z) exp(dgamma_exponential_pred_log(z, 5.5, 6)), 0, Inf)$value, 1,
               tolerance = 1e-6)
})

test_that("Normal model with known variance: posterior mean and variance", {
  y <- c(4.1, 5.3, 3.8, 6.0, 5.1); s2 <- 1.5^2
  post <- normal_known_var_update(y, s2, m0 = 0, v0 = 10)
  g <- seq(-5, 15, length.out = 8001)
  lp <- sapply(g, function(m) sum(dnorm(y, m, 1.5, log = TRUE))) + dnorm(g, 0, sqrt(10), log = TRUE)
  d <- grid_posterior(lp, g)$density
  m_grid <- simpson(g * d, g)
  expect_equal(post$mean, m_grid, tolerance = 1e-8)
  expect_equal(post$var, simpson((g - m_grid)^2 * d, g), tolerance = 1e-6)
})

test_that("Normal-Inverse-Gamma: marginal of mu and predictive match 2-D integration", {
  y <- c(2.1, 3.4, 1.9, 2.8, 3.0, 2.2)
  post <- nig_update(y, m0 = 0, k0 = 0.5, a0 = 2, b0 = 1)
  # joint log posterior on a grid of (mu, sigma2)
  mu_g <- seq(0, 5, length.out = 801); s2_g <- seq(0.005, 10, length.out = 1601)
  lj <- outer(mu_g, s2_g, Vectorize(function(m, s2)
    sum(dnorm(y, m, sqrt(s2), log = TRUE)) + dnorm(m, 0, sqrt(s2 / 0.5), log = TRUE) +
      dinvgamma_log(s2, 2, 1)))
  J <- exp(lj - max(lj))
  marg_mu <- apply(J, 1, function(r) simpson(r, s2_g))
  marg_mu <- marg_mu / simpson(marg_mu, mu_g)
  t_par <- nig_marginal_mu(post)
  expect_equal(marg_mu, dt_ls(mu_g, t_par$df, t_par$location, t_par$scale), tolerance = 2e-3)
  # predictive t: compare with Monte Carlo from rnig
  draws <- with_seed(3, rnig(2e5, post))
  ynew <- rnorm(2e5, draws[, "mu"], sqrt(draws[, "sigma2"]))
  pr <- nig_predictive(post)
  expect_equal(quantile(ynew, c(0.1, 0.5, 0.9), names = FALSE),
               qt_ls(c(0.1, 0.5, 0.9), pr$df, pr$location, pr$scale), tolerance = 0.02)
})

test_that("Dirichlet-Multinomial update, sampler and predictive", {
  alpha <- c(1, 2, 3); counts <- c(5, 0, 2)
  post <- dirichlet_multinomial_update(counts, alpha)
  expect_equal(post, c(6, 2, 5))
  p <- with_seed(4, rdirichlet(1e5, post))
  expect_equal(colMeans(p), post / sum(post), tolerance = 0.01)
  expect_true(all(abs(rowSums(p) - 1) < 1e-12))
  # predictive pmf over all compositions of 3 future observations sums to one
  comps <- expand.grid(a = 0:3, b = 0:3); comps$c <- 3 - comps$a - comps$b
  comps <- comps[comps$c >= 0, ]
  tot <- sum(apply(comps, 1, function(x) exp(ddirmult_log(x, post))))
  expect_equal(tot, 1)
})

test_that("inverse-gamma and location-scale t helpers are proper densities", {
  expect_equal(integrate(function(x) exp(dinvgamma_log(x, 3, 2)), 0, Inf)$value, 1, tolerance = 1e-6)
  expect_equal(integrate(function(x) dt_ls(x, 5, 2, 3), -Inf, Inf)$value, 1, tolerance = 1e-6)
  expect_equal(pt_ls(qt_ls(0.3, 5, 2, 3), 5, 2, 3), 0.3)
})
