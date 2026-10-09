test_that("log_sum_exp is exact on moderate inputs and stable on extreme ones", {
  x <- c(-1, 0.5, 2)
  expect_equal(log_sum_exp(x), log(sum(exp(x))))
  # exp(1000) overflows; the stable version does not
  expect_equal(log_sum_exp(c(1000, 1000)), 1000 + log(2))
  # exp(-1000) underflows to 0; the stable version keeps the information
  expect_equal(log_sum_exp(c(-1000, -1000)), -1000 + log(2))
  expect_equal(log_sum_exp(c(-Inf, -Inf)), -Inf)
  expect_equal(log_sum_exp(c(-Inf, 0)), 0)
  expect_error(log_sum_exp(numeric(0)))
})

test_that("normalise_log_weights returns a probability vector", {
  w <- normalise_log_weights(c(-1000, -1001, -1002))
  expect_equal(sum(w), 1)
  expect_equal(w[1] / w[2], exp(1))
  expect_error(normalise_log_weights(c(-Inf, -Inf)))
})

test_that("log1mexp and log1pexp are accurate in both tails", {
  a <- c(1e-20, 1e-8, 0.1, log(2), 1, 40, 800)
  ref <- log(-expm1(-a))
  ref[a > log(2)] <- log1p(-exp(-a[a > log(2)]))
  expect_equal(log1mexp(a), ref)
  expect_equal(log1mexp(1e-20), log(1e-20), tolerance = 1e-12)   # naive log(1 - exp(-a)) gives -Inf
  expect_equal(log1pexp(800), 800)
  expect_equal(log1pexp(-800), exp(-800))
  expect_equal(log1pexp(0), log(2))
})

test_that("Simpson's rule is exact for cubics and accurate for smooth integrands", {
  x <- seq(0, 2, length.out = 11)
  expect_equal(simpson(x^3 - x, x), 2^4 / 4 - 2^2 / 2)
  x <- seq(-8, 8, length.out = 2001)
  expect_equal(simpson(dnorm(x), x), 1, tolerance = 1e-10)
  expect_error(simpson(x[-1]^2, x[-1]))       # even number of points
})

test_that("grid_posterior normalises a Beta kernel correctly", {
  g <- seq(0, 1, length.out = 2001)
  lp <- dbeta(g, 3, 5, log = TRUE)
  gp <- grid_posterior(lp + 17, g)            # arbitrary additive constant
  expect_equal(gp$density, dbeta(g, 3, 5), tolerance = 1e-8)
  expect_equal(gp$log_normalising_constant, 17, tolerance = 1e-8)
})

test_that("numerical_gradient agrees with an analytic gradient", {
  f <- function(x) sum(sin(x) * x^2)
  g <- function(x) cos(x) * x^2 + 2 * x * sin(x)
  x <- c(-1.3, 0.2, 2.5)
  expect_equal(numerical_gradient(f, x), g(x), tolerance = 1e-8)
})

test_that("rmvnorm_chol and dmvnorm_log are consistent", {
  Sigma <- matrix(c(2, 0.6, 0.6, 1), 2)
  mu <- c(1, -1)
  x <- with_seed(1, rmvnorm_chol(20000, mu, Sigma))
  expect_equal(colMeans(x), mu, tolerance = 0.05)
  expect_equal(cov(x), Sigma, tolerance = 0.05)
  # density against a direct formula
  pt <- c(0.3, -0.2)
  r <- pt - mu
  ref <- -log(2 * pi) - 0.5 * log(det(Sigma)) - 0.5 * sum(r * solve(Sigma, r))
  expect_equal(dmvnorm_log(pt, mu, Sigma), ref)
  expect_error(rmvnorm_chol(10, mu, matrix(c(1, 2, 2, 1), 2)))   # not positive definite
})

test_that("with_seed reproduces draws and restores the RNG state", {
  set.seed(99); before <- runif(1); set.seed(99)
  a <- with_seed(1, rnorm(3)); b <- with_seed(1, rnorm(3))
  expect_identical(a, b)
  expect_identical(runif(1), before)          # outer stream undisturbed
})

test_that("validation helpers reject invalid input", {
  expect_error(check_count(2.5, "n"), "whole number")
  expect_error(check_count(0, "n"), "whole number")
  expect_error(check_scalar(c(1, 2), "x"), "single finite")
  expect_error(check_scalar(-1, "x", positive = TRUE), "strictly positive")
  expect_error(check_spd_matrix(matrix(c(1, 0, 1, 1), 2), "S"), "symmetric")
  expect_equal(eval_log_density(function(x) NaN, 1), -Inf)
  expect_equal(eval_log_density(function(x) NA_real_, 1), -Inf)
  expect_error(eval_log_density(function(x) Inf, 1), "\\+Inf")
  expect_error(eval_log_density(function(x) c(1, 2), 1), "single numeric")
})

test_that("multivariate t sampler and density are consistent", {
  Sigma <- matrix(c(1, 0.3, 0.3, 2), 2); mu <- c(1, 0); nu <- 5
  x <- with_seed(130, rmvt_chol(2e5, mu, Sigma, nu))
  expect_equal(colMeans(x), mu, tolerance = 0.02)
  expect_equal(cov(x), Sigma * nu / (nu - 2), tolerance = 0.05)
  # one-dimensional case agrees with dt()
  expect_equal(dmvt_log(c(0.7, 1.9), 0, matrix(4), 3), dt(c(0.7, 1.9) / 2, 3, log = TRUE) - log(2))
})

test_that("Gauss-Hermite rule integrates polynomials and smooth functions against phi", {
  gh <- gauss_hermite(20)
  expect_equal(sum(gh$weights), 1)
  expect_equal(sum(gh$weights * gh$nodes^2), 1)
  expect_equal(sum(gh$weights * gh$nodes^4), 3)
  expect_equal(sum(gh$weights * exp(0.5 * gh$nodes)), exp(0.125), tolerance = 1e-10)   # E[e^{aZ}] = e^{a^2/2}
})
