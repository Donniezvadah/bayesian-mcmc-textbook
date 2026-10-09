test_that("forward-mode AD reproduces analytic derivatives", {
  f <- function(x) exp(sin(x[1])) * x[1]^2 / sqrt(x[1])
  x0 <- 1.3
  analytic <- exp(sin(x0)) * (cos(x0) * x0^1.5 + 1.5 * x0^0.5)
  expect_equal(ad_gradient(f, x0), analytic, tolerance = 1e-12)
  g <- function(x) log1p(exp(x[1] * x[2])) - x[2] / x[1] + expm1(x[2]) - cos(x[1])
  x1 <- c(0.7, -0.4)
  expect_equal(ad_gradient(g, x1), numerical_gradient(g, x1), tolerance = 1e-6)
  expect_error(ad_gradient(function(x) tan(x[1]), 1), "not implemented")
})

test_that("AD gradient of a logistic log posterior matches the hand-coded gradient", {
  set.seed(1)
  x <- rnorm(50); y <- rbinom(50, 1, inv_logit(0.3 + 0.8 * x))
  m <- logistic_regression_model(cbind(1, x), y, prior_sd = 2.5)
  f <- function(b) {
    eta <- b[1] + b[2] * x
    sum(y * eta - log1p(exp(eta))) - (b[1]^2 + b[2]^2) / (2 * 2.5^2)
  }
  b0 <- c(0.2, 0.5)
  expect_equal(ad_gradient(f, b0), m$grad(b0), tolerance = 1e-10)
})

test_that("rejection ABC with a sufficient summary recovers the exact posterior", {
  set.seed(2)
  n <- 20; y <- rnorm(n, 1, 1); s0 <- 3
  post <- normal_known_var_update(y, 1, 0, s0^2)
  abc <- abc_rejection(100000, function(k) rnorm(k, 0, s0), function(th) rnorm(n, th, 1),
                       mean, mean(y), keep = 0.01)
  expect_length(abc$draws, 1000)
  expect_within_mcse(mean(abc$draws), post$mean, sqrt(post$var / 1000))
  expect_equal(sd(abc$draws), sqrt(post$var), tolerance = 0.1)
  expect_true(abc$eps < 0.05)
})

test_that("mean-field ADVI recovers a diagonal Gaussian and underestimates correlated marginals", {
  set.seed(3)
  t1 <- target_mvn(c(1, -2), diag(c(0.5, 2)^2))
  v1 <- advi_meanfield(t1$log_density, t1$grad, c(0, 0), n_iter = 4000)
  expect_equal(v1$mean, c(1, -2), tolerance = 0.05)
  expect_equal(v1$sd, c(0.5, 2), tolerance = 0.05)
  t2 <- target_correlated_gaussian(0.9)
  v2 <- advi_meanfield(t2$log_density, t2$grad, c(0.5, 0.5), n_iter = 4000)
  expect_equal(v2$sd, rep(sqrt(1 - 0.9^2), 2), tolerance = 0.05)   # 1 / Lambda_ii
  expect_equal(dim(v2$draw(10)), c(10L, 2L))
})

test_that("SGLD stationary variance matches the closed form for the Gaussian mean model", {
  set.seed(4)
  y <- rnorm(1000, 2, 1); s0 <- 10; b <- 10; step <- 2e-4
  th <- sgld_gaussian_stationary(y, s0, b, step)
  draws <- sgld_sampler(function(t) -t / s0^2, function(t, rows) sum(rows - t), y, th[["mean"]],
                        n_iter = 40000, batch_size = b, step = step)[, 1]
  expect_within_mcse(mean(draws), th[["mean"]], sqrt(th[["var"]] / ess_mean(draws)))
  rel_se <- sqrt(2 / ess_mean(draws))                 # approximate relative s.e. of a variance
  expect_true(abs(var(draws) / th[["var"]] - 1) < 4 * rel_se)
  expect_gt(th[["var"]], 1.4 * th[["posterior_var"]])  # minibatch noise inflates the spread
  expect_error(sgld_gaussian_stationary(y, s0, b, 1), "unstable")
})
