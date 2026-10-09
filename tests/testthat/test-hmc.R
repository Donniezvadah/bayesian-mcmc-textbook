test_that("model gradients agree with finite differences", {
  set.seed(30)
  targets <- list(target_mvn(c(1, -1), matrix(c(2, 0.5, 0.5, 1), 2)), target_funnel(5),
                  target_banana(0.05), target_mixture(c(0.3, 0.7), c(-2, 2), c(0.5, 1)),
                  eight_schools_centered(), eight_schools_noncentered())
  for (tg in targets) {
    x <- rnorm(tg$dim)
    expect_equal(tg$grad(x), numerical_gradient(tg$log_density, x), tolerance = 1e-6)
  }
  X <- cbind(1, rnorm(30), rnorm(30))
  yb <- rbinom(30, 1, 0.4); yc <- rpois(30, 2)
  for (m in list(logistic_regression_model(X, yb), poisson_regression_model(X, yc))) {
    b <- rnorm(3, 0, 0.3)
    expect_equal(m$grad(b), numerical_gradient(m$log_density, b), tolerance = 1e-6)
  }
})

test_that("leapfrog is time-reversible and has O(eps^2) energy error", {
  tg <- target_mvn(c(0, 0), matrix(c(1, 0.9, 0.9, 1), 2))
  q0 <- c(1, -0.5); p0 <- c(0.3, 0.8)
  fw <- leapfrog(q0, p0, tg$grad, 0.1, 25)
  bw <- leapfrog(fw$q, -fw$p, tg$grad, 0.1, 25)       # flip momentum and integrate back
  expect_equal(bw$q, q0, tolerance = 1e-10)
  expect_equal(-bw$p, p0, tolerance = 1e-10)
  H <- function(q, p) -tg$log_density(q) + sum(p^2) / 2
  err <- sapply(c(0.04, 0.02, 0.01), function(e) {
    tr <- leapfrog(q0, p0, tg$grad, e, round(1 / e))  # fixed integration time 1
    abs(H(tr$q, tr$p) - H(q0, p0))
  })
  slopes <- diff(log(err)) / diff(log(c(0.04, 0.02, 0.01)))
  expect_true(all(abs(slopes - 2) < 0.3))
})

test_that("HMC recovers moments of a correlated Gaussian and of N(3, 2^2)", {
  tg <- target_mvn(c(1, -2), matrix(c(1, 0.95, 0.95, 1), 2))
  out <- with_seed(31, hmc_sampler(tg$log_density, tg$grad, c(0, 0), 4000, step_size = 0.15, n_steps = 20))
  d <- out$draws[-(1:200), ]
  expect_within_mcse(mean(d[, 1]), 1, mcse_mean(d[, 1]))
  expect_within_mcse(mean(d[, 2]), -2, mcse_mean(d[, 2]))
  expect_within_mcse(mean((d[, 1] - 1) * (d[, 2] + 2)), 0.95, mcse_mean((d[, 1] - 1) * (d[, 2] + 2)))
  expect_false(any(out$divergent))
  t1 <- target_normal(3, 2)
  o1 <- with_seed(32, hmc_sampler(t1$log_density, t1$grad, 0, 4000, 0.5, 10))$draws[-(1:200), 1]
  expect_within_mcse(var(o1), 4, mcse_mean((o1 - mean(o1))^2))
  expect_error(hmc_sampler(t1$log_density, t1$grad, 0, 10, -0.1, 10), "strictly positive")
  expect_error(hmc_sampler(t1$log_density, t1$grad, 0, 10, 0.1, 0), "whole number")
})

test_that("HMC flags divergences in the funnel's neck when the step size is too large", {
  tg <- target_funnel(10)
  out <- with_seed(33, hmc_sampler(tg$log_density, tg$grad, c(-3, rep(0.1, 9)), 300, 1.2, 10))
  expect_gt(sum(out$divergent), 0)
})

test_that("dual averaging finds a step size with acceptance near the target", {
  tg <- target_mvn(rep(0, 5), diag(c(1, 4, 9, 0.25, 1)))
  wu <- with_seed(34, hmc_warmup(tg$log_density, tg$grad, rep(1, 5), 600, n_steps = 10,
                                 target_accept = 0.8))
  out <- with_seed(35, hmc_sampler(tg$log_density, tg$grad, wu$position, 1500, wu$step_size, 10, jitter = FALSE))
  expect_gt(out$mean_accept_prob, 0.65)
  expect_lt(out$mean_accept_prob, 0.95)
})
