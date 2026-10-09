test_that("Newton's method finds the logistic-regression MAP found by optim()", {
  set.seed(60)
  X <- cbind(1, rnorm(80)); y <- rbinom(80, 1, plogis(X %*% c(-0.5, 1.2)))
  m <- logistic_regression_model(X, y, prior_sd = 2.5)
  hess <- function(b) { p <- plogis(as.vector(X %*% b)); -crossprod(X * (p * (1 - p)), X) - diag(1 / 2.5^2, 2) }
  nt <- newton_maximise(m$log_density, m$grad, hess, c(0, 0))
  ref <- optim(c(0, 0), m$log_density, m$grad, method = "BFGS", control = list(fnscale = -1, reltol = 1e-14))
  expect_true(nt$converged)
  expect_lt(nt$iterations, 10)                      # quadratic convergence
  expect_equal(nt$par, ref$par, tolerance = 1e-6)
  # analytic Hessian agrees with finite differences of the analytic gradient
  H_fd <- t(sapply(1:2, function(j) numerical_gradient(function(b) m$grad(b)[j], nt$par)))
  expect_equal(hess(nt$par), H_fd, tolerance = 1e-6)
})

test_that("the Laplace approximation is exact for a Gaussian and close for a Gamma kernel", {
  Sigma <- matrix(c(2, 0.3, 0.3, 0.5), 2); P <- solve(Sigma); mu <- c(1, -1)
  lp <- function(x) -0.5 * sum((x - mu) * (P %*% (x - mu)))
  la <- laplace_approximation(lp, function(x) -as.vector(P %*% (x - mu)), function(x) -P, c(0, 0))
  expect_equal(la$mode, mu, tolerance = 1e-8)
  expect_equal(la$cov, Sigma, tolerance = 1e-8)
  expect_equal(la$log_Z, log(2 * pi) + 0.5 * log(det(Sigma)), tolerance = 1e-8)
  # Gamma(a, 1) kernel on the log scale: h(u) = a u - exp(u), Z = Gamma(a)
  a <- 50
  la2 <- laplace_approximation(function(u) a * u - exp(u), function(u) a - exp(u),
                               function(u) matrix(-exp(u)), 1)
  expect_equal(la2$log_Z, lgamma(a), tolerance = 1e-3)  # Stirling: relative error O(1/a)
})
