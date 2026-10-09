test_that("adaptive Metropolis learns the banana's scale and recovers its moments", {
  tg <- target_banana(0.03)
  out <- with_seed(50, adaptive_metropolis(tg$log_density, c(0, 0), 60000, adapt_start = 1000))
  d <- out$draws[-(1:10000), ]
  expect_within_mcse(mean(d[, 1]), 0, mcse_mean(d[, 1]))
  expect_within_mcse(mean(d[, 2]), 0, mcse_mean(d[, 2]))
  expect_equal(out$final_proposal_cov[1, 1] / (2.38^2 / 2), 100, tolerance = 0.3)
})

test_that("MALA recovers a 5-D Gaussian with heterogeneous scales", {
  tg <- target_mvn(rep(0, 5), diag(c(1, 2, 0.5, 1, 1.5)))
  out <- with_seed(51, mala_sampler(tg$log_density, tg$grad, rep(0, 5), 30000, step = 0.6))
  d <- out$draws[-(1:1000), ]
  for (j in 1:5) expect_within_mcse(mean(d[, j]^2), diag(tg$cov)[j], mcse_mean(d[, j]^2))
  expect_gt(out$acceptance_rate, 0.3)
})

test_that("slice sampler recovers Gamma(3, 2) including the boundary at 0", {
  tg <- target_gamma(3, 2)
  d <- with_seed(52, slice_sampler(tg$log_density, 1, 20000, w = 1))[-(1:200)]
  expect_true(all(d > 0))
  expect_within_mcse(mean(d), 1.5, mcse_mean(d))
  expect_within_mcse(mean(d^2), 1.5^2 + 0.75, mcse_mean(d^2))
})

test_that("elliptical slice sampling reproduces a conjugate Gaussian posterior", {
  # prior f ~ N(0, Sigma), likelihood y ~ N(f, s2 I): posterior is Gaussian.
  Sigma <- matrix(c(1, 0.5, 0.5, 1), 2); s2 <- 0.5; y <- c(1, -0.5)
  log_lik <- function(f) sum(dnorm(y, f, sqrt(s2), log = TRUE))
  V <- solve(solve(Sigma) + diag(2) / s2); m <- as.vector(V %*% (y / s2))
  d <- with_seed(53, elliptical_slice(log_lik, c(0, 0), 20000, Sigma))[-(1:200), ]
  expect_within_mcse(mean(d[, 1]), m[1], mcse_mean(d[, 1]))
  expect_within_mcse(mean(d[, 2]), m[2], mcse_mean(d[, 2]))
  expect_within_mcse(var(d[, 1]), V[1, 1], mcse_mean((d[, 1] - mean(d[, 1]))^2))
})

test_that("parallel tempering visits both well-separated modes in the right proportions", {
  tg <- target_mixture(c(0.3, 0.7), c(-8, 8), c(1, 1))
  pt <- with_seed(54, parallel_tempering(tg$log_density, 0, 30000,
                                         betas = c(1, 0.3, 0.1, 0.03), proposal_sd = c(2, 3.5, 6, 11)))
  d <- pt$draws[-(1:1000)]
  ind <- (d > 0) * 1
  expect_within_mcse(mean(ind), 0.7, mcse_mean(ind))
  # plain RWM from the left mode does not find the right one in the same budget
  rw <- with_seed(55, metropolis_hastings(tg$log_density, -8, 30000, 1.5))$draws
  expect_lt(mean(rw > 0), 0.05)
})

test_that("importance sampling estimates a mean and a normalising constant", {
  # target: unnormalised N(2, 1) kernel, Z = sqrt(2 pi); proposal t_5(0, 2)
  out <- with_seed(56, importance_sampling(
    f = function(x) x, log_target = function(x) -0.5 * (x - 2)^2,
    r_proposal = function(n) 2 * rt(n, 5), log_proposal = function(x) dt(x / 2, 5, log = TRUE) - log(2),
    n = 50000))
  expect_within_mcse(out$estimate, 2, out$se)
  expect_equal(out$log_Z, log(sqrt(2 * pi)), tolerance = 0.01)
  expect_equal(kish_ess(rep(1, 10)), 10)
  expect_equal(kish_ess(c(1, 0, 0)), 1)
})

test_that("tempering SMC recovers a conjugate posterior and its marginal likelihood", {
  y <- c(1.2, 0.8, 1.9, 1.4, 0.6, 1.1)
  # theta ~ N(0, 2^2), y_i ~ N(theta, 1): exact posterior and evidence
  post <- normal_known_var_update(y, 1, 0, 4)
  n <- length(y)
  log_ev <- dmvnorm_log(y, rep(0, n), diag(n) + 4)
  out <- with_seed(57, smc_tempering(
    log_lik = function(th) sum(dnorm(y, th, 1, log = TRUE)),
    log_prior = function(th) dnorm(th, 0, 2, log = TRUE),
    r_prior = function(n) rnorm(n, 0, 2), n = 4000))
  expect_equal(mean(out$draws), post$mean, tolerance = 0.03)
  expect_equal(var(as.vector(out$draws)), post$var, tolerance = 0.1)
  expect_equal(out$log_marginal_likelihood, log_ev, tolerance = 0.05)
  expect_equal(tail(out$temperatures, 1), 1)
})

test_that("conjugate linear regression matches the textbook formulas", {
  X <- cbind(1, c(-1, 0, 1, 2, 3)); y <- c(0.8, 2.1, 2.9, 4.2, 5.1)
  fit <- conjugate_linear_regression(X, y, V0 = diag(1e8, 2), a0 = 1e-8, b0 = 1e-8)
  expect_equal(fit$m, as.vector(solve(crossprod(X), crossprod(X, y))), tolerance = 1e-6)
  rss <- sum(lm.fit(X, y)$residuals^2)
  expect_equal(fit$b, rss / 2, tolerance = 1e-5)
  d <- with_seed(58, fit$draw(1e5))
  expect_equal(unname(colMeans(d[, 1:2])), fit$m, tolerance = 0.02)
})

test_that("Robbins-Monro scale adaptation reaches the target acceptance rate", {
  tg <- target_mvn(rep(0, 10), diag(10))
  out <- with_seed(59, rwm_adapt_scale(tg$log_density, rep(0, 10), 3000, 20000, target_accept = 0.234))
  expect_equal(out$acceptance_rate, 0.234, tolerance = 0.15)
  expect_equal(out$scale, 2.38 / sqrt(10), tolerance = 0.25)     # optimal scaling theory
  d <- out$draws
  expect_within_mcse(mean(d[, 1]^2), 1, mcse_mean(d[, 1]^2))
})

test_that("pCN reproduces a conjugate Gaussian posterior", {
  Sigma <- matrix(c(1, 0.5, 0.5, 1), 2); s2 <- 0.5; y <- c(1, -0.5)
  log_lik <- function(f) sum(dnorm(y, f, sqrt(s2), log = TRUE))
  V <- solve(solve(Sigma) + diag(2) / s2); m <- as.vector(V %*% (y / s2))
  out <- with_seed(62, pcn_sampler(log_lik, c(0, 0), 40000, Sigma, beta = 0.5))
  d <- out$draws[-(1:500), ]
  expect_within_mcse(mean(d[, 1]), m[1], mcse_mean(d[, 1]))
  expect_within_mcse(mean(d[, 2]), m[2], mcse_mean(d[, 2]))
})
