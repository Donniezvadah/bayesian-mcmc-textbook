test_that("metropolis_hastings validates its inputs", {
  lt <- function(x) dnorm(x, log = TRUE)
  expect_error(metropolis_hastings("not a function", 0, 10, 1), "function")
  expect_error(metropolis_hastings(lt, c(0, 1), 10, 1), "single finite")
  expect_error(metropolis_hastings(lt, NA_real_, 10, 1), "single finite")
  expect_error(metropolis_hastings(lt, 0, 0, 1), "whole number")
  expect_error(metropolis_hastings(lt, 0, 10.5, 1), "whole number")
  expect_error(metropolis_hastings(lt, 0, 10, -1), "strictly positive")
  # initial value outside the support of a Gamma target
  expect_error(metropolis_hastings(target_gamma()$log_density, -1, 10, 1), "finite log target")
  # a target returning +Inf is a bug and must stop the sampler
  # (the chain is confined to [-1, 1], so a proposal outside is certain to occur)
  expect_error(with_seed(1, metropolis_hastings(function(x) if (abs(x) > 1) Inf else 0, 0, 1000, 5)), "\\+Inf")
})

test_that("random-walk Metropolis recovers N(3, 2^2) moments within MC error", {
  tg <- target_normal(3, 2)
  chains <- with_seed(10, sapply(c(-5, 0, 5, 10), function(x0)
    metropolis_hastings(tg$log_density, x0, 6000, proposal_sd = 4.8)$draws[-(1:1000)]))
  expect_lt(rhat_rank(chains), 1.01)
  expect_within_mcse(mean(chains), 3, mcse_mean(chains))
  x2 <- (chains - 3)^2
  expect_within_mcse(mean(x2), 4, mcse_mean(x2))
  q <- quantile(chains, 0.9, names = FALSE)
  expect_within_mcse(q, qnorm(0.9, 3, 2), mcse_quantile(chains, 0.9))
})

test_that("non-finite proposals (outside the support) are rejected, never accepted", {
  tg <- target_gamma(shape = 2, rate = 1)
  out <- with_seed(11, metropolis_hastings(tg$log_density, 1, 20000, proposal_sd = 2))
  expect_true(all(out$draws > 0))
  # NaN-producing target (log of negative number) is treated as zero density
  lt_nan <- function(x) suppressWarnings(log(x) - x)        # Gamma(2, 1) kernel, NaN for x < 0
  out2 <- with_seed(12, metropolis_hastings(lt_nan, 1, 20000, proposal_sd = 2))
  expect_true(all(out2$draws > 0))
  d <- out2$draws[-(1:1000)]
  expect_within_mcse(mean(d), 2, mcse_mean(d))
})

test_that("acceptance rate responds to the proposal scale as theory predicts", {
  lt <- target_normal()$log_density
  acc <- with_seed(13, sapply(c(0.1, 2.4, 25), function(s)
    metropolis_hastings(lt, 0, 5000, s)$acceptance_rate))
  expect_true(acc[1] > 0.9 && acc[3] < 0.15)
  expect_true(acc[2] > 0.3 && acc[2] < 0.6)
})

test_that("multivariate RWM targets a correlated Gaussian", {
  tg <- target_correlated_gaussian(0.8)
  out <- with_seed(14, rw_metropolis(tg$log_density, c(0, 0), 30000,
                                     proposal_cov = 2.38^2 / 2 * tg$cov))
  d <- out$draws[-(1:2000), ]
  expect_within_mcse(mean(d[, 1]), 0, mcse_mean(d[, 1]))
  expect_within_mcse(mean(d[, 1] * d[, 2]), 0.8, mcse_mean(d[, 1] * d[, 2]))
  expect_error(rw_metropolis(tg$log_density, c(0, 0), 10, proposal_cov = diag(3)), "2 x 2")
})

test_that("an asymmetric proposal needs the Hastings correction", {
  # Target Exponential(1) on (0, Inf); multiplicative log-normal random walk
  # x' = x exp(s Z) has q(x' | x) = LN(x' | log x, s^2), which is asymmetric.
  lt <- function(x) if (x <= 0) -Inf else -x
  r_prop <- function(x) x * exp(0.8 * rnorm(1))
  log_q <- function(to, from) dlnorm(to, log(from), 0.8, log = TRUE)
  right <- with_seed(15, mh_sampler(lt, 1, 40000, r_prop, log_q))$draws[-(1:1000)]
  wrong <- with_seed(15, mh_sampler(lt, 1, 40000, r_prop, function(to, from) 0))$draws[-(1:1000)]
  expect_within_mcse(mean(right), 1, mcse_mean(right))
  # Ignoring the correction makes the chain reversible with respect to
  # exp(-x) / x instead of exp(-x) (because x' q(x' | x) = x q(x | x')). That
  # function is not integrable at 0, so the chain drifts towards 0 and the
  # sample mean is badly biased.
  expect_gt(abs(mean(wrong) - 1), 10 * mcse_mean(right))
})

test_that("independence sampler with heavy-tailed proposal recovers N(0, 1)", {
  out <- with_seed(16, independence_sampler(
    target_normal()$log_density, 0, 20000,
    r_g = function(n) rt(n, df = 3), log_g = function(x) dt(x, df = 3, log = TRUE)))
  expect_within_mcse(mean(out$draws^2), 1, mcse_mean(out$draws^2))
})

test_that("MH transition matrix satisfies detailed balance exactly", {
  pi_u <- c(1, 4, 2, 0.5, 3)                  # unnormalised target on 5 states
  pi_n <- pi_u / sum(pi_u)
  K <- length(pi_u)
  Q <- matrix(0, K, K)                          # asymmetric proposal on a cycle
  for (i in 1:K) { Q[i, i %% K + 1] <- 0.7; Q[i, (i - 2) %% K + 1] <- 0.3 }
  P <- mh_transition_matrix(pi_u, Q)
  expect_equal(rowSums(P), rep(1, K))
  flux <- pi_n * P                              # flux[i, j] = pi_i P_ij
  expect_equal(flux, t(flux))                   # detailed balance
  expect_equal(as.vector(pi_n %*% P), pi_n)     # invariance
  # convergence of the distribution from a point mass
  mu <- c(1, 0, 0, 0, 0)
  for (k in 1:200) mu <- mu %*% P
  expect_equal(as.vector(mu), pi_n, tolerance = 1e-8)
})
