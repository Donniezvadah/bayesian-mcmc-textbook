test_that("HPD from draws matches the exact HPD of a Gamma distribution", {
  # exact HPD [l, u] for Gamma(3, 1): dgamma(l) = dgamma(u) and coverage 0.9
  f <- function(l) { u <- uniroot(function(u) dgamma(u, 3) - dgamma(l, 3), c(2, 30))$root
                     pgamma(u, 3) - pgamma(l, 3) - 0.9 }
  l <- uniroot(f, c(0.05, 1.99))$root
  u <- uniroot(function(u) dgamma(u, 3) - dgamma(l, 3), c(2, 30))$root
  x <- with_seed(80, rgamma(2e5, 3))
  expect_equal(hpd_interval(x, 0.9), c(l, u), tolerance = 0.02)
  expect_lt(diff(hpd_interval(x, 0.9)), diff(credible_interval(x, 0.9)))   # shortest
  g <- seq(0, 20, length.out = 20001)
  hg <- hpd_grid(dgamma(g, 3), g, 0.9)
  expect_equal(unname(hg$intervals[1, ]), c(l, u), tolerance = 1e-3)
})

test_that("HPD region of a bimodal density has two pieces", {
  g <- seq(-8, 8, length.out = 8001)
  dens <- 0.5 * dnorm(g, -3) + 0.5 * dnorm(g, 3)
  expect_equal(nrow(hpd_grid(dens, g, 0.8)$intervals), 2L)
})

test_that("Bayes actions: squared error -> mean, absolute error -> median, pinball -> quantile", {
  x <- with_seed(81, rgamma(1e5, 2, 1))
  expect_equal(bayes_action(x, function(a, t) (a - t)^2), mean(x), tolerance = 1e-4)
  expect_equal(bayes_action(x, function(a, t) abs(a - t)), median(x), tolerance = 1e-3)
  tau <- 0.8
  pinball <- function(a, t) ifelse(t > a, tau * (t - a), (1 - tau) * (a - t))
  expect_equal(bayes_action(x, pinball), unname(quantile(x, tau)), tolerance = 1e-3)
})
