# MCMC output analysis (Chapters 17-19) ---------------------------------------------
#
# Educational implementations of autocorrelation, effective sample size,
# (rank-normalised split) R-hat, Monte Carlo standard errors, Geweke's
# diagnostic and simulation-based calibration. They follow the definitions in
# Geyer (1992), Vehtari et al. (2021) and Talts et al. (2018) and are tested
# against the `posterior` and `coda` packages in tests/testthat.
#
# Convention: `draws` is an iterations x chains matrix for one scalar quantity
# (a vector is treated as a single chain).

as_draws_matrix1 <- function(draws) {
  draws <- as.matrix(draws)
  if (!is.numeric(draws)) course_stop("`draws` must be numeric.")
  draws
}

# ---- Autocorrelation -------------------------------------------------------------

#' Biased (divide-by-n) autocovariance at lags 0..n-1, computed with the FFT in
#' O(n log n). Zero padding to length >= 2n prevents circular wrap-around.
#' The biased estimator is used because it yields a positive semi-definite
#' sequence (Geyer, 1992).
autocovariance <- function(x) {
  n <- length(x)
  xc <- x - mean(x)
  m <- stats::nextn(2L * n)
  f <- stats::fft(c(xc, numeric(m - n)))
  acov <- Re(stats::fft(Mod(f)^2, inverse = TRUE))[seq_len(n)] / m
  acov / n
}

#' Autocorrelation at lags 0..n-1.
autocorrelation <- function(x) {
  g <- autocovariance(x)
  if (g[1] == 0) return(c(1, rep(NA_real_, length(x) - 1L)))
  g / g[1]
}

# ---- Effective sample size ----------------------------------------------------------

#' Effective sample size of one or several chains for the mean of a scalar quantity.
#'
#' Combines within-chain autocovariances with the between-chain variance
#' (Vehtari et al., 2021, Section 3.2):
#'   rho_t = 1 - (W - mean_m acov_{m,t}) / var_plus,
#' then sums the autocorrelations with Geyer's initial monotone sequence
#' estimator: pair sums P_k = rho_{2k} + rho_{2k+1} are accumulated while
#' positive and forced to be non-increasing. ESS = M N / tau with
#' tau = -1 + 2 sum_k P_k.
ess_chains <- function(draws) {
  draws <- as_draws_matrix1(draws)
  n <- nrow(draws); m <- ncol(draws)
  if (n < 4L || any(!is.finite(draws))) return(NA_real_)
  if (all(apply(draws, 2, function(z) max(z) - min(z)) == 0)) return(NA_real_)   # constant chains
  acov <- apply(draws, 2, autocovariance)
  acov <- matrix(acov, nrow = n)
  W <- mean(acov[1, ]) * n / (n - 1)          # mean within-chain variance (unbiased)
  var_plus <- W * (n - 1) / n
  if (m > 1L) var_plus <- var_plus + stats::var(colMeans(draws))
  rho <- 1 - (W - rowMeans(acov)) / var_plus
  rho[1] <- 1
  # Geyer's initial positive sequence on pairs of lags (0,1), (2,3), ...
  P <- numeric(0); k <- 0L
  while (2L * k + 2L <= n) {
    pk <- rho[2L * k + 1L] + rho[2L * k + 2L]
    if (k > 0L && pk <= 0) break
    P <- c(P, pk); k <- k + 1L
  }
  P <- cummin(P)                                  # initial monotone sequence
  tau <- -1 + 2 * sum(P)
  tau <- max(tau, 1 / log10(n * m))               # cap ESS at N log10(N), as in Stan
  n * m / tau
}

#' Integrated autocorrelation time tau = 1 + 2 sum_{t>=1} rho_t of a single chain.
iat <- function(x) length(x) / ess_chains(x)

#' Split each chain in half (dropping the middle draw if n is odd), doubling the
#' number of chains. Non-stationarity within a chain then shows up as
#' disagreement between its two halves.
split_chains <- function(draws) {
  draws <- as_draws_matrix1(draws)
  n <- nrow(draws); half <- n %/% 2L
  cbind(draws[seq_len(half), , drop = FALSE], draws[(n - half + 1L):n, , drop = FALSE])
}

#' Rank-normalisation: replace each draw by Phi^{-1}((r - 3/8) / (S + 1/4)),
#' where r is its rank among all S draws pooled across chains (ties averaged).
z_scale <- function(draws) {
  draws <- as_draws_matrix1(draws)
  S <- length(draws)
  r <- rank(as.vector(draws), ties.method = "average")
  matrix(stats::qnorm((r - 3 / 8) / (S + 1 / 4)), nrow = nrow(draws))
}

#' Classic potential scale reduction factor (Gelman and Rubin, 1992; BDA3 form).
#' R-hat = sqrt(var_plus / W), var_plus = (n - 1)/n W + B/n,
#' W = mean within-chain variance, B = n * variance of chain means.
rhat_basic <- function(draws, split = TRUE) {
  draws <- as_draws_matrix1(draws)
  if (split) draws <- split_chains(draws)
  n <- nrow(draws)
  W <- mean(apply(draws, 2, stats::var))
  B <- n * stats::var(colMeans(draws))
  if (W == 0) return(NA_real_)
  sqrt(((n - 1) / n * W + B / n) / W)
}

#' Rank-normalised split R-hat (Vehtari et al., 2021): the maximum of the bulk
#' R-hat (on rank-normalised draws) and the tail/scale R-hat (on rank-normalised
#' absolute deviations from the median).
rhat_rank <- function(draws) {
  draws <- as_draws_matrix1(draws)
  sp <- split_chains(draws)
  bulk <- rhat_basic(z_scale(sp), split = FALSE)
  folded <- abs(sp - stats::median(sp))
  tail <- rhat_basic(z_scale(folded), split = FALSE)
  max(bulk, tail)
}

#' Bulk-ESS: ESS of rank-normalised split chains (efficiency for the centre of
#' the distribution; well-defined even without finite moments).
ess_bulk <- function(draws) ess_chains(z_scale(split_chains(draws)))

#' Tail-ESS: the minimum ESS of the 5% and 95% quantile indicators.
ess_tail <- function(draws) {
  draws <- as_draws_matrix1(draws)
  sp <- split_chains(draws)
  q <- stats::quantile(sp, c(0.05, 0.95), names = FALSE)
  min(ess_chains((sp <= q[1]) * 1), ess_chains((sp <= q[2]) * 1))
}

#' ESS for the mean (no rank normalisation), on split chains.
ess_mean <- function(draws) ess_chains(split_chains(draws))

# ---- Monte Carlo standard errors --------------------------------------------------

#' MCSE of the posterior-mean estimate: sd / sqrt(ESS).
mcse_mean <- function(draws) {
  draws <- as_draws_matrix1(draws)
  stats::sd(as.vector(draws)) / sqrt(ess_mean(draws))
}

#' MCSE of a quantile estimate (Vehtari et al., 2021, Section 4.3): the ESS of
#' the indicator I(x <= q) is plugged into a Beta approximation of the
#' sampling distribution of the empirical CDF at q, and mapped back to the
#' draw scale through order statistics.
mcse_quantile <- function(draws, prob) {
  draws <- as_draws_matrix1(draws)
  x <- as.vector(draws); N <- length(x)
  q <- stats::quantile(x, prob, names = FALSE)
  S <- ess_chains(split_chains((draws <= q) * 1))
  a <- stats::qbeta(stats::pnorm(-1), S * prob + 1, S * (1 - prob) + 1)
  b <- stats::qbeta(stats::pnorm(1), S * prob + 1, S * (1 - prob) + 1)
  xs <- sort(x)
  (xs[min(ceiling(b * N), N)] - xs[max(floor(a * N), 1L)]) / 2
}

#' Batch-means estimate of the MCSE of the mean of a single chain.
#' Splits the chain into a batches of size b; the variance of batch means,
#' scaled by b, estimates the asymptotic variance sigma^2. Consistent when
#' b -> Inf and a -> Inf (b = floor(sqrt(n)) is the common default; see
#' Flegal and Jones, 2010).
mcse_batch_means <- function(x, batch_size = floor(sqrt(length(x)))) {
  n <- length(x); b <- check_count(batch_size, "batch_size", min = 1L)
  a <- n %/% b
  if (a < 2L) course_stop("Need at least two batches.")
  x <- x[seq_len(a * b)]
  bm <- colMeans(matrix(x, nrow = b))
  sigma2 <- b * stats::var(bm)
  sqrt(sigma2 / (a * b))
}

#' Spectral density at frequency zero via a fitted autoregression (as in
#' coda::spectrum0.ar): S(0) = sigma_e^2 / (1 - sum phi_k)^2. For a stationary
#' chain, S(0) equals the asymptotic variance sigma^2 in the MCMC CLT.
spectrum0_ar <- function(x) {
  fit <- stats::ar(x, aic = TRUE)
  fit$var.pred / (1 - sum(fit$ar))^2
}

#' Geweke's (1992) diagnostic: z-score comparing the mean of the first
#' `first` fraction with the mean of the last `last` fraction of one chain,
#' with variances from spectral densities at zero.
geweke_z <- function(x, first = 0.1, last = 0.5) {
  n <- length(x)
  # window boundaries follow coda::geweke.diag so the two implementations agree
  xa <- x[seq_len(ceiling(1 + first * (n - 1)))]
  xb <- x[floor(n - last * (n - 1)):n]
  (mean(xa) - mean(xb)) / sqrt(spectrum0_ar(xa) / length(xa) + spectrum0_ar(xb) / length(xb))
}

#' One-row summary of a scalar quantity from several chains.
summarise_draws_basic <- function(draws, name = "x") {
  draws <- as_draws_matrix1(draws)
  x <- as.vector(draws)
  data.frame(variable = name, mean = mean(x), sd = stats::sd(x),
             q5 = stats::quantile(x, 0.05, names = FALSE),
             q95 = stats::quantile(x, 0.95, names = FALSE),
             rhat = rhat_rank(draws), ess_bulk = ess_bulk(draws),
             ess_tail = ess_tail(draws), mcse_mean = mcse_mean(draws))
}

# ---- Simulation-based calibration --------------------------------------------------

#' Simulation-based calibration (Talts et al., 2018).
#'
#' For s = 1..n_sims: draw theta* ~ prior, y ~ p(y | theta*), L posterior draws
#' theta_1..L given y, and record the rank r = #{l : theta_l < theta*} in
#' {0, ..., L}. If the posterior sampler is correct, the ranks are uniform on
#' {0, ..., L} (exactly, when the L draws are independent).
#'
#' @param r_prior  function() returning one parameter value (scalar).
#' @param r_data   function(theta) returning a simulated data set.
#' @param posterior_draws function(y, L) returning L posterior draws of theta.
sbc_ranks <- function(n_sims, L, r_prior, r_data, posterior_draws) {
  n_sims <- check_count(n_sims, "n_sims"); L <- check_count(L, "L")
  vapply(seq_len(n_sims), function(s) {
    theta <- r_prior()
    y <- r_data(theta)
    draws <- posterior_draws(y, L)
    sum(draws < theta)
  }, numeric(1))
}

#' Chi-square test of uniformity of SBC ranks over `n_bins` equal-width bins.
#' A small p-value indicates miscalibration; a large one is not proof of
#' correctness (the test has limited power against some errors).
sbc_uniformity_test <- function(ranks, L, n_bins = 10) {
  bins <- cut(ranks, breaks = seq(-0.5, L + 0.5, length.out = n_bins + 1L))
  counts <- as.vector(table(bins))
  expected <- length(ranks) / n_bins
  stat <- sum((counts - expected)^2 / expected)
  list(counts = counts, statistic = stat,
       p_value = stats::pchisq(stat, df = n_bins - 1L, lower.tail = FALSE))
}

#' Histogram of SBC ranks with a 99% band for bin counts under uniformity.
sbc_plot <- function(ranks, L, n_bins = 20, main = "") {
  breaks <- seq(-0.5, L + 0.5, length.out = n_bins + 1L)
  h <- graphics::hist(ranks, breaks = breaks, plot = FALSE)
  N <- length(ranks)
  band <- stats::qbinom(c(0.005, 0.995), N, 1 / n_bins)
  graphics::plot(NULL, xlim = range(breaks), ylim = c(0, max(h$counts, band[2]) * 1.1),
                 xlab = "rank statistic", ylab = "count", main = main)
  graphics::rect(min(breaks), band[1], max(breaks), band[2], col = alpha_col(course_muted, 0.15), border = NA)
  graphics::abline(h = N / n_bins, col = course_muted, lty = 2)
  graphics::rect(head(breaks, -1), 0, tail(breaks, -1), h$counts, col = course_fill, border = "white")
  invisible(h$counts)
}
