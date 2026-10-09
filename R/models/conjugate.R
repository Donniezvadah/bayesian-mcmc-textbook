# Conjugate Bayesian models -----------------------------------------------------
#
# Closed-form posterior updates derived in Chapter 5. Each function returns the
# posterior hyperparameters; predictive helpers return log probabilities or log
# densities. Parameterisations (fixed throughout the course):
#   Beta(a, b)                 density  x^(a-1) (1-x)^(b-1) / B(a, b)
#   Gamma(a, b)                shape a, RATE b: density b^a x^(a-1) e^(-b x) / Gamma(a)
#   InvGamma(a, b)             shape a, SCALE b: density b^a x^(-a-1) e^(-b/x) / Gamma(a)
#   N(m, v)                    mean m, VARIANCE v
#   Dirichlet(alpha)           concentration vector alpha
# These match R's dbeta(), dgamma(shape, rate), dnorm(mean, sd = sqrt(v)).

# ---- Beta-Binomial -----------------------------------------------------------

#' Posterior for theta in y | theta ~ Binomial(n, theta), theta ~ Beta(a, b).
#' Posterior: Beta(a + y, b + n - y).
beta_binomial_update <- function(y, n, a = 1, b = 1) {
  check_count(n, "n", min = 0L)
  if (!is.numeric(y) || length(y) != 1L || y < 0 || y > n || y != round(y)) {
    course_stop("`y` must be a whole number between 0 and n.")
  }
  check_scalar(a, "a", positive = TRUE); check_scalar(b, "b", positive = TRUE)
  list(a = a + y, b = b + n - y)
}

#' Log pmf of the Beta-Binomial distribution: the (prior or posterior)
#' predictive distribution of k successes in m future trials when
#' theta ~ Beta(a, b):  p(k) = choose(m, k) B(k + a, m - k + b) / B(a, b).
dbetabinom_log <- function(k, m, a, b) {
  out <- lchoose(m, k) + lbeta(k + a, m - k + b) - lbeta(a, b)
  out[k < 0 | k > m | k != round(k)] <- -Inf
  out
}

#' Random draws from the Beta-Binomial distribution by composition.
rbetabinom <- function(n_draws, m, a, b) {
  stats::rbinom(n_draws, size = m, prob = stats::rbeta(n_draws, a, b))
}

# ---- Gamma-Poisson -----------------------------------------------------------

#' y_i | lambda ~ Poisson(lambda * t_i), lambda ~ Gamma(a, b) (rate b).
#' Posterior: Gamma(a + sum y_i, b + sum t_i). `exposure` defaults to t_i = 1.
gamma_poisson_update <- function(y, a, b, exposure = rep(1, length(y))) {
  if (any(y < 0 | y != round(y))) course_stop("`y` must contain non-negative integers.")
  check_scalar(a, "a", positive = TRUE); check_scalar(b, "b", positive = TRUE)
  if (length(exposure) != length(y) || any(exposure <= 0)) {
    course_stop("`exposure` must be positive and have the same length as `y`.")
  }
  list(a = a + sum(y), b = b + sum(exposure))
}

#' Posterior predictive for a new count with exposure t under Gamma(a, b):
#' Negative-Binomial with size a and success probability b / (b + t).
dgamma_poisson_pred_log <- function(y_new, a, b, t = 1) {
  stats::dnbinom(y_new, size = a, prob = b / (b + t), log = TRUE)
}

# ---- Gamma-Exponential -------------------------------------------------------

#' y_i | lambda ~ Exponential(rate lambda), lambda ~ Gamma(a, b).
#' Posterior: Gamma(a + n, b + sum y_i).
gamma_exponential_update <- function(y, a, b) {
  if (any(!is.finite(y)) || any(y < 0)) course_stop("`y` must contain non-negative values.")
  check_scalar(a, "a", positive = TRUE); check_scalar(b, "b", positive = TRUE)
  list(a = a + length(y), b = b + sum(y))
}

#' Posterior predictive density of a new waiting time under Gamma(a, b):
#' p(y) = a b^a / (b + y)^(a + 1), y >= 0  (a Lomax / Pareto type II law).
dgamma_exponential_pred_log <- function(y_new, a, b) {
  out <- log(a) + a * log(b) - (a + 1) * log(b + y_new)
  out[y_new < 0] <- -Inf
  out
}

# ---- Normal, known variance ---------------------------------------------------

#' y_i | mu ~ N(mu, sigma2) with sigma2 known, mu ~ N(m0, v0).
#' Posterior precision is the sum of prior and data precisions; the posterior
#' mean is the precision-weighted average of m0 and ybar.
normal_known_var_update <- function(y, sigma2, m0, v0) {
  check_finite_vector(y, "y")
  check_scalar(sigma2, "sigma2", positive = TRUE)
  check_scalar(m0, "m0"); check_scalar(v0, "v0", positive = TRUE)
  n <- length(y)
  prec <- 1 / v0 + n / sigma2
  m <- (m0 / v0 + sum(y) / sigma2) / prec
  list(mean = m, var = 1 / prec)
}

# ---- Normal, unknown mean and variance (Normal-Inverse-Gamma) ----------------

#' y_i | mu, s2 ~ N(mu, s2);  mu | s2 ~ N(m0, s2 / k0);  s2 ~ InvGamma(a0, b0).
#' Returns the four posterior hyperparameters (m, k, a, b) of the
#' Normal-Inverse-Gamma posterior. Derivation in Chapter 5.
nig_update <- function(y, m0, k0, a0, b0) {
  check_finite_vector(y, "y")
  check_scalar(m0, "m0"); check_scalar(k0, "k0", positive = TRUE)
  check_scalar(a0, "a0", positive = TRUE); check_scalar(b0, "b0", positive = TRUE)
  n <- length(y); ybar <- mean(y)
  ss <- sum((y - ybar)^2)
  kn <- k0 + n
  list(m = (k0 * m0 + n * ybar) / kn,
       k = kn,
       a = a0 + n / 2,
       b = b0 + ss / 2 + k0 * n * (ybar - m0)^2 / (2 * kn))
}

#' Marginal posterior of mu under the NIG posterior: a location-scale Student t
#' with 2a degrees of freedom, location m and squared scale b / (a k).
nig_marginal_mu <- function(post) {
  list(df = 2 * post$a, location = post$m, scale = sqrt(post$b / (post$a * post$k)))
}

#' Posterior predictive of a new observation under the NIG posterior:
#' Student t with 2a degrees of freedom, location m, squared scale b (k + 1) / (a k).
nig_predictive <- function(post) {
  list(df = 2 * post$a, location = post$m,
       scale = sqrt(post$b * (post$k + 1) / (post$a * post$k)))
}

#' Draw (mu, s2) from a Normal-Inverse-Gamma distribution.
rnig <- function(n_draws, post) {
  s2 <- 1 / stats::rgamma(n_draws, shape = post$a, rate = post$b)
  mu <- stats::rnorm(n_draws, post$m, sqrt(s2 / post$k))
  cbind(mu = mu, sigma2 = s2)
}

# ---- Dirichlet-Multinomial -----------------------------------------------------

#' counts | p ~ Multinomial(N, p), p ~ Dirichlet(alpha). Posterior: Dirichlet(alpha + counts).
dirichlet_multinomial_update <- function(counts, alpha) {
  if (any(counts < 0 | counts != round(counts))) course_stop("`counts` must be non-negative integers.")
  if (length(alpha) != length(counts) || any(alpha <= 0)) {
    course_stop("`alpha` must be positive and have the same length as `counts`.")
  }
  alpha + counts
}

#' Random draws from Dirichlet(alpha) by normalising independent Gamma(alpha_k, 1) variables.
rdirichlet <- function(n_draws, alpha) {
  g <- matrix(stats::rgamma(n_draws * length(alpha), shape = alpha), ncol = length(alpha), byrow = TRUE)
  g / rowSums(g)
}

#' Log density of Dirichlet(alpha) at a point p of the simplex.
ddirichlet_log <- function(p, alpha) {
  lgamma(sum(alpha)) - sum(lgamma(alpha)) + sum((alpha - 1) * log(p))
}

#' Log pmf of the Dirichlet-Multinomial predictive for a future count vector x.
ddirmult_log <- function(x, alpha) {
  N <- sum(x); A <- sum(alpha)
  lgamma(N + 1) - sum(lgamma(x + 1)) + lgamma(A) - lgamma(N + A) +
    sum(lgamma(x + alpha) - lgamma(alpha))
}

# ---- Auxiliary densities --------------------------------------------------------

#' Inverse-gamma log density (shape a, scale b) and sampler.
dinvgamma_log <- function(x, a, b) {
  out <- a * log(b) - lgamma(a) - (a + 1) * log(x) - b / x
  out[x <= 0] <- -Inf
  out
}
rinvgamma <- function(n_draws, a, b) 1 / stats::rgamma(n_draws, shape = a, rate = b)

#' Location-scale Student t: density, distribution function, quantiles, draws.
dt_ls <- function(x, df, location = 0, scale = 1, log = FALSE) {
  out <- stats::dt((x - location) / scale, df = df, log = TRUE) - log(scale)
  if (log) out else exp(out)
}
pt_ls <- function(q, df, location = 0, scale = 1) stats::pt((q - location) / scale, df = df)
qt_ls <- function(p, df, location = 0, scale = 1) location + scale * stats::qt(p, df = df)
rt_ls <- function(n_draws, df, location = 0, scale = 1) location + scale * stats::rt(n_draws, df = df)
