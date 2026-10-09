# Hamiltonian Monte Carlo (Chapter 15) ----------------------------------------------
#
# Educational HMC with a diagonal mass matrix M = diag(1 / inv_mass), static
# trajectory length, and optional step-size jitter. Production samplers
# (Stan's NUTS) additionally adapt the step size and mass matrix and choose the
# trajectory length dynamically (Chapter 16).
#
# Hamiltonian:  H(q, p) = U(q) + K(p),  U(q) = -log pi(q),  K(p) = p' M^{-1} p / 2.

#' Leapfrog (Stoermer-Verlet) integration of Hamilton's equations
#'   dq/dt = M^{-1} p,   dp/dt = -grad U(q) = grad log pi(q).
#' The map is volume-preserving and time-reversible (with momentum flip);
#' its global energy error is O(step_size^2) on a fixed time interval.
#'
#' @return list(q, p) at the end of the trajectory, or with `path = TRUE` the
#'   whole sequence of positions (useful for plotting trajectories).
leapfrog <- function(q, p, grad_log_target, step_size, n_steps, inv_mass = 1, path = FALSE) {
  qs <- if (path) matrix(NA_real_, n_steps + 1L, length(q)) else NULL
  if (path) qs[1, ] <- q
  p <- p + 0.5 * step_size * grad_log_target(q)          # half step for momentum
  for (l in seq_len(n_steps)) {
    q <- q + step_size * inv_mass * p                    # full step for position
    g <- grad_log_target(q)
    if (any(!is.finite(g))) return(list(q = q, p = p, diverged = TRUE, path = qs))
    if (l < n_steps) p <- p + step_size * g              # full step for momentum
    if (path) qs[l + 1L, ] <- q
  }
  p <- p + 0.5 * step_size * g                           # final half step
  list(q = q, p = p, diverged = FALSE, path = qs)
}

#' Hamiltonian Monte Carlo with fixed integration time.
#'
#' @param log_target, grad_log_target  log pi(q) up to a constant, and its gradient.
#' @param step_size  leapfrog step size epsilon (> 0).
#' @param n_steps    number of leapfrog steps L per iteration.
#' @param inv_mass   diagonal of M^{-1} (scalar or length-d vector); ideally the
#'                   posterior variances, which makes the problem better conditioned.
#' @param jitter     if TRUE, draw each iteration's step size uniformly from
#'                   [0.8, 1.2] x step_size to avoid periodic trajectories.
#' @param divergence_threshold  an energy error larger than this marks the
#'                   transition as divergent (Stan uses 1000).
#' @return list with draws, acceptance_rate, energy_error (H_new - H_old per
#'   iteration), divergent flags, and mean acceptance probability.
#' Complexity: n_iter * n_steps gradient evaluations.
hmc_sampler <- function(log_target, grad_log_target, initial, n_iter, step_size,
                        n_steps, inv_mass = 1, jitter = TRUE, divergence_threshold = 1000) {
  check_function(log_target, "log_target"); check_function(grad_log_target, "grad_log_target")
  check_finite_vector(initial, "initial")
  n_iter <- check_count(n_iter, "n_iter"); n_steps <- check_count(n_steps, "n_steps")
  check_scalar(step_size, "step_size", positive = TRUE)
  d <- length(initial)
  inv_mass <- rep_len(inv_mass, d)
  if (any(!is.finite(inv_mass)) || any(inv_mass <= 0)) course_stop("`inv_mass` must be positive.")

  q <- initial
  logp <- eval_log_density(log_target, q)
  if (!is.finite(logp)) course_stop("Initial value must have finite log target density.")

  draws <- matrix(NA_real_, n_iter, d, dimnames = list(NULL, names(initial)))
  energy_error <- numeric(n_iter); divergent <- logical(n_iter); accept_prob <- numeric(n_iter)
  accepted <- 0L
  for (i in seq_len(n_iter)) {
    p0 <- stats::rnorm(d, 0, sqrt(1 / inv_mass))         # p ~ N(0, M)
    eps <- if (jitter) step_size * stats::runif(1, 0.8, 1.2) else step_size
    traj <- leapfrog(q, p0, grad_log_target, eps, n_steps, inv_mass)
    logp_new <- if (traj$diverged) -Inf else eval_log_density(log_target, traj$q)
    H0 <- -logp + 0.5 * sum(inv_mass * p0^2)
    H1 <- -logp_new + 0.5 * sum(inv_mass * traj$p^2)
    dH <- H1 - H0                                          # +Inf if the trajectory left the support
    energy_error[i] <- dH
    divergent[i] <- !is.finite(dH) || dH > divergence_threshold
    accept_prob[i] <- if (is.finite(dH)) min(1, exp(-dH)) else 0
    # Metropolis correction: accept with probability min(1, exp(H0 - H1)).
    if (is.finite(dH) && log(stats::runif(1)) < -dH) {
      q <- traj$q; logp <- logp_new; accepted <- accepted + 1L
    }
    draws[i, ] <- q
  }
  list(draws = draws, acceptance_rate = accepted / n_iter,
       mean_accept_prob = mean(accept_prob), energy_error = energy_error,
       divergent = divergent,
       settings = list(step_size = step_size, n_steps = n_steps, inv_mass = inv_mass))
}

#' Dual-averaging adaptation of the HMC step size (Nesterov, 2009; the variant of
#' Hoffman and Gelman, 2014, Section 3.2). During warm-up the log step size is
#' moved so that the average acceptance probability approaches `target_accept`.
#' The adapted step size is then frozen; sampling after warm-up uses a fixed
#' kernel, which is what guarantees invariance of the retained draws.
hmc_warmup <- function(log_target, grad_log_target, initial, n_warmup, n_steps,
                       inv_mass = 1, step_size = 0.1, target_accept = 0.8,
                       gamma = 0.05, t0 = 10, kappa = 0.75) {
  n_warmup <- check_count(n_warmup, "n_warmup")
  mu <- log(10 * step_size); log_eps_bar <- 0; h_bar <- 0
  q <- initial; eps <- step_size
  for (m in seq_len(n_warmup)) {
    out <- hmc_sampler(log_target, grad_log_target, q, 1L, eps, n_steps, inv_mass, jitter = FALSE)
    q <- out$draws[1, ]
    a <- out$mean_accept_prob
    w <- 1 / (m + t0)
    h_bar <- (1 - w) * h_bar + w * (target_accept - a)
    log_eps <- mu - sqrt(m) / gamma * h_bar
    eta <- m^(-kappa)
    log_eps_bar <- eta * log_eps + (1 - eta) * log_eps_bar
    eps <- exp(log_eps)
  }
  list(step_size = exp(log_eps_bar), position = q)
}
