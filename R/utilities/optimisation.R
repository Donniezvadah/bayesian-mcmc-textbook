# Optimisation and Laplace approximation (Chapter 2) ----------------------------

#' Maximise a smooth function with Newton's method and backtracking line search.
#'
#' Iteration: x_new = x - step * H(x)^{-1} g(x), with step halved until the
#' objective increases (Armijo-type safeguard). Near a non-degenerate maximum,
#' where the Hessian H is negative definite, convergence is quadratic.
#'
#' @param f, grad, hess objective, gradient and Hessian functions of x.
#' @param x0 starting point. @param tol stop when max |gradient| < tol.
#' @return list(par, value, gradient, hessian, iterations, converged)
newton_maximise <- function(f, grad, hess, x0, tol = 1e-8, max_iter = 100L) {
  check_function(f, "f"); check_function(grad, "grad"); check_function(hess, "hess")
  check_finite_vector(x0, "x0")
  x <- x0; fx <- f(x)
  for (it in seq_len(max_iter)) {
    g <- grad(x)
    if (max(abs(g)) < tol) {
      return(list(par = x, value = fx, gradient = g, hessian = hess(x),
                  iterations = it - 1L, converged = TRUE))
    }
    H <- hess(x)
    # Newton direction solves H d = -g; if H is not negative definite, fall back
    # to a gradient-ascent direction so that the step is still uphill.
    d <- tryCatch(-solve(H, g), error = function(e) g)
    if (sum(d * g) <= 0) d <- g
    step <- 1
    repeat {
      x_new <- x + step * d
      f_new <- f(x_new)
      if (is.finite(f_new) && f_new >= fx + 1e-4 * step * sum(d * g)) break
      step <- step / 2
      if (step < 1e-12) break
    }
    x <- x_new; fx <- f_new
  }
  list(par = x, value = fx, gradient = grad(x), hessian = hess(x),
       iterations = max_iter, converged = max(abs(grad(x))) < tol)
}

#' Laplace approximation to a posterior known up to a constant.
#'
#' With h(theta) = log of the unnormalised posterior, mode theta_hat and
#' A = -Hessian of h at theta_hat (positive definite):
#'   posterior  ~  N(theta_hat, A^{-1}),
#'   log Z     ~=  h(theta_hat) + (d/2) log(2 pi) - (1/2) log det A.
#' Exact when h is quadratic (Gaussian posterior); otherwise an asymptotic
#' approximation (Tierney and Kadane, 1986).
laplace_approximation <- function(log_post, grad, hess, x0, ...) {
  opt <- newton_maximise(log_post, grad, hess, x0, ...)
  if (!opt$converged) warning("Newton's method did not converge; the Laplace approximation is unreliable.")
  A <- -opt$hessian
  U <- tryCatch(chol(A), error = function(e) course_stop("Negative Hessian is not positive definite at the mode."))
  d <- length(x0)
  list(mode = opt$par, cov = chol2inv(U),
       log_Z = opt$value + 0.5 * d * log(2 * pi) - sum(log(diag(U))),
       iterations = opt$iterations)
}
