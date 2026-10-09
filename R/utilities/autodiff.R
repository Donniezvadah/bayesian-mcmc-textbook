# Forward-mode automatic differentiation with dual numbers (Chapter 27) ---------------
#
# A dual number a + b*delta (with delta^2 = 0) carries a value `v` and a directional
# derivative `d`. Every elementary operation propagates both by the chain rule, so
# evaluating f on dual(x, direction) returns f(x) and the directional derivative
# grad f(x) . direction, exact up to floating-point rounding. This is a teaching
# implementation: production systems (Stan, JAX, PyTorch) use reverse mode, which
# returns the whole gradient in one backward pass.
#
# Values may be vectors (v and d have the same length), so vectorised R code such as
# sum(y * eta - log1p(exp(eta))) works unchanged when eta is a dual vector.

#' Construct a dual number (or vector of dual numbers).
dual <- function(value, deriv = 0 * value) {
  if (!is.numeric(value) || !is.numeric(deriv)) course_stop("Dual numbers need numeric parts.")
  structure(list(v = value, d = rep_len(deriv, length(value))), class = "dual")
}

as_dual <- function(x) if (inherits(x, "dual")) x else dual(x, 0 * x)

print.dual <- function(x, ...) {
  cat("<dual>\n"); print(cbind(value = x$v, deriv = x$d), ...); invisible(x)
}

"[.dual" <- function(x, i) dual(x$v[i], x$d[i])

Ops.dual <- function(e1, e2) {
  if (missing(e2)) {
    a <- as_dual(e1)
    return(switch(.Generic, "-" = dual(-a$v, -a$d), "+" = a,
                  course_stop("Unary '%s' is not defined for dual numbers.", .Generic)))
  }
  a <- as_dual(e1); b <- as_dual(e2)
  switch(.Generic,
    "+" = dual(a$v + b$v, a$d + b$d),
    "-" = dual(a$v - b$v, a$d - b$d),
    "*" = dual(a$v * b$v, a$d * b$v + a$v * b$d),
    "/" = dual(a$v / b$v, (a$d * b$v - a$v * b$d) / b$v^2),
    "^" = {
      val <- a$v^b$v
      der <- b$v * a$v^(b$v - 1) * a$d
      if (any(b$d != 0)) der <- der + val * log(a$v) * b$d
      dual(val, der)
    },
    course_stop("Operator '%s' is not defined for dual numbers.", .Generic))
}

Math.dual <- function(x, ...) {
  v <- x$v; d <- x$d
  switch(.Generic,
    exp   = { e <- exp(v); dual(e, e * d) },
    log   = dual(log(v), d / v),
    sqrt  = { s <- sqrt(v); dual(s, d / (2 * s)) },
    log1p = dual(log1p(v), d / (1 + v)),
    expm1 = dual(expm1(v), exp(v) * d),
    sin   = dual(sin(v), cos(v) * d),
    cos   = dual(cos(v), -sin(v) * d),
    course_stop("Function '%s' is not implemented for dual numbers.", .Generic))
}

Summary.dual <- function(..., na.rm = FALSE) {
  args <- list(...)
  if (.Generic != "sum" || length(args) != 1L) course_stop("Only sum() of a single dual is implemented.")
  dual(sum(args[[1]]$v), sum(args[[1]]$d))
}

#' Gradient of a scalar function f: R^d -> R by forward-mode AD: d passes, each
#' seeding the direction with a unit vector e_k. `f` must use only operations
#' implemented above (arithmetic, ^, exp, log, sqrt, log1p, expm1, sin, cos, sum, [).
ad_gradient <- function(f, x) {
  check_function(f, "f"); check_finite_vector(x, "x")
  vapply(seq_along(x), function(k) {
    out <- f(dual(x, as.numeric(seq_along(x) == k)))
    if (!inherits(out, "dual") || length(out$v) != 1L) course_stop("`f` must return a scalar dual.")
    out$d
  }, numeric(1))
}
