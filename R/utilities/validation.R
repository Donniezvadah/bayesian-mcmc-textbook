# Input validation helpers ----------------------------------------------------
#
# Small, explicit checks used by every sampler and model function in the course.
# They fail early with an informative message instead of letting an invalid
# argument propagate into a silent numerical error several iterations later.

#' Stop with a formatted message, without the call in the error text.
course_stop <- function(...) stop(sprintf(...), call. = FALSE)

#' Check that `x` is a function.
check_function <- function(x, name = deparse(substitute(x))) {
  if (!is.function(x)) course_stop("`%s` must be a function.", name)
  invisible(x)
}

#' Check that `x` is a single finite number, optionally strictly positive.
check_scalar <- function(x, name = deparse(substitute(x)), positive = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x)) {
    course_stop("`%s` must be a single finite number.", name)
  }
  if (positive && x <= 0) course_stop("`%s` must be strictly positive.", name)
  invisible(x)
}

#' Check that `x` is a single positive whole number (an iteration count, a
#' sample size, a number of leapfrog steps, ...).
check_count <- function(x, name = deparse(substitute(x)), min = 1L) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x != round(x) || x < min) {
    course_stop("`%s` must be a whole number >= %d.", name, as.integer(min))
  }
  invisible(as.integer(x))
}

#' Check that `x` is a numeric vector of finite values (optionally of a given length).
check_finite_vector <- function(x, name = deparse(substitute(x)), len = NULL) {
  if (!is.numeric(x) || length(x) == 0L || any(!is.finite(x))) {
    course_stop("`%s` must be a non-empty numeric vector of finite values.", name)
  }
  if (!is.null(len) && length(x) != len) {
    course_stop("`%s` must have length %d (it has length %d).", name, len, length(x))
  }
  invisible(x)
}

#' Check that `x` is a number in [0, 1] (or a vector of such numbers).
check_probability <- function(x, name = deparse(substitute(x))) {
  if (!is.numeric(x) || any(!is.finite(x)) || any(x < 0 | x > 1)) {
    course_stop("`%s` must contain probabilities in [0, 1].", name)
  }
  invisible(x)
}

#' Check that `S` is a symmetric positive-definite matrix of dimension `d`.
#' Positive definiteness is checked by attempting a Cholesky factorisation,
#' which is also what the samplers use downstream.
check_spd_matrix <- function(S, name = deparse(substitute(S)), d = NULL) {
  if (!is.matrix(S) || !is.numeric(S) || nrow(S) != ncol(S) || any(!is.finite(S))) {
    course_stop("`%s` must be a finite square numeric matrix.", name)
  }
  if (!is.null(d) && nrow(S) != d) {
    course_stop("`%s` must be %d x %d.", name, d, d)
  }
  if (max(abs(S - t(S))) > 1e-10 * max(1, max(abs(S)))) {
    course_stop("`%s` must be symmetric.", name)
  }
  ok <- tryCatch({ chol(S); TRUE }, error = function(e) FALSE)
  if (!ok) course_stop("`%s` must be positive definite.", name)
  invisible(S)
}

#' Evaluate a user-supplied log density and normalise its return value.
#'
#' A log density may legitimately return -Inf (zero density, e.g. outside the
#' support). NaN/NA is treated as -Inf so that a proposal outside the domain of a
#' formula (for instance log of a negative number) is rejected rather than
#' crashing the sampler. +Inf or a non-scalar return value indicates a bug in
#' the target and stops execution.
eval_log_density <- function(log_density, x) {
  value <- log_density(x)
  if (!is.numeric(value) || length(value) != 1L) {
    course_stop("The log density must return a single numeric value.")
  }
  if (is.na(value)) return(-Inf)   # covers both NaN and NA_real_
  if (value == Inf) {
    course_stop("The log density returned +Inf; the target is not a valid density at this point.")
  }
  value
}
