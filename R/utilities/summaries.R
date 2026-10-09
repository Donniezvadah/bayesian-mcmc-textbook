# Posterior summaries (Chapter 6) ----------------------------------------------------

#' Central (equal-tailed) credible interval from draws.
credible_interval <- function(draws, prob = 0.95) {
  check_probability(prob, "prob")
  a <- (1 - prob) / 2
  stats::quantile(draws, c(a, 1 - a), names = FALSE)
}

#' Shortest interval containing a fraction `prob` of the draws (Chen and Shao, 1999).
#'
#' Among the intervals [x_(i), x_(i+m)] between order statistics that contain
#' m + 1 = ceiling(prob * S) draws, return the shortest. For a unimodal
#' posterior this estimates the highest posterior density (HPD) interval
#' consistently as S -> Inf. For multimodal posteriors the HPD region may be a
#' union of intervals, which this function cannot represent.
hpd_interval <- function(draws, prob = 0.95) {
  check_probability(prob, "prob")
  x <- sort(draws); S <- length(x)
  m <- ceiling(prob * S) - 1L
  if (m < 1L || m >= S) course_stop("Too few draws for the requested probability.")
  widths <- x[(m + 1L):S] - x[1:(S - m)]
  i <- which.min(widths)
  c(x[i], x[i + m])
}

#' HPD region of a one-dimensional density evaluated on an equally spaced grid.
#' Returns the density threshold k and the grid points whose density exceeds it;
#' the region {theta : p(theta) >= k} may consist of several intervals.
hpd_grid <- function(density, grid, prob = 0.95) {
  h <- diff(grid[1:2])
  o <- order(density, decreasing = TRUE)
  cum <- cumsum(density[o]) * h
  k <- density[o][which(cum >= prob)[1]]
  inside <- density >= k
  # identify contiguous runs of grid points inside the region
  runs <- rle(inside); ends <- cumsum(runs$lengths); starts <- ends - runs$lengths + 1
  intervals <- cbind(lower = grid[starts[runs$values]], upper = grid[ends[runs$values]])
  list(threshold = k, intervals = intervals, coverage = sum(density[inside]) * h)
}

#' Bayes action under a loss function, minimising Monte Carlo posterior expected loss
#' (1/S) sum_s L(a, theta_s) over a one-dimensional action a.
bayes_action <- function(draws, loss, interval = range(draws)) {
  check_function(loss, "loss")
  stats::optimize(function(a) mean(loss(a, draws)), interval)$minimum
}
