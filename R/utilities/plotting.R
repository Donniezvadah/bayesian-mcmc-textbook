# Base-graphics plotting helpers -----------------------------------------------
#
# Figures in the book use base R graphics so that every plot is reproducible
# with no extra packages. Colors follow a fixed categorical order (validated for
# colour-vision deficiency); series identity is always also given by a legend or
# a direct label, never by colour alone.

course_cols <- c(
  blue    = "#2a78d6",
  orange  = "#eb6834",
  aqua    = "#1baf7a",
  yellow  = "#eda100",
  magenta = "#e87ba4",
  green   = "#008300",
  violet  = "#4a3aa7",
  red     = "#e34948"
)
course_ink   <- "#1d1d1f"   # primary text / reference marks
course_muted <- "#6b6a66"   # secondary ink: grids, reference lines, annotations
course_fill  <- "#cde2fb"   # light blue fill for histograms and bands

#' Semi-transparent version of a colour (for overlapping chains / bands).
alpha_col <- function(col, alpha = 0.4) {
  rgb_vals <- grDevices::col2rgb(col) / 255
  grDevices::rgb(rgb_vals[1, ], rgb_vals[2, ], rgb_vals[3, ], alpha = alpha)
}

#' Course default graphical parameters: thin axes, horizontal labels, compact margins.
course_par <- function(...) {
  graphics::par(mar = c(4, 4.2, 1.6, 1), mgp = c(2.5, 0.7, 0), las = 1,
                bty = "l", tcl = -0.3, cex.axis = 0.85, cex.lab = 0.95,
                cex.main = 1, font.main = 1, col.axis = course_ink,
                col.lab = course_ink, fg = course_muted, ...)
}

#' Trace plot for one or several chains.
#'
#' @param draws A numeric vector (one chain) or an iterations x chains matrix.
#' @param ylab  Axis label.
#' @param truth Optional reference value drawn as a dashed horizontal line.
trace_plot <- function(draws, ylab = "value", truth = NULL, main = "", ...) {
  draws <- as.matrix(draws)
  n_chain <- ncol(draws)
  cols <- if (n_chain == 1L) course_cols[["blue"]] else
    alpha_col(course_cols[seq_len(min(n_chain, 8L))], 0.75)
  graphics::matplot(draws, type = "l", lty = 1, lwd = 0.7, col = cols,
                    xlab = "iteration", ylab = ylab, main = main, ...)
  if (!is.null(truth)) graphics::abline(h = truth, lty = 2, col = course_ink, lwd = 1)
  if (n_chain > 1L) {
    graphics::legend("topright", legend = paste("chain", seq_len(n_chain)),
                     col = cols, lty = 1, lwd = 1.5, bty = "n", cex = 0.75,
                     ncol = min(n_chain, 4L))
  }
  invisible(NULL)
}

#' Empirical autocorrelation function as thin vertical bars.
acf_plot <- function(x, max_lag = 50, main = "", col = course_cols[["blue"]], add_bounds = FALSE) {
  rho <- autocorrelation(x)[seq_len(min(max_lag, length(x) - 1L) + 1L)]
  lags <- seq_along(rho) - 1L
  graphics::plot(lags, rho, type = "h", lwd = 2, col = col, ylim = c(min(0, min(rho)), 1),
                 xlab = "lag", ylab = "autocorrelation", main = main)
  graphics::abline(h = 0, col = course_muted)
  if (add_bounds) {
    # Pointwise +/- 2/sqrt(n) bands are valid only under an i.i.d. null.
    graphics::abline(h = c(-2, 2) / sqrt(length(x)), lty = 3, col = course_muted)
  }
  invisible(rho)
}

#' Histogram of draws with an optional exact density overlaid.
hist_density <- function(draws, density_fun = NULL, breaks = 50, xlab = "value",
                         main = "", xlim = NULL, ...) {
  h <- graphics::hist(draws, breaks = breaks, plot = FALSE)
  if (is.null(xlim)) xlim <- range(h$breaks)
  ymax <- max(h$density)
  if (!is.null(density_fun)) {
    grid <- seq(xlim[1], xlim[2], length.out = 501)
    dens <- suppressWarnings(density_fun(grid))
    ymax <- max(ymax, min(max(dens[is.finite(dens)]), 3 * ymax))
  }
  graphics::plot(h, freq = FALSE, col = course_fill, border = "white", xlab = xlab,
                 main = main, xlim = xlim, ylim = c(0, 1.05 * ymax), ...)
  if (!is.null(density_fun)) {
    graphics::curve(density_fun(x), add = TRUE, lwd = 2, col = course_cols[["orange"]], n = 501)
  }
  invisible(NULL)
}

#' Rank plot (Vehtari et al., 2021): histograms of the ranks of each chain's
#' draws among all draws pooled. If all chains sample the same distribution,
#' every histogram is approximately uniform; a chain stuck in one region
#' produces a histogram concentrated at low or high ranks.
rank_plot <- function(draws, bins = 20, main = "") {
  draws <- as.matrix(draws)
  r <- matrix(rank(as.vector(draws), ties.method = "average"), nrow = nrow(draws))
  M <- ncol(draws); S <- length(draws)
  breaks <- seq(0, S, length.out = bins + 1)
  counts <- apply(r, 2, function(rk) graphics::hist(rk, breaks = breaks, plot = FALSE)$counts)
  expected <- nrow(draws) / bins
  old <- graphics::par(mfrow = c(1, M), mar = c(2, 2, 2, 0.5)); on.exit(graphics::par(old))
  for (m in seq_len(M)) {
    graphics::barplot(counts[, m], space = 0, col = course_fill, border = "white", main = paste("chain", m),
                      ylim = c(0, max(counts) * 1.1), axes = m == 1, cex.main = 0.85)
    graphics::abline(h = expected, lty = 2, col = course_muted)
  }
  invisible(counts)
}
