# Benchmarking helpers (Chapters 18 and 26, Project 6) ------------------------------

#' Run `expr` and return its value together with elapsed wall-clock seconds.
timed <- function(expr) {
  t0 <- proc.time()[["elapsed"]]
  value <- expr
  list(value = value, seconds = proc.time()[["elapsed"]] - t0)
}

#' Efficiency summary for one scalar quantity: bulk-ESS, ESS per second and
#' ESS per 1000 gradient (or density) evaluations. Comparing samplers only by
#' ESS per iteration ignores the cost of an iteration; comparing only by ESS
#' per second ties the result to a machine and an implementation.
efficiency_row <- function(draws, seconds, n_evals, label) {
  e <- ess_bulk(draws)
  data.frame(sampler = label, ess_bulk = e, seconds = seconds,
             ess_per_sec = e / seconds, ess_per_1000_evals = 1000 * e / n_evals)
}

#' Compact description of the computing environment for benchmark reports.
environment_report <- function() {
  si <- utils::sessionInfo()
  cpu <- tryCatch(system("sysctl -n machdep.cpu.brand_string", intern = TRUE, ignore.stderr = TRUE),
                  error = function(e) NA_character_, warning = function(w) NA_character_)
  if (length(cpu) == 0L) cpu <- NA_character_
  data.frame(
    item = c("R version", "platform", "operating system", "CPU", "BLAS", "RNG kind"),
    value = c(R.version.string, si$platform, si$running, cpu[1],
              if (is.null(si$BLAS)) NA_character_ else basename(si$BLAS),
              paste(RNGkind()[1:2], collapse = " / "))
  )
}
