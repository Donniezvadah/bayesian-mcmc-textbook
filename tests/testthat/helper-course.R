# Loaded automatically by testthat before the tests run: sources all course code.
course_root <- normalizePath(file.path(testthat::test_path(), "..", ".."))
for (f in sort(list.files(file.path(course_root, "R"), pattern = "[.]R$",
                          recursive = TRUE, full.names = TRUE))) {
  sys.source(f, envir = globalenv())
}

#' Expect a Monte Carlo estimate to agree with the truth within k Monte Carlo
#' standard errors. With k = 4 a correct sampler fails with probability about
#' 6e-5 per check under the CLT approximation, so failures signal real bugs.
expect_within_mcse <- function(estimate, truth, mcse, k = 4) {
  testthat::expect_true(
    abs(estimate - truth) <= k * mcse,
    label = sprintf("|%.5g - %.5g| = %.3g <= %g x MCSE (%.3g)",
                    estimate, truth, abs(estimate - truth), k, mcse)
  )
}
