# Run the course test suite from the project root:
#   Rscript tests/run-tests.R
# Exits with a non-zero status if any test fails (used by CI).
if (!requireNamespace("testthat", quietly = TRUE)) stop("Install the 'testthat' package first.")
Sys.setenv(TESTTHAT_EDITION = "3")   # modern expectations (waldo comparisons)
results <- testthat::test_dir("tests/testthat", reporter = testthat::ProgressReporter$new(show_praise = FALSE),
                              stop_on_failure = TRUE)
invisible(results)
