# Fit the Stan versions of the course's example models with cmdstanr.
#
# Usage (from the project root): Rscript scripts/run-stan-examples.R
# Requires the cmdstanr package and a working CmdStan installation
# (see https://mc-stan.org/cmdstanr/). These fits are NOT run when the book is
# rendered, so that the site builds without a C++ toolchain; the chapters show
# the code and compare it conceptually with the educational samplers.

if (!requireNamespace("cmdstanr", quietly = TRUE)) stop("Install cmdstanr first.")
library(cmdstanr)
source("_common.R")

fit_model <- function(file, data, ...) {
  mod <- cmdstan_model(file)
  mod$sample(data = data, chains = 4, parallel_chains = 4, seed = 2026,
             iter_warmup = 1000, iter_sampling = 1000, refresh = 0, ...)
}

bioassay_data <- list(K = 4, dose = c(-0.86, -0.30, -0.05, 0.73), n = rep(5L, 4),
                      deaths = c(0L, 1L, 3L, 5L), prior_sd = 10)
fit_b <- fit_model("stan/bioassay.stan", bioassay_data)
print(fit_b$summary(c("alpha", "beta", "ld50")))

es_data <- list(J = 8L, y = eight_schools$y, sigma = eight_schools$sigma)
fit_c <- fit_model("stan/eight_schools_centered.stan", es_data)
fit_n <- fit_model("stan/eight_schools_noncentered.stan", es_data)
cat("\nCentred: divergent transitions =", sum(fit_c$diagnostic_summary()$num_divergent), "\n")
print(fit_c$summary(c("mu", "tau")))
cat("\nNon-centred: divergent transitions =", sum(fit_n$diagnostic_summary()$num_divergent), "\n")
print(fit_n$summary(c("mu", "tau")))
