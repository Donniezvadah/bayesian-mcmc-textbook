# A small, reusable container for multi-chain MCMC output (Chapter 26) ----------------

#' Run several chains of any sampler with reproducible, independent random-number
#' streams, optionally in parallel.
#'
#' @param sampler function(chain_id) returning a draws matrix (iterations x parameters)
#'   with column names; all chains must return the same dimensions.
#' @param n_chains number of chains.
#' @param seed integer seed for the L'Ecuyer-CMRG generator; chain m uses the m-th
#'   stream, so results do not depend on whether chains run sequentially or in parallel.
#' @param cores number of processes (forking via parallel::mclapply; ignored on Windows).
#' @return an object of class "mcmc_fit": a 3-d array draws[iteration, chain, parameter]
#'   plus timing and environment information.
run_chains <- function(sampler, n_chains = 4L, seed = 1L, cores = 1L) {
  check_function(sampler, "sampler"); n_chains <- check_count(n_chains, "n_chains")
  old_kind <- RNGkind()[1]
  on.exit(RNGkind(old_kind), add = TRUE)
  RNGkind("L'Ecuyer-CMRG"); set.seed(seed)
  streams <- vector("list", n_chains); s <- .Random.seed
  for (m in seq_len(n_chains)) { streams[[m]] <- s; s <- parallel::nextRNGStream(s) }
  one_chain <- function(m) { assign(".Random.seed", streams[[m]], envir = globalenv()); sampler(m) }
  t0 <- proc.time()[["elapsed"]]
  res <- if (cores > 1L && .Platform$OS.type != "windows") {
    parallel::mclapply(seq_len(n_chains), one_chain, mc.cores = cores, mc.set.seed = FALSE)
  } else lapply(seq_len(n_chains), one_chain)
  elapsed <- proc.time()[["elapsed"]] - t0
  dims <- vapply(res, function(d) dim(as.matrix(d)), integer(2))
  if (any(dims[1, ] != dims[1, 1]) || any(dims[2, ] != dims[2, 1])) course_stop("Chains returned different dimensions.")
  pars <- colnames(as.matrix(res[[1]]))
  if (is.null(pars)) pars <- paste0("par", seq_len(dims[2, 1]))
  draws <- array(NA_real_, dim = c(dims[1, 1], n_chains, dims[2, 1]), dimnames = list(NULL, NULL, pars))
  for (m in seq_len(n_chains)) draws[, m, ] <- as.matrix(res[[m]])
  structure(list(draws = draws, seconds = elapsed, seed = seed, cores = cores), class = "mcmc_fit")
}

#' Summary table: mean, sd, 5% and 95% quantiles, rank-normalised R-hat, bulk/tail ESS, MCSE.
summary.mcmc_fit <- function(object, ...) {
  pars <- dimnames(object$draws)[[3]]
  out <- do.call(rbind, lapply(pars, function(p) summarise_draws_basic(object$draws[, , p], p)))
  rownames(out) <- NULL
  out
}

print.mcmc_fit <- function(x, digits = 3, ...) {
  d <- dim(x$draws)
  cat(sprintf("MCMC fit: %d chains x %d iterations, %d parameters (%.2f s, seed %d)\n",
              d[2], d[1], d[3], x$seconds, x$seed))
  s <- summary(x); num <- vapply(s, is.numeric, logical(1))
  s[num] <- lapply(s[num], signif, digits = digits)
  print(s, row.names = FALSE)
  invisible(x)
}
