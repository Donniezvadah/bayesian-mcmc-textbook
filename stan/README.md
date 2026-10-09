# Stan programs

Stan versions of example models from the course, for comparison with the educational samplers in `R/`. They are run by `scripts/run-stan-examples.R` (requires `cmdstanr` and CmdStan) and are **not** executed when the book is rendered, so the website builds without a C++ toolchain.

| File | Model | Chapters |
|:--|:--|:--|
| `bioassay.stan` | logistic dose–response | 2, 12, 15 |
| `eight_schools_centered.stan` | normal hierarchical model, centred | 13, 15, 21 |
| `eight_schools_noncentered.stan` | normal hierarchical model, non-centred | 13, 15, 21 |

## Reference run

Output of `Rscript scripts/run-stan-examples.R` on 9 October 2026 (CmdStan 2.39.0, cmdstanr 0.9.0, R 4.6.1, Apple M4, macOS 26.6.2). Four chains, 1000 warm-up and 1000 sampling iterations each, `seed = 2026`. Results will differ slightly with other versions and platforms; conclusions should not.

**Bioassay** (prior sd 10):

| variable | mean | sd | q5 | q95 | rhat | ess_bulk |
|:--|--:|--:|--:|--:|--:|--:|
| alpha | 0.942 | 0.923 | −0.465 | 2.51 | 1.00 | 1184 |
| beta | 8.81 | 3.82 | 3.52 | 15.6 | 1.00 | 1263 |
| ld50 | −0.0994 | 0.109 | −0.256 | 0.0834 | 1.00 | 1641 |

**Eight schools, centred:** 113 of 4000 transitions divergent; 2 of 4 chains with E-BFMI below 0.3.

| variable | mean | sd | q5 | q95 | rhat | ess_bulk | ess_tail |
|:--|--:|--:|--:|--:|--:|--:|--:|
| mu | 4.29 | 3.30 | −1.12 | 9.46 | 1.00 | 462 | 416 |
| tau | 3.50 | 2.36 | 0.920 | 8.35 | 1.02 | 188 | 166 |

**Eight schools, non-centred:** 0 divergent transitions.

| variable | mean | sd | q5 | q95 | rhat | ess_bulk | ess_tail |
|:--|--:|--:|--:|--:|--:|--:|--:|
| mu | 4.44 | 3.27 | −1.05 | 9.75 | 1.00 | 4023 | 2525 |
| tau | 3.32 | 2.52 | 0.238 | 8.09 | 1.00 | 2331 | 1572 |

Note how the centred fit's 5% quantile of `tau` (0.92) is far above the non-centred fit's (0.24): the centred sampler fails to explore the neck of the funnel at small `tau`, and the divergences are the warning sign (Chapters 15 and 21).
