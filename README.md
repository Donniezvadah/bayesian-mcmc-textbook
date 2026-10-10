# Bayesian Statistics and Markov Chain Monte Carlo

**From Foundations to Advanced Computational Inference**: an open textbook and computational course by Donnie.

The book develops Bayesian inference and Markov chain Monte Carlo from probability foundations to research-level topics, in 27 chapters, 7 cumulative checkpoints, 7 projects and a capstone assessment. Every algorithm is derived, implemented from scratch in R, tested against known answers and diagnosed. The book is built with [Quarto](https://quarto.org) and published as a static website.

> Status: version 0.1.0 (first complete draft). The content has been checked by its author and by automated tests, but has not yet had external peer review. Corrections are welcome (see [CONTRIBUTING.md](CONTRIBUTING.md)).

## Contents

| Part | Chapters |
|:--|:--|
| I. Mathematical and computational foundations | 1–3 |
| II. Bayesian statistics from first principles | 4–8 |
| III. Monte Carlo foundations | 9–11 |
| IV. Core MCMC algorithms | 12–16 |
| V. Convergence, diagnostics and reliability | 17–19 |
| VI. Bayesian modelling in practice | 20–23 |
| VII. Advanced and research-level topics | 24–27 |
| Projects and capstone | P1–P7, capstone with rubric |
| Appendices | distributions, R code index, datasets, environment |

## Repository layout

```
_quarto.yml            book configuration (chapters, theme, bibliography)
_common.R              shared set-up sourced by every chapter (loads R/, knitr options)
index.qmd, front/      home page, how to use the book, notation
chapters/              27 chapters
checkpoints/           7 cumulative checkpoints
projects/              projects P1-P7 and the capstone
appendices/            reference appendices
R/                     course code: utilities/, models/, samplers/, diagnostics/
tests/                 testthat suite (tests/run-tests.R runs it)
stan/                  Stan versions of three models (run locally, not in CI)
scripts/               reference verification, Stan runs, link checking
references.bib         verified bibliography (BibTeX)
assets/                SCSS themes (light and dark) and CSS
_freeze/               cached computational results (committed; see below)
.github/workflows/     GitHub Actions: test, build and deploy to GitHub Pages
```

## Building locally

Requirements: R (4.3 or later), Quarto (1.4 or later), and the R packages listed in `DESCRIPTION`. The exact versions used to build the published site are pinned: R 4.6.1 and Quarto 1.10.18 in `.github/workflows/publish.yml`, and every R package in `renv.lock`.

```bash
# 1. install R package dependencies (from the repository root)
Rscript -e 'install.packages("remotes"); remotes::install_deps(dependencies = TRUE)'

# 2. run the test suite
Rscript tests/run-tests.R

# 3. render the book to _book/
quarto render

# 4. preview with live reload while editing
quarto preview
```

The book uses `freeze: auto`. Computational results are stored in `_freeze/`, and a chapter is re-executed only when its source changes. Commit `_freeze/` so that the published site and continuous integration use the same results. To re-execute everything, delete `_freeze/` before rendering.

**Exact package versions (optional).** To install the versions recorded in `renv.lock` instead of the latest ones, without switching the project to renv, run
`Rscript -e 'install.packages("renv"); renv::restore(lockfile = "renv.lock", library = .libPaths()[1], prompt = FALSE)'`.
After upgrading packages deliberately, regenerate the lockfile with
`Rscript -e 'options(repos = c(CRAN = "https://packagemanager.posit.co/cran/latest")); renv::snapshot(lockfile = "renv.lock", packages = c("knitr", "rmarkdown", "posterior", "coda", "loo", "ggplot2", "MASS", "boot", "testthat", "jsonlite"), prompt = FALSE, force = TRUE)'`,
re-run the tests, re-render, and commit `renv.lock` together with `_freeze/`.

The Stan comparisons in `stan/` are optional and are not run during the build. To run them locally, install [CmdStan](https://mc-stan.org/users/interfaces/cmdstan) and `cmdstanr`, then run `Rscript scripts/run-stan-examples.R`.

To re-check the bibliography against Crossref and arXiv (requires internet access and the `jsonlite` package): `Rscript scripts/verify-references.R`, which writes `scripts/reference-verification-report.md`.

## Publishing on GitHub Pages

The workflow `.github/workflows/publish.yml`:

1. installs pinned versions of R (4.6.1), Quarto (1.10.18) and the R packages in `renv.lock`, then runs the test suite;
2. renders the book with the `ci` profile, which adds the repository and site URLs (derived from the GitHub context, so nothing is hard-coded);
3. checks internal links;
4. deploys `_book/` to GitHub Pages.

The book is published at <https://donniezvadah.github.io/bayesian-mcmc-textbook/> from the repository <https://github.com/Donniezvadah/bayesian-mcmc-textbook>. Every push to `main` re-runs the workflow; pull requests are tested and built but not deployed.

To publish a fork or copy elsewhere: push it to the `main` branch of a GitHub repository, set **Settings → Pages → Source** to **GitHub Actions**, and push again (or run the workflow from the **Actions** tab). The site appears at `https://<user>.github.io/<repository>/`.

All links inside the book are relative, so the site works under a project sub-path.

## Licences

- Text, figures and other content: [Creative Commons Attribution 4.0 International](LICENSE-CONTENT.md) (CC BY 4.0).
- Code (`R/`, `tests/`, `scripts/`, `stan/`, `_common.R`): [MIT](LICENSE).
- Data sets are not redistributed. They come from R packages under their own licences, or are small published tables cited in the text. See Appendix C.

## Citation

See [CITATION.cff](CITATION.cff). Suggested form: Donnie (2026). *Bayesian Statistics and Markov Chain Monte Carlo: From Foundations to Advanced Computational Inference*, version 0.1.0.
