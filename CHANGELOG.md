# Changelog

## Unreleased

- Continuous integration now pins R 4.6.1 and Quarto 1.10.18, and restores R packages at the versions recorded in the new `renv.lock`.
- Reviewer brief for the research-level chapters (24-27) in `REVIEW.md`.
- First deployment to GitHub Pages; CI runner image pinned to ubuntu-24.04.

## 0.1.0 (2026-10-09)

First complete draft.

- 27 chapters in seven parts, each with objectives, derivations and proofs, worked examples, R implementations, computational experiments, diagnostics, exercises and separate solutions.
- Seven cumulative checkpoints, seven projects (P1-P7) and a capstone assessment with rubric.
- Appendices: distributions, R code index, datasets and provenance, computational environment.
- Course R code in `R/` (samplers, models, diagnostics, utilities) with a `testthat` suite.
- Stan versions of three models with a recorded reference run.
- Verified bibliography with a verification script and report.
- Light and dark themes; GitHub Actions workflow for testing, building and deploying to GitHub Pages.
