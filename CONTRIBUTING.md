# Contributing

Thank you for helping improve this book. Corrections of mathematical errors, broken code, inaccurate citations and unclear explanations are especially valuable.

## Reporting a problem

Open an issue describing:

- the chapter and section (or the file and line);
- what is wrong, and what you believe is correct, with a derivation, a reference or a minimal reproducible R example;
- your R and package versions (`sessionInfo()`) if the problem is computational.

If a computed number in your build differs from the book, first check whether the difference is within Monte Carlo error (the book reports MCSEs). Numerical results can differ slightly across platforms and R versions.

## Proposing a change

1. Fork the repository and create a branch.
2. Make the change. Follow the existing structure of chapters (the 16-element pattern described in "How to use this book") and the notation page.
3. If you change code in `R/`, add or update tests in `tests/testthat/` and run `Rscript tests/run-tests.R`. Statistical tests should use MCSE-based tolerances (`expect_within_mcse()`), not arbitrary ones.
4. Render the affected chapter (`quarto render chapters/<file>.qmd`) and check that every statement in the text about a computed number agrees with the rendered output.
5. If you add a citation, add a complete BibTeX entry with a DOI or stable URL, and run `Rscript scripts/verify-references.R`. Never add a reference you have not checked.
6. Open a pull request explaining the change.

## Standards

- Mathematical results state their assumptions. Proofs are complete, or the gap is stated with a reference.
- Code is plain R, validates its inputs, works on the log scale, and is tested against known answers.
- Claims about real data are supported by documented provenance. Synthetic data are labelled as such.
- Be respectful and constructive in all discussions.

By contributing you agree that your contributions are licensed under the repository's licences: MIT for code and CC BY 4.0 for content.
