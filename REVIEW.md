# Reviewer brief: research-level chapters (24–27)

Thank you for agreeing to review part of *Bayesian Statistics and Markov Chain Monte Carlo: From Foundations to Advanced Computational Inference* (version 0.1.0). This brief explains the scope, what to check, and how to send comments. The chapters were checked by the author and by automated tests, but not yet by anyone else, so an independent expert reading is the most valuable check they can receive.

## Scope and time

| Chapter | Topic | Suggested reviewer background | Approx. reading time |
|:--|:--|:--|:--|
| 24 Advanced MCMC theory | general state spaces, minorisation and drift, ergodicity, CLTs, spectral gaps | Markov chain theory / probability | 3–4 h |
| 25 Specialised samplers | slice, elliptical slice, pCN, parallel tempering, SMC, data augmentation | computational statistics | 3 h |
| 26 Computational statistics and software engineering | interfaces, testing with controlled error rates, profiling, parallel RNG, simulation studies | research software engineering / simulation methodology | 2 h |
| 27 Research frontiers | automatic differentiation, ABC, variational inference, SG-MCMC, generalised Bayes | approximate / scalable Bayesian inference | 3 h |
| Checkpoint VII | cumulative problems and solution outlines | any of the above | 1 h |

One reviewer for 24–25 and one for 26–27 is a sensible split. A full review of all four chapters takes about 12 hours.

## How to read the chapters

- **Rendered:** on the published site, or locally with `quarto preview` (see `README.md`). Solutions and proofs are in collapsible boxes.
- **Source:** `chapters/24-advanced-mcmc-theory.qmd` … `chapters/27-research-frontiers.qmd`. The R code that the chapters display is in `R/`, and its tests are in `tests/testthat/`.

## What to check

Please prioritise **correctness** over style:

1. Every theorem, proposition and definition: are the assumptions sufficient (and stated), the statement correct, and the proof complete?
2. Results stated **with a citation instead of a proof**: is the statement faithful to the cited source, including theorem numbers?
3. Numerical claims in the text: do they follow from the code output shown directly above them? (Each was checked against the rendered output, but interpretation can still be wrong.)
4. Omissions: is anything important missing, outdated, or presented with misleading emphasis for a doctoral reader?

### Formal results

| Label | Result | Proof in text? |
|:--|:--|:--|
| `def-small-set`, `def-ergodicity-types`, `def-drift` | small sets and minorisation; rates of convergence; geometric drift | — |
| `thm-uniform-ergodicity` | Doeblin's condition implies uniform ergodicity | yes |
| `thm-drift-minorisation` | drift + minorisation ⇒ geometric ergodicity, with converse | **no**: cited to Meyn and Tweedie (2009, Thm 15.0.1) and Roberts and Rosenthal (2004, Thm 9), with a proof sketch |
| `prp-ar1-drift` | geometric ergodicity of the Gaussian AR(1) chain | yes |
| `prp-gap-variance` | the spectral gap bounds the asymptotic variance (reversible chains) | yes |
| `tbl-clt-conditions` | sufficient conditions for the Markov chain CLT | **no**: table of cited results |
| `def-auxiliary` | auxiliary-variable method | — |
| `prp-pcn` | pCN acceptance involves only the likelihood ratio | yes |
| `prp-pt-invariance` | parallel tempering leaves the product target invariant | yes |
| `def-correctness`, `prp-suite-error` | software vs statistical correctness; false failures in a test suite | yes |
| `def-front-ppl`, `def-front-abc`, `def-front-vi` | probabilistic program; ABC posterior; VI and ELBO | — |
| `prp-front-forward-ad` | forward-mode AD is exact | yes |
| `prp-front-abc` | target of rejection ABC and its $\varepsilon\to0$ limit | yes |
| `prp-front-elbo` | evidence decomposition | yes |
| `prp-front-meanfield` | mean-field Gaussian VI gives variances $1/\Lambda_{ii}\le\Sigma_{ii}$ | yes |
| `prp-sgld-variance` | exact stationary mean and variance of SGLD for a Gaussian mean | yes |
| `prp-front-power` | coverage of a power posterior under variance misspecification | yes |

### Points where the author would most value a second opinion

- **Ch 24.** The theorem numbers cited for `thm-drift-minorisation`, and the precise conditions in `tbl-clt-conditions` (in particular the rows for uniform and polynomial ergodicity and for reversible chains via Kipnis–Varadhan). Also whether the converse direction is stated with the right hypotheses.
- **Ch 24.** The claim, used in the experiment, that the chain targeting a $t_{1.2}$ distribution gives a non-Gaussian limit for the sample mean.
- **Ch 25.** The statement of pCN's dimension-robustness (absolute continuity with respect to the Gaussian prior), and the guidance on temperature ladders.
- **Ch 26.** The description of L'Ecuyer-CMRG stream separation (streams $2^{127}$ steps apart), and the independence assumption in `prp-suite-error`.
- **Ch 27.** The regularity conditions in part (ii) of `prp-front-abc` (dominated convergence); the attribution and interpretation of the Pareto-$\hat k$ thresholds for VI to Yao et al. (2018); the summary of optimal-scaling rates in §27.3.5; and whether the claims about reverse-mode cost are stated with appropriate care.
- **Checkpoint VII.** The solution outline for the independence sampler (uniform ergodicity iff $\pi/g$ is bounded, citing Mengersen and Tweedie, 1996), and the Pareto tail shape $k=3/4$ derived in Problem 4.

## How to send comments

Either:

- open a GitHub issue using the **Review comment** template (one issue per substantive point), or
- annotate the rendered pages (PDF print or comments in the source) and send them by email.

Please classify each comment:

- **Error**: mathematically or factually wrong; must be fixed.
- **Gap**: correct but incomplete (missing assumption, proof step or reference).
- **Clarity**: correct but likely to mislead or confuse.
- **Suggestion**: optional improvement.

With your permission, reviewers will be acknowledged by name in the book's front matter and in `CHANGELOG.md` when their comments are incorporated.
