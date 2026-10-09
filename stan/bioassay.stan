// Dose-response (bioassay) logistic regression, Chapters 2, 12 and 15.
// Data: Racine et al. (1986), as analysed in Gelman et al. (2013, Section 3.7).
data {
  int<lower=1> K;                    // number of dose groups
  vector[K] dose;                    // log dose (log g/ml)
  array[K] int<lower=0> n;           // animals per group
  array[K] int<lower=0> deaths;      // deaths per group
  real<lower=0> prior_sd;            // prior sd of alpha and beta
}
parameters {
  real alpha;
  real beta;
}
model {
  alpha ~ normal(0, prior_sd);
  beta ~ normal(0, prior_sd);
  deaths ~ binomial_logit(n, alpha + beta * dose);
}
generated quantities {
  real ld50 = -alpha / beta;         // dose at which P(death) = 1/2
}
