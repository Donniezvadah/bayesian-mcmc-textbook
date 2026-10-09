// Eight schools, centred parameterisation (Chapters 13, 15, 21).
data {
  int<lower=1> J;
  vector[J] y;
  vector<lower=0>[J] sigma;
}
parameters {
  real mu;
  real<lower=0> tau;
  vector[J] theta;
}
model {
  mu ~ normal(0, 5);
  tau ~ normal(0, 5);                // half-normal because of the lower bound
  theta ~ normal(mu, tau);
  y ~ normal(theta, sigma);
}
