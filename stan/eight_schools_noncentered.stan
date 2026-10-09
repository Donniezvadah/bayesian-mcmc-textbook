// Eight schools, non-centred parameterisation (Chapters 13, 15, 21).
data {
  int<lower=1> J;
  vector[J] y;
  vector<lower=0>[J] sigma;
}
parameters {
  real mu;
  real<lower=0> tau;
  vector[J] eta;
}
transformed parameters {
  vector[J] theta = mu + tau * eta;
}
model {
  mu ~ normal(0, 5);
  tau ~ normal(0, 5);
  eta ~ std_normal();
  y ~ normal(theta, sigma);
}
