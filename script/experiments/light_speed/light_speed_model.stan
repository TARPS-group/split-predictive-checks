data {
  int<lower=1> N_obs;
  vector[N_obs] y_obs;
}
parameters {
  real mu;
  real<lower=0> sigma;
}
model {
  for(i in 1:N_obs)
    target += normal_lpdf(y_obs[i] | mu, sigma);
  target += -log(sigma);
}
