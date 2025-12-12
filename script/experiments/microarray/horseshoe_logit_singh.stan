// horseshoe_logit_singh.stan
data {
  int<lower=1> N;          // number of samples
  int<lower=1> P;          // number of predictors (genes)
  matrix[N, P] X;          // design matrix
  int<lower=0, upper=1> y[N]; // binary outcomes
}

parameters {
  real alpha;                 // intercept
  vector[P] beta_raw;         // raw coefficient vector
  vector<lower=0>[P] lambda;  // local shrinkage parameters
  real<lower=0> tau_unif;     // sampled global scale if tau <= 0
}

transformed parameters {
  vector[P] beta;             // effective coefficients
  for (j in 1:P) {
    beta[j] = beta_raw[j] * lambda[j] * tau_unif;
  }
}

model {
  // Priors
  alpha ~ normal(0, 5);
  beta_raw ~ normal(0, 1);
  lambda ~ cauchy(0, 1);
  tau_unif ~ cauchy(0, 1);

  // Likelihood
  y ~ bernoulli_logit(X * beta + rep_vector(alpha, N));
}

generated quantities {
  int y_rep[N];   // posterior predictive draws
  vector[N] p;    // probabilities
  vector[N] y_prior_rep;
  for (i in 1:N) {
    p[i] = inv_logit(alpha + dot_product(row(X, i), beta));
    y_rep[i] = bernoulli_rng(p[i]);
  }
  for (n in 1:N)
    y_prior_rep[n] = bernoulli_logit_rng(dot_product(X[n], beta));
}
