functions {
#include gpbasisfun_functions1.stan
}
data {
  int<lower=1> N_obs;      // number of observations
  int<lower=1> N_new;     
  vector[N_obs] x_obs;         // univariate covariate
  vector[N_obs] y_obs;         // target variable
  vector[N_new] x_new;         
  real<lower=0> c_f1;  // factor c to determine the boundary value L
  int<lower=1> M_f1;   // number of basis functions for smooth function
}
transformed data {
  // Normalize data
  real xmean_obs = mean(x_obs); 
 // real xmean_new = mean(x_new); 
  real ymean_obs = mean(y_obs);
  real xsd_obs = sd(x_obs);
 // real xsd_new = sd(x_new);
  real ysd_obs = sd(y_obs);
  vector[N_obs] xn_obs = (x_obs - xmean_obs)/xsd_obs; 
  vector[N_new] xn_new = (x_new - xmean_obs)/xsd_obs; 
  vector[N_obs] yn_obs = (y_obs - ymean_obs)/ysd_obs;
  // Basis functions for f1
  real L_f1 = c_f1*fmax(max(xn_obs), max(xn_new));
  matrix[N_obs,M_f1] PHI_f1 = PHI_EQ(N_obs, M_f1, L_f1, xn_obs);
  matrix[N_new,M_f1] PHI_f1_pred = PHI_EQ(N_new, M_f1, L_f1, xn_new);
}
parameters {
  real intercept;               // 
  vector[M_f1] beta_f1;         // the basis functions coefficients
  real<lower=0> lengthscale_f1; // lengthscale of f1
  real<lower=0> sigma_f1;       // scale of f1
  real<lower=0> sigma;          // residual scale
}
model {
  // spectral densities for f1
  vector[M_f1] diagSPD_f1 = diagSPD_EQ(sigma_f1, lengthscale_f1, L_f1, M_f1);
  // priors
  intercept ~ normal(0, 1);
  beta_f1 ~ normal(0, 1);
  lengthscale_f1 ~ lognormal(log(700/xsd_obs), 1);
  sigma_f1 ~ normal(0, .5);
  sigma ~ normal(0, .5);
  // model
  yn_obs ~ normal_id_glm(PHI_f1, intercept, diagSPD_f1 .* beta_f1, sigma); 
}
generated quantities {
  vector[N_obs] f;
  vector[N_new] fpred;
  vector[N_new] yrep;
  vector[N_obs] log_lik;
  {
    // spectral densities for f1
    vector[M_f1] diagSPD_f1 = diagSPD_EQ(sigma_f1, lengthscale_f1, L_f1, M_f1);
    // function scaled back to original scale
    f = (intercept + PHI_f1 * (diagSPD_f1 .* beta_f1))*ysd_obs + ymean_obs;
    fpred = (intercept + PHI_f1_pred * (diagSPD_f1 .* beta_f1))*ysd_obs + ymean_obs ;
    for (n in 1:N_new) yrep[n] = exp(normal_rng(fpred[n], sigma*ysd_obs));
    // log_liks for loo
    for (n in 1:N_obs) log_lik[n] = normal_lpdf(y_obs[n] | f[n], sigma*ysd_obs);
  }
}
