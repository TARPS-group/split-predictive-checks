functions {
#include gpbasisfun_functions.stan
}
data {
  int<lower=1> N_obs;      // number of observations
  vector[N_obs] x_obs;         // univariate covariate
  vector[N_obs] y_obs;         // target variable
  int day_of_week[N_obs];  // 
  int day_of_year[N_obs];  // 
  int<lower=1> N_new;      // number of observations
  vector[N_new] x_new;         // univariate covariate   
  real<lower=0> c_f1;  // factor c to determine the boundary value L
  int<lower=1> M_f1;   // number of basis functions for smooth function
  int<lower=1> J_f2;   // number of cos and sin functions for periodic
}
transformed data {
  // Normalize data
  real xmean_obs = mean(x_obs);
  real ymean_obs = mean(y_obs);
  real xsd_obs = sd(x_obs);
  real ysd_obs = sd(y_obs);
  vector[N_new] xn_new = (x_new - xmean_obs)/xsd_obs;
  vector[N_obs] xn_obs = (x_obs - xmean_obs)/xsd_obs;
  vector[N_obs] yn_obs = (y_obs - ymean_obs)/ysd_obs;
  // Basis functions for f1
  real L_f1 = c_f1*fmax(max(xn_obs), max(xn_new));
  matrix[N_obs,M_f1] PHI_f1 = PHI_EQ(N_obs, M_f1, L_f1, xn_obs);
  // Basis functions for f2
  real period_year = 365.25/xsd_obs;
  matrix[N_obs,2*J_f2] PHI_f2 = PHI_periodic(N_obs, J_f2, 2*pi()/period_year, xn_obs);
  // Concatenated basis functions for f1 and f2
  matrix[N_obs,M_f1+2*J_f2] PHI_f = append_col(PHI_f1, PHI_f2);
  matrix[N_new,M_f1] PHI_f1_pred = PHI_EQ(N_new, M_f1, L_f1, xn_new);
  // Basis functions for f2
  matrix[N_new,2*J_f2] PHI_f2_pred = PHI_periodic(N_new, J_f2, 2*pi()/period_year, xn_new);
  // Concatenated basis functions for f1 and f2
}
parameters {
  real intercept0;
  vector[M_f1] beta_f1;         // the basis functions coefficients for f1
  vector[2*J_f2] beta_f2;       // the basis functions coefficients for f2
  vector[6] beta_f3;            // day of week effect
  vector[366] beta_f4;          // day of year effect
  real<lower=0> lengthscale_f1; //
  real<lower=0> lengthscale_f2; //
  real<lower=0> sigma_f1;       // scale of f1
  real<lower=0> sigma_f2;       // scale of f2
  real<lower=0> sigma_f4;       // scale of day of year effect
  real<lower=0> sigma;          // residual scale
}
model {
  // spectral densities for f1 and f2
  vector[M_f1] diagSPD_f1 = diagSPD_EQ(sigma_f1, lengthscale_f1, L_f1, M_f1);
  vector[2*J_f2] diagSPD_f2 = diagSPD_periodic(sigma_f2, lengthscale_f2, J_f2);
  // day of week and day of year effects
  vector[7] f_day_of_week = append_row(0, beta_f3);
  vector[N_obs] intercept = intercept0 + f_day_of_week[day_of_week] + beta_f4[day_of_year];
  // priors
  intercept0 ~ normal(0, 1);
  beta_f1 ~ normal(0, 1);
  beta_f2 ~ normal(0, 1);
  beta_f3 ~ normal(0, 1);
  beta_f4 ~ normal(0, sigma_f4);
  lengthscale_f1 ~ lognormal(log(700/xsd_obs), 1);
  lengthscale_f2 ~ normal(0, .1);
  sigma_f1 ~ normal(0, 1);
  sigma_f2 ~ normal(0, 1);
  sigma_f4 ~ normal(0, 0.1);
  sigma ~ normal(0, 0.5);
  // model
  yn_obs ~ normal_id_glm(PHI_f,
		     intercept,
		     append_row(diagSPD_f1 .* beta_f1, diagSPD_f2 .* beta_f2),
		     sigma);
}
generated quantities {
  vector[N_obs] f1;
  vector[N_obs] f2;
  vector[N_obs] f;
  vector[N_new] f1_pred;
  vector[N_new] f2_pred;
  vector[N_new] fpred;
  vector[7] f_day_of_week;
  vector[N_obs] log_lik;
  {
    // spectral densities for f1 and f2
    vector[M_f1] diagSPD_f1 = diagSPD_EQ(sigma_f1, lengthscale_f1, L_f1, M_f1);
    vector[2*J_f2] diagSPD_f2 = diagSPD_periodic(sigma_f2, lengthscale_f2, J_f2);
  // day of week and day of year effects
    vector[7] f_day_of_week_n = append_row(0, beta_f3);
    vector[N_obs] intercept = intercept0 + f_day_of_week_n[day_of_week] + beta_f4[day_of_year];
    f_day_of_week = f_day_of_week_n*ysd_obs;
    // functions scaled back to original scale
    f1 = (intercept0 + PHI_f1 * (diagSPD_f1 .* beta_f1))*ysd_obs;
    f2 = (PHI_f2 * (diagSPD_f2 .* beta_f2))*ysd_obs;
    f = f1 + f2 + (intercept-intercept0)*ysd_obs + ymean_obs;
    // predictions
    f1_pred = (intercept0 + PHI_f1_pred * (diagSPD_f1 .* beta_f1))*ysd_obs;
    f2_pred = (PHI_f2_pred * (diagSPD_f2 .* beta_f2))*ysd_obs;
    fpred = f1_pred + f2_pred + (intercept-intercept0)*ysd_obs + ymean_obs;
    // log_liks for loo
    for (n in 1:N_obs) log_lik[n] = normal_lpdf(y_obs[n] | f[n], sigma*ysd_obs);
  }
}
