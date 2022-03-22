library("rprojroot")
library(cmdstanr)
library(posterior)
options(pillar.neg = FALSE, pillar.subtle=FALSE, pillar.sigfig=2)

model6 <- cmdstan_model(stan_file = root("script/experiments/birthdays", "gpbf6.stan"), include_paths = root("script/experiments/birthdays"))

para_mod <- NULL

stan_data6 <- function(data_obs, data_new){
  stan_data <- list(x_obs=data_obs$id,
                    y_obs=log(data_obs$births_relative100),
                    N_obs=length(data_obs$id),
                    x_new=data_new$id,
                    N_new=length(data_new$id),
                    c_f1=1.5, # factor c of basis functions for GP for f1
                    M_f1=10, # number of basis functions for GP for f1
                    J_f2=20, # number of basis functions for periodic f2
                    day_of_week=data_obs$day_of_week,
                    day_of_year=data_obs$day_of_year2)
                  #  day_of_week_new=data_new$day_of_week,
                  #  day_of_year_new=data_new$day_of_year2) # 1st March = 61 every year
  return(stan_data)
}

generate_yrep_with_par <- function(stan_data){
  opt6 <- model6$optimize(data=stan_data, init=0, algorithm='lbfgs',history=100, tol_obj=10)
  
  
  
  
  opt1 <- model1$optimize(data = stan_data, init = 0, algorithm='bfgs')
  odraws1 <- opt1$draws()
  init1 <- sapply(c('intercept','sigma_f1','lengthscale_f1','beta_f1','sigma'),
                  function(variable) {as.numeric(subset(odraws1, variable=variable))})
  fit1 <- model1$sample(data = stan_data, iter_warmup=100, iter_sampling=100,
                        chains= 3, parallel_chains = 3,
                        init=function() { init1 })
  draws1 <- fit1$draws()
  draws1 <- as_draws_matrix(draws1)
  yrep <- exp(apply(subset(draws1, variable='fpred'), 2, median))
  yobs_fit <- exp(apply(subset(draws1, variable='f'), 2, median))
  return(list("Ef_new" = yrep, "Ef_obs" = yobs_fit))
}