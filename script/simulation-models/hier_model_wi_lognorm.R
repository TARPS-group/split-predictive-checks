# Hierachical Gaussian example: misspecified cross observations by sampling from log normal
# Assume each group has the same num of individuals
# Assumed model: Xij ~ N(theta_i, 4), i = 1,...,I   j = 1, ..., J
#                theta_i ~ N(mu0, tau0), i = 1,...,I  
#                p(mu0, tau0) \propto 1/tau0
# Generative model:
#   well-specified case: Xij ~ N(theta_i, 4), i = 1,...,I   j = 1, ..., J
#                        theta_i ~ N(0,1)
#   misspecified case:   Xij ~ log-normal(theta_i, 4), i = 1,...,I   j = 1, ..., J
#                        theta_i ~ N(0,1)

source("generate_pvals.r")
library(MASS)
library(rstan)
rstan_options(auto_write=TRUE)
para_mod <- list("sigma" = 2)
# bayes_model <- stan_model("examples/hierarchical.stan")
bayes_model <- stan_model("examples/hierarchical_jeffrey.stan")

simulate_data_well <- function(num_of_groups, num_of_individuals){
  true_thetas <- rnorm(num_of_groups)
  X <- matrix(rnorm(num_of_groups * num_of_individuals, true_thetas, sd = 2),  num_of_groups, num_of_individuals)
  return(list("X_obs" = X, "para" = true_thetas))
}



simulate_data_mis <- function(num_of_groups, num_of_individuals){
  true_thetas <- rnorm(num_of_groups)
  X <- matrix(exp(rnorm(num_of_groups * num_of_individuals, true_thetas, sd = 2)),  num_of_groups, num_of_individuals)
  return(list("X_obs" = X, "para" = true_thetas))
}


get_stan_data <- function(X_obs){
  I_obs <- dim(X_obs)[1]
  J_obs <- dim(X_obs)[2]
  list("I" = I_obs, "J" = J_obs ,"X" = X_obs)
}

draw_mcmc_samples <- function(object, stan_data, R) {
  sampling_result <- sampling(object, data = stan_data, chains = 1, iter = R * 2)
  hier_params <- as.matrix(sampling_result, pars = c("mu0","tau0"))
  thetas <-  as.matrix(sampling_result, pars = c("theta"))
  samples <- list("hier_params" = hier_params, "thetas" = thetas)
  return(samples)
}


get_sample_from_posterior <- function(X_obs, R){
  stan_data <- get_stan_data(X_obs)
  samples <- draw_mcmc_samples(bayes_model, stan_data, R)
  return(samples)
}


draw_theta_from_hierarchy <- function(hier_params, I_new){
  mu0 <-  hier_params[, "mu0"]
  tau0 <- hier_params[, "tau0"]
  R <- length(mu0)
  theta <- map2(mu0, tau0, function(x, y) rnorm(I_new, x, y))
  theta <-  matrix(unlist(theta), I_new, R)
  return(theta)
}

generate_yrep_with_par <- function(theta, I_new, J_new, para_mod) {
  sigma0 <- para_mod$sigma
  X <- matrix(rnorm(I_new * J_new, theta, sd = sigma0),  I_new, J_new)
  return(X)
  
}





