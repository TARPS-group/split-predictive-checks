
# Model 1: Gaussian location model
# mu unknown, sigma known
# Assumued model:   X ~ N(mu, 1), 
#                   mu ~ N(mu0, sigma0)
# Generative Model: 
#   Well-specified: X ~ N(0,1), i.e., sigma = 1
#   Misspecified:  X ~ N(0,15^2), i.e., sigma = 15



# gaussian <- list("true_pars" = list("mu" = 0, "sigma" = 1), "mis_par" = 3)    
source("generate_pvals.r")

para_mod <- list("mu" = NA, "sigma" = 1)

simulate_data_well <- function(sample_size){
  mu <- 0
  sigma <- 1
  X <- rnorm(sample_size, mean = mu, sd = sigma)
  para <- list("mu" = mu, "sigma" = sigma)
  data <- list("X_obs" = X, "para" = para)
  return(data)
}


simulate_data_mis <- function(sample_size){
  mu <- 0
  sigma <- 15
  X <- rnorm(sample_size, mean = mu, sd = sigma)
  para <- list("mu" = mu, "sigma" = sigma)
  data <- list("X_obs" = X, "para" = para)
  return(data)
}


get_sample_from_posterior <- function(X, num_of_samples, para, mu0, sigma0){
  N <- length(X)
  sigma <- para$sigma
  muN <-  1/(1/sigma0^2 + N/sigma^2) * (mu0/sigma0^2 + sum(X)/sigma^2)
  sigmaN <- sqrt(1/(1/sigma0^2 + N/sigma^2))
  theta <- rnorm(num_of_samples, mean = muN, sd = sigmaN)
  return(theta)
}


generate_yrep_with_par <- function(theta, size, para){
  sigma <- para$sigma
  return(rnorm(size, mean = theta, sd = sigma))
}




run_gaussian_model <- function(sizes, data_path, discr_name, iter = 200){
  load(data_path)
  dat <- updated_data
  for(i in 1: length(sizes)){
    dat <- run_experiment(N = sizes[i], iter = iter, data = dat,
                                 discr_fun_name = discr_name, file_name = data_path,
                                 misspecified = FALSE)(mu0 = 0, sigma0 = 10000)
    dat <- run_experiment(N = sizes[i], iter = iter, data = dat,
                                 discr_fun_name = discr_name, file_name = data_path,
                                 misspecified = TRUE)(mu0 = 0, sigma0 = 10000)
  }
}



