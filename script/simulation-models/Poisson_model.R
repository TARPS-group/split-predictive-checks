# Simple model I: Poisson Model
# Assumed model: X ~ Possion(theta)
#                theta ~ Gamma(a0, b0)
# Generative model:
#     misspecified:  X ~ negbin(mu = 2, theta = 0.01)   # mu = 2, var = mu + mu^2/ theta = 402
#     well_specified:  X ~ pois(2)  # mu = var = 2


library(MASS)

para_mod <- NULL
 
simulate_data_well <- function(sample_size){
  para <- list("theta" = 2)
  X <- rpois(sample_size, para$theta)
  data <- list("X_obs" = X, "para" = para)
  return(data)
}

# P_0 = NegBin(2, 0.01) over-dispersed
simulate_data_mis <- function(sample_size){
  para <- list("mu" = 2, "theta" = 0.01)
  X <- rnegbin(sample_size, mu = para$mu, theta = para$theta)
  data <- list("X_obs" = X, "para" = para)
  return(data)
}

# # P_0 = Binom(2, 0.1) under-dispersed
# simulate_data_mis <- function(sample_size){
#   para <- list("n" = 2, "p" = 0.1)
#   X <- rbinom(sample_size, size = para$n, prob = para$p)
#   data <- list("X_obs" = X, "para" = para)
#   return(data)
# }

get_sample_from_posterior <- function(X, num_of_samples, para, a0, b0){
  N <- length(X)
  aN <- a0 + sum(X)
  bN <- b0 + N
  return(rgamma(num_of_samples, shape = aN, rate = bN))
}


generate_yrep_with_par <- function(theta, size, para){
  return(rpois(size, theta))
}


run_poisson_model <- function(sizes, data_path, discr_name, iter = 200){
  load(data_path)
  dat <- updated_data
  for(i in 1: length(sizes)){
    dat <- run_experiment(N = sizes[i], iter = iter, data = dat,
                                   discr_fun_name = discr_name,file_path = data_path,
                                   misspecified = FALSE)(a0 = 0.1, b0 = 0.2)
    dat <- run_experiment(N = sizes[i], iter = iter, data = dat, 
                                   discr_fun_name = discr_name,file_path = data_path,
                                   misspecified = TRUE)(a0 = 0.1, b0 = 0.2)
  }
}
