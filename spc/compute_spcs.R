
# Compute PPCs, Pop-PCs, single and divided SPCs p-values  

## Simple exponential family models

## Input: 
## - X_obs: original dataset saved as a vector
## - para: parameter information in generative distribution, mainly used in Pop-PC functions
## - dist_fun: type of p-values in use, left-tailed, right-tailed or two-sided 
## - discr_fun: statistics or discrepancy functions specified by user
## - R: number of iterations to estimate the p-value
## - get_sample_from_posterior, generate_yrep_with_par are passed directly into the function 
##   given specific models
## Output: one p-value


compute_ppc <- function(X_obs, para, dist_fun, discr_fun, R){
  function(...){
    N <- length(X_obs)
    PPC <- 0
    for(r in 1: R){
      theta <- get_sample_from_posterior(X_obs, 1, para, ...)
      X_rep <- generate_yrep_with_par(theta, N, para)
      
      d_obs <- discr_fun(X_obs, theta)
      d_rep <- discr_fun(X_rep, theta)
      PPC <- PPC + dist_fun(d_obs, d_rep)
    }
    return(PPC / R)
  }
}




compute_pop_pc_ideal <- function(X_obs, para, dist_fun, discr_fun, R){
  function(...){
    N <- length(X_obs)
    
    POP_PC_ideal <- 0
    new_data <- simulate_data_well(N)
    X_new <- new_data$X_obs
    
    for(r in 1:R){
      theta <- get_sample_from_posterior(X_obs, 1, para, ...)
      X_rep <- generate_yrep_with_par(theta, N, para)
      
      d_new <- discr_fun(X_new, theta)
      d_rep <- discr_fun(X_rep, theta)
      POP_PC_ideal <- POP_PC_ideal + dist_fun(d_new, d_rep)
    }
    return(POP_PC_ideal/ R)
  }
  
}



split_data <- function(X, spc_prop){
  N <- length(X)
  train_size <- floor(N * spc_prop)
  X_obs <- X[1: train_size]
  X_new <- X[-(1: train_size)]
  return(list("obs" = X_obs, "new" = X_new))
}



compute_singleSPC <- function(X_obs, para, dist_fun, discr_fun, 
                              R, spc_prop = 0.5){ 
  function(...){
    spc <- 0
    splitted_data <- split_data(X_obs, spc_prop)
    
    X_split_new <- splitted_data$new
    X_split_obs <- splitted_data$obs
    N_split <- length(X_split_new)
    
    
    for(r in 1:R){
      
      theta <- get_sample_from_posterior(X_split_obs, 1, para,...)
      X_rep <- generate_yrep_with_par(theta, N_split, para)
      
      d_new <- discr_fun(X_split_new, theta)
      d_rep <- discr_fun(X_rep, theta)
      # print(var(X_rep))
      # print(d_rep)
      spc <- spc + dist_fun(d_new, d_rep)
    }
    return(spc/ R)
  }
  
}



divide_data <- function(X, rate){
  N <- length(X)
  K <- floor(N^(rate))
  N_split <- floor(N/K)
  data <- list()
  for(k in 1:K){
    data[[k]] <- X[((k-1) * N_split + 1) : (k * N_split)]
  }
  return(data)
}





compute_divided_SPC <- function(X_obs, para, rate, dist_fun = ind_dis, discr_fun, 
                                R = 500, spc_prop = 0.5, metric = "Kolmogorov"){  
  
  function(...){
    N <- length(X_obs)
    K <- floor(N^(rate))
    divided_data <- divide_data(X_obs, rate)
    pvals <- rep(0, K)
    for(k in 1:K){
      pvals[k] <- compute_singleSPC(divided_data[[k]], para, dist_fun, discr_fun, R, spc_prop)(...)
    }
    if(metric == "2-Wass"){
      dspc_pval <- as.numeric(wass_test(pvals, runif(100), p = 2)["P-Value"])
    }else if(metric == "Kolmogorov"){
      dspc_pval <- ks.test(pvals + rnorm(K, 0, 0.0001), "punif", 0, 1)$p.value
    }else{
      stop('invalid uniformity test')
    }
    return(dspc_pval)
  }
  
}