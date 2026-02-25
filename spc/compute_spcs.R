# ------------------------------------------------------------
# Compute PPCs, Pop-PCs, single and divided SPC p-values
# for simple exponential family models
#
# Input:
# - X_obs: observed dataset (vector)
# - para: parameter information for generative distribution
# - dist_fun: function defining p-value type (left/right/two-sided)
# - discr_fun: user-specified discrepancy/statistic function
# - R: Monte Carlo iterations for estimating p-value
# - get_sample_from_posterior, generate_yrep_with_par:
#     model-specific functions supplied externally
#
# Output:
# - A function that returns a single p-value
# ------------------------------------------------------------


## ------------------------------------------------------------
## Posterior Predictive Check (PPC)
## Returns a function(...) to allow additional model arguments
## ------------------------------------------------------------
compute_ppc <- function(X_obs, para, dist_fun, discr_fun, R){
  
  function(...){
    
    N <- length(X_obs)
    PPC <- 0
    
    for(r in 1:R){
      
      # Draw one posterior sample of parameter
      theta <- get_sample_from_posterior(X_obs, 1, para, ...)
      
      # Generate replicated dataset from sampled parameter
      X_rep <- generate_yrep_with_par(theta, N, para)
      
      # Compute discrepancy for observed and replicated data
      # (may depend on theta in exponential family setting)
      d_obs <- discr_fun(X_obs, theta)
      d_rep <- discr_fun(X_rep, theta)
      
      # Accumulate tail comparison
      PPC <- PPC + dist_fun(d_obs, d_rep)
    }
    
    # Monte Carlo estimate of posterior predictive p-value
    return(PPC / R)
  }
}


## ------------------------------------------------------------
## Split data for single SPC
## First portion used as training set,
## remaining portion used as held-out data
## ------------------------------------------------------------
split_data <- function(X, spc_prop){
  
  N <- length(X)
  train_size <- floor(N * spc_prop)
  
  X_obs <- X[1:train_size]
  X_new <- X[-(1:train_size)]
  
  return(list("obs" = X_obs, "new" = X_new))
}


## ------------------------------------------------------------
## Single Split Predictive Check (SPC)
## 1. Split data into training and test
## 2. Fit posterior on training data
## 3. Compare replicated data to held-out data
##
## Returns a function(...) for flexible model arguments
## ------------------------------------------------------------
compute_singleSPC <- function(X_obs, para, dist_fun, discr_fun,
                              R, spc_prop = 0.5){
  
  function(...){
    
    spc <- 0
    
    # Split dataset
    splitted_data <- split_data(X_obs, spc_prop)
    X_split_obs <- splitted_data$obs
    X_split_new <- splitted_data$new
    N_split <- length(X_split_new)
    
    for(r in 1:R){
      
      # Fit posterior using training subset
      theta <- get_sample_from_posterior(X_split_obs, 1, para, ...)
      
      # Generate replicated data with same size as held-out set
      X_rep <- generate_yrep_with_par(theta, N_split, para)
      
      # Compute discrepancy on held-out and replicated data
      d_new <- discr_fun(X_split_new, theta)
      d_rep <- discr_fun(X_rep, theta)
      
      spc <- spc + dist_fun(d_new, d_rep)
    }
    
    # Monte Carlo estimate of SPC p-value
    return(spc / R)
  }
}


## ------------------------------------------------------------
## Divide dataset into K blocks
## K = floor(N^rate)
## Each block has approximately equal size
## ------------------------------------------------------------
divide_data <- function(X, rate){
  
  N <- length(X)
  K <- floor(N^(rate))
  N_split <- floor(N / K)
  
  data <- list()
  
  for(k in 1:K){
    data[[k]] <- X[((k - 1) * N_split + 1):(k * N_split)]
  }
  
  return(data)
}


## ------------------------------------------------------------
## Divided SPC (dSPC)
## 1. Divide data into K subsets
## 2. Compute single-SPC p-value on each subset
## 3. Test whether resulting p-values follow Uniform(0,1)
##
## metric:
##   - "Kolmogorov": KS test against Uniform(0,1)
##   - "2-Wass":     2-Wasserstein distance test
##
## Returns a function(...) for flexible model arguments
## ------------------------------------------------------------
compute_divided_SPC <- function(X_obs, para, rate,
                                dist_fun = ind_dis,
                                discr_fun,
                                R = 500,
                                spc_prop = 0.5,
                                metric = "Kolmogorov"){
  
  function(...){
    
    N <- length(X_obs)
    K <- floor(N^(rate))
    
    # Divide dataset
    divided_data <- divide_data(X_obs, rate)
    
    # Compute SPC p-values for each subset
    pvals <- rep(0, K)
    for(k in 1:K){
      pvals[k] <- compute_singleSPC(
        divided_data[[k]],
        para,
        dist_fun,
        discr_fun,
        R,
        spc_prop
      )(...)
    }
    
    # Uniformity test of p-values
    if(metric == "2-Wass"){
      
      # Compare empirical distribution to Uniform(0,1)
      dspc_pval <- as.numeric(
        wass_test(pvals, runif(100), p = 2)["P-Value"]
      )
      
    }else if(metric == "Kolmogorov"){
      
      # Add small jitter to avoid ties
      dspc_pval <- ks.test(
        pvals + rnorm(K, 0, 0.0001),
        "punif", 0, 1
      )$p.value
      
    }else{
      stop("invalid uniformity test")
    }
    
    return(dspc_pval)
  }
}
