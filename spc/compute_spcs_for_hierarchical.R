## ------------------------------------------------------------
## Posterior Predictive Check (PPC) for two-level hierarchical model
## ------------------------------------------------------------
compute_ppc_hier <- function(X_obs, para, dist_fun, discr_fun, R){
  
  # Dimensions: I = number of groups, J = number of observations per group
  I <- dim(X_obs)[1]  
  J <- dim(X_obs)[2] 
  
  PPC <- 0
  
  # Draw R posterior samples of parameters given observed data
  sampling <- get_sample_from_posterior(X_obs, R)
  thetas <- sampling$thetas
  
  for(r in 1:R){
    
    # Generate replicated dataset from posterior draw r
    X_rep <- generate_yrep_with_par(thetas[r, ], I, J, para)
    
    # Compute discrepancy statistic for observed and replicated data
    d_obs <- discr_fun(X_obs)
    d_rep <- discr_fun(X_rep)
    
    # Accumulate indicator / distance comparison
    PPC <- PPC + dist_fun(d_rep, d_obs)
  }
  
  # Monte Carlo estimate of posterior predictive p-value
  return(PPC / R)
}



## ------------------------------------------------------------
## Split data into training (observed) and held-out (new) portions
## cross_split = TRUE  -> split across groups (rows)
## cross_split = FALSE -> split across items   (columns)
## ------------------------------------------------------------
split_data_hier <- function(X, train_prop = 0.1, cross_split = TRUE){
  
  I <- dim(X)[1]
  J <- dim(X)[2]
  
  if(cross_split){
    # Hold out a subset of groups
    I_new <- ceiling(I * (1 - train_prop))
    index_new <- sample(I, I_new, replace = FALSE) 
    X_new <- X[index_new, ]
    X_obs <- X[-(index_new), ]
    
  }else{
    # Hold out a subset of items
    J_new <- ceiling(J * (1 - train_prop))
    index_new <- sample(J, J_new, replace = FALSE) 
    X_new <- X[, index_new]
    X_obs <- X[, -(index_new)]
  }
  
  return(list("obs" = X_obs, "new" = X_new))
}


## ------------------------------------------------------------
## Single split predictive check (SPC)
## ------------------------------------------------------------
compute_singleSPC_hier <- function(X_obs, para, discr_fun, dist_fun,
                                   train_prop = 0.5, R = 1000,
                                   cross_split = TRUE){
  
  # Split into training and held-out sets
  splitted_data <- split_data_hier(X_obs, train_prop,
                                   cross_split = cross_split)
  X_split_obs <- splitted_data$obs
  X_split_new <- splitted_data$new
  
  I_new <- dim(X_split_new)[1]
  J_new <- dim(X_split_new)[2]
  
  # Discrepancy on held-out data
  d_new <- discr_fun(X_split_new)
  
  POP_PC <- 0
  
  # Fit posterior using training data
  pos_samples <- get_sample_from_posterior(X_split_obs, R)
  
  if(cross_split){
    # When splitting across groups:
    # draw new group-level parameters from inferred hierarchy
    mu0_sigma0 <- pos_samples$hier_params
    thetas <- draw_theta_from_hierarchy(mu0_sigma0, I_new)
    
  }else{
    # When splitting across items:
    # reuse posterior draws of group-level parameters
    thetas <- t(pos_samples$thetas)
  }
  
  for(r in 1:R){
    
    # Generate replicated held-out data
    X_rep <- generate_yrep_with_par(thetas[, r], I_new, J_new, para)
    d_rep <- discr_fun(X_rep)
    
    # Compare discrepancy to held-out data
    POP_PC <- dist_fun(d_rep, d_new) + POP_PC
  }
  
  POP_PC <- POP_PC / R
  return(POP_PC)
}


## ------------------------------------------------------------
## Divide data into K roughly equal blocks
## K ≈ sqrt(I) if dividing by groups
## K ≈ sqrt(J) if dividing by items
## ------------------------------------------------------------
divide_data_hier <- function(X, cross_divide = FALSE){
  
  I <- dim(X)[1]
  J <- dim(X)[2]
  data <- list()
  
  if(cross_divide){
    
    # Divide across groups
    K <- floor(sqrt(I))
    I_size <- floor(I / K)
    
    for(k in 1:(K - 1)){
      data[[k]] <- X[((k - 1) * I_size + 1):(k * I_size), ]
    }
    
    data[[K]] <- X[-(1:((K - 1) * I_size)), ]
    
  }else{
    
    # Divide across items
    K <- floor(sqrt(J))
    J_size <- floor(J / K)
    
    for(k in 1:(K - 1)){
      data[[k]] <- X[, ((k - 1) * J_size + 1):(k * J_size)]
    }
    
    data[[K]] <- X[, -(1:((K - 1) * J_size))]
  }
  
  return(list("data" = data, "K" = K))
}


## ------------------------------------------------------------
## Divided SPC:
## 1. Divide data into K subsets
## 2. Compute SPC p-value on each subset
## 3. Test uniformity of p-values using KS test
## ------------------------------------------------------------
compute_dividedSPC_hier <- function(X_obs, para,
                                    dist_fun = ind_dis,
                                    discr_fun,
                                    R = 500,
                                    train_prop = 0.5,
                                    cross_divide = FALSE,
                                    cross_split = FALSE){
  
  divided_data <- divide_data_hier(X_obs,
                                   cross_divide = cross_divide)
  X_div <- divided_data$data
  K <- divided_data$K
  
  pvals <- rep(0, K)
  
  for(k in 1:K){
    pvals[k] <- compute_singleSPC_hier(
      X_div[[k]],
      para,
      discr_fun,
      dist_fun,
      train_prop = train_prop,
      R = R,
      cross_split = cross_split
    )
  }
  
  # Small jitter avoids ties in KS test
  KSpval <- ks.test(
    pvals + rnorm(K, 0, 0.0001),
    "punif",
    0, 1
  )$p.value
  
  return(KSpval)
}
