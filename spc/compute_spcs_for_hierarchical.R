

## compute spcs for two-level hierarchical models

compute_ppc_hier <- function(X_obs, para, dist_fun, discr_fun, R){
  
  I <- dim(X_obs)[1]  
  J <- dim(X_obs)[2] 
  
  PPC <- 0
  sampling <- get_sample_from_posterior(X_obs,  R)
  thetas <- sampling$thetas
  for(r in 1: R){
    X_rep <- generate_yrep_with_par(thetas[r, ], I, J, para)
    d_obs <- discr_fun(X_obs)
    d_rep <- discr_fun(X_rep)
    PPC <- PPC + dist_fun(d_rep, d_obs)
  }
  return(PPC / R)
  
}




compute_pop_pc_ideal_hier <- function(X_obs, para, discr_fun, dist_fun, R = 500){
  
  I <- dim(X_obs)[1]
  J <- dim(X_obs)[2]
  
  new_data <- simulate_data_well(I, J)
  X_new <- new_data$X_obs
  d_new <- discr_fun(X_new)
  
  POP_PC <- 0
  samples <- get_sample_from_posterior(X_obs, R)
  thetas <- samples$thetas
  for (r in 1:R){
    
    X_rep <- generate_yrep_with_par(thetas[r, ], I, J, para)
    
    d_rep <- discr_fun(X_rep)
    
    POP_PC <- dist_fun(d_rep, d_new) + POP_PC
  }
  POP_PC <- POP_PC/ R
  return(POP_PC)
  
}





split_data_hier <- function(X, train_prop = 0.1, cross_split = TRUE){
  I <- dim(X)[1]
  J <- dim(X)[2]
  if(cross_split){
    I_new <- ceiling(I * (1-train_prop))
    index_new <- sample(I, I_new, replace = FALSE) 
    X_new <- X[index_new, ]
    X_obs <- X[-(index_new), ]
  }else{
    J_new <- ceiling(J * (1 - train_prop))
    index_new <- sample(J, J_new, replace = FALSE) 
    X_new <- X[, index_new]
    X_obs <- X[, -(index_new)]
  }
  return(list("obs" = X_obs, "new" = X_new))
  
}


compute_singleSPC_hier <- function(X_obs, para, discr_fun, dist_fun, train_prop = 0.5, R = 1000, cross_split = TRUE){
  
  splitted_data <- split_data_hier(X_obs, train_prop, cross_split = cross_split)
  X_split_obs <- splitted_data$obs
  X_split_new <- splitted_data$new
  
  I_new <- dim(X_split_new)[1]
  J_new <- dim(X_split_new)[2]
  
  d_new <- discr_fun(X_split_new)
  
  POP_PC <- 0
  pos_samples <- get_sample_from_posterior(X_split_obs, R)
  if(cross_split){
    mu0_sigma0 <- pos_samples$hier_params
    thetas <- draw_theta_from_hierarchy(mu0_sigma0, I_new)
  }else{
    thetas <- t(pos_samples$thetas)
  }
  
  for(r in 1:R){
    X_rep <- generate_yrep_with_par(thetas[, r], I_new, J_new, para)
    d_rep <- discr_fun(X_rep)
    POP_PC <- dist_fun(d_rep, d_new) + POP_PC
  }
  
  POP_PC <- POP_PC/ R
  return(POP_PC)
}


divide_data_hier <- function(X, cross_divide = FALSE){
  I <- dim(X)[1]
  J <- dim(X)[2]
  data <- list()
  if(cross_divide){
    K <- floor(sqrt(I))
    I_size <- floor(I/K)
    for(k in 1: (K - 1)){
      data[[k]] <- X[((k-1) * I_size + 1) : (k * I_size), ]
    }
    data[[K]] <- X[-(1 : ((K-1) * I_size)), ]
  }else{
    K <- floor(sqrt(J))
    J_size <- floor(J/K)
    for(k in 1: (K - 1)){
      data[[k]] <- X[, ((k-1) * J_size + 1) : (k * J_size)]
    }
    data[[K]] <- X[, -(1 : ((K-1) * J_size))]
  }
  return(list("data" = data, "K" = K))
}


compute_dividedSPC_hier <- function(X_obs, para, dist_fun = ind_dis, discr_fun, 
                                    R = 500, train_prop = 0.5,  
                                    cross_divide = FALSE, cross_split = FALSE){  # try R = 5000
  
  divided_data <- divide_data_hier(X_obs, cross_divide = cross_divide)
  X_div <- divided_data$data
  K <- divided_data$K
  
  pvals <- rep(0, K)
  for(k in 1:K){
    pvals[k] <- compute_singleSPC_hier(X_div[[k]], para, discr_fun, dist_fun, 
                                       train_prop = train_prop, R = R, cross_split = cross_split)
  }
  
  KSpval <- ks.test(pvals + rnorm(K, 0, 0.0001), "punif", 0, 1)$p.value
  return(KSpval)
}


