


get_posterior_draws <- function(obs, model, R = 10000){
  n_iter <- 2 * R
  stan_data <- list("N_obs" = length(obs), "y_obs" = obs)
  sf <- sampling(model, iter = n_iter, chains = 1, data = stan_data)
  
  # get posterior draws
  sims <- as.array(sf)
  mus <- sims[, "chain:1", parameters = "mu"]
  sigmas <- sims[, "chain:1", parameters = "sigma"]
  
  y_pred <- rep(NA, R)
  for(r in 1:R){
    y_pred[r] <- rnorm(1, mean = mus[r], sd = sigmas[r])
  }
  return(list("mu" = mus, "sigma" = sigmas, "ypred" = y_pred))
}


compute_spc_light <- function(dataset, model, discr_fun, spc_prop = NULL,  R = 1000){
  
  if(is.null(spc_prop)){ # PPC
    y_obs <- y_new <- dataset  
  }else { # SPC: split data
    splitted_data <- split_data(dataset, spc_prop)
    y_obs <- splitted_data$obs
    y_new <- splitted_data$new
  }
  n_iter <- 2 * R
  stan_data <- list("N_obs" = length(y_obs), "y_obs" = y_obs)
  sf <- sampling(model, iter = n_iter, chains = 1, data = stan_data)
  
  # get posterior draws
  sims <- as.array(sf)
  mus <- sims[, "chain:1", parameters = "mu"]
  sigmas <- sims[, "chain:1", parameters = "sigma"]
  
  
  # get predictive values
  N_new <- length(y_new)
  y_rep <- list()
  d_rep <- rep(NA, R)
  d_new <- rep(NA, R)
  for(r in 1:R){
    y_rep[[r]] <- rnorm(N_new, mean = mus[r], sd = sigmas[r])
    d_rep[r] <- discr_fun(y_rep[[r]], mus[r])
    d_new[r] <- discr_fun(y_new, mus[r])
  }
  
  pval <- 2 * min(mean(d_new > d_rep), 1 - mean(d_new > d_rep)) # two-sided p-value
  
  return(pval)
}


compute_dspc_pval_light <- function(dataset, model, discr_fun, spc_prop, dspc_rate, R){
  N <- length(dataset)
  K <- floor(N^{dspc_rate})
  divided_data <- divide_data(dataset, rate = dspc_rate)
  pvals <- rep(NA, K)
  for(k in 1:K){
    pvals[k] <- compute_spc_light(divided_data[[k]], model, discr_fun, spc_prop,  R)
  }
  print(pvals)
  dspc_pval <- ks.test(pvals + rnorm(K, 0, 0.0001), "punif", 0, 1)$p.value
  return(dspc_pval)
}


run_experiment_light <- function(data, model, discr_name, data_path, R, iter = 1000){
  discr_fun <- get(discr_name)
  PPC <- rep(NA, iter)
  spc_0.5 <- rep(NA, iter)
  dspc_0.5 <- rep(NA, iter)
  for (it in 1:iter){
    PPC[it] <- compute_spc_light(dataset = data, model = model, discr_fun = discr_fun, spc_prop = NULL, R = R)
    spc_0.5[it] <- compute_spc_light(dataset = data, model = model, discr_fun = discr_fun, spc_prop = 0.5, R = R)
    dspc_0.5[it] <- compute_dspc_pval_light(dataset = data, model = model, 
                                            discr_fun = discr_fun,spc_prop = 0.5, dspc_rate = 0.49, R = R)
  }
  new_results <- tibble(method = c(rep("PPC", iter),
                                   rep("sinlge 0.5-SPC", iter),
                                   rep("divided 0.5-SPC, k = N^0.49", iter)),
                        pvals = c(PPC,spc_0.5, dspc_0.5),
                        discr_fun = rep(discr_name, 3 * iter))
  
  load(data_path)
  results <- rbind(new_results, results)
  save(results, file = data_path)

}

# plot data, posterior and posterior predictive distribution
plot_hist_comparison <- function(data, post_draws, post_pred_draws){
  den_pred <- density(post_pred_draws)
  den_post_mu <- density(post_draws)
  hist(data, breaks = 40, xlab = "Deviations from 24,800 nanoseconds", 
       ylim = c(0,0.35),xlim = c(min(light), 60),freq = FALSE, col = rgb(0,0.8,0.2), 
       border = "white", main = "")
  polygon(den_pred, col = rgb(0.2,0.6,0.8 , alpha =0.5))
  lines(den_pred, col = rgb(0.2,0.6,0.8), lwd = 2)
  polygon(den_post_mu, col = rgb(1, 0.6, 0.2 , alpha =0.5))
  lines(den_post_mu, col = rgb(1, 0.6, 0.2 ), lwd = 2)
  abline(v = 33.0, col="#666666", lty = 2, lwd = 2)
  legend("topleft", 
         c("True value", "Standard posterior distribution", "Posterior predictive distribution"), 
         lty=c(2, 1, 1), 
         lwd = c(2, 2, 2),
         col=c("#666666",rgb(1, 0.6, 0.2 ),rgb(0.2,0.6,0.8)), 
         bty = "n")
  par(lend = 1.2)
  legend(-48.15, 0.298, "Data distribution", col = rgb(0,0.8,0.2), lty = 1, lwd = 10, bty = "n")
  
}


