


generate_yreps_given_obs <- function(data_obs, data_new, R){
  N_new <- nrow(data_new)
  fit <- stan_glm(Arr_delay ~ sin(month/12) + cos(month/12) + 
                    sin(DayOfWeek/7) + cos(DayOfWeek/7) + distance, 
                  data = data_obs,
                  family = neg_binomial_2(link = "log"), 
                  algorithm = "sampling",
                  QR = TRUE,
                  seed = 14124869, chains = 1, iter = 2 * R)
  
  y_rep <- matrix(posterior_predict(fit, newdata = data_new, draws = R), R, N_new)
  return(y_rep)
}



compute_singleSPC_glm <- function(data, spc_prop, discr_name, interpolation, R){
  discr_fun <- get(discr_name)
  N <- nrow(data)
  if(interpolation == TRUE){
    held_out <- seq(from = 1, to = N, by = ceiling(1/(1 - spc_prop) ))
    data_obs <- data %>% filter(!ID %in% held_out )
    data_new <- data %>% filter(ID %in% held_out) 
  }else if(interpolation == FALSE){
    held_out <- seq(from = floor(N * spc_prop), to = N, by = 1 ) 
    data_obs <- data %>% filter(!ID %in% held_out )
    data_new <- data %>% filter(ID %in% held_out) 
  }else if(interpolation == "N/A"){
    data_obs <- data_new <- data
  }
  
  N_new <- nrow(data_new)
  y_reps <- matrix(generate_yreps_given_obs(data_obs, data_new, R = R), R, N_new)
  y_new <- data_new$Arr_delay
  d_rep <- apply(y_reps, 1, function(y) discr_fun(y))
  d_new <- discr_fun(y_new)
  pval <- mean(d_rep < d_new)
  
  return(pval)
}


divide_data_glm <- function(data, rate, interpolation){
  N <- nrow(data)
  N_sub <- floor(N^(1-rate))
  K <- floor(N^(rate))
  Xsubs <- list()
  if(interpolation){
    data$intid <- c(rep(1:K, N_sub + 1), rep(NA,N%%K))
    data <- data[1:(K*(N_sub + 1)), ]
    for(k in 1:K){
      Xsubs[[k]] <- data[data$intid == k, ]
    }
  }else if(!interpolation){
    for(k in 1:(K-1)){
      Xsubs[[k]] <- data[((k-1) * N_sub + 1) : (k * N_sub), ]
    }
    Xsubs[[K]] <- data[-(1 : ((K-1) * N_sub)), ]
  }else{
    stop("Error: missing interpolation type")
  }
  return(Xsubs)
}


compute_divided_SPC_glm <- function(data, spc_prop, dspc_rate, discr_name, spc_interpolation, divide_interpolation, R){
  N <- nrow(data)
  Xsubs <- divide_data_glm(data, dspc_rate, divide_interpolation)
  K <- floor(N^(rate))
  pvals <- rep(NA, K)
  for(k in 1:K){
    pvals[k] <- compute_singleSPC_glm(Xsubs[[k]], spc_prop = 0.5, 
                                      discr_name = discr_name, 
                                      interpolation = spc_interpolation, R = R)
  }
  dspc_pval <- ks.test(pvals + rnorm(K, 0, 0.0001), "punif", 0, 1)$p.value
  return(dspc_pval)
}


run_experiment_negbin <- function(data, discr_name, result_path, R){
  ppc <- compute_singleSPC_glm(flights, spc_prop = NULL, discr_name, interpolation = "N/A", R = R)
  singleSPC_int <- compute_singleSPC_glm(flights, spc_prop = 0.5, discr_name , interpolation = TRUE, R = R)
  singleSPC_ext <- compute_singleSPC_glm(flights, spc_prop = 0.5, discr_name, interpolation = FALSE, R = R)
  intdivided_intSPC <- compute_divided_SPC_glm(flights, spc_prop = 0.5, dspc_rate = 0.49, discr_name, spc_interpolation = TRUE, divide_interpolation = TRUE, R = R)
  intdivided_extSPC <- compute_divided_SPC_glm(flights, spc_prop = 0.5, dspc_rate = 0.49, discr_name, spc_interpolation = FALSE, divide_interpolation = TRUE, R = R)
  load(result_path) 
  new_results <- tibble(method = c( 
    "PPC",
    "interpolated single 0.5-SPC",
    "extrapolated single 0.9-SPC",
    "double-interpolated divided 0.5-SPC, k = N^0.49",
    "interpolated divided extrapolated 0.5-SPC, k = N^0.49"),
    discr_fun = rep(discr_name, 5),
    pvals = c(ppc, singleSPC_int, singleSPC_ext,
              intdivided_intSPC,intdivided_extSPC))
  
  results <- rbind(new_results, results)
  save(results, file = result_path)
  
}

plot_hist_comparison <- function(data, y_rep_geom, y_rep_nb){
  plt_nyc <- hist(y_nyc, 
                  breaks = seq(min(y_nyc), max(y_nyc), length.out = 27), 
                  plot = FALSE)
  plt_rep_nb <- hist(y_rep_nb, 
                     breaks = seq(min(y_rep_nb), max(y_rep_nb), length.out = 65),  
                     plot = FALSE)
  plt_rep_geom <- hist(y_rep_geom, 
                       breaks = seq(min(y_rep_geom), max(y_rep_geom), length.out = 4),   
                       plot = FALSE)
  
  plot(plt_nyc,yaxt='n', ylim = c(0, 6), col = rgb(0,0.8,0.2, alpha = 0.8), 
       border = "white", ylab = "Count", xlab = "Arrival delays (in minutes)", main = "")
  plot(plt_rep_nb, col = rgb(0.2,0.6,0.8 , alpha = 0.6), border = "white", add = TRUE)
  plot(plt_rep_geom,col = rgb(1, 0.6, 0.2 , alpha =0.6), border = "white", add = TRUE)
  eaxis(side = 2, at = seq(0, 6,  length.out = 3), 
        labels = pretty10exp(10^seq(0, 6,  length.out = 3), drop.1=TRUE, sub10 = 0))
  par(lend = 1.2)
  legend(450, 6 , 
         c("Data distribution", 
           "Posterior predictive distribution - NB/GLM", 
           "Posterior predictive distribution - Geom"), 
         col = c(rgb(0,0.8,0.2, alpha = 0.8), rgb(0.2,0.6,0.8 , alpha = 0.8), 
                 rgb(1, 0.6, 0.2 , alpha =0.8)), lty = c(1,1), lwd = c(10,10), bty = "n")
  
}