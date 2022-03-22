


single_SPC_births <- function(data, model, stan_data_fun, train_prop = 0.5, lag = 30, R = 200, method, path_to_results, path_to_rep_data ){
  N <- as.numeric(count(data))
 
  if(method == "Interpolated single 0.5-SPC"){
    # extreme method
    held_out <- seq(from = 1, to = N-1, by = ceiling(1/(1 - train_prop) ))
    obs_id <- seq(from = 2, to = N, by = ceiling(1/(1 - train_prop) ))
    data_obs <- data %>% filter(id %in% obs_id )
    data_new <- data %>% filter(id %in% held_out)
  }else if (method == "Extrapolated single 0.5-SPC"){
    # consecutively split the data
    held_out <- seq(from = floor(N * train_prop)+1, to = N-1, by = 1 ) 
    obs_id <- seq(from = 1, to = floor(N * train_prop) , by = 1 ) 
    
    data_obs <- data %>% filter(id %in% obs_id )
    data_new <- data %>% filter(id %in% held_out)
  }else if(method == "PPC"){
    data_obs <- data_new <- data
 
  }


  y_new <- data_new$births_relative100
  N_new = length(data_new$id)
  
  stan_data <- stan_data_fun(data_obs, data_new)
  stan_data0 <- stan_data_fun(data, data)
  opt <- model$optimize(data = stan_data0, init = 0, algorithm='bfgs')
  odraws <- opt$draws()
  init <- sapply(c('intercept0','sigma_f1','lengthscale_f1','beta_f1','sigma'),
                 function(variable) {as.numeric(subset(odraws, variable=variable))})
  fit <- model$sample(data = stan_data, iter_warmup = 200, iter_sampling = R,
                      chains = 1, parallel_chains = 1,
                      init=function() { init })
  draws <- fit$draws()
  draws <- as_draws_matrix(draws)
  lag_num <- lag + 1
  
  pvals <- rep(0, lag_num)
  lag_new <- rep(0, lag_num)
  lag_new <- acf(y_new, lag.max = lag, type = "correlation", plot = FALSE)$acf
  lag_rep <- rep(0, lag_num)
  for(r in 1:R){
    (y_rep <- as.vector(subset(draws[r, ], variable='yrep') ))
    lag_rep <- acf(y_rep, lag.max = lag, type = "correlation", plot = FALSE)$acf
    pvals <- pvals + ind_dis(lag_new, lag_rep)
    
  }
  
  pvals <- pvals/R
  pvals <- sapply(pvals, function(x) 2 * min(x, 1-x))
  lag_labels <- sapply(seq(0, lag, by = 1), toString)
  
  # also store the last copy of fitted values
  new_results <- tibble(method = rep(method, lag_num) ,
                        lag_k = lag_labels,
                        pvals = pvals)
  
  new_results$lag_k <- factor(new_results$lag_k, levels = lag_labels)
  
  
  load(path_to_results)
  results <- rbind(new_results, results)
  save(results,  file = path_to_results)
  
  load(path_to_rep_data)
  new_df <- tibble(method = rep(method, N_new),
                   Efpred_draws = y_rep
  )
  df <- rbind(new_df, df)
  save(df,  file = path_to_rep_data)
  return(pvals)
}


plot_fit <- function(data, results, train_prop = 0.5, method_name){
  set1 <- RColorBrewer::brewer.pal(7, "Set1")
  N <- as.numeric(count(data))
  if(method_name == "Interpolated single 0.5-SPC"){
    # extreme method
    held_out <- seq(from = 1, to = N-1, by = ceiling(1/(1 - train_prop) ))
    obs_id <- seq(from = 2, to = N, by = ceiling(1/(1 - train_prop) ))
  }else if(method_name == "Extrapolated single 0.5-SPC"){
    # consecutively split the data
    held_out <- seq(from = floor(N * train_prop)+1, to = N-1, by = 1 ) 
    obs_id <- seq(from = 1, to = floor(N * train_prop) , by = 1 ) 
  }
  
  highlight_obs <- data %>% filter(id %in% obs_id)
  highlight_new <- data %>% filter(id %in% held_out) 
  N_new <- as.numeric(count(highlight_new))
  
  spc_results <- results %>% filter(method == method_name)
  y_rep <- spc_results$Efpred_draws
  highlight_new <- highlight_new %>% mutate(yrep = y_rep, held_out_id = held_out)
  fig <- data %>% 
    ggplot(aes(x = date)) +
    geom_point(data = highlight_obs, aes(y=births_relative100, color="observed data") , alpha = 0.2) +
    geom_point(data = highlight_new, aes(y=births_relative100,color="held-out data"), alpha = 0.2) +
    geom_point(data = highlight_new, aes(y = yrep, color = "predicted data") , alpha = 0.2) +
    geom_hline(yintercept=100, color='gray') +
    scale_color_manual(values=c("observed data"= set1[2], "held-out data"= set1[4], "predicted data"= set1[1])) +
    #  scale_alpha_manual(values=c(0.2, 0.2, 0.2)) +
    theme(legend.position = "bottom") +
    # geom_vline(xintercept = mid_date, color='gray') +
    labs(x="Date", y="Relative number of births", color = "")
  
  
  
  ppc_results <- results %>% filter(method == "PPC")
  data <- data %>% mutate(yrep = ppc_results$Efpred_draws)
  PPC_fig <- data %>% 
    ggplot(aes(x = date)) +
    geom_point(aes(y=births_relative100,color="observed/held-out data") , alpha = 0.2) +
    geom_point(aes(y= yrep, color= "predicted data"), alpha = 0.2) +
    geom_hline(yintercept=100, color='gray') +
    scale_color_manual(values=c("observed/held-out data"= set1[2], "predicted data"= set1[1])) +
    #  scale_alpha_manual(values=c(0.2, 0.2, 0.2)) +
    theme(legend.position = "bottom") +
    # geom_vline(xintercept = mid_date, color='gray') +
    labs(x="Date", y="Relative number of births", color = "")
  
  
  return(list("fig" = fig, "PPC_fig" = PPC_fig))
}


plot_pvals_lags_combined <- function(results){
  plt_legend <- ggplot() + 
    geom_point(data = results, aes(x = as.factor(lag_k), y = pvals, color = method), position = position_dodge(width = 1)) +
    geom_hline(yintercept = 0.05, color = "black", linetype = "dotted")  + 
    scale_color_manual(breaks = c("PPC",  "Interpolated single 0.5-SPC", "Extrapolated single 0.5-SPC","Interpolated single 0.9-SPC", "Extrapolated single 0.9-SPC"),values=c("#999999", "#0072B2", "#D55E00","#0072B2", "#D55E00")) +
    theme( axis.text.x = element_text(angle = 0),
           legend.position = c(0.75,0.8), 
           legend.title = element_blank(),
           legend.key.size = unit(0.32, "cm"),
           axis.line = element_line(size = 0.5),
           legend.background = element_rect(color = "grey", fill = "white"),
           axis.title = element_text(size = 14, color = "black"),
           axis.text = element_text(size = 11, color = "black"),
           axis.ticks = element_blank(),
           strip.text = element_text(size = 13),
           legend.text = element_text(size = 13, color = "black")) +  labs(x="lag", y="p-values")
  plt_nonlegend <- ggplot() + 
    geom_point(data = results, aes(x = as.factor(lag_k), y = pvals, color = method), position = position_dodge(width = 1)) +
    geom_hline(yintercept = 0.05, color = "black", linetype = "dotted")  + 
    scale_color_manual(breaks = c("PPC",  "Interpolated single 0.5-SPC", "Extrapolated single 0.5-SPC","Interpolated single 0.9-SPC", "Extrapolated single 0.9-SPC"),values=c("#999999", "#0072B2", "#D55E00","#0072B2", "#D55E00")) +
    theme( axis.text.x = element_text(angle = 0),
           legend.position = "none", 
           legend.title = element_blank(),
           legend.key.size = unit(0.32, "cm"),
           axis.line = element_line(size = 0.5),
           legend.background = element_rect(color = "grey", fill = "white"),
           axis.title = element_text(size = 14, color = "black"),
           axis.text = element_text(size = 11, color = "black"),
           axis.ticks = element_blank(),
           strip.text = element_text(size = 13),
           legend.text = element_text(size = 13, color = "black")) +  labs(x="lag", y="p-values")
  return(list(plt_nonlegend,plt_legend))
}