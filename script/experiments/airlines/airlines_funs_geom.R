### Script for computing SPC p-values and analyzing results for airline delays data ###
### Analysis includes tests on a geometric model and negative binomial glm model.




run_experiment_geom <- function(N_obs, data, discr_name, results_path, R){
  discr_fun = get(discr_name)
  N <- length(data)
  
  iter <- floor(N/N_obs)
  for(it in 1:iter){
    X_obs <- data[((it - 1) * N_obs + 1): (it * N_obs)]
    ppc <- as.numeric(compute_ppc(X_obs, para_mod,
                                  dist_fun = ind_dis, discr_fun = discr_fun, R = R))
    spc_0.5 <- as.numeric(compute_singleSPC(X_obs, para_mod,
                                            dist_fun = ind_dis, discr_fun = discr_fun,
                                            spc_prop = 0.5, R = R))
    spc_0.9 <- as.numeric(compute_singleSPC(X_obs, para_mod,
                                            dist_fun = ind_dis, discr_fun = discr_fun,
                                            spc_prop = 0.9, R = R))
    d_spc_0.5_.49 <- as.numeric(compute_divided_SPC(X_obs, para_mod, rate = 0.49,
                                                    dist_fun = ind_dis, discr_fun = discr_fun,
                                                    R = R, spc_prop = 0.5))
    d_spc_0.9_.49 <- as.numeric(compute_divided_SPC(X_obs, para_mod, rate = 0.49,
                                                    dist_fun = ind_dis, discr_fun = discr_fun,
                                                    R = R, spc_prop = 0.9))
    d_spc_0.9_.39 <- as.numeric(compute_divided_SPC(X_obs, para_mod, rate = 0.39,
                                                    dist_fun = ind_dis, discr_fun = discr_fun,
                                                    R = R, spc_prop = 0.9))
    d_spc_0.5_.39 <- as.numeric(compute_divided_SPC(X_obs, para_mod, rate = 0.39,
                                                    dist_fun = ind_dis, discr_fun = discr_fun,
                                                    R = R, spc_prop = 0.5))
    
    
    new_results <- tibble(method = c( 
      "PPC",
      "single 0.5-SPC",
      "single 0.9-SPC",
      "divided 0.5-SPC, k = N^0.49",
      "divided 0.9-SPC, k = N^0.49",
      "divided 0.5-SPC, k = N^0.39",
      "divided 0.9-SPC, k = N^0.39"),
      discr_fun = rep(discr_name, 7),
      num_of_obs = rep(N_obs, 7),
      pvals = c(ppc, spc_0.5, spc_0.9,
                d_spc_0.5_.49,d_spc_0.9_.39,
                d_spc_0.5_.39,d_spc_0.9_.39))
    load(results_path)
    results <- rbind(new_results, results)
    save(results, file = results_path)
  }
  
}



compute_singleSPC <- function(X_obs, para_mod, dist_fun, discr_fun, spc_prop = 0.5, R){ 
  
  set.seed(NULL)
  
  splitted_data <- split_data(X_obs, spc_prop)
  
  Xspc_new <- splitted_data$new
  Xspc_obs <- splitted_data$obs
  N_new <- length(Xspc_new)
  
  # print(d_new)
  posterior_draws <- generate_yreps_given_obs(Xspc_obs, R, N_new)
  X_reps <- posterior_draws$reps
  thetas <- posterior_draws$thetas
  d_new <- rep(NA, R)
  for(r in 1:R){
    d_new[r] <- discr_fun(Xspc_new, thetas[r])
    d_reps <- discr_fun(X_reps[r, ], thetas[r])
  }
  # print(d_reps)
  SPC <- mean(d_reps < d_new)
  
  return(SPC)
  
}




compute_ppc <- function(X_obs, para, dist_fun, discr_fun, R){
  set.seed(NULL)
  N <- length(X_obs)
  PPC <- 0
  posterior_draws <- generate_yreps_given_obs(X_obs, R, N)
  X_reps <- posterior_draws$reps
  thetas <- posterior_draws$thetas
  d_new <- rep(NA, R)
  for(r in 1:R){
    d_new[r] <- discr_fun(X_obs, thetas[r])
    d_reps <- discr_fun(X_reps[r, ], thetas[r])
  }
  
  PPC <- mean(d_reps < d_new)
  return(PPC)
}






compute_divided_SPC <- function(X_obs, para_mod, rate, dist_fun = ind_dis, discr_fun, 
                          R = 500, spc_prop = 0.5, metric = "Kolmogorov"){  # try R = 5000
  
  N <- length(X_obs)
  K <- floor(N^(rate))
  # print(K)
  #  K <- floor(sqrt(N))  # make it a little smaller than sqrt(N); N^(1/4)
  divided_data <- divide_data(X_obs, rate)
  pvals <- rep(0, K)
  for(k in 1:K){
    pvals[k] <- compute_singleSPC(divided_data[[k]], para_mod, dist_fun, discr_fun, spc_prop, R)
    # print(pvals[k])
  }
  
  dspc_pval <- ks.test(pvals + rnorm(K, 0, 0.0001), "punif", 0, 1)$p.value
  return(dspc_pval)
}



generate_yreps_given_obs <- function(X_obs, R, N_new, a0 = 0.1, b0 = 0.2){
  N_obs <- length(X_obs)
  aN <- a0 + N_obs
  bN <- b0 + sum(X_obs)
  theta <- rbeta(R, aN, bN)
  return( list("thetas" = theta, "reps" = matrix(rgeom(R * N_new , theta), R, N_new)))
}



power_plot_airlines <- function(res, discr){
  ts_and_power <- res %>% filter( discr_fun == discr) %>%  
    group_by(num_of_obs, method) %>% 
    mutate(ts_or_power = mean(pvals > 0.975 | pvals < 0.025), iter = floor(N/num_of_obs),
           ci_lower = ts_or_power - sqrt(ts_or_power * (1 - ts_or_power)) * 1.96/ sqrt(iter),
           ci_upper = ts_or_power + sqrt(ts_or_power * (1 - ts_or_power)) * 1.96/ sqrt(iter)) %>% 
    summarize(ts_or_power = mean(ts_or_power), ci_lower = mean(ci_lower), ci_upper = mean(ci_upper)) 
  pd <- position_dodge(0.1)
  logpower_plot_legend <- ts_and_power %>% 
    ggplot(aes(num_of_obs, ts_or_power, group = as.factor(method), 
               color = as.factor(method))) + 
    geom_point(position = pd) + geom_line(size = 0.8,position = pd) + 
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), position = pd, alpha = 0.3)+
    #  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.1, position = pd) +
    labs( y = "power", x = "number of observations", color = "method", linetype = "method")+ 
    scale_color_manual(breaks = c( "PPC", "single 0.9-SPC", "single 0.5-SPC", 
                                   "divided 0.5-SPC, k = N^0.49", 
                                   "divided 0.9-SPC, k = N^0.49",
                                   "divided 0.5-SPC, k = N^0.39", 
                                   "divided 0.9-SPC, k = N^0.39"),
                       values=c("#0000FF","#333333","#333333" ,
                                "#FF9900", "#FF9900", "#33CC99", "#33CC99")) + 
    theme(
      # legend.position = "topleft",
      legend.position = c(0.75,0.2),  # gaussian mse power
      legend.key.size = unit(0.32, "cm"),
      panel.background = element_rect(fill = "white", colour = NA),
      axis.line = element_line(size = 0.5),
      axis.title = element_text(size = 16, color = "black"),
      axis.text = element_text(size = 13, color = "black"),
      axis.ticks = element_blank(),
      legend.text = element_text(size = 12, color = "black"),
      legend.background = element_rect(color = "grey"),
      legend.title = element_blank()
    ) + 
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) + ylim(0,1)
  
  
  logpower_plot <- ts_and_power %>% 
    ggplot(aes(num_of_obs, ts_or_power, group = as.factor(method), 
               color = as.factor(method))) + 
    geom_point(position = pd) + geom_line(size = 0.8, position = pd) + 
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), position = pd, alpha = 0.3)+
    #geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.1, position = pd) +
    labs( y = "power", x = "number of observations", color = "method", linetype = "method")+ 
    scale_color_manual(breaks = c( "PPC", "single 0.9-SPC", "single 0.5-SPC", 
                                   "divided 0.5-SPC, k = N^0.49", 
                                   "divided 0.9-SPC, k = N^0.49",
                                   "divided 0.5-SPC, k = N^0.39", 
                                   "divided 0.9-SPC, k = N^0.39"),
                       values=c("#0000FF","#333333","#333333" ,
                                "#FF9900", "#FF9900", "#33CC99", "#33CC99")) + 
    theme(
      # legend.position = "topleft",
      legend.position = "none",  # gaussian mse power
      legend.key.size = unit(0.32, "cm"),
      panel.background = element_rect(fill = "white", colour = NA),
      axis.line = element_line(size = 0.5),
      axis.title = element_text(size = 16, color = "black"),
      axis.text = element_text(size = 13, color = "black"),
      axis.ticks = element_blank(),
      legend.text = element_text(size = 12, color = "black"),
      legend.background = element_rect(color = "grey"),
      legend.title = element_blank()
    ) + 
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)))+ ylim(0,1)
  
  return(list(logpower_plot_legend,logpower_plot, ts_and_power))
}

