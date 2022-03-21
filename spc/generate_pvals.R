
# This script computes p-values using POP_PC, PPC for different models
# All the results are automatically saved to .rdata file

library(tidyverse)
library(rstan)
library(stats)
rstan_options(auto_write=TRUE)

##### I. Functions to compute p-values for single EF models

# This function produces the p-values of PPC, Pop-PC, single and divided SPCs.
# Input: 
# - N: sample size of data
# - iter: number of repeated experiments to generate a p-value
# - data: data frame used to store the results
# - discr_fun_name: statistics
# - file_path: 

# About misspecified/well-specified case: 
# fix the assumed model and change generative model to get misspecified case

run_experiment <- function(N, iter, data = NULL, discr_fun_name, file_path, misspecified = FALSE){
   function(...){
      discr_fun <- get(discr_fun_name)
      if(misspecified){   
         simulate_data <- simulate_data_mis
         model_type = "misspecified"
      }else{
         simulate_data <- simulate_data_well
         model_type = "well-specified"
      }
      
         pop_pc <- rep(0, iter)
         ppc <- rep(0, iter)
         
         spc_0.1 <- rep(0, iter)
         spc_0.3 <- rep(0, iter)
         spc_0.7 <- rep(0, iter)
         spc_0.5 <- rep(0, iter)
         spc_0.9 <- rep(0, iter)
         
         d_spc_0.5_.4 <- rep(0, iter)
         d_spc_0.5_.49 <- rep(0, iter)
         d_spc_0.5_.6 <- rep(0, iter)
         d_spc_0.5_.8 <- rep(0, iter)

         d_spc_0.1 <- rep(0, iter)
         d_spc_0.3 <- rep(0, iter)
         d_spc_0.5 <- rep(0, iter)
         d_spc_0.7 <- rep(0, iter)
         d_spc_0.9 <- rep(0, iter)

         for(i in 1:iter){
            obs_data <- simulate_data(N)
            X_obs <- obs_data$X_obs
            
            pop_pc[i] <- as.numeric(compute_pop_pc_ideal(X_obs, para_mod,dist_fun = ind_dis,
                                                         discr_fun = discr_fun, R = 500)(...))
            ppc[i] <- as.numeric(compute_ppc(X_obs, para_mod, dist_fun = ind_dis,
                                             discr_fun = discr_fun, R = 500)(...))
            
            spc_0.1[i] <- as.numeric(compute_singleSPC(X_obs, para_mod, dist_fun = ind_dis,
                                                       discr_fun = discr_fun, R = 500, spc_prop = 0.1)(...))
            spc_0.3[i] <- as.numeric(compute_singleSPC(X_obs, para_mod, dist_fun = ind_dis,
                                                       discr_fun = discr_fun, R = 500, spc_prop = 0.3)(...))
            spc_0.5[i] <- as.numeric(compute_singleSPC(X_obs, para_mod, dist_fun = ind_dis,
                                                       discr_fun = discr_fun, R = 500, spc_prop = 0.5)(...))
            spc_0.7[i] <- as.numeric(compute_singleSPC(X_obs, para_mod, dist_fun = ind_dis,
                                                       discr_fun = discr_fun, R = 500, spc_prop = 0.7)(...))
            spc_0.9[i] <- as.numeric(compute_singleSPC(X_obs, para_mod, dist_fun = ind_dis,
                                                       discr_fun = discr_fun, R = 500, spc_prop = 0.9)(...))
           
            d_spc_0.5_.4[i] <- as.numeric(compute_divided_SPC(X_obs, para_mod, rate = 0.4, dist_fun = ind_dis, discr_fun = discr_fun,
                                                         R = 500, spc_prop = 0.5)(...))
            d_spc_0.5_.49[i] <- as.numeric(compute_divided_SPC(X_obs, para_mod, rate = 0.49, dist_fun = ind_dis, discr_fun = discr_fun,
                                                        R = 500, spc_prop = 0.5)(...))
            d_spc_0.5_.6[i] <- as.numeric(compute_divided_SPC(X_obs, para_mod, rate = 0.6, dist_fun = ind_dis, discr_fun = discr_fun,
                                                         R = 500, spc_prop = 0.5)(...))
            d_spc_0.5_.8[i] <- as.numeric(compute_divided_SPC(X_obs, para_mod, rate = 0.8, dist_fun = ind_dis, discr_fun = discr_fun,
                                                         R = 500, spc_prop = 0.5)(...))
            
            d_spc_0.1[i] <- as.numeric(compute_divided_SPC(X_obs, para_mod,  rate = 0.49, dist_fun = ind_dis, discr_fun = discr_fun,
                                                     R = 500, spc_prop = 0.1)(...))
            d_spc_0.3[i] <- as.numeric(compute_divided_SPC(X_obs, para_mod,  rate = 0.49, dist_fun = ind_dis, discr_fun = discr_fun,
                                                     R = 500, spc_prop = 0.3)(...))
            d_spc_0.5[i] <- as.numeric(compute_divided_SPC(X_obs, para_mod,  rate = 0.49, dist_fun = ind_dis, discr_fun = discr_fun,
                                                     R = 500, spc_prop = 0.5)(...))
            d_spc_0.7[i] <- as.numeric(compute_divided_SPC(X_obs, para_mod, rate = 0.49, dist_fun = ind_dis, discr_fun = discr_fun,
                                                     R = 500, spc_prop = 0.7)(...))
            d_spc_0.9[i] <- as.numeric(compute_divided_SPC(X_obs, para_mod, rate = 0.49, dist_fun = ind_dis, discr_fun = discr_fun,
                                                     R = 500, spc_prop = 0.9)(...))

         }
         
         size_of_data <- 16 * iter
         new_data_frame <- tibble(size = rep(N, size_of_data),
                                  method = c(rep("PPC", iter), rep("Pop-PC ideal", iter), 
                                             rep("single 0.1-SPC", iter), rep("single 0.3-SPC", iter), rep("single 0.5-SPC", iter), rep("single 0.7-SPC", iter), rep("single 0.9-SPC", iter),
                                             rep("divided 0.1-SPC", iter), rep("divided 0.3-SPC", iter), rep("divided 0.5-SPC", iter), rep("divided 0.7-SPC", iter), rep("divided 0.9-SPC", iter),
                                             rep("divided 0.5-SPC, k = N^(.4)", iter), rep("divided 0.5-SPC, k = N^(.49)", iter),
                                             rep("divided 0.5-SPC, k = N^(.6)", iter), rep("divided 0.5-SPC, k = N^(.8)", iter)),
                                   model = rep(model_type, size_of_data), 
                                   discr_fun = rep(discr_fun_name, size_of_data), 
                                   pvals = c(ppc, pop_pc, spc_0.1, spc_0.3, spc_0.5, spc_0.7, spc_0.9,
                                             d_spc_0.1, d_spc_0.3, d_spc_0.5, d_spc_0.7, d_spc_0.9,
                                             d_spc_0.5_.4, d_spc_0.5_.49, d_spc_0.5_.6, d_spc_0.5_.8))
   
         updated_data <- rbind(new_data_frame, data)
         save(updated_data,  file = file_path)
         return(updated_data)
      }
   
}


# Compute p-values given over- and under-dispersed misspecification
# which creates different values of rho
# over-dispersion: "NegBin(2, 0.01)"-6, "NegBin(2, 0.1)"-5, "NegBin(2, 0.5)"-4
# under-dispersion: "Binom(30, 0.8)"-1, "Binom(30, 0.5)"-2, "Binom(30, 0.1)"-3


run_experiment_rho <- function(N, data = NULL, discr_fun_name, file_name, dispersion_type, para, iter){
   function(...){
      discr_fun <- get(discr_fun_name)
      if(dispersion_type == "over-dispersed"){
         simulate_data <- function(sample_size){
            X <- rnegbin(sample_size, mu = para$s, theta = para$t)
            return(X)
         }
         rho_squared <- as.character(para$s + para$s^2/para$t)
         
      }else if (dispersion_type == "under-dispersed"){
         simulate_data <- function(sample_size){
            X <- rbinom(sample_size, size = para$s, prob = para$t)
            return(X)
         }
         rho_squared <- as.character(1 - para$t)
      }
      
   
   ppc <- rep(0, iter)
   pop_pc <- rep(0, iter)
   spc_0.5 <- rep(0, iter)
   dspc_0.5 <- rep(0, iter)
   
   for(i in 1:iter){
      X_obs <- simulate_data(N)
 
      pop_pc[i] <- as.numeric(compute_pop_pc_ideal(X_obs, para_mod,dist_fun = ind_dis,
                                                   discr_fun = discr_fun, R = 500)(...))
      ppc[i] <- as.numeric(compute_ppc(X_obs, para_mod, dist_fun = ind_dis,
                                       discr_fun = discr_fun, R = 500)(...))
      spc_0.5[i] <- as.numeric(compute_singleSPC(X_obs, para_mod, dist_fun = ind_dis, discr_fun = discr_fun,
                                                 R = 500, spc_prop = 0.5)(...))
      d_spc_0.5[i] <- as.numeric(compute_divided_SPC(X_obs, para_mod, dist_fun = ind_dis, discr_fun = discr_fun,
                                                     R = 500, spc_prop = 0.5, rate = 0.49)(...))
      
   }
   
   size_of_data <- 4 * iter
   new_data_frame <- tibble(size = rep(N, size_of_data),
                            method = c(rep("PPC", iter), rep("Pop-PC ideal", iter), 
                                       rep("single 0.5-SPC", iter), rep("divided 0.5-SPC", iter)),
                            discr_fun = rep(discr_fun_name, size_of_data), 
                            pvals =  c(ppc, pop_pc, spc_0.5, d_spc_0.5),
                            rho_squared = rep(rho_squared, size_of_data)) 
   
   updated_data <- rbind(new_data_frame, data)
   save(updated_data,  file = file_name)
   return(updated_data)
   }
}

run_mismatch_comparison_pois <- function(sizes, discr_name, dispersion_type, true_paras, iter = 200){
   for(i in 1: length(sizes)){
      for(j in 1:3){
         load(data_path)
         para = list("s" = 2, "t" = true_paras$t[j])
         run_experiment_rho(sizes[i], data = updated_data, discr_fun_name = discr_name, 
                            file_name = data_path, dispersion_type, para = para, iter)
      }
   }
   
   load(data_path)
   updated_data$group[updated_data$rho_squared == "0.2"] <- "1"
   updated_data$group[updated_data$rho_squared == "0.5"] <- "2"
   updated_data$group[updated_data$rho_squared == "0.9"] <- "3"
   updated_data$group[updated_data$rho_squared == "5"] <- "4"
   updated_data$group[updated_data$rho_squared == "21"] <- "5"
   updated_data$group[updated_data$rho_squared == "201"] <- "6"
   save(updated_data, file = data_path)
}



# Compute prior ess for poisson models:
# Fix N_obs = 50 and true value theta_o = 25
# Vary r to create different beta's
# Compute alpha such that theta_o lies at 95th quantile of the prior distribution



run_single_experiment_ess_pois <- function(N, iter, data = NULL, r_theory, theta_star, quantile = 0.99,  
                                    discr_fun_name, file_name){
   discr_fun <- get(discr_fun_name)
  
   simulate_data <- function(sample_size, true_theta){
      X <- rpois(sample_size, theta_star)
      return(X)
   }
   
  
   beta <- r_theory/(1 - r_theory) * N
   alpha = optimize(function(a) (quantile - pgamma(theta_star, a, beta))^2, 
                    lower = 0.01, upper = theta_star*beta, 
                    maximum = FALSE, tol = 1e-06)$minimum
   
   ppc <- rep(0, iter)
   pop_pc <- rep(0, iter)
   spc_0.5 <- rep(0, iter)
   dspc_0.5 <- rep(0, iter)
   
   for(i in 1:iter){
      X_obs_Nppc <- simulate_data(N, theta_star)
      
      Nspc <- 2*N
      X_obs_Nspc <- simulate_data(Nspc, theta_star)
      
      Ndspc <- floor((2*N)^(1/0.51))
      X_obs_Ndspc <- simulate_data(Ndspc, theta_star)
      
      pop_pc[i] <- as.numeric(compute_pop_pc_ideal(X_obs_Nppc, para_mod,dist_fun = ind_dis,
                                                   discr_fun = discr_fun, R = 500)(a0 = alpha, b0 = beta))
      ppc[i] <- as.numeric(compute_ppc(X_obs_Nppc, para_mod, dist_fun = ind_dis,
                                       discr_fun = discr_fun, R = 500)(a0 = alpha, b0 = beta))
      spc_0.5[i] <- as.numeric(compute_singleSPC(X_obs_Nspc, para_mod, dist_fun = ind_dis, discr_fun = discr_fun,
                                                 R = 500, spc_prop = 0.5)(a0 = alpha,b0 = beta))
      d_spc_0.5[i] <- as.numeric(compute_divided_SPC(X_obs_Ndspc, para_mod, dist_fun = ind_dis, discr_fun = discr_fun,
                                               R = 500, spc_prop = 0.5, rate = 0.49)(a0 = alpha, b0 = beta))
      
   }
   
   size_of_data <- 4 * iter
   new_data_frame <- tibble(size = rep(N, size_of_data),
                            prior_beta = rep(beta, size_of_data),
                            prior_beta = rep(alpha, size_of_data),
                            method = c(rep("PPC", iter), rep("Pop-PC ideal", iter), rep("single 0.5-SPC", iter), rep("divided 0.5-SPC", iter)),
                            discr_fun = rep(discr_fun_name, size_of_data), 
                            pvals =  c(ppc, pop_pc, spc_0.5, d_spc_0.5),
                            true_theta = rep(theta_star, size_of_data),
                            quantile = rep(quantile, size_of_data))
   
   updated_data <- rbind(new_data_frame, data)
   save(updated_data,  file = file_name)
   return(updated_data)
}



run_multiple_experiments_ess_pois <- function(r, discr_fun, true_theta, quantile, data_path){
   for(i in 1:length(r_adj)){
      load(data_path)
      run_single_experiment_ess_pois(N = 50, iter = 200, data = updated_data, r_theory = r_adj[i], theta_star = true_theta, 
                              quantile = quantile,  discr_fun_name = discr_fun, 
                              file_name = data_path)
   }
   
}

##### II. Functions to compute p-values for hierarchical models

run_experiment_hier <- function(I, J, iter, data = NULL, discr_fun_name, file_path, misspecified = FALSE, R = 500){
      
      discr_fun <- get(discr_fun_name)

      if(misspecified){   
         simulate_data <- simulate_data_mis
         model_type = "misspecified"
      }else{
         simulate_data <- simulate_data_well
         model_type = "well-specified"
      }

         pop_pc <- rep(0, iter)
         ppc <- rep(0, iter)
         spc_0.5cross <- rep(0, iter)
         spc_0.5within <- rep(0, iter)
         crdiv_0.5_crSPC <- rep(0, iter)
         widiv_0.5_crSPC <- rep(0, iter)
         widiv_0.5_wiSPC <- rep(0, iter)
         crdiv_0.5_wiSPC <- rep(0, iter)

         for(i in 1:iter){
            obs_data <- simulate_data(I, J)
            X_obs <- obs_data$X_obs
             
            pop_pc[i] <- as.numeric(compute_pop_pc_ideal_hier(X_obs, para_mod, dist_fun = ind_dis, discr_fun = discr_fun, R = R))
            ppc[i] <- as.numeric(compute_ppc_hier(X_obs, para_mod,dist_fun = ind_dis,  discr_fun = discr_fun, R = R))
           
            spc_0.5cross[i] <- as.numeric(compute_singleSPC_hier(X_obs, para_mod, discr_fun,dist_fun = ind_dis, spc_prop = 0.5, R = R, cross_split = TRUE))
            spc_0.5within[i] <- as.numeric(compute_singleSPC_hier(X_obs, para_mod, discr_fun,dist_fun = ind_dis, spc_prop = 0.5, R = R, cross_split = FALSE))
            
            crdiv_0.5_crSPC[i] <- compute_dividedSPC_hier(X_obs, para_mod, dist_fun = ind_dis, discr_fun = discr_fun, R = R, 
                                                          spc_prop = 0.5,  cross_divide = TRUE, cross_split = TRUE)
            widiv_0.5_crSPC[i] <- compute_dividedSPC_hier(X_obs, para_mod, dist_fun = ind_dis, discr_fun = discr_fun, R = R, 
                                                          spc_prop = 0.5,  cross_divide = FALSE, cross_split = TRUE)
            widiv_0.5_wiSPC[i] <- compute_dividedSPC_hier(X_obs, para_mod, dist_fun = ind_dis, discr_fun = discr_fun, R = R, 
                                                          spc_prop = 0.5,  cross_divide = FALSE, cross_split = FALSE)
            crdiv_0.5_wiSPC[i] <- compute_dividedSPC_hier(X_obs, para_mod, dist_fun = ind_dis, discr_fun = discr_fun, R = R, 
                                                          spc_prop = 0.5,  cross_divide = TRUE, cross_split = FALSE)
            
         }
         
         size_of_data <- 8 * iter
         new_data_frame <- tibble( num_of_groups = rep(I, size_of_data),
                                   num_of_indi = rep(J, size_of_data),
                                   method = c(rep("Pop-PC", iter), rep("PPC", iter), 
                                              rep("single 0.5-crossSPC", iter), rep("single 0.5-withinSPC", iter),
                                              rep("cross-divided 0.5-crossSPC", iter), 
                                              rep("within-divided 0.5-crossSPC", iter), 
                                              rep("within-divided 0.5-withinSPC", iter), 
                                              rep("cross-divided 0.5-withinSPC", iter)), 
                                   model = rep(model_type, size_of_data), 
                                   discr_fun = rep(discr_fun_name, size_of_data), 
                                   pvals = c(pop_pc, ppc, spc_0.5cross, spc_0.5within,
                                             crdiv_0.5_crSPC, widiv_0.5_crSPC, widiv_0.5_wiSPC , crdiv_0.5_wiSPC))
         
         updated_data <- rbind(new_data_frame, data)
         save(updated_data, file = file_path)
         return(updated_data)
   }



   



