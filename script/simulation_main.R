
# This script includes simulation-models study for 
# conjugate Poisson, conjugate Gaussian and two-level Gaussian hierarchical models.

library(tidyverse)
library(rstan)
library("rprojroot")
rstan_options(auto_write=TRUE)
root <- has_file("README.md")$make_fix_file()
source(root("spc/codes", "utils.r"))
source(root("spc/codes", "generate_pvals.r"))
source(root("spc/codes", "compute_spcs.r"))
source(root("spc/codes", "compute_spcs_hierarchical.r"))

####### Initialize data file paths to store results ####### 
updated_data <- NULL
save(updated_data, file = "poisson_model.rda")
save(updated_data, file = "poisson_model_ess.rda")
save(updated_data, file = "poisson_mismatch_comparison.rda")
save(updated_data, file = "gaussian_model.rda")
save(updated_data, file = "test.rda")


#######  I. Poisson model ####### 

source(root("spc/simulation-models", "Poisson_model.r"))
load("poisson_model.rda")
start <- Sys.time()
sizes <- c(50, 100, 1000, 3000, 5000)
run_poisson_model(sizes, data_path = "poisson_model.rda", discr_name = "emp_mean", iter = 1000)
run_poisson_model(sizes, data_path = "poisson_model.rda", discr_name = "sec_moment", iter = 1000)
run_poisson_model(sizes, data_path = "poisson_model.rda", discr_name = "third_moment", iter = 1000)
run_poisson_model(sizes, data_path = "poisson_model.rda", discr_name = "mse", iter = 1000)
end <- Sys.time()
end - start


source(root("spc/simulation-models", "Poisson_model.r"))
load("poisson_mismatch_comparison.rda")

sizes <- c(1000, 10000, 50000)
paras <- list("s" = 2, "t" = c(0.01, 0.1, 0.5)) # over-dispersed
run_mismatch_comparison_pois(sizes, discr_name = "emp_mean", dispersion_type = "over-dispersed", true_paras = paras)
run_mismatch_comparison_pois(sizes, discr_name = "sec_moment", dispersion_type = "over-dispersed", true_paras = paras)
paras <- list("s" = 30, "t" = c(0.1, 0.5, 0.8)) # under-dispersed
run_mismatch_comparison_pois(sizes, discr_name = "emp_mean", dispersion_type = "under-dispersed", true_paras = paras)
run_mismatch_comparison_pois(sizes, discr_name = "sec_moment", dispersion_type = "under-dispersed", true_paras = paras)



source(root("spc/simulation-models","Poisson_model.r"))
r_adj <- c(seq(0.5, 0.05, by = -0.05), 0.01, 0.005, 0.003, 0.001)
run_multiple_experiments_ess_pois(r = r_adj, discr_fun = "emp_mean", true_theta = 25, quantile = 0.95, data_path = "poisson_model_ess.rda")
run_multiple_experiments_ess_pois(r = r_adj, discr_fun = "emp_mean", true_theta = 25, quantile = 0.99, data_path = "poisson_model_ess.rda")
run_multiple_experiments_ess_pois(r = r_adj, discr_fun = "sec_moment", true_theta = 25, quantile = 0.95, data_path = "poisson_model_ess.rda")
run_multiple_experiments_ess_pois(r = r_adj, discr_fun = "sec_moment", true_theta = 25, quantile = 0.99, data_path = "poisson_model_ess.rda")



####### II. Gaussian location model (Appendix) ####### 

source(root("spc/simulation-models","Gaussian_model.r"))
load("gaussian_model.rda")
start <- Sys.time()
sizes <- c(50, 100, 1000, 3000, 5000)
run_gaussian_model(sizes, data_path = "gaussian_model.rda", discr_name = "emp_mean", iter = 1000)
run_gaussian_model(sizes, data_path = "gaussian_model.rda", discr_name = "sec_moment", iter = 1000)
run_gaussian_model(sizes, data_path = "gaussian_model.rda", discr_name = "quantile_0.75", iter = 1000)
run_gaussian_model(sizes, data_path = "gaussian_model.rda", discr_name = "mse", iter = 1000)
end <- Sys.time()
end - start




####### III. Two-level bierachical gaussian model ####### 


## Scenario 1: Xij ~ N(theta_i, 4), i = 1,...,I   j = 1, ..., J
##              theta_i ~ N(0, 1)
## We fix J = 8, run for I = 20, 50, 100, 150, 200;
## We fix I = 20, run for J = 8, 50, 100, 150, 200; 
source(root("spc/simulation-models","hier_model_cross_group.r"))
start <- Sys.time()
I <- 20
J <- 8
load("hier_well_specified.rda")
run_experiment_hier(I, J, iter = 200, data = updated_data, discr_fun_name = "q0.75_group_means",
                    file_path = "hier_well_specified.rda", misspecified = FALSE, R = 500)
load("hier_well_specified.rda")
run_experiment_hier(I, J, iter = 200, data = updated_data, discr_fun_name = "grand_mean",
                    file_path = "hier_well_specified.rda", misspecified = FALSE, R = 500)
load("hier_well_specified.rda")
run_experiment_hier(I, J, iter = 200, data = updated_data, discr_fun_name = "mean_group_q75",
                    file_path = "hier_well_specified.rda", misspecified = FALSE, R = 500)
end <- Sys.time()
end - start



##  Scenario 2: Xij ~ N(theta_i, 4), i = 1,...,I   j = 1, ..., J
##               theta_i ~ gamma(0.6, 0.2)
## We fix J = 8, run for I = 20, 50, 100, 150, 200;
## We fix I = 20, run for J = 8, 50, 100, 150, 200; 


source(root("spc/simulation-models", "hier_model_cross_group.r"))
start <- Sys.time()
I <- 20
J <- 8
load("hier_cr_jeff.rda")
run_experiment_hier(I, J, iter = 200, data = updated_data, discr_fun_name = "q0.75_group_means",
                            file_path = "hier_cr_jeff.rda", misspecified = TRUE, R = 500)
load("hier_cr_jeff.rda")
run_experiment_hier(I, J, iter = 200, data = updated_data, discr_fun_name = "grand_mean",
                            file_path = "hier_cr_jeff.rda", misspecified = TRUE, R = 500)
load("hier_cr_jeff.rda")
run_experiment_hier(I, J, iter = 200, data = updated_data, discr_fun_name = "mean_group_q75",
                            file_path = "hier_cr_jeff.rda", misspecified = TRUE, R = 500)
end <- Sys.time()
end - start



##  Scenario 3: Xij ~ N(theta_i, 8), i = 1,...,I   j = 1, ..., J
##              theta_i ~ N(0,1)
## We fix J = 8, run for I = 20, 50, 100, 150, 200;
## We fix I = 20, run for J = 8, 50, 100, 150, 200; 

source(root("spc/simulation-models", "hier_model_wi_gr_var.r"))
start <- Sys.time()
I <- 20
J <- 8
load("hier_wi_gr_var.rda")
run_experiment_hier(I, J, iter = 200, data = updated_data, discr_fun_name = "q0.75_group_means",
                    file_path = "hier_wi_gr_var.rda", misspecified = TRUE, R = 500)
load("hier_wi_gr_var.rda")
run_experiment_hier(I, J, iter = 200, data = updated_data, discr_fun_name = "grand_mean",
                    file_path = "hier_wi_gr_var.rda", misspecified = TRUE, R = 500)
load("hier_wi_gr_var.rda")
run_experiment_hier(I, J, iter = 200, data = updated_data, discr_fun_name = "mean_group_q75",
                    file_path = "hier_wi_gr_var.rda", misspecified = TRUE, R = 500)
end <- Sys.time()
end - start



##  Scenario 4: Xij ~ log-normal(theta_i, 4), i = 1,...,I   j = 1, ..., J
##              theta_i ~ N(0,1)
## We fix J = 8, run for I = 20, 50, 100, 150, 200;
## We fix I = 20, run for J = 8, 50, 100, 150, 200; 

source(root("spc/simulation-models","hier_model_wi_lognorm.r"))
start <- Sys.time()
I <- 20
J <- 8
load("hier_wi_lognorm.rda")
run_experiment_hier(I, J, iter = 200, data = updated_data, discr_fun_name = "q0.75_group_means",
                    file_path = "hier_wi_lognorm.rda", misspecified = TRUE, R = 500)
load("hier_wi_lognorm.rda")
run_experiment_hier(I, J, iter = 200, data = updated_data, discr_fun_name = "grand_mean",
                    file_path = "hier_wi_lognorm.rda", misspecified = TRUE, R = 500)
load("hier_wi_lognorm.rda")
run_experiment_hier(I, J, iter = 200, data = updated_data, discr_fun_name = "mean_group_q75",
                    file_path = "hier_wi_lognorm.rda", misspecified = TRUE, R = 500)
end <- Sys.time()
end - start
