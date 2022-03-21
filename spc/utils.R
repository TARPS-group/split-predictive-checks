

# This script includes statistics functions and
# functions used to compute p-values of PPC, Pop-PC ideal, SPCs 
# with simple exponential family models and two-level hierachical models.

library(rstan)
rstan_options(auto_write=TRUE)
library(loo)
library(twosamples)

# I. Distance functions for p-values
ind_dis <- function(a, b) as.numeric(a > b)


## II. Statistics or Discrepancies
## Simple EF models:
## input: vector X;  output: a real-valued statistic 
emp_mean <- function(X, theta = NULL)  mean(X)
test_max <- function(X, theta = NULL)  max(X)
quantile_0.75 <- function(X, theta = NULL)  quantile(X, probs = 0.75)
quantile_0.95 <- function(X, theta = NULL)  quantile(X, probs = 0.95)
quantile_0.9 <- function(X, theta = NULL) quantile(X, probs = 0.9)
quantile_0.05 <- function(X, theta = NULL) quantile(X, probs = 0.05)
quantile_0.01 <- function(X, theta = NULL) quantile(X, probs = 0.01)
sec_moment <- function(X, theta = NULL) mean(X^2)
third_moment <- function(X, theta = NULL) mean(X^3)
mse <- function(X, theta) mean((X - theta)^2)

## Hierachical models:
## input: matrix X;  output: a real-valued statistic 
max_group_means <- function(X) max(apply(X, 1, mean))
q0.75_group_means <- function(X) quantile((apply(X, 1, mean)), probs = 0.75)
mean_group_q75 <- function(X) mean(apply(X, 1, function(x) quantile(x, probs = 0.75)))
grand_mean <- function(X) mean(X)

## Time series models
## input: vector X; output: vector of (k+1) auto-correlation coefficients 
lag_k <- function(X, k = 1) acf(X, lag.max = k, type = "correlation", plot = FALSE)$acf
success_rate <- function(X, theta = NULL) length(X)/sum(X)




