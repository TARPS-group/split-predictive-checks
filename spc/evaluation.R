# ============================================================
# Analysis & Visualization Utilities
# ============================================================

library(tidyverse)
library(ggthemes)
library(RColorBrewer)
library(scales)

# ============================================================
# 1. Test Size / Power Evaluation (Simple Models)
#
# Input:
#   data        : data.frame containing p-values and metadata
#   discr_name  : discrepancy function used for evaluation
#
# Output:
#   - Log-scale test size plots (with and without legend)
#   - Log-scale power plots (with and without legend)
#   - Summary evaluation table
#
# Rejection rule: p-value < 0.025 or > 0.975
# ============================================================

eval_checks <- function(data, discr_name, model_name = ""){
  
  ts_and_power <- data %>%
    filter(discr_fun == discr_name) %>%
    group_by(size, method, model, discr_fun) %>%
    mutate(iter = n(),
           ts_or_power = mean(pvals > 0.975 | pvals < 0.025)) %>%
    summarize(
      ts_or_power = mean(ts_or_power),
      ci = sqrt(ts_or_power * (1 - ts_or_power)) * 1.96 / sqrt(length(iter))
    ) %>%
    arrange(desc(model), method)
  
  pd <- position_dodge(0.1)
  
  # ----------------------------
  # Power plot (misspecified)
  # ----------------------------
  
  base_plot_pw <- ts_and_power %>%
    filter(model == "misspecified") %>%
    ggplot(aes(size, ts_or_power,
               group = as.factor(method),
               color = as.factor(method),
               linetype = as.factor(method))) +
    geom_point(position = pd) +
    geom_line(position = pd, size = 0.8) +
    ylim(0, 1) +
    geom_errorbar(aes(ymin = ts_or_power - ci,
                      ymax = ts_or_power + ci),
                  width = .1,
                  position = pd) +
    labs(y = "power",
         x = "number of observations",
         color = "method",
         linetype = "method") +
    geom_hline(yintercept = 1,
               linetype = "dashed",
               color = "gray")
  
  power_plot <- base_plot_pw +
    scale_x_log10(
      breaks = trans_breaks("log10", function(x) 10^x),
      labels = trans_format("log10", math_format(10^.x))
    ) +
    theme(
      panel.background = element_rect(fill = "white"),
      axis.line = element_line(size = 0.5),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 13),
      axis.ticks = element_blank(),
      legend.position = "none"
    )
  
  power_plot_legend <- power_plot +
    theme(
      legend.position = c(0.78, 0.28),
      legend.background = element_rect(color = "grey",
                                       fill = "white")
    )
  
  # ----------------------------
  # Test size plot (well-specified)
  # ----------------------------
  
  base_plot_ts <- ts_and_power %>%
    filter(model == "well-specified") %>%
    ggplot(aes(size, ts_or_power,
               color = as.factor(method),
               linetype = as.factor(method),
               group = as.factor(method))) +
    geom_point(position = pd) +
    geom_line(position = pd, size = 0.8) +
    geom_errorbar(aes(ymin = ts_or_power - ci,
                      ymax = ts_or_power + ci),
                  width = .1,
                  position = pd) +
    geom_abline(slope = 0,
                intercept = 0.05,
                linetype = "dashed") +
    labs(y = "test size",
         x = "number of observations")
  
  log_test_size_plot <- base_plot_ts +
    scale_x_log10(
      breaks = trans_breaks("log10", function(x) 10^x),
      labels = trans_format("log10", math_format(10^.x))
    ) +
    theme(
      panel.background = element_rect(fill = "white"),
      axis.line = element_line(size = 0.5),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 13),
      axis.ticks = element_blank(),
      legend.position = "none"
    )
  
  log_test_size_plot_legend <- log_test_size_plot +
    theme(
      legend.position = "bottom",
      legend.background = element_rect(color = "grey",
                                       fill = "white")
    )
  
  return(list(
    log_test_size_plot,
    log_test_size_plot_legend,
    power_plot,
    power_plot_legend,
    evaluations = ts_and_power
  ))
}


# ============================================================
# 2. Q-Q Plot Comparison (Simple Models)
#
# Visual comparison of empirical p-values vs Uniform(0,1)
# under well-specified models.
# ============================================================

qq_plot_comparison <- function(data, discr_name,
                               model_name = "",
                               iter = 1000){
  
  data$method_f <- factor(
    data$method,
    levels = c(
      'PPC',
      'single 0.9-SPC',
      'single 0.5-SPC',
      'divided 0.9-SPC',
      'divided 0.7-SPC',
      'divided 0.5-SPC',
      'divided 0.3-SPC',
      'divided 0.1-SPC'
    )
  )
  
  base_plot <- data %>%
    filter(model == "well-specified",
           discr_fun == discr_name) %>%
    group_by(size, method) %>%
    mutate(
      unifcdf = 1:iter / (iter + 1),
      order_pvals = sort(pvals)
    ) %>%
    ggplot(aes(unifcdf,
               order_pvals,
               color = as.factor(size))) +
    geom_line(size = 0.5) +
    geom_abline(slope = 1,
                intercept = 0,
                color = "grey") +
    labs(color = "sample size",
         y = "p-values") +
    ylim(0, 1) +
    xlim(0, 1) +
    facet_wrap(~as.factor(method_f))
  
  qq_plot <- base_plot +
    theme(
      legend.position = "none",
      panel.background = element_rect(fill = "white"),
      axis.line = element_line(size = 0.25),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 12),
      axis.ticks = element_blank()
    )
  
  qq_plot_legend <- base_plot +
    theme(
      panel.background = element_rect(fill = "white"),
      axis.line = element_line(size = 0.25),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 12),
      legend.background = element_rect(color = "grey",
                                       fill = "white")
    )
  
  return(list(qq_plot, qq_plot_legend))
}
