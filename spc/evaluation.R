# Analysis & visualization

library(tidyverse)
library(ggthemes)
library("RColorBrewer")
library(scales)
## 1. Compute test size/ power for checks
## eval_checks: simple models
## Input: data.frame, iter(# of iterations for each setting)
##        specify one discr_fun used only for plots
## Output: test size and power under each setting
##         plot test size/ power versus sample size for each setting provided a specific discr_fun
## discr_names now include "emp_mean", "quantile_0.75", "sec_moment", "third_moment", "mse"
##


eval_checks <- function(data, discr_name, model_name = ""){
  
  ts_and_power <- data %>% filter( discr_fun == discr_name) %>% 
      group_by(size, method, model, discr_fun) %>% 
      mutate(iter = n(), ts_or_power = mean(pvals > 0.975 | pvals < 0.025)) %>% 
      summarize(ts_or_power = mean(ts_or_power), 
                ci =  sqrt(ts_or_power * (1 - ts_or_power)) * 1.96/ sqrt(length(iter))  )   %>% 
      arrange(desc(model), method)
 
  
  pd <- position_dodge(0.1) # move them .05 to the left and right
  
  base_plot_pw <- ts_and_power %>% filter(model == "misspecified", discr_fun == discr_name) %>% 
    ggplot(aes(size, ts_or_power, group = as.factor(method), color = as.factor(method), linetype = as.factor(method))) + 
    geom_point(position = pd) + geom_line(position = pd, size = 0.8)  + ylim(0,1) +
    geom_errorbar(aes(ymin=ts_or_power - ci, ymax=ts_or_power + ci), width=.1, position=pd) +
    labs( y = "power", x = "number of observations", color = "method", linetype = "method")+ 
    scale_color_manual(breaks = c("Pop-PC ideal", "PPC", 
                                  "single 0.1-SPC", "single 0.3-SPC", "single 0.5-SPC","single 0.7-SPC", "single 0.9-SPC", 
                                  "divided 0.1-SPC", "divided 0.3-SPC", "divided 0.5-SPC","divided 0.7-SPC", "divided 0.9-SPC", 
                                  "divided 0.9-SPC, k = N^(.4)", "divided 0.5-SPC, k = N^(.4)",
                                  "divided 0.9-SPC, k = N^(.49)", "divided 0.5-SPC, k = N^(.49)",
                                  "divided 0.9-SPC, k = N^(.6)", "divided 0.5-SPC, k = N^(.6)",
                                  "divided 0.9-SPC, k = N^(.8)", "divided 0.5-SPC, k = N^(.8)"
    ),
    values=c("#666666","#0000FF" , 
             "#009999", "#FF9900","#D55E00" , "#6699FF", "#666666",
             "#6699FF", "#FF9900", "#009999", "#666666", "#D55E00",
             "#666666","#666666",
             "#009E73","#009E73",
             "#6633FF","#6633FF",
             "#6699FF", "#6699FF"
    )) + 
    scale_linetype_manual(breaks = c("Pop-PC ideal", "PPC", 
                                     "single 0.1-SPC", "single 0.3-SPC", "single 0.5-SPC","single 0.7-SPC", "single 0.9-SPC", 
                                     "divided 0.1-SPC", "divided 0.3-SPC", "divided 0.5-SPC","divided 0.7-SPC", "divided 0.9-SPC", 
                                     "divided 0.9-SPC, k = N^(.4)", "divided 0.5-SPC, k = N^(.4)",
                                     "divided 0.9-SPC, k = N^(.49)", "divided 0.5-SPC, k = N^(.49)",
                                     "divided 0.9-SPC, k = N^(.6)", "divided 0.5-SPC, k = N^(.6)",
                                     "divided 0.9-SPC, k = N^(.8)", "divided 0.5-SPC, k = N^(.8)"
    ),
    values = c("solid", "solid", 
               "solid", "solid", "solid", "solid","solid",
               "solid", "solid", "solid", "solid","solid",
               "solid", "solid", "solid", "solid",
               "solid", "solid", "solid", "solid")) + geom_hline(yintercept = 1, linetype = "dashed", color = "gray")
  
  
  power_plot_legend <- base_plot_pw +
    theme(
      legend.position = c(0.78, 0.28), 
      legend.key.size = unit(0.32, "cm"),
      panel.background = element_rect(fill = "white", colour = NA),
      axis.line = element_line(size = 0.5),
      axis.title = element_text(size = 16, color = "black"),
      axis.text = element_text(size = 13, color = "black"),
      axis.ticks = element_blank(),
      legend.text = element_text(size = 12, color = "black"),
      legend.background = element_rect(color = "grey", fill = "white"),
      legend.title = element_blank(),
      legend.key=element_blank()
    ) 

  power_plot <- base_plot_pw + theme(
    panel.background = element_rect(fill = "white", colour = NA),
    axis.line = element_line(size = 0.5),
    axis.title = element_text(size = 16, color = "black"),
    axis.text = element_text(size = 13, color = "black"),
    axis.ticks = element_blank(),
    legend.position = "none"
  )
  
  log_power_plot_legend <-  power_plot_legend + 
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)))
  log_power_plot <-  power_plot + 
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)))
  
  
  base_plot_ts <- ts_and_power %>% filter(model == "well-specified", discr_fun == discr_name) %>% 
    ggplot(aes(size, ts_or_power, color = as.factor(method), linetype = as.factor(method), group = as.factor(method))) + 
    geom_point(position = pd) + geom_line(position = pd, size = 0.8)  +
    geom_errorbar(aes(ymin=ts_or_power - ci, ymax=ts_or_power + ci), width=.1, position=pd) +
    labs(y = "test size", x = "number of observations", color = "method", linetype = "method") + 
    geom_abline(slope = 0, intercept = 0.05, linetype = "dashed") +
    scale_color_manual(breaks = c("Pop-PC ideal", "PPC", 
                                  "single 0.1-SPC", "single 0.3-SPC", "single 0.5-SPC","single 0.7-SPC", "single 0.9-SPC", 
                                  "divided 0.1-SPC", "divided 0.3-SPC", "divided 0.5-SPC","divided 0.7-SPC", "divided 0.9-SPC", 
                                  "divided 0.9-SPC, k = N^(.4)", "divided 0.5-SPC, k = N^(.4)",
                                  "divided 0.9-SPC, k = N^(.49)", "divided 0.5-SPC, k = N^(.49)",
                                  "divided 0.9-SPC, k = N^(.6)", "divided 0.5-SPC, k = N^(.6)",
                                  "divided 0.9-SPC, k = N^(.8)", "divided 0.5-SPC, k = N^(.8)"
    ),
    values=c("#666666","#0000FF" , 
             "#009999", "#FF9900","#D55E00" , "#6699FF", "#666666",
             "#6699FF", "#FF9900", "#009999", "#666666", "#D55E00",
             "#666666","#666666",
             "#009E73","#009E73",
             "#6633FF","#6633FF",
             "#6699FF", "#6699FF"
    )) + 
    scale_linetype_manual(breaks = c("Pop-PC ideal", "PPC", 
                                     "single 0.1-SPC", "single 0.3-SPC", "single 0.5-SPC","single 0.7-SPC", "single 0.9-SPC", 
                                     "divided 0.1-SPC", "divided 0.3-SPC", "divided 0.5-SPC","divided 0.7-SPC", "divided 0.9-SPC", 
                                     "divided 0.9-SPC, k = N^(.4)", "divided 0.5-SPC, k = N^(.4)",
                                     "divided 0.9-SPC, k = N^(.49)", "divided 0.5-SPC, k = N^(.49)",
                                     "divided 0.9-SPC, k = N^(.6)", "divided 0.5-SPC, k = N^(.6)",
                                     "divided 0.9-SPC, k = N^(.8)", "divided 0.5-SPC, k = N^(.8)"
    ),
    values = c("solid", "solid", 
               "solid", "solid", "solid", "solid","solid",
               "solid", "solid", "solid", "solid","solid",
               "solid", "solid", "solid", "solid",
               "solid", "solid", "solid", "solid")) 
  
  
  log_test_size_plot_legend <- base_plot_ts +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                  labels = trans_format("log10", math_format(10^.x))) +
    theme(
      legend.position = "bottom",
      # legend.position = c(0.85, 0.6),  # gaussian mse power
      legend.key.size = unit(0.32, "cm"),
      panel.background = element_rect(fill = "white", colour = NA),
      axis.line = element_line(size = 0.5),
      axis.title = element_text(size = 16, color = "black"),
      axis.text = element_text(size = 13, color = "black"),
      axis.ticks = element_blank(),
      legend.text = element_text(size = 12, color = "black"),
      legend.background = element_rect(color = "grey", fill = "white"),
      legend.title = element_blank(),
      legend.key=element_blank()
    ) + guides(color=guide_legend(nrow=2,byrow=TRUE))
  
  log_test_size_plot <- base_plot_ts +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                  labels = trans_format("log10", math_format(10^.x))) + 
    theme(
      legend.position = "none",
      legend.key.size = unit(0.32, "cm"),
      panel.background = element_rect(fill = "white", colour = NA),
      axis.line = element_line(size = 0.5),
      axis.title = element_text(size = 16, color = "black"),
      axis.text = element_text(size = 13, color = "black"),
      axis.ticks = element_blank(),
      legend.text = element_text(size = 12, color = "black"),
      legend.background = element_rect(color = "grey", fill = "white"),
      legend.title = element_blank(),
      legend.key=element_blank()
    )
  
  power_table <- ts_and_power 
  return(list( log_test_size_plot, log_test_size_plot_legend,
               log_power_plot,   log_power_plot_legend, 
               "evaluations" = power_table))
  
}

eval_checks_clean_mismatch <- function(data, discr_name, model_name = "", iter = 1000){
  if(discr_name == "sec_moment"){
    ts_and_power <- data %>% filter(discr_fun == discr_name) %>% 
      group_by(size, method, group, discr_fun, rho, rho_squared, trueMod) %>% 
      mutate(iter = n(), ts_or_power = mean(pvals > 0.975 | pvals < 0.025)) %>% 
      summarize(ts_or_power = mean(ts_or_power), 
                true_pw = 1,
                ci =  sqrt(ts_or_power * (1 - ts_or_power)) * 1.96/ sqrt(length(iter))) %>% 
      arrange(method)
    
    df_seg <- data.frame(x = 1:6 + 0.15, y = rep(1,6),
                         xend = 1:6 + 0.45, yend = rep(1,6))
    
    x_lab_name <- expression(nu[o]^2 /nu(theta[0])^2)
    
  }else{
    ts_and_power <- data %>% filter( discr_fun == discr_name) %>%  
      group_by(size, method, group, discr_fun, rho_squared, trueMod) %>% 
      mutate(iter = n(), ts_or_power = mean(pvals > 0.975 | pvals < 0.025)) %>% 
      summarize(ts_or_power = mean(ts_or_power), 
                true_pw = 2*pnorm(qnorm(0.025)/rho),
                ci =  sqrt(ts_or_power * (1 - ts_or_power)) * 1.96/ sqrt(iter),
                MCerror = sqrt(ts_or_power * (1 - ts_or_power))/ sqrt(iter)) %>% 
      arrange(method)
    
    df_seg <- data.frame(x = 1:6 + 0.15, y = sort(as.numeric(unique(ts_and_power$true_pw)), FALSE),
                         xend = 1:6 + 0.45, yend = sort(as.numeric(unique(ts_and_power$true_pw)), FALSE))
    
    x_lab_name <- expression(rho^2)
  }
  
  pd <- position_dodge(0.8) 
  
  ts_and_power <- ts_and_power %>% filter(discr_fun == discr_name, size %in% c(1000,10000,50000,100000))
  
  
  
  base_plot_pw <-  ts_and_power %>% ggplot(aes(factor(group, level = c('1', '2', '3','4','5' ,'6')),
                                               ts_or_power, color = as.factor(method),shape = as.factor(size))) +
    geom_rect(aes(xmin = as.numeric(group) - 0.5, 
                  xmax = as.numeric(group) + 0.5, 
                  ymin = -Inf, ymax = Inf, fill = as.factor(group)), alpha = 0.5, color = NA) +
    # annotate("rect", ymin = -Inf, ymax = Inf,
    #          xmin = c("0.2", "0.9",'21'), xmax =  c("0.2", "0.9", "21"),
    #          alpha = 0.1, fill = c("gray", "gray", "gray")) +
    geom_point(position = pd) + 
    # geom_point(aes(factor(group), true_pw), color = "#FF6666", shape = 4, size = 2, position = pd) + 
    geom_errorbar(aes(ymin=ts_or_power - ci, ymax=ts_or_power + ci), width=.3, position = pd)+
    # labs( y = "power", x = expression(rho^2), color = "method", shape = "size") + 
    labs( y = "power", x = x_lab_name, color = "method", shape = "size") + 
    scale_color_manual(breaks = c("Pop-PC ideal", "PPC", 
                                  "single 0.5-SPC", "divided 0.5-SPC", "True single 0.5-SPC"),
                       values=c("#330000","#0000FF" , 
                                "#009999", "#FF9900", "red")) + 
    scale_shape_manual(values = c(17, 16, 15,18),
                       labels = c("1000", "10000", "50000", "100000")) +
    scale_fill_manual(breaks = c('1', '2', '3','4','5' ,'6'),
                      values=c("gray","white", "gray" , "white", "gray", "white")) +
    scale_x_discrete(labels=c('1' = "0.2", '2' = "0.5", '3' = "0.9", '4' = "5", '5' = "21", '6' = "201"))+
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray") +
    # geom_segment(aes(x = 0.5, y = 1.172645e-05, xend = 1.8, yend = 1.172645e-05, color = 'red', linetype = "dashed"), show.legend = FALSE) +
    # geom_segment(aes(x = 1.5, y = 5.574598e-03, xend = 2.8, yend = 5.574598e-03, color = 'red', linetype = "dashed"), show.legend = FALSE) +
    # geom_segment(aes(x = 2.5, y = 3.883004e-02, xend = 3.8, yend = 3.883004e-02, color = 'red', linetype = "dashed"), show.legend = FALSE)+
    # geom_segment(aes(x = 3.5, y = 3.807460e-01, xend = 4.8, yend = 3.807460e-01, color = 'red', linetype = "dashed"), show.legend = FALSE)+
    # geom_segment(aes(x = 4.5, y = 6.688701e-01, xend = 5.8, yend = 6.688701e-01, color = 'red', linetype = "dashed"), show.legend = FALSE)+
    # geom_segment(aes(x = 5.5, y = 8.900467e-01, xend = 6.8, yend = 8.900467e-01, color = 'red', linetype = "dashed"), show.legend = FALSE)
    geom_segment(data = df_seg, aes(x, y, xend=xend, yend=yend),  color = "red", size = 0.55, linetype = "solid", inherit.aes = FALSE , show.legend = FALSE) 
  
  # aes(factor(type, level = c('under-dispersed','subtle', 
  # 'subtle-to-moderate', 'moderate','moderate-to-major')), 
  
  
  
  power_plot_legend <- base_plot_pw +
    theme(
      #  legend.position = "bottom",
      legend.position = c(0.92, 0.35),  # gaussian mse power
      legend.key.size = unit(0.15, "cm"),
      panel.background = element_rect(fill = FALSE, colour = NA),
      axis.line = element_line(size = 0.5),
      axis.title = element_text(size = 16, color = "black"),
      axis.text = element_text(size = 13, color = "black"),
      axis.ticks = element_blank(),
      legend.text = element_text(size = 10.5, color = "black"),
      legend.background = element_rect(colour = NA, fill = FALSE)
    ) + guides(fill = "none")
  
  
  #+  guides(color=guide_legend(nrow=3,byrow=TRUE))
  power_plot <- base_plot_pw + theme(
    panel.background = element_rect(fill = "white", colour = NA),
    axis.line = element_line(size = 0.5),
    axis.title = element_text(size = 16, color = "black"),
    axis.text = element_text(size = 13, color = "black"),
    axis.ticks = element_blank(),
    legend.position = "none"
  )
  
  # log_power_plot_legend <-  power_plot_legend + 
  #   scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #                 labels = trans_format("log10", math_format(10^.x)))
  # log_power_plot <-  power_plot + 
  #   scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #                 labels = trans_format("log10", math_format(10^.x)))
  # 
  
  
  power_table <- ts_and_power 
  return(list( power_plot,   power_plot_legend, 
               "evaluations" = power_table))
  # return(list( test_size_plot,  power_plot))
  
}



## eval_checks_hier functions: hierachical models
## Input: data.frame, iter(# of iterations for each setting)
##        specify one num_of_groups or num_of_indi used only for plots
##        specify one discr_fun used only for plots
## Output: test size or power under each setting
##         plot test size/ power versus sample size for each setting provided a specific discr_fun and a fixed num_of_groups
## discr_names include "grand_mean", "max_group_means", "mean_group_q75"

size_plot_hier_fixJ_clean <- function(data, discr_name, iter = 200){
  ts_and_power <- data %>% 
    group_by(num_of_groups, num_of_indi, method, discr_fun) %>% 
    mutate(ts_or_power = mean(pvals > 0.975 | pvals < 0.025)) %>% 
    summarize(ts_or_power = mean(ts_or_power), 
              ci = sqrt(ts_or_power * (1 - ts_or_power)) * 1.96/ sqrt(iter)) %>% 
    arrange(method) %>% mutate(rename_J = paste("number of individuals =", paste0("",num_of_indi,""))) 
  
  pd <- position_dodge(0.6) # move them .05 to the left and right
  
  base_plot <- ts_and_power %>%  filter( discr_fun == discr_name) %>% 
    ggplot(aes(x = as.factor(num_of_groups), y = ts_or_power, 
               group = as.factor(method), color = as.factor(method), linetype = as.factor(method))) + 
    geom_point(position = pd) + geom_line(position = pd, size = 0.8)  +
    geom_abline(slope = 0, intercept = 0.05, linetype = "dashed") +
    geom_errorbar(aes(ymin=ts_or_power - ci, ymax=ts_or_power + ci), width=.8, position=pd) + facet_wrap(~as.factor(rename_J)) +
    labs(y = "test size", color = "method", x = "number of groups", linetype = "method") +
    # scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
    scale_color_manual(breaks = c("Pop-PC ideal", "PPC", 
                                  "single 0.3-crossSPC", "single 0.5-crossSPC",
                                  "single 0.3-withinSPC","single 0.5-withinSPC",
                                  "within-divided 0.5-crossSPC","cross-divided 0.5-crossSPC",
                                  "within-divided 0.5-withinSPC","cross-divided 0.5-withinSPC"),
                       values=c("#999999", "#E69F00", 
                                "#009E73", "#009E73",
                                "#0072B2", "#0072B2",
                                "#D55E00", "#D55E00",
                                "#CC79A7","#CC79A7")) + 
    scale_linetype_manual(breaks = c("Pop-PC ideal", "PPC", 
                                     "single 0.3-crossSPC", "single 0.5-crossSPC",
                                     "single 0.3-withinSPC","single 0.5-withinSPC",
                                     "within-divided 0.5-crossSPC","cross-divided 0.5-crossSPC",
                                     "within-divided 0.5-withinSPC","cross-divided 0.5-withinSPC"),
                          values = c("solid", "solid", "dotted", "solid", 
                                     "dotted", "solid", "dotted", "solid",
                                     "dotted", "solid"
                          ))
  
  test_size_plot_legend <- base_plot + theme(
      #legend.position = "none",
      legend.position = "bottom", 
      legend.key.size = unit(0.32, "cm"),
      legend.key = element_blank(),
      legend.key.width= unit(1, 'cm'),
      panel.background = element_rect(fill = "white", colour = NA),
      axis.line = element_line(size = 0.5),
      axis.title = element_text(size = 16, color = "black"),
      axis.text = element_text(size = 13, color = "black"),
      axis.ticks = element_blank(),
      legend.text = element_text(size = 12, color = "black"),
      legend.background = element_rect(color = "grey"),
      legend.title = element_blank(),
      strip.background =element_rect(fill="white"),
      strip.text = element_text(size=12)
    ) + guides(color=guide_legend(nrow=4,byrow=TRUE)) 
  # + scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #                labels = trans_format("log10", math_format(10^.x)))
  
  test_size_plot <- base_plot +
    theme(
      legend.position = "none",
      legend.key.size = unit(0.32, "cm"),
      legend.key = element_blank(),
      legend.key.width= unit(0.6, 'cm'),
      panel.background = element_rect(fill = "white", colour = NA),
      axis.line = element_line(size = 0.5),
      axis.title = element_text(size = 16, color = "black"),
      axis.text = element_text(size = 13, color = "black"),
      axis.ticks = element_blank(),
      legend.text = element_text(size = 12, color = "black"),
      legend.background = element_rect(color = "grey"),
      legend.title = element_blank(),
      strip.background =element_rect(fill="white"),
      strip.text = element_text(size=12)
    )
  # + scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #                    labels = trans_format("log10", math_format(10^.x)))
  return(list(test_size_plot, test_size_plot_legend))
}

eval_checks_hier_fixJ_clean <- function(data, discr_name, iter = 200){
  ts_and_power <- data %>% 
    group_by(num_of_groups, num_of_indi, method, model, discr_fun) %>% 
    mutate(ts_or_power = mean(pvals > 0.975 | pvals < 0.025)) %>% 
    summarize(ts_or_power = mean(ts_or_power), 
              ci =  sqrt(ts_or_power * (1 - ts_or_power)) * 1.96/ sqrt(iter)) %>% 
    arrange(desc(model), method) %>% 
    mutate(rename_J = paste("number of observations per group =", paste0("",num_of_indi,""))) 
  
  pd <- position_dodge(0.6) # move them .05 to the left and right
  
  base_plot <- ts_and_power %>% filter(model == "misspecified") %>% 
    ggplot(aes(x = as.factor(num_of_groups), y = ts_or_power, group = as.factor(method), 
               color = as.factor(method), 
               linetype = as.factor(method))) + geom_point(position = pd) + 
    geom_errorbar(aes(ymin = ts_or_power - ci, ymax = ts_or_power + ci), position = pd, width=1.2) + 
    geom_line(size = 0.6,position = pd) +  facet_wrap(~as.factor(discr_fun)) +
    scale_color_manual(breaks = c("Pop-PC ideal", "PPC", 
                                  "single 0.5-crossSPC",
                                  "single 0.5-withinSPC",
                                  "within-divided 0.5-crossSPC","cross-divided 0.5-crossSPC",
                                  "within-divided 0.5-withinSPC","cross-divided 0.5-withinSPC"),
                       values=c("#999999", "#E69F00", 
                                "#009E73",
                                "#0072B2",
                                "#D55E00", "#D55E00",
                                "#CC79A7","#CC79A7")) + 
    scale_linetype_manual(breaks = c("Pop-PC ideal", "PPC", 
                                     "single 0.5-crossSPC", "single 0.5-withinSPC",
                                     "within-divided 0.5-crossSPC","cross-divided 0.5-crossSPC",
                                     "within-divided 0.5-withinSPC","cross-divided 0.5-withinSPC"),
                          values = c("solid", "solid", "solid", "solid", 
                                     "dotted", "solid", "dotted", "solid"
                          ))
  
  base_legend <- base_plot +  theme(
    #legend.position = "none",
    legend.position = "bottom", 
    legend.key.size = unit(0.32, "cm"),
    legend.key = element_blank(),
    legend.key.width= unit(1, 'cm'),
    panel.background = element_rect(fill = "white", colour = NA),
    axis.line = element_line(size = 0.5),
    axis.title = element_text(size = 16, color = "black"),
    axis.text = element_text(size = 13, color = "black"),
    axis.ticks = element_blank(),
    legend.text = element_text(size = 12, color = "black"),
    legend.background = element_rect(color = "grey"),
    legend.title = element_blank(),
    strip.background =element_rect(fill="white"),
    strip.text = element_text(size=12)
  ) + guides(color=guide_legend(nrow=4,byrow=TRUE)) 
  
  base_non_legend <- base_plot + theme(
    legend.position = "none",
    legend.key.size = unit(0.32, "cm"),
    legend.key = element_blank(),
    legend.key.width= unit(0.6, 'cm'),
    panel.background = element_rect(fill = "white", colour = NA),
    axis.line = element_line(size = 0.5),
    axis.title = element_text(size = 16, color = "black"),
    axis.text = element_text(size = 13, color = "black"),
    axis.ticks = element_blank(),
    legend.text = element_text(size = 12, color = "black"),
    legend.background = element_rect(color = "grey"),
    legend.title = element_blank(),
    strip.background =element_rect(fill="white"),
    strip.text = element_text(size=12)
  )
  
  power_plot_legend <- base_legend +
    labs(y = "power", color = "method", x = "number of groups", linetype = "method") 
  
  
  
  power_plot<- base_non_legend + 
    labs(y = "power", color = "method", x = "number of groups", linetype = "method") 
  # +scale_x_continuous(breaks = seq(0, 500, length.out = 3)) 
  
  
  
  # log_power_plot_legend <-  base_legend +  
  #   labs(y = "power", color = "method", x = "number of groups", linetype = "method") +
  #   scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x, n = 4),
  #                 labels = trans_format("log10", math_format(10^.x))) 
  # log_power_plot <-  base_non_legend + 
  #   labs(y = "power", color = "method", x = "number of groups", linetype = "method") +
  #   scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x, n = 4),
  #                 labels = trans_format("log10", math_format(10^.x)))
  # return(list(  power_plot, power_plot_legend, log_power_plot, log_power_plot_legend,ts_and_power))
  return(list(  power_plot, power_plot_legend, ts_and_power))
  
}

eval_checks_hier_fixI_clean <- function(data, discr_name, iter = 200){
  ts_and_power <- data %>% 
    group_by(num_of_groups, num_of_indi, method, model, discr_fun) %>% 
    mutate(ts_or_power = mean(pvals > 0.975 | pvals < 0.025)) %>% 
    summarize(ts_or_power = mean(ts_or_power), 
              ci = sqrt(ts_or_power * (1 - ts_or_power)) * 1.96/ sqrt(iter)) %>% 
    arrange(desc(model), method) %>% 
    mutate(rename_J = paste("number of groups =", paste0("",num_of_groups,""))) 
  
  pd <- position_dodge(0.6) # move them .05 to the left and right
  
  base_plot <- ts_and_power %>% filter(model == "misspecified") %>% 
    ggplot(aes(x = as.factor(num_of_indi), y = ts_or_power, group = as.factor(method), 
               color = as.factor(method), 
               linetype = as.factor(method))) + geom_point(position = pd) + 
    geom_errorbar(aes(ymin = ts_or_power - ci, ymax = ts_or_power + ci), position = pd, width=1.2) + 
    geom_line(size = 0.6,position = pd) +  facet_wrap(~as.factor(discr_fun)) +
    scale_color_manual(breaks = c("Pop-PC ideal", "PPC", 
                                  "single 0.3-crossSPC", "single 0.5-crossSPC",
                                  "single 0.3-withinSPC","single 0.5-withinSPC",
                                  "within-divided 0.5-crossSPC","cross-divided 0.5-crossSPC",
                                  "within-divided 0.5-withinSPC","cross-divided 0.5-withinSPC"),
                       values=c("#999999", "#E69F00", 
                                "#009E73", "#009E73",
                                "#0072B2", "#0072B2",
                                "#D55E00", "#D55E00",
                                "#CC79A7","#CC79A7")) + 
    scale_linetype_manual(breaks = c("Pop-PC ideal", "PPC", 
                                     "single 0.3-crossSPC", "single 0.5-crossSPC",
                                     "single 0.3-withinSPC","single 0.5-withinSPC",
                                     "within-divided 0.5-crossSPC","cross-divided 0.5-crossSPC",
                                     "within-divided 0.5-withinSPC","cross-divided 0.5-withinSPC"),
                          values = c("solid", "solid", "solid", "solid", 
                                     "solid", "solid", "dotted", "solid",
                                     "dotted", "solid"
                          )) 
  #+ scale_x_continuous(n.breaks = 3)
  
  base_legend <- base_plot +  theme(
    #legend.position = "none",
    legend.position = "bottom", 
    legend.key.size = unit(0.32, "cm"),
    legend.key = element_blank(),
    legend.key.width= unit(1, 'cm'),
    panel.background = element_rect(fill = "white", colour = NA),
    axis.line = element_line(size = 0.5),
    axis.title = element_text(size = 16, color = "black"),
    axis.text = element_text(size = 13, color = "black"),
    axis.ticks = element_blank(),
    legend.text = element_text(size = 12, color = "black"),
    legend.background = element_rect(color = "grey"),
    legend.title = element_blank(),
    strip.background =element_rect(fill="white"),
    strip.text = element_text(size=12)
  ) + guides(color=guide_legend(nrow=4,byrow=TRUE)) 
  
  base_non_legend <- base_plot + theme(
    legend.position = "none",
    legend.key.size = unit(0.32, "cm"),
    legend.key = element_blank(),
    legend.key.width= unit(0.6, 'cm'),
    panel.background = element_rect(fill = "white", colour = NA),
    axis.line = element_line(size = 0.5),
    axis.title = element_text(size = 16, color = "black"),
    axis.text = element_text(size = 13, color = "black"),
    axis.ticks = element_blank(),
    legend.text = element_text(size = 12, color = "black"),
    legend.background = element_rect(color = "grey"),
    legend.title = element_blank(),
    strip.background =element_rect(fill="white"),
    strip.text = element_text(size=12)
  )
  
  power_plot_legend <- base_legend +
    labs(y = "power", color = "method", x = "number of observations per group", linetype = "method") 
  
  
  
  power_plot<- base_non_legend +
    labs(y = "power", color = "method", x = "number of observations per group", linetype = "method") 
  
  
  
  # log_power_plot_legend <-  base_legend +  
  #   labs(y = "power", color = "method", x = "number of observations per group", linetype = "method") +
  #   scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x, n = 4),
  #                 labels = trans_format("log10", math_format(10^.x)))
  # log_power_plot <-  base_non_legend + 
  #   labs(y = "power", color = "method", x = "number of observations per group", linetype = "method") +
  #   scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x, n = 4),
  #                 labels = trans_format("log10", math_format(10^.x)))
  # return(list(  power_plot, power_plot_legend, log_power_plot, log_power_plot_legend,ts_and_power))
  
  return(list(  power_plot, power_plot_legend,ts_and_power))
}


## 2. Q-Q plots under well-specified case
## qq_plot_comparison: simple models
## Input: data.frame, iteration, method_names for each setting
## Output: Q-Q plot for a specific method of all sample sizes 

          
qq_plot_comparison <- function(data, discr_name, model_name = "", iter = 1000){
  data$method_f= factor(data$method, 
                        levels=c('Pop-PC ideal', 'PPC','single 0.9-SPC','single 0.5-SPC',
                                 'divided 0.9-SPC','divided 0.7-SPC','divided 0.5-SPC',
                                 'divided 0.3-SPC','divided 0.1-SPC'))
  bs_plot <- data %>% 
    filter(model == "well-specified", discr_fun == discr_name) %>% 
    group_by(size, method) %>% 
    mutate(unifcdf = 1: iter/(iter + 1), order_pvals = sort(pvals)) %>% 
    ggplot(mapping = aes(unifcdf, order_pvals, color = as.factor(size))) + 
    geom_line(size = 0.5) + 
    geom_abline(slope = 1, intercept = 0, color = "grey") +
    labs(color = "sample size", y = "(KS) p-values") + ylim(0,1) + xlim(0,1) + 
    facet_wrap(~as.factor(method_f)) 
  
  qq_plot_legend <- bs_plot +
    theme(
      panel.background = element_rect(fill = "white", colour = NA),
      strip.background =element_rect(fill="white"),
      axis.line = element_line(size = 0.25),
      panel.border = element_blank(),
      axis.ticks = element_blank(),
      legend.key=element_blank(),
      legend.background = element_rect(color = "grey", fill = "white"))
  
  qq_plot <- bs_plot +
    theme(
      legend.position = "none",
      panel.background = element_rect(fill = "white", colour = NA),
      strip.background =element_rect(fill="white"),
      panel.border = element_blank(),
      axis.ticks = element_blank(),
      legend.key=element_blank(),
      legend.background = element_rect(fill = "white"),
    )
  
  
  return(list(qq_plot, qq_plot_legend))
}


## qq_plot_comparison_hier:  hierachical models
## discr_names include "grand_mean", "max_group_means", "mean_group_q75"

qq_plot_comparison_hier <- function(data, fixed_num_of_indi, discr_name, iter = 200){
   base <- data %>% 
    filter(model == "well-specified", discr_fun == discr_name, num_of_indi == fixed_num_of_indi) %>% 
    group_by(num_of_groups, method) %>% 
    mutate(unifcdf = 1: iter/(iter + 1), order_pvals = sort(pvals)) %>% 
    ggplot(mapping = aes(unifcdf, order_pvals, color = as.factor(num_of_groups))) + 
    geom_line(size = 1) + geom_abline(slope = 1, intercept = 0, color = "grey") +
    labs(color = "number of groups", y = "p-values") + ylim(0,1) + xlim(0,1) +
    facet_wrap(~as.factor(method)) 
   qq_plot <- base + 
     theme(
       legend.position = "none", 
       panel.background = element_rect(fill = "white", colour = NA),
       strip.background =element_rect(fill="white"),
       axis.line = element_line(size = 0.5),
       axis.title = element_text(size = 14, color = "black"),
       axis.text = element_text(size = 12, color = "black"),
       axis.ticks = element_blank(),
       legend.text = element_text(size = 12, color = "black"),
       legend.key=element_blank(),
       legend.background = element_rect(color = "grey", fill = "white")
     )  
  qq_plot_legend <- base + 
    theme(
      # legend.position = "bottom", 
      panel.background = element_rect(fill = "white", colour = NA),
      strip.background =element_rect(fill="white"),
      axis.line = element_line(size = 0.5),
      axis.title = element_text(size = 14, color = "black"),
      axis.text = element_text(size = 12, color = "black"),
      axis.ticks = element_blank(),
      legend.text = element_text(size = 12, color = "black"),
      legend.key=element_blank(),
      legend.background = element_rect(color = "grey", fill = "white")
    )  
 
  return(list(qq_plot, qq_plot_legend))
  
}

