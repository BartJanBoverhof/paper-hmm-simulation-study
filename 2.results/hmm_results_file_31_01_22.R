###~~~~~~~~~~~~~~~~~~~~~~~~~~###
###~~~~~ HMM Paper Code ~~~~~###
###~~~~~~~~~~~~~~~~~~~~~~~~~~###

## Please press ALT + O (Cmd + Option + O on Mac) to collapse the script. 

###~~~~~~~~~~~~~~~~~~~~~~~~###
### Load packages and data ###
###~~~~~~~~~~~~~~~~~~~~~~~~###

### Load packages.
## Check if package "data.table" can be required, if not install. 
if(!require(data.table)) install.packages("data.table")
## Check if package "reshape" can be required, if not install. 
if(!require(reshape)) install.packages("reshape")
## Check if package "ggh4x" can be required, if not install. 
if(!require(ggh4x)) install.packages("ggh4x")
## Check if package "tidyverse" can be required, if not install. 
if(!require(tidyverse)) install.packages("tidyverse")
## Check if package "gridExtra" can be required, if not install. 
if(!require(gridExtra)) install.packages("plyr")
## Check if package "Hmisc" can be required, if not install. 
if(!require(Hmisc)) install.packages("Hmisc")
## Check if package "cat" can be required, if not install. 
if(!require(cat)) install.packages("cat")
## Require packages.
packages <- c("data.table", "reshape", "ggh4x", "tidyverse", "gridExtra", "Hmisc", "cat")
lapply(packages, library, character.only = TRUE)

### Set working directory. 
setwd("C:\\Academia\\HMM paper\\Data\\Data")

### Function for assigning data sets to list.
rda2list <- function(file) {
  e <- new.env()
  load(file, envir = e)
  as.list(e)
}

### Loading and subsequently assigning separate data sets to a list. 
files <- list.files(pattern = ".rda$")
results <- Map(rda2list, file.path(files))
names(results) <- tools::file_path_sans_ext(files)

###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
### Functions & Constant object definitions ###
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###

### Supress "friendly" warnings. 
options(dplyr.summarise.inform = F)

### Define a function for calculating absolute bias, relative bias, standard deviation, and coverage of emission probabilities. 
calcR_theta <- function(true_probs, nobs, cond, condition, cond_string, varobs = F){
  
  ### Absolute bias. 
  abs_bias <- data.frame()
  for (i in 1: length(lengths)){
    abs_cond1 <- cbind(melt(sweep(eval(parse(text = paste0("results$sim_HMM_", cond[1] ,"_obs_", nobs[1], "_t_", lengths[i], "$out_sim$emiss_mean"))), 
                        2, unlist(true_probs[1]))), length = lengths[i], condition = condition[1])
    abs_cond2 <- cbind(melt(sweep(eval(parse(text = paste0("results$sim_HMM_", cond[2] ,"_obs_", nobs[2], "_t_", lengths[i], "$out_sim$emiss_mean"))), 
                        2, unlist(true_probs[2]))), length = lengths[i], condition = condition[2])
    abs_cond3 <- cbind(melt(sweep(eval(parse(text = paste0("results$sim_HMM_", cond[3] ,"_obs_", nobs[3] ,"_t_", lengths[i], "$out_sim$emiss_mean"))), 
                        2, unlist(true_probs[3]))), length = lengths[i], condition = condition[3])
    abs_bias <- rbind(abs_bias, abs_cond1, abs_cond2, abs_cond3)
    }
  ord <- 1 : (length(levels(abs_bias$X2)) / 3)
  abs_bias <- abs_bias %>%
    dplyr::rename(cats = X2, abs_bias = value) %>%
    mutate(length = factor(length, levels = c("250", "500", "1000", "2000", "4000", "8000")),
           cats = factor(cats, levels = c(paste0("M_S1_cat", ord), paste0("M_S2_cat", ord), paste0("M_S3_cat", ord)))) %>%
    group_by(cats, length, condition) %>%
    dplyr::summarise(mean_abs_bias = mean(abs_bias))
  abs_bias_complete_plot <- abs_bias %>% 
    ggplot(aes(x = length, y = mean_abs_bias , color = condition, group = condition)) +
    facet_wrap(facets = vars(cats), nrow = 3, ncol = nobs[3], labeller = label_parsed) +
    geom_point() + 
    geom_line() +
    xlab("Sequence length") + 
    ylab("Absolute bias") +
    labs(fill = cond_string) +
    ggtitle(paste0("Absolute Bias Emission Probabilities ", cond_string , " Condition - Complete")) +
    labs(color = cond_string) +
    geom_hline(yintercept = 0) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 
  
  ### Relative bias. 
  rel_bias <- data.frame()
  for (i in 1: length(lengths)){
    rel_cond1 <- cbind(melt(sweep(sweep(eval(parse(text = paste0("results$sim_HMM_", cond[1] ,"_obs_", nobs[1], "_t_", lengths[i], "$out_sim$emiss_mean"))), 
                              2, unlist(true_probs[1])), 2, unlist(true_probs[1]), FUN = '/')), length = lengths[i], condition = condition[1])
    rel_cond2 <- cbind(melt(sweep(sweep(eval(parse(text = paste0("results$sim_HMM_", cond[2] ,"_obs_", nobs[2], "_t_", lengths[i], "$out_sim$emiss_mean"))), 
                              2, unlist(true_probs[2])), 2, unlist(true_probs[2]), FUN = '/')), length = lengths[i], condition = condition[2])
    rel_cond3 <- cbind(melt(sweep(sweep(eval(parse(text = paste0("results$sim_HMM_", cond[3] ,"_obs_", nobs[3], "_t_", lengths[i], "$out_sim$emiss_mean"))), 
                             2, unlist(true_probs[3])), 2, unlist(true_probs[3]), FUN = '/')), length = lengths[i], condition = condition[3])
    rel_bias <- rbind(rel_bias, rel_cond1, rel_cond2, rel_cond3)
    }
  ord <- 1 : (length(levels(rel_bias$X2)) / 3)
  rel_bias <- rel_bias %>%
    dplyr::rename(cats = X2, rel_bias = value) %>%
    mutate(length = factor(length, levels = c("250", "500", "1000", "2000", "4000", "8000")),
           cats = factor(cats, levels = c(paste0("M_S1_cat", ord), paste0("M_S2_cat", ord), paste0("M_S3_cat", ord)))) %>%
    group_by(cats, length, condition) %>%
    dplyr::summarise(mean_rel_bias = mean(rel_bias)) 
  rel_bias_complete_plot <- rel_bias %>%
    ggplot(aes(x = length, y = mean_rel_bias , color = condition, group = condition)) +
    facet_wrap(facets = vars(cats), nrow = 3, ncol = nobs[3], labeller = label_parsed) +
    geom_point() + 
    geom_line() +
    geom_hline(yintercept = 0.1, linetype = "dashed", color = "black") +
    geom_hline(yintercept = -0.1, linetype = "dashed", color = "black") +
    xlab("Sequence length") + 
    ylab("Relative bias") +
    ggtitle(paste0("Relative Bias Emission Probabilities ", cond_string , " Condition - Complete")) +
    labs(color = cond_string) +
    geom_hline(yintercept = 0) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 
  
  ### Standard deviation. 
  sdev <- data.frame()
  for (i in 1: length(lengths)){
    sdev_cond1 <- cbind(melt(eval(parse(text = paste0("results$sim_HMM_", cond[1], "_obs_", nobs[1], "_t_", lengths[i], "$out_sim$emiss_sd")))), 
                        length = lengths[i], condition = condition[1])
    sdev_cond2 <- cbind(melt(eval(parse(text = paste0("results$sim_HMM_", cond[2], "_obs_", nobs[2], "_t_", lengths[i], "$out_sim$emiss_sd")))), 
                        length = lengths[i], condition = condition[2])
    sdev_cond3 <- cbind(melt(eval(parse(text = paste0("results$sim_HMM_", cond[3], "_obs_", nobs[3], "_t_", lengths[i], "$out_sim$emiss_sd")))), 
                        length = lengths[i], condition = condition[3])
    sdev <- rbind(sdev, sdev_cond1, sdev_cond2, sdev_cond3)
    }
  ord <- 1 : (length(levels(sdev$X2)) / 3)
  sdev <- sdev %>%
    dplyr::rename(cats = X2, sdev = value) %>%
    mutate(length = factor(length, levels = c("250", "500", "1000", "2000", "4000", "8000")),
           cats = factor(cats, levels = c(paste0("sd_S1_cat", ord), paste0("sd_S2_cat", ord), paste0("sd_S3_cat", ord)))) %>%
    group_by(cats, length, condition) %>%
    dplyr::summarise(mean_sdev = mean(sdev)) 
  sdev_complete_plot <- sdev %>%
    ggplot(aes(x = length, y = mean_sdev , color = condition, group = condition)) +
    facet_wrap(facets = vars(cats), nrow = 3, ncol = nobs[3], labeller = label_parsed) +
    geom_point() + 
    geom_line() +
    xlab("Sequence length") + 
    ylab("Standard deviation") +
    ggtitle(paste0("Standard Deviation Emission Probabilities ", cond_string , " Condition - Complete")) +
    labs(color = cond_string) +
    geom_hline(yintercept = 0) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 
  
  if (varobs == T){
    
    ### Coverage. 
    cov <- data.frame()
    for (i in 1: length(lengths)){
      cov_cond1 <- cbind(melt(sweep(eval(parse(text = paste0("results$sim_HMM_", cond[1], "_obs_", nobs[1], "_t_", lengths[i], "$out_sim$emiss_mean"))), 
                                    2, 0)), length = lengths[i], condition = condition[1], 
                         ci_low = as.vector(eval(parse(text = paste0("results$sim_HMM_", cond[1], "_obs_", nobs[1], "_t_", lengths[i], "$out_sim$emiss_low")))), 
                         ci_high = as.vector(eval(parse(text = paste0("results$sim_HMM_", cond[1], "_obs_", nobs[1], "_t_", lengths[i], "$out_sim$emiss_up")))))
      cov_cond2 <- cbind(melt(sweep(eval(parse(text = paste0("results$sim_HMM_", cond[2], "_obs_", nobs[2], "_t_", lengths[i], "$out_sim$emiss_mean"))), 
                                    2, 0)), length = lengths[i], condition = condition[2], 
                         ci_low = as.vector(eval(parse(text = paste0("results$sim_HMM_", cond[2], "_obs_", nobs[2], "_t_", lengths[i], "$out_sim$emiss_low")))), 
                         ci_high = as.vector(eval(parse(text = paste0("results$sim_HMM_", cond[2], "_obs_", nobs[2], "_t_", lengths[i], "$out_sim$emiss_up")))))
      cov_cond3 <- cbind(melt(sweep(eval(parse(text = paste0("results$sim_HMM_", cond[3], "_obs_", nobs[3], "_t_", lengths[i], "$out_sim$emiss_mean"))), 
                                    2, 0)), length = lengths[i], condition = condition[3], 
                         ci_low = as.vector(eval(parse(text = paste0("results$sim_HMM_", cond[3], "_obs_", nobs[3], "_t_", lengths[i], "$out_sim$emiss_low")))), 
                         ci_high = as.vector(eval(parse(text = paste0("results$sim_HMM_", cond[3], "_obs_", nobs[3], "_t_", lengths[i], "$out_sim$emiss_up")))))
      cov <- rbind(cov, cov_cond1, cov_cond2, cov_cond3)
      }
    ord <- 1 : (length(levels(cov$X2)) / 3)
    cov <- cov %>%
      dplyr::rename(cats = X2, cov = value) %>%
      mutate(length = factor(length, levels = c("250", "500", "1000", "2000", "4000", "8000")),
             cats = factor(cats, levels = c(paste0("M_S1_cat", ord), paste0("M_S2_cat", ord), paste0("M_S3_cat", ord)))) 
    foo1 <- cbind(expand.grid(c(paste0("M_S1_cat", 1:3), paste0("M_S2_cat", 1:3), paste0("M_S3_cat", 1:3)), levels(cov$condition)[1]), unlist(true_probs[1]))
    foo2 <- cbind(expand.grid(c(paste0("M_S1_cat", 1:5), paste0("M_S2_cat", 1:5), paste0("M_S3_cat", 1:5)), levels(cov$condition)[2]), unlist(true_probs[2]))
    foo3 <- cbind(expand.grid(c(paste0("M_S1_cat", 1:7), paste0("M_S2_cat", 1:7), paste0("M_S3_cat", 1:7)), levels(cov$condition)[3]), unlist(true_probs[3]))
    colnames(foo1) <- c("cats", "condition", "true")
    colnames(foo2) <- c("cats", "condition", "true")
    colnames(foo3) <- c("cats", "condition", "true")
    temp <- rbind(foo1, foo2, foo3)
    cov <- full_join(cov, temp, by = c("cats", "condition"))
    cov <- cov %>%
      select(-cov) %>%
      mutate(cov = ifelse((true > ci_low & true < ci_high), 1, 0)) %>%
      group_by(cats, length, condition) %>%
      dplyr::summarise(mean_cov = mean(cov))
    cov_complete_plot <- cov %>%
      ggplot(aes(x = length, y = mean_cov , color = condition, group = condition)) +
      facet_wrap(facets = vars(cats), nrow = 3, ncol = nobs[3], labeller = label_parsed) +
      geom_point() + 
      geom_line() +
      geom_hline(yintercept = 0.95, linetype = "dashed", color = "black") +
      xlab("Sequence length") + 
      ylab("Coverage") +
      ggtitle(paste0("Coverage Emission Probabilities ", cond_string , " Condition - Complete")) +
      labs(color = cond_string) +
      geom_hline(yintercept = 0) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 
    
    } else {
      
      ### Coverage. 
      cov <- data.frame()
      for (i in 1: length(lengths)){
        cov_cond1 <- cbind(melt(sweep(eval(parse(text = paste0("results$sim_HMM_", cond[1], "_obs_", nobs[1], "_t_", lengths[i], "$out_sim$emiss_mean"))), 
                                      2, 0)), length = lengths[i], condition = condition[1], 
                           ci_low = as.vector(eval(parse(text = paste0("results$sim_HMM_", cond[1], "_obs_", nobs[1], "_t_", lengths[i], "$out_sim$emiss_low")))), 
                           ci_high = as.vector(eval(parse(text = paste0("results$sim_HMM_", cond[1], "_obs_", nobs[1], "_t_", lengths[i], "$out_sim$emiss_up")))))
        cov_cond2 <- cbind(melt(sweep(eval(parse(text = paste0("results$sim_HMM_", cond[2], "_obs_", nobs[2], "_t_", lengths[i], "$out_sim$emiss_mean"))), 
                                      2, 0)), length = lengths[i], condition = condition[2], 
                           ci_low = as.vector(eval(parse(text = paste0("results$sim_HMM_", cond[2], "_obs_", nobs[2], "_t_", lengths[i], "$out_sim$emiss_low")))), 
                           ci_high = as.vector(eval(parse(text = paste0("results$sim_HMM_", cond[2], "_obs_", nobs[2], "_t_", lengths[i], "$out_sim$emiss_up")))))
        cov_cond3 <- cbind(melt(sweep(eval(parse(text = paste0("results$sim_HMM_", cond[3], "_obs_", nobs[3], "_t_", lengths[i], "$out_sim$emiss_mean"))), 
                                      2, 0)), length = lengths[i], condition = condition[3], 
                           ci_low = as.vector(eval(parse(text = paste0("results$sim_HMM_", cond[3], "_obs_", nobs[3], "_t_", lengths[i], "$out_sim$emiss_low")))), 
                           ci_high = as.vector(eval(parse(text = paste0("results$sim_HMM_", cond[3], "_obs_", nobs[3], "_t_", lengths[i], "$out_sim$emiss_up")))))
        cov <- rbind(cov, cov_cond1, cov_cond2, cov_cond3)
        }
      ord <- 1 : (length(levels(cov$X2)) / 3)
      cov <- cov %>%
        dplyr::rename(cats = X2, cov = value) %>%
        mutate(length = factor(length, levels = c("250", "500", "1000", "2000", "4000", "8000")),
               cats = factor(cats, levels = c(paste0("M_S1_cat", ord), paste0("M_S2_cat", ord), paste0("M_S3_cat", ord)))) 
      temp <- cbind(expand.grid(levels(cov$cats), levels(cov$condition)), unlist(true_probs))
      colnames(temp) <- c("cats", "condition", "true")
      cov <- full_join(cov, temp, by = c("cats", "condition"))
      cov <- cov %>%
        select(-cov) %>%
        mutate(cov = ifelse((true > ci_low & true < ci_high), 1, 0)) %>%
        group_by(cats, length, condition) %>%
        dplyr::summarise(mean_cov = mean(cov))
      cov_complete_plot <- cov %>%
        ggplot(aes(x = length, y = mean_cov , color = condition, group = condition)) +
        facet_wrap(facets = vars(cats), nrow = 3, ncol = nobs[3], labeller = label_parsed) +
        geom_point() + 
        geom_line() +
        geom_hline(yintercept = 0.95, linetype = "dashed", color = "black") +
        xlab("Sequence length") + 
        ylab("Coverage") +
        ggtitle(paste0("Coverage Emission Probabilities ", cond_string , " Condition - Complete")) +
        labs(color = cond_string) +
        geom_hline(yintercept = 0) +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 
      
        }
  
  return(list(abs_bias = list(abs_bias = abs_bias, abs_bias_complete_plot = abs_bias_complete_plot), 
              rel_bias = list(rel_bias = rel_bias, rel_bias_complete_plot = rel_bias_complete_plot), 
              sdev = list(sdev = sdev, sdev_complete_plot = sdev_complete_plot), 
              cov = list(cov = cov, cov_complete_plot = cov_complete_plot)))
  
  }

### Define a function for calculating absolute bias, relative bias, standard deviation, and coverage of transition probabilities. 
calcR_gamma <- function(true_probs, nobs, cond, condition, cond_string, varobs = F){
  
  ### Absolute bias. 
  abs_bias <- data.frame()
  for (i in 1: length(lengths)){
    abs_cond1 <- cbind(melt(sweep(eval(parse(text = paste0("results$sim_HMM_", cond[1], "_obs_", nobs[1], "_t_", lengths[i], "$out_sim$gamma_mean"))), 
                                  2, true_probs)), length = lengths[i], condition = condition[1])
    abs_cond2 <- cbind(melt(sweep(eval(parse(text = paste0("results$sim_HMM_", cond[2], "_obs_", nobs[2], "_t_", lengths[i], "$out_sim$gamma_mean"))), 
                                  2, true_probs)), length = lengths[i], condition = condition[2])
    abs_cond3 <- cbind(melt(sweep(eval(parse(text = paste0("results$sim_HMM_", cond[3], "_obs_", nobs[3] ,"_t_", lengths[i], "$out_sim$gamma_mean"))), 
                                  2, true_probs)), length = lengths[i], condition = condition[3])
    abs_bias <- rbind(abs_bias, abs_cond1, abs_cond2, abs_cond3)
    }
  abs_bias <- abs_bias %>%
    dplyr::rename(cats = X2, abs_bias = value) %>%
    mutate(factor(length, levels = c("250", "500", "1000", "2000", "4000", "8000"))) %>%
    group_by(cats, length, condition) %>%
    dplyr::summarise(mean_abs_bias = mean(abs_bias))
  abs_bias_complete_plot <- abs_bias %>%
    ggplot(aes(x = length, y = mean_abs_bias , color = condition, group = condition)) +
    facet_wrap(facets = vars(cats), nrow = 3, ncol = 3, labeller = label_parsed) +
    geom_point() + 
    geom_line() +
    xlab("Sequence length") + 
    ylab("Absolute bias") +
    ggtitle(paste0("Absolute Bias Transition Probabilities ", cond_string , " Condition - Complete")) +
    labs(color = cond_string) +
    geom_hline(yintercept = 0) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 
  
  ### Relative bias. 
  rel_bias <- data.frame()
  for (i in 1: length(lengths)){
    rel_cond1 <- cbind(melt(sweep(sweep(eval(parse(text = paste0("results$sim_HMM_", cond[1], "_obs_", nobs[1], "_t_", lengths[i], "$out_sim$gamma_mean"))), 
                                        2, true_probs), 2, true_probs, FUN = '/')), length = lengths[i], condition = condition[1])
    rel_cond2 <- cbind(melt(sweep(sweep(eval(parse(text = paste0("results$sim_HMM_", cond[2], "_obs_", nobs[2], "_t_", lengths[i], "$out_sim$gamma_mean"))), 
                                        2, true_probs), 2, true_probs, FUN = '/')), length = lengths[i], condition = condition[2])
    rel_cond3 <- cbind(melt(sweep(sweep(eval(parse(text = paste0("results$sim_HMM_", cond[3], "_obs_", nobs[3], "_t_", lengths[i], "$out_sim$gamma_mean"))), 
                                        2, true_probs), 2, true_probs, FUN = '/')), length = lengths[i], condition = condition[3])
    rel_bias <- rbind(rel_bias, rel_cond1, rel_cond2, rel_cond3)
    }
  rel_bias <- rel_bias %>%
    dplyr::rename(cats = X2, rel_bias = value) %>%
    mutate(factor(length, levels = c("250", "500", "1000", "2000", "4000", "8000"))) %>%
    group_by(cats, length, condition) %>%
    dplyr::summarise(mean_rel_bias = mean(rel_bias))
  rel_bias_complete_plot <- rel_bias %>% 
    ggplot(aes(x = length, y = mean_rel_bias , color = condition, group = condition)) +
    facet_wrap(facets = vars(cats), nrow = 3, ncol = 3, labeller = label_parsed) +
    geom_point() + 
    geom_line() +
    geom_hline(yintercept = 0.1, linetype = "dashed", color = "black") +
    geom_hline(yintercept = -0.1, linetype = "dashed", color = "black") +
    xlab("Sequence length") + 
    ylab("Relative bias") +
    ggtitle(paste0("Relative Bias Transition Probabilities ", cond_string , " Condition - Complete")) +
    labs(color = cond_string) +
    geom_hline(yintercept = 0) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 
  
  ### Standard deviation. 
  sdev <- data.frame()
  for (i in 1: length(lengths)){
    sdev_cond1 <- cbind(melt(eval(parse(text = paste0("results$sim_HMM_", cond[1],"_obs_", nobs[1], "_t_", lengths[i], "$out_sim$gamma_sd")))), 
                        length = lengths[i], condition = condition[1])
    sdev_cond2 <- cbind(melt(eval(parse(text = paste0("results$sim_HMM_", cond[2],"_obs_", nobs[2], "_t_", lengths[i], "$out_sim$gamma_sd")))), 
                        length = lengths[i], condition = condition[2])
    sdev_cond3 <- cbind(melt(eval(parse(text = paste0("results$sim_HMM_", cond[3],"_obs_", nobs[3], "_t_", lengths[i], "$out_sim$gamma_sd")))), 
                        length = lengths[i], condition = condition[3])
    sdev <- rbind(sdev, sdev_cond1, sdev_cond2, sdev_cond3)
    }
  sdev <- sdev %>%
    dplyr::rename(cats = X2, sdev = value) %>%
    mutate(factor(length, levels = c("250", "500", "1000", "2000", "4000", "8000"))) %>%
    group_by(cats, length, condition) %>%
    dplyr::summarise(mean_sdev = mean(sdev))
  sdev_complete_plot <- sdev %>% 
    ggplot(aes(x = length, y = mean_sdev , color = condition, group = condition)) +
    facet_wrap(facets = vars(cats), nrow = 3, ncol = 3, labeller = label_parsed) +
    geom_point() + 
    geom_line() +
    xlab("Sequence length") + 
    ylab("Standard deviation") +
    ggtitle(paste0("Standard Deviation Transition Probabilities ", cond_string , " Condition - Complete")) +
    labs(color = cond_string) +
    geom_hline(yintercept = 0) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 
  
  ### Coverage. 
  cov <- data.frame()
  for (i in 1: length(lengths)){
    cov_cond1 <- cbind(melt(sweep(eval(parse(text = paste0("results$sim_HMM_", cond[1], "_obs_", nobs[1], "_t_", lengths[i], "$out_sim$gamma_mean"))), 
                                  2, 0)), length = lengths[i], condition = condition[1], 
                       ci_low = as.vector(eval(parse(text = paste0("results$sim_HMM_", cond[1], "_obs_", nobs[1], "_t_", lengths[i], "$out_sim$gamma_low")))), 
                       ci_high = as.vector(eval(parse(text = paste0("results$sim_HMM_", cond[1], "_obs_", nobs[1], "_t_", lengths[i], "$out_sim$gamma_up")))))
    cov_cond2 <- cbind(melt(sweep(eval(parse(text = paste0("results$sim_HMM_", cond[2], "_obs_", nobs[2], "_t_", lengths[i], "$out_sim$gamma_mean"))), 
                                  2, 0)), length = lengths[i], condition = condition[2], 
                       ci_low = as.vector(eval(parse(text = paste0("results$sim_HMM_", cond[2], "_obs_", nobs[2], "_t_", lengths[i], "$out_sim$gamma_low")))), 
                       ci_high = as.vector(eval(parse(text = paste0("results$sim_HMM_", cond[2], "_obs_", nobs[2], "_t_", lengths[i], "$out_sim$gamma_up")))))
    cov_cond3 <- cbind(melt(sweep(eval(parse(text = paste0("results$sim_HMM_", cond[3], "_obs_", nobs[3], "_t_", lengths[i], "$out_sim$gamma_mean"))), 
                                  2, 0)), length = lengths[i], condition = condition[3], 
                       ci_low = as.vector(eval(parse(text = paste0("results$sim_HMM_", cond[3], "_obs_", nobs[3], "_t_", lengths[i], "$out_sim$gamma_low")))), 
                       ci_high = as.vector(eval(parse(text = paste0("results$sim_HMM_", cond[3], "_obs_", nobs[3], "_t_", lengths[i], "$out_sim$gamma_up")))))
    cov <- rbind(cov, cov_cond1, cov_cond2, cov_cond3)
    }
  cov <- cov %>%
    dplyr::rename(cats = X2, cov = value) %>%
    mutate(length = factor(length, levels = c("250", "500", "1000", "2000", "4000", "8000"))) 
  temp <- cbind(expand.grid(levels(cov$cats), levels(cov$condition)), unlist(as.vector(gamma_true)))
  colnames(temp) <- c("cats", "condition", "true")
  cov <- full_join(cov, temp, by = c("cats", "condition"))
  cov <- cov %>%
    select(-cov) %>%
    mutate(cov = ifelse((true > ci_low & true < ci_high), 1, 0)) %>%
    group_by(cats, length, condition) %>%
    dplyr::summarise(mean_cov = mean(cov))
  cov_complete_plot <- cov %>%
    ggplot(aes(x = length, y = mean_cov , color = condition, group = condition)) +
    facet_wrap(facets = vars(cats), nrow = 3, ncol = 3, labeller = label_parsed) +
    geom_point() + 
    geom_line() +
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "black") +
    xlab("Sequence length") + 
    ylab("Coverage") +
    ggtitle(paste0("Coverage Transition Probabilities ", cond_string , " Condition - Complete")) +
    labs(color = cond_string) +
    geom_hline(yintercept = 0) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 
  
  return(list(abs_bias = list(abs_bias = abs_bias, abs_bias_complete_plot = abs_bias_complete_plot), 
              rel_bias = list(rel_bias = rel_bias, rel_bias_complete_plot = rel_bias_complete_plot), 
              sdev = list(sdev = sdev, sdev_complete_plot = sdev_complete_plot), 
              cov = list(cov = cov, cov_complete_plot = cov_complete_plot)))
  
  }

### Define various lengths of the simulation chains, this object remains constant throughout the script.  
lengths <- c("250", "500", "1000", "2000", "4000", "8000")

### Define the population level transition probabilities, these are held constant over all the conditions. 
gamma_true <- array(c(c(0.80, 0.10, 0.10), c(0.10, 0.80, 0.10), c(0.10, 0.10, 0.80)), dim = c(1, 9, 1))

###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
### Absolute bias, relative bias, and standard error of the emission and transition probabilities clarity condition ###
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###

### Calculate absolute bias, relative bias, standard deviation, and coverage of emission probabilities of the clarity condition, 
### with associated plots. 

## Specify number of observations. 
nobs <- c(5, 5, 5)
## Specify condition.
cond <- c("theta_cl", "theta_modcl", "theta_uncl")
condition <- c("Clear", "Moderately clear", "Unclear")
cond_string <- c("Clarity")
## Specify true probabilities.
true_probs <- list(c(0.92, 0.02, 0.02, 0.02, 0.02, 0.02, 0.47, 0.02, 0.47, 0.02, 0.02, 0.02, 0.47, 0.02, 0.47), 
                   c(0.68, 0.08, 0.08, 0.08, 0.08, 0.08, 0.38, 0.08, 0.38, 0.08, 0.08, 0.08, 0.38, 0.08, 0.38), 
                   c(0.44, 0.14, 0.14, 0.14, 0.14, 0.14, 0.29, 0.14, 0.29, 0.14, 0.14, 0.14, 0.29, 0.14, 0.29))
## Calculate the various measures and the associated plots. 
clarity_theta <- calcR_theta(true_probs = true_probs, nobs = nobs, cond = cond, condition = condition, cond_string = cond_string, varobs = F)

### Reduced and combined plots of the absolute bias, relative bias, standard deviation, and coverage of emission probabilities of the 
### clarity condition. 
{
## Absolute bias. 
# Reduced plot. 
clarity_theta$abs_bias$abs_bias$cats <- car::recode(clarity_theta$abs_bias$abs_bias$cats, 
                                                    "c('M_S1_cat1') = 'One-indicator state'; 
                                                     c('M_S2_cat2', 'M_S2_cat4', 'M_S3_cat3', 'M_S3_cat5') = 'Two-indicator state';
                                                     else = 'Noise'")
clarity_theta$abs_bias$abs_bias$cats <- factor(clarity_theta$abs_bias$abs_bias$cats, 
                                               levels = c("One-indicator state", "Two-indicator state", "Noise"))
clarity_theta_abs_bias_reduced_plot <- clarity_theta$abs_bias$abs_bias %>%
  group_by(cats, length, condition) %>%
  dplyr::summarize(mean = mean(mean_abs_bias), sd = sd(mean_abs_bias)) %>% 
  ggplot(aes(x = length, y = mean, color = condition, group = condition)) +
  facet_wrap(facets = vars(cats)) +
  geom_point() + 
  geom_line() +  
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.1) +
  xlab("Sequence length") + 
  ylab("Absolute bias") +
  ggtitle("Absolute Bias Emission Probabilities Clarity Condition - Reduced") + 
  labs(color = "Clarity") +
  geom_hline(yintercept = 0) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 
# Combined plot. 
abs_bias_clarity_theta_combined_plot <- grid.arrange(clarity_theta$abs_bias$abs_bias_complete_plot, 
                                                     clarity_theta_abs_bias_reduced_plot, ncol = 2)

## Relative bias. 
# Reduced plot. 
clarity_theta$rel_bias$rel_bias$cats <- car::recode(clarity_theta$rel_bias$rel_bias$cats, 
                                                    "c('M_S1_cat1') = 'One-indicator state'; 
                                                     c('M_S2_cat2', 'M_S2_cat4', 'M_S3_cat3', 'M_S3_cat5') = 'Two-indicator state';
                                                     else = 'Noise'")
clarity_theta$rel_bias$rel_bias$cats <- factor(clarity_theta$rel_bias$rel_bias$cats, 
                                               levels = c("One-indicator state", "Two-indicator state", "Noise"))
clarity_theta_rel_bias_reduced_plot <- clarity_theta$rel_bias$rel_bias %>% 
  group_by(cats, length, condition) %>%
  dplyr::summarize(mean = mean(mean_rel_bias), sd = sd(mean_rel_bias)) %>% 
  ggplot(aes(x = length, y = mean, color = condition, group = condition)) +
  facet_wrap(facets = vars(cats)) +
  geom_point() + 
  geom_line() +  
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2) +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "black") +
  geom_hline(yintercept = -0.1, linetype = "dashed", color = "black") +
  xlab("Sequence length") + 
  ylab("Relative bias") +
  ggtitle("Relative Bias Emission Probabilities Clarity Condition - Reduced") + 
  labs(color = "Clarity") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 
# Combined plot. 
rel_bias_clarity_theta_combined_plot <- grid.arrange(clarity_theta$rel_bias$rel_bias_complete_plot, 
                                                     clarity_theta_rel_bias_reduced_plot, ncol = 2)

## Standard deviation. 
# Reduced plot. 
clarity_theta$sdev$sdev$cats <- car::recode(clarity_theta$sdev$sdev$cats, 
                                            "c('sd_S1_cat1') = 'One-indicator state'; 
                                             c('sd_S2_cat2', 'sd_S2_cat4', 'sd_S3_cat3', 'sd_S3_cat5') = 'Two-indicator state';
                                             else = 'Noise'")
clarity_theta$sdev$sdev$cats <- factor(clarity_theta$sdev$sdev$cats, 
                                               levels = c("One-indicator state", "Two-indicator state", "Noise"))
clarity_theta_sdev_reduced_plot <- clarity_theta$sdev$sdev %>% 
  group_by(cats, length, condition) %>%
  dplyr::summarize(mean = mean(mean_sdev), sd = sd(mean_sdev)) %>% 
  ggplot(aes(x = length, y = mean, color = condition, group = condition)) +
  facet_wrap(facets = vars(cats)) +
  geom_point() + 
  geom_line() +  
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.1) +
  xlab("Sequence length") + 
  ylab("Standard deviation") +
  ggtitle("Standard Deviation Emission Probabilities Clarity Condition - Reduced") + 
  labs(color = "Clarity") +
  geom_hline(yintercept = 0) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 
# Combined plot. 
sdev_clarity_theta_combined_plot <- grid.arrange(clarity_theta$sdev$sdev_complete_plot, 
                                                 clarity_theta_sdev_reduced_plot, ncol = 2)

## Coverage 
# Reduced plot. 
clarity_theta$cov$cov$cats <- car::recode(clarity_theta$cov$cov$cats, 
                                          "c('M_S1_cat1') = 'One-indicator state'; 
                                           c('M_S2_cat2', 'M_S2_cat4', 'M_S3_cat3', 'M_S3_cat5') = 'Two-indicator state';
                                           else = 'Noise'")
clarity_theta$cov$cov$cats <- factor(clarity_theta$cov$cov$cats, 
                                     levels = c("One-indicator state", "Two-indicator state", "Noise"))
clarity_theta_cov_reduced_plot <- clarity_theta$cov$cov %>% 
  group_by(cats, length, condition) %>%
  dplyr::summarize(mean = mean(mean_cov), sd = sd(mean_cov)) %>% 
  ggplot(aes(x = length, y = mean, color = condition, group = condition)) +
  facet_wrap(facets = vars(cats)) +
  geom_point() + 
  geom_line() +  
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.1) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "black") +
  xlab("Sequence length") + 
  ylab("Coverage") +
  ggtitle("Coverage Emission Probabilities Clarity Condition - Reduced") + 
  labs(color = "Clarity") +
  geom_hline(yintercept = 0) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 
# Combined plot. 
cov_clarity_theta_combined_plot <- grid.arrange(clarity_theta$cov$cov_complete_plot, 
                                                clarity_theta_cov_reduced_plot, ncol = 2)

  }

### Calculate absolute bias, relative bias, standard deviation, and coverage of transition probabilities of the clarity condition, 
### with associated plots. 

## Specify number of observations. 
nobs <- c(5, 5, 5)
## Specify condition.
cond <- c("theta_cl", "theta_modcl", "theta_uncl")
condition <- c("Clear", "Moderately clear", "Unclear")
cond_string <- c("Clarity")
## Calculate the various measures and the associated plots. 
clarity_gamma <- calcR_gamma(true_probs = gamma_true, nobs = nobs, cond = cond, condition = condition, cond_string = cond_string)

### Reduced and combined plots of the absolute bias, relative bias, standard deviation, and coverage of transition probabilities of the 
### clarity condition, 
{
## Absolute bias. 
# Reduced plot. 
clarity_gamma$abs_bias$abs_bias$cats <- car::recode(clarity_gamma$abs_bias$abs_bias$cats, 
                                                    "c('M_S1_to_S1', 'M_S2_to_S2', 'M_S3_to_S3') = 'Self-transition'; 
                                                     else = 'State-to-state transition'")
clarity_gamma$abs_bias$abs_bias$cats <- factor(clarity_gamma$abs_bias$abs_bias$cats, 
                                               levels = c("Self-transition", "State-to-state transition"))
clarity_gamma_abs_bias_reduced_plot <- clarity_gamma$abs_bias$abs_bias %>% 
  group_by(cats, length, condition) %>%
  dplyr::summarize(mean = mean(mean_abs_bias), sd = sd(mean_abs_bias)) %>% 
  ggplot(aes(x = length, y = mean, color = condition, group = condition)) +
  facet_wrap(facets = vars(cats)) +
  geom_point() + 
  geom_line() +  
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2) +
  xlab("Sequence length") + 
  ylab("Absolute bias") +
  ggtitle("Absolute Bias Transition Probabilities Clarity Condition - Reduced") + 
  labs(color = "Clarity") +
  geom_hline(yintercept = 0) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 
# Combined plot. 
clarity_gamma_abs_bias_combined_plot <- grid.arrange(clarity_gamma$abs_bias$abs_bias_complete_plot, 
                                                     clarity_gamma_abs_bias_reduced_plot, ncol = 2)

## Relative bias. 
# Reduced plot. 
clarity_gamma$rel_bias$rel_bias$cats <- car::recode(clarity_gamma$rel_bias$rel_bias$cats, 
                                                    "c('M_S1_to_S1', 'M_S2_to_S2', 'M_S3_to_S3') = 'Self-transition'; 
                                                     else = 'State-to-state transition'")
clarity_gamma$rel_bias$rel_bias$cats <- factor(clarity_gamma$rel_bias$rel_bias$cats, 
                                               levels = c("Self-transition", "State-to-state transition"))
clarity_gamma_rel_bias_reduced_plot <- clarity_gamma$rel_bias$rel_bias %>% 
  group_by(cats, length, condition) %>%
  dplyr::summarize(mean = mean(mean_rel_bias), sd = sd(mean_rel_bias)) %>% 
  ggplot(aes(x = length, y = mean, color = condition, group = condition)) +
  facet_wrap(facets = vars(cats)) +
  geom_point() + 
  geom_line() +  
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.1) +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "black") +
  geom_hline(yintercept = -0.1, linetype = "dashed", color = "black") +
  xlab("Sequence length") + 
  ylab("Relative bias") +
  ggtitle("Relative Bias Transition Probabilities Clarity Condition - Reduced") + 
  labs(color = "Clarity") + 
  geom_hline(yintercept = 0) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 
# Combined plot. 
clarity_gamma_rel_bias_combined_plot <- grid.arrange(clarity_gamma$rel_bias$rel_bias_complete_plot, 
                                                     clarity_gamma_rel_bias_reduced_plot, ncol = 2)

## Standard deviation. 
# Reduced plot. 
clarity_gamma$sdev$sdev$cats <- car::recode(clarity_gamma$sdev$sdev$cats, 
                                            "c('sd_S1_to_S1', 'sd_S2_to_S2', 'sd_S3_to_S3') = 'Self-transition'; 
                                             else = 'State-to-state transition'")
clarity_gamma$sdev$sdev$cats <- factor(clarity_gamma$sdev$sdev$cats, 
                                               levels = c("Self-transition", "State-to-state transition"))
clarity_gamma_sdev_reduced_plot <- clarity_gamma$sdev$sdev %>% 
  group_by(cats, length, condition) %>%
  dplyr::summarize(mean = mean(mean_sdev), sd = sd(mean_sdev)) %>% 
  ggplot(aes(x = length, y = mean, color = condition, group = condition)) +
  facet_wrap(facets = vars(cats)) +
  geom_point() + 
  geom_line() +  
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.1) +
  xlab("Sequence length") + 
  ylab("Standard deviation") +
  ggtitle("Standard Deviation Transition Probabilities Clarity Condition - Reduced") + 
  labs(color = "Clarity") +
  geom_hline(yintercept = 0) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 
# Combined plot. 
clarity_gamma_sdev_combined_plot <- grid.arrange(clarity_gamma$sdev$sdev_complete_plot, 
                                                 clarity_gamma_sdev_reduced_plot, ncol = 2)

## Coverage 
# Reduced plot. 
clarity_gamma$cov$cov$cats <- car::recode(clarity_gamma$cov$cov$cats, 
                                                    "c('M_S1_to_S1', 'M_S2_to_S2', 'M_S3_to_S3') = 'Self-transition'; 
                                                     else = 'State-to-state transition'")
clarity_gamma$cov$cov$cats <- factor(clarity_gamma$cov$cov$cats, 
                                               levels = c("Self-transition", "State-to-state transition"))
clarity_gamma_cov_reduced_plot <- clarity_gamma$cov$cov %>% 
  group_by(cats, length, condition) %>%
  dplyr::summarize(mean = mean(mean_cov), sd = sd(mean_cov)) %>% 
  ggplot(aes(x = length, y = mean, color = condition, group = condition)) +
  facet_wrap(facets = vars(cats)) +
  geom_point() + 
  geom_line() +  
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.1) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "black") +
  xlab("Sequence length") + 
  ylab("Coverage") +
  ggtitle("Coverage Transition Probabilities Clarity Condition - Reduced") + 
  labs(color = "Clarity") +
  geom_hline(yintercept = 0) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 
# Combined plot. 
clarity_gamma_cov_combined_plot <- grid.arrange(clarity_gamma$cov$cov_complete_plot, 
                                                clarity_gamma_cov_reduced_plot, ncol = 2)

  }

###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
### Absolute bias, relative bias, and standard error of the emission and transition probabilities overlap condition ###
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###

### Calculate absolute bias, relative bias, standard deviation, and coverage of emission probabilities of the overlap condition, 
### with associated plots. 

## Specify number of observations. 
nobs <- c(5, 5, 5)
## Specify condition.
cond <- c("theta_nooverl", "theta_modoverl", "theta_muchoverl")
condition <- c("None", "Moderate", "Much")
cond_string <- c("Overlap")
## Specify true probabilities.
true_probs <- list(c(0.84, 0.04, 0.04, 0.04, 0.04, 0.04, 0.44, 0.44, 0.04, 0.04, 0.04, 0.04, 0.04, 0.44, 0.44), 
                   c(0.59, 0.29, 0.04, 0.04, 0.04, 0.04, 0.29, 0.59, 0.04, 0.04, 0.04, 0.04, 0.20, 0.36, 0.36), 
                   c(0.44, 0.44, 0.04, 0.04, 0.04, 0.04, 0.44, 0.44, 0.04, 0.04, 0.04, 0.04, 0.30, 0.31, 0.31))
## Calculate the various measures and the associated plots. 
overlap_theta <- calcR_theta(true_probs = true_probs, nobs = nobs, cond = cond, condition = condition, cond_string = cond_string)

### Reduced and combined plots of the absolute bias, relative bias, standard deviation, and coverage of emission probabilities of the 
### overlap condition,
{
## Absolute bias. 
# Reduced plot. 
temp <- cbind(expand.grid(levels(overlap_theta$abs_bias$abs_bias$cats), levels(overlap_theta$abs_bias$abs_bias$condition)), unlist(true_probs)) 
temp <- cbind(temp, c("One-indicator state", "Noise", "Noise", "Noise", "Noise", "Noise", "Two-indicator state", "Two-indicator state",
                      "Noise", "Noise", "Noise", "Noise", "Noise", "Two-indicator state", "Two-indicator state", "Two-indicator primary state", 
                      "Two-indicator secondary state", "Noise", "Noise", "Noise", "Noise", "Two-indicator secondary state", "Two-indicator primary state", 
                      "Noise", "Noise", "Noise",  "Noise", "Three-indicator secondary state", "Three-indicator primary state", "Three-indicator primary state", 
                      "Two-indicator state", "Two-indicator state", "Noise", "Noise", "Noise", "Noise", "Two-indicator state", "Two-indicator state", 
                      "Noise", "Noise", "Noise", "Noise", "Three-indicator state", "Three-indicator state", "Three-indicator state"))
colnames(temp) <- c("cats", "condition", "true_probs", "cats_n")
overlap_theta$abs_bias$abs_bias <- full_join(overlap_theta$abs_bias$abs_bias, temp, by = c("cats", "condition"))
overlap_theta$abs_bias$abs_bias$cats_n <- factor(overlap_theta$abs_bias$abs_bias$cats_n, 
                                                 levels = c("One-indicator state", "Two-indicator state", "Three-indicator state", 
                                                            "Two-indicator primary state", "Two-indicator secondary state", 
                                                            "Three-indicator primary state", "Three-indicator secondary state",
                                                            "Noise"))
overlap_theta_abs_bias_reduced_plot <- overlap_theta$abs_bias$abs_bias %>%
  group_by(cats_n, length, condition) %>%
  dplyr::summarize(mean = mean(mean_abs_bias), sd = sd(mean_abs_bias)) %>% 
  ggplot(aes(x = length, y = mean, color = condition, group = condition)) +
  facet_wrap(facets = vars(cats_n)) +
  geom_point() + 
  geom_line() +  
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.1) +
  xlab("Sequence length") + 
  ylab("Absolute bias") +
  ggtitle("Absolute Bias Emission Probabilities Overlap Condition - Reduced") + 
  labs(color = "Overlap") +
  geom_hline(yintercept = 0) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 
# Combined plot. 
overlap_theta_abs_bias_combined_plot <- grid.arrange(overlap_theta$abs_bias$abs_bias_complete_plot, 
                                                     overlap_theta_abs_bias_reduced_plot, ncol = 2)

## Relative bias. 
# Reduced plot. 
temp <- cbind(expand.grid(levels(overlap_theta$rel_bias$rel_bias$cats), levels(overlap_theta$rel_bias$rel_bias$condition)), unlist(true_probs)) 
temp <- cbind(temp, c("One-indicator state", "Noise", "Noise", "Noise", "Noise", "Noise", "Two-indicator state", "Two-indicator state",
                      "Noise", "Noise", "Noise", "Noise", "Noise", "Two-indicator state", "Two-indicator state", "Two-indicator primary state", 
                      "Two-indicator secondary state", "Noise", "Noise", "Noise", "Noise", "Two-indicator secondary state", "Two-indicator primary state", 
                      "Noise", "Noise", "Noise",  "Noise", "Three-indicator secondary state", "Three-indicator primary state", "Three-indicator primary state", 
                      "Two-indicator state", "Two-indicator state", "Noise", "Noise", "Noise", "Noise", "Two-indicator state", "Two-indicator state", 
                      "Noise", "Noise", "Noise", "Noise", "Three-indicator state", "Three-indicator state", "Three-indicator state"))
colnames(temp) <- c("cats", "condition", "true_probs", "cats_n")
overlap_theta$rel_bias$rel_bias <- full_join(overlap_theta$rel_bias$rel_bias, temp, by = c("cats", "condition"))
overlap_theta$rel_bias$rel_bias$cats_n <- factor(overlap_theta$rel_bias$rel_bias$cats_n, 
                                                 levels = c("One-indicator state", "Two-indicator state", "Three-indicator state", 
                                                            "Two-indicator primary state", "Two-indicator secondary state", 
                                                            "Three-indicator primary state", "Three-indicator secondary state",
                                                            "Noise"))
overlap_theta_rel_bias_reduced_plot <- overlap_theta$rel_bias$rel_bias %>%
  group_by(cats_n, length, condition) %>%
  dplyr::summarize(mean = mean(mean_rel_bias), sd = sd(mean_rel_bias)) %>% 
  ggplot(aes(x = length, y = mean, color = condition, group = condition)) +
  facet_wrap(facets = vars(cats_n)) +
  geom_point() + 
  geom_line() +  
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.1) +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "black") +
  geom_hline(yintercept = -0.1, linetype = "dashed", color = "black") +
  xlab("Sequence length") + 
  ylab("Relative bias") +
  ggtitle("Relative Bias Emission Probabilities Overlap Condition - Reduced") + 
  labs(color = "Overlap") +
  geom_hline(yintercept = 0) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 
# Combined plot. 
overlap_theta_rel_bias_combined_plot <- grid.arrange(overlap_theta$rel_bias$rel_bias_complete_plot, 
                                                     overlap_theta_rel_bias_reduced_plot, ncol = 2)

## Standard deviation. 
# Reduced plot. 
temp <- cbind(expand.grid(levels(overlap_theta$sdev$sdev$cats), levels(overlap_theta$sdev$sdev$condition)), unlist(true_probs)) 
temp <- cbind(temp, c("One-indicator state", "Noise", "Noise", "Noise", "Noise", "Noise", "Two-indicator state", "Two-indicator state",
                      "Noise", "Noise", "Noise", "Noise", "Noise", "Two-indicator state", "Two-indicator state", "Two-indicator primary state", 
                      "Two-indicator secondary state", "Noise", "Noise", "Noise", "Noise", "Two-indicator secondary state", "Two-indicator primary state", 
                      "Noise", "Noise", "Noise",  "Noise", "Three-indicator secondary state", "Three-indicator primary state", "Three-indicator primary state", 
                      "Two-indicator state", "Two-indicator state", "Noise", "Noise", "Noise", "Noise", "Two-indicator state", "Two-indicator state", 
                      "Noise", "Noise", "Noise", "Noise", "Three-indicator state", "Three-indicator state", "Three-indicator state"))
colnames(temp) <- c("cats", "condition", "true_probs", "cats_n")
overlap_theta$sdev$sdev <- full_join(overlap_theta$sdev$sdev, temp, by = c("cats", "condition"))
overlap_theta$sdev$sdev$cats_n <- factor(overlap_theta$sdev$sdev$cats_n, 
                                                 levels = c("One-indicator state", "Two-indicator state", "Three-indicator state", 
                                                            "Two-indicator primary state", "Two-indicator secondary state", 
                                                            "Three-indicator primary state", "Three-indicator secondary state",
                                                            "Noise"))
overlap_theta_sdev_reduced_plot <- overlap_theta$sdev$sdev %>%
  group_by(cats_n, length, condition) %>%
  dplyr::summarize(mean = mean(mean_sdev), sd = sd(mean_sdev)) %>% 
  ggplot(aes(x = length, y = mean, color = condition, group = condition)) +
  facet_wrap(facets = vars(cats_n)) +
  geom_point() + 
  geom_line() +  
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.1) +
  xlab("Sequence length") + 
  ylab("Standard deviation") +
  ggtitle("Standard Deviation Emission Probabilities Overlap Condition - Reduced") + 
  labs(color = "Overlap") +
  geom_hline(yintercept = 0) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 
# Combined plot. 
overlap_theta_sdev_combined_plot <- grid.arrange(overlap_theta$sdev$sdev_complete_plot, 
                                                 overlap_theta_sdev_reduced_plot, ncol = 2)

## Coverage 
# Reduced plot. 
temp <- cbind(expand.grid(levels(overlap_theta$cov$cov$cats), levels(overlap_theta$cov$cov$condition)), unlist(true_probs)) 
temp <- cbind(temp, c("One-indicator state", "Noise", "Noise", "Noise", "Noise", "Noise", "Two-indicator state", "Two-indicator state",
                      "Noise", "Noise", "Noise", "Noise", "Noise", "Two-indicator state", "Two-indicator state", "Two-indicator primary state", 
                      "Two-indicator secondary state", "Noise", "Noise", "Noise", "Noise", "Two-indicator secondary state", "Two-indicator primary state", 
                      "Noise", "Noise", "Noise",  "Noise", "Three-indicator secondary state", "Three-indicator primary state", "Three-indicator primary state", 
                      "Two-indicator state", "Two-indicator state", "Noise", "Noise", "Noise", "Noise", "Two-indicator state", "Two-indicator state", 
                      "Noise", "Noise", "Noise", "Noise", "Three-indicator state", "Three-indicator state", "Three-indicator state"))
colnames(temp) <- c("cats", "condition", "true_probs", "cats_n")
overlap_theta$cov$cov <- full_join(overlap_theta$cov$cov, temp, by = c("cats", "condition"))
overlap_theta$cov$cov$cats_n <- factor(overlap_theta$cov$cov$cats_n, 
                                       levels = c("One-indicator state", "Two-indicator state", "Three-indicator state", 
                                                  "Two-indicator primary state", "Two-indicator secondary state", 
                                                  "Three-indicator primary state", "Three-indicator secondary state",
                                                  "Noise"))
overlap_theta_cov_reduced_plot <- overlap_theta$cov$cov %>% 
  group_by(cats_n, length, condition) %>%
  dplyr::summarize(mean = mean(mean_cov), sd = sd(mean_cov)) %>% 
  ggplot(aes(x = length, y = mean, color = condition, group = condition)) +
  facet_wrap(facets = vars(cats_n)) +
  geom_point() + 
  geom_line() +  
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.1) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "black") +
  xlab("Sequence length") + 
  ylab("Coverage") +
  ggtitle("Coverage Emission Probabilities Overlap Condition - Reduced") + 
  labs(color = "Overlap") +
  geom_hline(yintercept = 0) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 
# Combined plot. 
cov_clarity_theta_combined_plot <- grid.arrange(overlap_theta$cov$cov_complete_plot, 
                                                overlap_theta_cov_reduced_plot, ncol = 2)

  }

### Calculate absolute bias, relative bias, standard deviation, and coverage of transition probabilities of the overlap condition, 
### with associated plots. 

## Specify number of observations. 
nobs <- c(5, 5, 5)
## Specify condition.
cond <- c("theta_nooverl", "theta_modoverl", "theta_muchoverl")
condition <- c("None", "Moderate", "Much")
cond_string <- c("Overlap")
## Calculate the various measures and the associated plots. 
overlap_gamma <- calcR_gamma(true_probs = gamma_true, nobs = nobs, cond = cond, condition = condition, cond_string = cond_string)

### Reduced and combined plots of the absolute bias, relative bias, standard deviation, and coverage of transition probabilities of the 
### overlap condition,
{
## Absolute bias. 
# Reduced plot. 
overlap_gamma$abs_bias$abs_bias$cats <- car::recode(overlap_gamma$abs_bias$abs_bias$cats, 
                                                    "c('M_S1_to_S1', 'M_S2_to_S2', 'M_S3_to_S3') = 'Self-transition'; 
                                                     else = 'State-to-state transition'")
overlap_gamma$abs_bias$abs_bias$cats <- factor(overlap_gamma$abs_bias$abs_bias$cats, 
                                               levels = c("Self-transition", "State-to-state transition"))
overlap_gamma_abs_bias_reduced_plot <- overlap_gamma$abs_bias$abs_bias %>%
  group_by(cats, length, condition) %>%
  dplyr::summarize(mean = mean(mean_abs_bias), sd = sd(mean_abs_bias)) %>% 
  ggplot(aes(x = length, y = mean, color = condition, group = condition)) +
  facet_wrap(facets = vars(cats)) +
  geom_point() + 
  geom_line() +  
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.1) +
  xlab("Sequence length") + 
  ylab("Absolute bias") +
  ggtitle("Absolute Bias Transition Probabilities Overlap Condition - Reduced") + 
  labs(color = "Overlap") +
  geom_hline(yintercept = 0) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 
# Combined plot. 
overlap_gamma_abs_bias_combined_plot <- grid.arrange(overlap_gamma$abs_bias$abs_bias_complete_plot, 
                                                     overlap_gamma_abs_bias_reduced_plot, ncol = 2)

## Relative bias. 
# Reduced plot. 
overlap_gamma$rel_bias$rel_bias$cats <- car::recode(overlap_gamma$rel_bias$rel_bias$cats, 
                                                    "c('M_S1_to_S1', 'M_S2_to_S2', 'M_S3_to_S3') = 'Self-transition'; 
                                                     else = 'State-to-state transition'")
overlap_gamma$rel_bias$rel_bias$cats <- factor(overlap_gamma$rel_bias$rel_bias$cats, 
                                               levels = c("Self-transition", "State-to-state transition"))
overlap_gamma_rel_bias_reduced_plot <- overlap_gamma$rel_bias$rel_bias %>%
  group_by(cats, length, condition) %>%
  dplyr::summarize(mean = mean(mean_rel_bias), sd = sd(mean_rel_bias)) %>% 
  ggplot(aes(x = length, y = mean, color = condition, group = condition)) +
  facet_wrap(facets = vars(cats)) +
  geom_point() + 
  geom_line() +  
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.1) +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "black") +
  geom_hline(yintercept = -0.1, linetype = "dashed", color = "black") +
  xlab("Sequence length") + 
  ylab("Relative bias") +
  ggtitle("Relative Bias Transition Probabilities Overlap Condition - Reduced") + 
  labs(color = "Overlap") +
  geom_hline(yintercept = 0) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 
# Combined plot. 
overlap_gamma_rel_bias_combined_plot <- grid.arrange(overlap_gamma$rel_bias$rel_bias_complete_plot, 
                                                     overlap_gamma_rel_bias_reduced_plot, ncol = 2)

## Standard deviation.  
# Reduced plot. 
overlap_gamma$sdev$sdev$cats <- car::recode(overlap_gamma$sdev$sdev$cats, 
                                            "c('sd_S1_to_S1', 'sd_S2_to_S2', 'sd_S3_to_S3') = 'Self-transition'; 
                                                     else = 'State-to-state transition'")
overlap_gamma$sdev$sdev$cats <- factor(overlap_gamma$sdev$sdev$cats, 
                                       levels = c("Self-transition", "State-to-state transition"))
overlap_gamma_sdev_reduced_plot <- overlap_gamma$sdev$sdev %>%
  group_by(cats, length, condition) %>%
  dplyr::summarize(mean = mean(mean_sdev), sd = sd(mean_sdev)) %>% 
  ggplot(aes(x = length, y = mean, color = condition, group = condition)) +
  facet_wrap(facets = vars(cats)) +
  geom_point() + 
  geom_line() +  
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.1) +
  xlab("Sequence length") + 
  ylab("Standard deviation") +
  ggtitle("Standard Deviation Transition Probabilities Overlap Condition - Reduced") + 
  labs(color = "Overlap") +
  geom_hline(yintercept = 0) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 
# Combined plot. 
overlap_gamma_sdev_combined_plot <- grid.arrange(overlap_gamma$sdev$sdev_complete_plot, 
                                                 overlap_gamma_sdev_reduced_plot, ncol = 2)

## Coverage 
# Reduced plot. 
overlap_gamma$cov$cov$cats <- car::recode(overlap_gamma$cov$cov$cats, 
                                          "c('M_S1_to_S1', 'M_S2_to_S2', 'M_S3_to_S3') = 'Self-transition'; 
                                                     else = 'State-to-state transition'")
overlap_gamma$cov$cov$cats <- factor(overlap_gamma$cov$cov$cats, 
                                     levels = c("Self-transition", "State-to-state transition"))
overlap_gamma_cov_reduced_plot <- overlap_gamma$cov$cov %>% 
  group_by(cats, length, condition) %>%
  dplyr::summarize(mean = mean(mean_cov), sd = sd(mean_cov)) %>% 
  ggplot(aes(x = length, y = mean, color = condition, group = condition)) +
  facet_wrap(facets = vars(cats)) +
  geom_point() + 
  geom_line() +  
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.1) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "black") +
  xlab("Sequence length") + 
  ylab("Coverage") +
  ggtitle("Coverage Transition Probabilities Overlap Condition - Reduced") + 
  labs(color = "Overlap") +
  geom_hline(yintercept = 0) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 
# Combined plot. 
overlap_gamma_cov_combined_plot <- grid.arrange(overlap_gamma$cov$cov_complete_plot, 
                                                overlap_gamma_cov_reduced_plot, ncol = 2)

  }

###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
### Absolute bias, relative bias, and standard error of the emission and transition probabilities unclear number of observations condition ###
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###

### Calculate absolute bias, relative bias, standard deviation, and coverage of the emission probabilities for the unclear number 
### of observations condition, with associated plots.

## Specify number of observations. 
nobs <- c(3, 5, 7)
## Specify condition.
cond <- c("theta_uncl", "theta_uncl", "theta_uncl")
condition <- c("3", "5", "7")
cond_string <- c("Unclear Number of Observations")
## Specify true probabilities.
true_probs <- list(c(0.52, 0.24, 0.24, 0.24, 0.52, 0.24, 0.24, 0.24, 0.52), 
                   c(0.44, 0.14, 0.14, 0.14, 0.14, 0.14, 0.29, 0.14, 0.29, 0.14, 0.14, 0.14, 0.29, 0.14, 0.29), 
                   c(0.40, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.25, 0.10, 0.25, 0.10, 0.10, 0.10, 0.08, 0.08, 0.22, 0.09, 0.22, 0.09, 0.22))
## Calculate the various measures and the associated plots. 
nobs_uncl_theta <- calcR_theta(true_probs = true_probs, nobs = nobs, cond = cond, condition = condition, cond_string = cond_string, varobs = T)

### Reduced and combined plots of the absolute bias, relative bias, standard deviation, and coverage of the emission probabilities for the number 
### of observations condition with unclear emission probabilities, with associated plots. 
{
## Absolute bias. 
# Reduced plot. 
foo1 <- cbind(expand.grid(c(paste0("M_S1_cat", 1:3), paste0("M_S2_cat", 1:3), paste0("M_S3_cat", 1:3)), levels(nobs_uncl_theta$abs_bias$abs_bias$condition)[1]), unlist(true_probs[1]))
foo2 <- cbind(expand.grid(c(paste0("M_S1_cat", 1:5), paste0("M_S2_cat", 1:5), paste0("M_S3_cat", 1:5)), levels(nobs_uncl_theta$abs_bias$abs_bias$condition)[2]), unlist(true_probs[2]))
foo3 <- cbind(expand.grid(c(paste0("M_S1_cat", 1:7), paste0("M_S2_cat", 1:7), paste0("M_S3_cat", 1:7)), levels(nobs_uncl_theta$abs_bias$abs_bias$condition)[3]), unlist(true_probs[3]))
colnames(foo1) <- c("cats", "condition", "true")
colnames(foo2) <- c("cats", "condition", "true")
colnames(foo3) <- c("cats", "condition", "true")
temp <- rbind(foo1, foo2, foo3)
temp <- cbind(temp, c("One-indicator state", "Noise", "Noise", "Noise", "One-indicator state", "Noise", "Noise", "Noise", "One-indicator state",
                      "One-indicator state", "Noise", "Noise", "Noise", "Noise", "Noise", "Two-indicator state", "Noise", "Two-indicator state", "Noise",
                      "Noise", "Noise", "Two-indicator state", "Noise", "Two-indicator state", "One-indicator state", "Noise", "Noise", "Noise", "Noise", "Noise",
                      "Noise", "Noise", "Two-indicator state", "Noise", "Two-indicator state",  "Noise", "Noise",  "Noise", "Noise", "Noise", "Three-indicator state",
                      "Noise", "Three-indicator state", "Noise", "Three-indicator state"))
colnames(temp) <- c("cats", "condition", "true_probs", "cats_n")
nobs_uncl_theta$abs_bias$abs_bias <- full_join(nobs_uncl_theta$abs_bias$abs_bias, temp, by = c("cats", "condition"))
nobs_uncl_theta$abs_bias$abs_bias$cats_n <- factor(nobs_uncl_theta$abs_bias$abs_bias$cats_n, 
                                                   levels = c("One-indicator state", "Two-indicator state", "Three-indicator state", 
                                                              "Noise"))
nobs_uncl_theta_abs_bias_reduced_plot <- nobs_uncl_theta$abs_bias$abs_bias %>%
  group_by(cats_n, length, condition) %>%
  dplyr::summarize(mean = mean(mean_abs_bias), sd = sd(mean_abs_bias)) %>% 
  ggplot(aes(x = length, y = mean, color = condition, group = condition)) +
  facet_wrap(facets = vars(cats_n)) +
  geom_point() + 
  geom_line() +  
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.1) +
  xlab("Sequence length") + 
  ylab("Absolute bias") +
  ggtitle("Absolute Bias Emission Probabilities Unclear Number of Observations Condition - Reduced") + 
  labs(color = "Unclear Number of Observations") +
  geom_hline(yintercept = 0) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 
# Combined plot. 
nobs_uncl_theta_abs_bias_combined_plot <- grid.arrange(nobs_uncl_theta$abs_bias$abs_bias_complete_plot, 
                                                       nobs_uncl_theta_abs_bias_reduced_plot, ncol = 2)

## Relative bias. 
# Reduced plot. 
foo1 <- cbind(expand.grid(c(paste0("M_S1_cat", 1:3), paste0("M_S2_cat", 1:3), paste0("M_S3_cat", 1:3)), levels(nobs_uncl_theta$rel_bias$rel_bias$condition)[1]), unlist(true_probs[1]))
foo2 <- cbind(expand.grid(c(paste0("M_S1_cat", 1:5), paste0("M_S2_cat", 1:5), paste0("M_S3_cat", 1:5)), levels(nobs_uncl_theta$rel_bias$rel_bias$condition)[2]), unlist(true_probs[2]))
foo3 <- cbind(expand.grid(c(paste0("M_S1_cat", 1:7), paste0("M_S2_cat", 1:7), paste0("M_S3_cat", 1:7)), levels(nobs_uncl_theta$rel_bias$rel_bias$condition)[3]), unlist(true_probs[3]))
colnames(foo1) <- c("cats", "condition", "true")
colnames(foo2) <- c("cats", "condition", "true")
colnames(foo3) <- c("cats", "condition", "true")
temp <- rbind(foo1, foo2, foo3)
temp <- cbind(temp, c("One-indicator state", "Noise", "Noise", "Noise", "One-indicator state", "Noise", "Noise", "Noise", "One-indicator state",
                      "One-indicator state", "Noise", "Noise", "Noise", "Noise", "Noise", "Two-indicator state", "Noise", "Two-indicator state", "Noise",
                      "Noise", "Noise", "Two-indicator state", "Noise", "Two-indicator state", "One-indicator state", "Noise", "Noise", "Noise", "Noise", "Noise",
                      "Noise", "Noise", "Two-indicator state", "Noise", "Two-indicator state",  "Noise", "Noise",  "Noise", "Noise", "Noise", "Three-indicator state",
                      "Noise", "Three-indicator state", "Noise", "Three-indicator state"))
colnames(temp) <- c("cats", "condition", "true_probs", "cats_n")
nobs_uncl_theta$rel_bias$rel_bias <- full_join(nobs_uncl_theta$rel_bias$rel_bias, temp, by = c("cats", "condition"))
nobs_uncl_theta$rel_bias$rel_bias$cats_n <- factor(nobs_uncl_theta$rel_bias$rel_bias$cats_n, 
                                                   levels = c("One-indicator state", "Two-indicator state", "Three-indicator state", 
                                                              "Noise"))
nobs_uncl_theta_rel_bias_reduced_plot <- nobs_uncl_theta$rel_bias$rel_bias %>%
  group_by(cats_n, length, condition) %>%
  dplyr::summarize(mean = mean(mean_rel_bias), sd = sd(mean_rel_bias)) %>% 
  ggplot(aes(x = length, y = mean, color = condition, group = condition)) +
  facet_wrap(facets = vars(cats_n)) +
  geom_point() + 
  geom_line() +  
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.1) +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "black") +
  geom_hline(yintercept = -0.1, linetype = "dashed", color = "black") +
  xlab("Sequence length") + 
  ylab("Relative bias") +
  ggtitle("Relative Bias Emission Probabilities Unclear Number of Observations Condition - Reduced") + 
  labs(color = "Unclear Number of Observations") +
  geom_hline(yintercept = 0) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 
# Combined plot. 
nobs_uncl_theta_rel_bias_combined_plot <- grid.arrange(nobs_uncl_theta$rel_bias$rel_bias_complete_plot, 
                                                       nobs_uncl_theta_rel_bias_reduced_plot, ncol = 2)

## Standard deviation. 
# Reduced plot. 
foo1 <- cbind(expand.grid(c(paste0("sd_S1_cat", 1:3), paste0("sd_S2_cat", 1:3), paste0("sd_S3_cat", 1:3)), levels(nobs_uncl_theta$sdev$sdev$condition)[1]), unlist(true_probs[1]))
foo2 <- cbind(expand.grid(c(paste0("sd_S1_cat", 1:5), paste0("sd_S2_cat", 1:5), paste0("sd_S3_cat", 1:5)), levels(nobs_uncl_theta$sdev$sdev$condition)[2]), unlist(true_probs[2]))
foo3 <- cbind(expand.grid(c(paste0("sd_S1_cat", 1:7), paste0("sd_S2_cat", 1:7), paste0("sd_S3_cat", 1:7)), levels(nobs_uncl_theta$sdev$sdev$condition)[3]), unlist(true_probs[3]))
colnames(foo1) <- c("cats", "condition", "true")
colnames(foo2) <- c("cats", "condition", "true")
colnames(foo3) <- c("cats", "condition", "true")
temp <- rbind(foo1, foo2, foo3)
temp <- cbind(temp, c("One-indicator state", "Noise", "Noise", "Noise", "One-indicator state", "Noise", "Noise", "Noise", "One-indicator state",
                      "One-indicator state", "Noise", "Noise", "Noise", "Noise", "Noise", "Two-indicator state", "Noise", "Two-indicator state", "Noise",
                      "Noise", "Noise", "Two-indicator state", "Noise", "Two-indicator state", "One-indicator state", "Noise", "Noise", "Noise", "Noise", "Noise",
                      "Noise", "Noise", "Two-indicator state", "Noise", "Two-indicator state",  "Noise", "Noise",  "Noise", "Noise", "Noise", "Three-indicator state",
                      "Noise", "Three-indicator state", "Noise", "Three-indicator state"))
colnames(temp) <- c("cats", "condition", "true_probs", "cats_n")
nobs_uncl_theta$sdev$sdev <- full_join(nobs_uncl_theta$sdev$sdev, temp, by = c("cats", "condition"))
nobs_uncl_theta$sdev$sdev$cats_n <- factor(nobs_uncl_theta$sdev$sdev$cats_n, 
                                                   levels = c("One-indicator state", "Two-indicator state", "Three-indicator state", 
                                                              "Noise"))
nobs_uncl_theta_sdev_reduced_plot <- nobs_uncl_theta$sdev$sdev %>%
  group_by(cats_n, length, condition) %>%
  dplyr::summarize(mean = mean(mean_sdev), sd = sd(mean_sdev)) %>% 
  ggplot(aes(x = length, y = mean, color = condition, group = condition)) +
  facet_wrap(facets = vars(cats_n)) +
  geom_point() + 
  geom_line() +  
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.1) +
  xlab("Sequence length") + 
  ylab("Standard deviation") +
  ggtitle("Standard Deviation Emission Probabilities Unclear Number of Observations Condition - Reduced") + 
  labs(color = "Unclear Number of Observations") +
  geom_hline(yintercept = 0) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 
# Combined plot. 
nobs_uncl_theta_sdev_combined_plot <- grid.arrange(nobs_uncl_theta$sdev$sdev_complete_plot, 
                                                   nobs_uncl_theta_sdev_reduced_plot, ncol = 2)

## Coverage. 
# Reduced plot. 
foo1 <- cbind(expand.grid(c(paste0("M_S1_cat", 1:3), paste0("M_S2_cat", 1:3), paste0("M_S3_cat", 1:3)), levels(nobs_uncl_theta$cov$cov$condition)[1]), unlist(true_probs[1]))
foo2 <- cbind(expand.grid(c(paste0("M_S1_cat", 1:5), paste0("M_S2_cat", 1:5), paste0("M_S3_cat", 1:5)), levels(nobs_uncl_theta$cov$cov$condition)[2]), unlist(true_probs[2]))
foo3 <- cbind(expand.grid(c(paste0("M_S1_cat", 1:7), paste0("M_S2_cat", 1:7), paste0("M_S3_cat", 1:7)), levels(nobs_uncl_theta$cov$cov$condition)[3]), unlist(true_probs[3]))
colnames(foo1) <- c("cats", "condition", "true")
colnames(foo2) <- c("cats", "condition", "true")
colnames(foo3) <- c("cats", "condition", "true")
temp <- rbind(foo1, foo2, foo3)
temp <- cbind(temp, c("One-indicator state", "Noise", "Noise", "Noise", "One-indicator state", "Noise", "Noise", "Noise", "One-indicator state",
                      "One-indicator state", "Noise", "Noise", "Noise", "Noise", "Noise", "Two-indicator state", "Noise", "Two-indicator state", "Noise",
                      "Noise", "Noise", "Two-indicator state", "Noise", "Two-indicator state", "One-indicator state", "Noise", "Noise", "Noise", "Noise", "Noise",
                      "Noise", "Noise", "Two-indicator state", "Noise", "Two-indicator state",  "Noise", "Noise",  "Noise", "Noise", "Noise", "Three-indicator state",
                      "Noise", "Three-indicator state", "Noise", "Three-indicator state"))
colnames(temp) <- c("cats", "condition", "true_probs", "cats_n")
nobs_uncl_theta$cov$cov <- full_join(nobs_uncl_theta$cov$cov, temp, by = c("cats", "condition"))
nobs_uncl_theta$cov$cov$cats_n <- factor(nobs_uncl_theta$cov$cov$cats_n, 
                                                   levels = c("One-indicator state", "Two-indicator state", "Three-indicator state", 
                                                              "Noise"))
nobs_uncl_theta_cov_reduced_plot <- nobs_uncl_theta$cov$cov %>%
  group_by(cats_n, length, condition) %>%
  dplyr::summarize(mean = mean(mean_cov), sd = sd(mean_cov)) %>% 
  ggplot(aes(x = length, y = mean, color = condition, group = condition)) +
  facet_wrap(facets = vars(cats_n)) +
  geom_point() + 
  geom_line() +  
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.1) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "black") +
  xlab("Sequence length") + 
  ylab("Absolute bias") +
  ggtitle("Relative Bias Emission Probabilities Unclear Number of Observations Condition - Reduced") + 
  labs(color = "Unclear Number of Observations") +
  geom_hline(yintercept = 0) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 
# Combined plot. 
nobs_uncl_theta_cov_combined_plot <- grid.arrange(nobs_uncl_theta$cov$cov_complete_plot, 
                                                  nobs_uncl_theta_cov_reduced_plot, ncol = 2)

  }

### Calculate absolute bias, relative bias, standard deviation, and coverage of the transition probabilities for the unclear number 
### of observations condition, with associated plots. 

## Specify number of observations. 
nobs <- c(3, 5, 7)
## Specify condition.
cond <- c("theta_uncl", "theta_uncl", "theta_uncl")
condition <- c("3", "5", "7")
cond_string <- c("Unclear Number of Observations")
## Calculate the various measures and the associated plots. 
nobs_uncl_gamma <- calcR_gamma(true_probs = gamma_true, nobs = nobs, cond = cond, condition = condition, cond_string = cond_string, varobs = T)

### Reduced and combined plots of the absolute bias, relative bias, standard deviation, and coverage of the emission probabilities for the unclear  
### number of observations condition, with associated plots. 
{
## Absolute bias. 
# Reduced plot. 
nobs_uncl_gamma$abs_bias$abs_bias$cats <- car::recode(nobs_uncl_gamma$abs_bias$abs_bias$cats, 
                                                      "c('M_S1_to_S1', 'M_S2_to_S2', 'M_S3_to_S3') = 'Self-transition'; 
                                                       else = 'State-to-state transition'")
nobs_uncl_gamma$abs_bias$abs_bias$cats <- factor(nobs_uncl_gamma$abs_bias$abs_bias$cats, 
                                                 levels = c("Self-transition", "State-to-state transition"))
nobs_uncl_gamma_abs_bias_reduced_plot <- nobs_uncl_gamma$abs_bias$abs_bias %>%
  group_by(cats, length, condition) %>%
  dplyr::summarize(mean = mean(mean_abs_bias), sd = sd(mean_abs_bias)) %>% 
  ggplot(aes(x = length, y = mean, color = condition, group = condition)) +
  facet_wrap(facets = vars(cats)) +
  geom_point() + 
  geom_line() +  
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.1) +
  xlab("Sequence length") + 
  ylab("Absolute bias") +
  ggtitle("Absolute Bias Transition Probabilities Unclear Number of Observations Condition - Reduced") + 
  labs(color = "Unclear Number of Observations") +
  geom_hline(yintercept = 0) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 
# Combined plot. 
nobs_uncl_gamma_abs_bias_combined_plot <- grid.arrange(nobs_uncl_gamma$abs_bias$abs_bias_complete_plot, 
                                                       nobs_uncl_gamma_abs_bias_reduced_plot, ncol = 2)

## Relative bias. 
# Reduced plot. 
nobs_uncl_gamma$rel_bias$rel_bias$cats <- car::recode(nobs_uncl_gamma$rel_bias$rel_bias$cats, 
                                                      "c('M_S1_to_S1', 'M_S2_to_S2', 'M_S3_to_S3') = 'Self-transition'; 
                                                       else = 'State-to-state transition'")
nobs_uncl_gamma$rel_bias$rel_bias$cats <- factor(nobs_uncl_gamma$rel_bias$rel_bias$cats, 
                                                 levels = c("Self-transition", "State-to-state transition"))
nobs_uncl_gamma_rel_bias_reduced_plot <- nobs_uncl_gamma$rel_bias$rel_bias %>%
  group_by(cats, length, condition) %>%
  dplyr::summarize(mean = mean(mean_rel_bias), sd = sd(mean_rel_bias)) %>% 
  ggplot(aes(x = length, y = mean, color = condition, group = condition)) +
  facet_wrap(facets = vars(cats)) +
  geom_point() + 
  geom_line() +  
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.1) +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "black") +
  geom_hline(yintercept = -0.1, linetype = "dashed", color = "black") +
  xlab("Sequence length") + 
  ylab("Relative bias") +
  ggtitle("Relative Bias Transition Probabilities Unclear Number of Observations Condition - Reduced") + 
  labs(color = "Unclear Number of Observations") +
  geom_hline(yintercept = 0) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 
# Combined plot. 
nobs_uncl_gamma_rel_bias_combined_plot <- grid.arrange(nobs_uncl_gamma$rel_bias$rel_bias_complete_plot, 
                                                       nobs_uncl_gamma_rel_bias_reduced_plot, ncol = 2)

## Standard deviation. 
# Reduced plot. 
nobs_uncl_gamma$sdev$sdev$cats <- car::recode(nobs_uncl_gamma$sdev$sdev$cats, 
                                                      "c('sd_S1_to_S1', 'sd_S2_to_S2', 'sd_S3_to_S3') = 'Self-transition'; 
                                                       else = 'State-to-state transition'")
nobs_uncl_gamma$sdev$sdev$cats <- factor(nobs_uncl_gamma$sdev$sdev$cats, 
                                                 levels = c("Self-transition", "State-to-state transition"))
nobs_uncl_gamma_sdev_reduced_plot <- nobs_uncl_gamma$sdev$sdev %>%
  group_by(cats, length, condition) %>%
  dplyr::summarize(mean = mean(mean_sdev), sd = sd(mean_sdev)) %>% 
  ggplot(aes(x = length, y = mean, color = condition, group = condition)) +
  facet_wrap(facets = vars(cats)) +
  geom_point() + 
  geom_line() +  
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.1) +
  xlab("Sequence length") + 
  ylab("Standard deviation") +
  ggtitle("Standard Deviation Transition Probabilities Unclear Number of Observations Condition - Reduced") + 
  labs(color = "Unclear Number of Observations") +
  geom_hline(yintercept = 0) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 
# Combined plot. 
nobs_uncl_gamma_sdev_combined_plot <- grid.arrange(nobs_uncl_gamma$sdev$sdev_complete_plot, 
                                                       nobs_uncl_gamma_sdev_reduced_plot, ncol = 2)

## Coverage
# Reduced plot. 
nobs_uncl_gamma$cov$cov$cats <- car::recode(nobs_uncl_gamma$cov$cov$cats, 
                                                      "c('M_S1_to_S1', 'M_S2_to_S2', 'M_S3_to_S3') = 'Self-transition'; 
                                                       else = 'State-to-state transition'")
nobs_uncl_gamma$cov$cov$cats <- factor(nobs_uncl_gamma$cov$cov$cats, 
                                                 levels = c("Self-transition", "State-to-state transition"))
nobs_uncl_gamma_cov_reduced_plot <- nobs_uncl_gamma$cov$cov %>%
  group_by(cats, length, condition) %>%
  dplyr::summarize(mean = mean(mean_cov), sd = sd(mean_cov)) %>% 
  ggplot(aes(x = length, y = mean, color = condition, group = condition)) +
  facet_wrap(facets = vars(cats)) +
  geom_point() + 
  geom_line() +  
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.1) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "black") +
  xlab("Sequence length") + 
  ylab("Coverage") +
  ggtitle("Coverage Transition Probabilities Unclear Number of Observations Condition - Reduced") + 
  labs(color = "Unclear Number of Observations") +
  geom_hline(yintercept = 0) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 
# Combined plot. 
nobs_uncl_gamma_cov_combined_plot <- grid.arrange(nobs_uncl_gamma$cov$cov_complete_plot, 
                                                       nobs_uncl_gamma_cov_reduced_plot, ncol = 2)

  }

###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
### Absolute bias, relative bias, and standard error of the emission and transition probabilities clear number of observations condition ###
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###

### Calculate absolute bias, relative bias, standard deviation, and coverage of the emission probabilities for the clear number 
### of observations condition, with associated plots.

## Specify number of observations. 
nobs <- c(3, 5, 7)
## Specify condition.
cond <- c("varobs", "varobs", "varobs")
condition <- c("3", "5", "7")
cond_string <- c("Clear Number of Observations")
## Specify true probabilities.
true_probs <- list(c(0.80, 0.10, 0.10, 0.10, 0.80, 0.10, 0.10, 0.10, 0.80), 
                   c(0.92, 0.02, 0.02, 0.02, 0.02, 0.02, 0.47, 0.02, 0.47, 0.02, 0.02, 0.02, 0.62, 0.02, 0.32), 
                   c(0.76, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.02, 0.45, 0.02, 0.45, 0.02, 0.02, 0.02, 0.01, 0.02, 0.31, 0.02, 0.31, 0.02, 0.31))
## Calculate the various measures and the associated plots. 
nobs_cl_theta <- calcR_theta(true_probs = true_probs, nobs = nobs, cond = cond, condition = condition, cond_string = cond_string, varobs = T)

### Reduced and combined plot of the absolute bias, relative bias, standard deviation, and coverage of the emission probabilities for the clear number 
### of observations condition, with associated plots.
{
## Absolute bias. 
# Reduced plot. 
foo1 <- cbind(expand.grid(c(paste0("M_S1_cat", 1:3), paste0("M_S2_cat", 1:3), paste0("M_S3_cat", 1:3)), levels(nobs_cl_theta$abs_bias$abs_bias$condition)[1]), unlist(true_probs[1]))
foo2 <- cbind(expand.grid(c(paste0("M_S1_cat", 1:5), paste0("M_S2_cat", 1:5), paste0("M_S3_cat", 1:5)), levels(nobs_cl_theta$abs_bias$abs_bias$condition)[2]), unlist(true_probs[2]))
foo3 <- cbind(expand.grid(c(paste0("M_S1_cat", 1:7), paste0("M_S2_cat", 1:7), paste0("M_S3_cat", 1:7)), levels(nobs_cl_theta$abs_bias$abs_bias$condition)[3]), unlist(true_probs[3]))
colnames(foo1) <- c("cats", "condition", "true")
colnames(foo2) <- c("cats", "condition", "true")
colnames(foo3) <- c("cats", "condition", "true")
temp <- rbind(foo1, foo2, foo3)
temp <- cbind(temp, c("One-indicator state", "Noise", "Noise", "Noise", "One-indicator state", "Noise", "Noise", "Noise", "One-indicator state",
                      "One-indicator state", "Noise", "Noise", "Noise", "Noise", "Noise", "Two-indicator state", "Noise", "Two-indicator state", "Noise",
                      "Noise", "Noise", "Two-indicator primary state", "Noise", "Two-indicator secondary state", "One-indicator state", "Noise", "Noise", "Noise", "Noise", "Noise",
                      "Noise", "Noise", "Two-indicator state", "Noise", "Two-indicator state",  "Noise", "Noise",  "Noise", "Noise", "Noise", "Three-indicator state",
                      "Noise", "Three-indicator state", "Noise", "Three-indicator state"))
colnames(temp) <- c("cats", "condition", "true_probs", "cats_n")
nobs_cl_theta$abs_bias$abs_bias <- full_join(nobs_cl_theta$abs_bias$abs_bias, temp, by = c("cats", "condition"))
nobs_cl_theta$abs_bias$abs_bias$cats_n <- factor(nobs_cl_theta$abs_bias$abs_bias$cats_n, 
                                                   levels = c("One-indicator state", "Two-indicator state", "Three-indicator state", 
                                                              "Two-indicator primary state", "Two-indicator secondary state", "Noise"))
nobs_cl_theta_abs_bias_reduced_plot <- nobs_cl_theta$abs_bias$abs_bias %>%
  group_by(cats_n, length, condition) %>%
  dplyr::summarize(mean = mean(mean_abs_bias), sd = sd(mean_abs_bias)) %>% 
  ggplot(aes(x = length, y = mean, color = condition, group = condition)) +
  facet_wrap(facets = vars(cats_n)) +
  geom_point() + 
  geom_line() +  
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.1) +
  xlab("Sequence length") + 
  ylab("Absolute bias") +
  ggtitle("Absolute Bias Emission Probabilities Clear Number of Observations Condition - Reduced") + 
  labs(color = "Clear Number of Observations") +
  geom_hline(yintercept = 0) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 
# Combined plot. 
nobs_cl_theta_abs_bias_combined_plot <- grid.arrange(nobs_cl_theta$abs_bias$abs_bias_complete_plot, 
                                                       nobs_cl_theta_abs_bias_reduced_plot, ncol = 2)

## Relative bias. 
# Reduced plot. 
foo1 <- cbind(expand.grid(c(paste0("M_S1_cat", 1:3), paste0("M_S2_cat", 1:3), paste0("M_S3_cat", 1:3)), levels(nobs_cl_theta$rel_bias$rel_bias$condition)[1]), unlist(true_probs[1]))
foo2 <- cbind(expand.grid(c(paste0("M_S1_cat", 1:5), paste0("M_S2_cat", 1:5), paste0("M_S3_cat", 1:5)), levels(nobs_cl_theta$rel_bias$rel_bias$condition)[2]), unlist(true_probs[2]))
foo3 <- cbind(expand.grid(c(paste0("M_S1_cat", 1:7), paste0("M_S2_cat", 1:7), paste0("M_S3_cat", 1:7)), levels(nobs_cl_theta$rel_bias$rel_bias$condition)[3]), unlist(true_probs[3]))
colnames(foo1) <- c("cats", "condition", "true")
colnames(foo2) <- c("cats", "condition", "true")
colnames(foo3) <- c("cats", "condition", "true")
temp <- rbind(foo1, foo2, foo3)
temp <- cbind(temp, c("One-indicator state", "Noise", "Noise", "Noise", "One-indicator state", "Noise", "Noise", "Noise", "One-indicator state",
                      "One-indicator state", "Noise", "Noise", "Noise", "Noise", "Noise", "Two-indicator state", "Noise", "Two-indicator state", "Noise",
                      "Noise", "Noise", "Two-indicator primary state", "Noise", "Two-indicator secondary state", "One-indicator state", "Noise", "Noise", "Noise", "Noise", "Noise",
                      "Noise", "Noise", "Two-indicator state", "Noise", "Two-indicator state",  "Noise", "Noise",  "Noise", "Noise", "Noise", "Three-indicator state",
                      "Noise", "Three-indicator state", "Noise", "Three-indicator state"))
colnames(temp) <- c("cats", "condition", "true_probs", "cats_n")
nobs_cl_theta$rel_bias$rel_bias <- full_join(nobs_cl_theta$rel_bias$rel_bias, temp, by = c("cats", "condition"))
nobs_cl_theta$rel_bias$rel_bias$cats_n <- factor(nobs_cl_theta$rel_bias$rel_bias$cats_n, 
                                                 levels = c("One-indicator state", "Two-indicator state", "Three-indicator state", 
                                                            "Two-indicator primary state", "Two-indicator secondary state", "Noise"))
nobs_cl_theta_rel_bias_reduced_plot <- nobs_cl_theta$rel_bias$rel_bias %>%
  group_by(cats_n, length, condition) %>%
  dplyr::summarize(mean = mean(mean_rel_bias), sd = sd(mean_rel_bias)) %>% 
  ggplot(aes(x = length, y = mean, color = condition, group = condition)) +
  facet_wrap(facets = vars(cats_n)) +
  geom_point() + 
  geom_line() +  
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.1) +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "black") +
  geom_hline(yintercept = -0.1, linetype = "dashed", color = "black") +
  xlab("Sequence length") + 
  ylab("Relative bias") +
  ggtitle("Relative Bias Emission Probabilities Clear Number of Observations Condition - Reduced") + 
  labs(color = "Clear Number of Observations") +
  geom_hline(yintercept = 0) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 
# Combined plot. 
nobs_cl_theta_rel_bias_combined_plot <- grid.arrange(nobs_cl_theta$rel_bias$rel_bias_complete_plot, 
                                                       nobs_cl_theta_rel_bias_reduced_plot, ncol = 2)

## Standard deviation. 
# Reduced plot. 
foo1 <- cbind(expand.grid(c(paste0("sd_S1_cat", 1:3), paste0("sd_S2_cat", 1:3), paste0("sd_S3_cat", 1:3)), levels(nobs_cl_theta$sdev$sdev$condition)[1]), unlist(true_probs[1]))
foo2 <- cbind(expand.grid(c(paste0("sd_S1_cat", 1:5), paste0("sd_S2_cat", 1:5), paste0("sd_S3_cat", 1:5)), levels(nobs_cl_theta$sdev$sdev$condition)[2]), unlist(true_probs[2]))
foo3 <- cbind(expand.grid(c(paste0("sd_S1_cat", 1:7), paste0("sd_S2_cat", 1:7), paste0("sd_S3_cat", 1:7)), levels(nobs_cl_theta$sdev$sdev$condition)[3]), unlist(true_probs[3]))
colnames(foo1) <- c("cats", "condition", "true")
colnames(foo2) <- c("cats", "condition", "true")
colnames(foo3) <- c("cats", "condition", "true")
temp <- rbind(foo1, foo2, foo3)
temp <- cbind(temp, c("One-indicator state", "Noise", "Noise", "Noise", "One-indicator state", "Noise", "Noise", "Noise", "One-indicator state",
                      "One-indicator state", "Noise", "Noise", "Noise", "Noise", "Noise", "Two-indicator state", "Noise", "Two-indicator state", "Noise",
                      "Noise", "Noise", "Two-indicator primary state", "Noise", "Two-indicator secondary state", "One-indicator state", "Noise", "Noise", "Noise", "Noise", "Noise",
                      "Noise", "Noise", "Two-indicator state", "Noise", "Two-indicator state",  "Noise", "Noise",  "Noise", "Noise", "Noise", "Three-indicator state",
                      "Noise", "Three-indicator state", "Noise", "Three-indicator state"))
colnames(temp) <- c("cats", "condition", "true_probs", "cats_n")
nobs_cl_theta$sdev$sdev <- full_join(nobs_cl_theta$sdev$sdev, temp, by = c("cats", "condition"))
nobs_cl_theta$sdev$sdev$cats_n <- factor(nobs_cl_theta$sdev$sdev$cats_n, 
                                         levels = c("One-indicator state", "Two-indicator state", "Three-indicator state", 
                                                    "Two-indicator primary state", "Two-indicator secondary state", "Noise"))
nobs_cl_theta_sdev_reduced_plot <- nobs_cl_theta$sdev$sdev %>%
  group_by(cats_n, length, condition) %>%
  dplyr::summarize(mean = mean(mean_sdev), sd = sd(mean_sdev)) %>% 
  ggplot(aes(x = length, y = mean, color = condition, group = condition)) +
  facet_wrap(facets = vars(cats_n)) +
  geom_point() + 
  geom_line() +  
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.1) +
  xlab("Sequence length") + 
  ylab("Standard deviation") +
  ggtitle("Standard Deviation Emission Probabilities Clear Number of Observations Condition - Reduced") + 
  labs(color = "Clear Number of Observations") +
  geom_hline(yintercept = 0) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 
# Combined plot. 
nobs_cl_theta_sdev_combined_plot <- grid.arrange(nobs_cl_theta$sdev$sdev_complete_plot, 
                                                   nobs_cl_theta_sdev_reduced_plot, ncol = 2)

## Coverage. 
# Reduced plot. 
foo1 <- cbind(expand.grid(c(paste0("M_S1_cat", 1:3), paste0("M_S2_cat", 1:3), paste0("M_S3_cat", 1:3)), levels(nobs_cl_theta$cov$cov$condition)[1]), unlist(true_probs[1]))
foo2 <- cbind(expand.grid(c(paste0("M_S1_cat", 1:5), paste0("M_S2_cat", 1:5), paste0("M_S3_cat", 1:5)), levels(nobs_cl_theta$cov$cov$condition)[2]), unlist(true_probs[2]))
foo3 <- cbind(expand.grid(c(paste0("M_S1_cat", 1:7), paste0("M_S2_cat", 1:7), paste0("M_S3_cat", 1:7)), levels(nobs_cl_theta$cov$cov$condition)[3]), unlist(true_probs[3]))
colnames(foo1) <- c("cats", "condition", "true")
colnames(foo2) <- c("cats", "condition", "true")
colnames(foo3) <- c("cats", "condition", "true")
temp <- rbind(foo1, foo2, foo3)
temp <- cbind(temp, c("One-indicator state", "Noise", "Noise", "Noise", "One-indicator state", "Noise", "Noise", "Noise", "One-indicator state",
                      "One-indicator state", "Noise", "Noise", "Noise", "Noise", "Noise", "Two-indicator state", "Noise", "Two-indicator state", "Noise",
                      "Noise", "Noise", "Two-indicator primary state", "Noise", "Two-indicator secondary state", "One-indicator state", "Noise", "Noise", "Noise", "Noise", "Noise",
                      "Noise", "Noise", "Two-indicator state", "Noise", "Two-indicator state",  "Noise", "Noise",  "Noise", "Noise", "Noise", "Three-indicator state",
                      "Noise", "Three-indicator state", "Noise", "Three-indicator state"))
colnames(temp) <- c("cats", "condition", "true_probs", "cats_n")
nobs_cl_theta$cov$cov <- full_join(nobs_cl_theta$cov$cov, temp, by = c("cats", "condition"))
nobs_cl_theta$cov$cov$cats_n <- factor(nobs_cl_theta$cov$cov$cats_n, 
                                       levels = c("One-indicator state", "Two-indicator state", "Three-indicator state", 
                                                  "Two-indicator primary state", "Two-indicator secondary state", "Noise"))
nobs_cl_theta_cov_reduced_plot <- nobs_cl_theta$cov$cov %>%
  group_by(cats_n, length, condition) %>%
  dplyr::summarize(mean = mean(mean_cov), sd = sd(mean_cov)) %>% 
  ggplot(aes(x = length, y = mean, color = condition, group = condition)) +
  facet_wrap(facets = vars(cats_n)) +
  geom_point() + 
  geom_line() +  
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.1) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "black") +
  xlab("Sequence length") + 
  ylab("Coverage") +
  ggtitle("Coverage Emission Probabilities Clear Number of Observations Condition - Reduced") + 
  labs(color = "Clear Number of Observations") +
  geom_hline(yintercept = 0) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 
# Combined plot. 
nobs_cl_theta_cov_combined_plot <- grid.arrange(nobs_cl_theta$cov$cov_complete_plot, 
                                                nobs_cl_theta_cov_reduced_plot, ncol = 2)

  }

### Calculate absolute bias, relative bias, standard deviation, and coverage of the transition probabilities for the clear number 
### of observations condition, with associated plots.

## Specify number of observations. 
nobs <- c(3, 5, 7)
## Specify condition.
cond <- c("varobs", "varobs", "varobs")
condition <- c("3", "5", "7")
cond_string <- c("Clear Number of Observations")
## Calculate the various measures and the associated plots. 
nobs_cl_gamma <- calcR_gamma(true_probs = gamma_true, nobs = nobs, cond = cond, condition = condition, cond_string = cond_string, varobs = T)

### Reduced and combined plot of the absolute bias, relative bias, standard deviation, and coverage of the emission probabilities for the clear number 
### of observations condition, with associated plots.
{
## Absolute bias. 
# Reduced plot. 
nobs_cl_gamma$abs_bias$abs_bias$cats <- car::recode(nobs_cl_gamma$abs_bias$abs_bias$cats, 
                                                      "c('M_S1_to_S1', 'M_S2_to_S2', 'M_S3_to_S3') = 'Self-transition'; 
                                                     else = 'State-to-state transition'")
nobs_cl_gamma$abs_bias$abs_bias$cats <- factor(nobs_cl_gamma$abs_bias$abs_bias$cats, 
                                                 levels = c("Self-transition", "State-to-state transition"))
nobs_cl_gamma_abs_bias_reduced_plot <- nobs_cl_gamma$abs_bias$abs_bias %>%
  group_by(cats, length, condition) %>%
  dplyr::summarize(mean = mean(mean_abs_bias), sd = sd(mean_abs_bias)) %>% 
  ggplot(aes(x = length, y = mean, color = condition, group = condition)) +
  facet_wrap(facets = vars(cats)) +
  geom_point() + 
  geom_line() +  
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.1) +
  xlab("Sequence length") + 
  ylab("Absolute bias") +
  ggtitle("Absolute Bias Transition Probabilities Clear Number of Observations Condition - Reduced") + 
  labs(color = "Clear Number of Observations") +
  geom_hline(yintercept = 0) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 
# Combined plot. 
nobs_cl_gamma_abs_bias_combined_plot <- grid.arrange(nobs_cl_gamma$abs_bias$abs_bias_complete_plot, 
                                                     nobs_cl_gamma_abs_bias_reduced_plot, ncol = 2)

## Relative bias. 
# Reduced plot. 
nobs_cl_gamma$rel_bias$rel_bias$cats <- car::recode(nobs_cl_gamma$rel_bias$rel_bias$cats, 
                                                      "c('M_S1_to_S1', 'M_S2_to_S2', 'M_S3_to_S3') = 'Self-transition'; 
                                                     else = 'State-to-state transition'")
nobs_cl_gamma$rel_bias$rel_bias$cats <- factor(nobs_cl_gamma$rel_bias$rel_bias$cats, 
                                                 levels = c("Self-transition", "State-to-state transition"))
nobs_cl_gamma_rel_bias_reduced_plot <- nobs_cl_gamma$rel_bias$rel_bias %>%
  group_by(cats, length, condition) %>%
  dplyr::summarize(mean = mean(mean_rel_bias), sd = sd(mean_rel_bias)) %>% 
  ggplot(aes(x = length, y = mean, color = condition, group = condition)) +
  facet_wrap(facets = vars(cats)) +
  geom_point() + 
  geom_line() +  
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.1) +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "black") +
  geom_hline(yintercept = -0.1, linetype = "dashed", color = "black") +
  xlab("Sequence length") + 
  ylab("Relative bias") +
  ggtitle("Relative Bias Transition Probabilities Clear Number of Observations Condition - Reduced") + 
  labs(color = "Clear Number of Observations") +
  geom_hline(yintercept = 0) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 
# Combined plot. 
nobs_cl_gamma_rel_bias_combined_plot <- grid.arrange(nobs_cl_gamma$rel_bias$rel_bias_complete_plot, 
                                                     nobs_cl_gamma_rel_bias_reduced_plot, ncol = 2)

## Standard deviation. 
# Reduced plot. 
nobs_cl_gamma$sdev$sdev$cats <- car::recode(nobs_cl_gamma$sdev$sdev$cats, 
                                              "c('sd_S1_to_S1', 'sd_S2_to_S2', 'sd_S3_to_S3') = 'Self-transition'; 
                                               else = 'State-to-state transition'")
nobs_cl_gamma$sdev$sdev$cats <- factor(nobs_cl_gamma$sdev$sdev$cats, 
                                         levels = c("Self-transition", "State-to-state transition"))
nobs_cl_gamma_sdev_reduced_plot <- nobs_cl_gamma$sdev$sdev %>%
  group_by(cats, length, condition) %>%
  dplyr::summarize(mean = mean(mean_sdev), sd = sd(mean_sdev)) %>% 
  ggplot(aes(x = length, y = mean, color = condition, group = condition)) +
  facet_wrap(facets = vars(cats)) +
  geom_point() + 
  geom_line() +  
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.1) +
  xlab("Sequence length") + 
  ylab("Standard deviation") +
  ggtitle("Standard Deviation Transition Probabilities Clear Number of Observations Condition - Reduced") + 
  labs(color = "Clear Number of Observations") +
  geom_hline(yintercept = 0) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 
# Combined plot. 
nobs_cl_gamma_sdev_combined_plot <- grid.arrange(nobs_cl_gamma$sdev$sdev_complete_plot, 
                                                 nobs_cl_gamma_sdev_reduced_plot, ncol = 2)

## Coverage
# Reduced plot. 
nobs_cl_gamma$cov$cov$cats <- car::recode(nobs_cl_gamma$cov$cov$cats, 
                                            "c('M_S1_to_S1', 'M_S2_to_S2', 'M_S3_to_S3') = 'Self-transition'; 
                                                     else = 'State-to-state transition'")
nobs_cl_gamma$cov$cov$cats <- factor(nobs_cl_gamma$cov$cov$cats, 
                                       levels = c("Self-transition", "State-to-state transition"))
nobs_cl_gamma_cov_reduced_plot <- nobs_cl_gamma$cov$cov %>%
  group_by(cats, length, condition) %>%
  dplyr::summarize(mean = mean(mean_cov), sd = sd(mean_cov)) %>% 
  ggplot(aes(x = length, y = mean, color = condition, group = condition)) +
  facet_wrap(facets = vars(cats)) +
  geom_point() + 
  geom_line() +  
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.1) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "black") +
  xlab("Sequence length") + 
  ylab("Coverage") +
  ggtitle("Coverage Transition Probabilities Clear Number of Observations Condition - Reduced") + 
  labs(color = "Clear Number of Observations") +
  geom_hline(yintercept = 0) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 
# Combined plot. 
nobs_cl_gamma_cov_combined_plot <- grid.arrange(nobs_cl_gamma$cov$cov_complete_plot, 
                                                nobs_cl_gamma_cov_reduced_plot, ncol = 2)
  
  }
