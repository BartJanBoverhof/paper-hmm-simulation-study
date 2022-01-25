###########################
##### HMM Paper Code ######
###########################

##############################
### Load packages and data ###
##############################

### Load packages.
if(!require(data.table)) install.packages("data.table")
if(!require(reshape)) install.packages("reshape")
if(!require(ggplot2)) install.packages("ggplot2")
if(!require(ggh4x)) install.packages("ggh4x")
if(!require(plyr)) install.packages("plyr")
if(!require(plyr)) install.packages("tidyverse")

library(data.table)
library(reshape)
library(ggplot2)
library(ggh4x)
library(plyr)
library(tidyverse)

### Set working directory. 
setwd("~/Documents/GitHub/paper-hmm-simulation-study/data")
### Function for assigning datasets to list.
rda2list <- function(file) {
  e <- new.env()
  load(file, envir = e)
  as.list(e)
}

### Loading in and assigning datasets to list. 
files <- list.files(pattern = ".rda$")
results <- Map(rda2list, file.path(files))
names(results) <- tools::file_path_sans_ext(files)

##########################################################
### Theta / Gamma ; Clear / Moderately Clear / Unclear ###
##########################################################

length_strings <- c("250", "500", "1000", "2000", "4000", "8000")

### Theta.

## True emission probabilities.
# Clear.
theta_clear_true <- c(0.92, 0.02, 0.02, 0.02, 0.02,
                      0.02, 0.47, 0.02, 0.47, 0.02,
                      0.02, 0.02, 0.47, 0.02, 0.47)
theta_clear_true = data_frame(theta_clear_true) %>% mutate(c("M_S1_cat1", "M_S1_cat2",  "M_S1_cat3", "M_S1_cat4", "M_S1_cat5",
                                                             "M_S2_cat1", "M_S2_cat2",  "M_S2_cat3", "M_S2_cat4", "M_S2_cat5",
                                                             "M_S3_cat1", "M_S3_cat2",  "M_S3_cat3", "M_S3_cat4", "M_S3_cat5"))
colnames(theta_clear_true) <- c("true", "X2")

# Moderately clear. 
theta_mod_true <- c(0.68, 0.08, 0.08, 0.08, 0.08,
                    0.08, 0.38, 0.08, 0.38, 0.08,
                    0.08, 0.08, 0.38, 0.08, 0.38)

theta_mod_true = data_frame(theta_mod_true) %>% mutate(c("M_S1_cat1", "M_S1_cat2",  "M_S1_cat3", "M_S1_cat4", "M_S1_cat5",
                                                         "M_S2_cat1", "M_S2_cat2",  "M_S2_cat3", "M_S2_cat4", "M_S2_cat5",
                                                         "M_S3_cat1", "M_S3_cat2",  "M_S3_cat3", "M_S3_cat4", "M_S3_cat5"))
colnames(theta_mod_true) <- c("true", "X2")

# Unclear. 
theta_unclear_true <- c(0.44, 0.14, 0.14, 0.14, 0.14,
                        0.14, 0.29, 0.14, 0.29, 0.14,
                        0.14, 0.14, 0.29, 0.14, 0.29)

theta_unclear_true = data_frame(theta_unclear_true) %>% mutate(c("M_S1_cat1", "M_S1_cat2",  "M_S1_cat3", "M_S1_cat4", "M_S1_cat5",
                                                         "M_S2_cat1", "M_S2_cat2",  "M_S2_cat3", "M_S2_cat4", "M_S2_cat5",
                                                         "M_S3_cat1", "M_S3_cat2",  "M_S3_cat3", "M_S3_cat4", "M_S3_cat5"))
colnames(theta_unclear_true) <- c("true", "X2")


for (i in 1: length(length_strings)){
  
  assign(paste0("cov_theta_cl_obs_5_t_", length_strings[i]), 
         sweep(eval(parse(text = paste0("results$sim_HMM_theta_cl_obs_5_t_", length_strings[i], "$out_sim$emiss_mean"))), 
               2, 0)) 
  
  assign(paste0("cov_theta_cl_obs_5_t_", length_strings[i]), 
         melt(eval(parse(text = paste0("cov_theta_cl_obs_5_t_", length_strings[i])))))
  
  assign(paste0("cov_theta_cl_obs_5_t_", length_strings[i]), 
         cbind(eval(parse(text = paste0("cov_theta_cl_obs_5_t_", length_strings[i]))), 
               Length = length_strings[i], Clarity = "Clear", 
               CILow = as.vector(eval(parse(text = paste0("results$sim_HMM_theta_cl_obs_5_t_", length_strings[i], "$out_sim$emiss_low")))),
               CIHigh = as.vector(  eval(parse(text = paste0("results$sim_HMM_theta_cl_obs_5_t_", length_strings[i], "$out_sim$emiss_up"))))))



  assign(paste0("cov_theta_modcl_obs_5_t_", length_strings[i]), 
         sweep(eval(parse(text = paste0("results$sim_HMM_theta_modcl_obs_5_t_", length_strings[i], "$out_sim$emiss_mean"))), 
               2, 0)) 
  
  assign(paste0("cov_theta_modcl_obs_5_t_", length_strings[i]), 
         melt(eval(parse(text = paste0("cov_theta_modcl_obs_5_t_", length_strings[i])))))
  
  assign(paste0("cov_theta_modcl_obs_5_t_", length_strings[i]), 
         cbind(eval(parse(text = paste0("cov_theta_modcl_obs_5_t_", length_strings[i]))), 
               Length = length_strings[i], Clarity = "Moderate", 
               CILow = as.vector(eval(parse(text = paste0("results$sim_HMM_theta_modcl_obs_5_t_", length_strings[i], "$out_sim$emiss_low")))),
               CIHigh = as.vector(eval(parse(text = paste0("results$sim_HMM_theta_modcl_obs_5_t_", length_strings[i], "$out_sim$emiss_up"))))))
  
  
  
  assign(paste0("cov_theta_uncl_obs_5_t_", length_strings[i]), 
         sweep(eval(parse(text = paste0("results$sim_HMM_theta_uncl_obs_5_t_", length_strings[i], "$out_sim$emiss_mean"))), 
               2, 0)) 
  
  assign(paste0("cov_theta_uncl_obs_5_t_", length_strings[i]), 
         melt(eval(parse(text = paste0("cov_theta_uncl_obs_5_t_", length_strings[i])))))
  
  assign(paste0("cov_theta_uncl_obs_5_t_", length_strings[i]), 
         cbind(eval(parse(text = paste0("cov_theta_uncl_obs_5_t_", length_strings[i]))), 
               Length = length_strings[i], Clarity = "Unclear", 
               CILow = as.vector(eval(parse(text = paste0("results$sim_HMM_theta_uncl_obs_5_t_", length_strings[i], "$out_sim$emiss_low")))),
               CIHigh = as.vector(eval(parse(text = paste0("results$sim_HMM_theta_uncl_obs_5_t_", length_strings[i], "$out_sim$emiss_up"))))))
  
  
}

## Clear.
cl_obs_5_theta <- rbind(data.frame(cov_theta_cl_obs_5_t_250), data.frame(cov_theta_cl_obs_5_t_500),  
                        data.frame(cov_theta_cl_obs_5_t_1000), data.frame(cov_theta_cl_obs_5_t_2000), 
                        data.frame(cov_theta_cl_obs_5_t_4000), data.frame(cov_theta_cl_obs_5_t_8000))

cl_obs_5_theta <- full_join(cl_obs_5_theta, theta_clear_true, by = "X2")

cl_obs_5_theta <- cl_obs_5_theta %>% mutate(Coverage = ifelse((value > CILow & true < CIHigh),1,0))
colnames(cl_obs_5_theta) <- c("Id", "S_to_obs", "Estimate", "Length", "Clarity","CiLow", "CiHigh","true","Coverage")

cl_obs_5_theta$Length <- factor(cl_obs_5_theta$Length, levels = c("250", "500", "1000", "2000", "4000", "8000"))
cl_obs_5_theta <- aggregate(Coverage ~ Length + S_to_obs + Clarity, cl_obs_5_theta, mean)


## Moderately clear. 
modcl_obs_5_theta <- rbind(data.frame(cov_theta_modcl_obs_5_t_250), data.frame(cov_theta_modcl_obs_5_t_500),  
                        data.frame(cov_theta_modcl_obs_5_t_1000), data.frame(cov_theta_modcl_obs_5_t_2000), 
                        data.frame(cov_theta_modcl_obs_5_t_4000), data.frame(cov_theta_modcl_obs_5_t_8000))

modcl_obs_5_theta <- full_join(modcl_obs_5_theta, theta_mod_true, by = "X2")

modcl_obs_5_theta <- modcl_obs_5_theta %>% mutate(Coverage = ifelse((value > CILow & true < CIHigh),1,0))
colnames(modcl_obs_5_theta) <- c("Id", "S_to_obs", "Estimate", "Length", "Clarity","CiLow", "CiHigh","true" ,"Coverage")

modcl_obs_5_theta$Length <- factor(modcl_obs_5_theta$Length, levels = c("250", "500", "1000", "2000", "4000", "8000"))
modcl_obs_5_theta <- aggregate(Coverage ~ Length + S_to_obs + Clarity, modcl_obs_5_theta, mean)


## Unclear. 
uncl_obs_5_theta <- rbind(data.frame(cov_theta_uncl_obs_5_t_250), data.frame(cov_theta_uncl_obs_5_t_500),  
                        data.frame(cov_theta_uncl_obs_5_t_1000), data.frame(cov_theta_uncl_obs_5_t_2000), 
                        data.frame(cov_theta_uncl_obs_5_t_4000), data.frame(cov_theta_uncl_obs_5_t_8000))

uncl_obs_5_theta <- full_join(uncl_obs_5_theta, theta_unclear_true, by = "X2")

uncl_obs_5_theta <- uncl_obs_5_theta %>% mutate(Coverage = ifelse((true > CILow & true < CIHigh),1,0))
colnames(uncl_obs_5_theta) <- c("Id", "S_to_obs", "Estimate", "Length", "Clarity","CiLow", "CiHigh","true","Coverage")

uncl_obs_5_theta$Length <- factor(uncl_obs_5_theta$Length, levels = c("250", "500", "1000", "2000", "4000", "8000"))
uncl_obs_5_theta <- aggregate(Coverage ~ Length + S_to_obs + Clarity, uncl_obs_5_theta, mean)



## Plot of bias emission probabilities by clarity and sequence length.
obs_5_theta <- rbind(data.frame(cl_obs_5_theta), data.frame(modcl_obs_5_theta),  
                     data.frame(uncl_obs_5_theta))
obs_5_theta$S_to_obs <- mapvalues(obs_5_theta$S_to_obs, 
                                  from = c("M_S1_cat1", "M_S1_cat2", "M_S1_cat3", "M_S1_cat4", "M_S1_cat5", 
                                           "M_S2_cat1", "M_S2_cat2", "M_S2_cat3", "M_S2_cat4", "M_S2_cat5", 
                                           "M_S3_cat1", "M_S3_cat2", "M_S3_cat3", "M_S3_cat4", "M_S3_cat5"), 
                                  to = c(rep(paste0("theta[", 1, "][", 1:5, "]"), each = 1), 
                                         rep(paste0("theta[", 2, "][", 1:5, "]"), each = 1), 
                                         rep(paste0("theta[", 3, "][", 1:5, "]"), each = 1)))
obs_5_theta$Clarity <- factor(obs_5_theta$Clarity, levels = c("Clear", "Moderate", "Unclear"))
ggplot(obs_5_theta, aes(x = Length, y = Coverage, color = Clarity, group = Clarity)) +
  facet_wrap(facets = vars(S_to_obs), nrow = 3, ncol = 5, labeller = label_parsed) +
  geom_point() + 
  geom_line() +  
  xlab("Sequence length") + 
  ylab("Bias") +
  geom_hline(yintercept = 0) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  ylim(0.5, 1) 



### Gamma.

## True transition probabilities.
gamma_true <- c(0.80, 0.10, 0.10, 
                0.10, 0.80, 0.10, 
                0.10, 0.10, 0.80)

gamma_true = data_frame(gamma_true) %>% mutate(c("M_S1_to_S1", "M_S1_to_S2", "M_S1_to_S3",
                                                 "M_S2_to_S1", "M_S2_to_S2", "M_S2_to_S3", 
                                                 "M_S3_to_S1", "M_S3_to_S2", "M_S3_to_S3")) 
                                                 
colnames(gamma_true) <- c("true", "X2")

for (i in 1: length(length_strings)){
  
  assign(paste0("cov_gamma_cl_obs_5_t_", length_strings[i]), 
         sweep(eval(parse(text = paste0("results$sim_HMM_theta_cl_obs_5_t_", length_strings[i], "$out_sim$gamma_mean"))), 
               2, gamma_true)) 
  
  assign(paste0("cov_gamma_cl_obs_5_t_", length_strings[i]), 
         melt(eval(parse(text = paste0("cov_gamma_cl_obs_5_t_", length_strings[i])))))
  
  assign(paste0("cov_gamma_cl_obs_5_t_", length_strings[i]), 
         cbind(eval(parse(text = paste0("cov_gamma_cl_obs_5_t_", length_strings[i]))), 
               Length = length_strings[i], Clarity = "Clear",
               CILow = as.vector(eval(parse(text = paste0("results$sim_HMM_theta_cl_obs_5_t_", length_strings[i], "$out_sim$gamma_low")))),
               CIHigh = as.vector(eval(parse(text = paste0("results$sim_HMM_theta_cl_obs_5_t_", length_strings[i], "$out_sim$gamma_up"))))))
  
  
  assign(paste0("cov_gamma_modcl_obs_5_t_", length_strings[i]), 
         sweep(eval(parse(text = paste0("results$sim_HMM_theta_modcl_obs_5_t_", length_strings[i], "$out_sim$gamma_mean"))), 
               2, gamma_true)) 
  
  assign(paste0("cov_gamma_modcl_obs_5_t_", length_strings[i]), 
         melt(eval(parse(text = paste0("cov_gamma_modcl_obs_5_t_", length_strings[i])))))
  
  assign(paste0("cov_gamma_modcl_obs_5_t_", length_strings[i]), 
         cbind(eval(parse(text = paste0("cov_gamma_modcl_obs_5_t_", length_strings[i]))), 
               Length = length_strings[i], Clarity = "Moderate",
               CILow = as.vector(eval(parse(text = paste0("results$sim_HMM_theta_modcl_obs_5_t_", length_strings[i], "$out_sim$gamma_low")))),
               CIHigh = as.vector(eval(parse(text = paste0("results$sim_HMM_theta_modcl_obs_5_t_", length_strings[i], "$out_sim$gamma_up"))))))
  
  
  
  
  assign(paste0("cov_gamma_uncl_obs_5_t_", length_strings[i]), 
         sweep(eval(parse(text = paste0("results$sim_HMM_theta_uncl_obs_5_t_", length_strings[i], "$out_sim$gamma_mean"))), 
               2, gamma_true)) 
  
  assign(paste0("cov_gamma_uncl_obs_5_t_", length_strings[i]), 
         melt(eval(parse(text = paste0("cov_gamma_uncl_obs_5_t_", length_strings[i])))))
  
  assign(paste0("cov_gamma_uncl_obs_5_t_", length_strings[i]), 
         cbind(eval(parse(text = paste0("cov_gamma_uncl_obs_5_t_", length_strings[i]))), 
               Length = length_strings[i], Clarity = "Unclear",
               CILow = as.vector(eval(parse(text = paste0("results$sim_HMM_theta_cl_obs_5_t_", length_strings[i], "$out_sim$gamma_low")))),
               CIHigh = as.vector(eval(parse(text = paste0("results$sim_HMM_theta_cl_obs_5_t_", length_strings[i], "$out_sim$gamma_up"))))))
  
}

# Clear. 
cl_obs_5_gamma <- rbind(data.frame(cov_gamma_cl_obs_5_t_250), data.frame(cov_gamma_cl_obs_5_t_500),  
                        data.frame(cov_gamma_cl_obs_5_t_1000), data.frame(cov_gamma_cl_obs_5_t_2000), 
                        data.frame(cov_gamma_cl_obs_5_t_4000), data.frame(cov_gamma_cl_obs_5_t_8000))
cl_obs_5_gamma <- full_join(cl_obs_5_gamma, gamma_true, by = "X2")

cl_obs_5_gamma <- cl_obs_5_gamma %>% mutate(Coverage = ifelse((true > CILow & true < CIHigh),1,0))
colnames(cl_obs_5_gamma) <- c("Id", "S_to_s", "Estimate", "Length", "Clarity","CiLow", "CiHigh","true","Coverage")

cl_obs_5_gamma$Length <- factor(cl_obs_5_gamma$Length, levels = c("250", "500", "1000", "2000", "4000", "8000"))
cl_obs_5_gamma <- aggregate(Coverage ~ Length + S_to_s + Clarity, cl_obs_5_gamma, mean)

# Moderately clear. 
modcl_obs_5_gamma <- rbind(data.frame(cov_gamma_modcl_obs_5_t_250), data.frame(cov_gamma_modcl_obs_5_t_500),  
                        data.frame(cov_gamma_modcl_obs_5_t_1000), data.frame(cov_gamma_modcl_obs_5_t_2000), 
                        data.frame(cov_gamma_modcl_obs_5_t_4000), data.frame(cov_gamma_modcl_obs_5_t_8000))
modcl_obs_5_gamma <- full_join(modcl_obs_5_gamma, gamma_true, by = "X2")

modcl_obs_5_gamma <- modcl_obs_5_gamma %>% mutate(Coverage = ifelse((true > CILow & true < CIHigh),1,0))
colnames(modcl_obs_5_gamma) <- c("Id", "S_to_s", "Estimate", "Length", "Clarity","CiLow", "CiHigh","true","Coverage")

modcl_obs_5_gamma$Length <- factor(modcl_obs_5_gamma$Length, levels = c("250", "500", "1000", "2000", "4000", "8000"))
modcl_obs_5_gamma <- aggregate(Coverage ~ Length + S_to_s + Clarity, modcl_obs_5_gamma, mean)

# Unclear. 
uncl_obs_5_gamma <- rbind(data.frame(cov_gamma_uncl_obs_5_t_250), data.frame(cov_gamma_uncl_obs_5_t_500),  
                           data.frame(cov_gamma_uncl_obs_5_t_1000), data.frame(cov_gamma_uncl_obs_5_t_2000), 
                           data.frame(cov_gamma_uncl_obs_5_t_4000), data.frame(cov_gamma_uncl_obs_5_t_8000))
uncl_obs_5_gamma <- full_join(uncl_obs_5_gamma, gamma_true, by = "X2")

uncl_obs_5_gamma <- uncl_obs_5_gamma %>% mutate(Coverage = ifelse((true > CILow & true < CIHigh),1,0))
colnames(uncl_obs_5_gamma) <- c("Id", "S_to_s", "Estimate", "Length", "Clarity","CiLow", "CiHigh","true","Coverage")

uncl_obs_5_gamma$Length <- factor(uncl_obs_5_gamma$Length, levels = c("250", "500", "1000", "2000", "4000", "8000"))
uncl_obs_5_gamma <- aggregate(Coverage ~ Length + S_to_s + Clarity, uncl_obs_5_gamma, mean)


# Plot of bias transition probabilities by clarity and sequence length.
obs_5_gamma <- rbind(data.frame(cl_obs_5_gamma), data.frame(modcl_obs_5_gamma),  
                     data.frame(uncl_obs_5_gamma))
obs_5_gamma$S_to_s <- mapvalues(obs_5_gamma$S_to_s, 
                                from = c("M_S1_to_S1", "M_S1_to_S2", "M_S1_to_S3", 
                                         "M_S2_to_S1", "M_S2_to_S2", "M_S2_to_S3",
                                         "M_S3_to_S1", "M_S3_to_S2", "M_S3_to_S3"), 
                                to = c(rep(paste0("gamma[", 1, "][", 1:3, "]"), each = 1), 
                                       rep(paste0("gamma[", 2, "][", 1:3, "]"), each = 1), 
                                       rep(paste0("gamma[", 3, "][", 1:3, "]"), each = 1)))
obs_5_gamma$Clarity <- factor(obs_5_gamma$Clarity, levels = c("Clear", "Moderate", "Unclear"))
ggplot(obs_5_gamma, aes(x = Length, y = Coverage, color = Clarity, group = Clarity)) +
  facet_wrap(facets = vars(S_to_s), nrow = 3, ncol = 3, labeller = label_parsed) +
  geom_point() + geom_line() +  
  ggtitle("") +
  xlab("Sequence length") + 
  ylab("Bias") +
  geom_hline(yintercept = 0) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  ylim(0.5, 1)


#########################
### Overlap condition ###
#########################

### Theta.

## True emission probabilities.

# No overlap.
no_overlap_true <- c(0.84, 0.04, 0.04, 0.04, 0.04,
                     0.04, 0.44, 0.44, 0.04, 0.04,
                     0.04, 0.04, 0.04, 0.44, 0.44)  

no_overlap_true = data_frame(no_overlap_true) %>% mutate(c("M_S1_cat1", "M_S1_cat2",  "M_S1_cat3", "M_S1_cat4", "M_S1_cat5",
                                                             "M_S2_cat1", "M_S2_cat2",  "M_S2_cat3", "M_S2_cat4", "M_S2_cat5",
                                                             "M_S3_cat1", "M_S3_cat2",  "M_S3_cat3", "M_S3_cat4", "M_S3_cat5"))
colnames(no_overlap_true) <- c("true", "X2")
# Moderate overlap.
mod_overlap_true <- c(0.59, 0.29, 0.04, 0.04, 0.04,
                      0.04, 0.29, 0.59, 0.04, 0.04,
                      0.04, 0.04, 0.20, 0.36, 0.36)

mod_overlap_true = data_frame(mod_overlap_true) %>% mutate(c("M_S1_cat1", "M_S1_cat2",  "M_S1_cat3", "M_S1_cat4", "M_S1_cat5",
                                                           "M_S2_cat1", "M_S2_cat2",  "M_S2_cat3", "M_S2_cat4", "M_S2_cat5",
                                                           "M_S3_cat1", "M_S3_cat2",  "M_S3_cat3", "M_S3_cat4", "M_S3_cat5"))
colnames(mod_overlap_true) <- c("true", "X2")

# Much overlap.
much_overlap_true <- c(0.44, 0.44, 0.04, 0.04, 0.04,
                       0.04, 0.44, 0.44, 0.04, 0.04,
                       0.04, 0.04, 0.30, 0.31, 0.31)

much_overlap_true = data_frame(much_overlap_true) %>% mutate(c("M_S1_cat1", "M_S1_cat2",  "M_S1_cat3", "M_S1_cat4", "M_S1_cat5",
                                                           "M_S2_cat1", "M_S2_cat2",  "M_S2_cat3", "M_S2_cat4", "M_S2_cat5",
                                                           "M_S3_cat1", "M_S3_cat2",  "M_S3_cat3", "M_S3_cat4", "M_S3_cat5"))
colnames(much_overlap_true) <- c("true", "X2")


for (i in 1: length(length_strings)){
  
  assign(paste0("cov_theta_nooverl_obs_5_t_", length_strings[i]), 
         sweep(eval(parse(text = paste0("results$sim_HMM_theta_nooverl_obs_5_t_", length_strings[i], "$out_sim$emiss_mean"))), 
               2, 0)) 
  
  assign(paste0("cov_theta_nooverl_obs_5_t_", length_strings[i]), 
         melt(eval(parse(text = paste0("cov_theta_nooverl_obs_5_t_", length_strings[i])))))
  
  assign(paste0("cov_theta_nooverl_obs_5_t_", length_strings[i]), 
         cbind(eval(parse(text = paste0("cov_theta_nooverl_obs_5_t_", length_strings[i]))), 
               Length = length_strings[i], Overlap = "None",
               CILow = as.vector(eval(parse(text = paste0("results$sim_HMM_theta_nooverl_obs_5_t_", length_strings[i], "$out_sim$emiss_low")))),
               CIHigh = as.vector(eval(parse(text = paste0("results$sim_HMM_theta_nooverl_obs_5_t_", length_strings[i], "$out_sim$emiss_up"))))))
  
  
  assign(paste0("cov_theta_modoverl_obs_5_t_", length_strings[i]), 
         sweep(eval(parse(text = paste0("results$sim_HMM_theta_modoverl_obs_5_t_", length_strings[i], "$out_sim$emiss_mean"))), 
               2, 0)) 
  
  assign(paste0("cov_theta_modoverl_obs_5_t_", length_strings[i]), 
         melt(eval(parse(text = paste0("cov_theta_modoverl_obs_5_t_", length_strings[i])))))
  
  assign(paste0("cov_theta_modoverl_obs_5_t_", length_strings[i]), 
         cbind(eval(parse(text = paste0("cov_theta_modoverl_obs_5_t_", length_strings[i]))), 
               Length = length_strings[i], Overlap = "Moderate",
               CILow = as.vector(eval(parse(text = paste0("results$sim_HMM_theta_modoverl_obs_5_t_", length_strings[i], "$out_sim$emiss_low")))),
               CIHigh = as.vector(eval(parse(text = paste0("results$sim_HMM_theta_modoverl_obs_5_t_", length_strings[i], "$out_sim$emiss_up"))))))
  
  
  assign(paste0("cov_theta_muchoverl_obs_5_t_", length_strings[i]), 
         sweep(eval(parse(text = paste0("results$sim_HMM_theta_muchoverl_obs_5_t_", length_strings[i], "$out_sim$emiss_mean"))), 
               2, 0)) 
  
  assign(paste0("cov_theta_muchoverl_obs_5_t_", length_strings[i]), 
         melt(eval(parse(text = paste0("cov_theta_muchoverl_obs_5_t_", length_strings[i])))))
  
  assign(paste0("cov_theta_muchoverl_obs_5_t_", length_strings[i]), 
         cbind(eval(parse(text = paste0("cov_theta_muchoverl_obs_5_t_", length_strings[i]))), 
               Length = length_strings[i], Overlap = "Much",
               CILow = as.vector(eval(parse(text = paste0("results$sim_HMM_theta_muchoverl_obs_5_t_", length_strings[i], "$out_sim$emiss_low")))),
               CIHigh = as.vector(eval(parse(text = paste0("results$sim_HMM_theta_muchoverl_obs_5_t_", length_strings[i], "$out_sim$emiss_up"))))))
  
}

# No overlap.
nooverl_obs_5_theta <- rbind(data.frame(cov_theta_nooverl_obs_5_t_250), data.frame(cov_theta_nooverl_obs_5_t_500),  
                             data.frame(cov_theta_nooverl_obs_5_t_1000), data.frame(cov_theta_nooverl_obs_5_t_2000), 
                             data.frame(cov_theta_nooverl_obs_5_t_4000), data.frame(cov_theta_nooverl_obs_5_t_8000))

nooverl_obs_5_theta <- full_join(nooverl_obs_5_theta, no_overlap_true, by = "X2")

nooverl_obs_5_theta <- nooverl_obs_5_theta %>% mutate(Coverage = ifelse((value > CILow & true < CIHigh),1,0))
colnames(nooverl_obs_5_theta) <- c("Id", "S_to_obs", "Estimate", "Length", "Overlap","CiLow", "CiHigh","true","Coverage")

nooverl_obs_5_theta$Length <- factor(nooverl_obs_5_theta$Length, levels = c("250", "500", "1000", "2000", "4000", "8000"))
nooverl_obs_5_theta <- aggregate(Coverage ~ Length + S_to_obs + Overlap, nooverl_obs_5_theta, mean)

# Moderate overlap. 
modoverl_obs_5_theta <- rbind(data.frame(cov_theta_modoverl_obs_5_t_250), data.frame(cov_theta_modoverl_obs_5_t_500),  
                              data.frame(cov_theta_modoverl_obs_5_t_1000), data.frame(cov_theta_modoverl_obs_5_t_2000), 
                              data.frame(cov_theta_modoverl_obs_5_t_4000), data.frame(cov_theta_modoverl_obs_5_t_8000))

modoverl_obs_5_theta <- full_join(modoverl_obs_5_theta, mod_overlap_true, by = "X2")

modoverl_obs_5_theta <- modoverl_obs_5_theta %>% mutate(Coverage = ifelse((value > CILow & true < CIHigh),1,0))
colnames(modoverl_obs_5_theta) <- c("Id", "S_to_obs", "Estimate", "Length", "Overlap","CiLow", "CiHigh","true","Coverage")

modoverl_obs_5_theta$Length <- factor(modoverl_obs_5_theta$Length, levels = c("250", "500", "1000", "2000", "4000", "8000"))
modoverl_obs_5_theta <- aggregate(Coverage ~ Length + S_to_obs + Overlap, modoverl_obs_5_theta, mean)

# Much overlap. 
muchoverl_obs_5_theta <- rbind(data.frame(cov_theta_muchoverl_obs_5_t_250), data.frame(cov_theta_muchoverl_obs_5_t_500),  
                               data.frame(cov_theta_muchoverl_obs_5_t_1000), data.frame(cov_theta_muchoverl_obs_5_t_2000), 
                               data.frame(cov_theta_muchoverl_obs_5_t_4000), data.frame(cov_theta_muchoverl_obs_5_t_8000))

muchoverl_obs_5_theta <- full_join(muchoverl_obs_5_theta, much_overlap_true, by = "X2")

muchoverl_obs_5_theta <- muchoverl_obs_5_theta %>% mutate(Coverage = ifelse((value > CILow & true < CIHigh),1,0))
colnames(muchoverl_obs_5_theta) <- c("Id", "S_to_obs", "Estimate", "Length", "Overlap","CiLow", "CiHigh","true","Coverage")

muchoverl_obs_5_theta$Length <- factor(muchoverl_obs_5_theta$Length, levels = c("250", "500", "1000", "2000", "4000", "8000"))
muchoverl_obs_5_theta <- aggregate(Coverage ~ Length + S_to_obs + Overlap, muchoverl_obs_5_theta, mean)

# Plot of bias emission probabilities by overlap and sequence length
overl_5_theta <- rbind(data.frame(nooverl_obs_5_theta), data.frame(modoverl_obs_5_theta),
                       data.frame(muchoverl_obs_5_theta))
overl_5_theta$S_to_obs <- mapvalues(overl_5_theta$S_to_obs, 
                                    from = c("M_S1_cat1", "M_S1_cat2", "M_S1_cat3", "M_S1_cat4", "M_S1_cat5", 
                                             "M_S2_cat1", "M_S2_cat2", "M_S2_cat3", "M_S2_cat4", "M_S2_cat5", 
                                             "M_S3_cat1", "M_S3_cat2", "M_S3_cat3", "M_S3_cat4", "M_S3_cat5"), 
                                    to = c(rep(paste0("theta[", 1, "][", 1:5, "]"), each = 1), 
                                           rep(paste0("theta[", 2, "][", 1:5, "]"), each = 1), 
                                           rep(paste0("theta[", 3, "][", 1:5, "]"), each = 1)))
overl_5_theta$Overlap <- factor(overl_5_theta$Overlap, levels = c("None", "Moderate", "Much"))
ggplot(overl_5_theta, aes(x = Length, y = Coverage, color = Overlap, group = Overlap)) +
  facet_wrap(facets = vars(S_to_obs), nrow = 3, ncol = 5, labeller = label_parsed) +
  geom_point() + 
  geom_line() +
  xlab("Sequence length") + 
  ylab("Bias") +
  geom_hline(yintercept = 0) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  ylim(0, 1) 


### Gamma.
for (i in 1: length(length_strings)){
  
  assign(paste0("cov_gamma_nooverl_obs_5_t_", length_strings[i]), 
         sweep(eval(parse(text = paste0("results$sim_HMM_theta_nooverl_obs_5_t_", length_strings[i], "$out_sim$gamma_mean"))), 
               2, 0)) 
  
  assign(paste0("cov_gamma_nooverl_obs_5_t_", length_strings[i]), 
         melt(eval(parse(text = paste0("cov_gamma_nooverl_obs_5_t_", length_strings[i])))))
  
  assign(paste0("cov_gamma_nooverl_obs_5_t_", length_strings[i]), 
         cbind(eval(parse(text = paste0("cov_gamma_nooverl_obs_5_t_", length_strings[i]))), 
               Length = length_strings[i], Overlap = "None",
               CILow = as.vector(eval(parse(text = paste0("results$sim_HMM_theta_nooverl_obs_5_t_", length_strings[i], "$out_sim$gamma_low")))),
               CIHigh = as.vector(eval(parse(text = paste0("results$sim_HMM_theta_nooverl_obs_5_t_", length_strings[i], "$out_sim$gamma_up"))))))
  
  
  assign(paste0("cov_gamma_modoverl_obs_5_t_", length_strings[i]), 
         sweep(eval(parse(text = paste0("results$sim_HMM_theta_modoverl_obs_5_t_", length_strings[i], "$out_sim$gamma_mean"))), 
               2, 0)) 
  
  assign(paste0("cov_gamma_modoverl_obs_5_t_", length_strings[i]), 
         melt(eval(parse(text = paste0("cov_gamma_modoverl_obs_5_t_", length_strings[i])))))
  
  assign(paste0("cov_gamma_modoverl_obs_5_t_", length_strings[i]), 
         cbind(eval(parse(text = paste0("cov_gamma_modoverl_obs_5_t_", length_strings[i]))), 
               Length = length_strings[i], Overlap = "Moderate",
               CILow = as.vector(eval(parse(text = paste0("results$sim_HMM_theta_modoverl_obs_5_t_", length_strings[i], "$out_sim$gamma_low")))),
               CIHigh = as.vector(eval(parse(text = paste0("results$sim_HMM_theta_modoverl_obs_5_t_", length_strings[i], "$out_sim$gamma_up"))))))
  
  assign(paste0("cov_gamma_muchoverl_obs_5_t_", length_strings[i]), 
         sweep(eval(parse(text = paste0("results$sim_HMM_theta_muchoverl_obs_5_t_", length_strings[i], "$out_sim$gamma_mean"))), 
               2, 0)) 
  
  assign(paste0("cov_gamma_muchoverl_obs_5_t_", length_strings[i]), 
         melt(eval(parse(text = paste0("cov_gamma_muchoverl_obs_5_t_", length_strings[i])))))
  
  assign(paste0("cov_gamma_muchoverl_obs_5_t_", length_strings[i]), 
         cbind(eval(parse(text = paste0("cov_gamma_muchoverl_obs_5_t_", length_strings[i]))), 
               Length = length_strings[i], Overlap = "Much",
               CILow = as.vector(eval(parse(text = paste0("results$sim_HMM_theta_muchoverl_obs_5_t_", length_strings[i], "$out_sim$gamma_low")))),
               CIHigh = as.vector(eval(parse(text = paste0("results$sim_HMM_theta_muchoverl_obs_5_t_", length_strings[i], "$out_sim$gamma_up"))))))
  
}

# No overlap.
nooverl_obs_5_gamma <- rbind(data.frame(cov_gamma_nooverl_obs_5_t_250), data.frame(cov_gamma_nooverl_obs_5_t_500),  
                             data.frame(cov_gamma_nooverl_obs_5_t_1000), data.frame(cov_gamma_nooverl_obs_5_t_2000), 
                             data.frame(cov_gamma_nooverl_obs_5_t_4000), data.frame(cov_gamma_nooverl_obs_5_t_8000))

nooverl_obs_5_gamma <- full_join(nooverl_obs_5_gamma, gamma_true, by = "X2")
nooverl_obs_5_gamma <- nooverl_obs_5_gamma %>% mutate(Coverage = ifelse((true > CILow & true < CIHigh),1,0))
colnames(nooverl_obs_5_gamma) <- c("Id", "S_to_s", "Estimate", "Length", "Overlap","CiLow", "CiHigh","true","Coverage")

nooverl_obs_5_gamma$Length <- factor(nooverl_obs_5_gamma$Length, levels = c("250", "500", "1000", "2000", "4000", "8000"))
nooverl_obs_5_gamma <- aggregate(Coverage ~ Length + S_to_s + Overlap, nooverl_obs_5_gamma, mean)

# Moderate overlap. 
modoverl_obs_5_gamma <- rbind(data.frame(cov_gamma_modoverl_obs_5_t_250), data.frame(cov_gamma_modoverl_obs_5_t_500),  
                              data.frame(cov_gamma_modoverl_obs_5_t_1000), data.frame(cov_gamma_modoverl_obs_5_t_2000), 
                              data.frame(cov_gamma_modoverl_obs_5_t_4000), data.frame(cov_gamma_modoverl_obs_5_t_8000))

modoverl_obs_5_gamma <- full_join(modoverl_obs_5_gamma, gamma_true, by = "X2")
modoverl_obs_5_gamma <- modoverl_obs_5_gamma %>% mutate(Coverage = ifelse((true > CILow & true < CIHigh),1,0))
colnames(modoverl_obs_5_gamma) <- c("Id", "S_to_s", "Estimate", "Length", "Overlap","CiLow", "CiHigh","true","Coverage")

modoverl_obs_5_gamma$Length <- factor(modoverl_obs_5_gamma$Length, levels = c("250", "500", "1000", "2000", "4000", "8000"))
modoverl_obs_5_gamma <- aggregate(Coverage ~ Length + S_to_s + Overlap, modoverl_obs_5_gamma, mean)

# Much overlap. 
muchoverl_obs_5_gamma <- rbind(data.frame(cov_gamma_muchoverl_obs_5_t_250), data.frame(cov_gamma_muchoverl_obs_5_t_500),  
                               data.frame(cov_gamma_muchoverl_obs_5_t_1000), data.frame(cov_gamma_muchoverl_obs_5_t_2000), 
                               data.frame(cov_gamma_muchoverl_obs_5_t_4000), data.frame(cov_gamma_muchoverl_obs_5_t_8000))

muchoverl_obs_5_gamma <- full_join(muchoverl_obs_5_gamma, gamma_true, by = "X2")
muchoverl_obs_5_gamma <- muchoverl_obs_5_gamma %>% mutate(Coverage = ifelse((true > CILow & true < CIHigh),1,0))
colnames(muchoverl_obs_5_gamma) <- c("Id", "S_to_s", "Estimate", "Length", "Overlap","CiLow", "CiHigh","true","Coverage")

muchoverl_obs_5_gamma$Length <- factor(muchoverl_obs_5_gamma$Length, levels = c("250", "500", "1000", "2000", "4000", "8000"))
muchoverl_obs_5_gamma <- aggregate(Coverage ~ Length + S_to_s + Overlap, muchoverl_obs_5_gamma, mean)

# Plot of bias transition probabilities by overlap and sequence length
overl_5_gamma <- rbind(data.frame(nooverl_obs_5_gamma), data.frame(modoverl_obs_5_gamma),
                       data.frame(muchoverl_obs_5_gamma))
overl_5_gamma$S_to_s <- mapvalues(overl_5_gamma$S_to_s, 
                                  from = c("M_S1_to_S1", "M_S1_to_S2", "M_S1_to_S3", 
                                           "M_S2_to_S1", "M_S2_to_S2", "M_S2_to_S3",
                                           "M_S3_to_S1", "M_S3_to_S2", "M_S3_to_S3"), 
                                  to = c(rep(paste0("gamma[", 1, "][", 1:3, "]"), each = 1), 
                                         rep(paste0("gamma[", 2, "][", 1:3, "]"), each = 1), 
                                         rep(paste0("gamma[", 3, "][", 1:3, "]"), each = 1)))
overl_5_gamma$Overlap <- factor(overl_5_gamma$Overlap, levels = c("None", "Moderate", "Much"))
ggplot(overl_5_gamma, aes(x = Length, y = Coverage, color = Overlap, group = Overlap)) +
  facet_wrap(facets = vars(S_to_s), nrow = 3, ncol = 3, labeller = label_parsed) +
  geom_point() + 
  geom_line() +  
  xlab("Sequence length") + 
  ylab("Bias") +
  geom_hline(yintercept = 0) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  ylim(0.8, 1)




##############################
### Number of observations ###
##############################

### Theta.

## True emission probabilities.

# Three observations.
true_theta_three_obs <- c(0.52, 0.24, 0.24,
                          0.24, 0.52, 0.24,
                          0.24, 0.24, 0.52)

true_theta_three_obs = data_frame(true_theta_three_obs) %>% mutate(c("M_S1_cat1", "M_S1_cat2","M_S1_cat3",
                                                                     "M_S2_cat1", "M_S2_cat2","M_S2_cat3",
                                                                     "M_S3_cat1", "M_S3_cat2","M_S3_cat3"))
colnames(true_theta_three_obs) <- c("true", "X2")

# Five observations.
true_theta_five_obs <- c(0.44, 0.14, 0.14, 0.14, 0.14,
                         0.14, 0.29, 0.14, 0.29, 0.14,
                         0.14, 0.14, 0.29, 0.14, 0.29)

true_theta_five_obs = data_frame(true_theta_five_obs) %>% mutate(c("M_S1_cat1", "M_S1_cat2",  "M_S1_cat3", "M_S1_cat4", "M_S1_cat5",
                                                             "M_S2_cat1", "M_S2_cat2",  "M_S2_cat3", "M_S2_cat4", "M_S2_cat5",
                                                             "M_S3_cat1", "M_S3_cat2",  "M_S3_cat3", "M_S3_cat4", "M_S3_cat5"))
colnames(true_theta_five_obs) <- c("true", "X2")

# Seven observations.
true_theta_seven_obs <- c(0.40, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10,
                          0.10, 0.25, 0.10, 0.25, 0.10, 0.10, 0.10,
                          0.08, 0.08, 0.22, 0.09, 0.22, 0.09, 0.22)

true_theta_seven_obs = data_frame(true_theta_seven_obs) %>% mutate(c("M_S1_cat1", "M_S1_cat2",  "M_S1_cat3", "M_S1_cat4", "M_S1_cat5", "M_S1_cat6", "M_S1_cat7",
                                                                   "M_S2_cat1", "M_S2_cat2",  "M_S2_cat3", "M_S2_cat4", "M_S2_cat5", "M_S2_cat6", "M_S2_cat7",
                                                                   "M_S3_cat1", "M_S3_cat2",  "M_S3_cat3", "M_S3_cat4", "M_S3_cat5", "M_S3_cat6", "M_S3_cat7"))
colnames(true_theta_seven_obs) <- c("true", "X2")


for (i in 1: length(length_strings)){
  
  assign(paste0("cov_theta_uncl_obs_3_t_", length_strings[i]), 
         sweep(eval(parse(text = paste0("results$sim_HMM_theta_uncl_obs_3_t_", length_strings[i], "$out_sim$emiss_mean"))), 
               2, 0)) 
  
  assign(paste0("cov_theta_uncl_obs_3_t_", length_strings[i]), 
         melt(eval(parse(text = paste0("cov_theta_uncl_obs_3_t_", length_strings[i])))))
  
  assign(paste0("cov_theta_uncl_obs_3_t_", length_strings[i]), 
         cbind(eval(parse(text = paste0("cov_theta_uncl_obs_3_t_", length_strings[i]))), 
               Length = length_strings[i], Nobs = "3",
               CILow = as.vector(eval(parse(text = paste0("results$sim_HMM_theta_uncl_obs_3_t_", length_strings[i], "$out_sim$emiss_low")))),
               CIHigh = as.vector(eval(parse(text = paste0("results$sim_HMM_theta_uncl_obs_3_t_", length_strings[i], "$out_sim$emiss_up"))))))
  
  assign(paste0("cov_theta_uncl_obs_5_t_", length_strings[i]), 
         sweep(eval(parse(text = paste0("results$sim_HMM_theta_uncl_obs_5_t_", length_strings[i], "$out_sim$emiss_mean"))), 
               2, 0)) 
  
  assign(paste0("cov_theta_uncl_obs_5_t_", length_strings[i]), 
         melt(eval(parse(text = paste0("cov_theta_uncl_obs_5_t_", length_strings[i])))))
  
  assign(paste0("cov_theta_uncl_obs_5_t_", length_strings[i]), 
         cbind(eval(parse(text = paste0("cov_theta_uncl_obs_5_t_", length_strings[i]))), 
               Length = length_strings[i], Nobs = "5",
               CILow = as.vector(eval(parse(text = paste0("results$sim_HMM_theta_uncl_obs_5_t_", length_strings[i], "$out_sim$emiss_low")))),
               CIHigh = as.vector(  eval(parse(text = paste0("results$sim_HMM_theta_uncl_obs_5_t_", length_strings[i], "$out_sim$emiss_up"))))))
  
  assign(paste0("cov_theta_uncl_obs_7_t_", length_strings[i]), 
         sweep(eval(parse(text = paste0("results$sim_HMM_theta_uncl_obs_7_t_", length_strings[i], "$out_sim$emiss_mean"))), 
               2, 0)) 
  
  assign(paste0("cov_theta_uncl_obs_7_t_", length_strings[i]), 
         melt(eval(parse(text = paste0("cov_theta_uncl_obs_7_t_", length_strings[i])))))
  
  assign(paste0("cov_theta_uncl_obs_7_t_", length_strings[i]), 
         cbind(eval(parse(text = paste0("cov_theta_uncl_obs_7_t_", length_strings[i]))), 
               Length = length_strings[i], Nobs = "7", 
               CILow = as.vector(eval(parse(text = paste0("results$sim_HMM_theta_uncl_obs_7_t_", length_strings[i], "$out_sim$emiss_low")))),
               CIHigh = as.vector(eval(parse(text = paste0("results$sim_HMM_theta_uncl_obs_7_t_", length_strings[i], "$out_sim$emiss_up"))))))
  
}

### Three observations.
uncl_obs_3_theta <- rbind(data.frame(cov_theta_uncl_obs_3_t_250), data.frame(cov_theta_uncl_obs_3_t_500),  
                          data.frame(cov_theta_uncl_obs_3_t_1000), data.frame(cov_theta_uncl_obs_3_t_2000), 
                          data.frame(cov_theta_uncl_obs_3_t_4000), data.frame(cov_theta_uncl_obs_3_t_8000))

uncl_obs_3_theta <- full_join(uncl_obs_3_theta, true_theta_three_obs, by = "X2")
uncl_obs_3_theta <- uncl_obs_3_theta %>% mutate(Coverage = ifelse((value > CILow & true < CIHigh),1,0))
colnames(uncl_obs_3_theta) <- c("Id", "S_to_obs", "Estimate", "Length", "Nobs","CiLow", "CiHigh","true","Coverage")

uncl_obs_3_theta$Length <- factor(uncl_obs_3_theta$Length, levels = c("250", "500", "1000", "2000", "4000", "8000"))
uncl_obs_3_theta <- aggregate(Coverage ~ Length + S_to_obs + Nobs, uncl_obs_3_theta, mean)

### Five observations.
uncl_obs_5_theta <- rbind(data.frame(cov_theta_uncl_obs_5_t_250), data.frame(cov_theta_uncl_obs_5_t_500),  
                          data.frame(cov_theta_uncl_obs_5_t_1000), data.frame(cov_theta_uncl_obs_5_t_2000), 
                          data.frame(cov_theta_uncl_obs_5_t_4000), data.frame(cov_theta_uncl_obs_5_t_8000))

uncl_obs_5_theta <- full_join(uncl_obs_5_theta, true_theta_five_obs, by = "X2")
uncl_obs_5_theta <- uncl_obs_5_theta %>% mutate(Coverage = ifelse((value > CILow & true < CIHigh),1,0))
colnames(uncl_obs_5_theta) <- c("Id", "S_to_obs", "Estimate", "Length", "Nobs","CiLow", "CiHigh","true","Coverage")

uncl_obs_5_theta$Length <- factor(uncl_obs_5_theta$Length, levels = c("250", "500", "1000", "2000", "4000", "8000"))
uncl_obs_5_theta <- aggregate(Coverage ~ Length + S_to_obs + Nobs, uncl_obs_5_theta, mean)

### Seven observations.
uncl_obs_7_theta <- rbind(data.frame(cov_theta_uncl_obs_7_t_250), data.frame(cov_theta_uncl_obs_7_t_500),  
                          data.frame(cov_theta_uncl_obs_7_t_1000), data.frame(cov_theta_uncl_obs_7_t_2000), 
                          data.frame(cov_theta_uncl_obs_7_t_4000), data.frame(cov_theta_uncl_obs_7_t_8000))

uncl_obs_7_theta <- full_join(uncl_obs_7_theta, true_theta_seven_obs, by = "X2")
uncl_obs_7_theta <- uncl_obs_7_theta %>% mutate(Coverage = ifelse((value > CILow & true < CIHigh),1,0))
colnames(uncl_obs_7_theta) <- c("Id", "S_to_obs", "Estimate", "Length", "Nobs","CiLow", "CiHigh","true","Coverage")

uncl_obs_7_theta$Length <- factor(uncl_obs_7_theta$Length, levels = c("250", "500", "1000", "2000", "4000", "8000"))
uncl_obs_7_theta$Nobs <- factor(uncl_obs_7_theta$Nobs, levels = c("3", "5", "7"))
uncl_obs_7_theta <- aggregate(Coverage ~ Length + S_to_obs + Nobs, uncl_obs_7_theta, mean)

# Plot of bias emission probabilities by number of observations and sequence length.
uncl_obs_theta <- rbind(data.frame(uncl_obs_3_theta), data.frame(uncl_obs_5_theta),
                        data.frame(uncl_obs_7_theta))
uncl_obs_theta$S_to_obs <- mapvalues(uncl_obs_theta$S_to_obs, 
                                     from = c("M_S1_cat1", "M_S1_cat2", "M_S1_cat3", "M_S1_cat4", "M_S1_cat5", "M_S1_cat6", "M_S1_cat7",
                                              "M_S2_cat1", "M_S2_cat2", "M_S2_cat3", "M_S2_cat4", "M_S2_cat5", "M_S2_cat6", "M_S2_cat7",
                                              "M_S3_cat1", "M_S3_cat2", "M_S3_cat3", "M_S3_cat4", "M_S3_cat5", "M_S3_cat6", "M_S3_cat7"), 
                                     to = c(rep(paste0("theta[", 1, "][", 1:7, "]"), each = 1), 
                                            rep(paste0("theta[", 2, "][", 1:7, "]"), each = 1), 
                                            rep(paste0("theta[", 3, "][", 1:7, "]"), each = 1)))
uncl_obs_theta$S_to_obs <- factor(uncl_obs_theta$S_to_obs, levels = c("theta[1][1]", "theta[1][2]", "theta[1][3]", "theta[1][4]", 
                                                                      "theta[1][5]", "theta[1][6]", "theta[1][7]",
                                                                      "theta[2][1]", "theta[2][2]", "theta[2][3]", "theta[2][4]", 
                                                                      "theta[2][5]", "theta[2][6]", "theta[2][7]",
                                                                      "theta[3][1]", "theta[3][2]", "theta[3][3]", "theta[3][4]", 
                                                                      "theta[3][5]", "theta[3][6]", "theta[3][7]"))
ggplot(uncl_obs_theta, aes(x = Length, y = Coverage, color = Nobs, group = Nobs)) +
  facet_wrap(facets = vars(S_to_obs), nrow = 3, ncol = 7, labeller = label_parsed) +
  geom_point() + geom_line() +  
  xlab("Sequence length") + 
  ylab("Bias") +
  geom_hline(yintercept = 0) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  ylim(0, 1) 


### Gamma. 
for (i in 1: length(length_strings)){
  
  assign(paste0("cov_gamma_uncl_obs_3_t_", length_strings[i]), 
         sweep(eval(parse(text = paste0("results$sim_HMM_theta_uncl_obs_3_t_", length_strings[i], "$out_sim$gamma_mean"))), 
               2, 0)) 
  
  assign(paste0("cov_gamma_uncl_obs_3_t_", length_strings[i]), 
         melt(eval(parse(text = paste0("cov_gamma_uncl_obs_3_t_", length_strings[i])))))
  
  assign(paste0("cov_gamma_uncl_obs_3_t_", length_strings[i]), 
         cbind(eval(parse(text = paste0("cov_gamma_uncl_obs_3_t_", length_strings[i]))), 
               Length = length_strings[i], Nobs = "3", 
               CILow = as.vector(eval(parse(text = paste0("results$sim_HMM_theta_uncl_obs_3_t_", length_strings[i], "$out_sim$gamma_low")))),
               CIHigh = as.vector(eval(parse(text = paste0("results$sim_HMM_theta_uncl_obs_3_t_", length_strings[i], "$out_sim$gamma_up"))))))
               
  
  assign(paste0("cov_gamma_uncl_obs_5_t_", length_strings[i]), 
         sweep(eval(parse(text = paste0("results$sim_HMM_theta_uncl_obs_5_t_", length_strings[i], "$out_sim$gamma_mean"))), 
               2, 0)) 
  
  assign(paste0("cov_gamma_uncl_obs_5_t_", length_strings[i]), 
         melt(eval(parse(text = paste0("cov_gamma_uncl_obs_5_t_", length_strings[i])))))
  
  assign(paste0("cov_gamma_uncl_obs_5_t_", length_strings[i]), 
         cbind(eval(parse(text = paste0("cov_gamma_uncl_obs_5_t_", length_strings[i]))), 
               Length = length_strings[i], Nobs = "5", 
               CILow = as.vector(eval(parse(text = paste0("results$sim_HMM_theta_uncl_obs_5_t_", length_strings[i], "$out_sim$gamma_low")))),
               CIHigh = as.vector(eval(parse(text = paste0("results$sim_HMM_theta_uncl_obs_5_t_", length_strings[i], "$out_sim$gamma_up"))))))

  assign(paste0("cov_gamma_uncl_obs_7_t_", length_strings[i]), 
         sweep(eval(parse(text = paste0("results$sim_HMM_theta_uncl_obs_7_t_", length_strings[i], "$out_sim$gamma_mean"))), 
               2, 0)) 
  
  assign(paste0("cov_gamma_uncl_obs_7_t_", length_strings[i]), 
         melt(eval(parse(text = paste0("cov_gamma_uncl_obs_7_t_", length_strings[i])))))
  
  assign(paste0("cov_gamma_uncl_obs_7_t_", length_strings[i]), 
         cbind(eval(parse(text = paste0("cov_gamma_uncl_obs_7_t_", length_strings[i]))), 
               Length = length_strings[i], Nobs = "7", 
               CILow = as.vector(eval(parse(text = paste0("results$sim_HMM_theta_uncl_obs_7_t_", length_strings[i], "$out_sim$gamma_low")))),
               CIHigh = as.vector(eval(parse(text = paste0("results$sim_HMM_theta_uncl_obs_7_t_", length_strings[i], "$out_sim$gamma_up"))))))
  
  

}

# 3 observations.
uncl_obs_3_gamma <- rbind(data.frame(cov_gamma_uncl_obs_3_t_250), data.frame(cov_gamma_uncl_obs_3_t_500),  
                          data.frame(cov_gamma_uncl_obs_3_t_1000), data.frame(cov_gamma_uncl_obs_3_t_2000), 
                          data.frame(cov_gamma_uncl_obs_3_t_4000), data.frame(cov_gamma_uncl_obs_3_t_8000))

uncl_obs_3_gamma <- full_join(uncl_obs_3_gamma, gamma_true, by = "X2")
uncl_obs_3_gamma <- uncl_obs_3_gamma %>% mutate(Coverage = ifelse((value > CILow & true < CIHigh),1,0))
colnames(uncl_obs_3_gamma) <- c("Id", "S_to_obs", "Estimate", "Length", "Nobs","CiLow", "CiHigh","true","Coverage")

uncl_obs_3_gamma$Length <- factor(uncl_obs_3_gamma$Length, levels = c("250", "500", "1000", "2000", "4000", "8000"))
uncl_obs_3_gamma <- aggregate(Coverage ~ Length + S_to_obs + Nobs, uncl_obs_3_gamma, mean)

# 5 observations. 
uncl_obs_5_gamma <- rbind(data.frame(cov_gamma_uncl_obs_5_t_250), data.frame(cov_gamma_uncl_obs_5_t_500),  
                          data.frame(cov_gamma_uncl_obs_5_t_1000), data.frame(cov_gamma_uncl_obs_5_t_2000), 
                          data.frame(cov_gamma_uncl_obs_5_t_4000), data.frame(cov_gamma_uncl_obs_5_t_8000))

uncl_obs_5_gamma <- full_join(uncl_obs_5_gamma, gamma_true, by = "X2")
uncl_obs_5_gamma <- uncl_obs_5_gamma %>% mutate(Coverage = ifelse((value > CILow & true < CIHigh),1,0))
colnames(uncl_obs_5_gamma) <- c("Id", "S_to_obs", "Estimate", "Length", "Nobs","CiLow", "CiHigh","true","Coverage")

uncl_obs_5_gamma$Length <- factor(uncl_obs_5_gamma$Length, levels = c("250", "500", "1000", "2000", "4000", "8000"))
uncl_obs_5_gamma <- aggregate(Coverage ~ Length + S_to_obs + Nobs, uncl_obs_5_gamma, mean)

# 7 observations. 
uncl_obs_7_gamma <- rbind(data.frame(cov_gamma_uncl_obs_7_t_250), data.frame(cov_gamma_uncl_obs_7_t_500),  
                          data.frame(cov_gamma_uncl_obs_7_t_1000), data.frame(cov_gamma_uncl_obs_7_t_2000), 
                          data.frame(cov_gamma_uncl_obs_7_t_4000), data.frame(cov_gamma_uncl_obs_7_t_8000))

uncl_obs_7_gamma <- full_join(uncl_obs_7_gamma, gamma_true, by = "X2")
uncl_obs_7_gamma <- uncl_obs_7_gamma %>% mutate(Coverage = ifelse((value > CILow & true < CIHigh),1,0))
colnames(uncl_obs_7_gamma) <- c("Id", "S_to_obs", "Estimate", "Length", "Nobs","CiLow", "CiHigh","true","Coverage")

uncl_obs_7_gamma$Length <- factor(uncl_obs_7_gamma$Length, levels = c("250", "500", "1000", "2000", "4000", "8000"))
uncl_obs_7_gamma <- aggregate(Coverage ~ Length + S_to_obs + Nobs, uncl_obs_7_gamma, mean)

# Plot of bias transition probabilities by number of observations and sequence length,
num_obs_gamma <- rbind(data.frame(uncl_obs_3_gamma), data.frame(uncl_obs_5_gamma),
                       data.frame(uncl_obs_7_gamma))
num_obs_gamma$S_to_obs <- mapvalues(num_obs_gamma$S_to_obs, 
                                  from = c("M_S1_to_S1", "M_S1_to_S2", "M_S1_to_S3", 
                                           "M_S2_to_S1", "M_S2_to_S2", "M_S2_to_S3",
                                           "M_S3_to_S1", "M_S3_to_S2", "M_S3_to_S3"), 
                                  to = c(rep(paste0("gamma[", 1, "][", 1:3, "]"), each = 1), 
                                         rep(paste0("gamma[", 2, "][", 1:3, "]"), each = 1), 
                                         rep(paste0("gamma[", 3, "][", 1:3, "]"), each = 1)))
num_obs_gamma$Nobs <- factor(num_obs_gamma$Nobs, levels = c("3", "5", "7"))
ggplot(num_obs_gamma, aes(x = Length, y = Coverage, color = Nobs, group = Nobs)) +
  facet_wrap(facets = vars(S_to_obs), nrow = 3, ncol = 3, labeller = label_parsed) +
  geom_point() + geom_line() +  
  xlab("Sequence length") + 
  ylab("Coverage") +
  geom_hline(yintercept = 0) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  ylim(0, 1)


##############
### Varobs ###
##############

### Theta.

## True emission probabilities.

# Three observations.
true_theta_three_obs <- c(0.80, 0.10, 0.10,
                          0.10, 0.80, 0.10,
                          0.10, 0.10, 0.80)

true_theta_three_obs = data_frame(true_theta_three_obs) %>% mutate(c("M_S1_cat1", "M_S1_cat2",  "M_S1_cat3",
                                                                     "M_S2_cat1", "M_S2_cat2",  "M_S2_cat3", 
                                                                     "M_S3_cat1", "M_S3_cat2",  "M_S3_cat3"))
colnames(true_theta_three_obs) <- c("true", "X2")

# Five observations.
true_theta_five_obs <- c(0.92, 0.02, 0.02, 0.02, 0.02,
                         0.02, 0.47, 0.02, 0.47, 0.02,
                         0.02, 0.02, 0.62, 0.02, 0.32)

true_theta_five_obs = data_frame(true_theta_five_obs) %>% mutate(c("M_S1_cat1", "M_S1_cat2",  "M_S1_cat3", "M_S1_cat4", "M_S1_cat5",
                                                                     "M_S2_cat1", "M_S2_cat2",  "M_S2_cat3", "M_S2_cat4", "M_S2_cat5",
                                                                     "M_S3_cat1", "M_S3_cat2",  "M_S3_cat3","M_S3_cat4", "M_S3_cat5"))
colnames(true_theta_five_obs) <- c("true", "X2")

# Seven observations.
true_theta_seven_obs <- c(0.76, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04,
                          0.02, 0.45, 0.02, 0.45, 0.02, 0.02, 0.02,
                          0.01, 0.02, 0.31, 0.02, 0.31, 0.02, 0.31)

true_theta_seven_obs = data_frame(true_theta_seven_obs) %>% mutate(c("M_S1_cat1", "M_S1_cat2",  "M_S1_cat3", "M_S1_cat4", "M_S1_cat5", "M_S1_cat6", "M_S1_cat7",
                                                                     "M_S2_cat1", "M_S2_cat2",  "M_S2_cat3", "M_S2_cat4", "M_S2_cat5", "M_S2_cat6", "M_S2_cat7",
                                                                     "M_S3_cat1", "M_S3_cat2",  "M_S3_cat3", "M_S3_cat4", "M_S3_cat5", "M_S3_cat6", "M_S3_cat7"))
colnames(true_theta_seven_obs) <- c("true", "X2")

for (i in 1: length(length_strings)){
  
  assign(paste0("cov_theta_var_obs_3_t_", length_strings[i]), 
         sweep(eval(parse(text = paste0("results$sim_HMM_varobs_obs_3_t_", length_strings[i], "$out_sim$emiss_mean"))), 
               2, 0)) 
  
  assign(paste0("cov_theta_var_obs_3_t_", length_strings[i]), 
         melt(eval(parse(text = paste0("cov_theta_var_obs_3_t_", length_strings[i])))))
  
  assign(paste0("cov_theta_var_obs_3_t_", length_strings[i]), 
         cbind(eval(parse(text = paste0("cov_theta_var_obs_3_t_", length_strings[i]))), 
               Length = length_strings[i], Nobs = "3", 
               CILow = as.vector(eval(parse(text = paste0("results$sim_HMM_varobs_obs_3_t_", length_strings[i], "$out_sim$emiss_low")))),
               CIHigh = as.vector(  eval(parse(text = paste0("results$sim_HMM_varobs_obs_3_t_", length_strings[i], "$out_sim$emiss_up"))))))
  
  assign(paste0("cov_theta_var_obs_5_t_", length_strings[i]), 
         sweep(eval(parse(text = paste0("results$sim_HMM_varobs_obs_5_t_", length_strings[i], "$out_sim$emiss_mean"))), 
               2, 0)) 
  
  assign(paste0("cov_theta_var_obs_5_t_", length_strings[i]), 
         melt(eval(parse(text = paste0("cov_theta_var_obs_5_t_", length_strings[i])))))
  
  assign(paste0("cov_theta_var_obs_5_t_", length_strings[i]), 
         cbind(eval(parse(text = paste0("cov_theta_var_obs_5_t_", length_strings[i]))), 
               Length = length_strings[i], Nobs = "5", 
               CILow = as.vector(eval(parse(text = paste0("results$sim_HMM_varobs_obs_5_t_", length_strings[i], "$out_sim$emiss_low")))),
               CIHigh = as.vector(  eval(parse(text = paste0("results$sim_HMM_varobs_obs_5_t_", length_strings[i], "$out_sim$emiss_up"))))))
  
  assign(paste0("cov_theta_var_obs_7_t_", length_strings[i]), 
         sweep(eval(parse(text = paste0("results$sim_HMM_varobs_obs_7_t_", length_strings[i], "$out_sim$emiss_mean"))), 
               2, 0)) 
  
  assign(paste0("cov_theta_var_obs_7_t_", length_strings[i]), 
         melt(eval(parse(text = paste0("cov_theta_var_obs_7_t_", length_strings[i])))))
  
  assign(paste0("cov_theta_var_obs_7_t_", length_strings[i]), 
         cbind(eval(parse(text = paste0("cov_theta_var_obs_7_t_", length_strings[i]))), 
               Length = length_strings[i], Nobs = "7", 
               CILow = as.vector(eval(parse(text = paste0("results$sim_HMM_varobs_obs_7_t_", length_strings[i], "$out_sim$emiss_low")))),
               CIHigh = as.vector(eval(parse(text = paste0("results$sim_HMM_varobs_obs_7_t_", length_strings[i], "$out_sim$emiss_up"))))))
  
}

### Three observations.
var_obs_3_theta <- rbind(data.frame(cov_theta_var_obs_3_t_250), data.frame(cov_theta_var_obs_3_t_500),  
                         data.frame(cov_theta_var_obs_3_t_1000), data.frame(cov_theta_var_obs_3_t_2000), 
                         data.frame(cov_theta_var_obs_3_t_4000), data.frame(cov_theta_var_obs_3_t_8000))

var_obs_3_theta <- full_join(var_obs_3_theta, true_theta_three_obs, by = "X2")

var_obs_3_theta <- var_obs_3_theta %>% mutate(Coverage = ifelse((value > CILow & true < CIHigh),1,0))
colnames(var_obs_3_theta) <- c("Id", "S_to_obs", "Estimate", "Length", "Nobs","CiLow", "CiHigh","true","Coverage")

var_obs_3_theta$Length <- factor(var_obs_3_theta$Length, levels = c("250", "500", "1000", "2000", "4000", "8000"))
var_obs_3_theta <- aggregate(Coverage ~ Length + S_to_obs + Nobs, var_obs_3_theta, mean)

### Five observations.
var_obs_5_theta <- rbind(data.frame(cov_theta_var_obs_5_t_250), data.frame(cov_theta_var_obs_5_t_500),  
                         data.frame(cov_theta_var_obs_5_t_1000), data.frame(cov_theta_var_obs_5_t_2000), 
                         data.frame(cov_theta_var_obs_5_t_4000), data.frame(cov_theta_var_obs_5_t_8000))

var_obs_5_theta <- full_join(var_obs_5_theta, true_theta_five_obs, by = "X2")

var_obs_5_theta <- var_obs_5_theta %>% mutate(Coverage = ifelse((value > CILow & true < CIHigh),1,0))
colnames(var_obs_5_theta) <- c("Id", "S_to_obs", "Estimate", "Length", "Nobs","CiLow", "CiHigh","true","Coverage")

var_obs_5_theta$Length <- factor(var_obs_5_theta$Length, levels = c("250", "500", "1000", "2000", "4000", "8000"))
var_obs_5_theta <- aggregate(Coverage ~ Length + S_to_obs + Nobs, var_obs_5_theta, mean)

### Seven observations.
var_obs_7_theta <- rbind(data.frame(cov_theta_var_obs_7_t_250), data.frame(cov_theta_var_obs_7_t_500),  
                         data.frame(cov_theta_var_obs_7_t_1000), data.frame(cov_theta_var_obs_7_t_2000), 
                         data.frame(cov_theta_var_obs_7_t_4000), data.frame(cov_theta_var_obs_7_t_8000))

var_obs_7_theta <- full_join(var_obs_7_theta, true_theta_seven_obs, by = "X2")

var_obs_7_theta <- var_obs_7_theta %>% mutate(Coverage = ifelse((value > CILow & true < CIHigh),1,0))
colnames(var_obs_7_theta) <- c("Id", "S_to_obs", "Estimate", "Length", "Nobs","CiLow", "CiHigh","true","Coverage")

var_obs_7_theta$Length <- factor(var_obs_7_theta$Length, levels = c("250", "500", "1000", "2000", "4000", "8000"))
var_obs_7_theta <- aggregate(Coverage ~ Length + S_to_obs + Nobs, var_obs_7_theta, mean)

# Plot of bias emission probabilities by number of observations and sequence length.
var_obs_theta <- rbind(data.frame(var_obs_3_theta), data.frame(var_obs_5_theta),
                       data.frame(var_obs_7_theta))
var_obs_theta$S_to_obs <- mapvalues(var_obs_theta$S_to_obs, 
                                    from = c("M_S1_cat1", "M_S1_cat2", "M_S1_cat3", "M_S1_cat4", "M_S1_cat5", "M_S1_cat6", "M_S1_cat7",
                                             "M_S2_cat1", "M_S2_cat2", "M_S2_cat3", "M_S2_cat4", "M_S2_cat5", "M_S2_cat6", "M_S2_cat7",
                                             "M_S3_cat1", "M_S3_cat2", "M_S3_cat3", "M_S3_cat4", "M_S3_cat5", "M_S3_cat6", "M_S3_cat7"), 
                                    to = c(rep(paste0("theta[", 1, "][", 1:7, "]"), each = 1), 
                                           rep(paste0("theta[", 2, "][", 1:7, "]"), each = 1), 
                                           rep(paste0("theta[", 3, "][", 1:7, "]"), each = 1)))
var_obs_theta$S_to_obs <- factor(var_obs_theta$S_to_obs, levels = c("theta[1][1]", "theta[1][2]", "theta[1][3]", "theta[1][4]", 
                                                                    "theta[1][5]", "theta[1][6]", "theta[1][7]",
                                                                    "theta[2][1]", "theta[2][2]", "theta[2][3]", "theta[2][4]", 
                                                                    "theta[2][5]", "theta[2][6]", "theta[2][7]",
                                                                    "theta[3][1]", "theta[3][2]", "theta[3][3]", "theta[3][4]", 
                                                                    "theta[3][5]", "theta[3][6]", "theta[3][7]"))
ggplot(var_obs_theta, aes(x = Length, y = Coverage, color = Nobs, group = Nobs)) +
  facet_wrap(facets = vars(S_to_obs), nrow = 3, ncol = 7, labeller = label_parsed) +
  geom_point() + geom_line() +  
  xlab("Sequence length") + 
  ylab("Coverage") +
  geom_hline(yintercept = 0) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  ylim(0.5, 1) 


### Gamma. 

gamma_true <- c(0.80, 0.10, 0.10, 
                0.10, 0.80, 0.10, 
                0.10, 0.10, 0.80)

gamma_true = data_frame(gamma_true) %>% mutate(c("M_S1_to_S1", "M_S1_to_S2", "M_S1_to_S3",
                                                 "M_S2_to_S1", "M_S2_to_S2", "M_S2_to_S3", 
                                                 "M_S3_to_S1", "M_S3_to_S2", "M_S3_to_S3")) 

colnames(gamma_true) <- c("true", "X2")

for (i in 1: length(length_strings)){
  
  assign(paste0("cov_gamma_var_obs_3_t_", length_strings[i]), 
         sweep(eval(parse(text = paste0("results$sim_HMM_varobs_obs_3_t_", length_strings[i], "$out_sim$gamma_mean"))), 
               2, 0)) 
  
  assign(paste0("cov_gamma_var_obs_3_t_", length_strings[i]), 
         melt(eval(parse(text = paste0("cov_gamma_var_obs_3_t_", length_strings[i])))))
  
  assign(paste0("cov_gamma_var_obs_3_t_", length_strings[i]), 
         cbind(eval(parse(text = paste0("cov_gamma_var_obs_3_t_", length_strings[i]))), 
               Length = length_strings[i], Nobs = "3",
               CILow = as.vector(eval(parse(text = paste0("results$sim_HMM_varobs_obs_3_t_", length_strings[i], "$out_sim$gamma_low")))),
               CIHigh = as.vector(eval(parse(text = paste0("results$sim_HMM_varobs_obs_3_t_", length_strings[i], "$out_sim$gamma_up"))))))
 
   
  assign(paste0("cov_gamma_var_obs_5_t_", length_strings[i]), 
         sweep(eval(parse(text = paste0("results$sim_HMM_varobs_obs_5_t_", length_strings[i], "$out_sim$gamma_mean"))), 
               2, 0)) 
  
  assign(paste0("cov_gamma_var_obs_5_t_", length_strings[i]), 
         melt(eval(parse(text = paste0("cov_gamma_var_obs_5_t_", length_strings[i])))))
  
  assign(paste0("cov_gamma_var_obs_5_t_", length_strings[i]), 
         cbind(eval(parse(text = paste0("cov_gamma_var_obs_5_t_", length_strings[i]))), 
               Length = length_strings[i], Nobs = "5",
               CILow = as.vector(eval(parse(text = paste0("results$sim_HMM_varobs_obs_5_t_", length_strings[i], "$out_sim$gamma_low")))),
               CIHigh = as.vector(eval(parse(text = paste0("results$sim_HMM_varobs_obs_5_t_", length_strings[i], "$out_sim$gamma_up"))))))
  
  
  assign(paste0("cov_gamma_var_obs_7_t_", length_strings[i]), 
         sweep(eval(parse(text = paste0("results$sim_HMM_varobs_obs_7_t_", length_strings[i], "$out_sim$gamma_mean"))), 
               2, 0)) 
  
  assign(paste0("cov_gamma_var_obs_7_t_", length_strings[i]), 
         melt(eval(parse(text = paste0("cov_gamma_var_obs_7_t_", length_strings[i])))))
  
  assign(paste0("cov_gamma_var_obs_7_t_", length_strings[i]), 
         cbind(eval(parse(text = paste0("cov_gamma_var_obs_7_t_", length_strings[i]))), 
               Length = length_strings[i], Nobs = "7",
               CILow = as.vector(eval(parse(text = paste0("results$sim_HMM_varobs_obs_7_t_", length_strings[i], "$out_sim$gamma_low")))),
               CIHigh = as.vector(eval(parse(text = paste0("results$sim_HMM_varobs_obs_7_t_", length_strings[i], "$out_sim$gamma_up"))))))
  
}

# 3 observations.
var_obs_3_gamma <- rbind(data.frame(cov_gamma_var_obs_3_t_250), data.frame(cov_gamma_var_obs_3_t_500),  
                         data.frame(cov_gamma_var_obs_3_t_1000), data.frame(cov_gamma_var_obs_3_t_2000), 
                         data.frame(cov_gamma_var_obs_3_t_4000), data.frame(cov_gamma_var_obs_3_t_8000))

var_obs_3_gamma <- full_join(var_obs_3_gamma, gamma_true, by = "X2")

var_obs_3_gamma <- var_obs_3_gamma %>% mutate(Coverage = ifelse((true > CILow & true < CIHigh),1,0))
colnames(var_obs_3_gamma) <- c("Id", "S_to_s", "Estimate", "Length", "Nobs","CiLow", "CiHigh","true","Coverage")

var_obs_3_gamma$Length <- factor(var_obs_3_gamma$Length, levels = c("250", "500", "1000", "2000", "4000", "8000"))
var_obs_3_gamma <- aggregate(Coverage ~ Length + S_to_s + Nobs, var_obs_3_gamma, mean)

# 5 observations. 
var_obs_5_gamma <- rbind(data.frame(cov_gamma_var_obs_5_t_250), data.frame(cov_gamma_var_obs_5_t_500),  
                         data.frame(cov_gamma_var_obs_5_t_1000), data.frame(cov_gamma_var_obs_5_t_2000), 
                         data.frame(cov_gamma_var_obs_5_t_4000), data.frame(cov_gamma_var_obs_5_t_8000))

var_obs_5_gamma <- full_join(var_obs_5_gamma, gamma_true, by = "X2")

var_obs_5_gamma <- var_obs_5_gamma %>% mutate(Coverage = ifelse((true > CILow & true < CIHigh),1,0))
colnames(var_obs_5_gamma) <- c("Id", "S_to_s", "Estimate", "Length", "Nobs","CiLow", "CiHigh","true","Coverage")

var_obs_5_gamma$Length <- factor(var_obs_5_gamma$Length, levels = c("250", "500", "1000", "2000", "4000", "8000"))
var_obs_5_gamma <- aggregate(Coverage ~ Length + S_to_s + Nobs, var_obs_5_gamma, mean)

# 7 observations. 
var_obs_7_gamma <- rbind(data.frame(cov_gamma_var_obs_7_t_250), data.frame(cov_gamma_var_obs_7_t_500),  
                         data.frame(cov_gamma_var_obs_7_t_1000), data.frame(cov_gamma_var_obs_7_t_2000), 
                         data.frame(cov_gamma_var_obs_7_t_4000), data.frame(cov_gamma_var_obs_7_t_8000))

var_obs_7_gamma <- full_join(var_obs_7_gamma, gamma_true, by = "X2")

var_obs_7_gamma <- var_obs_7_gamma %>% mutate(Coverage = ifelse((true > CILow & true < CIHigh),1,0))
colnames(var_obs_7_gamma) <- c("Id", "S_to_s", "Estimate", "Length", "Nobs","CiLow", "CiHigh","true","Coverage")

var_obs_7_gamma$Length <- factor(var_obs_7_gamma$Length, levels = c("250", "500", "1000", "2000", "4000", "8000"))
var_obs_7_gamma <- aggregate(Coverage ~ Length + S_to_s + Nobs, var_obs_7_gamma, mean)

# Plot of bias transition probabilities by number of observations and sequence length.
num_obs_gamma <- rbind(data.frame(var_obs_3_gamma), data.frame(var_obs_5_gamma),
                       data.frame(var_obs_7_gamma))
num_obs_gamma$S_to_s <- mapvalues(num_obs_gamma$S_to_s, 
                                  from = c("M_S1_to_S1", "M_S1_to_S2", "M_S1_to_S3", 
                                           "M_S2_to_S1", "M_S2_to_S2", "M_S2_to_S3",
                                           "M_S3_to_S1", "M_S3_to_S2", "M_S3_to_S3"), 
                                  to = c(rep(paste0("gamma[", 1, "][", 1:3, "]"), each = 1), 
                                         rep(paste0("gamma[", 2, "][", 1:3, "]"), each = 1), 
                                         rep(paste0("gamma[", 3, "][", 1:3, "]"), each = 1)))
num_obs_gamma$Nobs <- factor(num_obs_gamma$Nobs, levels = c("3", "5", "7"))
ggplot(num_obs_gamma, aes(x = Length, y = Coverage, color = Nobs, group = Nobs)) +
  facet_wrap(facets = vars(S_to_s), nrow = 3, ncol = 3, labeller = label_parsed) +
  geom_point() + geom_line() +  
  xlab("Sequence length") + 
  ylab("Coverage") +
  geom_hline(yintercept = 0) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  ylim(0.5, 1)

### Save plots to dir. 
plots.dir.path <- list.files(tempdir(), pattern="rs-graphics", full.names = TRUE); 
plots.png.paths <- list.files(plots.dir.path, pattern=".png", full.names = TRUE)
file.copy(from=plots.png.paths, to="/Users/bart-jan/Documents/GitHub/paper-hmm-simulation-study/plots")

plots.png.detials <- file.info(plots.png.paths)
plots.png.detials <- plots.png.detials[order(plots.png.detials$mtime),]
sorted.png.names <- gsub(plots.dir.path, "/Users/bart-jan/Documents/GitHub/paper-hmm-simulation-study/plots", row.names(plots.png.detials), fixed=TRUE)
numbered.png.names <- paste0("/Users/bart-jan/Documents/GitHub/paper-hmm-simulation-study/plots", 1:length(sorted.png.names), ".png")

# Rename all the .png files as: 1.png, 2.png, 3.png, and so on.
file.rename(from=sorted.png.names, to=numbered.png.names)
