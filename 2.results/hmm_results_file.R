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
library(data.table)
library(reshape)
library(ggplot2)
library(ggh4x)
library(plyr)

### Set working directory. 
setwd("C:\\Academia\\HMM paper\\Data\\Data")

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
# Moderately clear. 
theta_mod_true <- c(0.68, 0.08, 0.08, 0.08, 0.08,
                    0.08, 0.38, 0.08, 0.38, 0.08,
                    0.08, 0.08, 0.38, 0.08, 0.38)
# Unclear. 
theta_unclear_true <- c(0.44, 0.14, 0.14, 0.14, 0.14,
                        0.14, 0.29, 0.14, 0.29, 0.14,
                        0.14, 0.14, 0.29, 0.14, 0.29)

for (i in 1: length(length_strings)){
  
  assign(paste0("bias_theta_cl_obs_5_t_", length_strings[i]), 
         sweep(eval(parse(text = paste0("results$sim_HMM_theta_cl_obs_5_t_", length_strings[i], "$out_sim$emiss_mean"))), 
               2, theta_clear_true)) 
  
  assign(paste0("bias_theta_cl_obs_5_t_", length_strings[i]), 
         melt(eval(parse(text = paste0("bias_theta_cl_obs_5_t_", length_strings[i])))))
  
  assign(paste0("bias_theta_cl_obs_5_t_", length_strings[i]), 
         cbind(eval(parse(text = paste0("bias_theta_cl_obs_5_t_", length_strings[i]))), 
               Length = length_strings[i], Clarity = "Clear"))
  
  assign(paste0("bias_theta_modcl_obs_5_t_", length_strings[i]), 
         sweep(eval(parse(text = paste0("results$sim_HMM_theta_modcl_obs_5_t_", length_strings[i], "$out_sim$emiss_mean"))), 
               2, theta_mod_true)) 
  
  assign(paste0("bias_theta_modcl_obs_5_t_", length_strings[i]), 
         melt(eval(parse(text = paste0("bias_theta_modcl_obs_5_t_", length_strings[i])))))
  
  assign(paste0("bias_theta_modcl_obs_5_t_", length_strings[i]), 
         cbind(eval(parse(text = paste0("bias_theta_modcl_obs_5_t_", length_strings[i]))), 
               Length = length_strings[i], Clarity = "Moderate"))
  
  assign(paste0("bias_theta_uncl_obs_5_t_", length_strings[i]), 
         sweep(eval(parse(text = paste0("results$sim_HMM_theta_uncl_obs_5_t_", length_strings[i], "$out_sim$emiss_mean"))), 
               2, theta_unclear_true)) 
  
  assign(paste0("bias_theta_uncl_obs_5_t_", length_strings[i]), 
         melt(eval(parse(text = paste0("bias_theta_uncl_obs_5_t_", length_strings[i])))))
  
  assign(paste0("bias_theta_uncl_obs_5_t_", length_strings[i]), 
         cbind(eval(parse(text = paste0("bias_theta_uncl_obs_5_t_", length_strings[i]))), 
               Length = length_strings[i], Clarity = "Unclear"))
  
}

## Clear.
cl_obs_5_theta <- rbind(data.frame(bias_theta_cl_obs_5_t_250), data.frame(bias_theta_cl_obs_5_t_500),  
                        data.frame(bias_theta_cl_obs_5_t_1000), data.frame(bias_theta_cl_obs_5_t_2000), 
                        data.frame(bias_theta_cl_obs_5_t_4000), data.frame(bias_theta_cl_obs_5_t_8000))
colnames(cl_obs_5_theta) <- c("Id", "S_to_obs", "Abs_bias", "Length", "Clarity")
cl_obs_5_theta$Length <- factor(cl_obs_5_theta$Length, levels = c("250", "500", "1000", "2000", "4000", "8000"))
cl_obs_5_theta <- aggregate(Abs_bias ~ Length + S_to_obs + Clarity, cl_obs_5_theta, mean)

## Moderately clear. 
modcl_obs_5_theta <- rbind(data.frame(bias_theta_modcl_obs_5_t_250), data.frame(bias_theta_modcl_obs_5_t_500),  
                           data.frame(bias_theta_modcl_obs_5_t_1000), data.frame(bias_theta_modcl_obs_5_t_2000), 
                           data.frame(bias_theta_modcl_obs_5_t_4000), data.frame(bias_theta_modcl_obs_5_t_8000))
colnames(modcl_obs_5_theta) <- c("Id", "S_to_obs", "Abs_bias", "Length", "Clarity")
modcl_obs_5_theta$Length <- factor(modcl_obs_5_theta$Length, levels = c("250", "500", "1000", "2000", "4000", "8000"))
modcl_obs_5_theta <- aggregate(Abs_bias ~ Length + S_to_obs + Clarity, modcl_obs_5_theta, mean)

## Unclear. 
uncl_obs_5_theta <- rbind(data.frame(bias_theta_uncl_obs_5_t_250), data.frame(bias_theta_uncl_obs_5_t_500),  
                          data.frame(bias_theta_uncl_obs_5_t_1000), data.frame(bias_theta_uncl_obs_5_t_2000), 
                          data.frame(bias_theta_uncl_obs_5_t_4000), data.frame(bias_theta_uncl_obs_5_t_8000))
colnames(uncl_obs_5_theta) <- c("Id", "S_to_obs", "Abs_bias", "Length", "Clarity")
uncl_obs_5_theta$Length <- factor(uncl_obs_5_theta$Length, levels = c("250", "500", "1000", "2000", "4000", "8000"))
uncl_obs_5_theta <- aggregate(Abs_bias ~ Length + S_to_obs + Clarity, uncl_obs_5_theta, mean)

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
ggplot(obs_5_theta, aes(x = Length, y = Abs_bias, color = Clarity, group = Clarity)) +
  facet_wrap(facets = vars(S_to_obs), nrow = 3, ncol = 5, labeller = label_parsed) +
  geom_point() + 
  geom_line() +  
  xlab("Sequence length") + 
  ylab("Bias") +
  geom_hline(yintercept = 0) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  ylim(-0.4, 0.25) 

## Standard deviation.
for (i in 1: length(length_strings)){
  
  assign(paste0("sd_theta_cl_obs_5_t_", length_strings[i]), 
         melt(eval(parse(text = paste0("results$sim_HMM_theta_cl_obs_5_t_", length_strings[i], "$out_sim$emiss_sd")))))
  
  assign(paste0("sd_theta_cl_obs_5_t_", length_strings[i]), 
         cbind(eval(parse(text = paste0("sd_theta_cl_obs_5_t_", length_strings[i]))), 
               Length = length_strings[i], Clarity = "Clear"))
  
  assign(paste0("sd_theta_modcl_obs_5_t_", length_strings[i]), 
         melt(eval(parse(text = paste0("results$sim_HMM_theta_modcl_obs_5_t_", length_strings[i], "$out_sim$emiss_sd")))))
  
  assign(paste0("sd_theta_modcl_obs_5_t_", length_strings[i]), 
         cbind(eval(parse(text = paste0("sd_theta_modcl_obs_5_t_", length_strings[i]))), 
               Length = length_strings[i], Clarity = "Moderate"))
  
  assign(paste0("sd_theta_uncl_obs_5_t_", length_strings[i]), 
         melt(eval(parse(text = paste0("results$sim_HMM_theta_uncl_obs_5_t_", length_strings[i], "$out_sim$emiss_sd")))))
  
  assign(paste0("sd_theta_uncl_obs_5_t_", length_strings[i]), 
         cbind(eval(parse(text = paste0("sd_theta_uncl_obs_5_t_", length_strings[i]))), 
               Length = length_strings[i], Clarity = "Unclear"))
  
}

# Clear.
sd_cl_obs_5_theta <- rbind(data.frame(sd_theta_cl_obs_5_t_250), data.frame(sd_theta_cl_obs_5_t_500),  
                           data.frame(sd_theta_cl_obs_5_t_1000), data.frame(sd_theta_cl_obs_5_t_2000), 
                           data.frame(sd_theta_cl_obs_5_t_4000), data.frame(sd_theta_cl_obs_5_t_8000))
colnames(sd_cl_obs_5_theta) <- c("Id", "S_to_obs", "Sd", "Length", "Clarity")
sd_cl_obs_5_theta$Length <- factor(sd_cl_obs_5_theta$Length, levels = c("250", "500", "1000", "2000", "4000", "8000"))
sd_cl_obs_5_theta <- aggregate(Sd ~ Length + S_to_obs + Clarity, sd_cl_obs_5_theta, mean)

# Moderately clear. 
sd_modcl_obs_5_theta <- rbind(data.frame(sd_theta_modcl_obs_5_t_250), data.frame(sd_theta_modcl_obs_5_t_500),  
                              data.frame(sd_theta_modcl_obs_5_t_1000), data.frame(sd_theta_modcl_obs_5_t_2000), 
                              data.frame(sd_theta_modcl_obs_5_t_4000), data.frame(sd_theta_modcl_obs_5_t_8000))
colnames(sd_modcl_obs_5_theta) <- c("Id", "S_to_obs", "Sd", "Length", "Clarity")
sd_modcl_obs_5_theta$Length <- factor(sd_modcl_obs_5_theta$Length, levels = c("250", "500", "1000", "2000", "4000", "8000"))
sd_modcl_obs_5_theta <- aggregate(Sd ~ Length + S_to_obs + Clarity, sd_modcl_obs_5_theta, mean)

# Unclear. 
sd_uncl_obs_5_theta <- rbind(data.frame(sd_theta_uncl_obs_5_t_250), data.frame(sd_theta_uncl_obs_5_t_500),  
                             data.frame(sd_theta_uncl_obs_5_t_1000), data.frame(sd_theta_uncl_obs_5_t_2000), 
                             data.frame(sd_theta_uncl_obs_5_t_4000), data.frame(sd_theta_uncl_obs_5_t_8000))
colnames(sd_uncl_obs_5_theta) <- c("Id", "S_to_obs", "Sd", "Length", "Clarity")
sd_uncl_obs_5_theta$Length <- factor(sd_uncl_obs_5_theta$Length, levels = c("250", "500", "1000", "2000", "4000", "8000"))
sd_uncl_obs_5_theta <- aggregate(Sd ~ Length + S_to_obs + Clarity, sd_uncl_obs_5_theta, mean)

# Plot of standard deviation emission probabilities by clarity and sequence length.
sd_obs_5_theta <- rbind(data.frame(sd_cl_obs_5_theta), data.frame(sd_modcl_obs_5_theta),  
                        data.frame(sd_uncl_obs_5_theta))
sd_obs_5_theta$S_to_obs <- mapvalues(sd_obs_5_theta$S_to_obs, 
                                     from = c("sd_S1_cat1", "sd_S1_cat2", "sd_S1_cat3", "sd_S1_cat4", "sd_S1_cat5", 
                                              "sd_S2_cat1", "sd_S2_cat2", "sd_S2_cat3", "sd_S2_cat4", "sd_S2_cat5", 
                                              "sd_S3_cat1", "sd_S3_cat2", "sd_S3_cat3", "sd_S3_cat4", "sd_S3_cat5"), 
                                     to = c(rep(paste0("theta[", 1, "][", 1:5, "]"), each = 1), 
                                            rep(paste0("theta[", 2, "][", 1:5, "]"), each = 1), 
                                            rep(paste0("theta[", 3, "][", 1:5, "]"), each = 1)))
sd_obs_5_theta$Clarity <- factor(sd_obs_5_theta$Clarity, levels = c("Clear", "Moderate", "Unclear"))
ggplot(sd_obs_5_theta, aes(x = Length, y = Sd, color = Clarity, group = Clarity)) +
  facet_wrap(facets = vars(S_to_obs), nrow = 3, ncol = 5, labeller = label_parsed) +
  geom_point() + 
  geom_line() +
  xlab("Sequence length") + 
  ylab("Standard deviation") +
  geom_hline(yintercept = 0) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  ylim(NA, 0.20)

### Gamma.

## True transition probabilities.
gamma_true <- c(0.80, 0.10, 0.10, 
                0.10, 0.80, 0.10, 
                0.10, 0.10, 0.80)

for (i in 1: length(length_strings)){
  
  assign(paste0("bias_gamma_cl_obs_5_t_", length_strings[i]), 
         sweep(eval(parse(text = paste0("results$sim_HMM_theta_cl_obs_5_t_", length_strings[i], "$out_sim$gamma_mean"))), 
               2, gamma_true)) 
  
  assign(paste0("bias_gamma_cl_obs_5_t_", length_strings[i]), 
         melt(eval(parse(text = paste0("bias_gamma_cl_obs_5_t_", length_strings[i])))))
  
  assign(paste0("bias_gamma_cl_obs_5_t_", length_strings[i]), 
         cbind(eval(parse(text = paste0("bias_gamma_cl_obs_5_t_", length_strings[i]))), 
               Length = length_strings[i], Clarity = "Clear"))
  
  assign(paste0("bias_gamma_modcl_obs_5_t_", length_strings[i]), 
         sweep(eval(parse(text = paste0("results$sim_HMM_theta_modcl_obs_5_t_", length_strings[i], "$out_sim$gamma_mean"))), 
               2, gamma_true)) 
  
  assign(paste0("bias_gamma_modcl_obs_5_t_", length_strings[i]), 
         melt(eval(parse(text = paste0("bias_gamma_modcl_obs_5_t_", length_strings[i])))))
  
  assign(paste0("bias_gamma_modcl_obs_5_t_", length_strings[i]), 
         cbind(eval(parse(text = paste0("bias_gamma_modcl_obs_5_t_", length_strings[i]))), 
               Length = length_strings[i], Clarity = "Moderate"))
  
  assign(paste0("bias_gamma_uncl_obs_5_t_", length_strings[i]), 
         sweep(eval(parse(text = paste0("results$sim_HMM_theta_uncl_obs_5_t_", length_strings[i], "$out_sim$gamma_mean"))), 
               2, gamma_true)) 
  
  assign(paste0("bias_gamma_uncl_obs_5_t_", length_strings[i]), 
         melt(eval(parse(text = paste0("bias_gamma_uncl_obs_5_t_", length_strings[i])))))
  
  assign(paste0("bias_gamma_uncl_obs_5_t_", length_strings[i]), 
         cbind(eval(parse(text = paste0("bias_gamma_uncl_obs_5_t_", length_strings[i]))), 
               Length = length_strings[i], Clarity = "Unclear"))
  
}

# Clear. 
cl_obs_5_gamma <- rbind(data.frame(bias_gamma_cl_obs_5_t_250), data.frame(bias_gamma_cl_obs_5_t_500),  
                        data.frame(bias_gamma_cl_obs_5_t_1000), data.frame(bias_gamma_cl_obs_5_t_2000), 
                        data.frame(bias_gamma_cl_obs_5_t_4000), data.frame(bias_gamma_cl_obs_5_t_8000))
colnames(cl_obs_5_gamma) <- c("Id", "S_to_s", "Abs_bias", "Length", "Clarity")
cl_obs_5_gamma$Length <- factor(cl_obs_5_gamma$Length, levels = c("250", "500", "1000", "2000", "4000", "8000"))
cl_obs_5_gamma <- aggregate(Abs_bias ~ Length + S_to_s + Clarity, cl_obs_5_gamma, mean)

# Moderately clear. 
modcl_obs_5_gamma <- rbind(data.frame(bias_gamma_modcl_obs_5_t_250), data.frame(bias_gamma_modcl_obs_5_t_500),  
                           data.frame(bias_gamma_modcl_obs_5_t_1000), data.frame(bias_gamma_modcl_obs_5_t_2000), 
                           data.frame(bias_gamma_modcl_obs_5_t_4000), data.frame(bias_gamma_modcl_obs_5_t_8000))
colnames(modcl_obs_5_gamma) <- c("Id", "S_to_s", "Abs_bias", "Length", "Clarity")
modcl_obs_5_gamma$Length <- factor(modcl_obs_5_gamma$Length, levels = c("250", "500", "1000", "2000", "4000", "8000"))
modcl_obs_5_gamma <- aggregate(Abs_bias ~ Length + S_to_s + Clarity, modcl_obs_5_gamma, mean)

# Unclear. 
uncl_obs_5_gamma <- rbind(data.frame(bias_gamma_uncl_obs_5_t_250), data.frame(bias_gamma_uncl_obs_5_t_500),  
                          data.frame(bias_gamma_uncl_obs_5_t_1000), data.frame(bias_gamma_uncl_obs_5_t_2000), 
                          data.frame(bias_gamma_uncl_obs_5_t_4000), data.frame(bias_gamma_uncl_obs_5_t_8000))
colnames(uncl_obs_5_gamma) <- c("Id", "S_to_s", "Abs_bias", "Length", "Clarity")
uncl_obs_5_gamma$Length <- factor(uncl_obs_5_gamma$Length, levels = c("250", "500", "1000", "2000", "4000", "8000"))
uncl_obs_5_gamma <- aggregate(Abs_bias ~ Length + S_to_s + Clarity, uncl_obs_5_gamma, mean)

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
ggplot(obs_5_gamma, aes(x = Length, y = Abs_bias, color = Clarity, group = Clarity)) +
  facet_wrap(facets = vars(S_to_s), nrow = 3, ncol = 3, labeller = label_parsed) +
  geom_point() + geom_line() +  
  ggtitle("") +
  xlab("Sequence length") + 
  ylab("Bias") +
  geom_hline(yintercept = 0) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  ylim(-0.50, 0.25)

## Standard deviation.
for (i in 1: length(length_strings)){
  
  assign(paste0("sd_gamma_cl_obs_5_t_", length_strings[i]), 
         melt(eval(parse(text = paste0("results$sim_HMM_theta_cl_obs_5_t_", length_strings[i], "$out_sim$gamma_sd")))))
  
  assign(paste0("sd_gamma_cl_obs_5_t_", length_strings[i]), 
         cbind(eval(parse(text = paste0("sd_gamma_cl_obs_5_t_", length_strings[i]))), 
               Length = length_strings[i], Clarity = "Clear"))
  
  assign(paste0("sd_gamma_modcl_obs_5_t_", length_strings[i]), 
         melt(eval(parse(text = paste0("results$sim_HMM_theta_modcl_obs_5_t_", length_strings[i], "$out_sim$gamma_sd")))))
  
  assign(paste0("sd_gamma_modcl_obs_5_t_", length_strings[i]), 
         cbind(eval(parse(text = paste0("sd_gamma_modcl_obs_5_t_", length_strings[i]))), 
               Length = length_strings[i], Clarity = "Moderate"))
  
  assign(paste0("sd_gamma_uncl_obs_5_t_", length_strings[i]), 
         melt(eval(parse(text = paste0("results$sim_HMM_theta_uncl_obs_5_t_", length_strings[i], "$out_sim$gamma_sd")))))
  
  assign(paste0("sd_gamma_uncl_obs_5_t_", length_strings[i]), 
         cbind(eval(parse(text = paste0("sd_gamma_uncl_obs_5_t_", length_strings[i]))), 
               Length = length_strings[i], Clarity = "Unclear"))
  
}

# Clear.
sd_cl_obs_5_gamma <- rbind(data.frame(sd_gamma_cl_obs_5_t_250), data.frame(sd_gamma_cl_obs_5_t_500),  
                           data.frame(sd_gamma_cl_obs_5_t_1000), data.frame(sd_gamma_cl_obs_5_t_2000), 
                           data.frame(sd_gamma_cl_obs_5_t_4000), data.frame(sd_gamma_cl_obs_5_t_8000))
colnames(sd_cl_obs_5_gamma) <- c("Id", "S_to_s", "Sd", "Length", "Clarity")
sd_cl_obs_5_gamma$Length <- factor(sd_cl_obs_5_gamma$Length, levels = c("250", "500", "1000", "2000", "4000", "8000"))
sd_cl_obs_5_gamma <- aggregate(Sd ~ Length + S_to_s + Clarity, sd_cl_obs_5_gamma, mean)

# Moderately clear. 
sd_modcl_obs_5_gamma <- rbind(data.frame(sd_gamma_modcl_obs_5_t_250), data.frame(sd_gamma_modcl_obs_5_t_500),  
                              data.frame(sd_gamma_modcl_obs_5_t_1000), data.frame(sd_gamma_modcl_obs_5_t_2000), 
                              data.frame(sd_gamma_modcl_obs_5_t_4000), data.frame(sd_gamma_modcl_obs_5_t_8000))
colnames(sd_modcl_obs_5_gamma) <- c("Id", "S_to_s", "Sd", "Length", "Clarity")
sd_modcl_obs_5_gamma$Length <- factor(sd_modcl_obs_5_gamma$Length, levels = c("250", "500", "1000", "2000", "4000", "8000"))
sd_modcl_obs_5_gamma <- aggregate(Sd ~ Length + S_to_s + Clarity, sd_modcl_obs_5_gamma, mean)

# Unclear. 
sd_uncl_obs_5_gamma <- rbind(data.frame(sd_gamma_uncl_obs_5_t_250), data.frame(sd_gamma_uncl_obs_5_t_500),  
                             data.frame(sd_gamma_uncl_obs_5_t_1000), data.frame(sd_gamma_uncl_obs_5_t_2000), 
                             data.frame(sd_gamma_uncl_obs_5_t_4000), data.frame(sd_gamma_uncl_obs_5_t_8000))
colnames(sd_uncl_obs_5_gamma) <- c("Id", "S_to_s", "Sd", "Length", "Clarity")
sd_uncl_obs_5_gamma$Length <- factor(sd_uncl_obs_5_gamma$Length, levels = c("250", "500", "1000", "2000", "4000", "8000"))
sd_uncl_obs_5_gamma <- aggregate(Sd ~ Length + S_to_s + Clarity, sd_uncl_obs_5_gamma, mean)

# Plot of standard deviation transition probabilities by clarity and sequence length.
sd_obs_5_gamma <- rbind(data.frame(sd_cl_obs_5_gamma), data.frame(sd_modcl_obs_5_gamma),  
                        data.frame(sd_uncl_obs_5_gamma))
sd_obs_5_gamma$S_to_s <- mapvalues(sd_obs_5_gamma$S_to_s, 
                                   from = c("sd_S1_to_S1", "sd_S1_to_S2", "sd_S1_to_S3", 
                                            "sd_S2_to_S1", "sd_S2_to_S2", "sd_S2_to_S3",
                                            "sd_S3_to_S1", "sd_S3_to_S2", "sd_S3_to_S3"), 
                                   to = c(rep(paste0("gamma[", 1, "][", 1:3, "]"), each = 1), 
                                          rep(paste0("gamma[", 2, "][", 1:3, "]"), each = 1), 
                                          rep(paste0("gamma[", 3, "][", 1:3, "]"), each = 1)))
sd_obs_5_gamma$Clarity <- factor(sd_obs_5_gamma$Clarity, levels = c("Clear", "Moderate", "Unclear"))
ggplot(sd_obs_5_gamma, aes(x = Length, y = Sd, color = Clarity, group = Clarity)) +
  facet_wrap(facets = vars(S_to_s), nrow = 3, ncol = 3, labeller = label_parsed) +
  geom_point() + 
  geom_line() +  
  xlab("Sequence length") + 
  ylab("Standard deviation") +
  geom_hline(yintercept = 0) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  ylim(NA, 0.25)

#########################
### Overlap condition ###
#########################

### Theta.

## True emission probabilities.

# No overlap.
no_overlap_true <- c(0.84, 0.04, 0.04, 0.04, 0.04,
                     0.04, 0.44, 0.44, 0.04, 0.04,
                     0.04, 0.04, 0.04, 0.44, 0.44)  
# Moderate overlap.
mod_overlap_true <- c(0.59, 0.29, 0.04, 0.04, 0.04,
                      0.04, 0.29, 0.59, 0.04, 0.04,
                      0.04, 0.04, 0.20, 0.36, 0.36)
# Much overlap.
much_overlap_true <- c(0.44, 0.44, 0.04, 0.04, 0.04,
                       0.04, 0.44, 0.44, 0.04, 0.04,
                       0.04, 0.04, 0.30, 0.31, 0.31)

for (i in 1: length(length_strings)){
  
  assign(paste0("bias_theta_nooverl_obs_5_t_", length_strings[i]), 
         sweep(eval(parse(text = paste0("results$sim_HMM_theta_nooverl_obs_5_t_", length_strings[i], "$out_sim$emiss_mean"))), 
               2, no_overlap_true)) 
  
  assign(paste0("bias_theta_nooverl_obs_5_t_", length_strings[i]), 
         melt(eval(parse(text = paste0("bias_theta_nooverl_obs_5_t_", length_strings[i])))))
  
  assign(paste0("bias_theta_nooverl_obs_5_t_", length_strings[i]), 
         cbind(eval(parse(text = paste0("bias_theta_nooverl_obs_5_t_", length_strings[i]))), 
               Length = length_strings[i], Overlap = "None"))
  
  assign(paste0("bias_theta_modoverl_obs_5_t_", length_strings[i]), 
         sweep(eval(parse(text = paste0("results$sim_HMM_theta_modoverl_obs_5_t_", length_strings[i], "$out_sim$emiss_mean"))), 
               2, mod_overlap_true)) 
  
  assign(paste0("bias_theta_modoverl_obs_5_t_", length_strings[i]), 
         melt(eval(parse(text = paste0("bias_theta_modoverl_obs_5_t_", length_strings[i])))))
  
  assign(paste0("bias_theta_modoverl_obs_5_t_", length_strings[i]), 
         cbind(eval(parse(text = paste0("bias_theta_modoverl_obs_5_t_", length_strings[i]))), 
               Length = length_strings[i], Overlap = "Moderate"))
  
  assign(paste0("bias_theta_muchoverl_obs_5_t_", length_strings[i]), 
         sweep(eval(parse(text = paste0("results$sim_HMM_theta_muchoverl_obs_5_t_", length_strings[i], "$out_sim$emiss_mean"))), 
               2, much_overlap_true)) 
  
  assign(paste0("bias_theta_muchoverl_obs_5_t_", length_strings[i]), 
         melt(eval(parse(text = paste0("bias_theta_muchoverl_obs_5_t_", length_strings[i])))))
  
  assign(paste0("bias_theta_muchoverl_obs_5_t_", length_strings[i]), 
         cbind(eval(parse(text = paste0("bias_theta_muchoverl_obs_5_t_", length_strings[i]))), 
               Length = length_strings[i], Overlap = "Much"))
  
}

# No overlap.
nooverl_obs_5_theta <- rbind(data.frame(bias_theta_nooverl_obs_5_t_250), data.frame(bias_theta_nooverl_obs_5_t_500),  
                             data.frame(bias_theta_nooverl_obs_5_t_1000), data.frame(bias_theta_nooverl_obs_5_t_2000), 
                             data.frame(bias_theta_nooverl_obs_5_t_4000), data.frame(bias_theta_nooverl_obs_5_t_8000))
colnames(nooverl_obs_5_theta) <- c("Id", "S_to_obs", "Abs_bias", "Length", "Overlap")
nooverl_obs_5_theta$Length <- factor(nooverl_obs_5_theta$Length, levels = c("250", "500", "1000", "2000", "4000", "8000"))
nooverl_obs_5_theta <- aggregate(Abs_bias ~ Length + S_to_obs + Overlap, nooverl_obs_5_theta, mean)

# Moderate overlap. 
modoverl_obs_5_theta <- rbind(data.frame(bias_theta_modoverl_obs_5_t_250), data.frame(bias_theta_modoverl_obs_5_t_500),  
                              data.frame(bias_theta_modoverl_obs_5_t_1000), data.frame(bias_theta_modoverl_obs_5_t_2000), 
                              data.frame(bias_theta_modoverl_obs_5_t_4000), data.frame(bias_theta_modoverl_obs_5_t_8000))
colnames(modoverl_obs_5_theta) <- c("Id", "S_to_obs", "Abs_bias", "Length", "Overlap")
modoverl_obs_5_theta$Length <- factor(modoverl_obs_5_theta$Length, levels = c("250", "500", "1000", "2000", "4000", "8000"))
modoverl_obs_5_theta <- aggregate(Abs_bias ~ Length + S_to_obs + Overlap, modoverl_obs_5_theta, mean)

# Much overlap. 
muchoverl_obs_5_theta <- rbind(data.frame(bias_theta_muchoverl_obs_5_t_250), data.frame(bias_theta_muchoverl_obs_5_t_500),  
                               data.frame(bias_theta_muchoverl_obs_5_t_1000), data.frame(bias_theta_muchoverl_obs_5_t_2000), 
                               data.frame(bias_theta_muchoverl_obs_5_t_4000), data.frame(bias_theta_muchoverl_obs_5_t_8000))
colnames(muchoverl_obs_5_theta) <- c("Id", "S_to_obs", "Abs_bias", "Length", "Overlap")
muchoverl_obs_5_theta$Length <- factor(muchoverl_obs_5_theta$Length, levels = c("250", "500", "1000", "2000", "4000", "8000"))
muchoverl_obs_5_theta <- aggregate(Abs_bias ~ Length + S_to_obs + Overlap, muchoverl_obs_5_theta, mean)

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
ggplot(overl_5_theta, aes(x = Length, y = Abs_bias, color = Overlap, group = Overlap)) +
  facet_wrap(facets = vars(S_to_obs), nrow = 3, ncol = 5, labeller = label_parsed) +
  geom_point() + 
  geom_line() +
  xlab("Sequence length") + 
  ylab("Bias") +
  geom_hline(yintercept = 0) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  ylim(-0.4, 0.25) 

## Standard deviation. 
for (i in 1: length(length_strings)){
  
  assign(paste0("sd_theta_nooverl_obs_5_t_", length_strings[i]), 
         melt(eval(parse(text = paste0("results$sim_HMM_theta_nooverl_obs_5_t_", length_strings[i], "$out_sim$emiss_sd")))))
  
  assign(paste0("sd_theta_nooverl_obs_5_t_", length_strings[i]), 
         cbind(eval(parse(text = paste0("sd_theta_nooverl_obs_5_t_", length_strings[i]))), 
               Length = length_strings[i], Overlap = "None"))
  
  assign(paste0("sd_theta_modoverl_obs_5_t_", length_strings[i]), 
         melt(eval(parse(text = paste0("results$sim_HMM_theta_modoverl_obs_5_t_", length_strings[i], "$out_sim$emiss_sd")))))
  
  assign(paste0("sd_theta_modoverl_obs_5_t_", length_strings[i]), 
         cbind(eval(parse(text = paste0("sd_theta_modoverl_obs_5_t_", length_strings[i]))), 
               Length = length_strings[i], Overlap = "Moderate"))
  
  assign(paste0("sd_theta_muchoverl_obs_5_t_", length_strings[i]), 
         melt(eval(parse(text = paste0("results$sim_HMM_theta_muchoverl_obs_5_t_", length_strings[i], "$out_sim$emiss_sd")))))
  
  assign(paste0("sd_theta_muchoverl_obs_5_t_", length_strings[i]), 
         cbind(eval(parse(text = paste0("sd_theta_muchoverl_obs_5_t_", length_strings[i]))), 
               Length = length_strings[i], Overlap = "Much"))
  
}

# No overlap.
sd_nooverl_obs_5_theta <- rbind(data.frame(sd_theta_nooverl_obs_5_t_250), data.frame(sd_theta_nooverl_obs_5_t_500),  
                                data.frame(sd_theta_nooverl_obs_5_t_1000), data.frame(sd_theta_nooverl_obs_5_t_2000), 
                                data.frame(sd_theta_nooverl_obs_5_t_4000), data.frame(sd_theta_nooverl_obs_5_t_8000))
colnames(sd_nooverl_obs_5_theta) <- c("Id", "S_to_obs", "Sd", "Length", "Overlap")
sd_nooverl_obs_5_theta$Length <- factor(sd_nooverl_obs_5_theta$Length, levels = c("250", "500", "1000", "2000", "4000", "8000"))
sd_nooverl_obs_5_theta$Overlap <- factor(sd_nooverl_obs_5_theta$Overlap, levels = c("None", "Moderate", "Much"))
sd_nooverl_obs_5_theta <- aggregate(Sd ~ Length + S_to_obs + Overlap, sd_nooverl_obs_5_theta, mean)

# Moderate overlap. 
sd_modoverl_obs_5_theta <- rbind(data.frame(sd_theta_modoverl_obs_5_t_250), data.frame(sd_theta_modoverl_obs_5_t_500),  
                                 data.frame(sd_theta_modoverl_obs_5_t_1000), data.frame(sd_theta_modoverl_obs_5_t_2000), 
                                 data.frame(sd_theta_modoverl_obs_5_t_4000), data.frame(sd_theta_modoverl_obs_5_t_8000))
colnames(sd_modoverl_obs_5_theta) <- c("Id", "S_to_obs", "Sd", "Length", "Overlap")
sd_modoverl_obs_5_theta$Length <- factor(sd_modoverl_obs_5_theta$Length, levels = c("250", "500", "1000", "2000", "4000", "8000"))
sd_modoverl_obs_5_theta <- aggregate(Sd ~ Length + S_to_obs + Overlap, sd_modoverl_obs_5_theta, mean)

# Much overlap. 
sd_muchoverl_obs_5_theta <- rbind(data.frame(sd_theta_muchoverl_obs_5_t_250), data.frame(sd_theta_muchoverl_obs_5_t_500), 
                                  data.frame(sd_theta_muchoverl_obs_5_t_1000), data.frame(sd_theta_muchoverl_obs_5_t_2000), 
                                  data.frame(sd_theta_muchoverl_obs_5_t_4000), data.frame(sd_theta_muchoverl_obs_5_t_8000))
colnames(sd_muchoverl_obs_5_theta) <- c("Id", "S_to_obs", "Sd", "Length", "Overlap")
sd_muchoverl_obs_5_theta$Length <- factor(sd_muchoverl_obs_5_theta$Length, levels = c("250", "500", "1000", "2000", "4000", "8000"))
sd_muchoverl_obs_5_theta <- aggregate(Sd ~ Length + S_to_obs + Overlap, sd_muchoverl_obs_5_theta, mean)

# Plot of standard deviation emission probabilities by overlap and sequence length
sd_obs_5_theta_overl <- rbind(data.frame(sd_nooverl_obs_5_theta), data.frame(sd_modoverl_obs_5_theta),  
                              data.frame(sd_muchoverl_obs_5_theta))
sd_obs_5_theta_overl$S_to_obs <- mapvalues(sd_obs_5_theta_overl$S_to_obs, 
                                           from = c("sd_S1_cat1", "sd_S1_cat2", "sd_S1_cat3", "sd_S1_cat4", "sd_S1_cat5", 
                                                    "sd_S2_cat1", "sd_S2_cat2", "sd_S2_cat3", "sd_S2_cat4", "sd_S2_cat5", 
                                                    "sd_S3_cat1", "sd_S3_cat2", "sd_S3_cat3", "sd_S3_cat4", "sd_S3_cat5"), 
                                           to = c(rep(paste0("theta[", 1, "][", 1:5, "]"), each = 1), 
                                                  rep(paste0("theta[", 2, "][", 1:5, "]"), each = 1), 
                                                  rep(paste0("theta[", 3, "][", 1:5, "]"), each = 1)))
sd_obs_5_theta_overl$Overlap <- factor(sd_obs_5_theta_overl$Overlap, levels = c("None", "Moderate", "Much"))
ggplot(sd_obs_5_theta_overl, aes(x = Length, y = Sd, color = Overlap, group = Overlap)) +
  facet_wrap(facets = vars(S_to_obs), nrow = 3, ncol = 5, labeller = label_parsed) +
  geom_point() + 
  geom_line() +  
  xlab("Sequence length") + 
  ylab("Standard deviation") +
  geom_hline(yintercept = 0) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  ylim(NA, 0.20)

### Gamma.
for (i in 1: length(length_strings)){
  
  assign(paste0("bias_gamma_nooverl_obs_5_t_", length_strings[i]), 
         sweep(eval(parse(text = paste0("results$sim_HMM_theta_nooverl_obs_5_t_", length_strings[i], "$out_sim$gamma_mean"))), 
               2, gamma_true)) 
  
  assign(paste0("bias_gamma_nooverl_obs_5_t_", length_strings[i]), 
         melt(eval(parse(text = paste0("bias_gamma_nooverl_obs_5_t_", length_strings[i])))))
  
  assign(paste0("bias_gamma_nooverl_obs_5_t_", length_strings[i]), 
         cbind(eval(parse(text = paste0("bias_gamma_nooverl_obs_5_t_", length_strings[i]))), 
               Length = length_strings[i], Overlap = "None"))
  
  assign(paste0("bias_gamma_modoverl_obs_5_t_", length_strings[i]), 
         sweep(eval(parse(text = paste0("results$sim_HMM_theta_modoverl_obs_5_t_", length_strings[i], "$out_sim$gamma_mean"))), 
               2, gamma_true)) 
  
  assign(paste0("bias_gamma_modoverl_obs_5_t_", length_strings[i]), 
         melt(eval(parse(text = paste0("bias_gamma_modoverl_obs_5_t_", length_strings[i])))))
  
  assign(paste0("bias_gamma_modoverl_obs_5_t_", length_strings[i]), 
         cbind(eval(parse(text = paste0("bias_gamma_modoverl_obs_5_t_", length_strings[i]))), 
               Length = length_strings[i], Overlap = "Moderate"))
  
  assign(paste0("bias_gamma_muchoverl_obs_5_t_", length_strings[i]), 
         sweep(eval(parse(text = paste0("results$sim_HMM_theta_muchoverl_obs_5_t_", length_strings[i], "$out_sim$gamma_mean"))), 
               2, gamma_true)) 
  
  assign(paste0("bias_gamma_muchoverl_obs_5_t_", length_strings[i]), 
         melt(eval(parse(text = paste0("bias_gamma_muchoverl_obs_5_t_", length_strings[i])))))
  
  assign(paste0("bias_gamma_muchoverl_obs_5_t_", length_strings[i]), 
         cbind(eval(parse(text = paste0("bias_gamma_muchoverl_obs_5_t_", length_strings[i]))), 
               Length = length_strings[i], Overlap = "Much"))
  
}

# No overlap.
nooverl_obs_5_gamma <- rbind(data.frame(bias_gamma_nooverl_obs_5_t_250), data.frame(bias_gamma_nooverl_obs_5_t_500),  
                             data.frame(bias_gamma_nooverl_obs_5_t_1000), data.frame(bias_gamma_nooverl_obs_5_t_2000), 
                             data.frame(bias_gamma_nooverl_obs_5_t_4000), data.frame(bias_gamma_nooverl_obs_5_t_8000))
colnames(nooverl_obs_5_gamma) <- c("Id", "S_to_s", "Abs_bias", "Length", "Overlap")
nooverl_obs_5_gamma$Length <- factor(nooverl_obs_5_gamma$Length, levels = c("250", "500", "1000", "2000", "4000", "8000"))
nooverl_obs_5_gamma <- aggregate(Abs_bias ~ Length + S_to_s + Overlap, nooverl_obs_5_gamma, mean)

# Moderate overlap. 
modoverl_obs_5_gamma <- rbind(data.frame(bias_gamma_modoverl_obs_5_t_250), data.frame(bias_gamma_modoverl_obs_5_t_500),  
                              data.frame(bias_gamma_modoverl_obs_5_t_1000), data.frame(bias_gamma_modoverl_obs_5_t_2000), 
                              data.frame(bias_gamma_modoverl_obs_5_t_4000), data.frame(bias_gamma_modoverl_obs_5_t_8000))
colnames(modoverl_obs_5_gamma) <- c("Id", "S_to_s", "Abs_bias", "Length", "Overlap")
modoverl_obs_5_gamma$Length <- factor(modoverl_obs_5_gamma$Length, levels = c("250", "500", "1000", "2000", "4000", "8000"))
modoverl_obs_5_gamma <- aggregate(Abs_bias ~ Length + S_to_s + Overlap, modoverl_obs_5_gamma, mean)

# Much overlap. 
muchoverl_obs_5_gamma <- rbind(data.frame(bias_gamma_muchoverl_obs_5_t_250), data.frame(bias_gamma_muchoverl_obs_5_t_500),  
                               data.frame(bias_gamma_muchoverl_obs_5_t_1000), data.frame(bias_gamma_muchoverl_obs_5_t_2000), 
                               data.frame(bias_gamma_muchoverl_obs_5_t_4000), data.frame(bias_gamma_muchoverl_obs_5_t_8000))
colnames(muchoverl_obs_5_gamma) <- c("Id", "S_to_s", "Abs_bias", "Length", "Overlap")
muchoverl_obs_5_gamma$Length <- factor(muchoverl_obs_5_gamma$Length, levels = c("250", "500", "1000", "2000", "4000", "8000"))
muchoverl_obs_5_gamma <- aggregate(Abs_bias ~ Length + S_to_s + Overlap, muchoverl_obs_5_gamma, mean)

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
ggplot(overl_5_gamma, aes(x = Length, y = Abs_bias, color = Overlap, group = Overlap)) +
  facet_wrap(facets = vars(S_to_s), nrow = 3, ncol = 3, labeller = label_parsed) +
  geom_point() + 
  geom_line() +  
  xlab("Sequence length") + 
  ylab("Bias") +
  geom_hline(yintercept = 0) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  ylim(-0.50, 0.25)

## Standard deviation overlap condition gamma.
for (i in 1: length(length_strings)){
  
  assign(paste0("sd_gamma_nooverl_obs_5_t_", length_strings[i]), 
         melt(eval(parse(text = paste0("results$sim_HMM_theta_nooverl_obs_5_t_", length_strings[i], "$out_sim$gamma_sd")))))
  
  assign(paste0("sd_gamma_nooverl_obs_5_t_", length_strings[i]), 
         cbind(eval(parse(text = paste0("sd_gamma_nooverl_obs_5_t_", length_strings[i]))), 
               Length = length_strings[i], Overlap = "None"))
  
  assign(paste0("sd_gamma_modoverl_obs_5_t_", length_strings[i]), 
         melt(eval(parse(text = paste0("results$sim_HMM_theta_modoverl_obs_5_t_", length_strings[i], "$out_sim$gamma_sd")))))
  
  assign(paste0("sd_gamma_modoverl_obs_5_t_", length_strings[i]), 
         cbind(eval(parse(text = paste0("sd_gamma_modoverl_obs_5_t_", length_strings[i]))), 
               Length = length_strings[i], Overlap = "Moderate"))
  
  assign(paste0("sd_gamma_muchoverl_obs_5_t_", length_strings[i]), 
         melt(eval(parse(text = paste0("results$sim_HMM_theta_muchoverl_obs_5_t_", length_strings[i], "$out_sim$gamma_sd")))))
  
  assign(paste0("sd_gamma_muchoverl_obs_5_t_", length_strings[i]), 
         cbind(eval(parse(text = paste0("sd_gamma_muchoverl_obs_5_t_", length_strings[i]))), 
               Length = length_strings[i], Overlap = "Much"))
  
}

# No overlap.
sd_nooverl_obs_5_gamma <- rbind(data.frame(sd_gamma_nooverl_obs_5_t_250), data.frame(sd_gamma_nooverl_obs_5_t_500),  
                                data.frame(sd_gamma_nooverl_obs_5_t_1000), data.frame(sd_gamma_nooverl_obs_5_t_2000), 
                                data.frame(sd_gamma_nooverl_obs_5_t_4000), data.frame(sd_gamma_nooverl_obs_5_t_8000))
colnames(sd_nooverl_obs_5_gamma) <- c("Id", "S_to_s", "Sd", "Length", "Overlap")
sd_nooverl_obs_5_gamma$Length <- factor(sd_nooverl_obs_5_gamma$Length, levels = c("250", "500", "1000", "2000", "4000", "8000"))
sd_nooverl_obs_5_gamma <- aggregate(Sd ~ Length + S_to_s + Overlap, sd_nooverl_obs_5_gamma, mean)

# Moderate overlap.
sd_modoverl_obs_5_gamma <- rbind(data.frame(sd_gamma_modoverl_obs_5_t_250), data.frame(sd_gamma_modoverl_obs_5_t_500),  
                                 data.frame(sd_gamma_modoverl_obs_5_t_1000), data.frame(sd_gamma_modoverl_obs_5_t_2000), 
                                 data.frame(sd_gamma_modoverl_obs_5_t_4000), data.frame(sd_gamma_modoverl_obs_5_t_8000))
colnames(sd_modoverl_obs_5_gamma) <- c("Id", "S_to_s", "Sd", "Length", "Overlap")
sd_modoverl_obs_5_gamma$Length <- factor(sd_modoverl_obs_5_gamma$Length, levels = c("250", "500", "1000", "2000", "4000", "8000"))
sd_modoverl_obs_5_gamma <- aggregate(Sd ~ Length + S_to_s + Overlap, sd_modoverl_obs_5_gamma, mean)

# Much overlap. 
sd_muchoverl_obs_5_gamma <- rbind(data.frame(sd_gamma_muchoverl_obs_5_t_250), data.frame(sd_gamma_muchoverl_obs_5_t_500),  
                                  data.frame(sd_gamma_muchoverl_obs_5_t_1000), data.frame(sd_gamma_muchoverl_obs_5_t_2000), 
                                  data.frame(sd_gamma_muchoverl_obs_5_t_4000), data.frame(sd_gamma_muchoverl_obs_5_t_8000))
colnames(sd_muchoverl_obs_5_gamma) <- c("Id", "S_to_s", "Sd", "Length", "Overlap")
sd_muchoverl_obs_5_gamma$Length <- factor(sd_muchoverl_obs_5_gamma$Length, levels = c("250", "500", "1000", "2000", "4000", "8000"))
sd_muchoverl_obs_5_gamma <- aggregate(Sd ~ Length + S_to_s + Overlap, sd_muchoverl_obs_5_gamma, mean)

# Plot of standard deviation transition probabilities by overlap and sequence length.
sd_obs_5_gamma_overl <- rbind(data.frame(sd_nooverl_obs_5_gamma), data.frame(sd_modoverl_obs_5_gamma),  
                              data.frame(sd_muchoverl_obs_5_gamma))
sd_obs_5_gamma_overl$S_to_s <- mapvalues(sd_obs_5_gamma_overl$S_to_s, 
                                         from = c("sd_S1_to_S1", "sd_S1_to_S2", "sd_S1_to_S3", 
                                                  "sd_S2_to_S1", "sd_S2_to_S2", "sd_S2_to_S3",
                                                  "sd_S3_to_S1", "sd_S3_to_S2", "sd_S3_to_S3"), 
                                         to = c(rep(paste0("gamma[", 1, "][", 1:3, "]"), each = 1), 
                                                rep(paste0("gamma[", 2, "][", 1:3, "]"), each = 1), 
                                                rep(paste0("gamma[", 3, "][", 1:3, "]"), each = 1)))
sd_obs_5_gamma_overl$Overlap <- factor(sd_obs_5_gamma_overl$Overlap, levels = c("None", "Moderate", "Much"))
ggplot(sd_obs_5_gamma_overl, aes(x = Length, y = Sd, color = Overlap, group = Overlap)) +
  facet_wrap(facets = vars(S_to_s), nrow = 3, ncol = 3, labeller = label_parsed) +
  geom_point() + 
  geom_line() + 
  xlab("Sequence length") + 
  ylab("Standard deviation") +
  geom_hline(yintercept = 0) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  ylim(NA, 0.25)

##############################
### Number of observations ###
##############################

### Theta.

## True emission probabilities.

# Three observations.
true_theta_three_obs <- c(0.52, 0.24, 0.24,
                          0.24, 0.52, 0.24,
                          0.24, 0.24, 0.52)
# Five observations.
true_theta_five_obs <- c(0.44, 0.14, 0.14, 0.14, 0.14,
                         0.14, 0.29, 0.14, 0.29, 0.14,
                         0.14, 0.14, 0.29, 0.14, 0.29)

# Seven observations.
true_theta_seven_obs <- c(0.40, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10,
                          0.10, 0.25, 0.10, 0.25, 0.10, 0.10, 0.10,
                          0.08, 0.08, 0.22, 0.09, 0.22, 0.09, 0.22)

for (i in 1: length(length_strings)){
  
  assign(paste0("bias_theta_uncl_obs_3_t_", length_strings[i]), 
         sweep(eval(parse(text = paste0("results$sim_HMM_theta_uncl_obs_3_t_", length_strings[i], "$out_sim$emiss_mean"))), 
               2, true_theta_three_obs)) 
  
  assign(paste0("bias_theta_uncl_obs_3_t_", length_strings[i]), 
         melt(eval(parse(text = paste0("bias_theta_uncl_obs_3_t_", length_strings[i])))))
  
  assign(paste0("bias_theta_uncl_obs_3_t_", length_strings[i]), 
         cbind(eval(parse(text = paste0("bias_theta_uncl_obs_3_t_", length_strings[i]))), 
               Length = length_strings[i], Nobs = "3"))
  
  assign(paste0("bias_theta_uncl_obs_5_t_", length_strings[i]), 
         sweep(eval(parse(text = paste0("results$sim_HMM_theta_uncl_obs_5_t_", length_strings[i], "$out_sim$emiss_mean"))), 
               2, true_theta_five_obs)) 
  
  assign(paste0("bias_theta_uncl_obs_5_t_", length_strings[i]), 
         melt(eval(parse(text = paste0("bias_theta_uncl_obs_5_t_", length_strings[i])))))
  
  assign(paste0("bias_theta_uncl_obs_5_t_", length_strings[i]), 
         cbind(eval(parse(text = paste0("bias_theta_uncl_obs_5_t_", length_strings[i]))), 
               Length = length_strings[i], Nobs = "5"))
  
  assign(paste0("bias_theta_uncl_obs_7_t_", length_strings[i]), 
         sweep(eval(parse(text = paste0("results$sim_HMM_theta_uncl_obs_7_t_", length_strings[i], "$out_sim$emiss_mean"))), 
               2, true_theta_seven_obs)) 
  
  assign(paste0("bias_theta_uncl_obs_7_t_", length_strings[i]), 
         melt(eval(parse(text = paste0("bias_theta_uncl_obs_7_t_", length_strings[i])))))
  
  assign(paste0("bias_theta_uncl_obs_7_t_", length_strings[i]), 
         cbind(eval(parse(text = paste0("bias_theta_uncl_obs_7_t_", length_strings[i]))), 
               Length = length_strings[i], Nobs = "7"))
  
}

### Three observations.
uncl_obs_3_theta <- rbind(data.frame(bias_theta_uncl_obs_3_t_250), data.frame(bias_theta_uncl_obs_3_t_500),  
                          data.frame(bias_theta_uncl_obs_3_t_1000), data.frame(bias_theta_uncl_obs_3_t_2000), 
                          data.frame(bias_theta_uncl_obs_3_t_4000), data.frame(bias_theta_uncl_obs_3_t_8000))
colnames(uncl_obs_3_theta) <- c("Id", "S_to_obs", "Abs_bias", "Length", "Nobs")
uncl_obs_3_theta$Length <- factor(uncl_obs_3_theta$Length, levels = c("250", "500", "1000", "2000", "4000", "8000"))
uncl_obs_3_theta <- aggregate(Abs_bias ~ Length + S_to_obs + Nobs, uncl_obs_3_theta, mean)

### Five observations.
uncl_obs_5_theta <- rbind(data.frame(bias_theta_uncl_obs_5_t_250), data.frame(bias_theta_uncl_obs_5_t_500),  
                          data.frame(bias_theta_uncl_obs_5_t_1000), data.frame(bias_theta_uncl_obs_5_t_2000), 
                          data.frame(bias_theta_uncl_obs_5_t_4000), data.frame(bias_theta_uncl_obs_5_t_8000))
colnames(uncl_obs_5_theta) <- c("Id", "S_to_obs", "Abs_bias", "Length", "Nobs")
uncl_obs_5_theta$Length <- factor(uncl_obs_5_theta$Length, levels = c("250", "500", "1000", "2000", "4000", "8000"))
uncl_obs_5_theta <- aggregate(Abs_bias ~ Length + S_to_obs + Nobs, uncl_obs_5_theta, mean)

### Seven observations.
uncl_obs_7_theta <- rbind(data.frame(bias_theta_uncl_obs_7_t_250), data.frame(bias_theta_uncl_obs_7_t_500),  
                          data.frame(bias_theta_uncl_obs_7_t_1000), data.frame(bias_theta_uncl_obs_7_t_2000), 
                          data.frame(bias_theta_uncl_obs_7_t_4000), data.frame(bias_theta_uncl_obs_7_t_8000))
colnames(uncl_obs_7_theta) <- c("Id", "S_to_obs", "Abs_bias", "Length", "Nobs")
uncl_obs_7_theta$Length <- factor(uncl_obs_7_theta$Length, levels = c("250", "500", "1000", "2000", "4000", "8000"))
uncl_obs_7_theta$Nobs <- factor(uncl_obs_7_theta$Nobs, levels = c("3", "5", "7"))
uncl_obs_7_theta <- aggregate(Abs_bias ~ Length + S_to_obs + Nobs, uncl_obs_7_theta, mean)

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
ggplot(uncl_obs_theta, aes(x = Length, y = Abs_bias, color = Nobs, group = Nobs)) +
  facet_wrap(facets = vars(S_to_obs), nrow = 3, ncol = 7, labeller = label_parsed) +
  geom_point() + geom_line() +  
  xlab("Sequence length") + 
  ylab("Bias") +
  geom_hline(yintercept = 0) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  ylim(-0.4, 0.25) 

## Standard deviation.
for (i in 1: length(length_strings)){
  
  assign(paste0("sd_theta_uncl_obs_3_t_", length_strings[i]), 
         melt(eval(parse(text = paste0("results$sim_HMM_theta_uncl_obs_3_t_", length_strings[i], "$out_sim$emiss_sd")))))
  
  assign(paste0("sd_theta_uncl_obs_3_t_", length_strings[i]), 
         cbind(eval(parse(text = paste0("sd_theta_uncl_obs_3_t_", length_strings[i]))), 
               Length = length_strings[i], Nobs = "3"))
  
  assign(paste0("sd_theta_uncl_obs_5_t_", length_strings[i]), 
         melt(eval(parse(text = paste0("results$sim_HMM_theta_uncl_obs_5_t_", length_strings[i], "$out_sim$emiss_sd")))))
  
  assign(paste0("sd_theta_uncl_obs_5_t_", length_strings[i]), 
         cbind(eval(parse(text = paste0("sd_theta_uncl_obs_5_t_", length_strings[i]))), 
               Length = length_strings[i], Nobs = "5"))
  
  assign(paste0("sd_theta_uncl_obs_7_t_", length_strings[i]), 
         melt(eval(parse(text = paste0("results$sim_HMM_theta_uncl_obs_7_t_", length_strings[i], "$out_sim$emiss_sd")))))
  
  assign(paste0("sd_theta_uncl_obs_7_t_", length_strings[i]), 
         cbind(eval(parse(text = paste0("sd_theta_uncl_obs_7_t_", length_strings[i]))), 
               Length = length_strings[i], Nobs = "7"))
  
}

# 3 observations.
sd_uncl_obs_3_theta <- rbind(data.frame(sd_theta_uncl_obs_3_t_250), data.frame(sd_theta_uncl_obs_3_t_500),  
                             data.frame(sd_theta_uncl_obs_3_t_1000), data.frame(sd_theta_uncl_obs_3_t_2000), 
                             data.frame(sd_theta_uncl_obs_3_t_4000), data.frame(sd_theta_uncl_obs_3_t_8000))
colnames(sd_uncl_obs_3_theta) <- c("Id", "S_to_obs", "Sd", "Length", "Nobs")
sd_uncl_obs_3_theta$Length <- factor(sd_uncl_obs_3_theta$Length, levels = c("250", "500", "1000", "2000", "4000", "8000"))
sd_uncl_obs_3_theta <- aggregate(Sd ~ Length + S_to_obs + Nobs, sd_uncl_obs_3_theta, mean)

# 5 observations. 
sd_uncl_obs_5_theta <- rbind(data.frame(sd_theta_uncl_obs_5_t_250), data.frame(sd_theta_uncl_obs_5_t_500),  
                             data.frame(sd_theta_uncl_obs_5_t_1000), data.frame(sd_theta_uncl_obs_5_t_2000), 
                             data.frame(sd_theta_uncl_obs_5_t_4000), data.frame(sd_theta_uncl_obs_5_t_8000))
colnames(sd_uncl_obs_5_theta) <- c("Id", "S_to_obs", "Sd", "Length", "Nobs")
sd_uncl_obs_5_theta$Length <- factor(sd_uncl_obs_5_theta$Length, levels = c("250", "500", "1000", "2000", "4000", "8000"))
sd_uncl_obs_5_theta <- aggregate(Sd ~ Length + S_to_obs + Nobs, sd_uncl_obs_5_theta, mean)

# 7 observations.  
sd_uncl_obs_7_theta <- rbind(data.frame(sd_theta_uncl_obs_7_t_250), data.frame(sd_theta_uncl_obs_7_t_500),  
                             data.frame(sd_theta_uncl_obs_7_t_1000), data.frame(sd_theta_uncl_obs_7_t_2000), 
                             data.frame(sd_theta_uncl_obs_7_t_4000), data.frame(sd_theta_uncl_obs_7_t_8000))
colnames(sd_uncl_obs_7_theta) <- c("Id", "S_to_obs", "Sd", "Length", "Nobs")
sd_uncl_obs_7_theta$Length <- factor(sd_uncl_obs_7_theta$Length, levels = c("250", "500", "1000", "2000", "4000", "8000"))
sd_uncl_obs_7_theta <- aggregate(Sd ~ Length + S_to_obs + Nobs, sd_uncl_obs_7_theta, mean)

# Plot of standard deviation emission probabilities by number of observations and sequence length.
sd_obs_theta <- rbind(data.frame(sd_uncl_obs_3_theta), data.frame(sd_uncl_obs_5_theta),  
                      data.frame(sd_uncl_obs_7_theta))
sd_obs_theta$S_to_obs <- mapvalues(sd_obs_theta$S_to_obs, 
                                   from = c("sd_S1_cat1", "sd_S1_cat2", "sd_S1_cat3", "sd_S1_cat4", "sd_S1_cat5", "sd_S1_cat6", "sd_S1_cat7",
                                            "sd_S2_cat1", "sd_S2_cat2", "sd_S2_cat3", "sd_S2_cat4", "sd_S2_cat5", "sd_S2_cat6", "sd_S2_cat7", 
                                            "sd_S3_cat1", "sd_S3_cat2", "sd_S3_cat3", "sd_S3_cat4", "sd_S3_cat5", "sd_S3_cat6", "sd_S3_cat7"), 
                                   to = c(rep(paste0("theta[", 1, "][", 1:7, "]"), each = 1), 
                                          rep(paste0("theta[", 2, "][", 1:7, "]"), each = 1), 
                                          rep(paste0("theta[", 3, "][", 1:7, "]"), each = 1)))
sd_obs_theta$S_to_obs <- factor(sd_obs_theta$S_to_obs, levels = c("theta[1][1]", "theta[1][2]", "theta[1][3]", "theta[1][4]", 
                                                                  "theta[1][5]", "theta[1][6]", "theta[1][7]",
                                                                  "theta[2][1]", "theta[2][2]", "theta[2][3]", "theta[2][4]", 
                                                                  "theta[2][5]", "theta[2][6]", "theta[2][7]",
                                                                  "theta[3][1]", "theta[3][2]", "theta[3][3]", "theta[3][4]", 
                                                                  "theta[3][5]", "theta[3][6]", "theta[3][7]"))
sd_obs_theta$Nobs <- factor(sd_obs_theta$Nobs, levels = c("3", "5", "7"))
ggplot(sd_obs_theta, aes(x = Length, y = Sd, color = Nobs, group = Nobs)) +
  facet_wrap(facets = vars(S_to_obs), nrow = 3, ncol = 7, labeller = label_parsed) +
  geom_point() + geom_line() +
  xlab("Sequence length") + 
  ylab("Standard deviation") +
  geom_hline(yintercept = 0) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  ylim(NA, 0.20)

### Gamma. 
for (i in 1: length(length_strings)){
  
  assign(paste0("bias_gamma_uncl_obs_3_t_", length_strings[i]), 
         sweep(eval(parse(text = paste0("results$sim_HMM_theta_uncl_obs_3_t_", length_strings[i], "$out_sim$gamma_mean"))), 
               2, gamma_true)) 
  
  assign(paste0("bias_gamma_uncl_obs_3_t_", length_strings[i]), 
         melt(eval(parse(text = paste0("bias_gamma_uncl_obs_3_t_", length_strings[i])))))
  
  assign(paste0("bias_gamma_uncl_obs_3_t_", length_strings[i]), 
         cbind(eval(parse(text = paste0("bias_gamma_uncl_obs_3_t_", length_strings[i]))), 
               Length = length_strings[i], Nobs = "3"))
  
  assign(paste0("bias_gamma_uncl_obs_5_t_", length_strings[i]), 
         sweep(eval(parse(text = paste0("results$sim_HMM_theta_uncl_obs_5_t_", length_strings[i], "$out_sim$gamma_mean"))), 
               2, gamma_true)) 
  
  assign(paste0("bias_gamma_uncl_obs_5_t_", length_strings[i]), 
         melt(eval(parse(text = paste0("bias_gamma_uncl_obs_5_t_", length_strings[i])))))
  
  assign(paste0("bias_gamma_uncl_obs_5_t_", length_strings[i]), 
         cbind(eval(parse(text = paste0("bias_gamma_uncl_obs_5_t_", length_strings[i]))), 
               Length = length_strings[i], Nobs = "5"))
  
  assign(paste0("bias_gamma_uncl_obs_7_t_", length_strings[i]), 
         sweep(eval(parse(text = paste0("results$sim_HMM_theta_uncl_obs_7_t_", length_strings[i], "$out_sim$gamma_mean"))), 
               2, gamma_true)) 
  
  assign(paste0("bias_gamma_uncl_obs_7_t_", length_strings[i]), 
         melt(eval(parse(text = paste0("bias_gamma_uncl_obs_7_t_", length_strings[i])))))
  
  assign(paste0("bias_gamma_uncl_obs_7_t_", length_strings[i]), 
         cbind(eval(parse(text = paste0("bias_gamma_uncl_obs_7_t_", length_strings[i]))), 
               Length = length_strings[i], Nobs = "7"))
  
}

# 3 observations.
uncl_obs_3_gamma <- rbind(data.frame(bias_gamma_uncl_obs_3_t_250), data.frame(bias_gamma_uncl_obs_3_t_500),  
                          data.frame(bias_gamma_uncl_obs_3_t_1000), data.frame(bias_gamma_uncl_obs_3_t_2000), 
                          data.frame(bias_gamma_uncl_obs_3_t_4000), data.frame(bias_gamma_uncl_obs_3_t_8000))
colnames(uncl_obs_3_gamma) <- c("Id", "S_to_s", "Abs_bias", "Length", "Nobs")
uncl_obs_3_gamma$Length <- factor(uncl_obs_3_gamma$Length, levels = c("250", "500", "1000", "2000", "4000", "8000"))
uncl_obs_3_gamma <- aggregate(Abs_bias ~ Length + S_to_s + Nobs, uncl_obs_3_gamma, mean)

# 5 observations. 
uncl_obs_5_gamma <- rbind(data.frame(bias_gamma_uncl_obs_5_t_250), data.frame(bias_gamma_uncl_obs_5_t_500),  
                          data.frame(bias_gamma_uncl_obs_5_t_1000), data.frame(bias_gamma_uncl_obs_5_t_2000), 
                          data.frame(bias_gamma_uncl_obs_5_t_4000), data.frame(bias_gamma_uncl_obs_5_t_8000))
colnames(uncl_obs_5_gamma) <- c("Id", "S_to_s", "Abs_bias", "Length", "Nobs")
uncl_obs_5_gamma$Length <- factor(uncl_obs_5_gamma$Length, levels = c("250", "500", "1000", "2000", "4000", "8000"))
uncl_obs_5_gamma <- aggregate(Abs_bias ~ Length + S_to_s + Nobs, uncl_obs_5_gamma, mean)

# 7 observations. 
uncl_obs_7_gamma <- rbind(data.frame(bias_gamma_uncl_obs_7_t_250), data.frame(bias_gamma_uncl_obs_7_t_500),  
                          data.frame(bias_gamma_uncl_obs_7_t_1000), data.frame(bias_gamma_uncl_obs_7_t_2000), 
                          data.frame(bias_gamma_uncl_obs_7_t_4000), data.frame(bias_gamma_uncl_obs_7_t_8000))
colnames(uncl_obs_7_gamma) <- c("Id", "S_to_s", "Abs_bias", "Length", "Nobs")
uncl_obs_7_gamma$Length <- factor(uncl_obs_7_gamma$Length, levels = c("250", "500", "1000", "2000", "4000", "8000"))
uncl_obs_7_gamma <- aggregate(Abs_bias ~ Length + S_to_s + Nobs, uncl_obs_7_gamma, mean)

# Plot of bias transition probabilities by number of observations and sequence length,
num_obs_gamma <- rbind(data.frame(uncl_obs_3_gamma), data.frame(uncl_obs_5_gamma),
                       data.frame(uncl_obs_7_gamma))
num_obs_gamma$S_to_s <- mapvalues(num_obs_gamma$S_to_s, 
                                  from = c("M_S1_to_S1", "M_S1_to_S2", "M_S1_to_S3", 
                                           "M_S2_to_S1", "M_S2_to_S2", "M_S2_to_S3",
                                           "M_S3_to_S1", "M_S3_to_S2", "M_S3_to_S3"), 
                                  to = c(rep(paste0("gamma[", 1, "][", 1:3, "]"), each = 1), 
                                         rep(paste0("gamma[", 2, "][", 1:3, "]"), each = 1), 
                                         rep(paste0("gamma[", 3, "][", 1:3, "]"), each = 1)))
num_obs_gamma$Nobs <- factor(num_obs_gamma$Nobs, levels = c("3", "5", "7"))
ggplot(num_obs_gamma, aes(x = Length, y = Abs_bias, color = Nobs, group = Nobs)) +
  facet_wrap(facets = vars(S_to_s), nrow = 3, ncol = 3, labeller = label_parsed) +
  geom_point() + geom_line() +  
  xlab("Sequence length") + 
  ylab("Bias") +
  geom_hline(yintercept = 0) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  ylim(-0.50, 0.25)

## Standard deviation.
for (i in 1: length(length_strings)){
  
  assign(paste0("sd_gamma_uncl_obs_3_t_", length_strings[i]), 
         melt(eval(parse(text = paste0("results$sim_HMM_theta_uncl_obs_3_t_", length_strings[i], "$out_sim$gamma_sd")))))
  
  assign(paste0("sd_gamma_uncl_obs_3_t_", length_strings[i]), 
         cbind(eval(parse(text = paste0("sd_gamma_uncl_obs_3_t_", length_strings[i]))), 
               Length = length_strings[i], Nobs = "3"))
  
  assign(paste0("sd_gamma_uncl_obs_5_t_", length_strings[i]), 
         melt(eval(parse(text = paste0("results$sim_HMM_theta_uncl_obs_5_t_", length_strings[i], "$out_sim$gamma_sd")))))
  
  assign(paste0("sd_gamma_uncl_obs_5_t_", length_strings[i]), 
         cbind(eval(parse(text = paste0("sd_gamma_uncl_obs_5_t_", length_strings[i]))), 
               Length = length_strings[i], Nobs = "5"))
  
  assign(paste0("sd_gamma_uncl_obs_7_t_", length_strings[i]), 
         melt(eval(parse(text = paste0("results$sim_HMM_theta_uncl_obs_7_t_", length_strings[i], "$out_sim$gamma_sd")))))
  
  assign(paste0("sd_gamma_uncl_obs_7_t_", length_strings[i]), 
         cbind(eval(parse(text = paste0("sd_gamma_uncl_obs_7_t_", length_strings[i]))), 
               Length = length_strings[i], Nobs = "7"))
  
}

# 3 observations.
sd_uncl_obs_3_gamma <- rbind(data.frame(sd_gamma_uncl_obs_3_t_250), data.frame(sd_gamma_uncl_obs_3_t_500),  
                             data.frame(sd_gamma_uncl_obs_3_t_1000), data.frame(sd_gamma_uncl_obs_3_t_2000), 
                             data.frame(sd_gamma_uncl_obs_3_t_4000), data.frame(sd_gamma_uncl_obs_3_t_8000))
colnames(sd_uncl_obs_3_gamma) <- c("Id", "S_to_s", "Sd", "Length", "Nobs")
sd_uncl_obs_3_gamma$Length <- factor(sd_uncl_obs_3_gamma$Length, levels = c("250", "500", "1000", "2000", "4000", "8000"))
sd_uncl_obs_3_gamma <- aggregate(Sd ~ Length + S_to_s + Nobs, sd_uncl_obs_3_gamma, mean)

# 5 observations.
sd_uncl_obs_5_gamma <- rbind(data.frame(sd_gamma_uncl_obs_5_t_250), data.frame(sd_gamma_uncl_obs_5_t_500),  
                             data.frame(sd_gamma_uncl_obs_5_t_1000), data.frame(sd_gamma_uncl_obs_5_t_2000), 
                             data.frame(sd_gamma_uncl_obs_5_t_4000), data.frame(sd_gamma_uncl_obs_5_t_8000))
colnames(sd_uncl_obs_5_gamma) <- c("Id", "S_to_s", "Sd", "Length", "Nobs")
sd_uncl_obs_5_gamma$Length <- factor(sd_uncl_obs_5_gamma$Length, levels = c("250", "500", "1000", "2000", "4000", "8000"))
sd_uncl_obs_5_gamma <- aggregate(Sd ~ Length + S_to_s + Nobs, sd_uncl_obs_5_gamma, mean)

# 7 observations. 
sd_uncl_obs_7_gamma <- rbind(data.frame(sd_gamma_uncl_obs_7_t_250), data.frame(sd_gamma_uncl_obs_7_t_500),  
                             data.frame(sd_gamma_uncl_obs_7_t_1000), data.frame(sd_gamma_uncl_obs_7_t_2000), 
                             data.frame(sd_gamma_uncl_obs_7_t_4000), data.frame(sd_gamma_uncl_obs_7_t_8000))
colnames(sd_uncl_obs_7_gamma) <- c("Id", "S_to_s", "Sd", "Length", "Nobs")
sd_uncl_obs_7_gamma$Length <- factor(sd_uncl_obs_7_gamma$Length, levels = c("250", "500", "1000", "2000", "4000", "8000"))
sd_uncl_obs_7_gamma <- aggregate(Sd ~ Length + S_to_s + Nobs, sd_uncl_obs_7_gamma, mean)

# Plot of standard deviation transition probabilities by number of observations and sequence length.
sd_num_obs_gamma <- rbind(data.frame(sd_uncl_obs_3_gamma), data.frame(sd_uncl_obs_5_gamma),  
                          data.frame(sd_uncl_obs_7_gamma))
sd_num_obs_gamma$S_to_s <- mapvalues(sd_num_obs_gamma$S_to_s, 
                                     from = c("sd_S1_to_S1", "sd_S1_to_S2", "sd_S1_to_S3", 
                                              "sd_S2_to_S1", "sd_S2_to_S2", "sd_S2_to_S3",
                                              "sd_S3_to_S1", "sd_S3_to_S2", "sd_S3_to_S3"), 
                                     to = c(rep(paste0("gamma[", 1, "][", 1:3, "]"), each = 1), 
                                            rep(paste0("gamma[", 2, "][", 1:3, "]"), each = 1), 
                                            rep(paste0("gamma[", 3, "][", 1:3, "]"), each = 1)))
sd_num_obs_gamma$Nobs <- factor(sd_num_obs_gamma$Nobs, levels = c("3", "5", "7"))
ggplot(sd_num_obs_gamma, aes(x = Length, y = Sd, color = Nobs, group = Nobs)) +
  facet_wrap(facets = vars(S_to_s), nrow = 3, ncol = 3, labeller = label_parsed) +
  geom_point() + geom_line() +  
  xlab("Sequence length") + 
  ylab("Standard deviation") +
  geom_hline(yintercept = 0) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  ylim(NA, 0.25)

##############
### Varobs ###
##############

### Theta.

## True emission probabilities.

# Three observations.
true_theta_three_obs <- c(0.80, 0.10, 0.10,
                          0.10, 0.80, 0.10,
                          0.10, 0.10, 0.80)
# Five observations.
true_theta_five_obs <- c(0.92, 0.02, 0.02, 0.02, 0.02,
                         0.02, 0.47, 0.02, 0.47, 0.02,
                         0.02, 0.02, 0.62, 0.02, 0.32)

# Seven observations.
true_theta_seven_obs <- c(0.76, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04,
                          0.02, 0.45, 0.02, 0.45, 0.02, 0.02, 0.02,
                          0.01, 0.02, 0.31, 0.02, 0.31, 0.02, 0.31)

for (i in 1: length(length_strings)){
  
  assign(paste0("bias_theta_var_obs_3_t_", length_strings[i]), 
         sweep(eval(parse(text = paste0("results$sim_HMM_varobs_obs_3_t_", length_strings[i], "$out_sim$emiss_mean"))), 
               2, true_theta_three_obs)) 
  
  assign(paste0("bias_theta_var_obs_3_t_", length_strings[i]), 
         melt(eval(parse(text = paste0("bias_theta_var_obs_3_t_", length_strings[i])))))
  
  assign(paste0("bias_theta_var_obs_3_t_", length_strings[i]), 
         cbind(eval(parse(text = paste0("bias_theta_var_obs_3_t_", length_strings[i]))), 
               Length = length_strings[i], Nobs = "3"))
  
  assign(paste0("bias_theta_var_obs_5_t_", length_strings[i]), 
         sweep(eval(parse(text = paste0("results$sim_HMM_varobs_obs_5_t_", length_strings[i], "$out_sim$emiss_mean"))), 
               2, true_theta_five_obs)) 
  
  assign(paste0("bias_theta_var_obs_5_t_", length_strings[i]), 
         melt(eval(parse(text = paste0("bias_theta_var_obs_5_t_", length_strings[i])))))
  
  assign(paste0("bias_theta_var_obs_5_t_", length_strings[i]), 
         cbind(eval(parse(text = paste0("bias_theta_var_obs_5_t_", length_strings[i]))), 
               Length = length_strings[i], Nobs = "5"))
  
  assign(paste0("bias_theta_var_obs_7_t_", length_strings[i]), 
         sweep(eval(parse(text = paste0("results$sim_HMM_varobs_obs_7_t_", length_strings[i], "$out_sim$emiss_mean"))), 
               2, true_theta_seven_obs)) 
  
  assign(paste0("bias_theta_var_obs_7_t_", length_strings[i]), 
         melt(eval(parse(text = paste0("bias_theta_var_obs_7_t_", length_strings[i])))))
  
  assign(paste0("bias_theta_var_obs_7_t_", length_strings[i]), 
         cbind(eval(parse(text = paste0("bias_theta_var_obs_7_t_", length_strings[i]))), 
               Length = length_strings[i], Nobs = "7"))
  
}

### Three observations.
var_obs_3_theta <- rbind(data.frame(bias_theta_var_obs_3_t_250), data.frame(bias_theta_var_obs_3_t_500),  
                         data.frame(bias_theta_var_obs_3_t_1000), data.frame(bias_theta_var_obs_3_t_2000), 
                         data.frame(bias_theta_var_obs_3_t_4000), data.frame(bias_theta_var_obs_3_t_8000))
colnames(var_obs_3_theta) <- c("Id", "S_to_obs", "Abs_bias", "Length", "Nobs")
var_obs_3_theta$Length <- factor(var_obs_3_theta$Length, levels = c("250", "500", "1000", "2000", "4000", "8000"))
var_obs_3_theta <- aggregate(Abs_bias ~ Length + S_to_obs + Nobs, var_obs_3_theta, mean)

### Five observations.
var_obs_5_theta <- rbind(data.frame(bias_theta_var_obs_5_t_250), data.frame(bias_theta_var_obs_5_t_500),  
                         data.frame(bias_theta_var_obs_5_t_1000), data.frame(bias_theta_var_obs_5_t_2000), 
                         data.frame(bias_theta_var_obs_5_t_4000), data.frame(bias_theta_var_obs_5_t_8000))
colnames(var_obs_5_theta) <- c("Id", "S_to_obs", "Abs_bias", "Length", "Nobs")
var_obs_5_theta$Length <- factor(var_obs_5_theta$Length, levels = c("250", "500", "1000", "2000", "4000", "8000"))
var_obs_5_theta <- aggregate(Abs_bias ~ Length + S_to_obs + Nobs, var_obs_5_theta, mean)

### Seven observations.
var_obs_7_theta <- rbind(data.frame(bias_theta_var_obs_7_t_250), data.frame(bias_theta_var_obs_7_t_500),  
                         data.frame(bias_theta_var_obs_7_t_1000), data.frame(bias_theta_var_obs_7_t_2000), 
                         data.frame(bias_theta_var_obs_7_t_4000), data.frame(bias_theta_var_obs_7_t_8000))
colnames(var_obs_7_theta) <- c("Id", "S_to_obs", "Abs_bias", "Length", "Nobs")
var_obs_7_theta$Length <- factor(var_obs_7_theta$Length, levels = c("250", "500", "1000", "2000", "4000", "8000"))
var_obs_7_theta <- aggregate(Abs_bias ~ Length + S_to_obs + Nobs, var_obs_7_theta, mean)

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
ggplot(var_obs_theta, aes(x = Length, y = Abs_bias, color = Nobs, group = Nobs)) +
  facet_wrap(facets = vars(S_to_obs), nrow = 3, ncol = 7, labeller = label_parsed) +
  geom_point() + geom_line() +  
  xlab("Sequence length") + 
  ylab("Bias") +
  geom_hline(yintercept = 0) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  ylim(-0.4, 0.25) 

## Standard deviation.
for (i in 1: length(length_strings)){
  
  assign(paste0("sd_varobs_obs_3_t_", length_strings[i]), 
         melt(eval(parse(text = paste0("results$sim_HMM_varobs_obs_3_t_", length_strings[i], "$out_sim$emiss_sd")))))
  
  assign(paste0("sd_varobs_obs_3_t_", length_strings[i]), 
         cbind(eval(parse(text = paste0("sd_varobs_obs_3_t_", length_strings[i]))), 
               Length = length_strings[i], Nobs = "3"))
  
  assign(paste0("sd_varobs_obs_5_t_", length_strings[i]), 
         melt(eval(parse(text = paste0("results$sim_HMM_varobs_obs_5_t_", length_strings[i], "$out_sim$emiss_sd")))))
  
  assign(paste0("sd_varobs_obs_5_t_", length_strings[i]), 
         cbind(eval(parse(text = paste0("sd_varobs_obs_5_t_", length_strings[i]))), 
               Length = length_strings[i], Nobs = "5"))
  
  assign(paste0("sd_varobs_obs_7_t_", length_strings[i]), 
         melt(eval(parse(text = paste0("results$sim_HMM_varobs_obs_7_t_", length_strings[i], "$out_sim$emiss_sd")))))
  
  assign(paste0("sd_varobs_obs_7_t_", length_strings[i]), 
         cbind(eval(parse(text = paste0("sd_varobs_obs_7_t_", length_strings[i]))), 
               Length = length_strings[i], Nobs = "7"))
  
}

# 3 observations.
sd_var_obs_3_theta <- rbind(data.frame(sd_varobs_obs_3_t_250), data.frame(sd_varobs_obs_3_t_500),  
                            data.frame(sd_varobs_obs_3_t_1000), data.frame(sd_varobs_obs_3_t_2000), 
                            data.frame(sd_varobs_obs_3_t_4000), data.frame(sd_varobs_obs_3_t_8000))
colnames(sd_var_obs_3_theta) <- c("Id", "S_to_obs", "Sd", "Length", "Nobs")
sd_var_obs_3_theta$Length <- factor(sd_var_obs_3_theta$Length, levels = c("250", "500", "1000", "2000", "4000", "8000"))
sd_var_obs_3_theta <- aggregate(Sd ~ Length + S_to_obs + Nobs, sd_var_obs_3_theta, mean)

# 5 observations. 
sd_var_obs_5_theta <- rbind(data.frame(sd_varobs_obs_5_t_250), data.frame(sd_varobs_obs_5_t_500),  
                            data.frame(sd_varobs_obs_5_t_1000), data.frame(sd_varobs_obs_5_t_2000), 
                            data.frame(sd_varobs_obs_5_t_4000), data.frame(sd_varobs_obs_5_t_8000))
colnames(sd_var_obs_5_theta) <- c("Id", "S_to_obs", "Sd", "Length", "Nobs")
sd_var_obs_5_theta$Length <- factor(sd_var_obs_5_theta$Length, levels = c("250", "500", "1000", "2000", "4000", "8000"))
sd_var_obs_5_theta <- aggregate(Sd ~ Length + S_to_obs + Nobs, sd_var_obs_5_theta, mean)

# 7 observations.  
sd_var_obs_7_theta <- rbind(data.frame(sd_varobs_obs_7_t_250), data.frame(sd_varobs_obs_7_t_500),  
                            data.frame(sd_varobs_obs_7_t_1000), data.frame(sd_varobs_obs_7_t_2000), 
                            data.frame(sd_varobs_obs_7_t_4000), data.frame(sd_varobs_obs_7_t_8000))
colnames(sd_var_obs_7_theta) <- c("Id", "S_to_obs", "Sd", "Length", "Nobs")
sd_var_obs_7_theta$Length <- factor(sd_var_obs_7_theta$Length, levels = c("250", "500", "1000", "2000", "4000", "8000"))
sd_var_obs_7_theta <- aggregate(Sd ~ Length + S_to_obs + Nobs, sd_var_obs_7_theta, mean)

# Plot of standard deviation emission probabilities by number of observations and sequence length.
sd_obs_theta <- rbind(data.frame(sd_var_obs_3_theta), data.frame(sd_var_obs_5_theta),  
                      data.frame(sd_var_obs_7_theta))
sd_obs_theta$S_to_obs <- mapvalues(sd_obs_theta$S_to_obs, 
                                   from = c("sd_S1_cat1", "sd_S1_cat2", "sd_S1_cat3", "sd_S1_cat4", "sd_S1_cat5", "sd_S1_cat6", "sd_S1_cat7",
                                            "sd_S2_cat1", "sd_S2_cat2", "sd_S2_cat3", "sd_S2_cat4", "sd_S2_cat5", "sd_S2_cat6", "sd_S2_cat7", 
                                            "sd_S3_cat1", "sd_S3_cat2", "sd_S3_cat3", "sd_S3_cat4", "sd_S3_cat5", "sd_S3_cat6", "sd_S3_cat7"), 
                                   to = c(rep(paste0("theta[", 1, "][", 1:7, "]"), each = 1), 
                                          rep(paste0("theta[", 2, "][", 1:7, "]"), each = 1), 
                                          rep(paste0("theta[", 3, "][", 1:7, "]"), each = 1)))
sd_obs_theta$S_to_obs <- factor(sd_obs_theta$S_to_obs, levels = c("theta[1][1]", "theta[1][2]", "theta[1][3]", "theta[1][4]", 
                                                                  "theta[1][5]", "theta[1][6]", "theta[1][7]",
                                                                  "theta[2][1]", "theta[2][2]", "theta[2][3]", "theta[2][4]", 
                                                                  "theta[2][5]", "theta[2][6]", "theta[2][7]",
                                                                  "theta[3][1]", "theta[3][2]", "theta[3][3]", "theta[3][4]", 
                                                                  "theta[3][5]", "theta[3][6]", "theta[3][7]"))
sd_obs_theta$Nobs <- factor(sd_obs_theta$Nobs, levels = c("3", "5", "7"))
ggplot(sd_obs_theta, aes(x = Length, y = Sd, color = Nobs, group = Nobs)) +
  facet_wrap(facets = vars(S_to_obs), nrow = 3, ncol = 7, labeller = label_parsed) +
  geom_point() + geom_line() +  
  xlab("Sequence length") + 
  ylab("Standard deviation") +
  geom_hline(yintercept = 0) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  ylim(NA, 0.20)

### Gamma. 
for (i in 1: length(length_strings)){
  
  assign(paste0("bias_gamma_var_obs_3_t_", length_strings[i]), 
         sweep(eval(parse(text = paste0("results$sim_HMM_varobs_obs_3_t_", length_strings[i], "$out_sim$gamma_mean"))), 
               2, gamma_true)) 
  
  assign(paste0("bias_gamma_var_obs_3_t_", length_strings[i]), 
         melt(eval(parse(text = paste0("bias_gamma_var_obs_3_t_", length_strings[i])))))
  
  assign(paste0("bias_gamma_var_obs_3_t_", length_strings[i]), 
         cbind(eval(parse(text = paste0("bias_gamma_var_obs_3_t_", length_strings[i]))), 
               Length = length_strings[i], Nobs = "3"))
  
  assign(paste0("bias_gamma_var_obs_5_t_", length_strings[i]), 
         sweep(eval(parse(text = paste0("results$sim_HMM_varobs_obs_5_t_", length_strings[i], "$out_sim$gamma_mean"))), 
               2, gamma_true)) 
  
  assign(paste0("bias_gamma_var_obs_5_t_", length_strings[i]), 
         melt(eval(parse(text = paste0("bias_gamma_var_obs_5_t_", length_strings[i])))))
  
  assign(paste0("bias_gamma_var_obs_5_t_", length_strings[i]), 
         cbind(eval(parse(text = paste0("bias_gamma_var_obs_5_t_", length_strings[i]))), 
               Length = length_strings[i], Nobs = "5"))
  
  assign(paste0("bias_gamma_var_obs_7_t_", length_strings[i]), 
         sweep(eval(parse(text = paste0("results$sim_HMM_varobs_obs_7_t_", length_strings[i], "$out_sim$gamma_mean"))), 
               2, gamma_true)) 
  
  assign(paste0("bias_gamma_var_obs_7_t_", length_strings[i]), 
         melt(eval(parse(text = paste0("bias_gamma_var_obs_7_t_", length_strings[i])))))
  
  assign(paste0("bias_gamma_var_obs_7_t_", length_strings[i]), 
         cbind(eval(parse(text = paste0("bias_gamma_var_obs_7_t_", length_strings[i]))), 
               Length = length_strings[i], Nobs = "7"))
  
}

# 3 observations.
var_obs_3_gamma <- rbind(data.frame(bias_gamma_var_obs_3_t_250), data.frame(bias_gamma_var_obs_3_t_500),  
                         data.frame(bias_gamma_var_obs_3_t_1000), data.frame(bias_gamma_var_obs_3_t_2000), 
                         data.frame(bias_gamma_var_obs_3_t_4000), data.frame(bias_gamma_var_obs_3_t_8000))
colnames(var_obs_3_gamma) <- c("Id", "S_to_s", "Abs_bias", "Length", "Nobs")
var_obs_3_gamma$Length <- factor(var_obs_3_gamma$Length, levels = c("250", "500", "1000", "2000", "4000", "8000"))
var_obs_3_gamma <- aggregate(Abs_bias ~ Length + S_to_s + Nobs, var_obs_3_gamma, mean)

# 5 observations. 
var_obs_5_gamma <- rbind(data.frame(bias_gamma_var_obs_5_t_250), data.frame(bias_gamma_var_obs_5_t_500),  
                         data.frame(bias_gamma_var_obs_5_t_1000), data.frame(bias_gamma_var_obs_5_t_2000), 
                         data.frame(bias_gamma_var_obs_5_t_4000), data.frame(bias_gamma_var_obs_5_t_8000))
colnames(var_obs_5_gamma) <- c("Id", "S_to_s", "Abs_bias", "Length", "Nobs")
var_obs_5_gamma$Length <- factor(var_obs_5_gamma$Length, levels = c("250", "500", "1000", "2000", "4000", "8000"))
var_obs_5_gamma <- aggregate(Abs_bias ~ Length + S_to_s + Nobs, var_obs_5_gamma, mean)

# 7 observations. 
var_obs_7_gamma <- rbind(data.frame(bias_gamma_var_obs_7_t_250), data.frame(bias_gamma_var_obs_7_t_500),  
                         data.frame(bias_gamma_var_obs_7_t_1000), data.frame(bias_gamma_var_obs_7_t_2000), 
                         data.frame(bias_gamma_var_obs_7_t_4000), data.frame(bias_gamma_var_obs_7_t_8000))
colnames(var_obs_7_gamma) <- c("Id", "S_to_s", "Abs_bias", "Length", "Nobs")
var_obs_7_gamma$Length <- factor(var_obs_7_gamma$Length, levels = c("250", "500", "1000", "2000", "4000", "8000"))
var_obs_7_gamma <- aggregate(Abs_bias ~ Length + S_to_s + Nobs, var_obs_7_gamma, mean)

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
ggplot(num_obs_gamma, aes(x = Length, y = Abs_bias, color = Nobs, group = Nobs)) +
  facet_wrap(facets = vars(S_to_s), nrow = 3, ncol = 3, labeller = label_parsed) +
  geom_point() + geom_line() +  
  xlab("Sequence length") + 
  ylab("Bias") +
  geom_hline(yintercept = 0) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  ylim(-0.50, 0.25)

## Standard deviation.
for (i in 1: length(length_strings)){
  
  assign(paste0("sd_gamma_var_obs_3_t_", length_strings[i]), 
         melt(eval(parse(text = paste0("results$sim_HMM_varobs_obs_3_t_", length_strings[i], "$out_sim$gamma_sd")))))
  
  assign(paste0("sd_gamma_var_obs_3_t_", length_strings[i]), 
         cbind(eval(parse(text = paste0("sd_gamma_var_obs_3_t_", length_strings[i]))), 
               Length = length_strings[i], Nobs = "3"))
  
  assign(paste0("sd_gamma_var_obs_5_t_", length_strings[i]), 
         melt(eval(parse(text = paste0("results$sim_HMM_varobs_obs_5_t_", length_strings[i], "$out_sim$gamma_sd")))))
  
  assign(paste0("sd_gamma_var_obs_5_t_", length_strings[i]), 
         cbind(eval(parse(text = paste0("sd_gamma_var_obs_5_t_", length_strings[i]))), 
               Length = length_strings[i], Nobs = "5"))
  
  assign(paste0("sd_gamma_var_obs_7_t_", length_strings[i]), 
         melt(eval(parse(text = paste0("results$sim_HMM_varobs_obs_7_t_", length_strings[i], "$out_sim$gamma_sd")))))
  
  assign(paste0("sd_gamma_var_obs_7_t_", length_strings[i]), 
         cbind(eval(parse(text = paste0("sd_gamma_var_obs_7_t_", length_strings[i]))), 
               Length = length_strings[i], Nobs = "7"))
  
}

# 3 observations.
sd_var_obs_3_gamma <- rbind(data.frame(sd_gamma_var_obs_3_t_250), data.frame(sd_gamma_var_obs_3_t_500),  
                            data.frame(sd_gamma_var_obs_3_t_1000), data.frame(sd_gamma_var_obs_3_t_2000), 
                            data.frame(sd_gamma_var_obs_3_t_4000), data.frame(sd_gamma_var_obs_3_t_8000))
colnames(sd_var_obs_3_gamma) <- c("Id", "S_to_s", "Sd", "Length", "Nobs")
sd_var_obs_3_gamma$Length <- factor(sd_var_obs_3_gamma$Length, levels = c("250", "500", "1000", "2000", "4000", "8000"))
sd_var_obs_3_gamma <- aggregate(Sd ~ Length + S_to_s + Nobs, sd_var_obs_3_gamma, mean)

# 5 observations.
sd_var_obs_5_gamma <- rbind(data.frame(sd_gamma_var_obs_5_t_250), data.frame(sd_gamma_var_obs_5_t_500),  
                            data.frame(sd_gamma_var_obs_5_t_1000), data.frame(sd_gamma_var_obs_5_t_2000), 
                            data.frame(sd_gamma_var_obs_5_t_4000), data.frame(sd_gamma_var_obs_5_t_8000))
colnames(sd_var_obs_5_gamma) <- c("Id", "S_to_s", "Sd", "Length", "Nobs")
sd_var_obs_5_gamma$Length <- factor(sd_var_obs_5_gamma$Length, levels = c("250", "500", "1000", "2000", "4000", "8000"))
sd_var_obs_5_gamma <- aggregate(Sd ~ Length + S_to_s + Nobs, sd_var_obs_5_gamma, mean)

# 7 observations. 
sd_var_obs_7_gamma <- rbind(data.frame(sd_gamma_var_obs_7_t_250), data.frame(sd_gamma_var_obs_7_t_500),  
                            data.frame(sd_gamma_var_obs_7_t_1000), data.frame(sd_gamma_var_obs_7_t_2000), 
                            data.frame(sd_gamma_var_obs_7_t_4000), data.frame(sd_gamma_var_obs_7_t_8000))
colnames(sd_var_obs_7_gamma) <- c("Id", "S_to_s", "Sd", "Length", "Nobs")
sd_var_obs_7_gamma$Length <- factor(sd_var_obs_7_gamma$Length, levels = c("250", "500", "1000", "2000", "4000", "8000"))
sd_var_obs_7_gamma <- aggregate(Sd ~ Length + S_to_s + Nobs, sd_var_obs_7_gamma, mean)

# Plot of standard deviation transition probabilities by number of observations and sequence length.
sd_num_obs_gamma <- rbind(data.frame(sd_var_obs_3_gamma), data.frame(sd_var_obs_5_gamma),  
                          data.frame(sd_var_obs_7_gamma))
sd_num_obs_gamma$S_to_s <- mapvalues(sd_num_obs_gamma$S_to_s, 
                                     from = c("sd_S1_to_S1", "sd_S1_to_S2", "sd_S1_to_S3", 
                                              "sd_S2_to_S1", "sd_S2_to_S2", "sd_S2_to_S3",
                                              "sd_S3_to_S1", "sd_S3_to_S2", "sd_S3_to_S3"), 
                                     to = c(rep(paste0("gamma[", 1, "][", 1:3, "]"), each = 1), 
                                            rep(paste0("gamma[", 2, "][", 1:3, "]"), each = 1), 
                                            rep(paste0("gamma[", 3, "][", 1:3, "]"), each = 1)))
sd_num_obs_gamma$Nobs <- factor(sd_num_obs_gamma$Nobs, levels = c("3", "5", "7"))
ggplot(sd_num_obs_gamma, aes(x = Length, y = Sd, color = Nobs, group = Nobs)) +
  facet_wrap(facets = vars(S_to_s), nrow = 3, ncol = 3, labeller = label_parsed) +
  geom_point() + geom_line() +  
  xlab("Sequence length") + 
  ylab("Standard deviation") +
  geom_hline(yintercept = 0) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  ylim(NA, 0.25)

### Save plots to dir. 
plots.dir.path <- list.files(tempdir(), pattern="rs-graphics", full.names = TRUE); 
plots.png.paths <- list.files(plots.dir.path, pattern=".png", full.names = TRUE)
file.copy(from=plots.png.paths, to="C:/Academia/HMM paper/Plots")

plots.png.detials <- file.info(plots.png.paths)
plots.png.detials <- plots.png.detials[order(plots.png.detials$mtime),]
sorted.png.names <- gsub(plots.dir.path, "C:/Academia/HMM paper/Plots", row.names(plots.png.detials), fixed=TRUE)
numbered.png.names <- paste0("C:/Academia/HMM paper/Plots/", 1:length(sorted.png.names), ".png")

# Rename all the .png files as: 1.png, 2.png, 3.png, and so on.
file.rename(from=sorted.png.names, to=numbered.png.names)
