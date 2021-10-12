library(mHMMbayes)

sim_type <- "theta_cl" 

# simulation runs, iterations and burn_in:
n_sim <- 500
J <- 5000
burn_in <- 1000

# Set n subjects with each n_t observations:
n <- 1
n_t <- 1000

# (baseline) model properties
m <- 3
q_emiss <- 5
n_dep <- 1
gamma <-  matrix(c(0.80, 0.10, 0.10, 
                   0.10, 0.80, 0.10, 
                   0.10, 0.10, 0.80), byrow = TRUE, nrow  = m)

theta.clear <- matrix(c(0.92, 0.02, 0.02, 0.02, 0.02,
                        0.02, 0.47, 0.02, 0.47, 0.02,
                        0.02, 0.02, 0.47, 0.02, 0.47), byrow = TRUE, nrow = m)

# simulate data
data1 <- rep(list(NULL), n_sim)
set.seed(2145)

for(sim_i in 1:n_sim){
  data1[[sim_i]] <- sim_mHMM(n_t = n_t, n = n, m = m, q_emiss = q_emiss,
                             gamma = gamma, emiss_distr = list(theta.clear)) 
}

## settings for analysis
# Starting values
st_gamma <- matrix(c(0.70, 0.10, 0.20, 
                     0.10, 0.70, 0.20, 
                     0.10, 0.20, 0.70), byrow = TRUE, nrow  = m) 

st_emiss_distr <-  matrix(c(0.60, 0.10, 0.10, 0.10, 0.10,
                            0.10, 0.35, 0.10, 0.35, 0.10,
                            0.10, 0.10, 0.35, 0.10, 0.35), byrow = TRUE, nrow = m)

emiss_med             <- matrix(, nrow = n_sim, ncol = m * q_emiss)
colnames(emiss_med)   <- paste("S", rep(c(1:m), each = q_emiss), "_cat", rep(1:q_emiss, m), sep = "")
emiss_low             <- matrix(, nrow = n_sim, ncol = m * q_emiss)
colnames(emiss_low)   <- paste("L_S", rep(c(1:m), each = q_emiss), "_cat", rep(1:q_emiss, m), sep = "")
emiss_up              <- matrix(, nrow = n_sim, ncol = m * q_emiss)
colnames(emiss_up)    <- paste("U_S", rep(c(1:m), each = q_emiss), "_cat", rep(1:q_emiss, m), sep = "")
emiss_IQR             <- matrix(, nrow = n_sim, ncol = m * q_emiss)
colnames(emiss_IQR)   <- paste("IQR_S", rep(c(1:m), each = q_emiss), "_cat", rep(1:q_emiss, m), sep = "")
emiss_mean            <- matrix(, nrow = n_sim, ncol = m * q_emiss)
colnames(emiss_mean)  <- paste("M_S", rep(c(1:m), each = q_emiss), "_cat", rep(1:q_emiss, m), sep = "")
emiss_sd              <- matrix(, nrow = n_sim, ncol = m * q_emiss)
colnames(emiss_sd)    <- paste("sd_S", rep(c(1:m), each = q_emiss), "_cat", rep(1:q_emiss, m), sep = "")

gamma_med             <- matrix(, nrow = n_sim, ncol = m * m)
colnames(gamma_med)   <- paste("S", rep(c(1:m), each = m), "_to_S", rep(1:m, m), sep = "")
gamma_low             <- matrix(, nrow = n_sim, ncol = m * m)
colnames(gamma_low)   <- paste("L_S", rep(c(1:m), each = m), "_to_S", rep(1:m, m), sep = "")
gamma_up              <- matrix(, nrow = n_sim, ncol = m * m)
colnames(gamma_up)    <- paste("U_S", rep(c(1:m), each = m), "_to_S", rep(1:m, m), sep = "")
gamma_IQR             <- matrix(, nrow = n_sim, ncol = m * m)
colnames(gamma_IQR)   <- paste("IQR_S", rep(c(1:m), each = m), "_to_S", rep(1:m, m), sep = "")
gamma_mean            <- matrix(, nrow = n_sim, ncol = m * m)
colnames(gamma_mean)  <- paste("M_S", rep(c(1:m), each = m), "_to_S", rep(1:m, m), sep = "")
gamma_sd              <- matrix(, nrow = n_sim, ncol = m * m)
colnames(gamma_sd)    <- paste("sd_S", rep(c(1:m), each = m), "_to_S", rep(1:m, m), sep = "")


# Run the model on the simulated data:
for(sim_i in 1:n_sim){
  set.seed(10103 + sim_i * 3)
  out_HMM_single_sim <- HMM(s_data = data1[[sim_i]]$obs, 
                            gen = list(m = m, n_dep = n_dep, q_emiss = q_emiss),
                            start_val = list(st_gamma, st_emiss_distr),
                            mcmc = list(J = J, burn_in = burn_in), show_progress = FALSE) 
  
  emiss_med[sim_i, ]  <- apply(out_HMM_single_sim$PD[burn_in:J, 1:sum(q_emiss*m)], 2, median)
  emiss_low[sim_i, ]  <- apply(out_HMM_single_sim$PD[burn_in:J, 1:sum(q_emiss*m)], 2, quantile, probs = 0.025)
  emiss_up[sim_i, ]   <- apply(out_HMM_single_sim$PD[burn_in:J, 1:sum(q_emiss*m)], 2, quantile, probs = 0.975)
  emiss_IQR[sim_i, ]  <- apply(out_HMM_single_sim$PD[burn_in:J, 1:sum(q_emiss*m)], 2, IQR)
  emiss_mean[sim_i, ] <- apply(out_HMM_single_sim$PD[burn_in:J, 1:sum(q_emiss*m)], 2, mean)
  emiss_sd[sim_i, ]   <- sqrt(apply(out_HMM_single_sim$PD[burn_in:J, 1:sum(q_emiss*m)], 2, var))
  
  gamma_med[sim_i, ]  <- apply(out_HMM_single_sim$PD[burn_in:J, (sum(q_emiss*m) + 1):(sum(q_emiss*m) + m*m)], 2, median)
  gamma_low[sim_i, ]  <- apply(out_HMM_single_sim$PD[burn_in:J, (sum(q_emiss*m) + 1):(sum(q_emiss*m) + m*m)], 2, quantile, probs = 0.025)
  gamma_up[sim_i, ]   <- apply(out_HMM_single_sim$PD[burn_in:J, (sum(q_emiss*m) + 1):(sum(q_emiss*m) + m*m)], 2, quantile, probs = 0.975)
  gamma_IQR[sim_i, ]  <- apply(out_HMM_single_sim$PD[burn_in:J, (sum(q_emiss*m) + 1):(sum(q_emiss*m) + m*m)], 2, IQR)
  gamma_mean[sim_i, ] <- apply(out_HMM_single_sim$PD[burn_in:J, (sum(q_emiss*m) + 1):(sum(q_emiss*m) + m*m)], 2, mean)
  gamma_sd[sim_i, ]   <- sqrt(apply(out_HMM_single_sim$PD[burn_in:J, (sum(q_emiss*m) + 1):(sum(q_emiss*m) + m*m)], 2, var))
  
}

out_sim <- list(emiss_med = emiss_med,
                emiss_low = emiss_low,
                emiss_up = emiss_up,
                emiss_IQR = emiss_IQR,
                emiss_mean = emiss_mean,
                emiss_sd = emiss_sd,
                
                gamma_med = gamma_med, 
                gamma_low = gamma_low,
                gamma_up = gamma_up,
                gamma_IQR = gamma_IQR,
                gamma_mean = gamma_mean,
                gamma_sd = gamma_sd
)


save(out_sim,file= paste(paste("sim_HMM", sim_type, "obs", q_emiss,"t", n_t, sep = "_"), ".rda", sep = ""))
