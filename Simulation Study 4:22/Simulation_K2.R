rm(list = ls())
library(MASS)
library(mclust)
library(ggplot2)
library(dplyr)
library(tidyr)
library(foreach)
library(doParallel)
library(parallelly)
library(patchwork)

source("CLEAR.R")

# =====================================================================
# 1. Generalized DGP Function
# =====================================================================
generate_scalable_dgp <- function(M = 5, d = 5) {
  n_k <- sample(300:600, M, replace = TRUE)
  site_data_list <- list()
  
  for (k in 1:M) {
    rho <- runif(1, 0.2, 0.7)
    R <- matrix(rho, nrow = d, ncol = d)
    diag(R) <- 1
    
    Z <- mvrnorm(n = n_k[k], mu = rep(0, d), Sigma = R)
    U <- pnorm(Z)
    
    df_sim <- data.frame(matrix(ncol = d, nrow = n_k[k]))
    for(v in 1:d) {
      if (v %% 4 == 1) {
        df_sim[, v] <- qnorm(U[, v], mean = 50 + k*2, sd = 10)
      } else if (v %% 4 == 2) {
        df_sim[, v] <- qgamma(U[, v], shape = 2, rate = 0.5 + k*0.1)
      } else if (v %% 4 == 3) {
        df_sim[, v] <- qweibull(U[, v], shape = 1.5, scale = 15)
      } else {
        df_sim[, v] <- qlnorm(U[, v], meanlog = 3, sdlog = 0.5)
      }
    }
    colnames(df_sim) <- paste0("X", 1:d)
    site_data_list[[k]] <- df_sim
  }
  return(site_data_list)
}

# =====================================================================
# 2. Parallel Cluster Setup
# =====================================================================
n_cores <- parallelly::availableCores() - 1
cat(sprintf("Initializing Parallel Cluster with %d cores...\n", n_cores))

cl <- parallel::makeCluster(n_cores)
doParallel::registerDoParallel(cl)

# =====================================================================
# 3. Parallel Simulation Execution
# =====================================================================
# We vary the number of GMM components while keeping other parameters fixed
k_values <- c(2, 3, 5, 7)
fixed_d <- 5 
fixed_n0 <- 10000 
fixed_inf <- 1.5
n_sims <- 500  
target_site <- 1

# Flatten the loops into a single grid of tasks
sim_grid <- expand.grid(KComp = k_values, Simulation = 1:n_sims)
total_tasks <- nrow(sim_grid)

cat(sprintf("Starting %d total simulation tasks for varying K components...\n", total_tasks))

sim_results <- foreach(
  i = 1:total_tasks, 
  .combine = rbind,
  .packages = c("MASS", "mclust")
) %dopar% {
  
  # Extract current task parameters
  current_k <- sim_grid$KComp[i]
  sim <- sim_grid$Simulation[i]
  
  # 1. Generate Data (Fixed dimension)
  site_data_list <- generate_scalable_dgp(M = 5, d = fixed_d)
  truth_site <- calc_empirical_truth(site_data_list[[target_site]])
  
  # 2. Run CLEAR with the varying K_comp parameter
  est_site <- tryCatch({
    res <- CLEAR(site_data_list, K_comp = current_k, n_0 = fixed_n0, inf_factor = fixed_inf)
    res[[paste0("Site_", target_site)]]
  }, warning = function(w) {
    res <- suppressWarnings(CLEAR(site_data_list, K_comp = current_k, n_0 = fixed_n0, inf_factor = fixed_inf))
    res[[paste0("Site_", target_site)]]
  }, error = function(e) {
    return(NULL) 
  })
  
  # 3. Format and return raw estimates
  if (!is.null(est_site) && !is.null(est_site$Mean)) {
    df_out <- data.frame()
    
    # Means
    df_out <- rbind(df_out, data.frame(Metric = "Mean", Variable = 1:fixed_d, Est = est_site$Mean, Truth = truth_site$Mean))
    
    # Variances
    df_out <- rbind(df_out, data.frame(Metric = "Variance", Variable = 1:fixed_d, Est = est_site$Variance, Truth = truth_site$Variance))
    
    # Quantiles
    q_names <- c("Quantile_05", "Quantile_25", "Quantile_50", "Quantile_75", "Quantile_95")
    for(j in 1:5) {
      df_out <- rbind(df_out, data.frame(Metric = q_names[j], Variable = 1:fixed_d, Est = est_site$Quantiles[, j], Truth = truth_site$Quantiles[, j]))
    }
    
    # Covariances (flattened matrix)
    df_out <- rbind(df_out, data.frame(Metric = "Covariance", Variable = 1:(fixed_d^2), Est = as.vector(est_site$Covariance), Truth = as.vector(truth_site$Covariance)))
    
    # Append iteration keys
    df_out$KComp <- current_k
    df_out$Simulation <- sim
    
    return(df_out)
  } else {
    return(NULL)
  }
}

parallel::stopCluster(cl)

# =====================================================================
# 4. Save Raw Results
# =====================================================================
saveRDS(sim_results, "Results/Simulation_K_Results.rds")

# =====================================================================
# 5. Calculate Bias, Variance, and Percentage Bias
# =====================================================================
summary_metrics <- sim_results %>%
  group_by(KComp, Metric, Variable) %>%
  summarize(
    Bias_v = mean(Est - Truth, na.rm = TRUE),
    Var_v = var(Est, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(KComp, Metric) %>%
  summarize(
    MAE = abs(mean(Bias_v, na.rm = TRUE)),
    Avg_Bias = mean(Bias_v, na.rm = TRUE),
    Avg_Variance = mean(Var_v, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # Convert Metric to a factored categorical layout for logical legend ordering
  mutate(Metric = factor(Metric, levels = c("Mean", "Variance", "Covariance", 
                                            "Quantile_05", "Quantile_25", "Quantile_50", 
                                            "Quantile_75", "Quantile_95")))

# =====================================================================
# 6. Plotting
# =====================================================================

# Plot 1: Average Bias
p_bias <- ggplot(summary_metrics, aes(x = KComp, y = Avg_Bias, color = Metric, group = Metric)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_color_brewer(palette = "Set2") + 
  scale_x_continuous(breaks = k_values) + # Force strict integer breaks matching components
  theme_bw(base_size = 14) +
  labs(title = "Average Bias of Estimators", 
       x = "Number of GMM Components (K)", 
       y = "Average Bias")

# Plot 2: Variance
p_var <- ggplot(summary_metrics, aes(x = KComp, y = Avg_Variance, color = Metric, group = Metric)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_color_brewer(palette = "Set2") + 
  scale_y_log10() + 
  scale_x_continuous(breaks = k_values) +
  theme_bw(base_size = 14) +
  labs(title = "Variance of Estimators (Log Scale)", 
       x = "Number of GMM Components (K)", 
       y = "Variance")

# Combine the plots using patchwork
combined_plot <- p_bias + p_var + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")

# Display the combined plot
print(combined_plot)

# Save plots
ggsave("Results/Simulation_Plot_KComp.png", combined_plot, width = 12, height = 6, dpi = 300)