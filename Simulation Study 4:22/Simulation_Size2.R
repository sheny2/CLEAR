rm(list = ls())
library(MASS)
library(mclust)
library(ggplot2)
library(dplyr)
library(tidyr)
library(foreach)
library(doParallel)
library(parallelly)

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
n0_values <- c(2000, 5000, 10000, 20000)
fixed_d <- 5 
n_sims <- 500
target_site <- 1

sim_grid <- expand.grid(N0 = n0_values, Simulation = 1:n_sims)
total_tasks <- nrow(sim_grid)

cat(sprintf("Starting %d total simulation tasks for varying n_0...\n", total_tasks))

sim_results <- foreach(
  i = 1:total_tasks, 
  .combine = rbind,
  .packages = c("MASS", "mclust")
) %dopar% {
  
  current_n0 <- sim_grid$N0[i]
  sim <- sim_grid$Simulation[i]
  
  site_data_list <- generate_scalable_dgp(M = 5, d = fixed_d)
  truth_site <- calc_empirical_truth(site_data_list[[target_site]])
  
  est_site <- tryCatch({
    res <- CLEAR(site_data_list, K_comp = 3, n_0 = current_n0, inf_factor = 1.5)
    res[[paste0("Site_", target_site)]]
  }, warning = function(w) {
    res <- suppressWarnings(CLEAR(site_data_list, K_comp = 3, n_0 = current_n0, inf_factor = 1.5))
    res[[paste0("Site_", target_site)]]
  }, error = function(e) {
    return(NULL) 
  })
  
  if (!is.null(est_site) && !is.null(est_site$Mean)) {
    # Initialize a local dataframe to hold long-format raw estimates
    df_out <- data.frame()
    
    # 1. Means
    df_out <- rbind(df_out, data.frame(Metric = "Mean", Variable = 1:fixed_d, Est = est_site$Mean, Truth = truth_site$Mean))
    
    # 2. Variances
    df_out <- rbind(df_out, data.frame(Metric = "Variance", Variable = 1:fixed_d, Est = est_site$Variance, Truth = truth_site$Variance))
    
    # 3. Quantiles
    q_names <- c("Quantile_05", "Quantile_25", "Quantile_50", "Quantile_75", "Quantile_95")
    for(j in 1:5) {
      df_out <- rbind(df_out, data.frame(Metric = q_names[j], Variable = 1:fixed_d, Est = est_site$Quantiles[, j], Truth = truth_site$Quantiles[, j]))
    }
    
    # 4. Covariances (flattened matrix)
    df_out <- rbind(df_out, data.frame(Metric = "Covariance", Variable = 1:(fixed_d^2), Est = as.vector(est_site$Covariance), Truth = as.vector(truth_site$Covariance)))
    
    # Append iteration keys
    df_out$N0 <- current_n0
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
# Save the exact estimates and truth values before doing aggregate math
saveRDS(sim_results, "Results/Simulation_N0_Results.rds")

# =====================================================================
# 5. Calculate Bias, Variance, and Percentage Bias
# =====================================================================
summary_metrics <- sim_results %>%
  group_by(N0, Metric, Variable) %>%
  summarize(
    Bias_v = mean(Est - Truth, na.rm = TRUE),
    Var_v = var(Est, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(N0, Metric) %>%
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
p_bias <- ggplot(summary_metrics, aes(x = N0, y = Avg_Bias, color = Metric, group = Metric)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_color_brewer(palette = "Set2") + 
  theme_bw(base_size = 14) +
  labs(title = "Average Bias of Estimators", 
       x = expression("Size of Reference Data (" * n[0] * ")"), 
       y = "Average Bias")

# Plot 2: Variance
p_var <- ggplot(summary_metrics, aes(x = N0, y = Avg_Variance, color = Metric, group = Metric)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_color_brewer(palette = "Set2") + 
  scale_y_log10() + 
  theme_bw(base_size = 14) +
  labs(title = "Variance of Estimators (Log Scale)", 
       x = expression("Size of Reference Data (" * n[0] * ")"), 
       y = "Variance")

library(patchwork)
combined_plot <- p_bias + p_var + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")

# Display the combined plot
print(combined_plot)

# Save plots
ggsave("Results/Simulation_Plot_N0.png", combined_plot, width = 12, height = 6, dpi = 300)
