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

n_cores <- parallelly::availableCores() - 1
cat(sprintf("Initializing Parallel Cluster with %d cores...\n", n_cores))

cl <- parallel::makeCluster(n_cores)
doParallel::registerDoParallel(cl)

inf_values <- c(1.0, 1.5, 2.0, 3.0, 5.0)
fixed_d <- 5 
fixed_n0 <- 10000 
n_sims <- 500
target_site <- 1

sim_grid <- expand.grid(InfFactor = inf_values, Simulation = 1:n_sims)
total_tasks <- nrow(sim_grid)

cat(sprintf("Starting %d total tasks for varying Inflation Factor...\n", total_tasks))

sim_results <- foreach(
  i = 1:total_tasks, 
  .combine = rbind,
  .packages = c("MASS", "mclust")
) %dopar% {
  
  current_inf <- sim_grid$InfFactor[i]
  sim <- sim_grid$Simulation[i]
  
  site_data_list <- generate_scalable_dgp(M = 5, d = fixed_d)
  truth_site <- calc_empirical_truth(site_data_list[[target_site]])
  
  est_site <- tryCatch({
    res <- CLEAR(site_data_list, K_comp = 3, n_0 = fixed_n0, inf_factor = current_inf)
    res[[paste0("Site_", target_site)]]
  }, warning = function(w) {
    res <- suppressWarnings(CLEAR(site_data_list, K_comp = 3, n_0 = fixed_n0, inf_factor = current_inf))
    res[[paste0("Site_", target_site)]]
  }, error = function(e) {
    return(NULL) 
  })
  
  if (!is.null(est_site) && !is.null(est_site$Mean)) {
    df_out <- data.frame()
    
    df_out <- rbind(df_out, data.frame(Metric = "Mean", Variable = 1:fixed_d, Est = est_site$Mean, Truth = truth_site$Mean))
    df_out <- rbind(df_out, data.frame(Metric = "Variance", Variable = 1:fixed_d, Est = est_site$Variance, Truth = truth_site$Variance))
    
    q_names <- c("Quantile_05", "Quantile_25", "Quantile_50", "Quantile_75", "Quantile_95")
    for(j in 1:5) {
      df_out <- rbind(df_out, data.frame(Metric = q_names[j], Variable = 1:fixed_d, Est = est_site$Quantiles[, j], Truth = truth_site$Quantiles[, j]))
    }
    
    df_out <- rbind(df_out, data.frame(Metric = "Covariance", Variable = 1:(fixed_d^2), Est = as.vector(est_site$Covariance), Truth = as.vector(truth_site$Covariance)))
    
    df_out$InfFactor <- current_inf
    df_out$Simulation <- sim
    
    return(df_out)
  } else {
    return(NULL)
  }
}

parallel::stopCluster(cl)

saveRDS(sim_results, "Results/Simulation_InfFactor_Results.rds")

summary_metrics <- sim_results %>%
  group_by(InfFactor, Metric, Variable) %>%
  summarize(
    Bias_v = mean(Est - Truth, na.rm = TRUE),
    Var_v = var(Est, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(InfFactor, Metric) %>%
  summarize(
    MAE = abs(mean(Bias_v, na.rm = TRUE)),
    Avg_Bias = mean(Bias_v, na.rm = TRUE),
    Avg_Variance = mean(Var_v, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Metric = factor(Metric, levels = c("Mean", "Variance", "Covariance", 
                                            "Quantile_05", "Quantile_25", "Quantile_50", 
                                            "Quantile_75", "Quantile_95")))

p_bias <- ggplot(summary_metrics, aes(x = InfFactor, y = Avg_Bias, color = Metric, group = Metric)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_color_brewer(palette = "Set2") + 
  scale_x_continuous(breaks = inf_values) +
  theme_bw(base_size = 14) +
  labs(title = "Average Bias of Estimators", 
       x = "Covariance Inflation Factor (inf_factor)", 
       y = "Average Bias")

p_var <- ggplot(summary_metrics, aes(x = InfFactor, y = Avg_Variance, color = Metric, group = Metric)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_color_brewer(palette = "Set2") + 
  scale_y_log10() + 
  scale_x_continuous(breaks = inf_values) +
  theme_bw(base_size = 14) +
  labs(title = "Variance of Estimators (Log Scale)", 
       x = "Covariance Inflation Factor (inf_factor)", 
       y = "Variance")

combined_plot <- p_bias + p_var + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")

print(combined_plot)
ggsave("Results/Simulation_Plot_InfFactor.png", combined_plot, width = 12, height = 6, dpi = 300)