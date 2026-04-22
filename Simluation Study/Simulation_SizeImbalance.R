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
# 1. Imbalance-Aware DGP Function
# =====================================================================
# We modify the DGP to accept an exact vector of sample sizes (n_k_vec)
generate_imbalanced_dgp <- function(n_k_vec, d = 5) {
  M <- length(n_k_vec)
  site_data_list <- list()
  
  for (k in 1:M) {
    # Generate distinct geometries for each site
    rho <- runif(1, 0.2, 0.7)
    R <- matrix(rho, nrow = d, ncol = d)
    diag(R) <- 1
    
    Z <- mvrnorm(n = n_k_vec[k], mu = rep(0, d), Sigma = R)
    U <- pnorm(Z)
    
    df_sim <- data.frame(matrix(ncol = d, nrow = n_k_vec[k]))
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
# 3. Define Imbalance Scenarios
# =====================================================================
# Total network N = 2000 across 5 sites.
# We progressively shift data from the satellite sites to Site 1.
scenarios <- list(
  "1_Balanced" = c(400, 400, 400, 400, 400),
  "2_Slight"   = c(800, 400, 300, 300, 200),
  "3_Moderate" = c(1200, 300, 200, 150, 150),
  "4_Extreme"  = c(1600, 100, 100, 100, 100),
  "5_Severe"   = c(1800, 50, 50, 50, 50)
)

fixed_d <- 5 
n_sims <- 100
# CRITICAL: We evaluate the LAST site (Site 5), which is always the smallest.
# We want to see if the massive Site 1 destroys the estimates for the tiny Site 5.
target_site <- 5 

# Flatten the loops into a single grid of tasks
sim_grid <- expand.grid(Scenario = names(scenarios), Simulation = 1:n_sims)
total_tasks <- nrow(sim_grid)

cat(sprintf("Starting %d total simulation tasks for Sample Size Imbalance...\n", total_tasks))

sim_results <- foreach(
  i = 1:total_tasks, 
  .combine = rbind,
  .packages = c("MASS", "mclust", "class")
) %dopar% {
  
  # Extract current task parameters
  scen_name <- as.character(sim_grid$Scenario[i])
  n_k_vec <- scenarios[[scen_name]]
  sim <- sim_grid$Simulation[i]
  
  # 1. Generate Data 
  site_data_list <- generate_imbalanced_dgp(n_k_vec = n_k_vec, d = fixed_d)
  truth_site <- calc_empirical_truth(site_data_list[[target_site]])
  
  # 2. Run CLEAR 
  est_site <- tryCatch({
    res <- CLEAR(site_data_list, K_comp = 2, n_0 = 10000, inf_factor = 1.5)
    res[[paste0("Site_", target_site)]]
  }, warning = function(w) {
    res <- suppressWarnings(CLEAR(site_data_list, K_comp = 2, n_0 = 10000, inf_factor = 1.5))
    res[[paste0("Site_", target_site)]]
  }, error = function(e) {
    return(NULL) 
  })
  
  # 3. Calculate Errors and Return as Dataframe row
  if (!is.null(est_site) && !is.null(est_site$Mean)) {
    mae_mean  <- mean(abs(est_site$Mean - truth_site$Mean))
    mae_var   <- mean(abs(est_site$Variance - truth_site$Variance))
    mae_quant <- mean(abs(est_site$Quantiles - truth_site$Quantiles))
    mae_cov   <- mean(abs(est_site$Covariance - truth_site$Covariance))
    
    return(data.frame(
      Scenario = scen_name,
      Target_N = n_k_vec[target_site], # Track the sample size of the site being evaluated
      Simulation = sim,
      Mean_Error = mae_mean,
      Var_Error = mae_var,
      Quantile_Error = mae_quant,
      Covariance_Error = mae_cov
    ))
  } else {
    return(NULL)
  }
}

parallel::stopCluster(cl)

# =====================================================================
# 4. Data Transformation and Plotting
# =====================================================================
# Format scenario names for cleaner plotting
results_long <- sim_results %>%
  mutate(
    Scenario_Label = case_when(
      Scenario == "1_Balanced" ~ "Balanced\n(Target N=400)",
      Scenario == "2_Slight"   ~ "Slight\n(Target N=200)",
      Scenario == "3_Moderate" ~ "Moderate\n(Target N=150)",
      Scenario == "4_Extreme"  ~ "Extreme\n(Target N=100)",
      Scenario == "5_Severe"   ~ "Severe\n(Target N=50)"
    )
  ) %>%
  mutate(Scenario_Label = factor(Scenario_Label, levels = c(
    "Balanced\n(Target N=400)", 
    "Slight\n(Target N=200)", 
    "Moderate\n(Target N=150)", 
    "Extreme\n(Target N=100)", 
    "Severe\n(Target N=50)"
  ))) %>%
  pivot_longer(cols = ends_with("Error"), names_to = "Metric", values_to = "MAE") %>%
  mutate(Metric = factor(Metric, levels = c("Mean_Error", "Var_Error", "Quantile_Error", "Covariance_Error")))

eval_plot <- ggplot(results_long, aes(x = Scenario_Label, y = MAE, fill = Scenario_Label)) +
  geom_boxplot(alpha = 0.7, outlier.size = 1, outlier.alpha = 0.5) +
  facet_wrap(~ Metric, scales = "free_y", ncol = 2) +
  theme_bw(base_size = 14) +
  labs(
    title = "CLEAR Robustness to Network Sample Size Imbalance",
    subtitle = "Evaluating performance on the smallest satellite site across increasingly imbalanced networks",
    x = "Network Imbalance Scenario",
    y = "Mean Absolute Error (MAE)"
  ) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "#f0f0f0"),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# save the plot
ggsave("Results/Simulation_Imbalance.png", eval_plot, width = 12, height = 8, dpi = 300)