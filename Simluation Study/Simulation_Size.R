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
# Instead of varying dimensions, we vary n_0 and lock the dimension
n0_values <- c(2000, 5000, 10000, 20000)
fixed_d <- 5 # Lock the dimension (5 is stable and fast for GLM)
n_sims <- 100
target_site <- 1

# Flatten the loops into a single grid of tasks
sim_grid <- expand.grid(N0 = n0_values, Simulation = 1:n_sims)
total_tasks <- nrow(sim_grid)

cat(sprintf("Starting %d total simulation tasks for varying n_0...\n", total_tasks))

sim_results <- foreach(
  i = 1:total_tasks, 
  .combine = rbind,
  .packages = c("MASS", "mclust", "class")
) %dopar% {
  
  # Extract current task parameters
  current_n0 <- sim_grid$N0[i]
  sim <- sim_grid$Simulation[i]
  
  # 1. Generate Data (Fixed dimension)
  site_data_list <- generate_scalable_dgp(M = 5, d = fixed_d)
  truth_site <- calc_empirical_truth(site_data_list[[target_site]])
  
  # 2. Run CLEAR with the varying n_0 parameter
  est_site <- tryCatch({
    res <- CLEAR(site_data_list, K_comp = 2, n_0 = current_n0, inf_factor = 1.5)
    res[[paste0("Site_", target_site)]]
  }, warning = function(w) {
    res <- suppressWarnings(CLEAR(site_data_list, K_comp = 2, n_0 = current_n0, inf_factor = 1.5))
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
      N0 = current_n0,
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
results_long <- sim_results %>%
  pivot_longer(cols = ends_with("Error"), 
               names_to = "Metric", 
               values_to = "MAE") %>%
  mutate(Metric = factor(Metric, levels = c("Mean_Error", "Var_Error", "Quantile_Error", "Covariance_Error")),
         N0 = as.factor(N0)) # Convert N0 to factor for categorical plotting

eval_plot <- ggplot(results_long, aes(x = N0, y = MAE, fill = N0)) +
  geom_boxplot(alpha = 0.7, outlier.size = 1, outlier.alpha = 0.5) +
  facet_wrap(~ Metric, scales = "free_y", ncol = 2) +
  theme_bw(base_size = 14) +
  labs(
    title = bquote("CLEAR Algorithm Performance by Reference Data Size ("*n[0]*")"),
    subtitle = paste("Mean Absolute Error (Target Site", target_site, "| d =", fixed_d, ")"),
    x = expression("Size of Reference Data (" * n[0] * ")"),
    y = "Mean Absolute Error (MAE)",
    fill = expression(n[0])
  ) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "#f0f0f0"),
    strip.text = element_text(face = "bold")
  )

# save the plot
ggsave("Results/Simulation_N0.png", eval_plot, width = 12, height = 8, dpi = 300)