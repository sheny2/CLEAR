library(MASS)
library(mclust)
library(ks) # Kernel Smoothing package for multivariate KDE

#' Estimates IPW weights via Kernel Density Estimation (KDE) 
#' to recover federated statistics from a reference distribution.
CLEAR_KDE <- function(site_data_list, K_comp = 3, n_0 = 10000, inf_factor = 2) {
  H <- length(site_data_list)
  var_names <- colnames(site_data_list[[1]])
  d <- length(var_names)
  
  n_h_vec <- sapply(site_data_list, nrow)
  prob_h <- n_h_vec / sum(n_h_vec)
  
  # --- Step 1 & 2: Local GMMs and Reference Data Generation ---
  # (Keeping your original logic for GMM generation to ensure consistency)
  site_models <- list()
  for (h in 1:H) {
    df_h <- site_data_list[[h]]
    fit <- Mclust(df_h, G = K_comp, modelNames = "VVV", verbose = FALSE)
    if (is.null(fit)) fit <- Mclust(df_h, G = 1:K_comp, verbose = FALSE)
    site_models[[h]] <- list(pi = fit$parameters$pro, 
                             mu = as.matrix(fit$parameters$mean), 
                             Sigma = fit$parameters$variance$sigma)
  }
  
  Z_0 <- matrix(0, nrow = n_0, ncol = d)
  for (i in 1:n_0) {
    h <- sample(1:H, 1, prob = prob_h)
    m <- site_models[[h]]
    k <- if (is.null(m$pi)) 1 else sample(1:length(m$pi), 1, prob = m$pi)
    mu_hk <- m$mu[, k]
    Sigma_hk <- if (length(dim(m$Sigma)) == 3) m$Sigma[,, k] else m$Sigma
    Z_0[i, ] <- mvrnorm(n = 1, mu = mu_hk, Sigma = Sigma_hk * inf_factor)
  }
  Z_0 <- as.data.frame(Z_0)
  colnames(Z_0) <- var_names
  
  # --- Step 3: Density Ratio Estimation via KDE ---
  results_list <- list()
  
  for (h in 1:H) {
    df_h <- as.matrix(site_data_list[[h]])
    n_h <- nrow(df_h)
    Z0_mat <- as.matrix(Z_0)
    
    # 1. Estimate Bandwidth Matrix (H_mat) for Local Site
    # Hpi is a plug-in bandwidth selector; robust for multivariate data
    H_local <- Hpi(x = df_h)
    
    # 2. Estimate Bandwidth Matrix (H_mat) for Reference Z_0
    H_ref <- Hpi(x = Z0_mat)
    
    # 3. Compute Densities for each point in Z_0
    # f_local(z) = density of local data evaluated at reference points
    f_local <- kde(x = df_h, H = H_local, eval.points = Z0_mat)$estimate
    
    # f_ref(z) = density of reference data evaluated at reference points
    f_ref <- kde(x = Z0_mat, H = H_ref, eval.points = Z0_mat)$estimate
    
    # 4. Calculate Density Ratio weights
    # Note: No need for the (n_0/n_h) term here as f(x) is a true PDF
    omega <- f_local / f_ref
    
    # --- Step 4: Weighting and Stats ---
    # Guard against division by zero or extreme outliers
    omega[is.na(omega) | is.infinite(omega)] <- 0
    threshold <- quantile(omega, 0.99, na.rm = TRUE)
    w_h <- pmin(omega, threshold)
    w_norm <- w_h / sum(w_h)
    
    
    site_stats <- list(
      Mean = numeric(d), Variance = numeric(d), 
      Skewness = numeric(d), Kurtosis = numeric(d), 
      Covariance = matrix(0, nrow=d, ncol=d)
    )
    
    names(site_stats$Mean) <- names(site_stats$Variance) <- names(site_stats$Skewness) <- names(site_stats$Kurtosis) <- var_names
    rownames(site_stats$Covariance) <- colnames(site_stats$Covariance) <- var_names
    
    for (v in 1:d) {
      x <- Z_0[, v]
      mu_v <- sum(w_norm * x)
      var_v <- sum(w_norm * (x - mu_v)^2)
      sd_v <- sqrt(var_v)
      
      site_stats$Mean[v] <- mu_v
      site_stats$Variance[v] <- var_v
      site_stats$Skewness[v] <- sum(w_norm * (x - mu_v)^3) / (sd_v^3)
      site_stats$Kurtosis[v] <- sum(w_norm * (x - mu_v)^4) / (sd_v^4)
    }
    for (u in 1:d) {
      for (v in 1:d) {
        site_stats$Covariance[u, v] <- sum(w_norm * (Z_0[, u] - site_stats$Mean[u]) * (Z_0[, v] - site_stats$Mean[v]))
      }
    }
    site_stats$Weights <- w_norm
    results_list[[paste0("Site_", h)]] <- site_stats
  }
  
  results_list$Z_0 <- Z_0 
  return(results_list)
}





library(future.apply) # dependency for parallelization

CLEAR_KDE <- function(site_data_list, K_comp = 3, n_0 = 10000, inf_cfactor = 2, workers = parallelly::availableCores() - 1) {
  
  # --- Setup Parallel Backend ---
  plan(multisession, workers = workers)
  
  H <- length(site_data_list)
  var_names <- colnames(site_data_list[[1]])
  d <- length(var_names)
  
  n_h_vec <- sapply(site_data_list, nrow)
  prob_h <- n_h_vec / sum(n_h_vec)
  
  # --- Step 1 & 2: Local GMMs and Reference Data Generation ---
  # (Keeping GMM sequential as it's usually fast compared to KDE)
  site_models <- list()
  for (h in 1:H) {
    df_h <- site_data_list[[h]]
    fit <- Mclust(df_h, G = K_comp, modelNames = "VVV", verbose = FALSE)
    if (is.null(fit)) fit <- Mclust(df_h, G = 1:K_comp, verbose = FALSE)
    site_models[[h]] <- list(pi = fit$parameters$pro, 
                             mu = as.matrix(fit$parameters$mean), 
                             Sigma = fit$parameters$variance$sigma)
  }
  
  # Generate Z_0 (Reference points)
  Z_0_mat <- matrix(0, nrow = n_0, ncol = d)
  for (i in 1:n_0) {
    h_idx <- sample(1:H, 1, prob = prob_h)
    m <- site_models[[h_idx]]
    k <- if (is.null(m$pi)) 1 else sample(1:length(m$pi), 1, prob = m$pi)
    mu_hk <- m$mu[, k]
    Sigma_hk <- if (length(dim(m$Sigma)) == 3) m$Sigma[,, k] else m$Sigma
    Z_0_mat[i, ] <- mvrnorm(n = 1, mu = mu_hk, Sigma = Sigma_hk * inf_cfactor)
  }
  colnames(Z_0_mat) <- var_names
  
  # --- Step 3: Parallelized Density Ratio Estimation ---
  results_list <- future_lapply(1:H, function(h) {
    
    df_h <- as.matrix(site_data_list[[h]])
    
    # 1. Bandwidth Selection
    H_local <- ks::Hpi(x = df_h)
    H_ref   <- ks::Hpi(x = Z_0_mat)
    
    # 2. KDE Densities
    f_local <- ks::kde(x = df_h, H = H_local, eval.points = Z_0_mat)$estimate
    f_ref   <- ks::kde(x = Z_0_mat, H = H_ref, eval.points = Z_0_mat)$estimate
    
    # 3. Density Ratio weights
    omega <- f_local / f_ref
    omega[is.na(omega) | is.infinite(omega)] <- 0
    threshold <- quantile(omega, 0.99, na.rm = TRUE)
    w_h <- pmin(omega, threshold)
    w_norm <- w_h / sum(w_h)
    
    # --- Step 4: Stats ---
    means <- colSums(w_norm * Z_0_mat)
    Z_centered <- sweep(Z_0_mat, 2, means, "-")
    
    vars <- colSums(w_norm * Z_centered^2)
    sds  <- sqrt(vars)
    
    skewness <- colSums(w_norm * Z_centered^3) / (sds^3)
    kurtosis <- colSums(w_norm * Z_centered^4) / (sds^4)
    
    # Covariance Matrix: t(Z*w) %*% Z
    cov_mat <- t(Z_centered * w_norm) %*% Z_centered
    
    # Clean up names
    names(means) <- names(vars) <- names(skewness) <- names(kurtosis) <- var_names
    dimnames(cov_mat) <- list(var_names, var_names)
    
    return(list(
      Mean = means, 
      Variance = vars, 
      Skewness = skewness, 
      Kurtosis = kurtosis, 
      Covariance = cov_mat,
      Weights = w_norm
    ))
  }, future.seed = TRUE) # Ensure random stability
  
  # Cleanup names and add Z_0 to return
  names(results_list) <- paste0("Site_", 1:H)
  results_list$Z_0 <- as.data.frame(Z_0_mat)
  
  # Reset parallel plan
  plan(sequential)
  
  return(results_list)
}



#' Helper to calculate exact empirical moments of the raw data (Truth)
calc_empirical_truth <- function(df) {
  n <- nrow(df)
  means <- colMeans(df)
  # Use population variance formula (div by n) to match IPW formula
  vars <- apply(df, 2, function(x) sum((x - mean(x))^2) / n) 
  sds <- sqrt(vars)
  
  skew <- numeric(ncol(df)); kurt <- numeric(ncol(df))
  for(v in 1:ncol(df)) {
    x <- df[, v]
    skew[v] <- sum((x - means[v])^3 / n) / (sds[v]^3)
    kurt[v] <- sum((x - means[v])^4 / n) / (sds[v]^4)
  }
  names(skew) <- names(kurt) <- colnames(df)
  
  # Population covariance matrix
  cov_mat <- cov(df) * ((n - 1) / n)
  
  list(Mean = means, Variance = vars, Skewness = skew, Kurtosis = kurt, Covariance = cov_mat)
}



