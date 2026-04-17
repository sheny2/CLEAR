library(MASS)
library(mclust)
library(caret)  # For KNN and data partitioning

#' Applies local GMMs, creates a global reference distribution, and estimates IPW 
#' weights via KNN Classification to recover federated statistics.
CLEAR_KNN <- function(site_data_list, K_comp = 3, n_0 = 10000, inf_cfactor = 2, k_neighbors = 25) {
  H <- length(site_data_list)
  var_names <- colnames(site_data_list[[1]])
  d <- length(var_names)
  
  n_h_vec <- sapply(site_data_list, nrow)
  prob_h <- n_h_vec / sum(n_h_vec)
  
  # --- Step 1 & 2: Local GMMs and Reference Data Generation ---
  # (Keeping your original logic for GMM generation)
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
    Z_0[i, ] <- mvrnorm(n = 1, mu = mu_hk, Sigma = Sigma_hk * inf_cfactor)
  }
  Z_0 <- as.data.frame(Z_0)
  colnames(Z_0) <- var_names
  
  # --- Step 3: Density Ratio Estimation via KNN ---
  results_list <- list()
  
  for (h in 1:H) {
    df_h <- site_data_list[[h]]
    n_h <- nrow(df_h)
    
    # Prepare training data: Site h (Class 1) and Reference Z_0 (Class 0)
    # Note: KNN is sensitive to scale; standardization is mandatory here
    X_train <- rbind(df_h, Z_0)
    Y_train <- as.factor(c(rep("Local", n_h), rep("Ref", n_0)))
    
    # KNN Prediction for Z_0
    # We use knn() from the 'class' package with prob = TRUE
    # This returns the proportion of the votes for the winning class
    library(class)
    
    # Pre-scaling both datasets based on combined data
    X_train_scaled <- scale(X_train)
    Z_0_scaled <- X_train_scaled[(n_h + 1):(n_h + n_0), , drop = FALSE]
    Site_h_scaled <- X_train_scaled[1:n_h, , drop = FALSE]
    
    # Run KNN to get probabilities that Z_0 belongs to the 'Local' class
    # We predict classes for every point in Z_0
    knn_fit <- knn(train = X_train_scaled, test = Z_0_scaled, 
                   cl = Y_train, k = k_neighbors, prob = TRUE)
    
    # Extract the probability of the winning class
    winning_prob <- attr(knn_fit, "prob")
    
    # Convert 'winning prob' to 'probability of being Local'
    # If knn predicted 'Ref', the prob of 'Local' is (1 - winning_prob)
    prob <- ifelse(knn_fit == "Local", winning_prob, 1 - winning_prob)
    
    # --- Step 4: Weighting and Stats ---
    # Bound probabilities to avoid infinity in odds ratio
    prob <- pmin(pmax(prob, 1e-6), 1 - 1e-6)
    
    omega <- (n_0 / n_h) * (prob / (1 - prob)) 
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