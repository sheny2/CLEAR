library(MASS)
library(mclust)
library(caret) 
library(class)

#' Applies local GMMs, creates a global reference distribution, and estimates IPW 
#' weights via KNN Classification to recover federated statistics.
CLEAR_KNN <- function(site_data_list, K_comp = 3, n_0 = 10000, inf_cfactor = 2, k_neighbors = 15) {
  H <- length(site_data_list)
  var_names <- colnames(site_data_list[[1]])
  d <- length(var_names)
  
  n_h_vec <- sapply(site_data_list, nrow)
  prob_h <- n_h_vec / sum(n_h_vec)
  
  # --- Step 1 & 2: Local GMMs and Reference Data Generation ---
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
  q_probs <- c(0.05, 0.25, 0.5, 0.75, 0.95) # Probabilities for Quantiles
  
  for (h in 1:H) {
    df_h <- site_data_list[[h]]
    n_h <- nrow(df_h)
    
    X_train <- rbind(df_h, Z_0)
    Y_train <- as.factor(c(rep("Local", n_h), rep("Ref", n_0)))
    
    X_train_scaled <- scale(X_train)
    Z_0_scaled <- X_train_scaled[(n_h + 1):(n_h + n_0), , drop = FALSE]
    
    knn_fit <- knn(train = X_train_scaled, test = Z_0_scaled, 
                   cl = Y_train, k = k_neighbors, prob = TRUE)
    
    winning_prob <- attr(knn_fit, "prob")
    prob <- ifelse(knn_fit == "Local", winning_prob, 1 - winning_prob)
    
    # --- Step 4: Weighting and Stats ---
    prob <- pmin(pmax(prob, 1e-6), 1 - 1e-6)
    
    omega <- (n_0 / n_h) * (prob / (1 - prob)) 
    threshold <- quantile(omega, 0.99, na.rm = TRUE)
    w_h <- pmin(omega, threshold)
    w_norm <- w_h / sum(w_h)
    
    site_stats <- list(
      Mean = numeric(d), 
      Variance = numeric(d), 
      Quantiles = matrix(0, nrow=d, ncol=length(q_probs)),
      Covariance = matrix(0, nrow=d, ncol=d)
    )
    
    names(site_stats$Mean) <- names(site_stats$Variance) <- var_names
    rownames(site_stats$Quantiles) <- rownames(site_stats$Covariance) <- colnames(site_stats$Covariance) <- var_names
    colnames(site_stats$Quantiles) <- c("5%", "25%", "50%", "75%", "95%")
    
    for (v in 1:d) {
      x <- Z_0[, v]
      
      # Mean and Variance
      mu_v <- sum(w_norm * x)
      var_v <- sum(w_norm * (x - mu_v)^2)
      
      site_stats$Mean[v] <- mu_v
      site_stats$Variance[v] <- var_v
      
      # Weighted Quantiles (Base R implementation)
      ord <- order(x)
      x_ord <- x[ord]
      w_ord <- w_norm[ord]
      cum_w <- cumsum(w_ord) 
      site_stats$Quantiles[v, ] <- sapply(q_probs, function(p) x_ord[which(cum_w >= p)[1]])
    }
    
    # Vectorized Covariance Calculation
    Z_centered <- sweep(as.matrix(Z_0), 2, site_stats$Mean, FUN = "-")
    site_stats$Covariance <- t(Z_centered) %*% (Z_centered * w_norm)
    
    site_stats$Weights <- w_norm
    results_list[[paste0("Site_", h)]] <- site_stats
  }
  
  results_list$Z_0 <- Z_0 
  return(results_list)
}


#' Helper to calculate exact empirical moments and quantiles of the raw data (Truth)
calc_empirical_truth <- function(df) {
  n <- nrow(df)
  
  # 1. Means
  means <- colMeans(df)
  
  # 2. Population Variance (div by n instead of n-1 to match IPW mathematical bounds)
  vars <- apply(df, 2, function(x) sum((x - mean(x))^2) / n) 
  
  # 3. Empirical Quantiles (Replacing Skewness and Kurtosis)
  q_probs <- c(0.05, 0.25, 0.50, 0.75, 0.95)
  
  # apply() returns quantiles as rows and variables as columns. 
  # We transpose t() it so rows are variables, matching the CLEAR algorithm output.
  quant <- t(apply(df, 2, quantile, probs = q_probs))
  colnames(quant) <- c("5%", "25%", "50%", "75%", "95%")
  
  # 4. Population Covariance Matrix
  cov_mat <- cov(df) * ((n - 1) / n)
  
  list(
    Mean = means, 
    Variance = vars, 
    Quantiles = quant, 
    Covariance = cov_mat
  )
}