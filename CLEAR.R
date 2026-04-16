library(MASS)
library(mclust)
# library(glmnet)

#' Applies local GMMs, creates a global reference distribution, and estimates IPW 
#' weights via quadratic logistic regression to recover EDA statistics.
CLEAR <- function(site_data_list, K_comp = 2, n_0 = 10000, inf_factor = 2) {
  H <- length(site_data_list)
  var_names <- colnames(site_data_list[[1]])
  d <- length(var_names)
  
  n_h_vec <- sapply(site_data_list, nrow)
  prob_h <- n_h_vec / sum(n_h_vec)
  
  # Step 1: Local GMMs
  site_models <- list()
  for (h in 1:H) {
    df_h <- site_data_list[[h]]
    fit <- Mclust(df_h, G = K_comp, modelNames = "VVV", verbose = FALSE) # EM for GMM
    if (is.null(fit)) fit <- Mclust(df_h, G = 1:K_comp, verbose = FALSE)
    # if (is.null(fit)) fit <- Mclust(df_h, G = 1, verbose = FALSE)
    
    site_models[[h]] <- list(pi = fit$parameters$pro, 
                             mu = as.matrix(fit$parameters$mean), 
                             Sigma = fit$parameters$variance$sigma)
  }
  
  # Step 2: Reference Data
  Z_0 <- matrix(0, nrow = n_0, ncol = d)
  for (i in 1:n_0) {
    h <- sample(1:H, 1, prob = prob_h)
    m <- site_models[[h]]
    
    if (is.null(m$pi)) { k <- 1 } else { k <- sample(1:length(m$pi), 1, prob = m$pi) }
    mu_hk <- m$mu[, k]
    if (length(dim(m$Sigma)) == 3) { Sigma_hk <- m$Sigma[,, k] } else { Sigma_hk <- m$Sigma }
    
    Sigma_inflated <- Sigma_hk * inf_factor
    Z_0[i, ] <- mvrnorm(n = 1, mu = mu_hk, Sigma = Sigma_inflated)
  }
  Z_0 <- as.data.frame(Z_0)
  colnames(Z_0) <- var_names
  
  # Step 3 & 4: IPW and Get EDA stats
  quad_terms <- paste0("I(", var_names, "^2)", collapse = " + ")
  poly_form <- as.formula(paste("~ .^2 +", quad_terms))

  results_list <- list()
  for (h in 1:H) {
    df_h <- site_data_list[[h]]
    n_h <- nrow(df_h)
    Z_comb <- rbind(df_h, Z_0)
    Y_comb <- c(rep(1, n_h), rep(0, n_0))
    Z_poly <- as.data.frame(model.matrix(poly_form, data = Z_comb))
    Z_poly$Y_comb <- Y_comb
  
    # print(h)
    glm_fit <- suppressWarnings(glm(Y_comb ~ . -1, data = Z_poly, family = binomial))
    # check glm_fit convergence
    if (!glm_fit$converged) {warning(sprintf("GLM did not converge for Site %d", h))}
  
    Z0_poly <- Z_poly[(n_h + 1):nrow(Z_poly), ]
    prob <- predict(glm_fit, newdata = Z0_poly, type = "response")
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