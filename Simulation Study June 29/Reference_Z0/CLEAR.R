# =============================================================================
# CLEAR (logistic-regression IPW)
# -----------------------------------------------------------------------------
# Design:
#   * Sites only ever export fitted GMM parameters (pi, mu, Sigma) -- no raw
#     data and no synthetic points leave the site.
#   * The analysis center pools the shared GMMs into a global reference
#     distribution, then recovers each site's statistics via density-ratio
#     importance weighting (quadratic logistic regression).
#   * A simulation wrapper drives the full pipeline end to end.
# =============================================================================

library(MASS)     # mvrnorm
library(mclust)   # Mclust


# =============================================================================
# SITE-SIDE FUNCTIONS
# Run locally at each data-holding site. Only the returned object is shared.
# =============================================================================

#' Fit a local GMM and return only the parameters (the shareable object).
#'
#' @param df          A site's data: numeric data.frame / matrix, rows = obs.
#' @param K_comp      Number of mixture components to request (VVV).
#' @param fallback    If TRUE, fall back to BIC selection over 1:K_comp when the
#'                    requested VVV fit fails.
#' @return A list with: pi, mu (d x K), Sigma (d x d x K), var_names, n.
#'         This is the ONLY thing that leaves the site.
site_fit_gmm <- function(df, K_comp = 3, fallback = TRUE) {
  var_names <- colnames(df)
  if (is.null(var_names)) stop("`df` must have column names (variable names).")
  
  fit <- Mclust(df, G = K_comp, modelNames = "VVV", verbose = FALSE)
  if (is.null(fit) && fallback) {
    warning(sprintf("VVV (G=%d) failed; falling back to BIC over G=1:%d.",
                    K_comp, K_comp))
    fit <- Mclust(df, G = 1:K_comp, verbose = FALSE)
  }
  if (is.null(fit)) {
    stop("Mclust failed. Check for issues like near-zero variance columns or too few rows.")
  }
  
  Sigma <- fit$parameters$variance$sigma
  # Ensure Sigma is always a d x d x K array for downstream consistency.
  if (length(dim(Sigma)) != 3) {
    d <- nrow(as.matrix(fit$parameters$mean))
    Sigma <- array(Sigma, dim = c(d, d, 1))
  }
  
  list(
    pi        = fit$parameters$pro,
    mu        = as.matrix(fit$parameters$mean),  # d x K
    Sigma     = Sigma,                           # d x d x K
    var_names = var_names,
    n         = nrow(df)
  )
}


# =============================================================================
# ANALYSIS-CENTER FUNCTIONS
# Run at the central node using only the GMM parameters shared by sites.
# =============================================================================

#' Build a global reference distribution by sampling from the pooled site GMMs.
#'
#' Each draw: pick a site (prob proportional to n_h), pick a component within
#' that site (prob = pi), then sample from that Gaussian.
#'
#' @param site_models   List of objects returned by site_fit_gmm().
#' @param n_0           Number of reference points to generate.
#' @param inflate       Optional variance inflation factor (default 1 = off).
#'                      Kept as an opt-in knob; >1 widens the reference cloud.
#' @return A data.frame Z_0 (n_0 x d) with the shared variable names.
center_make_reference <- function(site_models, n_0 = 10000, inflate = 1) {
  H         <- length(site_models)
  var_names <- site_models[[1]]$var_names
  d         <- length(var_names)
  
  n_h_vec <- sapply(site_models, function(m) m$n)
  prob_h  <- n_h_vec / sum(n_h_vec)
  
  Z_0 <- matrix(0, nrow = n_0, ncol = d)
  for (i in seq_len(n_0)) {
    h <- sample(seq_len(H), 1, prob = prob_h)
    m <- site_models[[h]]
    
    k        <- if (is.null(m$pi)) 1L else sample(seq_along(m$pi), 1, prob = m$pi)
    mu_hk    <- m$mu[, k]
    Sigma_hk <- m$Sigma[, , k]
    
    Z_0[i, ] <- mvrnorm(n = 1, mu = mu_hk, Sigma = Sigma_hk * inflate)
  }
  
  Z_0 <- as.data.frame(Z_0)
  colnames(Z_0) <- var_names
  Z_0
}


#' Estimate IPW (density-ratio) weights for one site against the reference.
#'
#' Classifies reference points vs. that site's data with a quadratic logistic
#' regression, then converts predicted probabilities into odds-ratio weights.
#'
#' @param df_h    The site's raw data (used only when this runs where data
#'                lives, e.g. in simulation; the center never sees it in a
#'                real deployment -- see note below).
#' @param Z_0     The reference distribution from center_make_reference().
#' @return Normalized weights (length n_0, summing to 1).
#'
#' NOTE: density-ratio IPW requires the site's actual points to fit the
#' classifier. In a real deployment this step is therefore run AT THE SITE
#' (the site receives Z_0, fits the GLM locally, and returns only w_norm),
#' or replaced by a GMM-density-ratio computed analytically from shared params.
#' It is exposed here as a center/site-shared routine for flexibility.
estimate_ipw_weights <- function(df_h, Z_0) {
  var_names <- colnames(Z_0)
  n_h <- nrow(df_h)
  n_0 <- nrow(Z_0)
  
  quad_terms <- paste0("I(", var_names, "^2)", collapse = " + ")
  poly_form  <- as.formula(paste("~ .^2 +", quad_terms))
  
  Z_comb <- rbind(df_h, Z_0)
  Y_comb <- c(rep(1, n_h), rep(0, n_0))
  Z_poly <- as.data.frame(model.matrix(poly_form, data = Z_comb))
  Z_poly$Y_comb <- Y_comb
  
  glm_fit <- suppressWarnings(
    glm(Y_comb ~ . - 1, data = Z_poly, family = binomial))
  if (!glm_fit$converged) warning("GLM did not converge.")
  
  Z0_poly <- Z_poly[(n_h + 1):nrow(Z_poly), ]
  prob    <- predict(glm_fit, newdata = Z0_poly, type = "response")
  
  prob   <- pmin(pmax(prob, 1e-6), 1 - 1e-6)
  omega  <- (n_0 / n_h) * (prob / (1 - prob))   # density-ratio / odds weights
  omega / sum(omega)                            # normalize (no clipping)
}


#' Compute weighted EDA statistics on the reference points.
#'
#' @param Z_0      Reference data.frame (n_0 x d).
#' @param w_norm   Normalized importance weights (length n_0).
#' @return A list: Mean, Variance, Quantiles, Covariance, Weights.
weighted_eda_stats <- function(Z_0, w_norm) {
  var_names <- colnames(Z_0)
  d         <- length(var_names)
  q_probs   <- c(0.05, 0.25, 0.50, 0.75, 0.95)
  
  stats <- list(
    Mean       = numeric(d),
    Variance   = numeric(d),
    Quantiles  = matrix(0, nrow = d, ncol = length(q_probs)),
    Covariance = matrix(0, nrow = d, ncol = d)
  )
  names(stats$Mean) <- names(stats$Variance) <- var_names
  rownames(stats$Quantiles)  <- var_names
  colnames(stats$Quantiles)  <- c("5%", "25%", "50%", "75%", "95%")
  rownames(stats$Covariance) <- colnames(stats$Covariance) <- var_names
  
  for (v in seq_len(d)) {
    x    <- Z_0[, v]
    mu_v <- sum(w_norm * x)
    stats$Mean[v]     <- mu_v
    stats$Variance[v] <- sum(w_norm * (x - mu_v)^2)
    
    # Weighted quantiles: sort, accumulate weight, take first crossing.
    ord   <- order(x)
    x_ord <- x[ord]
    cum_w <- cumsum(w_norm[ord])
    stats$Quantiles[v, ] <- sapply(q_probs, function(p) x_ord[which(cum_w >= p)[1]])
  }
  
  Z_centered       <- sweep(as.matrix(Z_0), 2, stats$Mean, FUN = "-")
  stats$Covariance <- t(Z_centered) %*% (Z_centered * w_norm)
  stats$Weights    <- w_norm
  stats
}


# =============================================================================
# SIMULATION WRAPPER
# One call that runs the full pipeline. Convenient for experiments where all
# site data is available in one process.
# =============================================================================

#' Full CLEAR pipeline for simulation.
#'
#' @param site_data_list  List of site data.frames (shared column names).
#' @param K_comp          Components per site GMM.
#' @param n_0             Reference sample size.
#' @param inflate         Optional variance inflation (default 1 = off).
#' @return List with Site_h stats for each site, plus Z_0.
CLEAR_sim <- function(site_data_list, K_comp = 3, n_0 = 10000, inflate = 1) {
  H <- length(site_data_list)
  
  # Step 1: each site fits and "shares" its GMM.
  site_models <- lapply(site_data_list, site_fit_gmm, K_comp = K_comp)
  
  # Step 2: center builds the reference distribution.
  Z_0 <- center_make_reference(site_models, n_0 = n_0, inflate = inflate)
  
  # Steps 3 & 4: per-site IPW weights + weighted statistics.
  results_list <- list()
  for (h in seq_len(H)) {
    w_norm <- estimate_ipw_weights(site_data_list[[h]], Z_0)
    results_list[[paste0("Site_", h)]] <- weighted_eda_stats(Z_0, w_norm)
  }
  
  results_list$Z_0 <- Z_0
  results_list
}


# =============================================================================
# VALIDATION HELPER (simulation use)
# Exact pooled-truth statistics, using population (div-by-n) moments to match
# the mathematical form of the IPW estimates.
# =============================================================================

calc_empirical_truth <- function(df) {
  n     <- nrow(df)
  means <- colMeans(df)
  vars  <- apply(df, 2, function(x) sum((x - mean(x))^2) / n)
  
  q_probs <- c(0.05, 0.25, 0.50, 0.75, 0.95)
  quant   <- t(apply(df, 2, quantile, probs = q_probs))
  colnames(quant) <- c("5%", "25%", "50%", "75%", "95%")
  
  cov_mat <- cov(df) * ((n - 1) / n)
  
  list(Mean = means, Variance = vars, Quantiles = quant, Covariance = cov_mat)
}