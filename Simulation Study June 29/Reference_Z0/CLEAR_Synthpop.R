# =============================================================================
# CLEAR + synthpop  (hybrid federated method)
# -----------------------------------------------------------------------------
# IDEA
#   Keep CLEAR's analysis-center machinery -- density-ratio IPW weighting and
#   weighted EDA -- but BUILD THE REFERENCE DISTRIBUTION with synthpop instead
#   of by sampling from shared GMM parameters.
#
#   Original CLEAR reference:  each site shares GMM params (pi, mu, Sigma);
#                              center samples Z_0 from the pooled mixture.
#   Hybrid reference:          each site fits a synthpop CART model and emits a
#                              synthetic sample; the pooled synthetic sets ARE
#                              the reference Z_0.
#
#   Everything after the reference is IDENTICAL to CLEAR:
#     for each site h -> estimate_ipw_weights(df_h, Z_0) -> weighted_eda_stats().
#
# WHY THIS IS A VALID SWAP
#   CLEAR's downstream stages only require a reference data.frame Z_0 with the
#   shared column names. They are agnostic to how Z_0 was generated. So we
#   replace center_make_reference() with a synthpop-based builder and reuse
#   estimate_ipw_weights() and weighted_eda_stats() verbatim.
#
# PRIVACY NOTE
#   In CLEAR, only GMM parameters leave a site. In this hybrid, each site emits
#   a SYNTHETIC SAMPLE (not raw data). Synthetic microdata is a weaker privacy
#   guarantee than summary parameters, but still avoids sharing real records.
#   Choose per your threat model.
#
# Requires: CLEAR.R (for estimate_ipw_weights, weighted_eda_stats,
#           calc_empirical_truth) and the synthpop package.
# =============================================================================

library(MASS)
library(mclust)
library(synthpop)

# CLEAR.R provides: site_fit_gmm, center_make_reference, estimate_ipw_weights,
#                   weighted_eda_stats, CLEAR_sim, calc_empirical_truth.
source("CLEAR.R")   


# =============================================================================
# SITE-SIDE (hybrid): synthesize a local sample with synthpop.
# This REPLACES site_fit_gmm as the shareable object. What leaves the site is a
# synthetic data.frame (n_syn x d), not GMM parameters.
# =============================================================================

#' Fit a local synthpop CART model and return a synthetic sample.
#'
#' @param df        A site's data: numeric data.frame, rows = obs, named cols.
#' @param n_syn     Number of synthetic rows to emit. Default = nrow(df), which
#'                  preserves the site's sample-size weight when sets are pooled.
#' @param method    synthpop synthesis method (default "cart").
#' @param seed      Optional integer seed for reproducible synthesis.
#' @return A list with: syn (synthetic data.frame), var_names, n (real n),
#'         n_syn (synthetic n). The `syn` element is what "leaves" the site.
site_synthesize <- function(df, n_syn = NULL, method = "cart", seed = NULL) {
  var_names <- colnames(df)
  if (is.null(var_names)) stop("`df` must have column names (variable names).")
  if (is.null(n_syn)) n_syn <- nrow(df)
  if (is.null(seed))  seed  <- sample.int(.Machine$integer.max, 1)
  
  s <- syn(df, method = method, m = 1, k = n_syn, seed = seed, print.flag = FALSE)
  syn_df <- s$syn
  # synthpop preserves column names/order; guard anyway.
  syn_df <- syn_df[, var_names, drop = FALSE]
  
  list(syn = syn_df, var_names = var_names, n = nrow(df), n_syn = n_syn)
}


# =============================================================================
# ANALYSIS-CENTER (hybrid): build the reference by pooling synthetic samples.
# Mirrors center_make_reference()'s output contract: a data.frame Z_0 (n_0 x d).
# =============================================================================

#' Build a global reference distribution from pooled site synthetic samples.
#'
#' Sites are represented in proportion to their emitted synthetic n. If you set
#' each site's n_syn = its real n (the default in site_synthesize), the pooled
#' reference reflects the real site-size mixing, matching CLEAR's n_h weighting.
#'
#' Optionally rescale total reference size to a target n_0 while preserving the
#' per-site proportions (so this hybrid can be swept over n_0 like CLEAR).
#'
#' @param syn_models   List of objects returned by site_synthesize().
#' @param n_0          Optional target total reference size. If NULL, use the
#'                     pooled synthetic sets as-is (sum of each site's n_syn).
#' @return A data.frame Z_0 (n_0 x d) with the shared variable names.
center_make_reference_syn <- function(syn_models, n_0 = NULL) {
  var_names <- syn_models[[1]]$var_names
  H         <- length(syn_models)
  
  # Pool the synthetic sets directly.
  Z_pool <- do.call(rbind, lapply(syn_models, function(m) m$syn[, var_names, drop = FALSE]))
  
  if (is.null(n_0)) {
    Z_0 <- Z_pool
  } else {
    # Resample to a target size n_0, preserving each site's proportion
    # (proportional allocation by emitted synthetic n).
    n_syn_vec <- sapply(syn_models, function(m) m$n_syn)
    prob_h    <- n_syn_vec / sum(n_syn_vec)
    alloc     <- as.vector(rmultinom(1, size = n_0, prob = prob_h))
    
    parts <- lapply(seq_len(H), function(h) {
      src <- syn_models[[h]]$syn[, var_names, drop = FALSE]
      if (alloc[h] == 0) return(src[0, , drop = FALSE])
      idx <- sample(seq_len(nrow(src)), alloc[h], replace = TRUE)
      src[idx, , drop = FALSE]
    })
    Z_0 <- do.call(rbind, parts)
  }
  
  Z_0 <- as.data.frame(Z_0)
  colnames(Z_0) <- var_names
  rownames(Z_0) <- NULL
  Z_0
}


# =============================================================================
# SIMULATION WRAPPER (hybrid): full CLEAR+synthpop pipeline.
# Same return contract as CLEAR_sim(): a list of Site_h weighted stats + Z_0.
# =============================================================================

#' Full CLEAR+synthpop pipeline for simulation.
#'
#' Step 1 (site)   : each site synthesizes locally via synthpop CART.
#' Step 2 (center) : pool synthetic sets into the reference Z_0.
#' Steps 3 & 4     : ORIGINAL CLEAR -- per-site density-ratio IPW weights on
#'                   Z_0, then weighted EDA statistics.
#'
#' @param site_data_list  List of site data.frames (shared column names).
#' @param n_0             Optional target reference size. NULL = use pooled
#'                        synthetic sets as-is (each site emits its real n).
#' @param n_syn_each      Optional vector/length-1 of synthetic n per site.
#'                        NULL = each site emits its own real n.
#' @param method          synthpop method (default "cart").
#' @param base_seed       Optional base seed; site h uses base_seed + h.
#' @return List with Site_h stats for each site, plus Z_0.
CLEAR_synthpop_sim <- function(site_data_list, n_0 = NULL, n_syn_each = NULL,
                               method = "cart", base_seed = NULL) {
  H <- length(site_data_list)
  
  # Step 1: each site synthesizes and "shares" a synthetic sample.
  syn_models <- lapply(seq_len(H), function(h) {
    n_syn_h <- if (is.null(n_syn_each)) NULL
    else if (length(n_syn_each) == 1) n_syn_each
    else n_syn_each[h]
    seed_h  <- if (is.null(base_seed)) NULL else base_seed + h
    site_synthesize(site_data_list[[h]], n_syn = n_syn_h, method = method, seed = seed_h)
  })
  
  # Step 2: center builds the reference distribution from pooled synthetics.
  Z_0 <- center_make_reference_syn(syn_models, n_0 = n_0)
  
  # Steps 3 & 4: ORIGINAL CLEAR IPW + weighted EDA, unchanged.
  results_list <- list()
  for (h in seq_len(H)) {
    w_norm <- estimate_ipw_weights(site_data_list[[h]], Z_0)   # from CLEAR.R
    results_list[[paste0("Site_", h)]] <- weighted_eda_stats(Z_0, w_norm)  # from CLEAR.R
  }
  
  results_list$Z_0 <- Z_0
  results_list
}


# =============================================================================
# CONVENIENCE: pool per-site weighted stats into GLOBAL estimates.
# Same n_h-weighted mixing used in the CLEAR simulation scripts, provided here
# so the hybrid is self-contained.
# =============================================================================

#' Combine per-site CLEAR(+synthpop) stats into global pooled estimates.
#'
#' @param clear_out  Output of CLEAR_synthpop_sim() (or CLEAR_sim()).
#' @param n_h        Vector of real site sample sizes (for n_h weighting).
#' @return List: Mean, Variance, Covariance, Quantiles, w_global.
combine_global_stats <- function(clear_out, n_h) {
  H         <- length(n_h)
  site_keys <- paste0("Site_", seq_len(H))
  w_site    <- n_h / sum(n_h)
  Z_0       <- clear_out$Z_0
  var_names <- colnames(Z_0)
  D         <- length(var_names)
  q_probs   <- c(0.05, 0.25, 0.50, 0.75, 0.95)
  
  w_global <- Reduce(`+`, Map(function(k, w) w * clear_out[[k]]$Weights, site_keys, w_site))
  w_global <- w_global / sum(w_global)
  
  g_mean <- Reduce(`+`, Map(function(k, w) w * clear_out[[k]]$Mean, site_keys, w_site))
  g_ex2  <- Reduce(`+`, Map(function(k, w) w * (clear_out[[k]]$Variance + clear_out[[k]]$Mean^2),
                            site_keys, w_site))
  g_var  <- g_ex2 - g_mean^2
  g_ecc  <- Reduce(`+`, Map(function(k, w) {
    mu_h <- clear_out[[k]]$Mean; w * (clear_out[[k]]$Covariance + outer(mu_h, mu_h))
  }, site_keys, w_site))
  g_cov  <- g_ecc - outer(g_mean, g_mean)
  
  g_quant <- t(sapply(seq_len(D), function(v) {
    x <- Z_0[, v]; ord <- order(x); xo <- x[ord]; cw <- cumsum(w_global[ord])
    sapply(q_probs, function(p) xo[which(cw >= p)[1]])
  }))
  dimnames(g_quant) <- list(var_names, c("5%", "25%", "50%", "75%", "95%"))
  
  list(Mean = g_mean, Variance = g_var, Covariance = g_cov,
       Quantiles = g_quant, w_global = w_global)
}


# =============================================================================
# QUICK DEMO  (uncomment to run; needs CLEAR.R sourced)
# -----------------------------------------------------------------------------
# source("CLEAR.R")
# set.seed(1)
# mk <- function(n, shift) { X <- MASS::mvrnorm(n, c(0,0,0)+shift, diag(3))
#                            df <- as.data.frame(X); colnames(df) <- c("X1","X2","X3"); df }
# sites <- list(mk(400, 0), mk(600, 1), mk(300, -1))
# n_h   <- sapply(sites, nrow)
# pooled <- do.call(rbind, sites)
#
# truth <- calc_empirical_truth(pooled)
# hyb   <- CLEAR_synthpop_sim(sites, n_0 = 10000, base_seed = 100)
# gh    <- combine_global_stats(hyb, n_h)
#
# print(round(rbind(truth = truth$Mean, hybrid = gh$Mean), 4))
# print(round(rbind(truth = truth$Variance, hybrid = gh$Variance), 4))