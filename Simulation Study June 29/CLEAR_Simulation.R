# =============================================================================
# CLEAR simulation study: CLEAR_sim vs. calc_empirical_truth
# -----------------------------------------------------------------------------
# Setup:
#   * 5 sites, 3 dimensions, non-Gaussian (per-site Gaussian MIXTURE) truth.
#   * Pooled truth is multimodal/skewed -> a real stress test for the
#     GMM-reference + IPW recovery.
#   * 500 iterations run in parallel over parallel::detectCores() workers.
#   * Per iteration: compute CLEAR estimates AND pooled empirical truth on the
#     SAME data, store the (estimate - truth) errors for every statistic.
#   * Output: bias + RMSE tables per statistic, plus boxplots of error.
#
# Requires CLEAR_v2.R (site_fit_gmm, center_make_reference,
# estimate_ipw_weights, weighted_eda_stats, CLEAR_sim, calc_empirical_truth).
# =============================================================================

library(parallel)
library(MASS)
library(mclust)
library(ggplot2)
library(reshape2)

# Adjust path as needed.
source("CLEAR.R")


# -----------------------------------------------------------------------------
# 0. Global configuration
# -----------------------------------------------------------------------------
N_ITER   <- 500
N_SITES  <- 5
D        <- 3
VAR_NMS  <- c("X1", "X2", "X3")
K_COMP   <- 3        # components requested per site GMM in CLEAR
N_0      <- 10000    # reference sample size
Q_LABS   <- c("5%", "25%", "50%", "75%", "95%")

# Per-site sample sizes (heterogeneous).
SITE_N <- c(800, 1200, 600, 1000, 900)


# -----------------------------------------------------------------------------
# 1. Data-generating process: non-Gaussian, heterogeneous per-site mixtures
# -----------------------------------------------------------------------------
# Each site is a 2-component Gaussian mixture in 3D. Component 2 is shifted and
# given a different covariance to induce skew/multimodality. Sites differ in
# means, covariances, and mixing weights so the pooled truth is messy.

make_site_params <- function() {
  base_means <- list(
    c( 0,  0,  0), c( 2, -1,  1), c(-2,  1, -1),
    c( 1,  2,  0), c(-1, -2,  2)
  )
  base_cov <- function(s, rho) {
    M <- matrix(rho, D, D); diag(M) <- 1
    (s^2) * M
  }
  lapply(seq_len(N_SITES), function(h) {
    mu1 <- base_means[[h]]
    mu2 <- base_means[[h]] + c(3, 3, -2)         # shifted 2nd mode
    list(
      pi    = c(0.65, 0.35),
      mu    = list(mu1, mu2),
      Sigma = list(base_cov(1.0, 0.3 + 0.05 * h),
                   base_cov(1.6, -0.2)),         # different shape -> skew
      n     = SITE_N[h]
    )
  })
}

SITE_PARAMS <- make_site_params()

#' Draw one site's data from its 2-component mixture.
draw_site <- function(p) {
  comp <- sample(1:2, p$n, replace = TRUE, prob = p$pi)
  X <- matrix(0, p$n, D)
  for (k in 1:2) {
    idx <- which(comp == k)
    if (length(idx) > 0) {
      X[idx, ] <- mvrnorm(length(idx), mu = p$mu[[k]], Sigma = p$Sigma[[k]])
    }
  }
  df <- as.data.frame(X); colnames(df) <- VAR_NMS; df
}


# # Check dist
# site_data <- lapply(SITE_PARAMS, draw_site)
# # plot density plot of each variable for all sites
# par(mfrow = c(5, 3))
# for (h in seq_len(N_SITES)) {
#   for (v in VAR_NMS) {
#     plot(density(site_data[[h]][[v]]), main = paste("Site", h, v), xlab = v)
#   }
# }

# -----------------------------------------------------------------------------
# 2. One iteration: returns error vectors (estimate - truth) per statistic
# -----------------------------------------------------------------------------
# CLEAR yields PER-SITE statistics; the pooled empirical truth is one global
# target. To compare on equal footing we aggregate CLEAR's per-site estimates
# into a global estimate using the same n_h weighting the reference uses:
#   global_mean = sum_h (n_h/N) * mean_h, etc. Variance/covariance pooled via
# the law of total variance so the comparison is apples-to-apples.

run_one_iter <- function(iter_seed) {
  set.seed(iter_seed)
  
  # Generate this iteration's site data.
  site_data <- lapply(SITE_PARAMS, draw_site)
  pooled    <- do.call(rbind, site_data)
  n_h       <- sapply(site_data, nrow)
  N         <- sum(n_h)
  w_site    <- n_h / N
  
  # --- Truth on pooled raw data ---
  truth <- calc_empirical_truth(pooled)
  
  # --- CLEAR estimates (per site) ---
  clear <- tryCatch(
    CLEAR_sim(site_data, K_comp = K_COMP, n_0 = N_0, inflate = 1),
    error = function(e) NULL
  )
  if (is.null(clear)) return(NULL)
  
  site_keys <- paste0("Site_", seq_len(N_SITES))
  
  # Aggregate per-site CLEAR stats -> global estimate.
  # Mean: weighted average.
  est_mean <- Reduce(`+`, Map(function(k, w) w * clear[[k]]$Mean,
                              site_keys, w_site))
  
  # Variance (law of total variance):
  #   Var = sum_h w_h * (Var_h + Mean_h^2) - GlobalMean^2
  ex2 <- Reduce(`+`, Map(function(k, w) w * (clear[[k]]$Variance + clear[[k]]$Mean^2),
                         site_keys, w_site))
  est_var <- ex2 - est_mean^2
  
  # Covariance (law of total covariance):
  #   Cov = sum_h w_h * (Cov_h + mu_h mu_h^T) - GlobalMean GlobalMean^T
  ecc <- Reduce(`+`, Map(function(k, w) {
    mu_h <- clear[[k]]$Mean
    w * (clear[[k]]$Covariance + outer(mu_h, mu_h))
  }, site_keys, w_site))
  est_cov <- ecc - outer(est_mean, est_mean)
  
  # Quantiles: per-site weighted quantiles don't pool linearly. Recover a
  # GLOBAL weighted quantile directly on the shared reference Z_0, using the
  # n_h-weighted mixture of per-site IPW weights as the global weight vector.
  Z_0 <- clear$Z_0
  w_global <- Reduce(`+`, Map(function(k, w) w * clear[[k]]$Weights,
                              site_keys, w_site))
  w_global <- w_global / sum(w_global)
  est_quant <- t(sapply(seq_len(D), function(v) {
    x <- Z_0[, v]; ord <- order(x); xo <- x[ord]; cw <- cumsum(w_global[ord])
    sapply(c(0.05, 0.25, 0.50, 0.75, 0.95), function(p) xo[which(cw >= p)[1]])
  }))
  
  # --- Errors (estimate - truth) ---
  cov_lt <- lower.tri(truth$Covariance, diag = TRUE)
  list(
    mean_err = est_mean - truth$Mean,
    var_err  = est_var  - truth$Variance,
    quant_err = as.vector(est_quant - truth$Quantiles),
    cov_err  = est_cov[cov_lt] - truth$Covariance[cov_lt]
  )
}


# -----------------------------------------------------------------------------
# 3. Run 500 iterations in parallel
# -----------------------------------------------------------------------------
n_cores <- parallelly::availableCores()
message(sprintf("Running %d iterations on %d cores...", N_ITER, n_cores))

cl <- makeCluster(n_cores)
on.exit(stopCluster(cl), add = TRUE)

# Export everything workers need.
clusterExport(cl, c("SITE_PARAMS", "N_SITES", "D", "VAR_NMS", "K_COMP", "N_0",
                    "SITE_N", "draw_site", "run_one_iter", "calc_empirical_truth",
                    "site_fit_gmm", "center_make_reference", "estimate_ipw_weights",
                    "weighted_eda_stats", "CLEAR_sim"))
clusterEvalQ(cl, { library(MASS); library(mclust) })

seeds <- seq_len(N_ITER) + 1000L
results <- parLapply(cl, seeds, run_one_iter)
stopCluster(cl)

# Drop failed iterations.
ok <- !sapply(results, is.null)
message(sprintf("%d / %d iterations succeeded.", sum(ok), N_ITER))
results <- results[ok]


# -----------------------------------------------------------------------------
# 4. Aggregate: bias + RMSE tables per statistic
# -----------------------------------------------------------------------------
stack_errs <- function(field, col_names) {
  M <- t(sapply(results, function(r) r[[field]]))
  colnames(M) <- col_names
  M
}

# Column labels.
mean_lab  <- VAR_NMS
var_lab   <- VAR_NMS
quant_lab <- as.vector(outer(VAR_NMS, Q_LABS, paste, sep = "_"))
cov_idx   <- which(lower.tri(matrix(0, D, D), diag = TRUE), arr.ind = TRUE)
cov_lab   <- apply(cov_idx, 1, function(ij) paste(VAR_NMS[ij[1]], VAR_NMS[ij[2]], sep = "-"))

mean_M  <- stack_errs("mean_err",  mean_lab)
var_M   <- stack_errs("var_err",   var_lab)
quant_M <- stack_errs("quant_err", quant_lab)
cov_M   <- stack_errs("cov_err",   cov_lab)

summ_tbl <- function(M, stat_name) {
  data.frame(
    Statistic = stat_name,
    Quantity  = colnames(M),
    Bias      = colMeans(M),
    RMSE      = sqrt(colMeans(M^2)),
    SD        = apply(M, 2, sd),
    row.names = NULL
  )
}

summary_table <- rbind(
  summ_tbl(mean_M,  "Mean"),
  summ_tbl(var_M,   "Variance"),
  summ_tbl(quant_M, "Quantile"),
  summ_tbl(cov_M,   "Covariance")
)

cat("\n=========== CLEAR vs. Empirical Truth: Bias / RMSE ===========\n")
print(summary_table, digits = 4, row.names = FALSE)
write.csv(summary_table, "CLEAR_sim_summary.csv", row.names = FALSE)


# -----------------------------------------------------------------------------
# 5. Plots: boxplots of error by statistic
# -----------------------------------------------------------------------------
make_long <- function(M, stat_name) {
  df <- as.data.frame(M)
  long <- melt(df, variable.name = "Quantity", value.name = "Error")
  long$Statistic <- stat_name
  long
}

long_all <- rbind(
  make_long(mean_M,  "Mean"),
  make_long(var_M,   "Variance"),
  make_long(quant_M, "Quantile"),
  make_long(cov_M,   "Covariance")
)

p <- ggplot(long_all, aes(x = Quantity, y = Error)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_boxplot(outlier.size = 0.5, fill = "grey85") +
  facet_wrap(~ Statistic, scales = "free", ncol = 2) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "CLEAR estimation error vs. empirical truth (500 iterations)",
       subtitle = "Non-Gaussian mixture DGP, 5 sites x 3 dims",
       y = "Estimate - Truth", x = NULL)

ggsave("CLEAR_sim_boxplots.png", p, width = 11, height = 8, dpi = 150)
cat("\nSaved: CLEAR_sim_summary.csv, CLEAR_sim_boxplots.png\n")
print(p)