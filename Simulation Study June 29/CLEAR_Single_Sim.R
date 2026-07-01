# =============================================================================
# CLEAR single-run
# -----------------------------------------------------------------------------
# The DGP is a 4-component Gaussian mixture (shared proportions across sites).
# Can choose the 4-vector PROPORTIONS below. Components with weight 0 are
# dropped from data generation, so we can mimic fewer modes:
#
# Four comparison blocks (same as the prior diagnostic):
#   (A) Descriptive statistics vs. pooled truth
#   (B) Per-site density recovery
#   (C) Global density recovery
#   (D) Lasso: pooled raw vs. CLEAR weighted reference (X1 ~ X2 + X3)
# =============================================================================
rm(list = ls())
library(MASS)
library(mclust)
library(ggplot2)
library(reshape2)
library(glmnet)

source("CLEAR.R")

# set.seed(2026)

# -----------------------------------------------------------------------------
# 0. CONFIG -- edit PROPORTIONS to change the DGP shape
# -----------------------------------------------------------------------------
PROPORTIONS <- c(0.5, 0.4, 0.1, 0)   

stopifnot(length(PROPORTIONS) == 4, all(PROPORTIONS >= 0), sum(PROPORTIONS) > 0)
PROPORTIONS <- PROPORTIONS / sum(PROPORTIONS)   # normalize

N_SITES <- 5
D       <- 3
VAR_NMS <- c("X1", "X2", "X3")
K_COMP  <- 4                 
N_0     <- 10000
SITE_N  <- seq(200, 800, length.out = N_SITES)  # per-site sample sizes
Q_PROBS <- c(0.05, 0.25, 0.50, 0.75, 0.95)
Q_LABS  <- c("5%", "25%", "50%", "75%", "95%")


# -----------------------------------------------------------------------------
# 1. DGP: shared proportions across sites
# -----------------------------------------------------------------------------
# 4 component centers + covariances per site. Sites differ by a per-site location offset so the pooled distribution is heterogeneous.
COMP_CENTERS <- list(c( 0,  0,  0), c( 4,  4, -3),
                     c(-4,  3,  3), c( 3, -4,  4))
comp_cov <- function(s, rho) { M <- matrix(rho, D, D); diag(M) <- 1; (s^2) * M }
COMP_COVS <- list(comp_cov(1.0,  0.3), comp_cov(1.4, -0.2),
                  comp_cov(1.1,  0.1), comp_cov(1.3, -0.3))

site_offset <- function(h) c((h - 3) * 1.5, (h - 3) * -1.0, (h - 3) * 0.8)

draw_site <- function(h, n) {
  off  <- site_offset(h)
  comp <- sample(seq_len(4), n, replace = TRUE, prob = PROPORTIONS)
  X <- matrix(0, n, D)
  for (k in unique(comp)) {
    idx <- which(comp == k)
    X[idx, ] <- mvrnorm(length(idx), mu = COMP_CENTERS[[k]] + off, Sigma = COMP_COVS[[k]])
  }
  df <- as.data.frame(X); colnames(df) <- VAR_NMS; df
}


# -----------------------------------------------------------------------------
# 1b. Generate ONCE, run CLEAR ONCE
# -----------------------------------------------------------------------------
site_data <- Map(draw_site, seq_len(N_SITES), SITE_N)
pooled    <- do.call(rbind, site_data)
n_h       <- sapply(site_data, nrow)
N         <- sum(n_h)
w_site    <- n_h / N


clear <- CLEAR_sim(site_data, K_comp = K_COMP, n_0 = N_0, inflate = 1)
Z_0   <- clear$Z_0
site_keys <- paste0("Site_", seq_len(N_SITES))

w_global <- Reduce(`+`, Map(function(k, w) w * clear[[k]]$Weights, site_keys, w_site)) # barW
w_global <- w_global / sum(w_global)



# =============================================================================
# (A) DESCRIPTIVE STATISTICS
# =============================================================================
truth <- calc_empirical_truth(pooled)

est_mean <- Reduce(`+`, Map(function(k, w) w * clear[[k]]$Mean, site_keys, w_site))
ex2 <- Reduce(`+`, Map(function(k, w) w * (clear[[k]]$Variance + clear[[k]]$Mean^2),
                       site_keys, w_site))
est_var <- ex2 - est_mean^2
ecc <- Reduce(`+`, Map(function(k, w) {
  mu_h <- clear[[k]]$Mean; w * (clear[[k]]$Covariance + outer(mu_h, mu_h))
}, site_keys, w_site))
est_cov <- ecc - outer(est_mean, est_mean)

est_quant <- t(sapply(seq_len(D), function(v) {
  x <- Z_0[, v]; ord <- order(x); xo <- x[ord]; cw <- cumsum(w_global[ord])
  sapply(Q_PROBS, function(p) xo[which(cw >= p)[1]])
}))
dimnames(est_quant) <- list(VAR_NMS, Q_LABS)

cat("\n================= (A) DESCRIPTIVE STATISTICS =================\n")
cat("\n--- Mean ---\n")
print(round(data.frame(CLEAR = est_mean, Truth = truth$Mean, Diff = est_mean - truth$Mean), 4))
cat("\n--- Variance ---\n")
print(round(data.frame(CLEAR = est_var, Truth = truth$Variance, Diff = est_var - truth$Variance), 4))
cat("\n--- Quantiles (CLEAR) ---\n"); print(round(est_quant, 4))
cat("\n--- Quantiles (Truth) ---\n"); print(round(truth$Quantiles, 4))
cat("\n--- Covariance (CLEAR) ---\n"); print(round(est_cov, 4))
cat("\n--- Covariance (Truth) ---\n"); print(round(truth$Covariance, 4))


# =============================================================================
# (B) PER-SITE DENSITY RECOVERY
# =============================================================================
weighted_density_df <- function(z, w, var_name, label, n_draw = 20000) {
  idx <- sample(seq_along(z), n_draw, replace = TRUE, prob = w)
  data.frame(value = z[idx], Variable = var_name, Source = label)
}

for (h in seq_len(N_SITES)) {
  w_h <- clear[[site_keys[h]]]$Weights
  long_h <- do.call(rbind, lapply(seq_len(D), function(v) {
    rbind(
      data.frame(value = site_data[[h]][, v], Variable = VAR_NMS[v], Source = "Raw site data"),
      weighted_density_df(Z_0[, v], w_h, VAR_NMS[v], "CLEAR weighted ref")
    )
  }))
  ph <- ggplot(long_h, aes(value, color = Source, fill = Source)) +
    geom_density(alpha = 0.15) +
    facet_wrap(~ Variable, scales = "free", ncol = 3) +
    theme_bw(base_size = 10) +
    labs(title = paste0("Site ", h, ": density"), x = NULL, y = "Density")
  print(ph)
  # ggsave(sprintf("CLEAR_density_site%d.png", h), ph, width = 10, height = 3.2, dpi = 150)
}


# =============================================================================
# (C) GLOBAL DENSITY RECOVERY
# =============================================================================
long_global <- do.call(rbind, lapply(seq_len(D), function(v) {
  rbind(
    data.frame(value = pooled[, v], Variable = VAR_NMS[v], Source = "Pooled raw data"),
    weighted_density_df(Z_0[, v], w_global, VAR_NMS[v], "CLEAR global ref")
  )
}))
p_global <- ggplot(long_global, aes(value, color = Source, fill = Source)) +
  geom_density(alpha = 0.15) +
  facet_wrap(~ Variable, scales = "free", ncol = 3) +
  theme_bw(base_size = 11) +
  labs(title = "Global density recovery: pooled raw vs. CLEAR weighted reference",
       subtitle = sprintf("Proportions = (%s)", paste(round(PROPORTIONS, 2), collapse = ", ")),
       x = NULL, y = "Density")
print(p_global)
# ggsave("CLEAR_density_global.png", p_global, width = 11, height = 3.8, dpi = 150)




# =============================================================================
# (B2) PER-SITE JOINT DENSITY RECOVERY (pairwise 2D heatmaps)
# -----------------------------------------------------------------------------
# Filled 2D KDE heatmaps. Two rows: true site data (top) vs. weight-resampled
# CLEAR reference (bottom). Columns = variable pairs (X1-X2, X1-X3, X2-X3).
# =============================================================================

var_pairs <- combn(VAR_NMS, 2, simplify = FALSE)

# Weighted resample of the full reference (preserves joint structure).
weighted_resample <- function(Z, w, n_draw = 20000) {
  idx <- sample(seq_len(nrow(Z)), n_draw, replace = TRUE, prob = w)
  Z[idx, , drop = FALSE]
}

for (h in seq_len(N_SITES)) {
  w_h    <- clear[[site_keys[h]]]$Weights
  ref_rs <- weighted_resample(Z_0, w_h)
  
  long_pairs <- do.call(rbind, lapply(var_pairs, function(pr) {
    lab <- paste(pr[1], pr[2], sep = " vs ")
    rbind(
      data.frame(xv = site_data[[h]][, pr[1]], yv = site_data[[h]][, pr[2]],
                 Pair = lab, Source = "Raw site data"),
      data.frame(xv = ref_rs[, pr[1]], yv = ref_rs[, pr[2]],
                 Pair = lab, Source = "CLEAR weighted ref")
    )
  }))
  long_pairs$Source <- factor(long_pairs$Source,
                              levels = c("Raw site data", "CLEAR weighted ref"))
  
  ph_joint <- ggplot(long_pairs, aes(xv, yv)) +
    stat_density_2d(aes(fill = after_stat(density)),
                    geom = "raster", contour = FALSE) +
    facet_grid(Source ~ Pair, scales = "free") +   # 2 rows x pairs columns
    scale_fill_viridis_c(option = "magma") +
    theme_bw(base_size = 10) +
    theme(legend.position = "right") +
    labs(title = paste0("Site ", h, ": joint density heatmap"),
         x = NULL, y = NULL, fill = "Density")
  print(ph_joint)
  # ggsave(sprintf("CLEAR_joint_site%d.png", h), ph_joint, width = 10, height = 6, dpi = 150)
}


# =============================================================================
# (C2) GLOBAL JOINT DENSITY RECOVERY (pairwise 2D heatmaps)
# -----------------------------------------------------------------------------
# Pooled raw data (top row) vs. reference resampled by n_h-mixed global weights
# (bottom row).
# =============================================================================
ref_global_rs <- weighted_resample(Z_0, w_global)

long_pairs_global <- do.call(rbind, lapply(var_pairs, function(pr) {
  lab <- paste(pr[1], pr[2], sep = " vs ")
  rbind(
    data.frame(xv = pooled[, pr[1]], yv = pooled[, pr[2]],
               Pair = lab, Source = "Pooled raw data"),
    data.frame(xv = ref_global_rs[, pr[1]], yv = ref_global_rs[, pr[2]],
               Pair = lab, Source = "CLEAR global ref")
  )
}))
long_pairs_global$Source <- factor(long_pairs_global$Source,
                                   levels = c("Pooled raw data", "CLEAR global ref"))

p_joint_global <- ggplot(long_pairs_global, aes(xv, yv)) +
  stat_density_2d(aes(fill = after_stat(density)),
                  geom = "raster", contour = FALSE) +
  facet_grid(Source ~ Pair, scales = "free") +     # 2 rows x pairs columns
  scale_fill_viridis_c(option = "magma") +
  theme_bw(base_size = 11) +
  theme(legend.position = "right") +
  labs(title = "Global joint density heatmap: pooled raw vs. CLEAR weighted reference",
       subtitle = sprintf("Proportions = (%s)", paste(round(PROPORTIONS, 2), collapse = ", ")),
       x = NULL, y = NULL, fill = "Density")
print(p_joint_global)
# ggsave("CLEAR_joint_global.png", p_joint_global, width = 11, height = 6, dpi = 150)



# =============================================================================
# (D) LASSO: pooled raw vs. CLEAR weighted reference (X1 ~ X2 + X3)
# =============================================================================
outcome <- "X1"; predictors <- c("X2", "X3")

Xp <- as.matrix(pooled[, predictors]); yp <- pooled[[outcome]]
set.seed(1); cv_pool <- cv.glmnet(Xp, yp, alpha = 1, nfolds = 10)
coef_pool <- as.vector(coef(cv_pool, s = "lambda.min"))

Xr <- as.matrix(Z_0[, predictors]); yr <- Z_0[[outcome]]
set.seed(1); cv_ref <- cv.glmnet(Xr, yr, alpha = 1, weights = w_global, nfolds = 10)
coef_ref <- as.vector(coef(cv_ref, s = "lambda.min"))

lasso_tbl <- data.frame(
  Term           = c("(Intercept)", predictors),
  Pooled_raw     = round(coef_pool, 4),
  CLEAR_weighted = round(coef_ref, 4),
  Diff           = round(coef_ref - coef_pool, 4)
)
cat("\n================= (D) LASSO COEFFICIENTS =================\n")
cat(sprintf("Model: %s ~ %s  (lambda.min)\n", outcome, paste(predictors, collapse = " + ")))
print(lasso_tbl, row.names = FALSE)
write.csv(lasso_tbl, "CLEAR4_lasso_comparison.csv", row.names = FALSE)

path_df <- function(fit, label) {
  B <- as.matrix(fit$glmnet.fit$beta)
  data.frame(loglambda = rep(log(fit$glmnet.fit$lambda), each = nrow(B)),
             Coef = rep(rownames(B), times = ncol(B)),
             Value = as.vector(B), Source = label)
}
paths <- rbind(path_df(cv_pool, "Pooled raw"), path_df(cv_ref, "CLEAR weighted"))
p_lasso <- ggplot(paths, aes(loglambda, Value, color = Coef, linetype = Source)) +
  geom_line(linewidth = 0.7) + theme_bw(base_size = 11) +
  labs(title = paste0("Lasso paths: ", outcome, " ~ ", paste(predictors, collapse = " + ")),
       x = "log(lambda)", y = "Coefficient")
# ggsave("CLEAR_lasso_paths.png", p_lasso, width = 9, height = 5, dpi = 150)
print(p_lasso)

