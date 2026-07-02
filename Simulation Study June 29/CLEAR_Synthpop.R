# =============================================================================
# CLEAR  vs.  federated synthpop  vs.  TRUTH
# -----------------------------------------------------------------------------
# Same multi-site Gaussian-mixture DGP as CLEAR_Single_Sim.R. We add an EXISTING
# synthetic-data method -- synthpop::syn() (CART) -- as a competitor to CLEAR,
# in the FEDERATED/DISTRIBUTED setting only: each site synthesizes LOCALLY and
# the synthetic sets are pooled (no raw data ever leaves a site), mirroring
# CLEAR's privacy constraint. Everything is scored against pooled TRUTH.
#
#   (1) TRUTH   -- pooled empirical statistics (the target)
#   (2) CLEAR   -- federated GMM-share + density-ratio IPW reweighting
#   (3) synFED  -- federated synthpop: per-site synthesis, synthetic sets pooled
#

# =============================================================================
rm(list = ls())
library(MASS)
library(mclust)
library(ggplot2)
library(reshape2)
library(glmnet)
library(synthpop)

source("CLEAR.R")

# set.seed(2026)

# -----------------------------------------------------------------------------
# 0. CONFIG  (identical to CLEAR_Single_Sim.R)
# -----------------------------------------------------------------------------
PROPORTIONS <- c(0.5, 0.4, 0.1, 0)
stopifnot(length(PROPORTIONS) == 4, all(PROPORTIONS >= 0), sum(PROPORTIONS) > 0)
PROPORTIONS <- PROPORTIONS / sum(PROPORTIONS)

N_SITES <- 5
D       <- 3
VAR_NMS <- c("X1", "X2", "X3")
K_COMP  <- 4
N_0     <- 10000
SITE_N  <- seq(200, 800, length.out = N_SITES)
Q_PROBS <- c(0.05, 0.25, 0.50, 0.75, 0.95)
Q_LABS  <- c("5%", "25%", "50%", "75%", "95%")


# -----------------------------------------------------------------------------
# 1. DGP  (identical to CLEAR_Single_Sim.R)
# -----------------------------------------------------------------------------
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
# 2. Generate ONCE
# -----------------------------------------------------------------------------
site_data <- Map(draw_site, seq_len(N_SITES), SITE_N)
pooled    <- do.call(rbind, site_data)
n_h       <- sapply(site_data, nrow)
N         <- sum(n_h)
w_site    <- n_h / N

truth <- calc_empirical_truth(pooled)


# -----------------------------------------------------------------------------
# 3. CLEAR arm  (identical pipeline to CLEAR_Single_Sim.R)
# -----------------------------------------------------------------------------
clear <- CLEAR_sim(site_data, K_comp = K_COMP, n_0 = N_0, inflate = 1)
Z_0   <- clear$Z_0
site_keys <- paste0("Site_", seq_len(N_SITES))

w_global <- Reduce(`+`, Map(function(k, w) w * clear[[k]]$Weights, site_keys, w_site))
w_global <- w_global / sum(w_global)

clear_mean <- Reduce(`+`, Map(function(k, w) w * clear[[k]]$Mean, site_keys, w_site))
clear_ex2  <- Reduce(`+`, Map(function(k, w) w * (clear[[k]]$Variance + clear[[k]]$Mean^2),
                              site_keys, w_site))
clear_var  <- clear_ex2 - clear_mean^2
clear_ecc  <- Reduce(`+`, Map(function(k, w) {
  mu_h <- clear[[k]]$Mean; w * (clear[[k]]$Covariance + outer(mu_h, mu_h))
}, site_keys, w_site))
clear_cov  <- clear_ecc - outer(clear_mean, clear_mean)

clear_quant <- t(sapply(seq_len(D), function(v) {
  x <- Z_0[, v]; ord <- order(x); xo <- x[ord]; cw <- cumsum(w_global[ord])
  sapply(Q_PROBS, function(p) xo[which(cw >= p)[1]])
}))
dimnames(clear_quant) <- list(VAR_NMS, Q_LABS)


# -----------------------------------------------------------------------------
# 4. Federated synthpop arm  (the EXISTING method)
#    Synthesize AT EACH SITE, then pool synthetic sets. Each synthetic set is
#    sized to that site's real n so pooled synthetic N == real N and the site
#    mixing is preserved. No raw data leaves a site -- same privacy setting as
#    CLEAR.
# -----------------------------------------------------------------------------
syn_fed_list <- lapply(seq_len(N_SITES), function(h) {
  s <- syn(site_data[[h]], method = "cart", m = 1,
           k = n_h[h], seed = 1000 + h, print.flag = FALSE)
  s$syn
})
syn_fed <- do.call(rbind, syn_fed_list)

syn_stats <- function(df) {
  n     <- nrow(df)
  means <- colMeans(df)
  vars  <- apply(df, 2, function(x) sum((x - mean(x))^2) / n)
  quant <- t(apply(df, 2, quantile, probs = Q_PROBS)); colnames(quant) <- Q_LABS
  covm  <- cov(df) * ((n - 1) / n)
  list(Mean = means, Variance = vars, Quantiles = quant, Covariance = covm)
}
sfed <- syn_stats(syn_fed)


# =============================================================================
# (A) DESCRIPTIVE STATISTICS  --  Truth vs CLEAR vs synFED
# =============================================================================
cat("\n================= (A) DESCRIPTIVE STATISTICS =================\n")

cat("\n--- Mean ---\n")
print(round(data.frame(
  Truth   = truth$Mean,
  CLEAR   = clear_mean,
  synFED  = sfed$Mean,
  dCLEAR  = clear_mean - truth$Mean,
  dsynFED = sfed$Mean  - truth$Mean), 4))

cat("\n--- Variance ---\n")
print(round(data.frame(
  Truth   = truth$Variance,
  CLEAR   = clear_var,
  synFED  = sfed$Variance,
  dCLEAR  = clear_var     - truth$Variance,
  dsynFED = sfed$Variance - truth$Variance), 4))

cat("\n--- Quantiles: Truth ---\n");        print(round(truth$Quantiles, 4))
cat("\n--- Quantiles: CLEAR ---\n");        print(round(clear_quant, 4))
cat("\n--- Quantiles: synthpop FED ---\n"); print(round(sfed$Quantiles, 4))

cat("\n--- Covariance: Truth ---\n");        print(round(truth$Covariance, 4))
cat("\n--- Covariance: CLEAR ---\n");        print(round(clear_cov, 4))
cat("\n--- Covariance: synthpop FED ---\n"); print(round(sfed$Covariance, 4))


# =============================================================================
# (B) PER-SITE MARGINAL DENSITY RECOVERY
# =============================================================================
weighted_density_df <- function(z, w, var_name, label, n_draw = 20000) {
  idx <- sample(seq_along(z), n_draw, replace = TRUE, prob = w)
  data.frame(value = z[idx], Variable = var_name, Source = label)
}

for (h in seq_len(N_SITES)) {
  w_h   <- clear[[site_keys[h]]]$Weights
  syn_h <- syn_fed_list[[h]]
  long_h <- do.call(rbind, lapply(seq_len(D), function(v) {
    rbind(
      data.frame(value = site_data[[h]][, v], Variable = VAR_NMS[v], Source = "Raw site data"),
      weighted_density_df(Z_0[, v], w_h, VAR_NMS[v], "CLEAR weighted ref"),
      data.frame(value = syn_h[, v],          Variable = VAR_NMS[v], Source = "synthpop (site)")
    )
  }))
  ph <- ggplot(long_h, aes(value, color = Source, fill = Source)) +
    geom_density(alpha = 0.12) +
    facet_wrap(~ Variable, scales = "free", ncol = 3) +
    theme_bw(base_size = 10) +
    labs(title = paste0("Site ", h, ": density  (Raw vs CLEAR vs synthpop)"),
         x = NULL, y = "Density")
  print(ph)
}


# =============================================================================
# (C) GLOBAL MARGINAL DENSITY RECOVERY
# =============================================================================
long_global <- do.call(rbind, lapply(seq_len(D), function(v) {
  rbind(
    data.frame(value = pooled[, v],  Variable = VAR_NMS[v], Source = "Pooled raw data"),
    weighted_density_df(Z_0[, v], w_global, VAR_NMS[v], "CLEAR global ref"),
    data.frame(value = syn_fed[, v], Variable = VAR_NMS[v], Source = "synthpop FED")
  )
}))
p_global <- ggplot(long_global, aes(value, color = Source, fill = Source)) +
  geom_density(alpha = 0.12) +
  facet_wrap(~ Variable, scales = "free", ncol = 3) +
  theme_bw(base_size = 11) +
  labs(title = "Global density recovery: Truth vs CLEAR vs synthpop (FED)",
       subtitle = sprintf("Proportions = (%s)", paste(round(PROPORTIONS, 2), collapse = ", ")),
       x = NULL, y = "Density")
print(p_global)


# =============================================================================
# (B2) PER-SITE JOINT DENSITY RECOVERY (pairwise 2D heatmaps)
# -----------------------------------------------------------------------------
# Three rows: raw site data / CLEAR weight-resampled reference / synthpop (site).
# Columns = variable pairs (X1-X2, X1-X3, X2-X3).
# =============================================================================
var_pairs <- combn(VAR_NMS, 2, simplify = FALSE)

weighted_resample <- function(Z, w, n_draw = 20000) {
  idx <- sample(seq_len(nrow(Z)), n_draw, replace = TRUE, prob = w)
  Z[idx, , drop = FALSE]
}

for (h in seq_len(N_SITES)) {
  w_h    <- clear[[site_keys[h]]]$Weights
  ref_rs <- weighted_resample(Z_0, w_h)
  syn_h  <- syn_fed_list[[h]]
  
  long_pairs <- do.call(rbind, lapply(var_pairs, function(pr) {
    lab <- paste(pr[1], pr[2], sep = " vs ")
    rbind(
      data.frame(xv = site_data[[h]][, pr[1]], yv = site_data[[h]][, pr[2]],
                 Pair = lab, Source = "Raw site data"),
      data.frame(xv = ref_rs[, pr[1]],         yv = ref_rs[, pr[2]],
                 Pair = lab, Source = "CLEAR weighted ref"),
      data.frame(xv = syn_h[, pr[1]],          yv = syn_h[, pr[2]],
                 Pair = lab, Source = "synthpop (site)")
    )
  }))
  long_pairs$Source <- factor(long_pairs$Source,
                              levels = c("Raw site data", "CLEAR weighted ref", "synthpop (site)"))
  
  ph_joint <- ggplot(long_pairs, aes(xv, yv)) +
    stat_density_2d(aes(fill = after_stat(density)),
                    geom = "raster", contour = FALSE) +
    facet_grid(Source ~ Pair, scales = "free") +   # 3 rows x pairs columns
    scale_fill_viridis_c(option = "magma") +
    theme_bw(base_size = 10) +
    theme(legend.position = "right") +
    labs(title = paste0("Site ", h, ": joint density heatmap  (Raw / CLEAR / synthpop)"),
         x = NULL, y = NULL, fill = "Density")
  print(ph_joint)
}


# =============================================================================
# (C2) GLOBAL JOINT DENSITY RECOVERY (pairwise 2D heatmaps)
# -----------------------------------------------------------------------------
# Pooled raw / CLEAR global-weight resample / pooled synthpop-FED.
# =============================================================================
ref_global_rs <- weighted_resample(Z_0, w_global)

long_pairs_global <- do.call(rbind, lapply(var_pairs, function(pr) {
  lab <- paste(pr[1], pr[2], sep = " vs ")
  rbind(
    data.frame(xv = pooled[, pr[1]],        yv = pooled[, pr[2]],
               Pair = lab, Source = "Pooled raw data"),
    data.frame(xv = ref_global_rs[, pr[1]], yv = ref_global_rs[, pr[2]],
               Pair = lab, Source = "CLEAR global ref"),
    data.frame(xv = syn_fed[, pr[1]],       yv = syn_fed[, pr[2]],
               Pair = lab, Source = "synthpop FED")
  )
}))
long_pairs_global$Source <- factor(long_pairs_global$Source,
                                   levels = c("Pooled raw data", "CLEAR global ref", "synthpop FED"))

p_joint_global <- ggplot(long_pairs_global, aes(xv, yv)) +
  stat_density_2d(aes(fill = after_stat(density)),
                  geom = "raster", contour = FALSE) +
  facet_grid(Source ~ Pair, scales = "free") +     # 3 rows x pairs columns
  scale_fill_viridis_c(option = "magma") +
  theme_bw(base_size = 11) +
  theme(legend.position = "right") +
  labs(title = "Global joint density heatmap: pooled raw vs CLEAR vs synthpop (FED)",
       subtitle = sprintf("Proportions = (%s)", paste(round(PROPORTIONS, 2), collapse = ", ")),
       x = NULL, y = NULL, fill = "Density")
print(p_joint_global)


# =============================================================================
# (D) LASSO:  X1 ~ X2 + X3   coefficients, CLEAR & synFED vs pooled-raw truth
# =============================================================================
outcome <- "X1"; predictors <- c("X2", "X3")

fit_lasso <- function(X, y, w = NULL) {
  set.seed(1)
  cv <- if (is.null(w)) cv.glmnet(X, y, alpha = 1, nfolds = 10)
  else            cv.glmnet(X, y, alpha = 1, weights = w, nfolds = 10)
  as.vector(coef(cv, s = "lambda.min"))
}

coef_truth <- fit_lasso(as.matrix(pooled[, predictors]),  pooled[[outcome]])
coef_clear <- fit_lasso(as.matrix(Z_0[, predictors]),     Z_0[[outcome]], w = w_global)
coef_syn   <- fit_lasso(as.matrix(syn_fed[, predictors]), syn_fed[[outcome]])

lasso_tbl <- data.frame(
  Term       = c("(Intercept)", predictors),
  Truth_pool = round(coef_truth, 4),
  CLEAR      = round(coef_clear, 4),
  synFED     = round(coef_syn,   4),
  dCLEAR     = round(coef_clear - coef_truth, 4),
  dsynFED    = round(coef_syn   - coef_truth, 4)
)
cat("\n================= (D) LASSO COEFFICIENTS =================\n")
cat(sprintf("Model: %s ~ %s  (lambda.min); Truth = pooled raw\n",
            outcome, paste(predictors, collapse = " + ")))
print(lasso_tbl, row.names = FALSE)


# =============================================================================
# (E) ERROR SUMMARY  --  CLEAR vs synFED, distance to TRUTH (lower = better)
# =============================================================================
rmse <- function(a, b) sqrt(mean((a - b)^2))
frob <- function(A, B) sqrt(sum((A - B)^2))

err <- data.frame(
  Metric = c("Mean_RMSE", "Var_RMSE", "Quantile_RMSE", "Cov_Frobenius", "Lasso_RMSE"),
  CLEAR = c(
    rmse(clear_mean, truth$Mean),
    rmse(clear_var,  truth$Variance),
    rmse(as.vector(clear_quant), as.vector(truth$Quantiles)),
    frob(clear_cov,  truth$Covariance),
    rmse(coef_clear, coef_truth)
  ),
  synFED = c(
    rmse(sfed$Mean,      truth$Mean),
    rmse(sfed$Variance,  truth$Variance),
    rmse(as.vector(sfed$Quantiles), as.vector(truth$Quantiles)),
    frob(sfed$Covariance, truth$Covariance),
    rmse(coef_syn,       coef_truth)
  )
)
err[ , -1] <- round(err[ , -1], 4)
err$Winner <- ifelse(err$CLEAR < err$synFED, "CLEAR", "synFED")

cat("\n================= (E) ERROR SUMMARY (lower = better) =================\n")
print(err, row.names = FALSE)


