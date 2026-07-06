# =============================================================================
# CLEAR parallel simulation study -- performance across REFERENCE SAMPLE SIZE
# -----------------------------------------------------------------------------
# K is FIXED at 3 for BOTH the DGP and CLEAR. The DGP uses the 3-mixture
# proportions c(1/3, 1/3, 1/3, 0) throughout. Instead of varying DGP
# complexity, we vary the reference-distribution sample size n_0 to see how
# it affects the recovery of pooled statistics.
#
# Settings (only n_0 varies):
#     n_0 = 1000, 2500, 5000, 10000, 20000, 40000
#
# 500 iterations per setting, parallel over available cores.
# Each iteration: generate data -> run CLEAR -> compute global estimates ->
# compare to pooled empirical truth for: Mean, Variance, Covariance,
# Quantiles, and Lasso coefficients (X1 ~ X2 + X3).
#
# Per-element summary metric per statistic = SIGNED BIAS (estimate - truth).
#
# Output: per-statistic figures (x = n_0) plus a CSV summary.
# Requires CLEAR.R.
# =============================================================================

rm(list = ls())
library(parallel)
library(MASS)
library(mclust)
library(glmnet)
library(ggplot2)
source("CLEAR.R")


# -----------------------------------------------------------------------------
# 0. Fixed configuration
# -----------------------------------------------------------------------------
N_ITER  <- 500
N_SITES <- 10
D       <- 3
VAR_NMS <- c("X1", "X2", "X3")
K_COMP  <- 3                              # FIXED for CLEAR
SITE_N  <- round(seq(200, 800, length.out = N_SITES))
Q_PROBS <- c(0.05, 0.25, 0.50, 0.75, 0.95)
EPS     <- 1e-8

# DGP proportions are FIXED at the 3-mixture (K = 3) throughout.
DGP_PROPS <- c(1/3, 1/3, 1/3, 0)

# The settings we now sweep over: reference sample size n_0.
N0_GRID <- c(1000, 2500, 5000, 10000, 20000, 40000)
SETTINGS <- setNames(as.list(N0_GRID), paste0("n0=", N0_GRID))

# Fixed DGP component bank (same as before).
COMP_CENTERS <- list(c( 0,  0,  0), c( 4,  4, -3),
                     c(-4,  3,  3), c( 3, -4,  4))
comp_cov <- function(s, rho) { M <- matrix(rho, D, D); diag(M) <- 1; (s^2) * M }
COMP_COVS <- list(comp_cov(1.0,  0.3), comp_cov(1.4, -0.2),
                  comp_cov(1.1,  0.1), comp_cov(1.3, -0.3))
site_offset <- function(h) c((h - 3) * 1.5, (h - 3) * -1.0, (h - 3) * 0.8)

draw_site <- function(h, n, props) {
  off  <- site_offset(h)
  comp <- sample(seq_len(4), n, replace = TRUE, prob = props)
  X <- matrix(0, n, D)
  for (k in unique(comp)) {
    idx <- which(comp == k)
    X[idx, ] <- mvrnorm(length(idx), mu = COMP_CENTERS[[k]] + off, Sigma = COMP_COVS[[k]])
  }
  df <- as.data.frame(X); colnames(df) <- VAR_NMS; df
}

rel_rmse <- function(est, truth) {
  sqrt(mean((est - truth)^2)) / (sqrt(mean(truth^2)) + EPS)
}


# -----------------------------------------------------------------------------
# 1. One iteration for a given n_0
# -----------------------------------------------------------------------------
run_one <- function(iter_seed, n0) {
  set.seed(iter_seed)
  props <- DGP_PROPS / sum(DGP_PROPS)

  site_data <- Map(function(h, n) draw_site(h, n, props), seq_len(N_SITES), SITE_N)
  pooled    <- do.call(rbind, site_data)
  n_h       <- sapply(site_data, nrow)
  w_site    <- n_h / sum(n_h)

  truth <- calc_empirical_truth(pooled)

  clear <- tryCatch(CLEAR_sim(site_data, K_comp = K_COMP, n_0 = n0, inflate = 1),
                    error = function(e) NULL)
  if (is.null(clear)) return(NULL)

  Z_0       <- clear$Z_0
  site_keys <- paste0("Site_", seq_len(N_SITES))
  w_global  <- Reduce(`+`, Map(function(k, w) w * clear[[k]]$Weights, site_keys, w_site))
  w_global  <- w_global / sum(w_global)

  # --- Global CLEAR estimates ---
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

  # --- Lasso: X1 ~ X2 + X3, pooled vs. weighted reference ---
  lasso_bias <- tryCatch({
    Xp <- as.matrix(pooled[, c("X2", "X3")]); yp <- pooled[["X1"]]
    Xr <- as.matrix(Z_0[, c("X2", "X3")]);    yr <- Z_0[["X1"]]
    set.seed(1); cp <- as.vector(coef(cv.glmnet(Xp, yp, alpha = 1, nfolds = 5), s = "lambda.min"))
    set.seed(1); cr <- as.vector(coef(cv.glmnet(Xr, yr, alpha = 1, weights = w_global, nfolds = 5),
                                       s = "lambda.min"))
    cr - cp   # signed bias per coefficient; truth = pooled-raw lasso coefs
  }, error = function(e) rep(NA_real_, 3))

  # --- Per-element SIGNED BIAS (estimate - truth) ---
  cov_lt   <- lower.tri(truth$Covariance, diag = TRUE)
  cov_idx  <- which(cov_lt, arr.ind = TRUE)
  cov_nms  <- apply(cov_idx, 1, function(ij) paste(VAR_NMS[ij[1]], VAR_NMS[ij[2]], sep = "-"))
  q_nms    <- as.vector(outer(VAR_NMS, paste0("q", Q_PROBS * 100), paste, sep = "_"))

  out <- c(
    setNames(est_mean - truth$Mean,                          paste0("Mean__",  VAR_NMS)),
    setNames(est_var  - truth$Variance,                      paste0("Var__",   VAR_NMS)),
    setNames(est_cov[cov_lt] - truth$Covariance[cov_lt],     paste0("Cov__",   cov_nms)),
    setNames(as.vector(t(est_quant)) - as.vector(t(truth$Quantiles)),
             paste0("Quantile__", q_nms)),
    setNames(lasso_bias,                                     paste0("Lasso__",
             c("(Intercept)", "X2", "X3")))
  )
  out
}


# -----------------------------------------------------------------------------
# 2. Run all settings in parallel
# -----------------------------------------------------------------------------
n_cores <- parallelly::availableCores()
message(sprintf("Running %d n_0 settings x %d iters on %d cores...",
                length(SETTINGS), N_ITER, n_cores))

cl <- makeCluster(n_cores)
clusterExport(cl, c("N_SITES", "D", "VAR_NMS", "K_COMP", "SITE_N", "Q_PROBS",
                    "EPS", "DGP_PROPS", "COMP_CENTERS", "COMP_COVS", "comp_cov",
                    "site_offset", "draw_site", "rel_rmse", "run_one",
                    "CLEAR_sim", "calc_empirical_truth", "site_fit_gmm",
                    "center_make_reference", "estimate_ipw_weights", "weighted_eda_stats"))
clusterEvalQ(cl, { library(MASS); library(mclust); library(glmnet) })

all_rows <- list()
for (s in names(SETTINGS)) {
  message("  setting: ", s)
  seeds <- seq_len(N_ITER) + 1000L * match(s, names(SETTINGS))
  res   <- parLapply(cl, seeds, run_one, n0 = SETTINGS[[s]])
  res   <- res[!sapply(res, is.null)]
  nm    <- names(res[[1]])
  M     <- do.call(rbind, lapply(res, function(r) r[nm]))
  colnames(M) <- nm
  df    <- as.data.frame(M, check.names = FALSE)
  df$Setting <- s
  all_rows[[s]] <- df
  message(sprintf("    %d / %d iterations succeeded.", nrow(df), N_ITER))
}
stopCluster(cl)

results <- do.call(rbind, all_rows)


# -----------------------------------------------------------------------------
# 3. Reshape to long + summary CSV
# -----------------------------------------------------------------------------
quantity_cols <- setdiff(colnames(results), "Setting")

long <- do.call(rbind, lapply(quantity_cols, function(q) {
  parts <- strsplit(q, "__", fixed = TRUE)[[1]]
  data.frame(Setting   = results$Setting,
             Statistic = parts[1],
             Quantity  = parts[2],
             Bias      = results[[q]],
             stringsAsFactors = FALSE)
}))

stat_levels <- c("Mean", "Var", "Cov", "Quantile", "Lasso")
long$Setting   <- factor(long$Setting, levels = names(SETTINGS))
long$Statistic <- factor(long$Statistic, levels = stat_levels)
long$Panel <- factor(paste0(long$Statistic, ": ", long$Quantity),
                      levels = unique(paste0(long$Statistic, ": ", long$Quantity)))
long <- long[is.finite(long$Bias), ]

summary_tbl <- aggregate(Bias ~ Setting + Statistic + Quantity, data = long,
                         FUN = function(x) c(bias_mean = mean(x), bias_median = median(x),
                                             sd = sd(x)))
summary_tbl <- do.call(data.frame, summary_tbl)
print(summary_tbl, digits = 4, row.names = FALSE)
write.csv(summary_tbl, "CLEAR_refsize_bias_summary.csv", row.names = FALSE)


# -----------------------------------------------------------------------------
# 4. Per-statistic figures, x = n_0, shared y-axis across panels
# -----------------------------------------------------------------------------
stat_titles <- c(
  Mean     = "CLEAR bias vs reference size: Mean",
  Var      = "CLEAR bias vs reference size: Variance",
  Cov      = "CLEAR bias vs reference size: Covariance",
  Quantile = "CLEAR bias vs reference size: Quantiles",
  Lasso    = "CLEAR bias vs reference size: Lasso coefficients"
)

make_stat_plot <- function(st) {
  d <- long[long$Statistic == st, ]
  d$Quantity <- factor(d$Quantity, levels = unique(d$Quantity))
  n_q <- nlevels(d$Quantity)

  ggplot(d, aes(x = Setting, y = Bias, fill = Setting)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.3) +
    geom_boxplot(outlier.size = 0.3, alpha = 0.85, linewidth = 0.3) +
    facet_wrap(~ Quantity, nrow = 1) +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 7),
          strip.text  = element_text(size = 8),
          legend.position = "none") +
    labs(title = stat_titles[[st]],
         x = "Reference sample size (n_0)", y = "Bias (estimate - truth)")
}

for (st in levels(long$Statistic)) {
  if (!any(long$Statistic == st)) next
  p_st <- make_stat_plot(st)
  n_q  <- length(unique(long$Quantity[long$Statistic == st]))
  ggsave(sprintf("CLEAR_refsize_bias_%s.png", st), p_st,
         width = max(6, n_q * 2.2), height = 4.2, dpi = 150, limitsize = FALSE)
  print(p_st)
}

cat("\nSaved: CLEAR_refsize_bias_summary.csv and per-statistic plots ",
    "(CLEAR_refsize_bias_Mean.png, _Var.png, _Cov.png, _Quantile.png, _Lasso.png)\n", sep = "")