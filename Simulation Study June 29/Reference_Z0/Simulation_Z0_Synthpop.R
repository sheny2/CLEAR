# =============================================================================
# CLEAR parallel simulation study -- performance across REFERENCE SAMPLE SIZE
#                                    + federated synthpop as a comparator
# -----------------------------------------------------------------------------
# K is FIXED at 3 for BOTH the DGP and CLEAR. The DGP uses the 3-mixture
# proportions c(1/3, 1/3, 1/3, 0) throughout. We vary the reference-distribution
# sample size n_0 to see how it affects recovery of pooled statistics, and we
# add an EXISTING synthetic-data method -- federated synthpop (CART) -- as an
# extra category on every plot.
#
# Settings on the x-axis of every figure:
#     n0=1000, n0=2500, n0=5000, n0=10000, n0=20000, n0=40000, synthpop
#
# synthpop is federated (each site synthesizes locally, synthetic sets pooled;
# no raw data leaves a site -- same privacy setting as CLEAR) and does NOT
# depend on n_0, so it is computed ONCE per iteration and shown as a single
# extra category appended after the n_0 grid.
#
# 500 iterations per setting, parallel over available cores.
# Each iteration: generate data -> run CLEAR (per n_0) and federated synthpop ->
# compute global estimates -> compare to pooled empirical truth for: Mean,
# Variance, Covariance, Quantiles, and Lasso coefficients (X1 ~ X2 + X3).
#
# Per-element summary metric per statistic = SIGNED BIAS (estimate - truth).
#
# Output: per-statistic figures (x = setting) plus a CSV summary.
# Requires CLEAR.R and the synthpop package.
# =============================================================================

rm(list = ls())
library(parallel)
library(MASS)
library(mclust)
library(glmnet)
library(ggplot2)
library(synthpop)

source("CLEAR.R")


# -----------------------------------------------------------------------------
# 0. Fixed configuration
# -----------------------------------------------------------------------------
N_ITER  <- 5
N_SITES <- 10
D       <- 3
VAR_NMS <- c("X1", "X2", "X3")
K_COMP  <- 3                              # FIXED for CLEAR
SITE_N  <- round(seq(200, 800, length.out = N_SITES))
Q_PROBS <- c(0.05, 0.25, 0.50, 0.75, 0.95)
EPS     <- 1e-8

# DGP proportions are FIXED at the 3-mixture (K = 3) throughout.
DGP_PROPS <- c(1/3, 1/3, 1/3, 0)

# The n_0 settings we sweep over, plus a "synthpop" category (n_0 not used).
N0_GRID   <- c(1000, 2500, 5000, 10000, 20000, 40000)
CLEAR_SETTINGS <- setNames(as.list(N0_GRID), paste0("n0=", N0_GRID))
SYN_LABEL <- "synthpop"
# Full x-axis ordering: all CLEAR n_0 settings first, then synthpop last.
SETTING_LEVELS <- c(names(CLEAR_SETTINGS), SYN_LABEL)

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
# 0b. Bias-vector helper: given global estimates, return the same named
#     signed-bias vector used by both CLEAR and synthpop arms.
# -----------------------------------------------------------------------------
make_bias_vec <- function(est_mean, est_var, est_cov, est_quant, lasso_bias, truth) {
  cov_lt  <- lower.tri(truth$Covariance, diag = TRUE)
  cov_idx <- which(cov_lt, arr.ind = TRUE)
  cov_nms <- apply(cov_idx, 1, function(ij) paste(VAR_NMS[ij[1]], VAR_NMS[ij[2]], sep = "-"))
  q_nms   <- as.vector(outer(VAR_NMS, paste0("q", Q_PROBS * 100), paste, sep = "_"))
  
  c(
    setNames(est_mean - truth$Mean,                      paste0("Mean__",  VAR_NMS)),
    setNames(est_var  - truth$Variance,                  paste0("Var__",   VAR_NMS)),
    setNames(est_cov[cov_lt] - truth$Covariance[cov_lt], paste0("Cov__",   cov_nms)),
    setNames(as.vector(t(est_quant)) - as.vector(t(truth$Quantiles)),
             paste0("Quantile__", q_nms)),
    setNames(lasso_bias, paste0("Lasso__", c("(Intercept)", "X2", "X3")))
  )
}


# -----------------------------------------------------------------------------
# 1. One CLEAR iteration for a given n_0
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
  
  make_bias_vec(est_mean, est_var, est_cov, est_quant, lasso_bias, truth)
}


# -----------------------------------------------------------------------------
# 1b. One federated-synthpop iteration (no n_0 dependence).
#     Each site synthesizes locally via CART; synthetic sets are pooled (sized
#     to each site's real n). Plain unweighted stats on the pooled synthetic
#     data, scored against the same pooled empirical truth.
# -----------------------------------------------------------------------------
run_one_syn <- function(iter_seed) {
  set.seed(iter_seed)
  props <- DGP_PROPS / sum(DGP_PROPS)
  
  site_data <- Map(function(h, n) draw_site(h, n, props), seq_len(N_SITES), SITE_N)
  pooled    <- do.call(rbind, site_data)
  n_h       <- sapply(site_data, nrow)
  
  truth <- calc_empirical_truth(pooled)
  
  syn_pooled <- tryCatch({
    syn_list <- lapply(seq_len(N_SITES), function(h) {
      s <- syn(site_data[[h]], method = "cart", m = 1, k = n_h[h],
               seed = iter_seed * 100L + h, print.flag = FALSE)
      s$syn
    })
    do.call(rbind, syn_list)
  }, error = function(e) NULL)
  if (is.null(syn_pooled)) return(NULL)
  
  n_s <- nrow(syn_pooled)
  
  # --- Global synthpop estimates (plain empirical, population moments) ---
  est_mean <- colMeans(syn_pooled)
  est_var  <- apply(syn_pooled, 2, function(x) sum((x - mean(x))^2) / n_s)
  est_cov  <- cov(syn_pooled) * ((n_s - 1) / n_s)
  est_quant <- t(apply(syn_pooled, 2, quantile, probs = Q_PROBS))
  
  # --- Lasso: X1 ~ X2 + X3, pooled-raw truth vs. synthetic-data fit ---
  lasso_bias <- tryCatch({
    Xp <- as.matrix(pooled[, c("X2", "X3")]);     yp <- pooled[["X1"]]
    Xs <- as.matrix(syn_pooled[, c("X2", "X3")]); ys <- syn_pooled[["X1"]]
    set.seed(1); cp <- as.vector(coef(cv.glmnet(Xp, yp, alpha = 1, nfolds = 5), s = "lambda.min"))
    set.seed(1); cs <- as.vector(coef(cv.glmnet(Xs, ys, alpha = 1, nfolds = 5), s = "lambda.min"))
    cs - cp
  }, error = function(e) rep(NA_real_, 3))
  
  make_bias_vec(est_mean, est_var, est_cov, est_quant, lasso_bias, truth)
}


# -----------------------------------------------------------------------------
# 2. Run all settings in parallel
# -----------------------------------------------------------------------------
n_cores <- parallelly::availableCores()
message(sprintf("Running %d CLEAR n_0 settings + synthpop, x %d iters on %d cores...",
                length(CLEAR_SETTINGS), N_ITER, n_cores))

cl <- makeCluster(n_cores)
clusterExport(cl, c("N_SITES", "D", "VAR_NMS", "K_COMP", "SITE_N", "Q_PROBS",
                    "EPS", "DGP_PROPS", "COMP_CENTERS", "COMP_COVS", "comp_cov",
                    "site_offset", "draw_site", "rel_rmse", "make_bias_vec",
                    "run_one", "run_one_syn",
                    "CLEAR_sim", "calc_empirical_truth", "site_fit_gmm",
                    "center_make_reference", "estimate_ipw_weights", "weighted_eda_stats"))
clusterEvalQ(cl, { library(MASS); library(mclust); library(glmnet); library(synthpop) })

all_rows <- list()

# --- CLEAR settings (vary n_0) ---
for (s in names(CLEAR_SETTINGS)) {
  message("  CLEAR setting: ", s)
  seeds <- seq_len(N_ITER) + 1000L * match(s, names(CLEAR_SETTINGS))
  res   <- parLapply(cl, seeds, run_one, n0 = CLEAR_SETTINGS[[s]])
  res   <- res[!sapply(res, is.null)]
  nm    <- names(res[[1]])
  M     <- do.call(rbind, lapply(res, function(r) r[nm]))
  colnames(M) <- nm
  df    <- as.data.frame(M, check.names = FALSE)
  df$Setting <- s
  all_rows[[s]] <- df
  message(sprintf("    %d / %d iterations succeeded.", nrow(df), N_ITER))
}

# --- synthpop setting (single extra category; reuse the same seed offsets so
#     it sees the same simulated datasets as one CLEAR block for fairness) ---
message("  setting: ", SYN_LABEL)
seeds_syn <- seq_len(N_ITER) + 1000L * (length(CLEAR_SETTINGS) + 1L)
res_syn <- parLapply(cl, seeds_syn, run_one_syn)
res_syn <- res_syn[!sapply(res_syn, is.null)]
nm_syn  <- names(res_syn[[1]])
M_syn   <- do.call(rbind, lapply(res_syn, function(r) r[nm_syn]))
colnames(M_syn) <- nm_syn
df_syn  <- as.data.frame(M_syn, check.names = FALSE)
df_syn$Setting <- SYN_LABEL
all_rows[[SYN_LABEL]] <- df_syn
message(sprintf("    %d / %d iterations succeeded.", nrow(df_syn), N_ITER))

stopCluster(cl)

# Align columns across all blocks (CLEAR and synthpop share identical names).
common_cols <- Reduce(intersect, lapply(all_rows, colnames))
all_rows <- lapply(all_rows, function(df) df[, common_cols, drop = FALSE])
results  <- do.call(rbind, all_rows)


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
long$Setting   <- factor(long$Setting, levels = SETTING_LEVELS)   # synthpop appears last
long$Statistic <- factor(long$Statistic, levels = stat_levels)
long$Panel <- factor(paste0(long$Statistic, ": ", long$Quantity),
                     levels = unique(paste0(long$Statistic, ": ", long$Quantity)))
long <- long[is.finite(long$Bias), ]

summary_tbl <- aggregate(Bias ~ Setting + Statistic + Quantity, data = long,
                         FUN = function(x) c(bias_mean = mean(x), bias_median = median(x),
                                             sd = sd(x)))
summary_tbl <- do.call(data.frame, summary_tbl)
print(summary_tbl, digits = 4, row.names = FALSE)
# write.csv(summary_tbl, "CLEAR_refsize_bias_summary.csv", row.names = FALSE)


# -----------------------------------------------------------------------------
# 4. Per-statistic figures, x = setting (n_0 grid + synthpop), shared y-axis
# -----------------------------------------------------------------------------
stat_titles <- c(
  Mean     = "CLEAR (by n_0) vs synthpop: Mean bias",
  Var      = "CLEAR (by n_0) vs synthpop: Variance bias",
  Cov      = "CLEAR (by n_0) vs synthpop: Covariance bias",
  Quantile = "CLEAR (by n_0) vs synthpop: Quantile bias",
  Lasso    = "CLEAR (by n_0) vs synthpop: Lasso coefficient bias"
)

# Color synthpop distinctly from the n_0 gradient so it reads as "the comparator".
setting_fills <- c(
  setNames(scales::seq_gradient_pal("#c6dbef", "#08306b")(
    seq(0, 1, length.out = length(names(CLEAR_SETTINGS)))),
    names(CLEAR_SETTINGS)),
  setNames("#e6550d", SYN_LABEL)
)

make_stat_plot <- function(st) {
  d <- long[long$Statistic == st, ]
  d$Quantity <- factor(d$Quantity, levels = unique(d$Quantity))
  
  ggplot(d, aes(x = Setting, y = Bias, fill = Setting)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.3) +
    geom_boxplot(outlier.size = 0.3, alpha = 0.85, linewidth = 0.3) +
    facet_wrap(~ Quantity, nrow = 1) +
    scale_fill_manual(values = setting_fills) +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 7),
          strip.text  = element_text(size = 8),
          legend.position = "none") +
    labs(title = stat_titles[[st]],
         x = "Setting (CLEAR reference size n_0, then synthpop)",
         y = "Bias (estimate - truth)")
}

for (st in levels(long$Statistic)) {
  if (!any(long$Statistic == st)) next
  p_st <- make_stat_plot(st)
  n_q  <- length(unique(long$Quantity[long$Statistic == st]))
  ggsave(sprintf("CLEAR_refsize_synthpop_bias_%s.png", st), p_st,
         width = max(6, n_q * 2.4), height = 4.2, dpi = 150, limitsize = FALSE)
  print(p_st)
}

cat("\nSaved: CLEAR_refsize_bias_summary.csv and per-statistic plots ",
    "(CLEAR_refsize_bias_Mean.png, _Var.png, _Cov.png, _Quantile.png, _Lasso.png).\n",
    "Each plot's x-axis is: n0=1000, 2500, 5000, 10000, 20000, 40000, then synthpop.\n", sep = "")
