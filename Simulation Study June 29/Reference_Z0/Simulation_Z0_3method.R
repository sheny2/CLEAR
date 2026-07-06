# =============================================================================
# Simulation study: CLEAR vs synthpop vs CLEAR_synthpop
#                    across THREE reference-sample-size settings
# -----------------------------------------------------------------------------
# THREE METHODS
#   CLEAR           : sites share GMM params -> pooled GMM reference Z_0
#                     -> per-site density-ratio IPW -> weighted EDA.
#   synthpop        : federated CART synthesis per site -> pooled synthetic set
#                     is the estimate; plain (unweighted) stats on it.
#   CLEAR_synthpop  : sites synthesize via CART -> pooled synthetic set is the
#                     reference Z_0 -> ORIGINAL CLEAR IPW + weighted EDA on it.
#
# THREE REFERENCE-SIZE SETTINGS  (relative to total pooled site N)
#   half   = round(0.5 * N_total)
#   equal  =        N_total
#   double =        2   * N_total
#   where N_total = sum of all site sample sizes.
#   * For CLEAR and CLEAR_synthpop this is n_0 (reference size).
#   * For synthpop this is the total number of synthetic rows generated
#     (allocated across sites in proportion to their real n), so all three
#     methods genuinely respond to the size setting.
#
# => 3 methods x 3 sizes = 9 settings on the x-axis of every figure.
#
# Each iteration: generate multi-site data -> run all 9 method/size combos ->
# global estimates -> SIGNED BIAS (estimate - truth) for Mean, Variance,
# Covariance, Quantiles, Lasso (X1 ~ X2 + X3).
#
# Output: per-statistic figures (x = setting, grouped/colored by method) plus a
# CSV summary, in the same format as the existing reference-size study.
#
# Requires CLEAR.R, CLEAR_Synthpop.R, and the synthpop package.
# =============================================================================

rm(list = ls())
library(parallel)
library(MASS)
library(mclust)
library(glmnet)
library(ggplot2)
library(synthpop)

source("CLEAR.R")            # site_fit_gmm, center_make_reference,
# estimate_ipw_weights, weighted_eda_stats,
# CLEAR_sim, calc_empirical_truth
source("CLEAR_Synthpop.R")   # site_synthesize, center_make_reference_syn,
# CLEAR_synthpop_sim, combine_global_stats
# (CLEAR_Synthpop.R itself sources CLEAR.R; that is
#  harmless -- functions are simply redefined.)


# -----------------------------------------------------------------------------
# 0. Fixed configuration
# -----------------------------------------------------------------------------
N_ITER  <- 500
N_SITES <- 10
D       <- 3
VAR_NMS <- c("X1", "X2", "X3")
K_COMP  <- 3                              # FIXED for CLEAR / CLEAR_synthpop GMM
SITE_N  <- round(seq(200, 800, length.out = N_SITES))
Q_PROBS <- c(0.05, 0.25, 0.50, 0.75, 0.95)
EPS     <- 1e-8

N_TOTAL <- sum(SITE_N)                    # total pooled site sample size (fixed)

# DGP proportions are FIXED at the 3-mixture (K = 3) throughout.
DGP_PROPS <- c(1/3, 1/3, 1/3, 0)

# Three reference-size settings relative to total pooled N.
SIZE_FACTORS <- c(half = 0.5, equal = 1, double = 2)
SIZE_N0      <- round(SIZE_FACTORS * N_TOTAL)   # named: half / equal / double
METHODS      <- c("CLEAR", "synthpop", "CLEAR_synthpop")

# Full setting grid: one row per method x size. Label = "Method (size)".
SETTING_GRID <- expand.grid(Method = METHODS, Size = names(SIZE_FACTORS),
                            stringsAsFactors = FALSE)
SETTING_GRID$n0    <- SIZE_N0[SETTING_GRID$Size]
SETTING_GRID$Label <- paste0(SETTING_GRID$Method, " (", SETTING_GRID$Size, ")")

# x-axis ordering: group by size (half, equal, double), method within each.
SIZE_ORDER   <- names(SIZE_FACTORS)
SETTING_GRID <- SETTING_GRID[order(match(SETTING_GRID$Size, SIZE_ORDER),
                                   match(SETTING_GRID$Method, METHODS)), ]
SETTING_LEVELS <- SETTING_GRID$Label

# Fixed DGP component bank.
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


# -----------------------------------------------------------------------------
# 0b. Shared bias-vector helper (identical to the existing study).
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

# Global estimates from a CLEAR-style per-site weighted-stats object.
# (Works for both CLEAR_sim and CLEAR_synthpop_sim output.)
global_from_clear <- function(clear_out, n_h) {
  site_keys <- paste0("Site_", seq_along(n_h))
  w_site    <- n_h / sum(n_h)
  Z_0       <- clear_out$Z_0
  
  w_global <- Reduce(`+`, Map(function(k, w) w * clear_out[[k]]$Weights, site_keys, w_site))
  w_global <- w_global / sum(w_global)
  
  est_mean <- Reduce(`+`, Map(function(k, w) w * clear_out[[k]]$Mean, site_keys, w_site))
  ex2 <- Reduce(`+`, Map(function(k, w) w * (clear_out[[k]]$Variance + clear_out[[k]]$Mean^2),
                         site_keys, w_site))
  est_var <- ex2 - est_mean^2
  ecc <- Reduce(`+`, Map(function(k, w) {
    mu_h <- clear_out[[k]]$Mean; w * (clear_out[[k]]$Covariance + outer(mu_h, mu_h))
  }, site_keys, w_site))
  est_cov <- ecc - outer(est_mean, est_mean)
  
  est_quant <- t(sapply(seq_len(D), function(v) {
    x <- Z_0[, v]; ord <- order(x); xo <- x[ord]; cw <- cumsum(w_global[ord])
    sapply(Q_PROBS, function(p) xo[which(cw >= p)[1]])
  }))
  
  list(mean = est_mean, var = est_var, cov = est_cov, quant = est_quant,
       Z_0 = Z_0, w_global = w_global)
}


# -----------------------------------------------------------------------------
# 1. ONE ITERATION: run ALL 9 method/size combos on the SAME simulated data.
#    Returns a named list of bias vectors keyed by setting Label.
# -----------------------------------------------------------------------------
run_one_all <- function(iter_seed) {
  set.seed(iter_seed)
  props <- DGP_PROPS / sum(DGP_PROPS)
  
  site_data <- Map(function(h, n) draw_site(h, n, props), seq_len(N_SITES), SITE_N)
  pooled    <- do.call(rbind, site_data)
  n_h       <- sapply(site_data, nrow)
  
  truth <- calc_empirical_truth(pooled)
  
  # Pooled-raw Lasso coefs = truth for the Lasso comparison (computed once).
  cp <- tryCatch({
    Xp <- as.matrix(pooled[, c("X2", "X3")]); yp <- pooled[["X1"]]
    set.seed(1); as.vector(coef(cv.glmnet(Xp, yp, alpha = 1, nfolds = 5), s = "lambda.min"))
  }, error = function(e) rep(NA_real_, 3))
  
  out <- list()
  
  for (i in seq_len(nrow(SETTING_GRID))) {
    method <- SETTING_GRID$Method[i]
    n0     <- SETTING_GRID$n0[i]
    label  <- SETTING_GRID$Label[i]
    
    biasv <- tryCatch({
      
      if (method == "CLEAR") {
        clear <- CLEAR_sim(site_data, K_comp = K_COMP, n_0 = n0, inflate = 1)
        g <- global_from_clear(clear, n_h)
        lasso_bias <- tryCatch({
          Xr <- as.matrix(g$Z_0[, c("X2", "X3")]); yr <- g$Z_0[["X1"]]
          set.seed(1); cr <- as.vector(coef(cv.glmnet(Xr, yr, alpha = 1,
                                                      weights = g$w_global, nfolds = 5), s = "lambda.min"))
          cr - cp
        }, error = function(e) rep(NA_real_, 3))
        make_bias_vec(g$mean, g$var, g$cov, g$quant, lasso_bias, truth)
        
      } else if (method == "CLEAR_synthpop") {
        hyb <- CLEAR_synthpop_sim(site_data, n_0 = n0,
                                  method = "cart", base_seed = iter_seed * 10L)
        g <- global_from_clear(hyb, n_h)
        lasso_bias <- tryCatch({
          Xr <- as.matrix(g$Z_0[, c("X2", "X3")]); yr <- g$Z_0[["X1"]]
          set.seed(1); cr <- as.vector(coef(cv.glmnet(Xr, yr, alpha = 1,
                                                      weights = g$w_global, nfolds = 5), s = "lambda.min"))
          cr - cp
        }, error = function(e) rep(NA_real_, 3))
        make_bias_vec(g$mean, g$var, g$cov, g$quant, lasso_bias, truth)
        
      } else {  # ---- plain federated synthpop ----
        # Generate n0 synthetic rows total, allocated across sites in
        # proportion to real n so synthpop responds to the size setting.
        prob_h <- n_h / sum(n_h)
        alloc  <- as.vector(rmultinom(1, size = n0, prob = prob_h))
        syn_list <- lapply(seq_len(N_SITES), function(h) {
          k_h <- max(alloc[h], 1L)
          s <- syn(site_data[[h]], method = "cart", m = 1, k = k_h,
                   seed = iter_seed * 100L + h, print.flag = FALSE)
          s$syn
        })
        syn_pooled <- do.call(rbind, syn_list)
        n_s <- nrow(syn_pooled)
        
        est_mean  <- colMeans(syn_pooled)
        est_var   <- apply(syn_pooled, 2, function(x) sum((x - mean(x))^2) / n_s)
        est_cov   <- cov(syn_pooled) * ((n_s - 1) / n_s)
        est_quant <- t(apply(syn_pooled, 2, quantile, probs = Q_PROBS))
        
        lasso_bias <- tryCatch({
          Xs <- as.matrix(syn_pooled[, c("X2", "X3")]); ys <- syn_pooled[["X1"]]
          set.seed(1); cs <- as.vector(coef(cv.glmnet(Xs, ys, alpha = 1, nfolds = 5),
                                            s = "lambda.min"))
          cs - cp
        }, error = function(e) rep(NA_real_, 3))
        make_bias_vec(est_mean, est_var, est_cov, est_quant, lasso_bias, truth)
      }
      
    }, error = function(e) NULL)
    
    out[[label]] <- biasv
  }
  
  out
}


# -----------------------------------------------------------------------------
# 2. Run in parallel over iterations
# -----------------------------------------------------------------------------
n_cores <- parallelly::availableCores()
message(sprintf("Running %d method/size settings x %d iters on %d cores (N_total=%d)...",
                nrow(SETTING_GRID), N_ITER, n_cores, N_TOTAL))

cl <- makeCluster(n_cores)
clusterExport(cl, c("N_SITES", "D", "VAR_NMS", "K_COMP", "SITE_N", "Q_PROBS",
                    "EPS", "N_TOTAL", "DGP_PROPS", "SIZE_FACTORS", "SIZE_N0",
                    "METHODS", "SETTING_GRID", "SETTING_LEVELS",
                    "COMP_CENTERS", "COMP_COVS", "comp_cov", "site_offset",
                    "draw_site", "make_bias_vec", "global_from_clear",
                    "run_one_all",
                    # CLEAR.R:
                    "CLEAR_sim", "calc_empirical_truth", "site_fit_gmm",
                    "center_make_reference", "estimate_ipw_weights",
                    "weighted_eda_stats",
                    # CLEAR_Synthpop.R:
                    "CLEAR_synthpop_sim", "site_synthesize",
                    "center_make_reference_syn", "combine_global_stats"))
clusterEvalQ(cl, { library(MASS); library(mclust); library(glmnet); library(synthpop) })

seeds <- seq_len(N_ITER) + 7000L
res_list <- parLapply(cl, seeds, run_one_all)
stopCluster(cl)


# -----------------------------------------------------------------------------
# 3. Assemble results: one data.frame per setting Label, then bind.
# -----------------------------------------------------------------------------
all_rows <- list()
for (label in SETTING_LEVELS) {
  rows <- lapply(res_list, function(it) it[[label]])
  rows <- rows[!sapply(rows, is.null)]
  if (length(rows) == 0) next
  nm <- names(rows[[1]])
  M  <- do.call(rbind, lapply(rows, function(r) r[nm]))
  colnames(M) <- nm
  df <- as.data.frame(M, check.names = FALSE)
  df$Setting <- label
  all_rows[[label]] <- df
  message(sprintf("  %-24s : %d / %d iterations succeeded.", label, nrow(df), N_ITER))
}

common_cols <- Reduce(intersect, lapply(all_rows, colnames))
all_rows <- lapply(all_rows, function(df) df[, common_cols, drop = FALSE])
results  <- do.call(rbind, all_rows)


# -----------------------------------------------------------------------------
# 4. Reshape to long + summary CSV (same format as existing study)
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

# Attach Method / Size back onto each row (for coloring + summary columns).
lab2method <- setNames(SETTING_GRID$Method, SETTING_GRID$Label)
lab2size   <- setNames(SETTING_GRID$Size,   SETTING_GRID$Label)
long$Method <- factor(lab2method[long$Setting], levels = METHODS)
long$Size   <- factor(lab2size[long$Setting],   levels = names(SIZE_FACTORS))

stat_levels <- c("Mean", "Var", "Cov", "Quantile", "Lasso")
long$Setting   <- factor(long$Setting, levels = SETTING_LEVELS)
long$Statistic <- factor(long$Statistic, levels = stat_levels)
long <- long[is.finite(long$Bias), ]

summary_tbl <- aggregate(Bias ~ Setting + Method + Size + Statistic + Quantity,
                         data = long,
                         FUN = function(x) c(bias_mean = mean(x),
                                             bias_median = median(x),
                                             sd = sd(x),
                                             rmse = sqrt(mean(x^2))))
summary_tbl <- do.call(data.frame, summary_tbl)
print(summary_tbl, digits = 4, row.names = FALSE)
write.csv(summary_tbl, "CLEAR_3method_3size_bias_summary.csv", row.names = FALSE)


# -----------------------------------------------------------------------------
# 5. Per-statistic figures: x = setting (9 combos), fill = Method
# -----------------------------------------------------------------------------
stat_titles <- c(
  Mean     = "Mean bias: CLEAR vs synthpop vs CLEAR_synthpop (by reference size)",
  Var      = "Variance bias: CLEAR vs synthpop vs CLEAR_synthpop (by reference size)",
  Cov      = "Covariance bias: CLEAR vs synthpop vs CLEAR_synthpop (by reference size)",
  Quantile = "Quantile bias: CLEAR vs synthpop vs CLEAR_synthpop (by reference size)",
  Lasso    = "Lasso coefficient bias: CLEAR vs synthpop vs CLEAR_synthpop (by reference size)"
)

method_fills <- c(CLEAR = "#08519c", synthpop = "#e6550d", CLEAR_synthpop = "#31a354")

make_stat_plot <- function(st) {
  d <- long[long$Statistic == st, ]
  d$Quantity <- factor(d$Quantity, levels = unique(d$Quantity))
  
  ggplot(d, aes(x = Setting, y = Bias, fill = Method)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.3) +
    geom_boxplot(outlier.size = 0.3, alpha = 0.85, linewidth = 0.3) +
    facet_wrap(~ Quantity, nrow = 1) +
    scale_fill_manual(values = method_fills) +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 7),
          strip.text  = element_text(size = 8),
          legend.position = "top") +
    labs(title = stat_titles[[st]],
         x = "Setting: Method (reference size = half / equal / double of pooled N)",
         y = "Bias (estimate - truth)")
}

for (st in levels(long$Statistic)) {
  if (!any(long$Statistic == st)) next
  p_st <- make_stat_plot(st)
  n_q  <- length(unique(long$Quantity[long$Statistic == st]))
  ggsave(sprintf("CLEAR_3method_3size_bias_%s.png", st), p_st,
         width = max(8, n_q * 2.8), height = 4.6, dpi = 150, limitsize = FALSE)
  print(p_st)
}




# -----------------------------------------------------------------------------
# 6. Aggregate bias / SD / RMSE across quantities within each Statistic x Setting
#    so each metric is comparable on one panel per statistic.


summary_tbl <- read.csv("CLEAR_3method_3size_bias_summary.csv", stringsAsFactors = FALSE)
metric_tbl <- aggregate(
  cbind(AbsBias = abs(Bias.bias_mean), SD = Bias.sd, RMSE = Bias.rmse) ~
    Setting + Method + Size + Statistic,
  data = summary_tbl,
  FUN  = mean                      # average over the quantities within a statistic
)

# Long form for faceting: one row per (Setting, Statistic, Metric).
metric_long <- reshape(
  metric_tbl,
  varying   = c("AbsBias", "SD", "RMSE"),
  v.names   = "Value",
  timevar   = "Metric",
  times     = c("AbsBias", "SD", "RMSE"),
  direction = "long"
)
metric_long$Metric  <- factor(metric_long$Metric, levels = c("AbsBias", "SD", "RMSE"))
metric_long$Setting <- factor(metric_long$Setting, levels = SETTING_LEVELS)
metric_long$Method  <- factor(metric_long$Method,  levels = METHODS)

# -----------------------------------------------------------------------------
# 7. One figure per statistic: rows = metric (|bias| / SD / RMSE), x = setting
# -----------------------------------------------------------------------------
make_metric_plot <- function(st) {
  d <- metric_long[metric_long$Statistic == st, ]
  ggplot(d, aes(x = Setting, y = Value, fill = Method)) +
    geom_col(width = 0.7, alpha = 0.9) +
    facet_wrap(~ Metric, nrow = 1, scales = "free_y") +
    scale_fill_manual(values = method_fills) +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 7),
          legend.position = "top") +
    labs(title = sprintf("%s: |mean bias|, SD, and RMSE by method x reference size", st),
         x = "Setting: Method (reference size = half / equal / double of pooled N)",
         y = "Value (averaged over quantities)")
}

for (st in levels(long$Statistic)) {
  if (!any(metric_long$Statistic == st)) next
  p_m <- make_metric_plot(st)
  ggsave(sprintf("CLEAR_3method_3size_metrics_%s.png", st), p_m,
         width = 11, height = 4.2, dpi = 150, limitsize = FALSE)
  print(p_m)
}

# -----------------------------------------------------------------------------
# 8. Compact wide table: rows = Statistic x Setting, cols = the three metrics.
#    Good for dropping straight into a paper/appendix.
# -----------------------------------------------------------------------------
metric_wide <- reshape(
  metric_tbl[, c("Statistic", "Method", "Size", "Setting", "AbsBias", "SD", "RMSE")],
  idvar   = c("Statistic", "Setting", "Method", "Size"),
  timevar = NULL, direction = "wide"
)
metric_wide <- metric_wide[order(metric_wide$Statistic,
                                 match(metric_wide$Size, names(SIZE_FACTORS)),
                                 match(metric_wide$Method, METHODS)), ]

cat("\n===== |bias| / SD / RMSE (averaged over quantities) =====\n")
print(metric_wide, digits = 4, row.names = FALSE)
write.csv(metric_wide, "CLEAR_3method_3size_metric_summary.csv", row.names = FALSE)