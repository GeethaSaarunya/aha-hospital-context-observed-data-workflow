# ============================================================
# 14_rebuild_masked_smc_network_v4.R
#
# PURPOSE
# -------
# Build a stable masked Simple Matching Coefficient (SMC) network
# from the v4 weighted/masked bundle, explicitly addressing the
# main potatoes observed in earlier runs:
#
# 1) ultra-sparse variables were destabilizing pairwise overlap
# 2) pairwise overlap could be ~0, yielding NaN / Inf summaries
# 3) weak-overlap pairs should not drive the similarity graph
# 4) module derivation should use a stable variable subset only
#
# INPUT
# -----
# - derived_data/aha_2023_masked_weight_bundle_v4.rds
#
# OUTPUTS
# -------
# - derived_data/aha_2023_masked_smc_bundle_v4.rds
# - tables/masked_smc_similarity_summary_2023_v4.csv
# - tables/masked_smc_overlap_summary_2023_v4.csv
# - tables/masked_smc_variable_filtering_2023_v4.csv
# - tables/masked_smc_module_assignments_2023_v4.csv
# - tables/masked_smc_module_selection_2023_v4.csv
# - logs/14_rebuild_masked_smc_network_v4_<timestamp>.txt
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999)

setwd("~/Desktop/AHA")

CONFIG_PATH <- file.path(getwd(), "metadata", "project_config.rds")
if (!file.exists(CONFIG_PATH)) stop("Run 00_setup_paths.R first.")
cfg <- readRDS(CONFIG_PATH)
DIRS <- cfg$DIRS

# ----------------------------
# Parameters
# ----------------------------
MIN_VAR_OBSERVED_PCT <- 0.20      # drop ultra-sparse variables
MIN_VAR_POSITIVE_N   <- 20L       # drop nearly all-zero variables
MIN_VAR_ZERO_N       <- 20L       # drop nearly all-one variables

MIN_PAIRWISE_OVERLAP_N <- 250L    # require enough jointly observed hospitals
OVERLAP_SHRINK_TARGET  <- 1000L   # shrink weak-overlap pairs toward 0
MIN_SMC_FOR_EDGE       <- 0.55    # only keep reasonably similar pairs as graph edges

MODULE_K_GRID <- 2:10             # candidate cutree sizes
SEED <- 20260606L

# ----------------------------
# Inputs
# ----------------------------
bundle_path <- file.path(DIRS$derived_data, "aha_2023_masked_weight_bundle_v4.rds")
if (!file.exists(bundle_path)) {
  stop("Missing v4 masked-weight bundle. Run 13_rebuild_weight_and_mask_missingness_v4.R first.")
}

bundle <- readRDS(bundle_path)
masked_binary <- bundle$masked_binary
var_audit <- bundle$variable_audit

if (!("ID" %in% names(masked_binary))) stop("Expected ID in masked_binary.")
if (!("ipw_response" %in% names(masked_binary))) stop("Expected ipw_response in masked_binary.")

# ----------------------------
# Helpers
# ----------------------------
weighted_mean_safe <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & !is.na(x) & !is.na(w)
  if (!any(ok)) return(NA_real_)
  sum(x[ok] * w[ok]) / sum(w[ok])
}

summarize_vec <- function(x) {
  x <- x[is.finite(x) & !is.na(x)]
  if (length(x) == 0L) {
    return(c(n = 0, min = NA, p01 = NA, p25 = NA, median = NA, mean = NA, p75 = NA, p99 = NA, max = NA))
  }
  c(
    n = length(x),
    min = min(x),
    p01 = unname(quantile(x, 0.01)),
    p25 = unname(quantile(x, 0.25)),
    median = median(x),
    mean = mean(x),
    p75 = unname(quantile(x, 0.75)),
    p99 = unname(quantile(x, 0.99)),
    max = max(x)
  )
}

# pairwise masked weighted SMC for binary variables
pairwise_masked_weighted_smc <- function(x, y, w,
                                         min_overlap_n = MIN_PAIRWISE_OVERLAP_N,
                                         shrink_target = OVERLAP_SHRINK_TARGET) {
  ok <- !is.na(x) & !is.na(y) & is.finite(w) & !is.na(w)
  overlap_n <- sum(ok)

  if (overlap_n < min_overlap_n) {
    return(c(overlap_n = overlap_n,
             overlap_w = ifelse(overlap_n == 0, NA_real_, sum(w[ok])),
             smc_raw = NA_real_,
             smc_shrunk = NA_real_))
  }

  x_ok <- x[ok]
  y_ok <- y[ok]
  w_ok <- w[ok]

  # restrict to binary 0/1 observed pairs only
  valid_bin <- x_ok %in% c(0, 1) & y_ok %in% c(0, 1)
  if (!any(valid_bin)) {
    return(c(overlap_n = overlap_n,
             overlap_w = sum(w_ok),
             smc_raw = NA_real_,
             smc_shrunk = NA_real_))
  }

  x_ok <- x_ok[valid_bin]
  y_ok <- y_ok[valid_bin]
  w_ok <- w_ok[valid_bin]

  overlap_n_bin <- length(x_ok)
  overlap_w <- sum(w_ok)

  if (overlap_n_bin < min_overlap_n || overlap_w <= 0) {
    return(c(overlap_n = overlap_n_bin,
             overlap_w = overlap_w,
             smc_raw = NA_real_,
             smc_shrunk = NA_real_))
  }

  matches <- as.numeric(x_ok == y_ok)
  smc_raw <- sum(matches * w_ok) / overlap_w

  shrink_factor <- min(1, overlap_n_bin / shrink_target)
  smc_shrunk <- smc_raw * shrink_factor

  c(
    overlap_n = overlap_n_bin,
    overlap_w = overlap_w,
    smc_raw = smc_raw,
    smc_shrunk = smc_shrunk
  )
}

compute_average_silhouette <- function(dist_obj, membership) {
  if (length(unique(membership)) < 2L) return(NA_real_)
  sil <- try(cluster::silhouette(membership, dist_obj), silent = TRUE)
  if (inherits(sil, "try-error")) return(NA_real_)
  mean(sil[, "sil_width"])
}

# ----------------------------
# Stable variable subset
# ----------------------------
candidate_vars <- setdiff(names(masked_binary), c("ID", "responded_it", "ipw_response"))
candidate_vars <- intersect(candidate_vars, var_audit$var)

vf <- merge(
  data.frame(var = candidate_vars, stringsAsFactors = FALSE),
  var_audit,
  by = "var",
  all.x = TRUE,
  sort = FALSE
)

vf$keep_obs <- vf$pct_observed >= MIN_VAR_OBSERVED_PCT
vf$keep_pos <- vf$n_positive >= MIN_VAR_POSITIVE_N
vf$keep_zero <- vf$n_zero >= MIN_VAR_ZERO_N
vf$keep_all <- vf$keep_obs & vf$keep_pos & vf$keep_zero

stable_vars <- vf$var[vf$keep_all]

if (length(stable_vars) < 10L) {
  stop("Too few stable variables survived filtering. Loosen thresholds or inspect variable audit.")
}

vf$status <- ifelse(vf$keep_all, "kept", "dropped")
vf$drop_reason <- ""
vf$drop_reason[!vf$keep_obs] <- paste0(vf$drop_reason[!vf$keep_obs], "low_obs;")
vf$drop_reason[!vf$keep_pos] <- paste0(vf$drop_reason[!vf$keep_pos], "low_pos;")
vf$drop_reason[!vf$keep_zero] <- paste0(vf$drop_reason[!vf$keep_zero], "low_zero;")

write.csv(vf, file.path(DIRS$tables, "masked_smc_variable_filtering_2023_v4.csv"), row.names = FALSE)

# ----------------------------
# Pairwise masked SMC
# ----------------------------
X <- masked_binary[, stable_vars, drop = FALSE]
W <- masked_binary$ipw_response

p <- ncol(X)
smc_raw_mat <- matrix(NA_real_, p, p, dimnames = list(stable_vars, stable_vars))
smc_shrunk_mat <- matrix(NA_real_, p, p, dimnames = list(stable_vars, stable_vars))
overlap_n_mat <- matrix(0, p, p, dimnames = list(stable_vars, stable_vars))
overlap_w_mat <- matrix(NA_real_, p, p, dimnames = list(stable_vars, stable_vars))

for (i in seq_len(p)) {
  xi <- X[[i]]
  for (j in i:p) {
    yj <- X[[j]]
    out <- pairwise_masked_weighted_smc(xi, yj, W)

    overlap_n_mat[i, j] <- overlap_n_mat[j, i] <- out["overlap_n"]
    overlap_w_mat[i, j] <- overlap_w_mat[j, i] <- out["overlap_w"]
    smc_raw_mat[i, j] <- smc_raw_mat[j, i] <- out["smc_raw"]
    smc_shrunk_mat[i, j] <- smc_shrunk_mat[j, i] <- out["smc_shrunk"]
  }
}

diag(smc_raw_mat) <- 1
diag(smc_shrunk_mat) <- 1

# replace any remaining non-finite off-diagonals with 0 similarity
bad_offdiag <- !is.finite(smc_shrunk_mat)
diag(bad_offdiag) <- FALSE
smc_shrunk_mat[bad_offdiag] <- 0

bad_offdiag_raw <- !is.finite(smc_raw_mat)
diag(bad_offdiag_raw) <- FALSE
smc_raw_mat[bad_offdiag_raw] <- 0

# ----------------------------
# Summaries
# ----------------------------
upper_idx <- upper.tri(smc_shrunk_mat, diag = FALSE)

smc_similarity_summary <- data.frame(
  metric = names(summarize_vec(smc_shrunk_mat[upper_idx])),
  value = as.numeric(summarize_vec(smc_shrunk_mat[upper_idx])),
  stringsAsFactors = FALSE
)

smc_overlap_summary <- data.frame(
  metric = names(summarize_vec(overlap_n_mat[upper_idx])),
  value = as.numeric(summarize_vec(overlap_n_mat[upper_idx])),
  stringsAsFactors = FALSE
)

write.csv(
  smc_similarity_summary,
  file.path(DIRS$tables, "masked_smc_similarity_summary_2023_v4.csv"),
  row.names = FALSE
)

write.csv(
  smc_overlap_summary,
  file.path(DIRS$tables, "masked_smc_overlap_summary_2023_v4.csv"),
  row.names = FALSE
)

# ----------------------------
# Graph pruning and module derivation
# ----------------------------
# For clustering, convert similarity to a bounded distance.
# Keep only stable similarities above threshold for graph summaries,
# but use full shrunk matrix for hierarchical clustering.
graph_edge_count <- sum(smc_shrunk_mat[upper_idx] >= MIN_SMC_FOR_EDGE, na.rm = TRUE)

dist_mat <- 1 - smc_shrunk_mat
diag(dist_mat) <- 0
dist_obj <- as.dist(dist_mat)

hc <- hclust(dist_obj, method = "average")

module_sel <- data.frame(
  k = integer(0),
  avg_silhouette = numeric(0),
  min_module_size = integer(0),
  max_module_size = integer(0),
  n_singletons = integer(0),
  stringsAsFactors = FALSE
)

for (k in MODULE_K_GRID) {
  mem <- cutree(hc, k = k)
  tab <- table(mem)
  module_sel <- rbind(
    module_sel,
    data.frame(
      k = k,
      avg_silhouette = compute_average_silhouette(dist_obj, mem),
      min_module_size = min(tab),
      max_module_size = max(tab),
      n_singletons = sum(tab == 1),
      stringsAsFactors = FALSE
    )
  )
}

# choose k with highest silhouette, then fewer singletons, then larger min module
ord <- order(-module_sel$avg_silhouette, module_sel$n_singletons, -module_sel$min_module_size, module_sel$k)
best_k <- module_sel$k[ord[1]]

module_membership <- cutree(hc, k = best_k)

module_assignments <- data.frame(
  var = stable_vars,
  module = paste0("M", module_membership),
  stringsAsFactors = FALSE
)

write.csv(
  module_assignments,
  file.path(DIRS$tables, "masked_smc_module_assignments_2023_v4.csv"),
  row.names = FALSE
)

write.csv(
  module_sel,
  file.path(DIRS$tables, "masked_smc_module_selection_2023_v4.csv"),
  row.names = FALSE
)

# ----------------------------
# Save bundle
# ----------------------------
smc_bundle <- list(
  stable_vars = stable_vars,
  variable_filtering = vf,
  smc_raw = smc_raw_mat,
  smc_shrunk = smc_shrunk_mat,
  overlap_n = overlap_n_mat,
  overlap_w = overlap_w_mat,
  hclust = hc,
  module_selection = module_sel,
  best_k = best_k,
  module_assignments = module_assignments,
  graph_edge_count = graph_edge_count,
  min_pairwise_overlap_n = MIN_PAIRWISE_OVERLAP_N,
  overlap_shrink_target = OVERLAP_SHRINK_TARGET,
  min_smc_for_edge = MIN_SMC_FOR_EDGE
)

saveRDS(
  smc_bundle,
  file.path(DIRS$derived_data, "aha_2023_masked_smc_bundle_v4.rds")
)

# ----------------------------
# Log
# ----------------------------
log_file <- file.path(
  DIRS$logs,
  paste0("14_rebuild_masked_smc_network_v4_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt")
)

sink(log_file)
cat("14_rebuild_masked_smc_network_v4.R log\n")
cat("=======================================\n")
cat("Run time: ", as.character(Sys.time()), "\n\n")
cat("Stable vars kept: ", length(stable_vars), " of ", length(candidate_vars), "\n")
cat("Minimum variable observed pct: ", MIN_VAR_OBSERVED_PCT, "\n")
cat("Minimum pairwise overlap N: ", MIN_PAIRWISE_OVERLAP_N, "\n")
cat("Overlap shrink target: ", OVERLAP_SHRINK_TARGET, "\n")
cat("Minimum SMC edge threshold: ", MIN_SMC_FOR_EDGE, "\n")
cat("Graph edge count (thresholded): ", graph_edge_count, "\n")
cat("Selected module count (best_k): ", best_k, "\n\n")
cat("Similarity summary:\n")
print(smc_similarity_summary)
cat("\nOverlap summary:\n")
print(smc_overlap_summary)
cat("\nModule selection table:\n")
print(module_sel)
cat("\nTop dropped variables:\n")
print(utils::head(vf[vf$status == "dropped", c("var", "pct_observed", "n_positive", "n_zero", "drop_reason")], 20))
cat("\nSession info:\n")
print(sessionInfo())
sink()

cat("Script 14 v4 complete.\n")
