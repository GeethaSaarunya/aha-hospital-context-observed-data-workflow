# ============================================================
# 15_rebuild_weighted_masked_clustering_v5.R
#
# PURPOSE
# -------
# Rebuild hospital grouping from the v4 masked-SMC module structure.
#
# THIS VERSION FIXES
# ------------------
# Error seen in v4:
#   Error in if (all(cl == cl_prev)) break :
#     missing value where TRUE/FALSE needed
#
# Why that happened:
# - some rows in x or centers could still contain NA / non-finite values
# - some distance rows therefore became all-NA / all-Inf
# - cluster assignment then produced NA values
#
# Fixes here:
# 1) explicit sanitization of module-score matrix before clustering
# 2) drop zero-variance score columns before scaling
# 3) replace any non-finite scaled values with 0
# 4) robust distance handling inside weighted k-means
# 5) deterministic fallback when a row has invalid distances
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999)

suppressPackageStartupMessages({
  library(cluster)
})

setwd("~/Desktop/AHA")

CONFIG_PATH <- file.path(getwd(), "metadata", "project_config.rds")
if (!file.exists(CONFIG_PATH)) stop("Run 00_setup_paths.R first.")
cfg <- readRDS(CONFIG_PATH)
DIRS <- cfg$DIRS

# ----------------------------
# Parameters
# ----------------------------
MIN_USABLE_MODULE_SIZE <- 3L
MIN_MODULES_OBSERVED_FOR_ELIGIBILITY <- 2L
MIN_MEAN_COMPLETENESS_FOR_ELIGIBILITY <- 0.25
K_GRID <- 2:5
SEED <- 20260606L
MAX_KMEANS_ITERS <- 100L
N_STARTS <- 10L

# ----------------------------
# Inputs
# ----------------------------
mask_bundle_path <- file.path(DIRS$derived_data, "aha_2023_masked_weight_bundle_v4.rds")
smc_bundle_path  <- file.path(DIRS$derived_data, "aha_2023_masked_smc_bundle_v4.rds")

if (!file.exists(mask_bundle_path)) stop("Run Script 13 v4 first.")
if (!file.exists(smc_bundle_path))  stop("Run Script 14 v4 first.")

mask_bundle <- readRDS(mask_bundle_path)
smc_bundle  <- readRDS(smc_bundle_path)

masked_binary <- mask_bundle$masked_binary
smc_mat <- smc_bundle$smc_shrunk
module_assignments <- smc_bundle$module_assignments

if (!("ID" %in% names(masked_binary))) stop("Expected ID in masked_binary.")
if (!("ipw_response" %in% names(masked_binary))) stop("Expected ipw_response in masked_binary.")

# ----------------------------
# Helpers
# ----------------------------
weighted_median <- function(x, w) {
  ok <- is.finite(x) & !is.na(x) & is.finite(w) & !is.na(w)
  x <- x[ok]
  w <- w[ok]
  if (length(x) == 0L) return(NA_real_)
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  cw <- cumsum(w) / sum(w)
  x[which(cw >= 0.5)[1]]
}

weighted_mean_safe <- function(x, w = NULL) {
  ok <- !is.na(x) & is.finite(x)
  if (!any(ok)) return(NA_real_)
  if (is.null(w)) return(mean(x[ok]))
  wok <- ok & !is.na(w) & is.finite(w)
  if (!any(wok)) return(mean(x[ok]))
  sum(x[wok] * w[wok]) / sum(w[wok])
}

compute_avg_silhouette <- function(xmat, membership) {
  if (length(unique(membership)) < 2L) return(NA_real_)
  d <- try(dist(xmat), silent = TRUE)
  if (inherits(d, "try-error")) return(NA_real_)
  sil <- try(cluster::silhouette(membership, d), silent = TRUE)
  if (inherits(sil, "try-error")) return(NA_real_)
  mean(sil[, "sil_width"])
}

size_balance_index <- function(tab) {
  tab <- as.numeric(tab)
  if (length(tab) == 0L || sum(tab) == 0L) return(NA_real_)
  min(tab) / max(tab)
}

sanitize_matrix_for_clustering <- function(df, weights) {
  x <- as.data.frame(df)

  # replace non-finite / NA by weighted medians column-wise
  for (nm in names(x)) {
    vec <- as.numeric(x[[nm]])
    med <- weighted_median(vec, weights)
    if (!is.finite(med) || is.na(med)) med <- 0
    bad <- !is.finite(vec) | is.na(vec)
    vec[bad] <- med
    x[[nm]] <- vec
  }

  # drop zero-variance columns
  sds <- vapply(x, function(z) stats::sd(z, na.rm = TRUE), numeric(1))
  keep <- is.finite(sds) & !is.na(sds) & sds > 0
  x <- x[, keep, drop = FALSE]

  if (ncol(x) < 2L) {
    stop("Too few non-degenerate score columns remain after sanitization.")
  }

  # scale and replace any non-finite remnants with 0
  xs <- scale(as.matrix(x))
  xs[!is.finite(xs)] <- 0

  list(
    xmat = xs,
    kept_cols = colnames(x)
  )
}

weighted_kmeans_once <- function(x, w, centers_init, max_iter = MAX_KMEANS_ITERS) {
  n <- nrow(x)
  k <- nrow(centers_init)

  centers <- centers_init
  centers[!is.finite(centers)] <- 0

  cl_prev <- rep(-1L, n)

  for (iter in seq_len(max_iter)) {
    dmat <- matrix(Inf, n, k)

    for (j in seq_len(k)) {
      dif <- sweep(x, 2, centers[j, ], FUN = "-")
      d <- rowSums(dif^2)
      d[!is.finite(d)] <- Inf
      dmat[, j] <- d
    }

    # robust assignment: if all distances invalid, assign to largest-weight cluster seed 1
    cl <- apply(dmat, 1, function(z) {
      z[!is.finite(z)] <- Inf
      if (all(is.infinite(z))) {
        return(1L)
      } else {
        return(which.min(z))
      }
    })
    cl <- as.integer(cl)

    if (all(cl == cl_prev)) break
    cl_prev <- cl

    for (j in seq_len(k)) {
      idx <- which(cl == j)
      if (length(idx) == 0L) {
        repl <- sample.int(n, 1, prob = w)
        centers[j, ] <- x[repl, ]
      } else {
        ww <- w[idx]
        if (!all(is.finite(ww)) || sum(ww, na.rm = TRUE) <= 0) ww <- rep(1, length(idx))
        centers[j, ] <- colSums(x[idx, , drop = FALSE] * ww) / sum(ww)
      }
    }

    centers[!is.finite(centers)] <- 0
  }

  d_final <- numeric(n)
  for (i in seq_len(n)) {
    dif <- x[i, ] - centers[cl[i], ]
    dval <- sum(dif^2)
    if (!is.finite(dval) || is.na(dval)) dval <- 0
    d_final[i] <- dval
  }
  tot_withinss_w <- sum(w * d_final)

  list(cluster = cl, centers = centers, tot_withinss_w = tot_withinss_w, iterations = iter)
}

weighted_kmeans_multi <- function(x, w, k, nstart = N_STARTS, seed = SEED) {
  set.seed(seed + k)
  n <- nrow(x)

  best <- NULL
  for (s in seq_len(nstart)) {
    init_idx <- sample.int(n, k, prob = w, replace = FALSE)
    fit <- weighted_kmeans_once(x, w, centers_init = x[init_idx, , drop = FALSE])
    if (is.null(best) || fit$tot_withinss_w < best$tot_withinss_w) {
      best <- fit
    }
  }
  best
}

# ----------------------------
# Resolve tiny modules
# ----------------------------
module_counts <- sort(table(module_assignments$module), decreasing = TRUE)
usable_modules <- names(module_counts)[module_counts >= MIN_USABLE_MODULE_SIZE]
tiny_modules   <- names(module_counts)[module_counts < MIN_USABLE_MODULE_SIZE]

module_assignments$module_initial <- module_assignments$module
module_assignments$module_final <- module_assignments$module

if (length(usable_modules) == 0L) {
  stop("No usable modules found after size filtering. Inspect Script 14 output.")
}

if (length(tiny_modules) > 0L) {
  for (tm in tiny_modules) {
    tiny_vars <- module_assignments$var[module_assignments$module_initial == tm]
    if (length(tiny_vars) == 0L) next

    target_scores <- sapply(usable_modules, function(um) {
      big_vars <- module_assignments$var[module_assignments$module_initial == um]
      if (length(big_vars) == 0L) return(-Inf)
      sims <- as.numeric(smc_mat[tiny_vars, big_vars, drop = FALSE])
      sims <- sims[is.finite(sims) & !is.na(sims)]
      if (length(sims) == 0L) return(-Inf)
      mean(sims)
    })

    best_um <- names(which.max(target_scores))
    if (length(best_um) == 0L || !is.finite(max(target_scores))) {
      module_assignments$module_final[module_assignments$module_initial == tm] <- "M0_residual"
    } else {
      module_assignments$module_final[module_assignments$module_initial == tm] <- best_um
    }
  }
}

usable_final_modules <- sort(unique(module_assignments$module_final[module_assignments$module_final != "M0_residual"]))

write.csv(
  module_assignments[, c("var", "module_initial", "module_final")],
  file.path(DIRS$tables, "masked_module_assignments_usable_2023_v5.csv"),
  row.names = FALSE
)

# ----------------------------
# Build hospital-level module scores
# ----------------------------
score_df <- data.frame(
  ID = masked_binary$ID,
  ipw_response = masked_binary$ipw_response,
  responded_it = masked_binary$responded_it,
  stringsAsFactors = FALSE
)

comp_long <- list()

for (mod in usable_final_modules) {
  vars_mod <- module_assignments$var[module_assignments$module_final == mod]
  vars_mod <- intersect(vars_mod, names(masked_binary))
  if (length(vars_mod) == 0L) next

  xmat_mod <- masked_binary[, vars_mod, drop = FALSE]
  obs_count <- rowSums(!is.na(xmat_mod))
  mod_size <- ncol(xmat_mod)
  completeness <- obs_count / mod_size

  score <- apply(xmat_mod, 1, function(z) {
    ok <- !is.na(z)
    if (!any(ok)) return(NA_real_)
    mean(as.numeric(z[ok]))
  })

  score_df[[paste0(mod, "_score")]] <- score
  score_df[[paste0(mod, "_obs_n")]] <- obs_count
  score_df[[paste0(mod, "_complete")]] <- completeness

  comp_long[[mod]] <- data.frame(
    ID = masked_binary$ID,
    module = mod,
    module_size = mod_size,
    obs_n = obs_count,
    completeness = completeness,
    score = score,
    stringsAsFactors = FALSE
  )
}

comp_df <- do.call(rbind, comp_long)
write.csv(
  comp_df,
  file.path(DIRS$tables, "masked_module_score_completeness_2023_v5.csv"),
  row.names = FALSE
)

score_cols <- grep("_score$", names(score_df), value = TRUE)
complete_cols <- grep("_complete$", names(score_df), value = TRUE)

if (length(score_cols) < 2L) {
  stop("Too few usable module score columns were created for clustering.")
}

score_df$n_modules_observed <- rowSums(!is.na(score_df[, score_cols, drop = FALSE]))
score_df$mean_completeness <- apply(score_df[, complete_cols, drop = FALSE], 1, function(z) {
  if (all(is.na(z))) return(NA_real_)
  mean(z, na.rm = TRUE)
})

score_df$cluster_eligible <- with(
  score_df,
  n_modules_observed >= MIN_MODULES_OBSERVED_FOR_ELIGIBILITY &
    mean_completeness >= MIN_MEAN_COMPLETENESS_FOR_ELIGIBILITY
)

eligible_df <- score_df[score_df$cluster_eligible, c("ID", "ipw_response", score_cols), drop = FALSE]
if (nrow(eligible_df) < 50L) {
  stop("Too few eligible hospitals for clustering after eligibility filtering.")
}

# sanitize before clustering
san <- sanitize_matrix_for_clustering(eligible_df[, score_cols, drop = FALSE], eligible_df$ipw_response)
xmat <- san$xmat
kept_score_cols <- san$kept_cols
w <- eligible_df$ipw_response
w[!is.finite(w) | is.na(w) | w <= 0] <- 1

# ----------------------------
# Weighted clustering
# ----------------------------
fit_tbl <- data.frame(
  k = integer(0),
  n_eligible = integer(0),
  n_assigned = integer(0),
  n_score_cols_used = integer(0),
  weighted_withinss = numeric(0),
  avg_silhouette = numeric(0),
  min_cluster_size = integer(0),
  max_cluster_size = integer(0),
  size_balance = numeric(0),
  stringsAsFactors = FALSE
)

fit_objs <- list()

for (k in K_GRID) {
  fit <- weighted_kmeans_multi(xmat, w, k = k, nstart = N_STARTS, seed = SEED)
  tab <- table(fit$cluster)

  fit_tbl <- rbind(
    fit_tbl,
    data.frame(
      k = k,
      n_eligible = nrow(eligible_df),
      n_assigned = length(fit$cluster),
      n_score_cols_used = ncol(xmat),
      weighted_withinss = fit$tot_withinss_w,
      avg_silhouette = compute_avg_silhouette(xmat, fit$cluster),
      min_cluster_size = min(tab),
      max_cluster_size = max(tab),
      size_balance = size_balance_index(tab),
      stringsAsFactors = FALSE
    )
  )

  fit_objs[[as.character(k)]] <- fit
}

ord <- order(-fit_tbl$avg_silhouette, fit_tbl$weighted_withinss, -fit_tbl$size_balance, fit_tbl$k)
best_k <- fit_tbl$k[ord[1]]
best_fit <- fit_objs[[as.character(best_k)]]

write.csv(
  fit_tbl,
  file.path(DIRS$tables, "masked_weighted_clustering_fit_2023_v5.csv"),
  row.names = FALSE
)

# ----------------------------
# Assign profiles back to full frame
# ----------------------------
profile_assign <- data.frame(
  ID = score_df$ID,
  cluster_eligible = score_df$cluster_eligible,
  n_modules_observed = score_df$n_modules_observed,
  mean_completeness = score_df$mean_completeness,
  profile = NA_character_,
  stringsAsFactors = FALSE
)

profile_assign$profile[match(eligible_df$ID, profile_assign$ID)] <- paste0("P", best_fit$cluster)

profile_assign <- merge(
  profile_assign,
  score_df[, c("ID", score_cols, complete_cols), drop = FALSE],
  by = "ID", all.x = TRUE, sort = FALSE
)

profile_sizes <- as.data.frame(table(profile_assign$profile), stringsAsFactors = FALSE)
names(profile_sizes) <- c("profile", "n_hospitals")
profile_sizes <- profile_sizes[!is.na(profile_sizes$profile) & profile_sizes$profile != "", , drop = FALSE]

write.csv(
  profile_assign,
  file.path(DIRS$tables, "masked_profile_assignments_2023_v5.csv"),
  row.names = FALSE
)

write.csv(
  profile_sizes,
  file.path(DIRS$tables, "masked_profile_sizes_2023_v5.csv"),
  row.names = FALSE
)

profile_score_means <- NULL
if (nrow(profile_sizes) > 0L) {
  plev <- sort(unique(na.omit(profile_assign$profile)))
  profile_score_means <- do.call(rbind, lapply(plev, function(pp) {
    d <- profile_assign[profile_assign$profile == pp, , drop = FALSE]
    out <- data.frame(profile = pp, n_hospitals = nrow(d), stringsAsFactors = FALSE)
    for (cc in score_cols) {
      out[[cc]] <- weighted_mean_safe(d[[cc]])
    }
    out
  }))
}

clust_bundle <- list(
  module_assignments_usable = module_assignments,
  usable_final_modules = usable_final_modules,
  module_score_df = score_df,
  kept_score_cols = kept_score_cols,
  eligible_df = eligible_df,
  fit_table = fit_tbl,
  best_k = best_k,
  best_fit = best_fit,
  profile_assignments = profile_assign,
  profile_sizes = profile_sizes,
  profile_score_means = profile_score_means
)

saveRDS(
  clust_bundle,
  file.path(DIRS$derived_data, "aha_2023_masked_clustering_bundle_v5.rds")
)

# ----------------------------
# Log
# ----------------------------
log_file <- file.path(
  DIRS$logs,
  paste0("15_rebuild_weighted_masked_clustering_v5_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt")
)

sink(log_file)
cat("15_rebuild_weighted_masked_clustering_v5.R log\n")
cat("=======================================\n")
cat("Run time: ", as.character(Sys.time()), "\n\n")
cat("Initial modules:\n")
print(module_counts)
cat("\nFinal module table after tiny-module handling:\n")
print(table(module_assignments$module_final))
cat("\nScore columns kept for clustering:\n")
print(kept_score_cols)
cat("\nEligible hospitals: ", nrow(eligible_df), " of ", nrow(score_df), "\n")
cat("Selected K: ", best_k, "\n\n")
cat("Clustering fit table:\n")
print(fit_tbl)
cat("\nProfile sizes:\n")
print(profile_sizes)
cat("\nSession info:\n")
print(sessionInfo())
sink()

cat("Script 15 v5 complete.\n")
