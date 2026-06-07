# ============================================================
# 13_rebuild_weight_and_mask_missingness_v4.R
#
# PURPOSE
# -------
# Stable, no-imputation primary preprocessing for the 2023 masked workflow.
#
# FIXES IN THIS VERSION
# ---------------------
# 1) Never assumes aux_candidates exist in work_df.
# 2) Creates responded_it before response-model construction.
# 3) Screens response predictors to existing columns only.
# 4) Falls back to a minimal response model if auxiliary predictors are sparse.
# 5) Writes a clean weighted/masked bundle for downstream scripts.
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999)

suppressPackageStartupMessages({
  library(DBI)
  library(duckdb)
})

setwd("~/Desktop/AHA")

CONFIG_PATH <- file.path(getwd(), "metadata", "project_config.rds")
if (!file.exists(CONFIG_PATH)) stop("Run 00_setup_paths.R first.")
cfg <- readRDS(CONFIG_PATH)
DIRS <- cfg$DIRS

SEED <- 20260606L
TRIM_LOWER <- 0.01
TRIM_UPPER <- 0.99
MIN_FACTOR_COUNT <- 25L
MAX_FACTOR_LEVELS <- 12L
MAX_NUMERIC_PREDICTORS <- 12L
MIN_PCT_OBSERVED_KEEP <- 0.02

left_path    <- file.path(DIRS$derived_data, "aha_2023_integrated_left_asdb.rds")
profile_path <- file.path(DIRS$derived_data, "aha_2023_candidate_profile_matrix_left.rds")
asdb_sel_p   <- file.path(DIRS$tables, "selected_asdb_profile_vars.csv")
it_sel_p     <- file.path(DIRS$tables, "selected_it_profile_vars.csv")

req <- c(left_path, profile_path, asdb_sel_p, it_sel_p)
miss <- req[!file.exists(req)]
if (length(miss) > 0) stop("Missing required inputs:\n", paste(miss, collapse = "\n"))

left_df  <- readRDS(left_path)
prof_df  <- readRDS(profile_path)
asdb_sel <- read.csv(asdb_sel_p, stringsAsFactors = FALSE, check.names = FALSE)
it_sel   <- read.csv(it_sel_p, stringsAsFactors = FALSE, check.names = FALSE)

safe_lower <- function(x) tolower(ifelse(is.na(x), "", x))

normalize_name <- function(x) {
  x <- trimws(x)
  x <- gsub("\\.+", "_", x)
  x <- gsub("\\s+", "_", x)
  x <- gsub("[^A-Za-z0-9_]", "", x)
  safe_lower(x)
}

find_candidates <- function(df, patterns) {
  nm <- names(df)
  nmn <- normalize_name(nm)
  hits <- rep(FALSE, length(nm))
  for (p in patterns) hits <- hits | grepl(p, nmn)
  nm[hits]
}

collapse_factor_levels <- function(x, min_count = MIN_FACTOR_COUNT, max_levels = MAX_FACTOR_LEVELS) {
  x <- as.character(x)
  x[is.na(x) | trimws(x) == ""] <- "Missing"
  tab <- sort(table(x), decreasing = TRUE)
  keep <- names(tab)[tab >= min_count]
  if (length(keep) == 0L) {
    keep <- names(tab)[seq_len(min(length(tab), max_levels - 1L))]
  }
  keep <- keep[seq_len(min(length(keep), max_levels - 1L))]
  x[!(x %in% keep)] <- "Other"
  factor(x)
}

prep_ipw_design <- function(df, predictors) {
  if (length(predictors) == 0L) return(df[, integer(0), drop = FALSE])
  out <- df[, predictors, drop = FALSE]

  for (nm in names(out)) {
    x <- out[[nm]]

    if (is.character(x) || is.factor(x)) {
      out[[nm]] <- collapse_factor_levels(x)
    } else if (is.logical(x)) {
      out[[nm]] <- factor(ifelse(is.na(x), "Missing", ifelse(x, "Yes", "No")))
    } else if (is.numeric(x) || is.integer(x)) {
      x[!is.finite(x)] <- NA_real_
      med <- suppressWarnings(median(x, na.rm = TRUE))
      if (!is.finite(med)) med <- 0
      x[is.na(x)] <- med
      out[[nm]] <- as.numeric(x)
    } else {
      out[[nm]] <- collapse_factor_levels(as.character(x))
    }
  }

  keep1 <- vapply(out, function(x) length(unique(x)) > 1L, logical(1))
  out <- out[, keep1, drop = FALSE]

  if (ncol(out) > 1L) {
    dup <- rep(FALSE, ncol(out))
    for (j in seq_len(ncol(out) - 1L)) {
      if (dup[j]) next
      for (k in (j + 1L):ncol(out)) {
        if (identical(out[[j]], out[[k]])) dup[k] <- TRUE
      }
    }
    out <- out[, !dup, drop = FALSE]
  }

  is_num <- vapply(out, is.numeric, logical(1))
  if (sum(is_num) > MAX_NUMERIC_PREDICTORS) {
    num_names <- names(out)[is_num]
    vars <- vapply(out[num_names], function(z) {
      v <- suppressWarnings(var(z, na.rm = TRUE))
      if (!is.finite(v)) 0 else v
    }, numeric(1))
    keep_num <- names(sort(vars, decreasing = TRUE))[seq_len(MAX_NUMERIC_PREDICTORS)]
    out <- out[, c(names(out)[!is_num], keep_num), drop = FALSE]
  }

  out
}

fit_ipw_model_safe <- function(resp_df, outcome = "responded_it") {
  predictors <- setdiff(names(resp_df), outcome)

  if (length(predictors) == 0L) {
    p_hat <- rep(mean(resp_df[[outcome]] == 1, na.rm = TRUE), nrow(resp_df))
    return(list(engine = "intercept_only", fit = NULL, formula = as.formula(paste(outcome, "~ 1")), model_df = resp_df, p_hat = p_hat))
  }

  design_df <- prep_ipw_design(resp_df, predictors)

  if (ncol(design_df) == 0L) {
    p_hat <- rep(mean(resp_df[[outcome]] == 1, na.rm = TRUE), nrow(resp_df))
    return(list(engine = "intercept_only", fit = NULL, formula = as.formula(paste(outcome, "~ 1")), model_df = resp_df, p_hat = p_hat))
  }

  model_df <- cbind(resp_df[, outcome, drop = FALSE], design_df)
  form <- as.formula(paste(outcome, "~", paste(names(design_df), collapse = " + ")))

  fit <- try(
    glm(form, data = model_df, family = binomial(),
        control = glm.control(maxit = 30, epsilon = 1e-8, trace = FALSE)),
    silent = TRUE
  )

  if (!inherits(fit, "try-error")) {
    p_hat <- as.numeric(stats::predict(fit, type = "response"))
    return(list(engine = "glm", fit = fit, formula = form, model_df = model_df, p_hat = p_hat))
  }

  if (!requireNamespace("glmnet", quietly = TRUE)) {
    p_hat <- rep(mean(resp_df[[outcome]] == 1, na.rm = TRUE), nrow(resp_df))
    return(list(engine = "intercept_only_after_glm_fail", fit = NULL, formula = as.formula(paste(outcome, "~ 1")), model_df = resp_df, p_hat = p_hat))
  }

  mm <- stats::model.matrix(form, data = model_df)
  if (ncol(mm) <= 1L) {
    p_hat <- rep(mean(resp_df[[outcome]] == 1, na.rm = TRUE), nrow(resp_df))
    return(list(engine = "intercept_only_after_mm_fail", fit = NULL, formula = as.formula(paste(outcome, "~ 1")), model_df = resp_df, p_hat = p_hat))
  }

  x <- mm[, -1, drop = FALSE]
  y <- model_df[[outcome]]

  cvfit <- glmnet::cv.glmnet(
    x = x,
    y = y,
    family = "binomial",
    alpha = 0,
    nfolds = 5,
    type.measure = "deviance"
  )

  p_hat <- as.numeric(stats::predict(cvfit, newx = x, s = "lambda.1se", type = "response"))
  list(engine = "glmnet_ridge", fit = cvfit, formula = form, model_df = model_df, p_hat = p_hat)
}

selected_asdb <- unique(asdb_sel$var_original[asdb_sel$var_original %in% names(prof_df)])
selected_it   <- unique(it_sel$var_original[it_sel$var_original %in% names(prof_df)])

selected_vars <- unique(c(selected_asdb, selected_it))
selected_vars <- selected_vars[selected_vars %in% names(left_df)]

if (!("ID" %in% names(left_df))) stop("Expected ID column not found in left_df.")

obs_pct <- vapply(selected_vars, function(v) mean(!is.na(left_df[[v]])), numeric(1))
selected_vars <- selected_vars[obs_pct >= MIN_PCT_OBSERVED_KEEP]

work_df <- left_df[, c("ID", selected_vars), drop = FALSE]

it_cols_present <- intersect(selected_it, names(work_df))
if (length(it_cols_present) == 0L) {
  stop("No selected IT variables found in working frame; cannot define responded_it.")
}
work_df$responded_it <- as.integer(rowSums(!is.na(work_df[, it_cols_present, drop = FALSE])) > 0)

aux_candidates <- unique(c(
  find_candidates(left_df, c("^mstate$", "^state$", "region", "division")),
  find_candidates(left_df, c("bed", "hospbd", "bdtot")),
  find_candidates(left_df, c("teach", "teaching", "resid")),
  find_candidates(left_df, c("owner", "ownership", "control")),
  find_candidates(left_df, c("rural", "urban", "metro", "micropolitan", "cbsa")),
  find_candidates(left_df, c("system", "network", "member", "^sys"))
))

aux_candidates <- setdiff(aux_candidates, c("ID", selected_vars, "responded_it"))
aux_candidates <- aux_candidates[aux_candidates %in% names(left_df)]

if (length(aux_candidates) == 0L) {
  fallback_geo <- find_candidates(left_df, c("^mstate$", "^state$", "region"))
  aux_candidates <- fallback_geo[fallback_geo %in% names(left_df)]
}

if (length(aux_candidates) > 0L) {
  for (nm in aux_candidates) work_df[[nm]] <- left_df[[nm]]
}

resp_predictors <- intersect(aux_candidates, names(work_df))
if (length(resp_predictors) == 0L) {
  warning("No auxiliary response predictors available after screening; using intercept-only response model.")
}

resp_df <- work_df[, c("responded_it", resp_predictors), drop = FALSE]

set.seed(SEED)
ipw_obj <- fit_ipw_model_safe(resp_df, outcome = "responded_it")

p_obs <- mean(resp_df$responded_it == 1, na.rm = TRUE)
p_hat <- pmin(pmax(ipw_obj$p_hat, 1e-6), 1 - 1e-6)

raw_ipw <- ifelse(
  resp_df$responded_it == 1,
  p_obs / p_hat,
  (1 - p_obs) / (1 - p_hat)
)

q_lo <- unname(quantile(raw_ipw, probs = TRIM_LOWER, na.rm = TRUE))
q_hi <- unname(quantile(raw_ipw, probs = TRIM_UPPER, na.rm = TRUE))
ipw_trim <- pmin(pmax(raw_ipw, q_lo), q_hi)

work_df$ipw_response <- ipw_trim

variable_audit <- do.call(rbind, lapply(selected_vars, function(v) {
  x <- work_df[[v]]
  obs <- !is.na(x)
  pos <- sum(x == 1, na.rm = TRUE)
  zero <- sum(x == 0, na.rm = TRUE)
  data.frame(
    var = v,
    source = ifelse(v %in% selected_it, "IT", "ASDB"),
    n_total = nrow(work_df),
    n_observed = sum(obs),
    pct_observed = mean(obs),
    n_positive = pos,
    n_zero = zero,
    stringsAsFactors = FALSE
  )
}))

masked_binary <- work_df[, c("ID", "responded_it", "ipw_response", selected_vars), drop = FALSE]
obs_flags <- as.data.frame(lapply(masked_binary[, selected_vars, drop = FALSE], function(x) as.integer(!is.na(x))))
names(obs_flags) <- paste0(selected_vars, "__obs")

bundle <- list(
  id_var = "ID",
  full_frame_n = nrow(work_df),
  selected_vars = selected_vars,
  selected_asdb = intersect(selected_asdb, selected_vars),
  selected_it = intersect(selected_it, selected_vars),
  response_predictors = resp_predictors,
  response_engine = ipw_obj$engine,
  response_formula = deparse(ipw_obj$formula),
  response_fit = ipw_obj$fit,
  work_df = work_df[, c("ID", "responded_it", "ipw_response", resp_predictors, selected_vars), drop = FALSE],
  masked_binary = masked_binary,
  observed_flags = cbind(ID = work_df$ID, obs_flags),
  variable_audit = variable_audit
)

saveRDS(bundle, file.path(DIRS$derived_data, "aha_2023_masked_weight_bundle_v4.rds"))

write.csv(variable_audit, file.path(DIRS$tables, "masked_variable_audit_2023_v4.csv"), row.names = FALSE)

w_dist <- data.frame(
  metric = c("n", "mean_weight", "sd_weight", "min_weight", "p01", "median", "p99", "max_weight"),
  value = c(
    length(ipw_trim),
    mean(ipw_trim),
    sd(ipw_trim),
    min(ipw_trim),
    quantile(ipw_trim, 0.01),
    median(ipw_trim),
    quantile(ipw_trim, 0.99),
    max(ipw_trim)
  ),
  stringsAsFactors = FALSE
)

write.csv(w_dist, file.path(DIRS$tables, "ipw_weight_distribution_2023_masked_v4.csv"), row.names = FALSE)

if (!is.null(ipw_obj$fit) && inherits(ipw_obj$fit, "glm")) {
  coef_tbl <- data.frame(
    term = rownames(summary(ipw_obj$fit)$coefficients),
    summary(ipw_obj$fit)$coefficients,
    stringsAsFactors = FALSE
  )
} else if (!is.null(ipw_obj$fit) && "cv.glmnet" %in% class(ipw_obj$fit)) {
  beta <- as.matrix(stats::coef(ipw_obj$fit, s = "lambda.1se"))
  coef_tbl <- data.frame(
    term = rownames(beta),
    estimate = as.numeric(beta[, 1]),
    stringsAsFactors = FALSE
  )
} else {
  coef_tbl <- data.frame(
    term = "(Intercept)",
    estimate = qlogis(mean(resp_df$responded_it == 1, na.rm = TRUE)),
    stringsAsFactors = FALSE
  )
}

write.csv(coef_tbl, file.path(DIRS$tables, "ipw_response_model_summary_2023_masked_v4.csv"), row.names = FALSE)

duckdb_path <- file.path(DIRS$derived_data, "aha_masked_pipeline_cache_v4.duckdb")
con <- dbConnect(duckdb::duckdb(), duckdb_path)
dbWriteTable(con, "masked_work_df_v4", bundle$work_df, overwrite = TRUE)
dbWriteTable(con, "masked_variable_audit_2023_v4", variable_audit, overwrite = TRUE)
dbDisconnect(con, shutdown = TRUE)

log_file <- file.path(
  DIRS$logs,
  paste0("13_rebuild_weight_and_mask_missingness_v4_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt")
)

sink(log_file)
cat("13_rebuild_weight_and_mask_missingness_v4.R log\n")
cat("=======================================\n")
cat("Run time: ", as.character(Sys.time()), "\n\n")
cat("Full frame N: ", nrow(work_df), "\n")
cat("Selected vars kept: ", length(selected_vars), "\n")
cat("Selected IT vars present: ", length(it_cols_present), "\n")
cat("Auxiliary response predictors kept: ", length(resp_predictors), "\n")
cat("Response engine: ", ipw_obj$engine, "\n\n")
cat("Response formula:\n")
print(ipw_obj$formula)
cat("\nWeight distribution:\n")
print(w_dist)
cat("\nTop variable audit rows:\n")
print(utils::head(variable_audit, 20))
cat("\nSession info:\n")
print(sessionInfo())
sink()

cat("Script 13 v4 complete.\n")
