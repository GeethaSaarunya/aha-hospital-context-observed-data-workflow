# ============================================================
# 16_rebuild_system_nesting_and_validation_masked_v7.R
#
# PURPOSE
# -------
# Post-clustering evaluation for the v5 masked weighted pipeline.
#
# WHY v7 EXISTS
# -------------
# v6 incorrectly selected SYSTEM_PRIMARY_CONTACT as the "system" field,
# producing person-name groupings instead of health-system groupings.
#
# v7 fixes this by:
#   1) EXCLUDING obvious contact/name/person fields
#   2) PRIORITIZING true system-like ID/member variables
#   3) SKIPPING system nesting entirely if no believable system field exists
#   4) STILL writing clean outputs and optional outcome-validation placeholders
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
OUTCOME_TEMPLATE_NAME <- "external_outcomes_2023_template_v7.csv"
MIN_GROUPS_FOR_ICC <- 2L
MIN_OBS_FOR_ICC <- 10L

# ----------------------------
# Inputs
# ----------------------------
left_path <- file.path(DIRS$derived_data, "aha_2023_integrated_left_asdb.rds")
clust_bundle_path <- file.path(DIRS$derived_data, "aha_2023_masked_clustering_bundle_v5.rds")

if (!file.exists(left_path)) stop("Missing aha_2023_integrated_left_asdb.rds")
if (!file.exists(clust_bundle_path)) stop("Run Script 15 v5 first.")

left_df <- readRDS(left_path)
clust_bundle <- readRDS(clust_bundle_path)

core_df <- clust_bundle$profile_assignments
if (is.null(core_df) || nrow(core_df) == 0L) stop("No profile assignments found in clustering bundle.")
if (!("ID" %in% names(core_df))) stop("Expected ID column in profile assignments.")
if (!("profile" %in% names(core_df))) stop("Expected profile column in profile assignments.")

dir.create(DIRS$tables, recursive = TRUE, showWarnings = FALSE)
dir.create(DIRS$derived_data, recursive = TRUE, showWarnings = FALSE)
dir.create(DIRS$logs, recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# Helpers
# ----------------------------
normalize_name <- function(x) {
  x <- trimws(x)
  x <- gsub("\\.+", "_", x)
  x <- gsub("\\s+", "_", x)
  x <- gsub("[^A-Za-z0-9_]", "", x)
  tolower(x)
}

safe_profile_tab <- function(x) {
  x <- x[!is.na(x) & trimws(x) != ""]
  if (length(x) == 0L) return(integer(0))
  table(x)
}

icc_like_one <- function(score, group) {
  ok <- is.finite(score) & !is.na(score) & !is.na(group) & trimws(group) != ""
  score <- score[ok]
  group <- as.character(group[ok])

  if (length(score) < MIN_OBS_FOR_ICC || length(unique(group)) < MIN_GROUPS_FOR_ICC) {
    return(c(
      n_obs = length(score),
      n_groups = length(unique(group)),
      grand_mean = ifelse(length(score) == 0L, NA_real_, mean(score)),
      between_var = NA_real_,
      within_var = NA_real_,
      icc_like = NA_real_
    ))
  }

  grand_mean <- mean(score)
  grp <- split(score, group)
  group_means <- vapply(grp, mean, numeric(1))
  group_ns <- vapply(grp, length, integer(1))
  between_var <- weighted.mean((group_means - grand_mean)^2, w = group_ns)

  within_var <- mean(vapply(grp, function(z) {
    if (length(z) <= 1L) return(0)
    var(z)
  }, numeric(1)), na.rm = TRUE)

  icc_like <- between_var / (between_var + within_var)

  c(
    n_obs = length(score),
    n_groups = length(unique(group)),
    grand_mean = grand_mean,
    between_var = between_var,
    within_var = within_var,
    icc_like = icc_like
  )
}

is_believable_system_field <- function(df, colname) {
  nmn <- normalize_name(colname)

  # hard exclusions: likely person / contact / title fields
  bad_patterns <- c(
    "contact", "primary_contact", "name$", "_name$", "firstname", "lastname",
    "ceo", "president", "officer", "director", "chair", "executive",
    "mr_", "mrs_", "ms_", "dr_", "email", "phone", "fax", "address"
  )
  if (any(grepl(paste(bad_patterns, collapse = "|"), nmn))) return(FALSE)

  x <- df[[colname]]
  x_chr <- as.character(x)
  x_chr[is.na(x_chr) | trimws(x_chr) == ""] <- NA_character_
  ux <- unique(x_chr[!is.na(x_chr)])
  if (length(ux) < 2L) return(FALSE)

  # reject fields that look like mostly personal names
  sampled <- head(ux, 100)
  has_title <- mean(grepl("^(mr|mrs|ms|dr|prof)\\b", tolower(sampled))) > 0.20
  has_person_like <- mean(grepl("^[A-Za-z]+\\.?\\s+[A-Za-z]", sampled)) > 0.50
  if (has_title || has_person_like) return(FALSE)

  TRUE
}

score_system_candidate <- function(df, colname) {
  nmn <- normalize_name(colname)
  score <- 0

  # positive signals
  if (grepl("system_id|sys_id|network_id|mhs_id|multihospital", nmn)) score <- score + 10
  if (grepl("system_member|member_system|sys_member|network_member", nmn)) score <- score + 8
  if (grepl("^mhs$|_mhs$|mhs_", nmn)) score <- score + 7
  if (grepl("system|network|sys", nmn)) score <- score + 5

  # negative signals
  if (grepl("contact|name|ceo|president|officer|director|executive", nmn)) score <- score - 20
  if (grepl("email|phone|fax|address", nmn)) score <- score - 10

  # prefer moderate-cardinality system variables
  x <- as.character(df[[colname]])
  x[is.na(x) | trimws(x) == ""] <- NA_character_
  n_unique <- length(unique(x[!is.na(x)]))
  score <- score + ifelse(n_unique >= 2 & n_unique <= 2000, 2, 0)

  score
}

# ----------------------------
# Detect believable system field
# ----------------------------
nm <- names(left_df)
system_like <- nm[grepl("system|network|sys|mhs|multihospital|member", normalize_name(nm))]
system_like <- unique(system_like)

system_like <- system_like[vapply(system_like, function(v) is_believable_system_field(left_df, v), logical(1))]

system_var <- character(0)
system_detection_note <- "No believable system variable detected; system nesting skipped."

if (length(system_like) > 0L) {
  scores <- vapply(system_like, function(v) score_system_candidate(left_df, v), numeric(1))
  system_like <- system_like[order(scores, decreasing = TRUE)]
  chosen <- system_like[1]

  # sanity guard: require at least two non-missing distinct groups
  x <- as.character(left_df[[chosen]])
  x[is.na(x) | trimws(x) == ""] <- NA_character_
  if (length(unique(x[!is.na(x)])) >= 2L) {
    system_var <- chosen
    system_detection_note <- paste0("Selected believable system variable: ", system_var)
  }
}

# ----------------------------
# Merge grouping variable
# ----------------------------
if (length(system_var) == 0L) {
  system_df <- data.frame(ID = left_df$ID, .system_group_raw = NA_character_, stringsAsFactors = FALSE)
} else {
  system_df <- left_df[, c("ID", system_var), drop = FALSE]
  names(system_df) <- c("ID", ".system_group_raw")
}

core_df <- merge(core_df, system_df, by = "ID", all.x = TRUE, sort = FALSE)

sys <- as.character(core_df$.system_group_raw)
sys[is.na(sys) | trimws(sys) == ""] <- "MissingSystem"
core_df$.system_group <- sys

assigned_df <- core_df[!is.na(core_df$profile) & trimws(core_df$profile) != "", , drop = FALSE]
unassigned_df <- core_df[is.na(core_df$profile) | trimws(core_df$profile) == "", , drop = FALSE]

# ----------------------------
# System nesting summary: assigned
# ----------------------------
if (length(system_var) == 0L) {
  sys_tbl <- data.frame(
    note = system_detection_note,
    stringsAsFactors = FALSE
  )
} else if (nrow(assigned_df) > 0L) {
  assigned_split <- split(assigned_df, assigned_df$.system_group)

  sys_list <- lapply(assigned_split, function(d) {
    tab <- safe_profile_tab(d$profile)
    if (length(tab) == 0L) return(NULL)

    data.frame(
      system_group = unique(d$.system_group)[1],
      n_hospitals = nrow(d),
      n_profiles_present = sum(tab > 0),
      dominant_profile = names(tab)[which.max(tab)],
      dominant_profile_share = max(tab) / sum(tab),
      stringsAsFactors = FALSE
    )
  })
  sys_list <- sys_list[!vapply(sys_list, is.null, logical(1))]

  if (length(sys_list) == 0L) {
    sys_tbl <- data.frame(
      note = "No assigned hospitals with valid profile counts found for selected system variable.",
      stringsAsFactors = FALSE
    )
  } else {
    sys_tbl <- do.call(rbind, sys_list)
    sys_tbl <- sys_tbl[order(-sys_tbl$n_hospitals, sys_tbl$system_group), , drop = FALSE]
  }
} else {
  sys_tbl <- data.frame(
    note = "No assigned hospitals found in profile assignments.",
    stringsAsFactors = FALSE
  )
}

write.csv(
  sys_tbl,
  file.path(DIRS$tables, "system_nesting_summary_masked_2023_v7.csv"),
  row.names = FALSE
)

# ----------------------------
# System nesting summary: unassigned
# ----------------------------
if (length(system_var) == 0L) {
  unassigned_tbl <- data.frame(
    note = system_detection_note,
    stringsAsFactors = FALSE
  )
} else if (nrow(unassigned_df) > 0L) {
  unassigned_tbl <- aggregate(
    x = list(n_unassigned = rep(1L, nrow(unassigned_df))),
    by = list(system_group = unassigned_df$.system_group),
    FUN = sum
  )
  unassigned_tbl <- unassigned_tbl[order(-unassigned_tbl$n_unassigned, unassigned_tbl$system_group), , drop = FALSE]
} else {
  unassigned_tbl <- data.frame(
    note = "No unassigned hospitals.",
    stringsAsFactors = FALSE
  )
}

write.csv(
  unassigned_tbl,
  file.path(DIRS$tables, "system_nesting_unassigned_summary_masked_2023_v7.csv"),
  row.names = FALSE
)

# ----------------------------
# ICC-like module score nesting
# ----------------------------
score_df <- clust_bundle$module_score_df
score_cols <- grep("_score$", names(score_df), value = TRUE)

if (length(system_var) == 0L) {
  icc_tbl <- data.frame(
    note = system_detection_note,
    stringsAsFactors = FALSE
  )
} else if (!is.null(score_df) && length(score_cols) > 0L) {
  score_df2 <- merge(score_df, system_df, by = "ID", all.x = TRUE, sort = FALSE)
  sys2 <- as.character(score_df2$.system_group_raw)
  sys2[is.na(sys2) | trimws(sys2) == ""] <- "MissingSystem"
  score_df2$.system_group <- sys2

  icc_tbl <- do.call(rbind, lapply(score_cols, function(sc) {
    out <- icc_like_one(score_df2[[sc]], score_df2$.system_group)
    data.frame(module_score = sc, t(out), stringsAsFactors = FALSE)
  }))
  rownames(icc_tbl) <- NULL
} else {
  icc_tbl <- data.frame(
    note = "No module score columns available for ICC-like summary.",
    stringsAsFactors = FALSE
  )
}

write.csv(
  icc_tbl,
  file.path(DIRS$tables, "system_icc_module_scores_masked_2023_v7.csv"),
  row.names = FALSE
)

# ----------------------------
# Optional outcome validation
# ----------------------------
outcome_candidates <- c(
  file.path(DIRS$derived_data, "external_outcomes_2023.csv"),
  file.path(DIRS$tables, "external_outcomes_2023.csv"),
  file.path(getwd(), "external_outcomes_2023.csv")
)

outcome_path <- outcome_candidates[file.exists(outcome_candidates)]
have_real_outcomes <- length(outcome_path) > 0L

outcome_summary <- data.frame(
  note = "No real external outcomes file detected. Validation not run.",
  stringsAsFactors = FALSE
)
outcome_incremental <- data.frame(
  note = "No real external outcomes file detected. Incremental validation not run.",
  stringsAsFactors = FALSE
)

if (have_real_outcomes) {
  outcome_df <- read.csv(outcome_path[1], stringsAsFactors = FALSE, check.names = FALSE)

  if ("ID" %in% names(outcome_df)) {
    valid_df <- merge(
      assigned_df[, c("ID", "profile"), drop = FALSE],
      outcome_df, by = "ID", all = FALSE, sort = FALSE
    )

    outcome_vars <- setdiff(names(valid_df), c("ID", "profile"))
    numeric_outcomes <- outcome_vars[vapply(valid_df[outcome_vars], is.numeric, logical(1))]

    if (length(numeric_outcomes) > 0L) {
      outcome_summary <- do.call(rbind, lapply(numeric_outcomes, function(yy) {
        d <- valid_df[is.finite(valid_df[[yy]]) & !is.na(valid_df[[yy]]), c("profile", yy), drop = FALSE]
        if (nrow(d) < 20L || length(unique(d$profile)) < 2L) {
          return(data.frame(outcome = yy, n = nrow(d), p_value = NA_real_, stringsAsFactors = FALSE))
        }
        fit <- try(stats::lm(stats::as.formula(paste(yy, "~ profile")), data = d), silent = TRUE)
        if (inherits(fit, "try-error")) {
          return(data.frame(outcome = yy, n = nrow(d), p_value = NA_real_, stringsAsFactors = FALSE))
        }
        a <- anova(fit)
        pv <- if ("profile" %in% rownames(a)) a["profile", "Pr(>F)"] else NA_real_
        data.frame(outcome = yy, n = nrow(d), p_value = pv, stringsAsFactors = FALSE)
      }))
      outcome_incremental <- outcome_summary
      names(outcome_incremental)[names(outcome_incremental) == "p_value"] <- "profile_model_p_value"
    } else {
      outcome_summary <- data.frame(note = "Outcomes file found, but no numeric outcomes were available.", stringsAsFactors = FALSE)
      outcome_incremental <- outcome_summary
    }
  } else {
    outcome_summary <- data.frame(note = "Outcomes file found, but no ID column was present.", stringsAsFactors = FALSE)
    outcome_incremental <- outcome_summary
  }
} else {
  template_path <- file.path(DIRS$tables, OUTCOME_TEMPLATE_NAME)
  if (!file.exists(template_path)) {
    template <- data.frame(
      ID = character(0),
      outcome_1 = numeric(0),
      outcome_2 = numeric(0),
      stringsAsFactors = FALSE
    )
    write.csv(template, template_path, row.names = FALSE)
  }
}

write.csv(
  outcome_summary,
  file.path(DIRS$tables, "outcome_validation_summary_masked_2023_v7.csv"),
  row.names = FALSE
)

write.csv(
  outcome_incremental,
  file.path(DIRS$tables, "outcome_incremental_value_summary_masked_2023_v7.csv"),
  row.names = FALSE
)

# ----------------------------
# Save summary bundle
# ----------------------------
eval_bundle <- list(
  system_var = ifelse(length(system_var) == 0L, NA_character_, system_var),
  system_detection_note = system_detection_note,
  assigned_n = nrow(assigned_df),
  unassigned_n = nrow(unassigned_df),
  system_nesting_summary = sys_tbl,
  system_unassigned_summary = unassigned_tbl,
  system_icc_module_scores = icc_tbl,
  outcome_validation_summary = outcome_summary,
  outcome_incremental_value_summary = outcome_incremental
)

saveRDS(
  eval_bundle,
  file.path(DIRS$derived_data, "aha_2023_masked_validation_bundle_v7.rds")
)

# ----------------------------
# Log
# ----------------------------
log_file <- file.path(
  DIRS$logs,
  paste0("16_rebuild_system_nesting_and_validation_masked_v7_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt")
)

sink(log_file)
cat("16_rebuild_system_nesting_and_validation_masked_v7.R log\n")
cat("=======================================\n")
cat("Run time: ", as.character(Sys.time()), "\n\n")
cat("System detection note: ", system_detection_note, "\n")
cat("Assigned hospitals: ", nrow(assigned_df), "\n")
cat("Unassigned hospitals: ", nrow(unassigned_df), "\n")
cat("Real outcomes detected: ", have_real_outcomes, "\n\n")
cat("Top system nesting summary rows:\n")
print(utils::head(sys_tbl, 20))
cat("\nTop unassigned system rows:\n")
print(utils::head(unassigned_tbl, 20))
cat("\nICC-like module score summary:\n")
print(icc_tbl)
cat("\nOutcome validation summary:\n")
print(outcome_summary)
cat("\nSession info:\n")
print(sessionInfo())
sink()

cat("Script 16 v7 complete.\n")
