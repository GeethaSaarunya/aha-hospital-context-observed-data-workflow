# ============================================================
# 17_update_manuscript_tables_and_figures_v3.R
#
# PURPOSE
# -------
# Journal-polished reporting/update script for the remediated pipeline.
#
# This version keeps the exact same analytic content as v2, but uses a
# more restrained, publication-style visual design:
#   - smaller titles and subtitles
#   - quieter grid lines
#   - panel tags (A-D)
#   - tighter margins
#   - color-blind-friendly Okabe-Ito palette with lower visual noise
#   - figure captions saved separately
#
# IMPORTANT
# ---------
# This script does NOT run new analyses. It only reads the final outputs
# from Scripts 13 v4, 14 v4, 15 v5, and 16 v7 and updates manuscript-
# ready numbers, tables, and figures.
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999)

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
  library(gridExtra)
})

# ------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------
rescale01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (!all(is.finite(rng)) || diff(rng) == 0) {
    return(rep(0.5, length(x)))
  }
  (x - rng[1]) / diff(rng)
}

fmt_int <- function(x) format(round(x), big.mark = ",", trim = TRUE, scientific = FALSE)
fmt_num <- function(x, digits = 2) format(round(x, digits), nsmall = digits, trim = TRUE, scientific = FALSE)

# Okabe-Ito palette
okabe <- c(
  orange = "#E69F00",
  sky = "#56B4E9",
  bluishgreen = "#009E73",
  yellow = "#F0E442",
  blue = "#0072B2",
  vermillion = "#D55E00",
  reddishpurple = "#CC79A7",
  black = "#000000",
  gray = "#999999",
  lightgray = "#E6E6E6",
  darkgray = "#4D4D4D"
)

panel_tag_theme <- function() {
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "#DDDDDD", linewidth = 0.35),
    axis.title = element_text(face = "bold", size = 10),
    axis.text = element_text(size = 9, color = okabe["darkgray"]),
    axis.text.x = element_text(size = 9),
    plot.title = element_text(face = "bold", size = 12, hjust = 0),
    plot.subtitle = element_text(size = 9.5, hjust = 0, color = okabe["darkgray"]),
    plot.tag = element_text(face = "bold", size = 13),
    plot.margin = margin(6, 8, 6, 8),
    legend.title = element_text(face = "bold", size = 9.5),
    legend.text = element_text(size = 8.5)
  )
}

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------
setwd("~/Desktop/AHA")
CONFIG_PATH <- file.path(getwd(), "metadata", "project_config.rds")
if (!file.exists(CONFIG_PATH)) stop("Run 00_setup_paths.R first.")

cfg <- readRDS(CONFIG_PATH)
DIRS <- cfg$DIRS

dir.create(DIRS$tables, recursive = TRUE, showWarnings = FALSE)
dir.create(DIRS$figures, recursive = TRUE, showWarnings = FALSE)
dir.create(DIRS$derived_data, recursive = TRUE, showWarnings = FALSE)
dir.create(DIRS$logs, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Inputs
# ------------------------------------------------------------
bundle13_path <- file.path(DIRS$derived_data, "aha_2023_masked_weight_bundle_v4.rds")
bundle14_path <- file.path(DIRS$derived_data, "aha_2023_masked_smc_bundle_v4.rds")
bundle15_path <- file.path(DIRS$derived_data, "aha_2023_masked_clustering_bundle_v5.rds")
bundle16_path <- file.path(DIRS$derived_data, "aha_2023_masked_validation_bundle_v7.rds")

required <- c(bundle13_path, bundle14_path, bundle15_path, bundle16_path)
missing_files <- required[!file.exists(required)]
if (length(missing_files) > 0) stop("Missing required bundles:\n", paste(missing_files, collapse = "\n"))

b13 <- readRDS(bundle13_path)
b14 <- readRDS(bundle14_path)
b15 <- readRDS(bundle15_path)
b16 <- readRDS(bundle16_path)

# ------------------------------------------------------------
# Extract key numbers
# ------------------------------------------------------------
full_frame_n <- if (!is.null(b13$full_frame_n)) b13$full_frame_n else nrow(b13$work_df)
selected_vars_n <- length(b13$selected_vars)
selected_it_n <- length(b13$selected_it)
stable_vars_n <- length(b14$stable_vars)

module_final_tab <- sort(table(b15$module_assignments_usable$module_final), decreasing = TRUE)
usable_modules <- names(module_final_tab)[names(module_final_tab) != "M0_residual"]
usable_module_count <- length(usable_modules)
score_cols <- b15$kept_score_cols
score_cols_n <- length(score_cols)

eligible_n <- if (!is.null(b15$eligible_df)) nrow(b15$eligible_df) else sum(b15$module_score_df$cluster_eligible, na.rm = TRUE)
assigned_n <- if (!is.null(b15$profile_assignments)) sum(!is.na(b15$profile_assignments$profile) & trimws(b15$profile_assignments$profile) != "") else 0L
unassigned_n <- full_frame_n - assigned_n

fit_tbl <- b15$fit_table
best_k <- b15$best_k
profile_sizes <- b15$profile_sizes
if (is.null(profile_sizes) || nrow(profile_sizes) == 0L) stop("No profile sizes found in v5 clustering bundle.")
profile_sizes <- profile_sizes[order(profile_sizes$profile), , drop = FALSE]
assigned_tab <- setNames(profile_sizes$n_hospitals, profile_sizes$profile)

system_note <- if (!is.null(b16$system_detection_note)) b16$system_detection_note else NA_character_
system_var <- if (!is.null(b16$system_var)) b16$system_var else NA_character_

# Profile means
profile_assign <- b15$profile_assignments
score_cols_all <- grep("_score$", names(profile_assign), value = TRUE)

profile_mean_rows <- list()
for (pp in sort(unique(na.omit(profile_assign$profile)))) {
  d <- profile_assign[profile_assign$profile == pp, , drop = FALSE]
  out <- data.frame(profile = pp, stringsAsFactors = FALSE)
  for (sc in score_cols_all) {
    vals <- d[[sc]]
    out[[sc]] <- if (all(is.na(vals))) NA_real_ else mean(vals, na.rm = TRUE)
  }
  profile_mean_rows[[pp]] <- out
}
profile_means <- do.call(rbind, profile_mean_rows)

# ------------------------------------------------------------
# Save manuscript numbers and table
# ------------------------------------------------------------
numbers_df <- data.frame(
  metric = c(
    "full_frame_n",
    "selected_vars_n",
    "selected_it_n",
    "stable_vars_n",
    "usable_module_count",
    "score_cols_n",
    "eligible_n",
    "assigned_n",
    "unassigned_n",
    "best_k"
  ),
  value = c(
    full_frame_n, selected_vars_n, selected_it_n, stable_vars_n,
    usable_module_count, score_cols_n, eligible_n, assigned_n,
    unassigned_n, best_k
  ),
  stringsAsFactors = FALSE
)

profile_num_df <- data.frame(
  metric = paste0("profile_size_", profile_sizes$profile),
  value = profile_sizes$n_hospitals,
  stringsAsFactors = FALSE
)

numbers_out <- rbind(numbers_df, profile_num_df)
write.csv(numbers_out, file.path(DIRS$tables, "manuscript_numbers_v17_v3.csv"), row.names = FALSE)

table1 <- data.frame(
  item = c(
    "Full linked hospital frame",
    "Selected variables retained",
    "Selected IT variables retained",
    "Stable variables retained after masked filtering",
    "Usable modules after tiny-module absorption",
    "Usable module-score dimensions for clustering",
    "Cluster-eligible hospitals",
    "Assigned hospitals",
    "Unassigned hospitals",
    "Selected profile solution (K)",
    "Profile P1 size",
    "Profile P2 size",
    "Profile P3 size"
  ),
  value = c(
    full_frame_n, selected_vars_n, selected_it_n, stable_vars_n,
    usable_module_count, score_cols_n, eligible_n, assigned_n,
    unassigned_n, best_k,
    ifelse("P1" %in% names(assigned_tab), assigned_tab["P1"], NA),
    ifelse("P2" %in% names(assigned_tab), assigned_tab["P2"], NA),
    ifelse("P3" %in% names(assigned_tab), assigned_tab["P3"], NA)
  ),
  interpretation = c(
    "Weighted observed-data analysis frame",
    "Variables retained from ASDB + IT selection",
    "IT subset used to define response/masking",
    "Variables surviving observedness and binary-support filters",
    "Modules remaining after tiny-module absorption",
    "Final dimensions used for hospital clustering",
    "Hospitals meeting the module-score eligibility rule",
    "Hospitals receiving profile assignment",
    "Hospitals not assigned under the current eligibility rule",
    "Practical clustering solution selected from the fit table",
    "Largest coarse profile",
    "Intermediate coarse profile",
    "Smaller high-deviation tail profile"
  ),
  stringsAsFactors = FALSE
)
write.csv(table1, file.path(DIRS$tables, "table1_manuscript_summary_v17_v3.csv"), row.names = FALSE)

# ------------------------------------------------------------
# Figure 1 data prep
# ------------------------------------------------------------
hosp_df <- data.frame(
  stage = factor(
    c("Full linked frame", "Cluster-eligible hospitals", "Assigned hospitals"),
    levels = c("Full linked frame", "Cluster-eligible hospitals", "Assigned hospitals")
  ),
  n = c(full_frame_n, eligible_n, assigned_n),
  fill = c(okabe["blue"], okabe["sky"], okabe["bluishgreen"]),
  stringsAsFactors = FALSE
)

var_df <- data.frame(
  stage = factor(
    c("Selected variables", "Selected IT variables", "Stable variables", "Usable modules", "Score dimensions"),
    levels = c("Selected variables", "Selected IT variables", "Stable variables", "Usable modules", "Score dimensions")
  ),
  n = c(selected_vars_n, selected_it_n, stable_vars_n, usable_module_count, score_cols_n),
  fill = c(okabe["orange"], okabe["yellow"], okabe["vermillion"], okabe["reddishpurple"], okabe["black"]),
  stringsAsFactors = FALSE
)

module_plot_df <- data.frame(
  module = factor(names(module_final_tab), levels = names(module_final_tab)),
  n_vars = as.numeric(module_final_tab),
  stringsAsFactors = FALSE
)
module_fill_map <- c(M1 = okabe["blue"], M2 = okabe["orange"], M0_residual = okabe["gray"])
module_plot_df$fill <- unname(module_fill_map[as.character(module_plot_df$module)])
module_plot_df$fill[is.na(module_plot_df$fill)] <- okabe["gray"]

profile_plot_df <- profile_sizes
profile_plot_df$profile <- factor(profile_plot_df$profile, levels = profile_plot_df$profile)
profile_fill_map <- c(P1 = okabe["blue"], P2 = okabe["orange"], P3 = okabe["bluishgreen"], P4 = okabe["reddishpurple"], P5 = okabe["vermillion"])
profile_plot_df$fill <- unname(profile_fill_map[as.character(profile_plot_df$profile)])
profile_plot_df$fill[is.na(profile_plot_df$fill)] <- okabe["gray"]

base_theme <- theme_minimal(base_size = 11) + panel_tag_theme()

p_hosp <- ggplot(hosp_df, aes(x = stage, y = n, fill = stage)) +
  geom_col(width = 0.62, show.legend = FALSE) +
  geom_text(aes(label = fmt_int(n)), vjust = -0.25, size = 3.4, color = okabe["black"]) +
  scale_fill_manual(values = setNames(hosp_df$fill, hosp_df$stage)) +
  labs(
    tag = "A",
    title = "Hospital analysis flow",
    subtitle = "From full linked frame to assigned hospitals",
    x = NULL, y = "Hospitals"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.10))) +
  base_theme +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

p_vars <- ggplot(var_df, aes(x = stage, y = n, fill = stage)) +
  geom_col(width = 0.62, show.legend = FALSE) +
  geom_text(aes(label = fmt_int(n)), vjust = -0.25, size = 3.4, color = okabe["black"]) +
  scale_fill_manual(values = setNames(var_df$fill, var_df$stage)) +
  labs(
    tag = "B",
    title = "Variable filtering summary",
    subtitle = "From selected inputs to final score dimensions",
    x = NULL, y = "Count"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.10))) +
  base_theme +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

p_modules <- ggplot(module_plot_df, aes(x = module, y = n_vars, fill = module)) +
  geom_col(width = 0.62, show.legend = FALSE) +
  geom_text(aes(label = fmt_int(n_vars)), vjust = -0.25, size = 3.4, color = okabe["black"]) +
  scale_fill_manual(values = setNames(module_plot_df$fill, module_plot_df$module)) +
  labs(
    tag = "C",
    title = "Final module sizes",
    subtitle = "After tiny-module absorption",
    x = NULL, y = "Variables"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.10))) +
  base_theme

p_profiles <- ggplot(profile_plot_df, aes(x = profile, y = n_hospitals, fill = profile)) +
  geom_col(width = 0.62, show.legend = FALSE) +
  geom_text(aes(label = fmt_int(n_hospitals)), vjust = -0.25, size = 3.4, color = okabe["black"]) +
  scale_fill_manual(values = setNames(profile_plot_df$fill, profile_plot_df$profile)) +
  labs(
    tag = "D",
    title = "Hospital profile sizes",
    subtitle = paste0("Selected clustering solution: K = ", best_k),
    x = NULL, y = "Hospitals"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.10))) +
  base_theme

fig1_path_png <- file.path(DIRS$figures, "figure1_pipeline_modules_profiles_v17_v3.png")
fig1_path_pdf <- file.path(DIRS$figures, "figure1_pipeline_modules_profiles_v17_v3.pdf")

png(fig1_path_png, width = 2000, height = 1400, res = 180, bg = "white")
grid.arrange(p_hosp, p_vars, p_modules, p_profiles, ncol = 2)
dev.off()

pdf(fig1_path_pdf, width = 12.5, height = 8.8, bg = "white")
grid.arrange(p_hosp, p_vars, p_modules, p_profiles, ncol = 2)
dev.off()

# ------------------------------------------------------------
# Figure 2: journal-polished heatmap
# ------------------------------------------------------------
fig2_path_png <- NA_character_
fig2_path_pdf <- NA_character_

if (!is.null(profile_means) && nrow(profile_means) > 0 && length(score_cols_all) > 0) {
  heat_list <- list()
  for (sc in score_cols_all) {
    heat_list[[sc]] <- data.frame(
      profile = profile_means$profile,
      module_score = sc,
      raw_value = profile_means[[sc]],
      stringsAsFactors = FALSE
    )
  }
  heat_df <- do.call(rbind, heat_list)
  heat_df$module_score_clean <- gsub("_score$", "", heat_df$module_score)
  heat_df$module_score_clean <- factor(heat_df$module_score_clean, levels = unique(gsub("_score$", "", score_cols_all)))
  heat_df$profile <- factor(heat_df$profile, levels = sort(unique(profile_means$profile)))

  # Rescale fill within each score dimension for visual comparability
  heat_df$fill_scaled <- ave(heat_df$raw_value, heat_df$module_score, FUN = rescale01)
  heat_df$label <- ifelse(is.na(heat_df$raw_value), "NA", fmt_num(heat_df$raw_value, 2))
  heat_df$label_col <- ifelse(is.na(heat_df$fill_scaled), okabe["black"], ifelse(heat_df$fill_scaled >= 0.60, "white", okabe["black"]))

  p_heat <- ggplot(heat_df, aes(x = module_score_clean, y = profile, fill = fill_scaled)) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(aes(label = label, color = label_col), size = 3.8, show.legend = FALSE) +
    scale_color_identity() +
    scale_fill_gradientn(
      colours = c("#F7F7F7", okabe["sky"], okabe["blue"]),
      limits = c(0, 1),
      breaks = c(0, 0.5, 1),
      labels = c("Low", "Mid", "High"),
      name = "Relative level\nwithin score"
    ) +
    labs(
      title = "Profile mean module scores",
      subtitle = "Tile shading is rescaled within each score dimension; labels show raw means",
      x = "Module score dimension", y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.title = element_text(face = "bold", size = 10),
      axis.text = element_text(size = 9.5, color = okabe["darkgray"]),
      plot.title = element_text(face = "bold", size = 12, hjust = 0),
      plot.subtitle = element_text(size = 9.5, hjust = 0, color = okabe["darkgray"]),
      legend.title = element_text(face = "bold", size = 9.5),
      legend.text = element_text(size = 8.5),
      plot.margin = margin(6, 8, 6, 8)
    )

  fig2_path_png <- file.path(DIRS$figures, "figure2_profile_module_heatmap_v17_v3.png")
  fig2_path_pdf <- file.path(DIRS$figures, "figure2_profile_module_heatmap_v17_v3.pdf")

  ggsave(fig2_path_png, p_heat, width = 6.8, height = 4.6, dpi = 200, bg = "white")
  ggsave(fig2_path_pdf, p_heat, width = 6.8, height = 4.6, bg = "white")
}

# ------------------------------------------------------------
# Captions
# ------------------------------------------------------------
captions <- data.frame(
  figure = c("Figure 1", "Figure 2"),
  caption = c(
    "Summary of the remediated weighted observed-data pipeline. Panel A shows the hospital analysis flow from the full linked frame to assigned hospitals. Panel B shows the variable-filtering sequence from selected inputs to final score dimensions. Panel C shows final module sizes after tiny-module absorption. Panel D shows hospital profile sizes under the selected K = 3 solution.",
    "Heatmap of mean module scores by assigned hospital profile. Tile shading is rescaled within each score dimension to improve visual comparability, while overlaid labels display the raw mean scores. The final clustering solution should be interpreted as a coarse descriptive grouping rather than a detailed national taxonomy."
  ),
  stringsAsFactors = FALSE
)
write.csv(captions, file.path(DIRS$tables, "figure_captions_v17_v3.csv"), row.names = FALSE)

# ------------------------------------------------------------
# Log
# ------------------------------------------------------------
log_file <- file.path(
  DIRS$logs,
  paste0("17_update_manuscript_tables_and_figures_v3_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt")
)

sink(log_file)
cat("17_update_manuscript_tables_and_figures_v3.R log\n")
cat("================================================\n")
cat("Run time: ", as.character(Sys.time()), "\n\n")
cat("Full frame N: ", full_frame_n, "\n")
cat("Selected vars retained: ", selected_vars_n, "\n")
cat("Selected IT vars retained: ", selected_it_n, "\n")
cat("Stable vars retained: ", stable_vars_n, "\n")
cat("Usable modules: ", usable_module_count, "\n")
cat("Score dimensions: ", score_cols_n, "\n")
cat("Eligible hospitals: ", eligible_n, "\n")
cat("Assigned hospitals: ", assigned_n, "\n")
cat("Unassigned hospitals: ", unassigned_n, "\n")
cat("Best K: ", best_k, "\n\n")
cat("Profile sizes:\n")
print(profile_sizes)
cat("\nSystem note:\n")
print(system_note)
cat("\nSelected system variable:\n")
print(system_var)
cat("\nSession info:\n")
print(sessionInfo())
sink()

cat("Script 17 v3 complete.\n")
