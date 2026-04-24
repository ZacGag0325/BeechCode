#!/usr/bin/env Rscript

# ============================================================================
# nearest_neighbor_sampling_check.R
# ----------------------------------------------------------------------------
# Standalone script to evaluate nearest-neighbor distances among sampled beech
# stems within each site.
#
# What this script does:
#   1) Reads stem-level data from a user-specified file.
#   2) Resolves required columns (site ID, sample ID, x, y).
#   3) Calculates nearest sampled neighbor distances within each site.
#   4) Writes per-sample and per-site summary tables to disk.
#   5) Creates and saves two publication-style figures.
#
# NOTE:
#   - This script is standalone and does not modify any project pipeline.
#   - Edit the USER SETTINGS section before running.
# ============================================================================

# ----------------------------- USER SETTINGS ---------------------------------
# Set your input file path here (CSV/TSV supported by extension).
input_file <- "data/stem_level_data.csv"

# Optional output directory (will be created if needed).
output_dir <- "outputs/nearest_neighbor_sampling_check"

# REQUIRED COLUMN MAPPINGS:
# - If you know the exact column names, set them as strings.
# - If left as NULL, the script will attempt to auto-detect from common names.
site_col   <- NULL  # e.g., "site", "site_id"
sample_col <- NULL  # e.g., "sample_id", "individual_id", "stem_id"
x_col      <- NULL  # e.g., "x", "x_m", "utm_x"
y_col      <- NULL  # e.g., "y", "y_m", "utm_y"

# Distance thresholds (meters) used in site-level summary proportions.
thresholds_m <- c(8, 10, 12)
# ----------------------------------------------------------------------------

# ----------------------------- PACKAGE CHECK ---------------------------------
required_pkgs <- c("dplyr", "readr", "stringr", "ggplot2", "tidyr", "purrr")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_pkgs) > 0) {
  stop(
    paste0(
      "Missing required packages: ", paste(missing_pkgs, collapse = ", "),
      "\nPlease install them before running, e.g.:\n",
      "install.packages(c(",
      paste(sprintf('"%s"', missing_pkgs), collapse = ", "),
      "))"
    ),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(tidyr)
  library(purrr)
})
# ----------------------------------------------------------------------------

# ----------------------------- HELPER FUNCTIONS ------------------------------
read_input_data <- function(path) {
  if (!file.exists(path)) {
    xlsx_candidates <- list.files(
      path = ".",
      pattern = "\\.(xlsx|xls)$",
      recursive = TRUE,
      full.names = TRUE,
      ignore.case = TRUE
    )
    
    stop(
      paste0(
        "Input file not found: ", path,
        "\nPlease update 'input_file' in the USER SETTINGS section.",
        "\n\nIf you intended to run the adjacent-number Excel analysis, use:",
        "\n  source('scripts/nearest_neighbour_adjacent_numbers_analysis.R')",
        "\n(or the US spelling alias if present).",
        if (length(xlsx_candidates) > 0) {
          paste0(
            "\n\nExcel files detected in this project (examples):\n  - ",
            paste(utils::head(xlsx_candidates, 10), collapse = "\n  - ")
          )
        } else {
          "\n\nNo Excel files were detected in the current project tree."
        }
      ),
      call. = FALSE
    )
  }
  
  ext <- tolower(tools::file_ext(path))
  
  if (ext %in% c("csv")) {
    readr::read_csv(path, show_col_types = FALSE)
  } else if (ext %in% c("tsv", "txt")) {
    readr::read_tsv(path, show_col_types = FALSE)
  } else {
    stop(
      paste0(
        "Unsupported file extension '.", ext, "'.\n",
        "Please provide a .csv, .tsv, or .txt file."
      ),
      call. = FALSE
    )
  }
}

choose_col <- function(data, user_col, candidate_names, col_label) {
  nms <- names(data)
  
  # If user explicitly provides a column name, validate it.
  if (!is.null(user_col)) {
    if (!user_col %in% nms) {
      stop(
        paste0(
          "Configured ", col_label, " column '", user_col, "' not found in data.\n",
          "Available columns: ", paste(nms, collapse = ", ")
        ),
        call. = FALSE
      )
    }
    return(user_col)
  }
  
  # Otherwise, attempt auto-detection (case-insensitive exact match first).
  lower_nms <- tolower(nms)
  lower_candidates <- tolower(candidate_names)
  
  exact_hits_idx <- which(lower_nms %in% lower_candidates)
  if (length(exact_hits_idx) >= 1) {
    return(nms[exact_hits_idx[1]])
  }
  
  # Then attempt partial match as fallback.
  partial_hits_idx <- unique(unlist(lapply(lower_candidates, function(cn) which(str_detect(lower_nms, fixed(cn))))))
  if (length(partial_hits_idx) >= 1) {
    return(nms[partial_hits_idx[1]])
  }
  
  stop(
    paste0(
      "Could not auto-detect ", col_label, " column.\n",
      "Please set '", gsub(" ", "_", col_label), "' in USER SETTINGS.\n",
      "Available columns: ", paste(nms, collapse = ", ")
    ),
    call. = FALSE
  )
}

compute_nn_within_site <- function(df_site) {
  # Expected columns: site, sample_id, x, y
  n <- nrow(df_site)
  
  if (n < 2) {
    # Graceful handling for sites with insufficient points.
    return(
      df_site %>%
        mutate(
          nearest_neighbor_id = NA_character_,
          nearest_neighbor_distance_m = NA_real_
        )
    )
  }
  
  coords <- as.matrix(df_site[, c("x", "y")])
  
  # Euclidean distance matrix.
  dmat <- as.matrix(stats::dist(coords, method = "euclidean"))
  
  # Exclude self-distance from nearest-neighbor search.
  diag(dmat) <- Inf
  
  nn_index <- max.col(-dmat, ties.method = "first")
  nn_distance <- dmat[cbind(seq_len(n), nn_index)]
  
  df_site %>%
    mutate(
      nearest_neighbor_id = sample_id[nn_index],
      nearest_neighbor_distance_m = nn_distance
    )
}
# ----------------------------------------------------------------------------

# ----------------------------- MAIN WORKFLOW ---------------------------------
# 1) Read data.
stems_raw <- read_input_data(input_file)

# 2) Resolve required columns.
resolved_site_col <- choose_col(
  data = stems_raw,
  user_col = site_col,
  candidate_names = c("site", "site_id", "siteid", "population", "plot", "stand"),
  col_label = "site"
)

resolved_sample_col <- choose_col(
  data = stems_raw,
  user_col = sample_col,
  candidate_names = c("sample_id", "individual_id", "id", "stem_id", "tree_id", "sample", "individual"),
  col_label = "sample id"
)

resolved_x_col <- choose_col(
  data = stems_raw,
  user_col = x_col,
  candidate_names = c("x", "x_m", "xcoord", "x_coord", "easting", "utm_x", "xpos"),
  col_label = "x coordinate"
)

resolved_y_col <- choose_col(
  data = stems_raw,
  user_col = y_col,
  candidate_names = c("y", "y_m", "ycoord", "y_coord", "northing", "utm_y", "ypos"),
  col_label = "y coordinate"
)

# 3) Standardize and validate.
stems <- stems_raw %>%
  transmute(
    site = as.character(.data[[resolved_site_col]]),
    sample_id = as.character(.data[[resolved_sample_col]]),
    x = suppressWarnings(as.numeric(.data[[resolved_x_col]])),
    y = suppressWarnings(as.numeric(.data[[resolved_y_col]]))
  )

if (any(is.na(stems$site) | stems$site == "")) {
  stop("Site column contains missing/blank values. Please fix input data.", call. = FALSE)
}

if (any(is.na(stems$sample_id) | stems$sample_id == "")) {
  stop("Sample ID column contains missing/blank values. Please fix input data.", call. = FALSE)
}

if (any(is.na(stems$x) | is.na(stems$y))) {
  bad_n <- sum(is.na(stems$x) | is.na(stems$y))
  stop(
    paste0(
      "Coordinate columns contain non-numeric or missing values in ", bad_n, " row(s).\n",
      "Please clean coordinates or update x/y column mappings in USER SETTINGS."
    ),
    call. = FALSE
  )
}

# Handle duplicated sample IDs within a site.
dupes <- stems %>%
  count(site, sample_id, name = "n") %>%
  filter(n > 1)

if (nrow(dupes) > 0) {
  stop(
    paste0(
      "Found duplicated sample IDs within site. Please ensure unique sample IDs per site.\n",
      "Examples: ",
      paste0(head(paste(dupes$site, dupes$sample_id, sep = ":"), 5), collapse = ", ")
    ),
    call. = FALSE
  )
}

# 4) Compute nearest-neighbor distances within each site.
insufficient_sites <- stems %>%
  count(site, name = "n") %>%
  filter(n < 2)

if (nrow(insufficient_sites) > 0) {
  warning(
    paste0(
      "Some sites have fewer than 2 sampled individuals; nearest-neighbor metrics will be NA for those sites: ",
      paste(insufficient_sites$site, collapse = ", ")
    ),
    call. = FALSE
  )
}

nn_table <- stems %>%
  arrange(site, sample_id) %>%
  group_by(site) %>%
  group_modify(~ compute_nn_within_site(.x)) %>%
  ungroup() %>%
  select(site, sample_id, nearest_neighbor_id, nearest_neighbor_distance_m)

# 5) Build site-level summary.
site_summary <- nn_table %>%
  group_by(site) %>%
  summarise(
    n_sampled_individuals = n(),
    mean_nearest_neighbor_distance_m = ifelse(all(is.na(nearest_neighbor_distance_m)), NA_real_, mean(nearest_neighbor_distance_m, na.rm = TRUE)),
    median_nearest_neighbor_distance_m = ifelse(all(is.na(nearest_neighbor_distance_m)), NA_real_, median(nearest_neighbor_distance_m, na.rm = TRUE)),
    min_nearest_neighbor_distance_m = ifelse(all(is.na(nearest_neighbor_distance_m)), NA_real_, min(nearest_neighbor_distance_m, na.rm = TRUE)),
    max_nearest_neighbor_distance_m = ifelse(all(is.na(nearest_neighbor_distance_m)), NA_real_, max(nearest_neighbor_distance_m, na.rm = TRUE)),
    sd_nearest_neighbor_distance_m = ifelse(all(is.na(nearest_neighbor_distance_m)), NA_real_, sd(nearest_neighbor_distance_m, na.rm = TRUE)),
    prop_nn_le_8m = ifelse(all(is.na(nearest_neighbor_distance_m)), NA_real_, mean(nearest_neighbor_distance_m <= 8, na.rm = TRUE)),
    prop_nn_le_10m = ifelse(all(is.na(nearest_neighbor_distance_m)), NA_real_, mean(nearest_neighbor_distance_m <= 10, na.rm = TRUE)),
    prop_nn_le_12m = ifelse(all(is.na(nearest_neighbor_distance_m)), NA_real_, mean(nearest_neighbor_distance_m <= 12, na.rm = TRUE)),
    .groups = "drop"
  )

# 6) Save outputs.
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

table_out_path <- file.path(output_dir, "nearest_neighbor_per_sample.csv")
summary_out_path <- file.path(output_dir, "nearest_neighbor_site_summary.csv")

readr::write_csv(nn_table, table_out_path)
readr::write_csv(site_summary, summary_out_path)

# 7) Create and save figures.
plot_theme <- theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95", color = "grey70"),
    strip.text = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )

hist_plot <- nn_table %>%
  filter(!is.na(nearest_neighbor_distance_m)) %>%
  ggplot(aes(x = nearest_neighbor_distance_m)) +
  geom_histogram(bins = 20, color = "white", fill = "#2C7FB8") +
  geom_vline(xintercept = 8, linetype = "dashed", linewidth = 0.8, color = "#D95F02") +
  facet_wrap(~ site, scales = "free_y") +
  labs(
    title = "Nearest-neighbor distances among sampled stems",
    x = "Nearest-neighbor distance (m)",
    y = "Count"
  ) +
  plot_theme

box_plot <- nn_table %>%
  filter(!is.na(nearest_neighbor_distance_m)) %>%
  ggplot(aes(x = site, y = nearest_neighbor_distance_m)) +
  geom_boxplot(fill = "#74A9CF", color = "#1F78B4", outlier.alpha = 0.6) +
  geom_hline(yintercept = 8, linetype = "dashed", linewidth = 0.8, color = "#D95F02") +
  labs(
    title = "Nearest-neighbor distance by site",
    x = "Site",
    y = "Nearest-neighbor distance (m)"
  ) +
  plot_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

hist_out_path <- file.path(output_dir, "nearest_neighbor_histogram_by_site.png")
box_out_path <- file.path(output_dir, "nearest_neighbor_boxplot_by_site.png")

ggsave(filename = hist_out_path, plot = hist_plot, width = 10, height = 7, dpi = 400)
ggsave(filename = box_out_path, plot = box_plot, width = 9, height = 6, dpi = 400)

# 8) Final console message.
message("Nearest-neighbor sampling check complete.")
message("Resolved columns:")
message("  site:      ", resolved_site_col)
message("  sample_id: ", resolved_sample_col)
message("  x:         ", resolved_x_col)
message("  y:         ", resolved_y_col)
message("Outputs saved to: ", normalizePath(output_dir, winslash = "/", mustWork = FALSE))
message("  - ", normalizePath(table_out_path, winslash = "/", mustWork = FALSE))
message("  - ", normalizePath(summary_out_path, winslash = "/", mustWork = FALSE))
message("  - ", normalizePath(hist_out_path, winslash = "/", mustWork = FALSE))
message("  - ", normalizePath(box_out_path, winslash = "/", mustWork = FALSE))