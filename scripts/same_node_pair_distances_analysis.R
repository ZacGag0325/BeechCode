#!/usr/bin/env Rscript

# ============================================================================
# same_node_pair_distances_analysis.R
# ----------------------------------------------------------------------------
# Standalone script to compute distances only between stems sampled at the
# same node/grid point (within site), with special focus on nodes containing
# exactly two stems.
#
# Outputs (saved to outputs/same_node_pair_distances/):
#   - same_node_pair_distances.csv
#   - same_node_summary_by_site.csv
#   - same_node_threshold_summary_by_site.csv
#   - same_node_distance_distribution_by_site.png
#   - same_node_threshold_percentages.png
# ============================================================================

# ----------------------------- USER SETTINGS ---------------------------------
input_file <- "donnees_modifiees_west_summer2024 copie.xlsx"
output_dir <- "outputs/same_node_pair_distances"

# Optional manual overrides (exact names in the input table)
sheet_name_override <- NULL
site_col_override <- NULL
sample_col_override <- NULL
node_col_override <- NULL
x_col_override <- NULL
y_col_override <- NULL
lat_col_override <- NULL
lon_col_override <- NULL

# Distance thresholds (meters)
thresholds_m <- c(0.5, 1, 2, 5, 8)
# ----------------------------------------------------------------------------

required_pkgs <- c("readxl", "dplyr", "stringr", "ggplot2", "tidyr", "purrr", "readr")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    paste0(
      "Missing packages: ", paste(missing_pkgs, collapse = ", "), "\n",
      "Install with: install.packages(c(",
      paste(sprintf('"%s"', missing_pkgs), collapse = ", "),
      "))"
    ),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(tidyr)
  library(purrr)
  library(readr)
})

resolve_input_file <- function(path_in) {
  if (file.exists(path_in)) return(path_in)
  
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
  script_dir <- if (length(file_arg) > 0) {
    dirname(normalizePath(file_arg[1], winslash = "/", mustWork = FALSE))
  } else {
    getwd()
  }
  
  candidates <- unique(c(
    path_in,
    file.path(getwd(), path_in),
    file.path(script_dir, path_in),
    file.path(script_dir, "..", path_in),
    file.path(script_dir, "..", "data", path_in),
    file.path(script_dir, "..", "data", "raw", path_in),
    file.path(script_dir, "..", "inputs", path_in)
  ))
  
  hit <- candidates[file.exists(candidates)]
  if (length(hit) >= 1) return(normalizePath(hit[1], winslash = "/", mustWork = TRUE))
  
  stop(
    paste0(
      "Input file not found from '", path_in, "'.\n",
      "Checked:\n  - ", paste(candidates, collapse = "\n  - "), "\n",
      "Set 'input_file' to the correct path before running."
    ),
    call. = FALSE
  )
}

normalize_names <- function(x) {
  x %>%
    stringr::str_to_lower() %>%
    stringr::str_replace_all("[^a-z0-9]+", "_") %>%
    stringr::str_replace_all("^_+|_+$", "")
}

as_numeric_safely <- function(x) {
  if (is.numeric(x)) return(as.numeric(x))
  x_chr <- as.character(x)
  x_chr <- stringr::str_replace_all(x_chr, ",", ".")
  suppressWarnings(as.numeric(x_chr))
}

pick_unique_column <- function(df, override, candidates, label, sheet_name) {
  nms <- names(df)
  
  if (!is.null(override)) {
    if (!(override %in% nms)) {
      stop(paste0("Configured ", label, " override '", override, "' not found."), call. = FALSE)
    }
    return(override)
  }
  
  nms_norm <- normalize_names(nms)
  candidates_norm <- normalize_names(candidates)
  
  exact_hits <- nms[nms_norm %in% candidates_norm]
  if (length(exact_hits) == 1) return(exact_hits)
  if (length(exact_hits) > 1) {
    stop(
      paste0(
        "Ambiguous ", label, " column detection in sheet '", sheet_name, "': ",
        paste(exact_hits, collapse = ", "), "\n",
        "Please set the corresponding *_override setting explicitly."
      ),
      call. = FALSE
    )
  }
  
  partial_idx <- unique(unlist(lapply(candidates_norm, function(cd) {
    which(stringr::str_detect(nms_norm, stringr::fixed(cd)))
  })))
  partial_hits <- nms[partial_idx]
  
  if (length(partial_hits) == 1) return(partial_hits)
  if (length(partial_hits) > 1) {
    stop(
      paste0(
        "Ambiguous partial match for ", label, " in sheet '", sheet_name, "': ",
        paste(partial_hits, collapse = ", "), "\n",
        "Please set the corresponding *_override setting explicitly."
      ),
      call. = FALSE
    )
  }
  
  NA_character_
}

haversine_m <- function(lat1, lon1, lat2, lon2) {
  r <- 6371000
  to_rad <- pi / 180
  dlat <- (lat2 - lat1) * to_rad
  dlon <- (lon2 - lon1) * to_rad
  a <- sin(dlat / 2)^2 + cos(lat1 * to_rad) * cos(lat2 * to_rad) * sin(dlon / 2)^2
  2 * r * atan2(sqrt(a), sqrt(1 - a))
}

pair_distance_m <- function(coord_type, x1, y1, x2, y2) {
  if (coord_type == "latlon") {
    haversine_m(lat1 = y1, lon1 = x1, lat2 = y2, lon2 = x2)
  } else {
    sqrt((x2 - x1)^2 + (y2 - y1)^2)
  }
}

# ------------------------------ MAIN -----------------------------------------
input_file <- resolve_input_file(input_file)

sheets <- readxl::excel_sheets(input_file)
if (length(sheets) == 0) {
  stop("No sheets found in the input workbook.", call. = FALSE)
}

choose_sheet <- function() {
  if (!is.null(sheet_name_override)) {
    if (!(sheet_name_override %in% sheets)) {
      stop(paste0("sheet_name_override '", sheet_name_override, "' not found."), call. = FALSE)
    }
    return(sheet_name_override)
  }
  
  for (sh in sheets) {
    df_try <- suppressMessages(readxl::read_excel(input_file, sheet = sh))
    if (ncol(df_try) < 5) next
    
    site_try <- pick_unique_column(df_try, NULL,
                                   c("site", "site_id", "population", "plot", "stand", "location"),
                                   "site", sh)
    sample_try <- pick_unique_column(df_try, NULL,
                                     c("sample", "sample_id", "sample_number", "individual", "stem", "id_tige", "tree", "numero", "numero_individu"),
                                     "sample", sh)
    node_try <- pick_unique_column(df_try, NULL,
                                   c("node", "grid_point", "point", "sampling_point", "quadrat_point", "pair_id", "subplot", "sous_parcelle", "second_node", "double_sample_point"),
                                   "node/grid point", sh)
    
    lat_try <- pick_unique_column(df_try, NULL, c("lat", "latitude", "y_wgs84"), "latitude", sh)
    lon_try <- pick_unique_column(df_try, NULL, c("lon", "long", "longitude", "x_wgs84"), "longitude", sh)
    x_try <- pick_unique_column(df_try, NULL, c("x", "easting", "utm_x", "coord_x"), "x", sh)
    y_try <- pick_unique_column(df_try, NULL, c("y", "northing", "utm_y", "coord_y"), "y", sh)
    
    has_required <- !is.na(site_try) && !is.na(sample_try) && !is.na(node_try)
    has_coords <- (!is.na(lat_try) && !is.na(lon_try)) || (!is.na(x_try) && !is.na(y_try))
    if (has_required && has_coords) return(sh)
  }
  
  stop(
    paste0(
      "Could not auto-detect a sheet containing site + sample + node + coordinates.\n",
      "Set sheet_name_override and column overrides manually."
    ),
    call. = FALSE
  )
}

sheet_use <- choose_sheet()
raw_df <- suppressMessages(readxl::read_excel(input_file, sheet = sheet_use))

site_col <- pick_unique_column(raw_df, site_col_override,
                               c("site", "site_id", "population", "plot", "stand", "location"),
                               "site", sheet_use)
sample_col <- pick_unique_column(raw_df, sample_col_override,
                                 c("sample", "sample_id", "sample_number", "individual", "stem", "id_tige", "tree", "numero", "numero_individu"),
                                 "sample", sheet_use)
node_col <- pick_unique_column(raw_df, node_col_override,
                               c("node", "grid_point", "point", "sampling_point", "quadrat_point", "pair_id", "subplot", "sous_parcelle", "second_node", "double_sample_point"),
                               "node/grid point", sheet_use)

lat_col <- pick_unique_column(raw_df, lat_col_override, c("lat", "latitude", "y_wgs84"), "latitude", sheet_use)
lon_col <- pick_unique_column(raw_df, lon_col_override, c("lon", "long", "longitude", "x_wgs84"), "longitude", sheet_use)
x_col <- pick_unique_column(raw_df, x_col_override, c("x", "easting", "utm_x", "coord_x"), "x", sheet_use)
y_col <- pick_unique_column(raw_df, y_col_override, c("y", "northing", "utm_y", "coord_y"), "y", sheet_use)

if (is.na(site_col) || is.na(sample_col) || is.na(node_col)) {
  stop(
    paste0(
      "Required columns not detected. Found: site=", site_col,
      ", sample=", sample_col, ", node=", node_col, "\n",
      "Please set *_override values in USER SETTINGS."
    ),
    call. = FALSE
  )
}

coord_type <- if (!is.na(lat_col) && !is.na(lon_col)) {
  "latlon"
} else if (!is.na(x_col) && !is.na(y_col)) {
  "xy"
} else {
  stop("Could not detect a coordinate pair (lat/lon or x/y).", call. = FALSE)
}

cleaned <- if (coord_type == "latlon") {
  raw_df %>%
    transmute(
      site = as.character(.data[[site_col]]),
      sample_id = as.character(.data[[sample_col]]),
      node_id = as.character(.data[[node_col]]),
      x = as_numeric_safely(.data[[lon_col]]),
      y = as_numeric_safely(.data[[lat_col]])
    )
} else {
  raw_df %>%
    transmute(
      site = as.character(.data[[site_col]]),
      sample_id = as.character(.data[[sample_col]]),
      node_id = as.character(.data[[node_col]]),
      x = as_numeric_safely(.data[[x_col]]),
      y = as_numeric_safely(.data[[y_col]])
    )
}

cleaned <- cleaned %>%
  filter(!is.na(site), site != "", !is.na(sample_id), sample_id != "", !is.na(node_id), node_id != "", !is.na(x), !is.na(y))

if (nrow(cleaned) == 0) {
  stop("No usable rows after cleaning.", call. = FALSE)
}

pair_tbl <- cleaned %>%
  group_by(site, node_id) %>%
  arrange(sample_id, .by_group = TRUE) %>%
  mutate(node_n_stems = n()) %>%
  filter(node_n_stems == 2) %>%
  summarise(
    sample_id_1 = sample_id[1],
    sample_id_2 = sample_id[2],
    x1 = x[1], y1 = y[1],
    x2 = x[2], y2 = y[2],
    distance_m = pair_distance_m(coord_type, x1, y1, x2, y2),
    .groups = "drop"
  )

summary_by_site <- pair_tbl %>%
  group_by(site) %>%
  summarise(
    n_same_node_pairs = n(),
    median_same_node_distance_m = median(distance_m, na.rm = TRUE),
    mean_same_node_distance_m = mean(distance_m, na.rm = TRUE),
    min_same_node_distance_m = min(distance_m, na.rm = TRUE),
    max_same_node_distance_m = max(distance_m, na.rm = TRUE),
    pct_within_1m = mean(distance_m <= 1, na.rm = TRUE) * 100,
    pct_within_2m = mean(distance_m <= 2, na.rm = TRUE) * 100,
    pct_within_5m = mean(distance_m <= 5, na.rm = TRUE) * 100,
    pct_within_8m = mean(distance_m <= 8, na.rm = TRUE) * 100,
    .groups = "drop"
  )

threshold_summary <- pair_tbl %>%
  tidyr::crossing(threshold_m = thresholds_m) %>%
  group_by(site, threshold_m) %>%
  summarise(
    n_pairs = n(),
    n_within_threshold = sum(distance_m <= threshold_m, na.rm = TRUE),
    pct_within_threshold = mean(distance_m <= threshold_m, na.rm = TRUE) * 100,
    .groups = "drop"
  )

# Include sites with zero same-node pairs in summaries when possible
all_sites <- cleaned %>% distinct(site)
summary_by_site <- all_sites %>%
  left_join(summary_by_site, by = "site") %>%
  mutate(
    n_same_node_pairs = ifelse(is.na(n_same_node_pairs), 0L, as.integer(n_same_node_pairs))
  )

threshold_summary <- all_sites %>%
  tidyr::crossing(threshold_m = thresholds_m) %>%
  left_join(threshold_summary, by = c("site", "threshold_m")) %>%
  mutate(
    n_pairs = ifelse(is.na(n_pairs), 0L, as.integer(n_pairs)),
    n_within_threshold = ifelse(is.na(n_within_threshold), 0L, as.integer(n_within_threshold)),
    pct_within_threshold = ifelse(is.na(pct_within_threshold), NA_real_, pct_within_threshold)
  )

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

pair_csv <- file.path(output_dir, "same_node_pair_distances.csv")
summary_csv <- file.path(output_dir, "same_node_summary_by_site.csv")
threshold_csv <- file.path(output_dir, "same_node_threshold_summary_by_site.csv")
plot_dist_png <- file.path(output_dir, "same_node_distance_distribution_by_site.png")
plot_thr_png <- file.path(output_dir, "same_node_threshold_percentages.png")

readr::write_csv(pair_tbl, pair_csv)
readr::write_csv(summary_by_site, summary_csv)
readr::write_csv(threshold_summary, threshold_csv)

plot_theme <- theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

p_dist <- pair_tbl %>%
  ggplot(aes(x = site, y = distance_m)) +
  geom_boxplot(fill = "#9ecae1", color = "#2171b5", outlier.alpha = 0.6, na.rm = TRUE) +
  labs(
    title = "Distance between two stems sampled at the same node",
    x = "Site",
    y = "Same-node pair distance (m)"
  ) +
  plot_theme

ggsave(plot_dist_png, p_dist, width = 10, height = 6, dpi = 320)

p_thr <- threshold_summary %>%
  ggplot(aes(x = factor(threshold_m), y = pct_within_threshold, fill = site)) +
  geom_col(position = position_dodge(), na.rm = TRUE) +
  labs(
    title = "Percent of same-node pairs within distance thresholds",
    x = "Distance threshold (m)",
    y = "Pairs within threshold (%)",
    fill = "Site"
  ) +
  plot_theme

ggsave(plot_thr_png, p_thr, width = 11, height = 6, dpi = 320)

# Required console messages
message("Same-node pair distance analysis complete.")
message("Input file: ", normalizePath(input_file, winslash = "/", mustWork = FALSE))
message("Sheet used: ", sheet_use)
message("Resolved columns: site=", site_col, ", sample=", sample_col, ", node=", node_col)
if (coord_type == "latlon") {
  message("Coordinates: latitude=", lat_col, ", longitude=", lon_col, " (Haversine distance, meters)")
} else {
  message("Coordinates: x=", x_col, ", y=", y_col, " (Euclidean distance, meters)")
}

per_site_console <- summary_by_site %>%
  arrange(site)

if (nrow(per_site_console) == 0) {
  message("No sites found after cleaning.")
} else {
  message("Per-site same-node pair summary:")
  for (i in seq_len(nrow(per_site_console))) {
    row <- per_site_console[i, ]
    message(
      "  Site ", row$site,
      ": pairs=", row$n_same_node_pairs,
      ", median_m=", ifelse(is.na(row$median_same_node_distance_m), "NA", sprintf("%.3f", row$median_same_node_distance_m)),
      ", <=1m=", ifelse(is.na(row$pct_within_1m), "NA", sprintf("%.1f%%", row$pct_within_1m)),
      ", <=2m=", ifelse(is.na(row$pct_within_2m), "NA", sprintf("%.1f%%", row$pct_within_2m)),
      ", <=5m=", ifelse(is.na(row$pct_within_5m), "NA", sprintf("%.1f%%", row$pct_within_5m)),
      ", <=8m=", ifelse(is.na(row$pct_within_8m), "NA", sprintf("%.1f%%", row$pct_within_8m))
    )
  }
}

message("Outputs saved to: ", normalizePath(output_dir, winslash = "/", mustWork = FALSE))
message("  - ", normalizePath(pair_csv, winslash = "/", mustWork = FALSE))
message("  - ", normalizePath(summary_csv, winslash = "/", mustWork = FALSE))
message("  - ", normalizePath(threshold_csv, winslash = "/", mustWork = FALSE))
message("  - ", normalizePath(plot_dist_png, winslash = "/", mustWork = FALSE))
message("  - ", normalizePath(plot_thr_png, winslash = "/", mustWork = FALSE))