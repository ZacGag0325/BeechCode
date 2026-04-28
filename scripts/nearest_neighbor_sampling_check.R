#!/usr/bin/env Rscript

# ============================================================================
# nearest_neighbor_sampling_check.R
# ----------------------------------------------------------------------------
# Combined standalone script for TRUE nearest-neighbour distance analysis.
#
# This script:
#   1) Reads stem-level coordinates from CSV/TSV/TXT or Excel.
#   2) Auto-detects site/sample/coordinate columns (or uses user overrides).
#   3) Computes TRUE nearest neighbour per stem within each site:
#      - all pairwise distances within site
#      - self excluded
#      - keep the minimum distance for each stem
#   4) Exports tables and figures for sampling-check interpretation.
#
# IMPORTANT:
#   - This script does NOT use adjacent sample-number pairs.
#   - It uses spatial coordinates to compute nearest neighbours.
# ============================================================================

# ----------------------------- USER SETTINGS ---------------------------------
input_file <- "data/stem_level_data.csv"
output_dir <- "outputs/nearest_neighbor_sampling_check"

# Optional Excel sheet override (only used for .xlsx/.xls)
sheet_name_override <- NULL

# Optional column overrides (set exact column names to force mapping)
site_col_override <- NULL
sample_col_override <- NULL
x_col_override <- NULL
y_col_override <- NULL
lat_col_override <- NULL
lon_col_override <- NULL

# Thresholds for console summaries
thresholds_m <- c(1, 2, 5, 8)
# ----------------------------------------------------------------------------

# ----------------------------- PACKAGE CHECK ---------------------------------
required_pkgs <- c("dplyr", "readr", "stringr", "ggplot2", "tidyr", "purrr", "readxl")
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
  library(readxl)
})
# ----------------------------------------------------------------------------

# ----------------------------- HELPERS ---------------------------------------
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

pick_unique_column <- function(df, override, candidates, label) {
  nms <- names(df)
  
  if (!is.null(override)) {
    if (!(override %in% nms)) {
      stop(
        paste0("Configured ", label, " override '", override, "' not found.\n",
               "Available columns: ", paste(nms, collapse = ", ")),
        call. = FALSE
      )
    }
    return(override)
  }
  
  nms_norm <- normalize_names(nms)
  candidates_norm <- normalize_names(candidates)
  
  exact_hits <- nms[nms_norm %in% candidates_norm]
  if (length(exact_hits) == 1) return(exact_hits)
  if (length(exact_hits) > 1) {
    stop(
      paste0("Ambiguous auto-detection for ", label, ": ", paste(exact_hits, collapse = ", "),
             "\nPlease set an override in USER SETTINGS."),
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
      paste0("Ambiguous partial match for ", label, ": ", paste(partial_hits, collapse = ", "),
             "\nPlease set an override in USER SETTINGS."),
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

pair_distance <- function(site_df, i, j) {
  if (site_df$coord_type[1] == "latlon") {
    haversine_m(site_df$lat[i], site_df$lon[i], site_df$lat[j], site_df$lon[j])
  } else {
    dx <- site_df$x[i] - site_df$x[j]
    dy <- site_df$y[i] - site_df$y[j]
    sqrt(dx^2 + dy^2)
  }
}

compute_true_nn_within_site <- function(site_df) {
  n <- nrow(site_df)
  
  if (n < 2) {
    return(
      site_df %>%
        mutate(
          nearest_neighbor_id = NA_character_,
          nearest_neighbor_distance_m = NA_real_
        )
    )
  }
  
  idx <- seq_len(n)
  
  nearest_idx <- map_int(idx, function(i) {
    dists <- map_dbl(idx, function(j) {
      if (i == j) return(Inf)
      pair_distance(site_df, i, j)
    })
    which.min(dists)
  })
  
  nearest_dist <- map2_dbl(idx, nearest_idx, ~ pair_distance(site_df, .x, .y))
  
  site_df %>%
    mutate(
      nearest_neighbor_id = sample_id[nearest_idx],
      nearest_neighbor_distance_m = nearest_dist
    )
}

resolve_input_data <- function(path_in, sheet_override = NULL) {
  if (!file.exists(path_in)) {
    stop(
      paste0("Input file not found: ", path_in,
             "\nPlease update 'input_file' in USER SETTINGS."),
      call. = FALSE
    )
  }
  
  ext <- tolower(tools::file_ext(path_in))
  
  if (ext == "csv") {
    data <- readr::read_csv(path_in, show_col_types = FALSE)
    return(list(data = data, source = basename(path_in), sheet = NA_character_))
  }
  
  if (ext %in% c("tsv", "txt")) {
    data <- readr::read_tsv(path_in, show_col_types = FALSE)
    return(list(data = data, source = basename(path_in), sheet = NA_character_))
  }
  
  if (ext %in% c("xlsx", "xls")) {
    sheets <- readxl::excel_sheets(path_in)
    if (length(sheets) == 0) stop("No sheets found in Excel file.", call. = FALSE)
    
    sheet_use <- if (!is.null(sheet_override)) {
      if (!(sheet_override %in% sheets)) {
        stop(
          paste0("sheet_name_override '", sheet_override, "' not found in workbook."),
          call. = FALSE
        )
      }
      sheet_override
    } else {
      sheets[1]
    }
    
    data <- readxl::read_excel(path_in, sheet = sheet_use)
    return(list(data = data, source = basename(path_in), sheet = sheet_use))
  }
  
  stop(
    paste0("Unsupported extension '.", ext, "'. Use csv/tsv/txt/xlsx/xls."),
    call. = FALSE
  )
}
# ----------------------------------------------------------------------------

# ----------------------------- MAIN ------------------------------------------
input_obj <- resolve_input_data(input_file, sheet_name_override)
raw_df <- input_obj$data

site_col <- pick_unique_column(
  raw_df,
  site_col_override,
  c("site", "site_id", "siteid", "population", "pop", "plot", "stand", "location"),
  "site"
)

sample_col <- pick_unique_column(
  raw_df,
  sample_col_override,
  c("sample", "sample_id", "sample_number", "sample_no", "individual", "individual_id", "ind", "id", "tree", "stem", "stem_id", "numero"),
  "sample"
)

lat_col <- pick_unique_column(
  raw_df,
  lat_col_override,
  c("lat", "latitude", "y_wgs84"),
  "latitude"
)

lon_col <- pick_unique_column(
  raw_df,
  lon_col_override,
  c("lon", "long", "longitude", "x_wgs84"),
  "longitude"
)

x_col <- pick_unique_column(
  raw_df,
  x_col_override,
  c("x", "x_m", "xcoord", "x_coord", "easting", "utm_x", "coord_x"),
  "x"
)

y_col <- pick_unique_column(
  raw_df,
  y_col_override,
  c("y", "y_m", "ycoord", "y_coord", "northing", "utm_y", "coord_y"),
  "y"
)

coord_type <- if (!is.na(lat_col) && !is.na(lon_col)) {
  "latlon"
} else if (!is.na(x_col) && !is.na(y_col)) {
  "xy"
} else {
  stop(
    paste0(
      "Could not detect coordinate columns.\n",
      "Need either lat/lon or x/y columns.\n",
      "Available columns: ", paste(names(raw_df), collapse = ", ")
    ),
    call. = FALSE
  )
}

stems <- if (coord_type == "latlon") {
  raw_df %>%
    transmute(
      site = as.character(.data[[site_col]]),
      sample_id = as.character(.data[[sample_col]]),
      lat = as_numeric_safely(.data[[lat_col]]),
      lon = as_numeric_safely(.data[[lon_col]]),
      coord_type = "latlon"
    ) %>%
    filter(!is.na(site), site != "", !is.na(sample_id), sample_id != "", !is.na(lat), !is.na(lon))
} else {
  raw_df %>%
    transmute(
      site = as.character(.data[[site_col]]),
      sample_id = as.character(.data[[sample_col]]),
      x = as_numeric_safely(.data[[x_col]]),
      y = as_numeric_safely(.data[[y_col]]),
      coord_type = "xy"
    ) %>%
    filter(!is.na(site), site != "", !is.na(sample_id), sample_id != "", !is.na(x), !is.na(y))
}

if (nrow(stems) == 0) {
  stop("No rows left after cleaning. Check column mapping and coordinate values.", call. = FALSE)
}

dupes <- stems %>%
  count(site, sample_id, name = "n") %>%
  filter(n > 1)

if (nrow(dupes) > 0) {
  stop(
    paste0(
      "Found duplicated sample IDs within site. Please ensure uniqueness.\n",
      "Examples: ", paste(head(paste(dupes$site, dupes$sample_id, sep = ":"), 10), collapse = ", ")
    ),
    call. = FALSE
  )
}

nn_table <- stems %>%
  arrange(site, sample_id) %>%
  group_by(site) %>%
  group_modify(~ compute_true_nn_within_site(.x)) %>%
  ungroup() %>%
  select(site, sample_id, nearest_neighbor_id, nearest_neighbor_distance_m)

site_summary <- nn_table %>%
  group_by(site) %>%
  summarise(
    n_stems = n(),
    mean_nearest_neighbor_distance_m = ifelse(all(is.na(nearest_neighbor_distance_m)), NA_real_, mean(nearest_neighbor_distance_m, na.rm = TRUE)),
    median_nearest_neighbor_distance_m = ifelse(all(is.na(nearest_neighbor_distance_m)), NA_real_, median(nearest_neighbor_distance_m, na.rm = TRUE)),
    n_within_1m = ifelse(all(is.na(nearest_neighbor_distance_m)), NA_real_, sum(nearest_neighbor_distance_m <= 1, na.rm = TRUE)),
    pct_within_1m = ifelse(all(is.na(nearest_neighbor_distance_m)), NA_real_, mean(nearest_neighbor_distance_m <= 1, na.rm = TRUE) * 100),
    n_within_2m = ifelse(all(is.na(nearest_neighbor_distance_m)), NA_real_, sum(nearest_neighbor_distance_m <= 2, na.rm = TRUE)),
    pct_within_2m = ifelse(all(is.na(nearest_neighbor_distance_m)), NA_real_, mean(nearest_neighbor_distance_m <= 2, na.rm = TRUE) * 100),
    n_within_5m = ifelse(all(is.na(nearest_neighbor_distance_m)), NA_real_, sum(nearest_neighbor_distance_m <= 5, na.rm = TRUE)),
    pct_within_5m = ifelse(all(is.na(nearest_neighbor_distance_m)), NA_real_, mean(nearest_neighbor_distance_m <= 5, na.rm = TRUE) * 100),
    n_within_8m = ifelse(all(is.na(nearest_neighbor_distance_m)), NA_real_, sum(nearest_neighbor_distance_m <= 8, na.rm = TRUE)),
    pct_within_8m = ifelse(all(is.na(nearest_neighbor_distance_m)), NA_real_, mean(nearest_neighbor_distance_m <= 8, na.rm = TRUE) * 100),
    .groups = "drop"
  )

# Save outputs
# Keep existing outputs and add required file names.
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

legacy_table_out <- file.path(output_dir, "nearest_neighbor_per_sample.csv")
legacy_summary_out <- file.path(output_dir, "nearest_neighbor_site_summary.csv")
true_nn_stem_out <- file.path(output_dir, "true_nearest_neighbour_distance_by_stem.csv")

readr::write_csv(nn_table, legacy_table_out)
readr::write_csv(site_summary, legacy_summary_out)
readr::write_csv(nn_table, true_nn_stem_out)

# Histogram-like class plot: 1 m bins up to 8 m + >8 m
bin_labels <- c("0–1 m", "1–2 m", "2–3 m", "3–4 m", "4–5 m", "5–6 m", "6–7 m", "7–8 m", ">8 m")

histogram_data <- nn_table %>%
  filter(!is.na(nearest_neighbor_distance_m), nearest_neighbor_distance_m >= 0) %>%
  mutate(
    distance_class = cut(
      nearest_neighbor_distance_m,
      breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, Inf),
      labels = bin_labels,
      include.lowest = TRUE,
      right = FALSE
    )
  ) %>%
  count(distance_class, name = "n_stems", .drop = FALSE) %>%
  mutate(distance_class = factor(distance_class, levels = bin_labels))

plot_theme <- theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold")
  )

hist_plot <- ggplot(histogram_data, aes(x = distance_class, y = n_stems)) +
  geom_col(fill = "#2C7FB8", color = "white") +
  labs(
    title = "Distribution of true nearest-neighbour distances among sampled stems",
    x = "True nearest-neighbour distance (m)",
    y = "Number of stems",
    caption = "Bin width = 1 m; final class aggregates stems with nearest-neighbour distance >8 m."
  ) +
  plot_theme +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

box_plot <- nn_table %>%
  filter(!is.na(nearest_neighbor_distance_m), nearest_neighbor_distance_m >= 0) %>%
  ggplot(aes(x = site, y = nearest_neighbor_distance_m)) +
  geom_boxplot(fill = "#74A9CF", color = "#1F78B4", outlier.alpha = 0.6) +
  labs(
    title = "True nearest-neighbour distance by site",
    x = "Site",
    y = "True nearest-neighbour distance (m)"
  ) +
  plot_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

hist_out <- file.path(output_dir, "true_nearest_neighbour_distance_histogram.png")
box_out <- file.path(output_dir, "true_nearest_neighbour_distance_by_site.png")

ggsave(filename = hist_out, plot = hist_plot, width = 10, height = 7, dpi = 400)
ggsave(filename = box_out, plot = box_plot, width = 9, height = 6, dpi = 400)

# Console summaries
nn_non_missing <- nn_table %>%
  filter(!is.na(nearest_neighbor_distance_m), nearest_neighbor_distance_m >= 0)

total_stems <- nrow(nn_non_missing)

count_within <- function(x, threshold) sum(x <= threshold, na.rm = TRUE)
fmt_count_pct <- function(count, total) {
  if (total == 0) return("0 (0.0%)")
  paste0(count, " (", sprintf("%.1f", (count / total) * 100), "%)")
}

message("Nearest-neighbour sampling check complete.")
message("Resolved source: ", input_obj$source)
if (!is.na(input_obj$sheet)) message("Resolved sheet: ", input_obj$sheet)
message("Resolved columns:")
message("  site: ", site_col)
message("  sample: ", sample_col)
if (coord_type == "latlon") {
  message("  latitude: ", lat_col)
  message("  longitude: ", lon_col)
  message("Distance method: geodesic (Haversine) meters for latitude/longitude coordinates.")
} else {
  message("  x: ", x_col)
  message("  y: ", y_col)
  message("Distance method: Euclidean meters for projected X/Y coordinates.")
}

message("Total number of stems: ", total_stems)
message("Median nearest-neighbour distance overall (m): ", sprintf("%.3f", median(nn_non_missing$nearest_neighbor_distance_m, na.rm = TRUE)))
message("Mean nearest-neighbour distance overall (m): ", sprintf("%.3f", mean(nn_non_missing$nearest_neighbor_distance_m, na.rm = TRUE)))
for (th in thresholds_m) {
  n_th <- count_within(nn_non_missing$nearest_neighbor_distance_m, th)
  message("Number and percent of stems within ", th, " m: ", fmt_count_pct(n_th, total_stems))
}

message("Site-level nearest-neighbour summaries:")
for (i in seq_len(nrow(site_summary))) {
  this_site <- site_summary$site[i]
  this_n <- site_summary$n_stems[i]
  this_med <- site_summary$median_nearest_neighbor_distance_m[i]
  this_mean <- site_summary$mean_nearest_neighbor_distance_m[i]
  
  message("  Site ", this_site, ":")
  message("    stems = ", this_n)
  message("    median (m) = ", sprintf("%.3f", this_med))
  message("    mean (m) = ", sprintf("%.3f", this_mean))
  message("    within 1 m = ", fmt_count_pct(site_summary$n_within_1m[i], this_n))
  message("    within 2 m = ", fmt_count_pct(site_summary$n_within_2m[i], this_n))
  message("    within 5 m = ", fmt_count_pct(site_summary$n_within_5m[i], this_n))
  message("    within 8 m = ", fmt_count_pct(site_summary$n_within_8m[i], this_n))
}

message("Outputs saved to: ", normalizePath(output_dir, winslash = "/", mustWork = FALSE))
message("  - ", normalizePath(true_nn_stem_out, winslash = "/", mustWork = FALSE))
message("  - ", normalizePath(hist_out, winslash = "/", mustWork = FALSE))
message("  - ", normalizePath(box_out, winslash = "/", mustWork = FALSE))
message("  - ", normalizePath(legacy_table_out, winslash = "/", mustWork = FALSE), " (kept)")
message("  - ", normalizePath(legacy_summary_out, winslash = "/", mustWork = FALSE), " (kept)")