#!/usr/bin/env Rscript

# ============================================================================
# nearest_neighbour_adjacent_numbers_analysis.R
# ----------------------------------------------------------------------------
# Standalone script for adjacent-sample and nearest-neighbour distance analysis
# from an Excel workbook.
#
# This script:
#   1) Reads an Excel file.
#   2) Finds the relevant sheet/columns (site, sample number, coordinates).
#   3) Cleans data and converts sample IDs / coordinates to numeric.
#   4) Computes adjacent-sample distances within each site.
#   5) Computes true nearest-neighbour distances within each site.
#   6) Summarizes site-level patterns and threshold proportions.
#   7) Exports CSV outputs and figures to a new output folder.
#
# IMPORTANT:
#   - Standalone script (not connected to any project pipeline).
#   - If auto-detection is ambiguous, the script prints column names and stops
#     with clear instructions for manual edits.
# ============================================================================

# ----------------------------- USER SETTINGS ---------------------------------
input_file <- "donnees_modifiees_west_summer2024 copie.xlsx"
output_dir <- "outputs/nearest_neighbour_adjacent_numbers"

# Optional manual overrides (set to exact column names if needed)
sheet_name_override <- NULL
site_col_override <- NULL
sample_col_override <- NULL
x_col_override <- NULL
y_col_override <- NULL
lat_col_override <- NULL
lon_col_override <- NULL

# Distance thresholds (meters)
thresholds_m <- c(1, 2, 5, 8, 10, 20)

# Optional per-site map/plot output (can be many files)
make_site_maps <- TRUE
# ----------------------------------------------------------------------------

# Resolve input path from several likely locations so the script can be run
# either from project root or from another working directory.
resolve_input_file <- function(path_in) {
  if (file.exists(path_in)) return(path_in)
  
  # Try relative to script directory and project-like subfolders.
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
    file.path(script_dir, "..", "inputs", path_in)
  ))
  
  hit <- candidates[file.exists(candidates)]
  if (length(hit) >= 1) return(normalizePath(hit[1], winslash = "/", mustWork = TRUE))
  
  # Fallback: search by basename in common folders
  base <- basename(path_in)
  fallback_dirs <- unique(c(getwd(), script_dir, file.path(script_dir, ".."), file.path(script_dir, "..", "data")))
  search_hits <- unlist(lapply(fallback_dirs, function(d) {
    if (!dir.exists(d)) return(character())
    list.files(d, pattern = paste0("^", gsub("([.()\\[\\]{}+*?^$|\\\\])", "\\\\\\1", base), "$"),
               full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
  }))
  search_hits <- unique(search_hits[file.exists(search_hits)])
  
  if (length(search_hits) >= 1) return(normalizePath(search_hits[1], winslash = "/", mustWork = TRUE))
  
  # Helpful error with discovered Excel files
  visible_excels <- unique(unlist(lapply(fallback_dirs, function(d) {
    if (!dir.exists(d)) return(character())
    list.files(d, pattern = "\\.(xlsx|xls)$", full.names = TRUE, recursive = FALSE, ignore.case = TRUE)
  })))
  
  stop(
    paste0(
      "Input Excel file not found from '", path_in, "'.\n",
      "Checked locations:\n  - ", paste(candidates, collapse = "\n  - "), "\n",
      if (length(visible_excels) > 0) {
        paste0("Excel files seen in common locations:\n  - ", paste(visible_excels, collapse = "\n  - "), "\n")
      } else {
        "No Excel files found in immediate common folders.\n"
      },
      "Please set 'input_file' to an absolute path or a correct relative path."
    ),
    call. = FALSE
  )
}

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

print_cols_and_stop <- function(sheet_name, cols, reason) {
  msg <- paste0(
    "\n", reason, "\n",
    "Sheet checked: '", sheet_name, "'\n",
    "Available columns:\n  - ", paste(cols, collapse = "\n  - "), "\n\n",
    "Please edit one or more of:\n",
    "  sheet_name_override, site_col_override, sample_col_override,\n",
    "  lat_col_override/lon_col_override OR x_col_override/y_col_override."
  )
  stop(msg, call. = FALSE)
}

pick_unique_column <- function(df, override, candidates, label, sheet_name) {
  nms <- names(df)
  if (!is.null(override)) {
    if (!(override %in% nms)) {
      stop(
        paste0("Configured ", label, " override '", override, "' not found."),
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
    print_cols_and_stop(
      sheet_name,
      nms,
      paste0("Ambiguous auto-detection for ", label, ": ", paste(exact_hits, collapse = ", "))
    )
  }
  
  partial_idx <- unique(unlist(lapply(candidates_norm, function(cd) {
    which(stringr::str_detect(nms_norm, stringr::fixed(cd)))
  })))
  
  partial_hits <- nms[partial_idx]
  if (length(partial_hits) == 1) return(partial_hits)
  if (length(partial_hits) > 1) {
    print_cols_and_stop(
      sheet_name,
      nms,
      paste0("Ambiguous partial match for ", label, ": ", paste(partial_hits, collapse = ", "))
    )
  }
  
  NA_character_
}

infer_coord_type <- function(df, lat_col, lon_col, x_col, y_col) {
  if (!is.na(lat_col) && !is.na(lon_col)) return("latlon")
  if (!is.na(x_col) && !is.na(y_col)) return("xy")
  stop("Could not determine coordinate type.", call. = FALSE)
}

haversine_m <- function(lat1, lon1, lat2, lon2) {
  # Geodesic distance in meters (Haversine, spherical Earth)
  r <- 6371000
  to_rad <- pi / 180
  dlat <- (lat2 - lat1) * to_rad
  dlon <- (lon2 - lon1) * to_rad
  a <- sin(dlat / 2)^2 + cos(lat1 * to_rad) * cos(lat2 * to_rad) * sin(dlon / 2)^2
  2 * r * atan2(sqrt(a), sqrt(1 - a))
}

pair_distance <- function(d, i, j) {
  if (nrow(d) == 0) return(NA_real_)
  if (d$coord_type[1] == "latlon") {
    haversine_m(d$lat[i], d$lon[i], d$lat[j], d$lon[j])
  } else {
    dx <- d$x[i] - d$x[j]
    dy <- d$y[i] - d$y[j]
    sqrt(dx^2 + dy^2)
  }
}

compute_site_outputs <- function(site_df) {
  site_df <- site_df %>% arrange(sample_num, sample_id)
  n <- nrow(site_df)
  
  # Adjacent-by-number pairs (1->2, 2->3, ... after sorting)
  adjacent_tbl <- if (n >= 2) {
    idx1 <- seq_len(n - 1)
    idx2 <- idx1 + 1
    tibble(
      site = site_df$site[1],
      sample_id_1 = site_df$sample_id[idx1],
      sample_id_2 = site_df$sample_id[idx2],
      sample_num_1 = site_df$sample_num[idx1],
      sample_num_2 = site_df$sample_num[idx2],
      adjacent_distance_m = map2_dbl(idx1, idx2, ~ pair_distance(site_df, .x, .y))
    )
  } else {
    tibble(
      site = character(), sample_id_1 = character(), sample_id_2 = character(),
      sample_num_1 = numeric(), sample_num_2 = numeric(), adjacent_distance_m = numeric()
    )
  }
  
  # True nearest neighbour for each individual
  nn_tbl <- if (n >= 2) {
    idx <- seq_len(n)
    nearest_idx <- map_int(idx, function(i) {
      dist_i <- map_dbl(idx, function(j) {
        if (i == j) return(Inf)
        pair_distance(site_df, i, j)
      })
      which.min(dist_i)
    })
    
    nearest_dist <- map2_dbl(idx, nearest_idx, ~ pair_distance(site_df, .x, .y))
    
    tibble(
      site = site_df$site,
      sample_id = site_df$sample_id,
      sample_num = site_df$sample_num,
      nearest_sample_id = site_df$sample_id[nearest_idx],
      nearest_sample_num = site_df$sample_num[nearest_idx],
      nearest_distance_m = nearest_dist
    )
  } else {
    tibble(
      site = site_df$site,
      sample_id = site_df$sample_id,
      sample_num = site_df$sample_num,
      nearest_sample_id = NA_character_,
      nearest_sample_num = NA_real_,
      nearest_distance_m = NA_real_
    )
  }
  
  # Is adjacent-number neighbour also true nearest neighbour?
  adjacent_match_tbl <- if (n >= 2) {
    tibble(i = seq_len(n)) %>%
      mutate(
        prev_i = i - 1,
        next_i = i + 1,
        prev_dist = ifelse(prev_i >= 1, map2_dbl(i, prev_i, ~ pair_distance(site_df, .x, .y)), Inf),
        next_dist = ifelse(next_i <= n, map2_dbl(i, next_i, ~ pair_distance(site_df, .x, .y)), Inf),
        adjacent_partner_i = ifelse(prev_dist <= next_dist, prev_i, next_i),
        adjacent_partner_i = ifelse(adjacent_partner_i < 1 | adjacent_partner_i > n, NA, adjacent_partner_i),
        adjacent_partner_id = ifelse(is.na(adjacent_partner_i), NA_character_, site_df$sample_id[adjacent_partner_i]),
        adjacent_partner_distance_m = ifelse(
          is.na(adjacent_partner_i),
          NA_real_,
          map2_dbl(i, adjacent_partner_i, ~ pair_distance(site_df, .x, .y))
        )
      ) %>%
      left_join(nn_tbl %>% select(sample_id, nearest_sample_id, nearest_distance_m), by = c("adjacent_partner_id" = "sample_id")) %>%
      transmute(
        site = site_df$site[1],
        sample_id = site_df$sample_id[i],
        sample_num = site_df$sample_num[i],
        adjacent_partner_id,
        adjacent_partner_distance_m,
        is_adjacent_partner_true_nearest = ifelse(is.na(adjacent_partner_id), NA, adjacent_partner_id == nearest_sample_id)
      )
  } else {
    tibble(
      site = site_df$site,
      sample_id = site_df$sample_id,
      sample_num = site_df$sample_num,
      adjacent_partner_id = NA_character_,
      adjacent_partner_distance_m = NA_real_,
      is_adjacent_partner_true_nearest = NA
    )
  }
  
  list(adjacent = adjacent_tbl, nearest = nn_tbl, adjacent_match = adjacent_match_tbl)
}

# ------------------------------ MAIN -----------------------------------------
input_file <- resolve_input_file(input_file)

sheets <- readxl::excel_sheets(input_file)
if (length(sheets) == 0) {
  stop("No sheets found in the Excel file.", call. = FALSE)
}

choose_sheet <- function() {
  if (!is.null(sheet_name_override)) {
    if (!(sheet_name_override %in% sheets)) {
      stop(
        paste0("sheet_name_override '", sheet_name_override, "' not found in workbook."),
        call. = FALSE
      )
    }
    return(sheet_name_override)
  }
  
  # Heuristic: pick first sheet where site+sample+coords can be detected.
  for (sh in sheets) {
    df_try <- suppressMessages(readxl::read_excel(input_file, sheet = sh))
    if (ncol(df_try) < 4) next
    
    site_try <- pick_unique_column(df_try, NULL,
                                   c("site", "population", "pop", "plot", "stand", "location"),
                                   "site", sh
    )
    sample_try <- pick_unique_column(df_try, NULL,
                                     c("sample", "sample_number", "sample_no", "individual", "ind", "id", "tree", "stem"),
                                     "sample", sh
    )
    
    lat_try <- pick_unique_column(df_try, NULL,
                                  c("lat", "latitude", "y_wgs84"), "latitude", sh
    )
    lon_try <- pick_unique_column(df_try, NULL,
                                  c("lon", "long", "longitude", "x_wgs84"), "longitude", sh
    )
    x_try <- pick_unique_column(df_try, NULL,
                                c("x", "easting", "utm_x", "coord_x"), "x", sh
    )
    y_try <- pick_unique_column(df_try, NULL,
                                c("y", "northing", "utm_y", "coord_y"), "y", sh
    )
    
    has_site_sample <- !is.na(site_try) && !is.na(sample_try)
    has_coord <- (!is.na(lat_try) && !is.na(lon_try)) || (!is.na(x_try) && !is.na(y_try))
    
    if (has_site_sample && has_coord) return(sh)
  }
  
  stop(
    paste0(
      "Could not auto-detect a sheet with site/sample/coordinate columns.\n",
      "Available sheets: ", paste(sheets, collapse = ", "), "\n",
      "Set sheet_name_override in USER SETTINGS."
    ),
    call. = FALSE
  )
}

sheet_use <- choose_sheet()
raw_df <- suppressMessages(readxl::read_excel(input_file, sheet = sheet_use))

site_col <- pick_unique_column(raw_df, site_col_override,
                               c("site", "population", "pop", "plot", "stand", "location"),
                               "site", sheet_use
)
sample_col <- pick_unique_column(raw_df, sample_col_override,
                                 c("sample", "sample_number", "sample_no", "individual", "ind", "id", "tree", "stem"),
                                 "sample", sheet_use
)

lat_col <- pick_unique_column(raw_df, lat_col_override,
                              c("lat", "latitude", "y_wgs84"), "latitude", sheet_use
)
lon_col <- pick_unique_column(raw_df, lon_col_override,
                              c("lon", "long", "longitude", "x_wgs84"), "longitude", sheet_use
)
x_col <- pick_unique_column(raw_df, x_col_override,
                            c("x", "easting", "utm_x", "coord_x"), "x", sheet_use
)
y_col <- pick_unique_column(raw_df, y_col_override,
                            c("y", "northing", "utm_y", "coord_y"), "y", sheet_use
)

coord_type <- infer_coord_type(raw_df, lat_col, lon_col, x_col, y_col)

if (is.na(site_col) || is.na(sample_col)) {
  print_cols_and_stop(sheet_use, names(raw_df), "Failed to detect required site/sample columns.")
}

if (coord_type == "latlon") {
  cleaned <- raw_df %>%
    transmute(
      site = as.character(.data[[site_col]]),
      sample_id = as.character(.data[[sample_col]]),
      sample_num = as_numeric_safely(.data[[sample_col]]),
      lat = as_numeric_safely(.data[[lat_col]]),
      lon = as_numeric_safely(.data[[lon_col]]),
      coord_type = "latlon"
    )
} else {
  cleaned <- raw_df %>%
    transmute(
      site = as.character(.data[[site_col]]),
      sample_id = as.character(.data[[sample_col]]),
      sample_num = as_numeric_safely(.data[[sample_col]]),
      x = as_numeric_safely(.data[[x_col]]),
      y = as_numeric_safely(.data[[y_col]]),
      coord_type = "xy"
    )
}

if (coord_type == "latlon") {
  cleaned <- cleaned %>%
    filter(!is.na(site), site != "", !is.na(sample_id), sample_id != "", !is.na(sample_num), !is.na(lat), !is.na(lon))
} else {
  cleaned <- cleaned %>%
    filter(!is.na(site), site != "", !is.na(sample_id), sample_id != "", !is.na(sample_num), !is.na(x), !is.na(y))
}

if (nrow(cleaned) == 0) {
  stop("No rows left after cleaning. Check column mappings and value formats.", call. = FALSE)
}

dup_check <- cleaned %>% count(site, sample_id) %>% filter(n > 1)
if (nrow(dup_check) > 0) {
  stop(
    paste0(
      "Duplicate sample_id values found within site after cleaning (first examples): ",
      paste(head(paste0(dup_check$site, ":", dup_check$sample_id), 10), collapse = ", "),
      "\nPlease ensure unique sample identifiers per site."
    ),
    call. = FALSE
  )
}

site_results <- cleaned %>%
  group_by(site) %>%
  group_split() %>%
  set_names(cleaned %>% distinct(site) %>% pull(site)) %>%
  map(compute_site_outputs)

adjacent_pairs <- bind_rows(map(site_results, "adjacent"))
nearest_results <- bind_rows(map(site_results, "nearest"))
adjacent_vs_nn <- bind_rows(map(site_results, "adjacent_match"))

summary_stats <- adjacent_pairs %>%
  group_by(site) %>%
  summarise(
    n_adjacent_pairs = n(),
    mean_adjacent_distance_m = mean(adjacent_distance_m, na.rm = TRUE),
    median_adjacent_distance_m = median(adjacent_distance_m, na.rm = TRUE),
    min_adjacent_distance_m = min(adjacent_distance_m, na.rm = TRUE),
    q25_adjacent_distance_m = quantile(adjacent_distance_m, probs = 0.25, na.rm = TRUE, names = FALSE),
    q75_adjacent_distance_m = quantile(adjacent_distance_m, probs = 0.75, na.rm = TRUE, names = FALSE),
    max_adjacent_distance_m = max(adjacent_distance_m, na.rm = TRUE),
    .groups = "drop"
  )

threshold_summary <- adjacent_pairs %>%
  tidyr::crossing(threshold_m = thresholds_m) %>%
  group_by(site, threshold_m) %>%
  summarise(prop_within_threshold = mean(adjacent_distance_m <= threshold_m, na.rm = TRUE), .groups = "drop")

adj_vs_nn_summary <- adjacent_vs_nn %>%
  group_by(site) %>%
  summarise(
    n_individuals = n(),
    prop_adjacent_partner_is_true_nearest = mean(is_adjacent_partner_true_nearest, na.rm = TRUE),
    .groups = "drop"
  )

site_summary <- summary_stats %>%
  left_join(adj_vs_nn_summary, by = "site")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

readr::write_csv(adjacent_pairs, file.path(output_dir, "adjacent_number_pair_distances.csv"))
readr::write_csv(nearest_results, file.path(output_dir, "true_nearest_neighbour_by_individual.csv"))
readr::write_csv(adjacent_vs_nn, file.path(output_dir, "adjacent_vs_true_nearest_by_individual.csv"))
readr::write_csv(site_summary, file.path(output_dir, "summary_by_site.csv"))
readr::write_csv(threshold_summary, file.path(output_dir, "threshold_summary_by_site.csv"))

p1 <- ggplot(adjacent_pairs, aes(x = site, y = adjacent_distance_m)) +
  geom_boxplot(fill = "#9ecae1", color = "#2171b5", outlier.alpha = 0.6) +
  theme_bw(base_size = 12) +
  labs(
    title = "Adjacent sample-number distances by site",
    x = "Site",
    y = "Distance between adjacent sample numbers (m)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  filename = file.path(output_dir, "adjacent_distance_distribution_by_site.png"),
  plot = p1,
  width = 10,
  height = 6,
  dpi = 320
)

p2 <- ggplot(threshold_summary, aes(x = factor(threshold_m), y = prop_within_threshold * 100, fill = site)) +
  geom_col(position = position_dodge()) +
  theme_bw(base_size = 12) +
  labs(
    title = "Percent of adjacent-number pairs within thresholds",
    x = "Distance threshold (m)",
    y = "Pairs within threshold (%)",
    fill = "Site"
  )

ggsave(
  filename = file.path(output_dir, "adjacent_pairs_threshold_percentages.png"),
  plot = p2,
  width = 11,
  height = 6,
  dpi = 320
)

if (make_site_maps) {
  dir.create(file.path(output_dir, "site_maps"), recursive = TRUE, showWarnings = FALSE)
  
  map_data <- if (coord_type == "latlon") {
    cleaned %>% mutate(px = lon, py = lat)
  } else {
    cleaned %>% mutate(px = x, py = y)
  }
  
  invisible(
    map_data %>%
      group_by(site) %>%
      group_split() %>%
      walk(function(d_site) {
        d_site <- d_site %>% arrange(sample_num)
        line_df <- if (nrow(d_site) >= 2) {
          tibble(
            x = head(d_site$px, -1),
            y = head(d_site$py, -1),
            xend = tail(d_site$px, -1),
            yend = tail(d_site$py, -1)
          )
        } else {
          tibble(x = numeric(), y = numeric(), xend = numeric(), yend = numeric())
        }
        
        p_site <- ggplot() +
          geom_segment(data = line_df, aes(x = x, y = y, xend = xend, yend = yend),
                       color = "grey55", linewidth = 0.5) +
          geom_point(data = d_site, aes(x = px, y = py), color = "#08519c", size = 2) +
          geom_text(data = d_site, aes(x = px, y = py, label = sample_num), vjust = -0.7, size = 2.8) +
          theme_bw(base_size = 11) +
          labs(
            title = paste0("Site ", d_site$site[1], ": adjacent-number connections"),
            x = ifelse(coord_type == "latlon", "Longitude", "X"),
            y = ifelse(coord_type == "latlon", "Latitude", "Y")
          )
        
        out <- file.path(output_dir, "site_maps", paste0("site_", make.names(d_site$site[1]), "_adjacent_map.png"))
        ggsave(out, p_site, width = 7, height = 6, dpi = 320)
      })
  )
}

message("Analysis complete.")
message("Input file: ", normalizePath(input_file, mustWork = FALSE))
message("Sheet used: ", sheet_use)
message("Resolved columns:")
message("  site: ", site_col)
message("  sample: ", sample_col)
if (coord_type == "latlon") {
  message("  latitude: ", lat_col)
  message("  longitude: ", lon_col)
  message("Distance method: geodesic (Haversine) in meters for latitude/longitude coordinates.")
} else {
  message("  x: ", x_col)
  message("  y: ", y_col)
  message("Distance method: Euclidean in meters for projected X/Y coordinates.")
}
message("Rows kept after cleaning: ", nrow(cleaned))
message("Output directory: ", normalizePath(output_dir, mustWork = FALSE))