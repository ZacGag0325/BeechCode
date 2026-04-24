#!/usr/bin/env Rscript

# =============================================================================
# nearest_neighbor_sampling_check.R
# -----------------------------------------------------------------------------
# Standalone script to evaluate whether adjacent sample numbers within each site
# were geographically close, and to compare adjacent-number neighbors against
# true nearest geographic neighbors.
#
# This script:
#   1) Reads a user-specified Excel file.
#   2) Auto-detects (or uses user overrides for) site, sample number, and
#      coordinate columns.
#   3) Computes distance between adjacent sample numbers (n to n+1) by site.
#   4) Computes true nearest geographic neighbor for each sampled individual.
#   5) Summarizes results by site and writes output tables/figures.
#
# IMPORTANT:
#   - Standalone: does NOT modify any run-all or master pipeline.
#   - Edit USER SETTINGS before running.
# =============================================================================

# ------------------------------- USER SETTINGS -------------------------------
# Input Excel file provided by user.
input_excel <- "donnees_modifiees_west_summer2024 copie.xlsx"

# Optional explicit sheet name. If NULL, script will auto-select the best sheet.
sheet_name <- NULL

# Output folder (new, clearly named).
output_dir <- "outputs/nearest_neighbour_adjacent_numbers"

# Optional explicit column mappings. Leave as NULL to auto-detect.
site_col <- NULL
sample_col <- NULL
lat_col <- NULL
lon_col <- NULL
x_col <- NULL
y_col <- NULL

# Set TRUE to also save optional per-site point plots with adjacent-number lines.
make_site_maps <- TRUE

# Thresholds (meters) for adjacent-pair proportion summaries.
thresholds_m <- c(1, 2, 5, 8, 10, 20)
# -----------------------------------------------------------------------------

# ------------------------------- PACKAGE CHECK -------------------------------
required_pkgs <- c(
  "dplyr", "readr", "readxl", "stringr", "ggplot2",
  "tidyr", "purrr", "geosphere", "scales"
)
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_pkgs) > 0) {
  stop(
    paste0(
      "Missing required packages: ", paste(missing_pkgs, collapse = ", "),
      "\nInstall with, for example:\n",
      "install.packages(c(", paste(sprintf('"%s"', missing_pkgs), collapse = ", "), "))"
    ),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(readxl)
  library(stringr)
  library(ggplot2)
  library(tidyr)
  library(purrr)
  library(geosphere)
  library(scales)
})
# -----------------------------------------------------------------------------

# ------------------------------- HELPERS -------------------------------------
clean_names <- function(x) {
  x %>%
    stringr::str_trim() %>%
    stringr::str_to_lower() %>%
    stringr::str_replace_all("[^a-z0-9]+", "_") %>%
    stringr::str_replace_all("^_+|_+$", "")
}

choose_col <- function(data, user_col, candidate_patterns, label) {
  nms <- names(data)
  
  if (!is.null(user_col)) {
    if (!user_col %in% nms) {
      stop(
        paste0(
          "Configured ", label, " column '", user_col, "' not found.\n",
          "Available columns: ", paste(nms, collapse = ", ")
        ),
        call. = FALSE
      )
    }
    return(user_col)
  }
  
  nms_clean <- clean_names(nms)
  
  exact_idx <- which(nms_clean %in% candidate_patterns)
  if (length(exact_idx) == 1) return(nms[exact_idx])
  if (length(exact_idx) > 1) {
    stop(
      paste0(
        "Ambiguous auto-detection for ", label, " column.\n",
        "Multiple matches: ", paste(nms[exact_idx], collapse = ", "), "\n",
        "Please set '", gsub(" ", "_", label), "' manually in USER SETTINGS.\n",
        "All columns: ", paste(nms, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
  partial_idx <- unique(unlist(lapply(candidate_patterns, function(p) which(str_detect(nms_clean, fixed(p))))))
  if (length(partial_idx) == 1) return(nms[partial_idx])
  if (length(partial_idx) > 1) {
    stop(
      paste0(
        "Ambiguous auto-detection for ", label, " column.\n",
        "Multiple partial matches: ", paste(nms[partial_idx], collapse = ", "), "\n",
        "Please set '", gsub(" ", "_", label), "' manually in USER SETTINGS.\n",
        "All columns: ", paste(nms, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
  stop(
    paste0(
      "Could not auto-detect ", label, " column.\n",
      "Please set '", gsub(" ", "_", label), "' manually in USER SETTINGS.\n",
      "All columns: ", paste(nms, collapse = ", ")
    ),
    call. = FALSE
  )
}

sheet_score <- function(df) {
  nms <- clean_names(names(df))
  
  site_hits <- sum(nms %in% c("site", "site_id", "population", "pop", "stand", "plot"))
  sample_hits <- sum(nms %in% c("sample", "sample_id", "individual", "individual_id", "tree_id", "stem_id", "numero", "number", "id"))
  lat_hits <- sum(nms %in% c("lat", "latitude", "y_lat", "gps_lat", "coord_lat"))
  lon_hits <- sum(nms %in% c("lon", "long", "longitude", "x_lon", "gps_lon", "coord_lon"))
  x_hits <- sum(nms %in% c("x", "x_m", "x_coord", "easting", "utm_x"))
  y_hits <- sum(nms %in% c("y", "y_m", "y_coord", "northing", "utm_y"))
  
  coord_pair_hits <- ifelse(lat_hits > 0 & lon_hits > 0, 2L, 0L) + ifelse(x_hits > 0 & y_hits > 0, 2L, 0L)
  site_hits + sample_hits + coord_pair_hits
}

get_distance_matrix <- function(df_site, coord_type) {
  if (nrow(df_site) < 2) return(matrix(NA_real_, nrow(df_site), nrow(df_site)))
  
  if (coord_type == "latlon") {
    # Geodesic distances in meters (WGS84) using Haversine.
    geosphere::distm(
      x = as.matrix(df_site[, c("lon", "lat")]),
      fun = geosphere::distHaversine
    )
  } else {
    # Euclidean distances in projected coordinate units (expected meters).
    as.matrix(stats::dist(as.matrix(df_site[, c("x", "y")]), method = "euclidean"))
  }
}

pair_distance <- function(a_lon, a_lat, b_lon, b_lat, a_x, a_y, b_x, b_y, coord_type) {
  if (coord_type == "latlon") {
    geosphere::distHaversine(cbind(a_lon, a_lat), cbind(b_lon, b_lat))
  } else {
    sqrt((a_x - b_x)^2 + (a_y - b_y)^2)
  }
}
# -----------------------------------------------------------------------------

# ------------------------------- READ DATA -----------------------------------
if (!file.exists(input_excel)) {
  stop(
    paste0(
      "Input Excel file not found: ", input_excel,
      "\nPlease update 'input_excel' in USER SETTINGS."
    ),
    call. = FALSE
  )
}

xlsx_sheets <- readxl::excel_sheets(input_excel)
if (length(xlsx_sheets) == 0) {
  stop("No sheets found in the Excel file.", call. = FALSE)
}

if (is.null(sheet_name)) {
  sheet_meta <- purrr::map_dfr(xlsx_sheets, function(s) {
    tmp <- suppressMessages(readxl::read_excel(input_excel, sheet = s, n_max = 0))
    tibble(sheet = s, score = sheet_score(tmp), n_cols = ncol(tmp))
  })
  
  best_score <- max(sheet_meta$score)
  best_sheets <- sheet_meta %>% filter(score == best_score) %>% pull(sheet)
  
  if (best_score <= 1) {
    stop(
      paste0(
        "Could not confidently detect a sheet with site/sample/coordinates.\n",
        "Sheets and scores:\n",
        paste(paste0("  - ", sheet_meta$sheet, ": ", sheet_meta$score), collapse = "\n"),
        "\nPlease set 'sheet_name' manually in USER SETTINGS."
      ),
      call. = FALSE
    )
  }
  
  if (length(best_sheets) > 1) {
    stop(
      paste0(
        "Multiple sheets tie for best match: ", paste(best_sheets, collapse = ", "), "\n",
        "Please set 'sheet_name' manually in USER SETTINGS."
      ),
      call. = FALSE
    )
  }
  
  sheet_name <- best_sheets[[1]]
}

raw <- suppressMessages(readxl::read_excel(input_excel, sheet = sheet_name))

if (nrow(raw) == 0) stop("Selected sheet has no rows.", call. = FALSE)
if (ncol(raw) == 0) stop("Selected sheet has no columns.", call. = FALSE)

# ------------------------------- COLUMN RESOLUTION ---------------------------
site_col_res <- choose_col(
  raw,
  site_col,
  candidate_patterns = c("site", "site_id", "population", "pop", "stand", "plot"),
  label = "site"
)

sample_col_res <- choose_col(
  raw,
  sample_col,
  candidate_patterns = c("sample", "sample_id", "individual", "individual_id", "tree_id", "stem_id", "numero", "number", "id"),
  label = "sample"
)

# Detect coordinate mode: lat/lon preferred if clearly available; else x/y.
coord_mode <- NULL

if (!is.null(lat_col) && !is.null(lon_col)) {
  lat_col_res <- choose_col(raw, lat_col, c("lat", "latitude"), "latitude")
  lon_col_res <- choose_col(raw, lon_col, c("lon", "long", "longitude"), "longitude")
  coord_mode <- "latlon"
} else if (!is.null(x_col) && !is.null(y_col)) {
  x_col_res <- choose_col(raw, x_col, c("x", "easting", "utm_x"), "x coordinate")
  y_col_res <- choose_col(raw, y_col, c("y", "northing", "utm_y"), "y coordinate")
  coord_mode <- "xy"
} else {
  # Try automatic lat/lon first.
  lat_try <- tryCatch(
    choose_col(raw, NULL, c("lat", "latitude", "gps_lat", "coord_lat", "y_lat"), "latitude"),
    error = function(e) NULL
  )
  lon_try <- tryCatch(
    choose_col(raw, NULL, c("lon", "long", "longitude", "gps_lon", "coord_lon", "x_lon"), "longitude"),
    error = function(e) NULL
  )
  
  if (!is.null(lat_try) && !is.null(lon_try)) {
    lat_col_res <- lat_try
    lon_col_res <- lon_try
    coord_mode <- "latlon"
  } else {
    x_try <- tryCatch(
      choose_col(raw, NULL, c("x", "x_m", "x_coord", "easting", "utm_x"), "x coordinate"),
      error = function(e) NULL
    )
    y_try <- tryCatch(
      choose_col(raw, NULL, c("y", "y_m", "y_coord", "northing", "utm_y"), "y coordinate"),
      error = function(e) NULL
    )
    
    if (!is.null(x_try) && !is.null(y_try)) {
      x_col_res <- x_try
      y_col_res <- y_try
      coord_mode <- "xy"
    }
  }
}

if (is.null(coord_mode)) {
  stop(
    paste0(
      "Could not detect coordinate columns as either lat/lon or x/y.\n",
      "Available columns:\n",
      paste0("  - ", names(raw), collapse = "\n"),
      "\nSet 'lat_col'+'lon_col' OR 'x_col'+'y_col' in USER SETTINGS."
    ),
    call. = FALSE
  )
}

# ------------------------------- CLEAN DATA ----------------------------------
if (coord_mode == "latlon") {
  dat <- raw %>%
    transmute(
      site = as.character(.data[[site_col_res]]),
      sample_id_raw = as.character(.data[[sample_col_res]]),
      sample_num = suppressWarnings(as.numeric(as.character(.data[[sample_col_res]]))),
      lat = suppressWarnings(as.numeric(.data[[lat_col_res]])),
      lon = suppressWarnings(as.numeric(.data[[lon_col_res]]))
    ) %>%
    filter(!is.na(site), site != "", !is.na(sample_id_raw), sample_id_raw != "") %>%
    filter(!is.na(sample_num), !is.na(lat), !is.na(lon))
  
  # Lat/lon sanity check.
  if (nrow(dat) > 0 && any(dat$lat < -90 | dat$lat > 90 | dat$lon < -180 | dat$lon > 180)) {
    stop(
      "Detected latitude/longitude columns but found values outside valid ranges. Check coordinate columns.",
      call. = FALSE
    )
  }
} else {
  dat <- raw %>%
    transmute(
      site = as.character(.data[[site_col_res]]),
      sample_id_raw = as.character(.data[[sample_col_res]]),
      sample_num = suppressWarnings(as.numeric(as.character(.data[[sample_col_res]]))),
      x = suppressWarnings(as.numeric(.data[[x_col_res]])),
      y = suppressWarnings(as.numeric(.data[[y_col_res]]))
    ) %>%
    filter(!is.na(site), site != "", !is.na(sample_id_raw), sample_id_raw != "") %>%
    filter(!is.na(sample_num), !is.na(x), !is.na(y))
}

if (nrow(dat) == 0) {
  stop(
    "After cleaning, no rows remained with non-missing site, sample number, and coordinates.",
    call. = FALSE
  )
}

# Ensure sample numbers are unique per site for adjacent-number pair logic.
dup_sample_nums <- dat %>%
  count(site, sample_num, name = "n") %>%
  filter(n > 1)

if (nrow(dup_sample_nums) > 0) {
  stop(
    paste0(
      "Duplicate sample numbers found within site. Adjacent-number analysis requires unique sample numbers per site.\n",
      "Examples: ", paste(head(paste0(dup_sample_nums$site, "#", dup_sample_nums$sample_num), 10), collapse = ", ")
    ),
    call. = FALSE
  )
}

# ------------------------------- DISTANCE CALCS ------------------------------
# 1) True nearest geographic neighbor per individual within each site.
nearest_neighbor_table <- dat %>%
  arrange(site, sample_num) %>%
  group_by(site) %>%
  group_modify(function(.x, .y) {
    n <- nrow(.x)
    if (n < 2) {
      out <- .x %>%
        mutate(
          nearest_neighbor_sample_num = NA_real_,
          nearest_neighbor_sample_id = NA_character_,
          nearest_neighbor_distance_m = NA_real_,
          nearest_is_adjacent_number = NA
        )
      return(out)
    }
    
    dmat <- get_distance_matrix(.x, coord_mode)
    diag(dmat) <- Inf
    
    nn_idx <- max.col(-dmat, ties.method = "first")
    
    .x %>%
      mutate(
        nearest_neighbor_sample_num = sample_num[nn_idx],
        nearest_neighbor_sample_id = sample_id_raw[nn_idx],
        nearest_neighbor_distance_m = dmat[cbind(seq_len(n), nn_idx)],
        nearest_is_adjacent_number = abs(sample_num - nearest_neighbor_sample_num) == 1
      )
  }) %>%
  ungroup() %>%
  select(
    site,
    sample_num,
    sample_id_raw,
    nearest_neighbor_sample_num,
    nearest_neighbor_sample_id,
    nearest_neighbor_distance_m,
    nearest_is_adjacent_number
  )

# 2) Adjacent sample-number pairs (n to n+1) within each site.
left_tbl <- dat %>% rename(sample_num_a = sample_num, sample_id_a = sample_id_raw)
right_tbl <- dat %>% rename(sample_num_b = sample_num, sample_id_b = sample_id_raw)

adjacent_pairs <- left_tbl %>%
  mutate(sample_num_b = sample_num_a + 1) %>%
  inner_join(right_tbl, by = c("site", "sample_num_b")) %>%
  mutate(
    adjacent_distance_m = pair_distance(
      a_lon = if (coord_mode == "latlon") lon.x else NA_real_,
      a_lat = if (coord_mode == "latlon") lat.x else NA_real_,
      b_lon = if (coord_mode == "latlon") lon.y else NA_real_,
      b_lat = if (coord_mode == "latlon") lat.y else NA_real_,
      a_x = if (coord_mode == "xy") x.x else NA_real_,
      a_y = if (coord_mode == "xy") y.x else NA_real_,
      b_x = if (coord_mode == "xy") x.y else NA_real_,
      b_y = if (coord_mode == "xy") y.y else NA_real_,
      coord_type = coord_mode
    )
  )

# Compare adjacency pair vs true nearest-neighbor relation.
nn_lookup <- nearest_neighbor_table %>%
  select(site, sample_num, nearest_neighbor_sample_num)

adjacent_pairs <- adjacent_pairs %>%
  left_join(
    nn_lookup %>%
      rename(sample_num_a = sample_num, nn_of_a = nearest_neighbor_sample_num),
    by = c("site", "sample_num_a")
  ) %>%
  left_join(
    nn_lookup %>%
      rename(sample_num_b = sample_num, nn_of_b = nearest_neighbor_sample_num),
    by = c("site", "sample_num_b")
  ) %>%
  mutate(
    a_has_b_as_true_nn = nn_of_a == sample_num_b,
    b_has_a_as_true_nn = nn_of_b == sample_num_a,
    pair_any_true_nn_match = a_has_b_as_true_nn | b_has_a_as_true_nn,
    pair_mutual_true_nn = a_has_b_as_true_nn & b_has_a_as_true_nn
  ) %>%
  select(
    site,
    sample_num_a,
    sample_id_a,
    sample_num_b,
    sample_id_b,
    adjacent_distance_m,
    a_has_b_as_true_nn,
    b_has_a_as_true_nn,
    pair_any_true_nn_match,
    pair_mutual_true_nn
  ) %>%
  arrange(site, sample_num_a)

# ------------------------------- SITE SUMMARY --------------------------------
summary_base <- adjacent_pairs %>%
  group_by(site) %>%
  summarise(
    n_adjacent_pairs = n(),
    mean_adjacent_distance_m = mean(adjacent_distance_m, na.rm = TRUE),
    median_adjacent_distance_m = median(adjacent_distance_m, na.rm = TRUE),
    min_adjacent_distance_m = min(adjacent_distance_m, na.rm = TRUE),
    q1_adjacent_distance_m = quantile(adjacent_distance_m, 0.25, na.rm = TRUE),
    q3_adjacent_distance_m = quantile(adjacent_distance_m, 0.75, na.rm = TRUE),
    max_adjacent_distance_m = max(adjacent_distance_m, na.rm = TRUE),
    prop_pair_any_true_nn_match = mean(pair_any_true_nn_match, na.rm = TRUE),
    prop_pair_mutual_true_nn = mean(pair_mutual_true_nn, na.rm = TRUE),
    .groups = "drop"
  )

threshold_summary <- purrr::map_dfr(thresholds_m, function(th) {
  adjacent_pairs %>%
    group_by(site) %>%
    summarise(
      threshold_m = th,
      prop_pairs_within_threshold = mean(adjacent_distance_m <= th, na.rm = TRUE),
      .groups = "drop"
    )
})

threshold_wide <- threshold_summary %>%
  mutate(col_name = paste0("prop_adjacent_le_", threshold_m, "m")) %>%
  select(site, col_name, prop_pairs_within_threshold) %>%
  tidyr::pivot_wider(names_from = col_name, values_from = prop_pairs_within_threshold)

site_summary <- summary_base %>%
  left_join(threshold_wide, by = "site") %>%
  arrange(site)

# ------------------------------- OUTPUTS -------------------------------------
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

aadj_path <- file.path(output_dir, "adjacent_number_pair_distances.csv")
nn_path <- file.path(output_dir, "true_nearest_neighbour_results.csv")
summary_path <- file.path(output_dir, "adjacent_number_distance_summary_by_site.csv")
threshold_long_path <- file.path(output_dir, "adjacent_number_threshold_summary_long.csv")

readr::write_csv(adjacent_pairs, aadj_path)
readr::write_csv(nearest_neighbor_table, nn_path)
readr::write_csv(site_summary, summary_path)
readr::write_csv(threshold_summary, threshold_long_path)

# ------------------------------- FIGURES -------------------------------------
plot_theme <- theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95", color = "grey75"),
    strip.text = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )

# Figure 1: Distribution of adjacent-number distances by site.
fig1 <- adjacent_pairs %>%
  ggplot(aes(x = adjacent_distance_m)) +
  geom_histogram(bins = 20, color = "white", fill = "#2C7FB8") +
  facet_wrap(~ site, scales = "free_y") +
  labs(
    title = "Distribution of adjacent sample-number distances by site",
    x = "Distance between adjacent sample numbers (m)",
    y = "Count"
  ) +
  plot_theme

fig1_path <- file.path(output_dir, "figure_adjacent_distance_distribution_by_site.png")
ggsave(fig1_path, fig1, width = 11, height = 7, dpi = 400)

# Figure 2: Percentage of adjacent-number pairs within thresholds.
fig2 <- threshold_summary %>%
  ggplot(aes(x = factor(threshold_m), y = prop_pairs_within_threshold, group = site, color = site)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Adjacent-pair proximity across distance thresholds",
    x = "Distance threshold (m)",
    y = "Pairs within threshold (%)",
    color = "Site"
  ) +
  plot_theme

fig2_path <- file.path(output_dir, "figure_adjacent_pair_threshold_percentages.png")
ggsave(fig2_path, fig2, width = 10, height = 6, dpi = 400)

# Optional site plots: points + labels + lines connecting adjacent sample numbers.
if (isTRUE(make_site_maps)) {
  map_dir <- file.path(output_dir, "site_maps_adjacent_connections")
  dir.create(map_dir, recursive = TRUE, showWarnings = FALSE)
  
  map_input <- if (coord_mode == "latlon") {
    dat %>% mutate(plot_x = lon, plot_y = lat)
  } else {
    dat %>% mutate(plot_x = x, plot_y = y)
  }
  
  map_lines <- adjacent_pairs %>%
    left_join(
      map_input %>%
        select(site, sample_num, plot_x, plot_y) %>%
        rename(sample_num_a = sample_num, x_a = plot_x, y_a = plot_y),
      by = c("site", "sample_num_a")
    ) %>%
    left_join(
      map_input %>%
        select(site, sample_num, plot_x, plot_y) %>%
        rename(sample_num_b = sample_num, x_b = plot_x, y_b = plot_y),
      by = c("site", "sample_num_b")
    )
  
  site_ids <- sort(unique(map_input$site))
  
  for (s in site_ids) {
    p_points <- map_input %>% filter(site == s)
    p_lines <- map_lines %>% filter(site == s)
    
    p <- ggplot() +
      geom_segment(
        data = p_lines,
        aes(x = x_a, y = y_a, xend = x_b, yend = y_b),
        linewidth = 0.7,
        color = "#D95F02",
        alpha = 0.8
      ) +
      geom_point(
        data = p_points,
        aes(x = plot_x, y = plot_y),
        color = "#1F78B4",
        size = 2.2
      ) +
      geom_text(
        data = p_points,
        aes(x = plot_x, y = plot_y, label = sample_num),
        nudge_y = 0.02 * (max(p_points$plot_y, na.rm = TRUE) - min(p_points$plot_y, na.rm = TRUE) + 1e-9),
        size = 3
      ) +
      labs(
        title = paste0("Site: ", s, " (adjacent sample-number connections)"),
        x = ifelse(coord_mode == "latlon", "Longitude", "X"),
        y = ifelse(coord_mode == "latlon", "Latitude", "Y")
      ) +
      plot_theme
    
    ggsave(
      filename = file.path(map_dir, paste0("site_map_adjacent_", make.names(s), ".png")),
      plot = p,
      width = 7,
      height = 6,
      dpi = 400
    )
  }
}

# ------------------------------- CONSOLE MESSAGE -----------------------------
message("Adjacent-number nearest-neighbour analysis complete.")
message("Input file: ", normalizePath(input_excel, winslash = "/", mustWork = FALSE))
message("Sheet used: ", sheet_name)
message("Resolved columns:")
message("  site:   ", site_col_res)
message("  sample: ", sample_col_res)

if (coord_mode == "latlon") {
  message("  latitude:  ", lat_col_res)
  message("  longitude: ", lon_col_res)
  message("Distance method: geodesic (Haversine) in meters via geosphere::distHaversine.")
} else {
  message("  x: ", x_col_res)
  message("  y: ", y_col_res)
  message("Distance method: Euclidean distance in projected coordinate units (expected meters).")
}

message("Outputs saved to: ", normalizePath(output_dir, winslash = "/", mustWork = FALSE))
message("  - ", normalizePath(aadj_path, winslash = "/", mustWork = FALSE))
message("  - ", normalizePath(nn_path, winslash = "/", mustWork = FALSE))
message("  - ", normalizePath(summary_path, winslash = "/", mustWork = FALSE))
message("  - ", normalizePath(threshold_long_path, winslash = "/", mustWork = FALSE))
message("  - ", normalizePath(fig1_path, winslash = "/", mustWork = FALSE))
message("  - ", normalizePath(fig2_path, winslash = "/", mustWork = FALSE))
if (isTRUE(make_site_maps)) {
  message("  - ", normalizePath(file.path(output_dir, "site_maps_adjacent_connections"), winslash = "/", mustWork = FALSE))
}