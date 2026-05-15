# scripts/sampled_stem_descriptive_results.R
############################################################
# Article-ready descriptive summaries for sampled American beech stems
#
# Components:
# 1) Spatial distances among sampled stems, including true nearest-neighbor
#    distances within each site and intended paired-node distances when a
#    paired-sample column can be detected.
# 2) Diameter / developmental-stage summaries for sampled stems in the full
#    genetic/clonality dataset (gi; ramets, not clone-corrected gi_mll).
#
# Main outputs:
# - outputs/tables/sampled_stem_distance_summary.csv
# - outputs/tables/sampled_stem_distance_by_site.csv
# - outputs/tables/sampled_stem_nearest_neighbor_pairs.csv
# - outputs/tables/sampled_stem_intended_pair_distances.csv (if pair metadata exists)
# - outputs/tables/sampled_stem_developmental_stage_summary.csv
# - outputs/tables/sampled_stem_developmental_stage_by_site.csv
# - outputs/tables/sampled_stem_dbh_distribution.csv
# - outputs/figures/sampled_stem_nearest_neighbor_distance_histogram.png
# - outputs/figures/sampled_stem_nearest_neighbor_distance_by_site.png
# - outputs/figures/sampled_stem_diameter_distribution.png
# - outputs/figures/sampled_stem_developmental_stage_by_site.png
############################################################

suppressPackageStartupMessages({
  library(adegenet)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readxl)
})

source("scripts/_load_objects.R")

SCRIPT_TAG <- "[sampled_stem_descriptive_results]"
SITE_LOOKUP_SHEET <- "site_lookup"
GENETIC_SHEET_HINT <- "genetique"
EXPECTED_SITE_LABELS <- c("S1", "S2", "S3", "S4", "S5", "S6", "N1", "N2", "N3", "N4", "N5", "N6")

DISTANCE_CLASS_LEVELS <- c("0-1 m", "1-2 m", "2-5 m", ">5 m")
DIAMETER_CLASS_LEVELS <- c("<1 cm", "1-2 cm", "2-5 cm", "5-10 cm", ">10 cm")

# Optional manual override. Leave NULL to auto-discover the workbook in data/raw
# or the project root. In this project, the workbook has typically been named
# something like "donnees_modifiees_west_summer2024 copie.xlsx".
INPUT_WORKBOOK_OVERRIDE <- NULL
COORDINATE_SHEET_OVERRIDE <- NULL

# ------------------------------------------------------------------
# Generic helpers
# ------------------------------------------------------------------
normalize_lookup_key <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("\uFEFF", "", x, fixed = TRUE)
  x <- gsub("[[:cntrl:]]", "", x)
  x[is.na(x)] <- ""
  x
}

normalize_ascii_key <- function(x) {
  x <- iconv(as.character(x), from = "", to = "ASCII//TRANSLIT")
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x[is.na(x)] <- ""
  x
}

coerce_numeric <- function(x) {
  if (is.numeric(x)) return(as.numeric(x))
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "N/A", "NaN", "NULL", "null", "na")] <- NA_character_
  x <- gsub(",", ".", x, fixed = TRUE)
  x <- gsub("[^0-9eE+.-]+", "", x)
  suppressWarnings(as.numeric(x))
}

safe_mean <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (length(x) == 0 || all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

safe_median <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (length(x) == 0 || all(is.na(x))) return(NA_real_)
  stats::median(x, na.rm = TRUE)
}

safe_min <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (length(x) == 0 || all(is.na(x))) return(NA_real_)
  min(x, na.rm = TRUE)
}

safe_max <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (length(x) == 0 || all(is.na(x))) return(NA_real_)
  max(x, na.rm = TRUE)
}

safe_sd <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[!is.na(x)]
  if (length(x) <= 1) return(NA_real_)
  stats::sd(x)
}

format_pct <- function(n, denom) {
  if (is.na(denom) || denom == 0) return(NA_real_)
  100 * n / denom
}

pick_column <- function(df, candidates, label = "column", required = FALSE) {
  nms <- names(df)
  nms_norm <- normalize_ascii_key(nms)
  candidate_norm <- normalize_ascii_key(candidates)
  
  idx <- match(candidate_norm, nms_norm, nomatch = 0)
  idx <- idx[idx > 0]
  if (length(idx) > 0) return(nms[idx[1]])
  
  if (required) {
    stop(
      SCRIPT_TAG, " Could not detect required ", label,
      ". Accepted names include: ", paste(candidates, collapse = ", "),
      "\nAvailable columns: ", paste(nms, collapse = ", "),
      call. = FALSE
    )
  }
  
  NA_character_
}

write_csv_msg <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.csv(df, path, row.names = FALSE, na = "")
  message(SCRIPT_TAG, " Saved: ", path)
  invisible(path)
}

site_order_tbl <- function() {
  data.frame(Site = EXPECTED_SITE_LABELS, Site_order = seq_along(EXPECTED_SITE_LABELS), stringsAsFactors = FALSE)
}

order_by_expected_site <- function(df) {
  if (!"Site" %in% names(df)) return(df)
  df %>%
    mutate(Site = factor(as.character(Site), levels = EXPECTED_SITE_LABELS)) %>%
    arrange(Site) %>%
    mutate(Site = as.character(Site))
}

# ------------------------------------------------------------------
# Workbook and site_lookup helpers
# ------------------------------------------------------------------
list_candidate_workbooks <- function() {
  search_dirs <- unique(c(
    file.path(PROJECT_ROOT, "data", "raw"),
    PROJECT_ROOT,
    file.path(PROJECT_ROOT, "data")
  ))
  search_dirs <- search_dirs[dir.exists(search_dirs)]
  unique(unlist(lapply(search_dirs, function(d) {
    list.files(d, pattern = "\\.(xlsx|xls)$", full.names = TRUE, recursive = FALSE, ignore.case = TRUE)
  })))
}

workbook_score <- function(path) {
  sheets <- tryCatch(readxl::excel_sheets(path), error = function(e) character(0))
  if (length(sheets) == 0) return(-Inf)
  score <- 0
  sheet_norm <- normalize_ascii_key(sheets)
  base_norm <- normalize_ascii_key(basename(path))
  if (normalize_ascii_key(SITE_LOOKUP_SHEET) %in% sheet_norm) score <- score + 100
  if (normalize_ascii_key(GENETIC_SHEET_HINT) %in% sheet_norm) score <- score + 50
  if (grepl("donnee|modifie|west|summer|copie", base_norm)) score <- score + 20
  score
}

find_input_workbook <- function() {
  if (!is.null(INPUT_WORKBOOK_OVERRIDE)) {
    if (!file.exists(INPUT_WORKBOOK_OVERRIDE)) {
      stop(SCRIPT_TAG, " INPUT_WORKBOOK_OVERRIDE does not exist: ", INPUT_WORKBOOK_OVERRIDE, call. = FALSE)
    }
    return(normalizePath(INPUT_WORKBOOK_OVERRIDE, winslash = "/", mustWork = TRUE))
  }
  
  candidates <- list_candidate_workbooks()
  if (length(candidates) == 0) {
    stop(
      SCRIPT_TAG, " No Excel workbook found in data/raw, project root, or data/. ",
      "Set INPUT_WORKBOOK_OVERRIDE at the top of this script.",
      call. = FALSE
    )
  }
  
  scores <- vapply(candidates, workbook_score, numeric(1))
  candidates <- candidates[order(scores, decreasing = TRUE)]
  scores <- scores[order(scores, decreasing = TRUE)]
  
  if (!is.finite(scores[1]) || scores[1] < 0) {
    stop(SCRIPT_TAG, " Excel workbooks were found, but none could be read.", call. = FALSE)
  }
  
  if (sum(scores == scores[1]) > 1) {
    message(SCRIPT_TAG, " Multiple equally scored workbooks found; using: ", basename(candidates[1]))
  } else {
    message(SCRIPT_TAG, " Using workbook: ", basename(candidates[1]))
  }
  
  normalizePath(candidates[1], winslash = "/", mustWork = TRUE)
}

load_site_lookup <- function(workbook) {
  sheets <- readxl::excel_sheets(workbook)
  sheet_idx <- match(normalize_ascii_key(SITE_LOOKUP_SHEET), normalize_ascii_key(sheets))
  if (is.na(sheet_idx)) {
    message(SCRIPT_TAG, " site_lookup sheet not found; retaining existing site labels where possible.")
    return(NULL)
  }
  
  lookup <- suppressMessages(readxl::read_excel(workbook, sheet = sheets[sheet_idx])) %>%
    as.data.frame(stringsAsFactors = FALSE)
  
  old_site_col <- pick_column(lookup, c("Site", "site", "site_code", "old_site", "code_site"), "old site code", required = TRUE)
  label_col <- pick_column(lookup, c("Site_label", "site_label", "new_site", "site_new", "label", "site_id"), "display site label")
  order_col <- pick_column(lookup, c("Site_order", "site_order", "order", "ordre", "sort", "south_north_order"), "site order")
  region_col <- pick_column(lookup, c("Region", "region"), "region")
  
  if (is.na(label_col)) label_col <- old_site_col
  
  out <- lookup %>%
    transmute(
      raw_site = normalize_lookup_key(.data[[old_site_col]]),
      Site = normalize_lookup_key(.data[[label_col]]),
      Site_order = if (!is.na(order_col)) coerce_numeric(.data[[order_col]]) else match(normalize_lookup_key(.data[[label_col]]), EXPECTED_SITE_LABELS),
      Region = if (!is.na(region_col)) normalize_lookup_key(.data[[region_col]]) else NA_character_
    ) %>%
    filter(nzchar(raw_site), nzchar(Site)) %>%
    distinct(raw_site, .keep_all = TRUE)
  
  if (nrow(out) == 0) return(NULL)
  message(SCRIPT_TAG, " Loaded site_lookup rows: ", nrow(out))
  out
}

map_site_labels <- function(site_values, lookup) {
  raw <- normalize_lookup_key(site_values)
  if (is.null(lookup)) return(raw)
  idx <- match(raw, lookup$raw_site)
  mapped <- lookup$Site[idx]
  out <- ifelse(is.na(mapped) | !nzchar(mapped), raw, mapped)
  out
}

validate_expected_sites <- function(site_values, context) {
  observed <- unique(as.character(site_values))
  if (!setequal(observed, EXPECTED_SITE_LABELS)) {
    stop(
      SCRIPT_TAG, " ", context, " must contain exactly the expected site labels.",
      "\nExpected: ", paste(EXPECTED_SITE_LABELS, collapse = ", "),
      "\nObserved: ", paste(sort(observed), collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

# ------------------------------------------------------------------
# Coordinate / nearest-neighbor distance analysis
# ------------------------------------------------------------------
site_candidates <- c("Site", "site", "site_code", "Code_site", "code_site", "population", "Population", "pop", "localite", "localité", "Numéro_Population", "Numero_Population")
sample_candidates <- c("ind_id", "Individual", "individual", "individual_id", "sample", "sample_id", "Sample_ID", "stem", "stem_id", "id_tige", "ID_tige", "numero", "num", "no", "numero_individu", "Nom_Labo_Échantillons", "Nom_Labo_Echantillons")
lat_candidates <- c("lat", "latitude", "Latitude", "LAT", "y_wgs84", "Y_WGS84")
lon_candidates <- c("lon", "long", "longitude", "Longitude", "LONG", "x_wgs84", "X_WGS84")
x_candidates <- c("x", "X", "x_m", "X_m", "xcoord", "x_coord", "coord_x", "Coord_X", "easting", "utm_x", "UTM_X")
y_candidates <- c("y", "Y", "y_m", "Y_m", "ycoord", "y_coord", "coord_y", "Coord_Y", "northing", "utm_y", "UTM_Y")
node_candidates <- c("node", "Node", "grid_node", "sampling_node", "point", "Point", "station", "Station", "noeud", "Noeud", "noeud_echantillonnage")
pair_role_candidates <- c("role", "Role", "sample_type", "type", "Type", "paired", "pair", "neighbor", "neighbour", "voisin", "Voisin", "nearest_neighbor", "nearest_neighbour")

sheet_coordinate_score <- function(df) {
  if (ncol(df) < 4) return(-Inf)
  site_col <- pick_column(df, site_candidates, required = FALSE)
  sample_col <- pick_column(df, sample_candidates, required = FALSE)
  lat_col <- pick_column(df, lat_candidates, required = FALSE)
  lon_col <- pick_column(df, lon_candidates, required = FALSE)
  x_col <- pick_column(df, x_candidates, required = FALSE)
  y_col <- pick_column(df, y_candidates, required = FALSE)
  has_site_sample <- !is.na(site_col) && !is.na(sample_col)
  has_latlon <- !is.na(lat_col) && !is.na(lon_col)
  has_xy <- !is.na(x_col) && !is.na(y_col)
  if (!has_site_sample || (!has_latlon && !has_xy)) return(-Inf)
  score <- 100
  if (has_latlon) score <- score + 20
  if (has_xy) score <- score + 20
  if (!is.na(pick_column(df, node_candidates, required = FALSE))) score <- score + 5
  if (!is.na(pick_column(df, pair_role_candidates, required = FALSE))) score <- score + 5
  score
}

read_coordinate_sheet <- function(workbook) {
  sheets <- readxl::excel_sheets(workbook)
  
  if (!is.null(COORDINATE_SHEET_OVERRIDE)) {
    if (!COORDINATE_SHEET_OVERRIDE %in% sheets) {
      stop(SCRIPT_TAG, " COORDINATE_SHEET_OVERRIDE not found: ", COORDINATE_SHEET_OVERRIDE, call. = FALSE)
    }
    df <- suppressMessages(readxl::read_excel(workbook, sheet = COORDINATE_SHEET_OVERRIDE)) %>% as.data.frame(stringsAsFactors = FALSE)
    return(list(sheet = COORDINATE_SHEET_OVERRIDE, data = df))
  }
  
  candidates <- lapply(sheets, function(sh) {
    df <- suppressMessages(readxl::read_excel(workbook, sheet = sh)) %>% as.data.frame(stringsAsFactors = FALSE)
    list(sheet = sh, data = df, score = sheet_coordinate_score(df))
  })
  
  scores <- vapply(candidates, function(x) x$score, numeric(1))
  if (all(!is.finite(scores))) {
    message(SCRIPT_TAG, " No sheet with site/sample/coordinate columns was detected; distance outputs will be skipped.")
    return(NULL)
  }
  
  best <- candidates[[which.max(scores)]]
  message(SCRIPT_TAG, " Coordinate source sheet: ", best$sheet)
  best
}

haversine_m <- function(lat1, lon1, lat2, lon2) {
  r <- 6371000
  to_rad <- pi / 180
  dlat <- (lat2 - lat1) * to_rad
  dlon <- (lon2 - lon1) * to_rad
  a <- sin(dlat / 2)^2 + cos(lat1 * to_rad) * cos(lat2 * to_rad) * sin(dlon / 2)^2
  2 * r * atan2(sqrt(a), sqrt(1 - a))
}

pair_distance_m <- function(site_df, i, j) {
  if (site_df$coordinate_type[1] == "latlon") {
    haversine_m(site_df$latitude[i], site_df$longitude[i], site_df$latitude[j], site_df$longitude[j])
  } else {
    sqrt((site_df$x_m[i] - site_df$x_m[j])^2 + (site_df$y_m[i] - site_df$y_m[j])^2)
  }
}

prepare_coordinate_data <- function(coord_source, lookup) {
  if (is.null(coord_source)) return(NULL)
  
  df <- coord_source$data
  site_col <- pick_column(df, site_candidates, "site", required = TRUE)
  sample_col <- pick_column(df, sample_candidates, "sample/stem ID", required = TRUE)
  lat_col <- pick_column(df, lat_candidates, "latitude", required = FALSE)
  lon_col <- pick_column(df, lon_candidates, "longitude", required = FALSE)
  x_col <- pick_column(df, x_candidates, "x coordinate", required = FALSE)
  y_col <- pick_column(df, y_candidates, "y coordinate", required = FALSE)
  node_col <- pick_column(df, node_candidates, "sampling node", required = FALSE)
  role_col <- pick_column(df, pair_role_candidates, "paired-sample role", required = FALSE)
  
  has_latlon <- !is.na(lat_col) && !is.na(lon_col)
  has_xy <- !is.na(x_col) && !is.na(y_col)
  
  if (!has_latlon && !has_xy) {
    message(SCRIPT_TAG, " Coordinate sheet lacks usable lat/lon or x/y columns; distance outputs will be skipped.")
    return(NULL)
  }
  
  use_latlon <- has_latlon
  if (use_latlon) {
    message(SCRIPT_TAG, " Distance coordinate columns used: latitude=", lat_col, ", longitude=", lon_col, " (Haversine meters).")
    out <- df %>%
      transmute(
        Site_raw = normalize_lookup_key(.data[[site_col]]),
        Site = map_site_labels(.data[[site_col]], lookup),
        sample_id = normalize_lookup_key(.data[[sample_col]]),
        latitude = coerce_numeric(.data[[lat_col]]),
        longitude = coerce_numeric(.data[[lon_col]]),
        x_m = NA_real_,
        y_m = NA_real_,
        coordinate_type = "latlon",
        sampling_node = if (!is.na(node_col)) normalize_lookup_key(.data[[node_col]]) else NA_character_,
        pair_role = if (!is.na(role_col)) normalize_lookup_key(.data[[role_col]]) else NA_character_
      ) %>%
      filter(nzchar(Site), nzchar(sample_id), !is.na(latitude), !is.na(longitude))
  } else {
    message(SCRIPT_TAG, " Distance coordinate columns used: x=", x_col, ", y=", y_col, " (Euclidean meters).")
    out <- df %>%
      transmute(
        Site_raw = normalize_lookup_key(.data[[site_col]]),
        Site = map_site_labels(.data[[site_col]], lookup),
        sample_id = normalize_lookup_key(.data[[sample_col]]),
        latitude = NA_real_,
        longitude = NA_real_,
        x_m = coerce_numeric(.data[[x_col]]),
        y_m = coerce_numeric(.data[[y_col]]),
        coordinate_type = "xy",
        sampling_node = if (!is.na(node_col)) normalize_lookup_key(.data[[node_col]]) else NA_character_,
        pair_role = if (!is.na(role_col)) normalize_lookup_key(.data[[role_col]]) else NA_character_
      ) %>%
      filter(nzchar(Site), nzchar(sample_id), !is.na(x_m), !is.na(y_m))
  }
  
  validate_expected_sites(out$Site, "coordinate data")
  
  dupes <- out %>% count(Site, sample_id, name = "n") %>% filter(n > 1)
  if (nrow(dupes) > 0) {
    stop(
      SCRIPT_TAG, " Duplicate sample IDs within site in coordinate data. Examples: ",
      paste(head(paste(dupes$Site, dupes$sample_id, sep = ":"), 10), collapse = ", "),
      call. = FALSE
    )
  }
  
  out %>% order_by_expected_site()
}

compute_nearest_neighbors <- function(stems) {
  if (is.null(stems) || nrow(stems) == 0) return(NULL)
  
  out_rows <- list()
  row_i <- 1L
  
  for (site in EXPECTED_SITE_LABELS) {
    site_df <- stems %>% filter(Site == site) %>% arrange(sample_id)
    n <- nrow(site_df)
    if (n == 0) next
    
    for (i in seq_len(n)) {
      if (n < 2) {
        nearest_j <- NA_integer_
        nearest_dist <- NA_real_
      } else {
        dists <- vapply(seq_len(n), function(j) if (i == j) Inf else pair_distance_m(site_df, i, j), numeric(1))
        nearest_j <- which.min(dists)
        nearest_dist <- dists[nearest_j]
      }
      
      out_rows[[row_i]] <- data.frame(
        Site = site,
        sample_id = site_df$sample_id[i],
        nearest_neighbor_id = if (!is.na(nearest_j)) site_df$sample_id[nearest_j] else NA_character_,
        nearest_neighbor_distance_m = nearest_dist,
        distance_class = cut(
          nearest_dist,
          breaks = c(0, 1, 2, 5, Inf),
          labels = DISTANCE_CLASS_LEVELS,
          include.lowest = TRUE,
          right = FALSE
        ),
        coordinate_type = site_df$coordinate_type[i],
        sampling_node = site_df$sampling_node[i],
        pair_role = site_df$pair_role[i],
        stringsAsFactors = FALSE
      )
      row_i <- row_i + 1L
    }
  }
  
  bind_rows(out_rows) %>%
    mutate(distance_class = as.character(distance_class)) %>%
    order_by_expected_site()
}

summarise_distance_by_site <- function(nn_tbl) {
  class_counts <- nn_tbl %>%
    mutate(distance_class = factor(distance_class, levels = DISTANCE_CLASS_LEVELS)) %>%
    count(Site, distance_class, name = "n") %>%
    complete(Site = EXPECTED_SITE_LABELS, distance_class = DISTANCE_CLASS_LEVELS, fill = list(n = 0L)) %>%
    group_by(Site) %>%
    mutate(pct = format_pct(n, sum(n))) %>%
    ungroup() %>%
    pivot_wider(
      names_from = distance_class,
      values_from = c(n, pct),
      names_glue = "{.value}_{distance_class}"
    )
  
  stats_tbl <- nn_tbl %>%
    group_by(Site) %>%
    summarise(
      N_sampled_stems = n(),
      min_nearest_neighbor_distance_m = safe_min(nearest_neighbor_distance_m),
      mean_nearest_neighbor_distance_m = safe_mean(nearest_neighbor_distance_m),
      median_nearest_neighbor_distance_m = safe_median(nearest_neighbor_distance_m),
      max_nearest_neighbor_distance_m = safe_max(nearest_neighbor_distance_m),
      sd_nearest_neighbor_distance_m = safe_sd(nearest_neighbor_distance_m),
      .groups = "drop"
    )
  
  stats_tbl %>%
    left_join(class_counts, by = "Site") %>%
    order_by_expected_site()
}

summarise_distance_overall <- function(nn_tbl) {
  total <- nrow(nn_tbl)
  class_tbl <- nn_tbl %>%
    mutate(distance_class = factor(distance_class, levels = DISTANCE_CLASS_LEVELS)) %>%
    count(distance_class, name = "n") %>%
    complete(distance_class = DISTANCE_CLASS_LEVELS, fill = list(n = 0L)) %>%
    mutate(pct = format_pct(n, total))
  
  out <- data.frame(
    N_sampled_stems = total,
    min_nearest_neighbor_distance_m = safe_min(nn_tbl$nearest_neighbor_distance_m),
    mean_nearest_neighbor_distance_m = safe_mean(nn_tbl$nearest_neighbor_distance_m),
    median_nearest_neighbor_distance_m = safe_median(nn_tbl$nearest_neighbor_distance_m),
    max_nearest_neighbor_distance_m = safe_max(nn_tbl$nearest_neighbor_distance_m),
    sd_nearest_neighbor_distance_m = safe_sd(nn_tbl$nearest_neighbor_distance_m),
    stringsAsFactors = FALSE
  )
  
  for (cls in DISTANCE_CLASS_LEVELS) {
    n_val <- class_tbl$n[class_tbl$distance_class == cls]
    pct_val <- class_tbl$pct[class_tbl$distance_class == cls]
    nm <- gsub("[^A-Za-z0-9]+", "_", cls)
    out[[paste0("n_", nm)]] <- n_val
    out[[paste0("pct_", nm)]] <- pct_val
  }
  
  out
}

compute_intended_pairs <- function(stems) {
  if (is.null(stems) || nrow(stems) == 0) return(NULL)
  if (all(is.na(stems$sampling_node)) || all(!nzchar(stems$sampling_node))) {
    message(SCRIPT_TAG, " No sampling-node column detected; intended paired-node distances will not be calculated separately.")
    return(NULL)
  }
  
  pairs <- stems %>%
    filter(!is.na(sampling_node), nzchar(sampling_node)) %>%
    group_by(Site, sampling_node) %>%
    filter(n() >= 2) %>%
    group_modify(function(.x, .y) {
      if (nrow(.x) < 2) return(data.frame())
      comb <- utils::combn(seq_len(nrow(.x)), 2)
      rows <- lapply(seq_len(ncol(comb)), function(k) {
        i <- comb[1, k]
        j <- comb[2, k]
        data.frame(
          sample_id_1 = .x$sample_id[i],
          sample_id_2 = .x$sample_id[j],
          intended_pair_distance_m = pair_distance_m(.x, i, j),
          pair_role_1 = .x$pair_role[i],
          pair_role_2 = .x$pair_role[j],
          coordinate_type = .x$coordinate_type[i],
          stringsAsFactors = FALSE
        )
      })
      bind_rows(rows)
    }) %>%
    ungroup() %>%
    order_by_expected_site()
  
  if (nrow(pairs) == 0) {
    message(SCRIPT_TAG, " Sampling-node metadata were present, but no nodes contained >=2 sampled stems.")
    return(NULL)
  }
  
  pairs
}

make_distance_figures <- function(nn_tbl, out_paths) {
  hist_plot <- nn_tbl %>%
    filter(!is.na(nearest_neighbor_distance_m), nearest_neighbor_distance_m >= 0) %>%
    ggplot(aes(x = nearest_neighbor_distance_m)) +
    geom_histogram(binwidth = 0.5, boundary = 0, fill = "#2E8B57", color = "white") +
    theme_bw(base_size = 12) +
    labs(
      title = "Nearest-neighbor distances among sampled stems",
      x = "Nearest-neighbor distance (m)",
      y = "Number of sampled stems"
    )
  
  site_plot <- nn_tbl %>%
    filter(!is.na(nearest_neighbor_distance_m), nearest_neighbor_distance_m >= 0) %>%
    mutate(Site = factor(Site, levels = EXPECTED_SITE_LABELS)) %>%
    ggplot(aes(x = Site, y = nearest_neighbor_distance_m)) +
    geom_boxplot(fill = "#2E8B57", outlier.alpha = 0.6) +
    theme_bw(base_size = 12) +
    labs(
      title = "Nearest-neighbor distance by site",
      x = "Site",
      y = "Nearest-neighbor distance (m)"
    )
  
  ggsave(out_paths$histogram, hist_plot, width = 8, height = 5, dpi = 320)
  ggsave(out_paths$by_site, site_plot, width = 8, height = 5, dpi = 320)
}

# ------------------------------------------------------------------
# Diameter / developmental-stage analysis
# ------------------------------------------------------------------
dbh_candidates <- c("DBH", "dbh", "Dbh", "dbh_cm", "DBH_cm", "DHP", "dhp", "Dhp", "DHP_cm", "dhp_cm", "Dhp_cm", "Dhp_tige", "dhp_tige", "DHP_tige", "Dhp_tige_cm", "dhp_tige_cm")
basal_candidates <- c("basal_diameter", "Basal_diameter", "basal_diameter_cm", "diametre_basal", "diamètre_basal", "diametre_collet", "diamètre_collet", "diameter_cm", "diametre_cm", "diamètre_cm")
stage_candidates <- c("stage", "Stage", "developmental_stage", "Developmental_stage", "stade", "Stade", "classe_developpement", "classe_développement")

normalize_id_for_match <- function(x) {
  x <- normalize_lookup_key(x)
  x[x %in% c("", "NA", "N/A", "NaN", "NULL", "null", "na")] <- NA_character_
  x
}

normalize_tree_number <- function(x) {
  x <- normalize_lookup_key(x)
  x[x %in% c("", "NA", "NaN", "NULL")] <- NA_character_
  x_num <- suppressWarnings(as.numeric(gsub(",", ".", x, fixed = TRUE)))
  ifelse(!is.na(x_num) & abs(x_num - round(x_num)) < .Machine$double.eps^0.5, as.character(as.integer(round(x_num))), x)
}

compact_key <- function(...) {
  raw <- paste(..., sep = "_")
  out <- toupper(gsub("[^A-Za-z0-9]+", "", raw))
  out[is.na(out)] <- ""
  out
}

find_tree_col <- function(df) {
  pick_column(
    df,
    c("tree", "Tree", "tree_id", "tree_number", "arbre", "no_arbre", "numero_arbre", "Numéro_Arbre", "id_tige", "stem", "stem_id", "numero_individu", "sample_number", "sample_no", "numero", "num", "no"),
    required = FALSE
  )
}

build_genetic_stem_table <- function() {
  aligned <- align_df_ids_to_genind(gi, df_ids, context = "[sampled_stem_descriptive_results gi]")
  out <- data.frame(
    Individual = adegenet::indNames(gi),
    Site_raw = normalize_lookup_key(aligned$Site),
    stringsAsFactors = FALSE
  ) %>%
    bind_cols(aligned %>% select(-ind_id, -Site))
  out
}

choose_diameter_column <- function(df) {
  dbh_col <- pick_column(df, dbh_candidates, required = FALSE)
  basal_col <- pick_column(df, basal_candidates, required = FALSE)
  stage_col <- pick_column(df, stage_candidates, required = FALSE)
  
  if (!is.na(dbh_col)) {
    return(list(diameter_col = dbh_col, diameter_type = "DBH_cm", stage_col = stage_col, diagnostic = paste0("Using DBH column '", dbh_col, "'.")))
  }
  
  if (!is.na(basal_col)) {
    return(list(diameter_col = basal_col, diameter_type = "basal_diameter_cm", stage_col = stage_col, diagnostic = paste0("No DBH column detected; using basal diameter column '", basal_col, "'.")))
  }
  
  list(diameter_col = NA_character_, diameter_type = NA_character_, stage_col = stage_col, diagnostic = "No DBH or basal diameter column detected.")
}

find_diameter_source <- function(workbook, stem_tbl) {
  direct <- choose_diameter_column(stem_tbl)
  if (!is.na(direct$diameter_col)) {
    out <- stem_tbl
    out$Diameter_cm <- coerce_numeric(stem_tbl[[direct$diameter_col]])
    out$Diameter_variable <- direct$diameter_type
    out$Field_stage <- if (!is.na(direct$stage_col)) normalize_lookup_key(stem_tbl[[direct$stage_col]]) else NA_character_
    return(list(stem_tbl = out, source = paste0("df_ids column '", direct$diameter_col, "'"), diagnostic = direct$diagnostic))
  }
  
  if (exists("meta", inherits = TRUE) && is.data.frame(meta)) {
    meta_choice <- choose_diameter_column(meta)
    if (!is.na(meta_choice$diameter_col)) {
      id_col <- pick_column(meta, c(sample_candidates, "ID", "Id", "id"), required = FALSE)
      if (!is.na(id_col)) {
        source_ids <- normalize_id_for_match(meta[[id_col]])
        retained_ids <- normalize_id_for_match(stem_tbl$Individual)
        idx <- match(retained_ids, source_ids)
        if (sum(!is.na(idx)) / nrow(stem_tbl) >= 0.95) {
          out <- stem_tbl
          out$Diameter_cm <- ifelse(!is.na(idx), coerce_numeric(meta[[meta_choice$diameter_col]])[idx], NA_real_)
          out$Diameter_variable <- meta_choice$diameter_type
          out$Field_stage <- if (!is.na(meta_choice$stage_col)) normalize_lookup_key(meta[[meta_choice$stage_col]])[idx] else NA_character_
          return(list(stem_tbl = out, source = paste0("meta column '", meta_choice$diameter_col, "' matched by '", id_col, "'"), diagnostic = meta_choice$diagnostic))
        }
      }
    }
  }
  
  sheets <- readxl::excel_sheets(workbook)
  sheet_order <- unique(c(sheets[normalize_ascii_key(sheets) == normalize_ascii_key(GENETIC_SHEET_HINT)], sheets))
  
  for (sheet in sheet_order) {
    raw_df <- suppressMessages(readxl::read_excel(workbook, sheet = sheet)) %>% as.data.frame(stringsAsFactors = FALSE)
    choice <- choose_diameter_column(raw_df)
    if (is.na(choice$diameter_col)) next
    
    raw_diameter <- coerce_numeric(raw_df[[choice$diameter_col]])
    retained_ids <- normalize_id_for_match(stem_tbl$Individual)
    
    id_cols <- names(raw_df)[normalize_ascii_key(names(raw_df)) %in% normalize_ascii_key(c(sample_candidates, "ID", "Id", "id"))]
    for (id_col in id_cols) {
      raw_ids <- normalize_id_for_match(raw_df[[id_col]])
      idx <- match(retained_ids, raw_ids)
      if (sum(!is.na(idx)) / nrow(stem_tbl) >= 0.95 && !anyDuplicated(raw_ids[!is.na(raw_ids) & raw_ids %in% retained_ids])) {
        out <- stem_tbl
        out$Diameter_cm <- ifelse(!is.na(idx), raw_diameter[idx], NA_real_)
        out$Diameter_variable <- choice$diameter_type
        out$Field_stage <- if (!is.na(choice$stage_col)) normalize_lookup_key(raw_df[[choice$stage_col]])[idx] else NA_character_
        return(list(stem_tbl = out, source = paste0(workbook, " | sheet '", sheet, "' | column '", choice$diameter_col, "' matched by ID column '", id_col, "'"), diagnostic = choice$diagnostic))
      }
    }
    
    raw_site_col <- pick_column(raw_df, site_candidates, required = FALSE)
    raw_tree_col <- find_tree_col(raw_df)
    stem_tree_col <- find_tree_col(stem_tbl)
    if (!is.na(raw_site_col) && !is.na(raw_tree_col)) {
      stem_tree <- if (!is.na(stem_tree_col)) normalize_tree_number(stem_tbl[[stem_tree_col]]) else normalize_tree_number(gsub(".*?([0-9]+)[^0-9]*$", "\\1", stem_tbl$Individual))
      raw_key <- compact_key(normalize_lookup_key(raw_df[[raw_site_col]]), normalize_tree_number(raw_df[[raw_tree_col]]))
      stem_key <- compact_key(stem_tbl$Site_raw, stem_tree)
      idx <- match(stem_key, raw_key)
      if (sum(!is.na(idx)) / nrow(stem_tbl) >= 0.95 && !anyDuplicated(raw_key[!is.na(raw_key) & raw_key %in% stem_key])) {
        out <- stem_tbl
        out$Diameter_cm <- ifelse(!is.na(idx), raw_diameter[idx], NA_real_)
        out$Diameter_variable <- choice$diameter_type
        out$Field_stage <- if (!is.na(choice$stage_col)) normalize_lookup_key(raw_df[[choice$stage_col]])[idx] else NA_character_
        return(list(stem_tbl = out, source = paste0(workbook, " | sheet '", sheet, "' | column '", choice$diameter_col, "' matched by site + tree/stem number"), diagnostic = choice$diagnostic))
      }
    }
  }
  
  stop(
    SCRIPT_TAG, " Could not find or match a DBH/basal diameter column for >=95% of sampled genetic stems.",
    " Accepted DBH columns include: ", paste(dbh_candidates, collapse = ", "),
    "; accepted basal diameter columns include: ", paste(basal_candidates, collapse = ", "),
    call. = FALSE
  )
}

classify_diameter <- function(x) {
  out <- cut(
    x,
    breaks = c(-Inf, 1, 2, 5, 10, Inf),
    labels = DIAMETER_CLASS_LEVELS,
    right = FALSE
  )
  factor(out, levels = DIAMETER_CLASS_LEVELS, ordered = TRUE)
}

classify_developmental_stage <- function(diameter_cm, field_stage) {
  field_stage <- normalize_lookup_key(field_stage)
  if (any(nzchar(field_stage), na.rm = TRUE)) {
    out <- field_stage
    out[!nzchar(out)] <- NA_character_
    return(out)
  }
  
  # Conservative, reproducible fallback: do not claim true field stage when no
  # stage column exists. Use measured diameter classes as diameter-derived stages.
  as.character(classify_diameter(diameter_cm))
}

prepare_diameter_data <- function(workbook, lookup) {
  stem_tbl <- build_genetic_stem_table()
  resolved <- find_diameter_source(workbook, stem_tbl)
  out <- resolved$stem_tbl
  
  out <- out %>%
    mutate(
      Site = map_site_labels(Site_raw, lookup),
      Diameter_cm = coerce_numeric(Diameter_cm),
      Diameter_class = classify_diameter(Diameter_cm),
      Developmental_stage = classify_developmental_stage(Diameter_cm, Field_stage),
      Developmental_stage_basis = ifelse(any(nzchar(normalize_lookup_key(Field_stage)), na.rm = TRUE), "field developmental-stage column", paste0("measured ", unique(Diameter_variable)[1], " class"))
    ) %>%
    select(Individual, Site, Site_raw, Diameter_variable, Diameter_cm, Diameter_class, Field_stage, Developmental_stage, Developmental_stage_basis, everything()) %>%
    order_by_expected_site()
  
  validate_expected_sites(out$Site, "diameter/development data")
  
  message(SCRIPT_TAG, " Diameter source: ", resolved$source)
  message(SCRIPT_TAG, " Diameter diagnostic: ", resolved$diagnostic)
  message(SCRIPT_TAG, " Developmental-stage basis: ", unique(out$Developmental_stage_basis)[1])
  
  out
}

summarise_diameter_distribution <- function(stems) {
  stems %>%
    mutate(Diameter_class = factor(Diameter_class, levels = DIAMETER_CLASS_LEVELS, ordered = TRUE)) %>%
    count(Diameter_variable, Diameter_class, name = "N_stems") %>%
    complete(Diameter_variable, Diameter_class = DIAMETER_CLASS_LEVELS, fill = list(N_stems = 0L)) %>%
    group_by(Diameter_variable) %>%
    mutate(Percent_stems = format_pct(N_stems, sum(N_stems))) %>%
    ungroup() %>%
    arrange(Diameter_variable, Diameter_class)
}

summarise_development_overall <- function(stems) {
  total <- sum(!is.na(stems$Developmental_stage))
  stems %>%
    filter(!is.na(Developmental_stage), nzchar(Developmental_stage)) %>%
    count(Developmental_stage_basis, Developmental_stage, name = "N_stems") %>%
    mutate(Percent_stems = format_pct(N_stems, total)) %>%
    arrange(Developmental_stage)
}

summarise_development_by_site <- function(stems) {
  stage_counts <- stems %>%
    filter(!is.na(Developmental_stage), nzchar(Developmental_stage)) %>%
    count(Site, Developmental_stage_basis, Developmental_stage, name = "N_stems") %>%
    group_by(Site) %>%
    mutate(Percent_stems = format_pct(N_stems, sum(N_stems))) %>%
    ungroup()
  
  diameter_stats <- stems %>%
    group_by(Site, Diameter_variable) %>%
    summarise(
      N_sampled_stems = n(),
      N_with_diameter = sum(!is.na(Diameter_cm)),
      mean_diameter_cm = safe_mean(Diameter_cm),
      median_diameter_cm = safe_median(Diameter_cm),
      min_diameter_cm = safe_min(Diameter_cm),
      max_diameter_cm = safe_max(Diameter_cm),
      sd_diameter_cm = safe_sd(Diameter_cm),
      .groups = "drop"
    )
  
  diameter_class_wide <- stems %>%
    mutate(Diameter_class = factor(Diameter_class, levels = DIAMETER_CLASS_LEVELS, ordered = TRUE)) %>%
    count(Site, Diameter_class, name = "n") %>%
    complete(Site = EXPECTED_SITE_LABELS, Diameter_class = DIAMETER_CLASS_LEVELS, fill = list(n = 0L)) %>%
    group_by(Site) %>%
    mutate(pct = format_pct(n, sum(n))) %>%
    ungroup() %>%
    pivot_wider(names_from = Diameter_class, values_from = c(n, pct), names_glue = "{.value}_{Diameter_class}")
  
  stage_wide <- stage_counts %>%
    select(Site, Developmental_stage, N_stems, Percent_stems) %>%
    pivot_wider(names_from = Developmental_stage, values_from = c(N_stems, Percent_stems), values_fill = 0, names_glue = "{.value}_{Developmental_stage}")
  
  diameter_stats %>%
    left_join(diameter_class_wide, by = "Site") %>%
    left_join(stage_wide, by = "Site") %>%
    order_by_expected_site()
}

make_diameter_figures <- function(stems, out_paths) {
  diameter_plot <- stems %>%
    filter(!is.na(Diameter_class)) %>%
    mutate(Diameter_class = factor(Diameter_class, levels = DIAMETER_CLASS_LEVELS, ordered = TRUE)) %>%
    ggplot(aes(x = Diameter_class)) +
    geom_bar(fill = "#4C78A8") +
    theme_bw(base_size = 12) +
    labs(
      title = "Diameter distribution of sampled American beech stems",
      x = "Diameter class",
      y = "Number of sampled stems"
    )
  
  stage_plot <- stems %>%
    filter(!is.na(Developmental_stage), nzchar(Developmental_stage)) %>%
    mutate(Site = factor(Site, levels = EXPECTED_SITE_LABELS)) %>%
    ggplot(aes(x = Site, fill = Developmental_stage)) +
    geom_bar(position = "fill") +
    scale_y_continuous(labels = function(x) paste0(round(100 * x), "%")) +
    theme_bw(base_size = 12) +
    labs(
      title = "Diameter-derived developmental/stem-size classes by site",
      x = "Site",
      y = "Percent of sampled stems",
      fill = "Class"
    )
  
  ggsave(out_paths$diameter, diameter_plot, width = 8, height = 5, dpi = 320)
  ggsave(out_paths$stage_by_site, stage_plot, width = 9, height = 5, dpi = 320)
}

# ------------------------------------------------------------------
# Main workflow
# ------------------------------------------------------------------
output_paths <- character(0)
skipped_outputs <- character(0)

input_workbook <- find_input_workbook()
site_lookup <- load_site_lookup(input_workbook)

# Distances among sampled stems
coord_source <- read_coordinate_sheet(input_workbook)
coordinate_stems <- prepare_coordinate_data(coord_source, site_lookup)

if (!is.null(coordinate_stems) && nrow(coordinate_stems) > 0) {
  nearest_neighbor_pairs <- compute_nearest_neighbors(coordinate_stems)
  distance_by_site <- summarise_distance_by_site(nearest_neighbor_pairs)
  distance_summary <- summarise_distance_overall(nearest_neighbor_pairs)
  intended_pairs <- compute_intended_pairs(coordinate_stems)
  
  distance_summary_path <- file.path(TABLES_DIR, "sampled_stem_distance_summary.csv")
  distance_by_site_path <- file.path(TABLES_DIR, "sampled_stem_distance_by_site.csv")
  nearest_pairs_path <- file.path(TABLES_DIR, "sampled_stem_nearest_neighbor_pairs.csv")
  intended_pairs_path <- file.path(TABLES_DIR, "sampled_stem_intended_pair_distances.csv")
  
  output_paths <- c(
    output_paths,
    write_csv_msg(distance_summary, distance_summary_path),
    write_csv_msg(distance_by_site, distance_by_site_path),
    write_csv_msg(nearest_neighbor_pairs, nearest_pairs_path)
  )
  
  if (!is.null(intended_pairs) && nrow(intended_pairs) > 0) {
    output_paths <- c(output_paths, write_csv_msg(intended_pairs, intended_pairs_path))
  } else {
    skipped_outputs <- c(skipped_outputs, intended_pairs_path)
  }
  
  distance_figure_paths <- list(
    histogram = file.path(FIGURES_DIR, "sampled_stem_nearest_neighbor_distance_histogram.png"),
    by_site = file.path(FIGURES_DIR, "sampled_stem_nearest_neighbor_distance_by_site.png")
  )
  make_distance_figures(nearest_neighbor_pairs, distance_figure_paths)
  output_paths <- c(output_paths, unlist(distance_figure_paths))
} else {
  skipped_outputs <- c(
    skipped_outputs,
    file.path(TABLES_DIR, "sampled_stem_distance_summary.csv"),
    file.path(TABLES_DIR, "sampled_stem_distance_by_site.csv"),
    file.path(TABLES_DIR, "sampled_stem_nearest_neighbor_pairs.csv"),
    file.path(FIGURES_DIR, "sampled_stem_nearest_neighbor_distance_histogram.png"),
    file.path(FIGURES_DIR, "sampled_stem_nearest_neighbor_distance_by_site.png")
  )
}

# Diameter and developmental/stem-size summaries
diameter_stems <- prepare_diameter_data(input_workbook, site_lookup)
dbh_distribution <- summarise_diameter_distribution(diameter_stems)
development_summary <- summarise_development_overall(diameter_stems)
development_by_site <- summarise_development_by_site(diameter_stems)

# Include the individual-level stem table because it is useful for audit/reanalysis.
stem_diameter_values <- diameter_stems %>%
  select(Individual, Site, Site_raw, Diameter_variable, Diameter_cm, Diameter_class, Field_stage, Developmental_stage, Developmental_stage_basis)

# Requested table paths.
dev_summary_path <- file.path(TABLES_DIR, "sampled_stem_developmental_stage_summary.csv")
dev_by_site_path <- file.path(TABLES_DIR, "sampled_stem_developmental_stage_by_site.csv")
dbh_distribution_path <- file.path(TABLES_DIR, "sampled_stem_dbh_distribution.csv")
stem_values_path <- file.path(TABLES_DIR, "sampled_stem_diameter_values.csv")

output_paths <- c(
  output_paths,
  write_csv_msg(development_summary, dev_summary_path),
  write_csv_msg(development_by_site, dev_by_site_path),
  write_csv_msg(dbh_distribution, dbh_distribution_path),
  write_csv_msg(stem_diameter_values, stem_values_path)
)

diameter_figure_paths <- list(
  diameter = file.path(FIGURES_DIR, "sampled_stem_diameter_distribution.png"),
  stage_by_site = file.path(FIGURES_DIR, "sampled_stem_developmental_stage_by_site.png")
)
make_diameter_figures(diameter_stems, diameter_figure_paths)
output_paths <- c(output_paths, unlist(diameter_figure_paths))

# ------------------------------------------------------------------
# Console summary
# ------------------------------------------------------------------
cat("\n", SCRIPT_TAG, " Article-ready sampled-stem descriptive outputs complete.\n", sep = "")
cat(SCRIPT_TAG, " Input workbook: ", input_workbook, "\n", sep = "")

if (exists("nearest_neighbor_pairs") && nrow(nearest_neighbor_pairs) > 0) {
  cat(SCRIPT_TAG, " Nearest-neighbor distances: N stems = ", nrow(nearest_neighbor_pairs),
      "; median = ", sprintf("%.3f", safe_median(nearest_neighbor_pairs$nearest_neighbor_distance_m)),
      " m; mean = ", sprintf("%.3f", safe_mean(nearest_neighbor_pairs$nearest_neighbor_distance_m)),
      " m; range = ", sprintf("%.3f", safe_min(nearest_neighbor_pairs$nearest_neighbor_distance_m)),
      "-", sprintf("%.3f", safe_max(nearest_neighbor_pairs$nearest_neighbor_distance_m)), " m.\n", sep = "")
} else {
  cat(SCRIPT_TAG, " Nearest-neighbor distances: skipped because no usable coordinate sheet was detected.\n", sep = "")
}

cat(SCRIPT_TAG, " Diameter data: N stems = ", nrow(diameter_stems),
    "; N with diameter = ", sum(!is.na(diameter_stems$Diameter_cm)),
    "; median = ", sprintf("%.3f", safe_median(diameter_stems$Diameter_cm)),
    " cm; mean = ", sprintf("%.3f", safe_mean(diameter_stems$Diameter_cm)),
    " cm; range = ", sprintf("%.3f", safe_min(diameter_stems$Diameter_cm)),
    "-", sprintf("%.3f", safe_max(diameter_stems$Diameter_cm)), " cm.\n", sep = "")

cat(SCRIPT_TAG, " Outputs created:\n", sep = "")
for (path in output_paths) cat("  - ", path, "\n", sep = "")

if (length(skipped_outputs) > 0) {
  cat(SCRIPT_TAG, " Outputs skipped because required metadata were not available:\n", sep = "")
  for (path in skipped_outputs) cat("  - ", path, "\n", sep = "")
}

cat(SCRIPT_TAG, " Send these files back to ChatGPT for writing the Results paragraph:\n", sep = "")
files_to_send <- c(
  file.path(TABLES_DIR, "sampled_stem_distance_summary.csv"),
  file.path(TABLES_DIR, "sampled_stem_distance_by_site.csv"),
  file.path(TABLES_DIR, "sampled_stem_nearest_neighbor_pairs.csv"),
  file.path(TABLES_DIR, "sampled_stem_developmental_stage_summary.csv"),
  file.path(TABLES_DIR, "sampled_stem_developmental_stage_by_site.csv"),
  file.path(TABLES_DIR, "sampled_stem_dbh_distribution.csv"),
  file.path(TABLES_DIR, "sampled_stem_diameter_values.csv")
)
for (path in files_to_send) cat("  - ", path, "\n", sep = "")