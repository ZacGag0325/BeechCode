# scripts/11_isolation_by_distance.R
############################################################
# Isolation by Distance (Mantel test)
# Genetic distance input: pairwise Jost's D (clone-corrected, site-level)
#
# Preconditions for a valid Mantel / IBD analysis in this workflow:
# 1) A site-by-site genetic distance matrix already exists
#    (here: outputs/matrices/pairwise_jostD.csv from gi_mll).
# 2) Exactly one valid latitude/longitude pair can be assigned to each site.
# 3) Site names in the genetic matrix and coordinate table either match
#    directly or can be safely normalized and matched.
#
# Outputs:
# - outputs/tables/mantel_test_results.csv
# - outputs/tables/site_coordinate_audit.csv
# - outputs/matrices/geographic_distance_matrix.csv
# - outputs/tables/ibd_points.csv
# - outputs/figures/isolation_by_distance.jpeg
############################################################

suppressPackageStartupMessages({
  library(vegan)
  library(geosphere)
  library(ggplot2)
  library(dplyr)
})

source("scripts/_load_objects.R")

cat("[11_isolation_by_distance] Running Mantel test: Jost's D vs geographic distance...\n")
cat("[11_isolation_by_distance] Genetic distance matrix expected from scripts/06_distance_matrices.R using gi_mll.\n")

normalize_site <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("\\uFEFF", "", x, fixed = TRUE)
  x <- gsub("[[:cntrl:]]", "", x)
  x <- gsub("\\s+", " ", x)
  toupper(x)
}

normalize_name <- function(x) {
  x <- iconv(as.character(x), from = "", to = "ASCII//TRANSLIT")
  x[is.na(x)] <- as.character(x[is.na(x)])
  x <- tolower(trimws(x))
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

is_numeric_like <- function(x) {
  x <- trimws(as.character(x))
  x[x == ""] <- NA_character_
  x <- gsub(",", ".", x, fixed = TRUE)
  num <- suppressWarnings(as.numeric(x))
  sum(!is.na(num)) > 0 && all(is.na(x) | !is.na(num))
}

coerce_numeric_safely <- function(x) {
  x <- trimws(as.character(x))
  x[x == ""] <- NA_character_
  x <- gsub(",", ".", x, fixed = TRUE)
  suppressWarnings(as.numeric(x))
}

safe_read_square_matrix <- function(path) {
  if (!file.exists(path)) {
    stop("Missing file: ", path)
  }
  
  cat("[11_isolation_by_distance] Loading genetic matrix from: ", path, "\n", sep = "")
  
  df <- tryCatch(
    read.csv(path, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE, strip.white = TRUE),
    error = function(e) NULL
  )
  
  if (is.null(df)) {
    raw_df <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, strip.white = TRUE)
    if (ncol(raw_df) < 2) {
      stop("Imported matrix has fewer than 2 columns: ", path)
    }
    row_id <- trimws(as.character(raw_df[[1]]))
    df <- raw_df[, -1, drop = FALSE]
    rownames(df) <- row_id
  }
  
  row_id <- trimws(rownames(df))
  colnames(df) <- trimws(colnames(df))
  
  keep <- nzchar(row_id)
  if (!all(keep)) {
    cat("[11_isolation_by_distance] Dropping ", sum(!keep), " rows with empty site names from imported matrix.\n", sep = "")
    df <- df[keep, , drop = FALSE]
    row_id <- row_id[keep]
  }
  rownames(df) <- row_id
  
  if (ncol(df) > 0) {
    first_col_chr <- trimws(as.character(df[[1]]))
    if (all(first_col_chr == rownames(df))) {
      cat("[11_isolation_by_distance] Dropping duplicated site-name column from imported matrix.\n")
      df <- df[, -1, drop = FALSE]
    }
  }
  
  if (anyDuplicated(rownames(df))) {
    dups <- unique(rownames(df)[duplicated(rownames(df))])
    stop("Duplicate row site names in genetic matrix: ", paste(dups, collapse = ", "))
  }
  
  bad_cols <- character(0)
  for (j in seq_len(ncol(df))) {
    x <- trimws(as.character(df[[j]]))
    x[x == ""] <- NA_character_
    if (any(grepl("^[+-]?[0-9]+,[0-9]+$", x[!is.na(x)]))) {
      x <- gsub(",", ".", x, fixed = TRUE)
    }
    num <- suppressWarnings(as.numeric(x))
    bad <- !is.na(x) & is.na(num)
    if (any(bad)) {
      bad_cols <- c(bad_cols, colnames(df)[j])
    }
    df[[j]] <- num
  }
  
  if (length(bad_cols) > 0) {
    stop(
      "Imported genetic matrix has non-numeric values in column(s): ",
      paste(unique(bad_cols), collapse = ", "),
      " (", path, ")"
    )
  }
  
  mat <- as.matrix(df)
  storage.mode(mat) <- "numeric"
  
  if (nrow(mat) != ncol(mat)) {
    stop("Imported genetic matrix is not square: ", nrow(mat), "x", ncol(mat), " (", path, ")")
  }
  if (!setequal(rownames(mat), colnames(mat))) {
    stop("Imported matrix row/column names do not match as sets: ", path)
  }
  
  mat <- mat[, rownames(mat), drop = FALSE]
  
  if (!identical(rownames(mat), colnames(mat))) {
    stop("Imported matrix row/column names could not be aligned identically: ", path)
  }
  
  diag_idx <- row(mat) == col(mat)
  finite_mask <- is.finite(mat)
  n_nonfinite_diag <- sum(!finite_mask[diag_idx], na.rm = TRUE)
  n_nonfinite_offdiag <- sum(!finite_mask[!diag_idx], na.rm = TRUE)
  
  cat("[11_isolation_by_distance] Non-finite entries on diagonal: ", n_nonfinite_diag, "\n", sep = "")
  cat("[11_isolation_by_distance] Non-finite entries off-diagonal: ", n_nonfinite_offdiag, "\n", sep = "")
  
  if (n_nonfinite_offdiag > 0) {
    stop("Imported genetic matrix contains non-finite OFF-diagonal values; cannot run Mantel.")
  }
  if (n_nonfinite_diag > 0) {
    mat[diag_idx & !finite_mask] <- 0
    cat("[11_isolation_by_distance] Replaced non-finite diagonal values with 0 for distance-matrix use.\n")
  }
  
  cat("[11_isolation_by_distance] Genetic matrix dimensions: ", nrow(mat), " x ", ncol(mat), "\n", sep = "")
  cat("[11_isolation_by_distance] Genetic matrix sites (first 10): ", paste(head(rownames(mat), 10), collapse = ", "), "\n", sep = "")
  
  mat
}

score_site_column <- function(df) {
  nms <- names(df)
  std <- normalize_name(nms)
  scores <- rep(0, length(nms))
  
  scores <- scores + ifelse(std %in% c("site", "population", "pop", "site_name", "population_name"), 100, 0)
  scores <- scores + ifelse(grepl("(^|_)site($|_)|(^|_)population($|_)|(^|_)pop($|_)", std), 40, 0)
  scores <- scores - ifelse(grepl("sample|individual|specimen|nom_labo|echant|barcode|clone|mll|mlg|id", std), 100, 0)
  
  for (i in seq_along(nms)) {
    values <- trimws(as.character(df[[i]]))
    values <- values[nzchar(values) & !is.na(values)]
    if (length(values) == 0) {
      scores[i] <- scores[i] - 50
      next
    }
    uniq_prop <- length(unique(values)) / length(values)
    if (uniq_prop >= 0.95) {
      scores[i] <- scores[i] - 40
    }
  }
  
  best_idx <- which.max(scores)
  if (length(best_idx) == 0 || scores[best_idx] <= 0) {
    return(NA_character_)
  }
  nms[best_idx]
}

score_coordinate_column <- function(df, type = c("lat", "lon")) {
  type <- match.arg(type)
  nms <- names(df)
  std <- normalize_name(nms)
  scores <- rep(-Inf, length(nms))
  
  for (i in seq_along(nms)) {
    nm <- std[i]
    raw <- df[[i]]
    num <- coerce_numeric_safely(raw)
    non_missing <- !is.na(num)
    n_non_missing <- sum(non_missing)
    
    if (n_non_missing == 0) {
      next
    }
    if (!is_numeric_like(raw)) {
      next
    }
    if (grepl("sample|individual|specimen|nom_labo|echant|barcode|clone|mll|mlg|id", nm)) {
      next
    }
    
    val <- num[non_missing]
    if (type == "lat") {
      if (any(val < -90 | val > 90)) {
        next
      }
      score <- 0
      if (nm %in% c("latitude", "lat", "decimal_latitude", "lat_dd")) score <- score + 100
      if (grepl("latitude|(^|_)lat($|_)", nm)) score <- score + 50
      if (grepl("longitude|(^|_)lon($|_)|(^|_)long($|_)", nm)) score <- score - 120
    } else {
      if (any(val < -180 | val > 180)) {
        next
      }
      score <- 0
      if (nm %in% c("longitude", "lon", "long", "decimal_longitude", "decimal_long", "lon_dd", "long_dd")) score <- score + 100
      if (grepl("longitude|(^|_)lon($|_)|(^|_)long($|_)", nm)) score <- score + 50
      if (grepl("latitude|(^|_)lat($|_)", nm)) score <- score - 120
    }
    
    if (length(unique(round(val, 4))) < 2) score <- score - 10
    if (mean(is.na(num)) > 0.5) score <- score - 15
    scores[i] <- score
  }
  
  best_idx <- which.max(scores)
  if (length(best_idx) == 0 || !is.finite(scores[best_idx]) || scores[best_idx] <= 0) {
    return(NA_character_)
  }
  nms[best_idx]
}

collapse_site_coordinates <- function(meta_df, site_col, lat_col, lon_col, source_name) {
  coords_raw <- data.frame(
    Site = trimws(as.character(meta_df[[site_col]])),
    Latitude = coerce_numeric_safely(meta_df[[lat_col]]),
    Longitude = coerce_numeric_safely(meta_df[[lon_col]]),
    stringsAsFactors = FALSE
  )
  
  coords_raw$Site_norm <- normalize_site(coords_raw$Site)
  coords_raw <- coords_raw %>%
    filter(nzchar(Site), nzchar(Site_norm), !is.na(Latitude), !is.na(Longitude))
  
  if (nrow(coords_raw) == 0) {
    stop(
      "No usable site coordinates remained after filtering in ", source_name,
      ". Required columns must contain non-missing numeric latitude and longitude values."
    )
  }
  
  if (any(coords_raw$Latitude < -90 | coords_raw$Latitude > 90)) {
    bad <- unique(coords_raw$Site[coords_raw$Latitude < -90 | coords_raw$Latitude > 90])
    stop("Latitude values outside [-90, 90] detected for site(s): ", paste(head(bad, 20), collapse = ", "))
  }
  if (any(coords_raw$Longitude < -180 | coords_raw$Longitude > 180)) {
    bad <- unique(coords_raw$Site[coords_raw$Longitude < -180 | coords_raw$Longitude > 180])
    stop("Longitude values outside [-180, 180] detected for site(s): ", paste(head(bad, 20), collapse = ", "))
  }
  
  coord_audit <- coords_raw %>%
    group_by(Site_norm) %>%
    summarise(
      Site = dplyr::first(Site),
      n_rows = n(),
      n_unique_lat = n_distinct(round(Latitude, 8)),
      n_unique_lon = n_distinct(round(Longitude, 8)),
      Latitude_mean = mean(Latitude, na.rm = TRUE),
      Longitude_mean = mean(Longitude, na.rm = TRUE),
      Latitude_min = min(Latitude, na.rm = TRUE),
      Latitude_max = max(Latitude, na.rm = TRUE),
      Longitude_min = min(Longitude, na.rm = TRUE),
      Longitude_max = max(Longitude, na.rm = TRUE),
      inconsistent_within_site = n_unique_lat > 1 | n_unique_lon > 1,
      .groups = "drop"
    ) %>%
    arrange(Site)
  
  inconsistent_sites <- coord_audit$Site[coord_audit$inconsistent_within_site]
  if (length(inconsistent_sites) > 0) {
    warning(
      "[11_isolation_by_distance] Inconsistent duplicate coordinates detected within site(s): ",
      paste(inconsistent_sites, collapse = ", "),
      ". Using the mean latitude/longitude per site as the site centroid for Mantel/IBD."
    )
  }
  
  coords <- coord_audit %>%
    transmute(
      Site = Site,
      Site_norm = Site_norm,
      Latitude = Latitude_mean,
      Longitude = Longitude_mean,
      n_rows_collapsed = n_rows,
      inconsistent_within_site = inconsistent_within_site
    )
  
  list(coords = coords, audit = coord_audit)
}

candidate_coordinate_files <- function() {
  candidates <- c(
    file.path(PROJECT_ROOT, "inputs", "site_metadata.csv"),
    file.path(PROJECT_ROOT, "inputs", "site_coordinates.csv"),
    file.path(PROJECT_ROOT, "inputs", "site_coords.csv"),
    file.path(PROJECT_ROOT, "data", "raw", "site_metadata.csv"),
    file.path(PROJECT_ROOT, "data", "raw", "site_coordinates.csv"),
    file.path(PROJECT_ROOT, "data", "raw", "site_coords.csv")
  )
  unique(candidates[file.exists(candidates)])
}

read_fallback_coordinate_table <- function() {
  files <- candidate_coordinate_files()
  if (length(files) == 0) {
    stop(
      "Could not extract clean site coordinates from meta, and no dedicated site-coordinate table was found.\n",
      "Provide one of these files: inputs/site_metadata.csv, inputs/site_coordinates.csv, inputs/site_coords.csv, data/raw/site_metadata.csv, data/raw/site_coordinates.csv, or data/raw/site_coords.csv.\n",
      "The table must contain one site column (e.g. Site) plus numeric latitude and longitude columns (e.g. Latitude and Longitude)."
    )
  }
  
  path <- files[1]
  cat("[11_isolation_by_distance] Falling back to dedicated coordinate table: ", path, "\n", sep = "")
  df <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  list(data = df, source_name = path)
}

build_site_coordinates <- function(meta_df) {
  if (!is.data.frame(meta_df) || ncol(meta_df) == 0) {
    stop("[11_isolation_by_distance] meta must be a non-empty data.frame.")
  }
  
  attempt_df <- meta_df
  attempt_source <- "meta object (outputs/v1/objects/meta.rds)"
  
  site_col <- score_site_column(attempt_df)
  lat_col <- score_coordinate_column(attempt_df, "lat")
  lon_col <- score_coordinate_column(attempt_df, "lon")
  
  if (is.na(site_col) || is.na(lat_col) || is.na(lon_col)) {
    cat("[11_isolation_by_distance] Clean coordinates not detected in meta; trying dedicated coordinate table fallback.\n")
    fallback <- read_fallback_coordinate_table()
    attempt_df <- fallback$data
    attempt_source <- fallback$source_name
    site_col <- score_site_column(attempt_df)
    lat_col <- score_coordinate_column(attempt_df, "lat")
    lon_col <- score_coordinate_column(attempt_df, "lon")
  }
  
  if (is.na(site_col) || is.na(lat_col) || is.na(lon_col)) {
    stop(
      "[11_isolation_by_distance] Could not identify a valid Site / Latitude / Longitude column set in ", attempt_source, ".\n",
      "Detected columns were insufficient or unsafe. A valid table must contain:\n",
      "- one site column (e.g. Site or Population)\n",
      "- one numeric latitude column with values in [-90, 90]\n",
      "- one numeric longitude column with values in [-180, 180]\n",
      "Sample-ID columns (e.g. Nom_Labo_Échantillons) are explicitly rejected as coordinate columns."
    )
  }
  
  cat("[11_isolation_by_distance] Detected Site column: ", site_col, "\n", sep = "")
  cat("[11_isolation_by_distance] Detected Latitude column: ", lat_col, "\n", sep = "")
  cat("[11_isolation_by_distance] Detected Longitude column: ", lon_col, "\n", sep = "")
  cat("[11_isolation_by_distance] Coordinate source: ", attempt_source, "\n", sep = "")
  
  collapsed <- collapse_site_coordinates(
    meta_df = attempt_df,
    site_col = site_col,
    lat_col = lat_col,
    lon_col = lon_col,
    source_name = attempt_source
  )
  
  if (nrow(collapsed$coords) < 3) {
    stop("[11_isolation_by_distance] Fewer than 3 sites have valid coordinates after collapse; Mantel test cannot be run.")
  }
  
  list(
    coords = collapsed$coords,
    audit = collapsed$audit,
    source_name = attempt_source,
    site_col = site_col,
    lat_col = lat_col,
    lon_col = lon_col
  )
}

# ----------------------------
# 1) Load genetic distance matrix (Jost's D)
# ----------------------------
jost_file <- file.path(MATRICES_DIR, "pairwise_jostD.csv")
pairwise_jostD <- safe_read_square_matrix(jost_file)

# ----------------------------
# 2) Load and validate site coordinates
# ----------------------------
coord_build <- build_site_coordinates(meta)
site_meta <- coord_build$coords
coord_audit <- coord_build$audit

audit_file <- file.path(TABLES_DIR, "site_coordinate_audit.csv")
write.csv(coord_audit, audit_file, row.names = FALSE)
cat("[11_isolation_by_distance] Saved coordinate audit: ", audit_file, "\n", sep = "")

# ----------------------------
# 3) Align shared sites and build geographic distance matrix
# ----------------------------
gen_sites_raw <- trimws(rownames(pairwise_jostD))
col_sites_raw <- trimws(colnames(pairwise_jostD))
gen_sites_norm <- normalize_site(gen_sites_raw)
col_sites_norm <- normalize_site(col_sites_raw)
meta_sites_norm <- normalize_site(site_meta$Site)

if (!identical(gen_sites_norm, col_sites_norm)) {
  stop("Genetic matrix row/column order mismatch after normalization; cannot build valid distance object.")
}

cat("[11_isolation_by_distance] Genetic sites before matching (first 10): ", paste(head(gen_sites_raw, 10), collapse = ", "), "\n", sep = "")
cat("[11_isolation_by_distance] Metadata sites before matching (first 10): ", paste(head(site_meta$Site, 10), collapse = ", "), "\n", sep = "")

shared_norm <- intersect(gen_sites_norm, meta_sites_norm)
cat("[11_isolation_by_distance] Shared site count after normalization: ", length(shared_norm), "\n", sep = "")

if (length(setdiff(gen_sites_norm, shared_norm)) > 0) {
  cat("[11_isolation_by_distance] Dropping genetic-only normalized sites: ",
      paste(setdiff(gen_sites_norm, shared_norm), collapse = ", "), "\n", sep = "")
}
if (length(setdiff(meta_sites_norm, shared_norm)) > 0) {
  cat("[11_isolation_by_distance] Dropping coordinate-only normalized sites: ",
      paste(setdiff(meta_sites_norm, shared_norm), collapse = ", "), "\n", sep = "")
}

if (length(shared_norm) < 3) {
  stop(
    "Need at least 3 shared sites between pairwise_jostD.csv and the validated coordinate table. Found: ",
    length(shared_norm),
    "\nGenetic normalized sites: ", paste(gen_sites_norm, collapse = ", "),
    "\nCoordinate normalized sites: ", paste(meta_sites_norm, collapse = ", ")
  )
}

meta_key <- setNames(site_meta$Site, meta_sites_norm)
shared_display <- unname(meta_key[shared_norm])
shared_display[is.na(shared_display)] <- shared_norm[is.na(shared_display)]

idx_gen <- match(shared_norm, gen_sites_norm)
pairwise_jostD <- pairwise_jostD[idx_gen, idx_gen, drop = FALSE]
rownames(pairwise_jostD) <- shared_display
colnames(pairwise_jostD) <- shared_display

idx_meta <- match(shared_norm, meta_sites_norm)
site_meta <- site_meta[idx_meta, , drop = FALSE]
site_meta$Site <- shared_display

if (anyDuplicated(site_meta$Site)) {
  stop("[11_isolation_by_distance] Duplicate site labels remain after alignment: ",
       paste(unique(site_meta$Site[duplicated(site_meta$Site)]), collapse = ", "))
}

coords <- as.matrix(site_meta[, c("Longitude", "Latitude")])
if (!is.matrix(coords) || nrow(coords) != nrow(site_meta)) {
  stop("[11_isolation_by_distance] Failed to build coordinate matrix for geographic distances.")
}
if (any(!is.finite(coords))) {
  stop("[11_isolation_by_distance] Coordinate matrix contains non-finite values after validation.")
}

geographic_km <- geosphere::distm(coords, fun = geosphere::distHaversine) / 1000
rownames(geographic_km) <- shared_display
colnames(geographic_km) <- shared_display

geo_file <- file.path(MATRICES_DIR, "geographic_distance_matrix.csv")
write.csv(geographic_km, geo_file, row.names = TRUE)
cat("[11_isolation_by_distance] Saved geographic matrix: ", geo_file, "\n", sep = "")

# ----------------------------
# 4) Mantel test
# ----------------------------
pairwise_jostD <- as.matrix(pairwise_jostD)
geographic_km <- as.matrix(geographic_km)
storage.mode(pairwise_jostD) <- "numeric"
storage.mode(geographic_km) <- "numeric"

if (!identical(rownames(pairwise_jostD), rownames(geographic_km))) {
  stop("Site-order mismatch between genetic and geographic matrices before Mantel.")
}
if (any(!is.finite(pairwise_jostD[row(pairwise_jostD) != col(pairwise_jostD)]))) {
  stop("Genetic matrix has non-finite off-diagonal values before Mantel.")
}
if (any(!is.finite(geographic_km[row(geographic_km) != col(geographic_km)]))) {
  stop("Geographic matrix has non-finite off-diagonal values before Mantel.")
}

diag(pairwise_jostD) <- 0
diag(geographic_km) <- 0

gen_dist <- as.dist(pairwise_jostD)
geo_dist <- as.dist(geographic_km)

if (!identical(attr(gen_dist, "Labels"), attr(geo_dist, "Labels"))) {
  stop("Mantel labels mismatch between genetic and geographic distance objects after alignment.")
}

mantel_fit <- vegan::mantel(
  gen_dist,
  geo_dist,
  method = "pearson",
  permutations = 9999
)

mantel_results <- data.frame(
  test = "Mantel_JostD_vs_Geographic",
  statistic_r = as.numeric(mantel_fit$statistic),
  p_value = as.numeric(mantel_fit$signif),
  permutations = as.integer(mantel_fit$permutations),
  n_shared_sites = length(shared_display),
  genetic_distance_input = "pairwise_jostD.csv",
  coordinate_source = coord_build$source_name,
  site_column = coord_build$site_col,
  latitude_column = coord_build$lat_col,
  longitude_column = coord_build$lon_col,
  stringsAsFactors = FALSE
)

mantel_file <- file.path(TABLES_DIR, "mantel_test_results.csv")
write.csv(mantel_results, mantel_file, row.names = FALSE)
cat("[11_isolation_by_distance] Saved Mantel results: ", mantel_file, "\n", sep = "")

# ----------------------------
# 5) Export IBD points and plot
# ----------------------------
upper_idx <- upper.tri(pairwise_jostD)

ibd_points <- data.frame(
  Site1 = rownames(pairwise_jostD)[row(pairwise_jostD)[upper_idx]],
  Site2 = colnames(pairwise_jostD)[col(pairwise_jostD)[upper_idx]],
  Geographic_km = as.numeric(geographic_km[upper_idx]),
  JostD = as.numeric(pairwise_jostD[upper_idx]),
  stringsAsFactors = FALSE
)

ibd_points_file <- file.path(TABLES_DIR, "ibd_points.csv")
write.csv(ibd_points, ibd_points_file, row.names = FALSE)
cat("[11_isolation_by_distance] Saved IBD points: ", ibd_points_file, "\n", sep = "")

ann_text <- sprintf(
  "Mantel r = %.3f\np = %.4f\nSites = %d",
  as.numeric(mantel_fit$statistic),
  as.numeric(mantel_fit$signif),
  length(shared_display)
)

ibd_plot <- ggplot(ibd_points, aes(Geographic_km, JostD)) +
  geom_point(size = 2.2, alpha = 0.8, color = "#2C3E50") +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
  annotate("text", x = Inf, y = Inf, label = ann_text, hjust = 1.05, vjust = 1.3, size = 4.1) +
  theme_bw(base_size = 12) +
  labs(
    title = "Isolation by distance",
    subtitle = "Clone-corrected pairwise Jost's D by site versus geographic distance",
    x = "Geographic distance (km)",
    y = "Pairwise Jost's D"
  )

ibd_plot_file <- file.path(FIGURES_DIR, "isolation_by_distance.jpeg")
ggsave(ibd_plot_file, plot = ibd_plot, width = 8, height = 6, dpi = 320)
cat("[11_isolation_by_distance] Saved figure: ", ibd_plot_file, "\n", sep = "")