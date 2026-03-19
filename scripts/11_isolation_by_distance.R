# scripts/11_isolation_by_distance.R
############################################################
# Isolation by Distance (Mantel test)
# Genetic distance input: pairwise Jost's D (clone-corrected, site-level)
#
# Biological rationale:
# - The primary IBD analysis in this pipeline is site-level because the main
#   differentiation summaries (Jost's D, FST, AMOVA) are all interpreted among
#   sampling localities, not among ramets.
# - Therefore, the primary Mantel test compares a site-by-site genetic
#   distance matrix (pairwise Jost's D from gi_mll) against a site-by-site
#   geographic distance matrix built from mean site centroids.
# - A secondary OPTIONAL individual-level IBD analysis is added here only as a
#   complementary diagnostic. It does not replace the main site-level result.
#
# Interpretation note:
# - Site-level IBD asks whether genetically differentiated sites tend to be
#   farther apart geographically.
# - Individual-level IBD, when available, asks whether pairwise genetic
#   differences among sampled individuals increase with pairwise geographic
#   distance. Because that mixes within- and among-site structure, it should be
#   interpreted as a secondary diagnostic rather than the main inference.
#
# Preconditions for the primary site-level Mantel / IBD analysis:
# 1) A site-by-site genetic distance matrix already exists
#    (outputs/matrices/pairwise_jostD.csv from gi_mll).
# 2) Exactly one valid centroid can be assigned to each site after collapsing
#    individual/site metadata to mean latitude and longitude.
# 3) Site names in the genetic matrix and coordinate table either match
#    directly or can be safely normalized and aligned.
#
# Outputs:
# - outputs/tables/mantel_test_results.csv
# - outputs/tables/site_coordinate_audit.csv
# - outputs/tables/site_coordinates_clean.csv
# - outputs/matrices/geographic_distance_matrix.csv
# - outputs/tables/ibd_points.csv
# - outputs/figures/isolation_by_distance.jpeg
# - outputs/tables/individual_coordinate_audit.csv (optional)
# - outputs/tables/individual_mantel_test_results.csv (optional)
# - outputs/tables/individual_ibd_points.csv (optional)
# - outputs/matrices/individual_geographic_distance_matrix.csv (optional)
# - outputs/matrices/individual_genetic_distance_matrix.csv (optional)
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
cat("[11_isolation_by_distance] Primary interpretation = among-site IBD using mean site centroids.\n")
cat("[11_isolation_by_distance] Secondary interpretation = among-individual IBD using per-individual GPS (optional; does not replace main result).\n")

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

coerce_numeric_safely <- function(x, column_name = "unknown", context = "[11_isolation_by_distance]") {
  raw <- trimws(as.character(x))
  raw[raw == ""] <- NA_character_
  raw <- gsub(",", ".", raw, fixed = TRUE)
  num <- suppressWarnings(as.numeric(raw))
  bad <- !is.na(raw) & is.na(num)
  if (any(bad)) {
    warning(
      context,
      " Non-numeric values detected in column '", column_name,
      "'; affected rows will be excluded. Example values: ",
      paste(head(unique(raw[bad]), 5), collapse = ", ")
    )
  }
  num
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

score_id_column <- function(df) {
  nms <- names(df)
  std <- normalize_name(nms)
  scores <- rep(-Inf, length(nms))
  
  for (i in seq_along(nms)) {
    nm <- std[i]
    values <- trimws(as.character(df[[i]]))
    values <- values[nzchar(values) & !is.na(values)]
    if (length(values) == 0) next
    
    score <- 0
    if (nm %in% c("ind", "individual", "sample", "sampleid", "id", "ind_id")) score <- score + 100
    if (grepl("individual|sample|specimen|barcode|nom_labo|(^|_)id($|_)|(^|_)ind($|_)", nm)) score <- score + 50
    if (grepl("latitude|longitude|site|population|pop", nm)) score <- score - 60
    
    uniq_prop <- length(unique(values)) / length(values)
    if (uniq_prop >= 0.7) score <- score + 20
    scores[i] <- score
  }
  
  best_idx <- which.max(scores)
  if (length(best_idx) == 0 || !is.finite(scores[best_idx]) || scores[best_idx] <= 0) {
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
    if (!is_numeric_like(raw)) next
    x <- trimws(as.character(raw))
    x[x == ""] <- NA_character_
    x <- gsub(",", ".", x, fixed = TRUE)
    num <- suppressWarnings(as.numeric(x))
    non_missing <- !is.na(num)
    n_non_missing <- sum(non_missing)
    
    if (n_non_missing == 0) next
    if (grepl("sample|individual|specimen|nom_labo|echant|barcode|clone|mll|mlg|id", nm)) next
    
    val <- num[non_missing]
    if (type == "lat") {
      if (any(val < -90 | val > 90)) next
      score <- 0
      if (nm %in% c("latitude", "lat", "decimal_latitude", "lat_dd")) score <- score + 100
      if (grepl("latitude|(^|_)lat($|_)", nm)) score <- score + 50
      if (grepl("longitude|(^|_)lon($|_)|(^|_)long($|_)", nm)) score <- score - 120
    } else {
      if (any(val < -180 | val > 180)) next
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

assert_site_name_consistency <- function(reference_sites, candidate_sites, context) {
  ref_norm <- unique(normalize_site(reference_sites))
  cand_norm <- unique(normalize_site(candidate_sites))
  missing_in_candidate <- setdiff(ref_norm, cand_norm)
  extra_in_candidate <- setdiff(cand_norm, ref_norm)
  
  if (length(missing_in_candidate) > 0) {
    warning(
      context,
      " Site names present in the primary object but absent from the comparison table: ",
      paste(head(missing_in_candidate, 10), collapse = ", ")
    )
  }
  if (length(extra_in_candidate) > 0) {
    warning(
      context,
      " Site names present in the comparison table but absent from the primary object: ",
      paste(head(extra_in_candidate, 10), collapse = ", ")
    )
  }
}

collapse_site_coordinates <- function(meta_df, site_col, lat_col, lon_col, source_name) {
  coords_raw <- data.frame(
    Site = trimws(as.character(meta_df[[site_col]])),
    Latitude = coerce_numeric_safely(meta_df[[lat_col]], column_name = lat_col),
    Longitude = coerce_numeric_safely(meta_df[[lon_col]], column_name = lon_col),
    stringsAsFactors = FALSE
  )
  
  coords_raw$Site_norm <- normalize_site(coords_raw$Site)
  
  missing_site <- sum(!nzchar(coords_raw$Site) | is.na(coords_raw$Site))
  missing_lat <- sum(is.na(coords_raw$Latitude))
  missing_lon <- sum(is.na(coords_raw$Longitude))
  if (missing_site > 0 || missing_lat > 0 || missing_lon > 0) {
    warning(
      "[11_isolation_by_distance] Coordinate table '", source_name,
      "' contains missing values (missing Site rows=", missing_site,
      ", missing Latitude rows=", missing_lat,
      ", missing Longitude rows=", missing_lon,
      "). Incomplete rows will be excluded before computing site centroids."
    )
  }
  
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
      site = dplyr::first(Site),
      n_records = n(),
      mean_latitude = mean(Latitude, na.rm = TRUE),
      mean_longitude = mean(Longitude, na.rm = TRUE),
      min_latitude = min(Latitude, na.rm = TRUE),
      max_latitude = max(Latitude, na.rm = TRUE),
      min_longitude = min(Longitude, na.rm = TRUE),
      max_longitude = max(Longitude, na.rm = TRUE),
      coordinates_inconsistent_within_site =
        n_distinct(round(Latitude, 8)) > 1 | n_distinct(round(Longitude, 8)) > 1,
      .groups = "drop"
    ) %>%
    arrange(mean_latitude, site)
  
  centroid_sites <- coord_audit$site[coord_audit$n_records > 1]
  if (length(centroid_sites) > 0) {
    message(
      "[11_isolation_by_distance] Multiple coordinate records were detected for site(s): ",
      paste(centroid_sites, collapse = ", "),
      ". Mean latitude/longitude are being used as site centroids for site-level Mantel/IBD."
    )
  } else {
    message("[11_isolation_by_distance] One coordinate record per site detected; those coordinates are used directly as site centroids for site-level Mantel/IBD.")
  }
  
  inconsistent_sites <- coord_audit$site[coord_audit$coordinates_inconsistent_within_site]
  if (length(inconsistent_sites) > 0) {
    warning(
      "[11_isolation_by_distance] Inconsistent duplicate coordinates detected within site(s): ",
      paste(inconsistent_sites, collapse = ", "),
      ". Mean latitude/longitude per site are being used as centroids; review site_coordinate_audit.csv for coordinate ranges."
    )
  }
  
  coords <- coord_audit %>%
    transmute(
      Site = site,
      Site_norm = Site_norm,
      Latitude = mean_latitude,
      Longitude = mean_longitude,
      n_rows_collapsed = n_records,
      inconsistent_within_site = coordinates_inconsistent_within_site
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

build_individual_coordinates <- function(meta_df, df_ids_tbl) {
  if (!is.data.frame(meta_df) || ncol(meta_df) == 0) {
    warning("[11_isolation_by_distance] Optional individual-level IBD skipped: meta is empty.")
    return(NULL)
  }
  
  id_col <- score_id_column(meta_df)
  site_col <- score_site_column(meta_df)
  lat_col <- score_coordinate_column(meta_df, "lat")
  lon_col <- score_coordinate_column(meta_df, "lon")
  
  if (is.na(id_col) || is.na(lat_col) || is.na(lon_col)) {
    warning(
      "[11_isolation_by_distance] Optional individual-level IBD skipped: could not identify an ID / Latitude / Longitude column set in meta."
    )
    return(NULL)
  }
  
  id_norm_meta <- normalize_id(meta_df[[id_col]])
  if (anyDuplicated(id_norm_meta[nzchar(id_norm_meta)])) {
    dup_ids <- unique(trimws(as.character(meta_df[[id_col]])[duplicated(id_norm_meta) & nzchar(id_norm_meta)]))
    warning(
      "[11_isolation_by_distance] Duplicate individual IDs detected in the metadata used for optional individual-level IBD. Keeping the first record per ID. Examples: ",
      paste(head(dup_ids, 10), collapse = ", ")
    )
  }
  
  indiv_raw <- data.frame(
    Individual = trimws(as.character(meta_df[[id_col]])),
    Site = if (!is.na(site_col)) trimws(as.character(meta_df[[site_col]])) else NA_character_,
    Latitude = coerce_numeric_safely(meta_df[[lat_col]], column_name = lat_col),
    Longitude = coerce_numeric_safely(meta_df[[lon_col]], column_name = lon_col),
    stringsAsFactors = FALSE
  )
  indiv_raw$Individual_norm <- normalize_id(indiv_raw$Individual)
  indiv_raw <- indiv_raw[!duplicated(indiv_raw$Individual_norm), , drop = FALSE]
  
  missing_rows <- sum(!nzchar(indiv_raw$Individual) | is.na(indiv_raw$Latitude) | is.na(indiv_raw$Longitude))
  if (missing_rows > 0) {
    warning(
      "[11_isolation_by_distance] Optional individual-level IBD: ", missing_rows,
      " metadata row(s) have missing individual ID or coordinates and will be excluded."
    )
  }
  
  indiv_raw <- indiv_raw %>%
    filter(nzchar(Individual), !is.na(Latitude), !is.na(Longitude))
  
  if (nrow(indiv_raw) == 0) {
    warning("[11_isolation_by_distance] Optional individual-level IBD skipped: no individuals retained after filtering coordinates.")
    return(NULL)
  }
  
  if (any(indiv_raw$Latitude < -90 | indiv_raw$Latitude > 90)) {
    bad_ids <- unique(indiv_raw$Individual[indiv_raw$Latitude < -90 | indiv_raw$Latitude > 90])
    warning(
      "[11_isolation_by_distance] Optional individual-level IBD skipped: latitude outside [-90, 90] for individual(s): ",
      paste(head(bad_ids, 10), collapse = ", ")
    )
    return(NULL)
  }
  if (any(indiv_raw$Longitude < -180 | indiv_raw$Longitude > 180)) {
    bad_ids <- unique(indiv_raw$Individual[indiv_raw$Longitude < -180 | indiv_raw$Longitude > 180])
    warning(
      "[11_isolation_by_distance] Optional individual-level IBD skipped: longitude outside [-180, 180] for individual(s): ",
      paste(head(bad_ids, 10), collapse = ", ")
    )
    return(NULL)
  }
  
  site_map <- setNames(as.character(df_ids_tbl$Site), normalize_id(df_ids_tbl$ind_id))
  matched_site <- site_map[indiv_raw$Individual_norm]
  missing_site_match <- sum(is.na(matched_site) | !nzchar(matched_site))
  if (missing_site_match > 0) {
    warning(
      "[11_isolation_by_distance] Optional individual-level IBD: ", missing_site_match,
      " metadata individual(s) with coordinates could not be matched to df_ids and will be excluded."
    )
  }
  
  indiv_raw$Site <- matched_site
  indiv_raw <- indiv_raw %>%
    filter(!is.na(Site), nzchar(Site))
  
  if (nrow(indiv_raw) < 3) {
    warning("[11_isolation_by_distance] Optional individual-level IBD skipped: fewer than 3 individuals remain after alignment to df_ids.")
    return(NULL)
  }
  
  missing_in_meta <- setdiff(normalize_id(df_ids_tbl$ind_id), indiv_raw$Individual_norm)
  if (length(missing_in_meta) > 0) {
    warning(
      "[11_isolation_by_distance] Optional individual-level IBD: ", length(missing_in_meta),
      " gi individuals do not have usable per-individual coordinates in meta and will be omitted from the optional analysis."
    )
  }
  
  audit <- indiv_raw %>%
    transmute(
      Individual = Individual,
      Site = Site,
      Latitude = Latitude,
      Longitude = Longitude,
      source_id_column = id_col,
      source_latitude_column = lat_col,
      source_longitude_column = lon_col
    )
  
  list(data = indiv_raw, audit = audit, id_col = id_col, lat_col = lat_col, lon_col = lon_col)
}

run_individual_ibd_optional <- function(individual_coords) {
  if (is.null(individual_coords)) return(invisible(NULL))
  
  indiv_tbl <- individual_coords$data
  keep_idx <- match(normalize_id(adegenet::indNames(gi)), indiv_tbl$Individual_norm)
  keep <- !is.na(keep_idx)
  
  if (sum(keep) < 3) {
    warning("[11_isolation_by_distance] Optional individual-level IBD skipped: fewer than 3 gi individuals matched to per-individual coordinates.")
    return(invisible(NULL))
  }
  
  gi_ind <- gi[keep, , drop = FALSE]
  indiv_tbl <- indiv_tbl[keep_idx[keep], , drop = FALSE]
  
  if (!all(normalize_id(adegenet::indNames(gi_ind)) == indiv_tbl$Individual_norm)) {
    warning("[11_isolation_by_distance] Optional individual-level IBD skipped: gi and coordinate table could not be aligned identically.")
    return(invisible(NULL))
  }
  
  if (nrow(indiv_tbl) < 3) {
    warning("[11_isolation_by_distance] Optional individual-level IBD skipped after alignment: fewer than 3 individuals remain.")
    return(invisible(NULL))
  }
  
  coords_ind <- as.matrix(indiv_tbl[, c("Longitude", "Latitude")])
  storage.mode(coords_ind) <- "numeric"
  if (any(!is.finite(coords_ind))) {
    warning("[11_isolation_by_distance] Optional individual-level IBD skipped: non-finite individual coordinates remain after alignment.")
    return(invisible(NULL))
  }
  
  geographic_ind_km <- geosphere::distm(coords_ind, fun = geosphere::distHaversine) / 1000
  rownames(geographic_ind_km) <- indiv_tbl$Individual
  colnames(geographic_ind_km) <- indiv_tbl$Individual
  
  genetic_ind <- as.matrix(poppr::prevosti.dist(gi_ind))
  storage.mode(genetic_ind) <- "numeric"
  rownames(genetic_ind) <- indiv_tbl$Individual
  colnames(genetic_ind) <- indiv_tbl$Individual
  
  if (!identical(rownames(genetic_ind), rownames(geographic_ind_km))) {
    warning("[11_isolation_by_distance] Optional individual-level IBD skipped: label mismatch between individual genetic and geographic matrices.")
    return(invisible(NULL))
  }
  
  diag(genetic_ind) <- 0
  diag(geographic_ind_km) <- 0
  
  if (any(!is.finite(genetic_ind[row(genetic_ind) != col(genetic_ind)]))) {
    warning("[11_isolation_by_distance] Optional individual-level IBD skipped: individual genetic distance matrix contains non-finite off-diagonal values.")
    return(invisible(NULL))
  }
  if (any(!is.finite(geographic_ind_km[row(geographic_ind_km) != col(geographic_ind_km)]))) {
    warning("[11_isolation_by_distance] Optional individual-level IBD skipped: individual geographic distance matrix contains non-finite off-diagonal values.")
    return(invisible(NULL))
  }
  
  mantel_ind <- vegan::mantel(
    as.dist(genetic_ind),
    as.dist(geographic_ind_km),
    method = "pearson",
    permutations = 9999
  )
  
  ind_geo_file <- file.path(MATRICES_DIR, "individual_geographic_distance_matrix.csv")
  write.csv(geographic_ind_km, ind_geo_file, row.names = TRUE)
  cat("[11_isolation_by_distance] Saved optional individual geographic matrix: ", ind_geo_file, "\n", sep = "")
  
  ind_gen_file <- file.path(MATRICES_DIR, "individual_genetic_distance_matrix.csv")
  write.csv(genetic_ind, ind_gen_file, row.names = TRUE)
  cat("[11_isolation_by_distance] Saved optional individual genetic matrix: ", ind_gen_file, "\n", sep = "")
  
  upper_idx <- upper.tri(genetic_ind)
  ind_points <- data.frame(
    Individual1 = rownames(genetic_ind)[row(genetic_ind)[upper_idx]],
    Individual2 = colnames(genetic_ind)[col(genetic_ind)[upper_idx]],
    Site1 = indiv_tbl$Site[row(genetic_ind)[upper_idx]],
    Site2 = indiv_tbl$Site[col(genetic_ind)[upper_idx]],
    Geographic_km = as.numeric(geographic_ind_km[upper_idx]),
    Genetic_distance_prevosti = as.numeric(genetic_ind[upper_idx]),
    stringsAsFactors = FALSE
  )
  
  ind_points_file <- file.path(TABLES_DIR, "individual_ibd_points.csv")
  write.csv(ind_points, ind_points_file, row.names = FALSE)
  cat("[11_isolation_by_distance] Saved optional individual IBD points: ", ind_points_file, "\n", sep = "")
  
  ind_results <- data.frame(
    test = "Mantel_Individual_Prevosti_vs_Geographic",
    statistic_r = as.numeric(mantel_ind$statistic),
    p_value = as.numeric(mantel_ind$signif),
    permutations = as.integer(mantel_ind$permutations),
    n_individuals = nrow(indiv_tbl),
    genotype_object = "gi",
    genetic_distance_input = "poppr::prevosti.dist(gi)",
    coordinate_source = "meta object (per-individual GPS)",
    id_column = individual_coords$id_col,
    latitude_column = individual_coords$lat_col,
    longitude_column = individual_coords$lon_col,
    interpretation_note = paste(
      "Optional secondary analysis.",
      "Unlike the primary site-level Mantel test, this individual-level analysis mixes within-site and among-site structure and should be interpreted cautiously."
    ),
    stringsAsFactors = FALSE
  )
  
  ind_results_file <- file.path(TABLES_DIR, "individual_mantel_test_results.csv")
  write.csv(ind_results, ind_results_file, row.names = FALSE)
  cat("[11_isolation_by_distance] Saved optional individual Mantel results: ", ind_results_file, "\n", sep = "")
  
  invisible(NULL)
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

assert_site_name_consistency(rownames(pairwise_jostD), site_meta$Site, "[11_isolation_by_distance]")
assert_site_name_consistency(unique(as.character(df_ids_mll$Site)), site_meta$Site, "[11_isolation_by_distance]")

audit_file <- file.path(TABLES_DIR, "site_coordinate_audit.csv")
write.csv(coord_audit, audit_file, row.names = FALSE)
cat("[11_isolation_by_distance] Saved coordinate audit: ", audit_file, "\n", sep = "")

site_coord_file <- file.path(TABLES_DIR, "site_coordinates_clean.csv")
write.csv(site_meta, site_coord_file, row.names = FALSE)
cat("[11_isolation_by_distance] Saved clean site coordinates: ", site_coord_file, "\n", sep = "")

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
  warning(
    "[11_isolation_by_distance] Dropping genetic-only normalized sites: ",
    paste(setdiff(gen_sites_norm, shared_norm), collapse = ", ")
  )
}
if (length(setdiff(meta_sites_norm, shared_norm)) > 0) {
  warning(
    "[11_isolation_by_distance] Dropping coordinate-only normalized sites: ",
    paste(setdiff(meta_sites_norm, shared_norm), collapse = ", ")
  )
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
# 4) Mantel test (primary site-level analysis)
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
  interpretation_note = paste(
    "Primary site-level IBD analysis.",
    "Geographic distances are calculated among mean site centroids, not among individual tree coordinates."
  ),
  stringsAsFactors = FALSE
)

mantel_file <- file.path(TABLES_DIR, "mantel_test_results.csv")
write.csv(mantel_results, mantel_file, row.names = FALSE)
cat("[11_isolation_by_distance] Saved Mantel results: ", mantel_file, "\n", sep = "")

# ----------------------------
# 5) Export IBD points and plot (primary site-level)
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
    subtitle = "Primary analysis: clone-corrected pairwise Jost's D by site versus geographic distance among mean site centroids",
    x = "Geographic distance (km)",
    y = "Pairwise Jost's D"
  )

ibd_plot_file <- file.path(FIGURES_DIR, "isolation_by_distance.jpeg")
ggsave(ibd_plot_file, plot = ibd_plot, width = 8, height = 6, dpi = 320)
cat("[11_isolation_by_distance] Saved figure: ", ibd_plot_file, "\n", sep = "")

# ----------------------------
# 6) OPTIONAL secondary individual-level IBD using per-individual GPS
# ----------------------------
individual_coords <- build_individual_coordinates(meta, df_ids)
if (!is.null(individual_coords)) {
  individual_audit_file <- file.path(TABLES_DIR, "individual_coordinate_audit.csv")
  write.csv(individual_coords$audit, individual_audit_file, row.names = FALSE)
  cat("[11_isolation_by_distance] Saved optional individual coordinate audit: ", individual_audit_file, "\n", sep = "")
}
run_individual_ibd_optional(individual_coords)