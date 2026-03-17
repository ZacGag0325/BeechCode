# scripts/11_isolation_by_distance.R
############################################################
# Isolation by Distance (Mantel test)
# Genetic distance input: pairwise Jost's D (clone-corrected, site-level)
# Outputs:
# - outputs/tables/mantel_test_results.csv
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

normalize_site <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("\\uFEFF", "", x, fixed = TRUE)
  x <- gsub("[[:cntrl:]]", "", x)
  x <- gsub("\\s+", " ", x)
  toupper(x)
}

safe_read_square_matrix <- function(path) {
  if (!file.exists(path)) {
    stop("Missing file: ", path)
  }
  
  cat("[11_isolation_by_distance] Loading genetic matrix from: ", path, "\n", sep = "")
  
  # First attempt: canonical CSV written with row.names = TRUE
  df <- tryCatch(
    read.csv(path, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE, strip.white = TRUE),
    error = function(e) NULL
  )
  
  # Fallback: permissive read and manual row-id extraction
  if (is.null(df)) {
    raw_df <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, strip.white = TRUE)
    if (ncol(raw_df) < 2) {
      stop("Imported matrix has fewer than 2 columns: ", path)
    }
    row_id <- trimws(as.character(raw_df[[1]]))
    df <- raw_df[, -1, drop = FALSE]
    rownames(df) <- row_id
  }
  
  # Trim row/column names and drop empty row names
  row_id <- trimws(rownames(df))
  colnames(df) <- trimws(colnames(df))
  
  keep <- nzchar(row_id)
  if (!all(keep)) {
    cat("[11_isolation_by_distance] Dropping ", sum(!keep), " rows with empty site names from imported matrix.\n", sep = "")
    df <- df[keep, , drop = FALSE]
    row_id <- row_id[keep]
  }
  rownames(df) <- row_id
  
  # Some exports accidentally include a duplicate ID column after row names.
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
  
  # Convert every column safely to numeric and report exactly which columns fail.
  bad_cols <- character(0)
  for (j in seq_len(ncol(df))) {
    x <- trimws(as.character(df[[j]]))
    x[x == ""] <- NA_character_
    
    # Handle decimal comma if present (e.g., 0,123)
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
  
  # Reorder columns to exactly match row order
  mat <- mat[, rownames(mat), drop = FALSE]
  
  if (!identical(rownames(mat), colnames(mat))) {
    stop("Imported matrix row/column names could not be aligned identically: ", path)
  }
  
  if (any(!is.finite(mat), na.rm = TRUE)) {
    stop("Imported genetic matrix contains non-finite values.")
  }
  
  cat("[11_isolation_by_distance] Genetic matrix dimensions: ", nrow(mat), " x ", ncol(mat), "\n", sep = "")
  cat("[11_isolation_by_distance] Genetic matrix sites (first 10): ", paste(head(rownames(mat), 10), collapse = ", "), "\n", sep = "")
  
  mat
}

resolve_col <- function(df, candidates, regex = NULL) {
  nms <- names(df)
  nms_low <- tolower(nms)
  
  if (!is.null(candidates)) {
    idx <- match(TRUE, nms_low %in% tolower(candidates), nomatch = 0)
    if (idx > 0) return(nms[idx])
  }
  
  if (!is.null(regex)) {
    idx <- grep(regex, nms_low)
    if (length(idx) > 0) return(nms[idx[1]])
  }
  
  NA_character_
}

build_site_coordinates <- function(meta_df, fallback_path) {
  site_col <- resolve_col(meta_df, c("site", "population", "pop"))
  lat_col <- resolve_col(meta_df, NULL, regex = "lat")
  lon_col <- resolve_col(meta_df, NULL, regex = "lon|long")
  
  source_name <- "meta object (outputs/v1/objects/meta.rds)"
  
  if (is.na(site_col) || is.na(lat_col) || is.na(lon_col)) {
    if (!file.exists(fallback_path)) {
      stop(
        "Could not find Site/Latitude/Longitude columns in meta, and fallback file is missing: ",
        fallback_path
      )
    }
    
    cat("[11_isolation_by_distance] Metadata coordinates not found in meta; using fallback: ", fallback_path, "\n", sep = "")
    meta_df <- read.csv(fallback_path, stringsAsFactors = FALSE, check.names = FALSE)
    site_col <- resolve_col(meta_df, c("site", "population", "pop"))
    lat_col <- resolve_col(meta_df, NULL, regex = "lat")
    lon_col <- resolve_col(meta_df, NULL, regex = "lon|long")
    source_name <- "inputs/site_metadata.csv"
  }
  
  if (is.na(site_col) || is.na(lat_col) || is.na(lon_col)) {
    stop("Metadata must contain site/pop, latitude, and longitude columns.")
  }
  
  cat("[11_isolation_by_distance] Detected metadata Site column: ", site_col, "\n", sep = "")
  cat("[11_isolation_by_distance] Detected metadata Latitude column: ", lat_col, "\n", sep = "")
  cat("[11_isolation_by_distance] Detected metadata Longitude column: ", lon_col, "\n", sep = "")
  
  coords <- data.frame(
    Site = trimws(as.character(meta_df[[site_col]])),
    Latitude = suppressWarnings(as.numeric(meta_df[[lat_col]])),
    Longitude = suppressWarnings(as.numeric(meta_df[[lon_col]])),
    stringsAsFactors = FALSE
  )
  
  coords <- coords %>%
    filter(nzchar(Site), !is.na(Latitude), !is.na(Longitude))
  
  if (nrow(coords) == 0) {
    stop("No valid site coordinates available after filtering missing values.")
  }
  
  if (anyDuplicated(coords$Site)) {
    coords <- coords %>%
      mutate(Site_norm = normalize_site(Site)) %>%
      group_by(Site_norm) %>%
      slice(1) %>%
      ungroup() %>%
      select(-Site_norm)
    cat("[11_isolation_by_distance] Duplicate site rows found in metadata; keeping first row per normalized site name.\n")
  }
  
  cat("[11_isolation_by_distance] Coordinate source: ", source_name, "\n", sep = "")
  cat("[11_isolation_by_distance] Coordinate table sites (first 10): ", paste(head(coords$Site, 10), collapse = ", "), "\n", sep = "")
  
  coords
}

# ----------------------------
# 1) Load genetic distance matrix (Jost's D)
# ----------------------------
jost_file <- file.path(MATRICES_DIR, "pairwise_jostD.csv")
pairwise_jostD <- safe_read_square_matrix(jost_file)

# ----------------------------
# 2) Load site coordinates from metadata
# ----------------------------
site_meta <- build_site_coordinates(
  meta_df = meta,
  fallback_path = file.path(PROJECT_ROOT, "inputs", "site_metadata.csv")
)

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
    "Need at least 3 shared sites between pairwise_jostD.csv and metadata coordinates. Found: ",
    length(shared_norm),
    "\nGenetic normalized sites: ", paste(gen_sites_norm, collapse = ", "),
    "\nMetadata normalized sites: ", paste(meta_sites_norm, collapse = ", ")
  )
}

# Keep original display names from metadata when available
meta_key <- setNames(site_meta$Site, meta_sites_norm)
shared_display <- unname(meta_key[shared_norm])
shared_display[is.na(shared_display)] <- shared_norm[is.na(shared_display)]

# subset and reorder genetic matrix by normalized keys
idx_gen <- match(shared_norm, gen_sites_norm)
pairwise_jostD <- pairwise_jostD[idx_gen, idx_gen, drop = FALSE]
rownames(pairwise_jostD) <- shared_display
colnames(pairwise_jostD) <- shared_display

# subset and reorder site metadata by normalized keys
idx_meta <- match(shared_norm, meta_sites_norm)
site_meta <- site_meta[idx_meta, , drop = FALSE]
site_meta$Site <- shared_display

coords <- as.matrix(site_meta[, c("Longitude", "Latitude")])
geographic_km <- geosphere::distm(coords, fun = geosphere::distHaversine) / 1000
rownames(geographic_km) <- shared_display
colnames(geographic_km) <- shared_display

geo_file <- file.path(MATRICES_DIR, "geographic_distance_matrix.csv")
write.csv(geographic_km, geo_file, row.names = TRUE)
cat("[11_isolation_by_distance] Saved geographic matrix: ", geo_file, "\n", sep = "")

# ----------------------------
# 4) Mantel test
# ----------------------------
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
    x = "Geographic distance (km)",
    y = "Pairwise Jost's D"
  )

ibd_plot_file <- file.path(FIGURES_DIR, "isolation_by_distance.jpeg")
ggsave(ibd_plot_file, plot = ibd_plot, width = 8, height = 6, dpi = 320)
cat("[11_isolation_by_distance] Saved figure: ", ibd_plot_file, "\n", sep = "")