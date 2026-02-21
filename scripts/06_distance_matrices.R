############################################################
# scripts/06_distance_matrices.R
# Build and save Nei genetic + geographic distance matrices
############################################################

options(repos = c(CRAN = "https://cloud.r-project.org"))
pkgs <- c("adegenet", "poppr", "geosphere", "readr", "dplyr")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)

suppressPackageStartupMessages({
  library(adegenet)
  library(poppr)
  library(geosphere)
  library(readr)
  library(dplyr)
})

find_project_root <- function() {
  candidates <- c(getwd(), normalizePath(file.path(getwd(), ".."), mustWork = FALSE))
  cmd_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
  if (length(cmd_file) > 0 && nzchar(cmd_file[1])) {
    candidates <- c(candidates, dirname(normalizePath(cmd_file[1], mustWork = FALSE)))
  }
  for (start in unique(candidates)) {
    cur <- normalizePath(start, mustWork = FALSE)
    repeat {
      if (file.exists(file.path(cur, "scripts", "_load_objects.R"))) return(cur)
      parent <- dirname(cur)
      if (identical(parent, cur)) break
      cur <- parent
    }
  }
  stop("Cannot find project root containing scripts/_load_objects.R")
}

normalize_site <- function(x) {
  x <- trimws(as.character(x))
  x[x == ""] <- NA_character_
  x
}

match_first_col <- function(df, patterns) {
  nms <- names(df)
  nms_low <- tolower(nms)
  for (pat in patterns) {
    hit <- which(grepl(pat, nms_low, perl = TRUE))
    if (length(hit) > 0) return(nms[hit[1]])
  }
  NA_character_
}

read_site_metadata <- function(meta_path) {
  meta <- readr::read_csv(meta_path, show_col_types = FALSE)
  if (!nrow(meta)) stop("site_metadata.csv is empty: ", meta_path)
  
  site_col <- match_first_col(meta, c("^site$", "site", "population", "pop"))
  if (is.na(site_col)) stop("Could not detect site column in ", meta_path)
  
  lat_col <- match_first_col(meta, c("^lat$", "latitude", "^y$", "^y_coord", "northing"))
  lon_col <- match_first_col(meta, c("^lon$", "^long$", "longitude", "^x$", "^x_coord", "easting"))
  
  out <- meta %>%
    mutate(
      Site = normalize_site(.data[[site_col]]),
      coord_a = suppressWarnings(as.numeric(.data[[lat_col]])),
      coord_b = suppressWarnings(as.numeric(.data[[lon_col]]))
    ) %>%
    filter(!is.na(Site)) %>%
    distinct(Site, .keep_all = TRUE)
  
  if (is.na(lat_col) || is.na(lon_col)) {
    warning("Could not detect coordinate columns in site metadata.")
    return(list(data = out, mode = "none", lat_col = lat_col, lon_col = lon_col))
  }
  
  lat_name <- tolower(lat_col)
  lon_name <- tolower(lon_col)
  likely_latlon <- any(grepl("lat|lon|long|latitude|longitude", c(lat_name, lon_name)))
  
  mode <- if (likely_latlon) "latlon" else "projected"
  list(data = out, mode = mode, lat_col = lat_col, lon_col = lon_col)
}

compute_nei_matrix <- function(gi) {
  if (!inherits(gi, "genind")) stop("gi must be a genind object")
  if (adegenet::nInd(gi) < 2 || adegenet::nPop(gi) < 2) {
    warning("Too few individuals/populations to compute Nei distance.")
    return(NULL)
  }
  
  gp <- adegenet::genind2genpop(gi)
  nei <- tryCatch(
    adegenet::dist.genpop(gp, method = 1),
    error = function(e) tryCatch(poppr::nei.dist(gp), error = function(e2) e2)
  )
  if (inherits(nei, "error")) {
    warning("Nei distance computation failed: ", conditionMessage(nei))
    return(NULL)
  }
  
  mat <- as.matrix(nei)
  mat <- (mat + t(mat)) / 2
  diag(mat) <- 0
  mat
}

compute_geo_matrix <- function(meta_obj, gi_sites = NULL) {
  dat <- meta_obj$data
  dat <- dat[!is.na(dat$coord_a) & !is.na(dat$coord_b), , drop = FALSE]
  if (!is.null(gi_sites)) dat <- dat[dat$Site %in% gi_sites, , drop = FALSE]
  if (nrow(dat) < 2) {
    warning("Too few sites with coordinates to compute geographic distances.")
    return(NULL)
  }
  
  dat <- dat[order(dat$Site), , drop = FALSE]
  
  if (meta_obj$mode == "latlon") {
    mat <- geosphere::distm(as.matrix(dat[, c("coord_b", "coord_a")]), fun = geosphere::distHaversine) / 1000
    message("Computed geographic matrix with Haversine distances (km) using columns: ", meta_obj$lat_col, ", ", meta_obj$lon_col)
  } else if (meta_obj$mode == "projected") {
    mat <- as.matrix(stats::dist(dat[, c("coord_b", "coord_a")])) / 1000
    message("Computed geographic matrix with Euclidean distances (km) using projected coordinates: ", meta_obj$lat_col, ", ", meta_obj$lon_col)
  } else {
    warning("No usable coordinate mode detected.")
    return(NULL)
  }
  
  rownames(mat) <- dat$Site
  colnames(mat) <- dat$Site
  mat <- (mat + t(mat)) / 2
  diag(mat) <- 0
  mat
}

save_matrix_outputs <- function(mat, csv_path, rds_path, label) {
  if (is.null(mat)) {
    message("Skipped saving ", label, " (matrix unavailable).")
    return(invisible(FALSE))
  }
  write.csv(mat, csv_path, row.names = TRUE)
  saveRDS(mat, rds_path)
  message("Created ", label, " CSV: ", csv_path)
  message("Created ", label, " RDS: ", rds_path)
  invisible(TRUE)
}

PROJECT_ROOT <- find_project_root()
setwd(PROJECT_ROOT)
source(file.path(PROJECT_ROOT, "scripts", "_load_objects.R"))

RUN_TAG <- if (exists("RUN_TAG", inherits = TRUE)) get("RUN_TAG", inherits = TRUE) else "v1"
RUN_OUT <- if (exists("RUN_OUT", inherits = TRUE)) get("RUN_OUT", inherits = TRUE) else file.path(PROJECT_ROOT, "outputs", RUN_TAG)
DIST_DIR <- file.path(PROJECT_ROOT, "outputs", RUN_TAG, "distances")
INPUTS_DIR <- file.path(PROJECT_ROOT, "inputs")
dir.create(file.path(PROJECT_ROOT, "outputs"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(PROJECT_ROOT, "outputs", RUN_TAG), showWarnings = FALSE, recursive = TRUE)
dir.create(DIST_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(INPUTS_DIR, showWarnings = FALSE, recursive = TRUE)

nei_csv <- file.path(DIST_DIR, "genetic_nei_dist.csv")
nei_rds <- file.path(DIST_DIR, "genetic_nei_dist.rds")
geo_csv <- file.path(DIST_DIR, "geographic_dist.csv")
geo_rds <- file.path(DIST_DIR, "geographic_dist.rds")

# Legacy compatibility paths expected by older downstream scripts
nei_csv_legacy <- file.path(RUN_OUT, "matrix_genetic_distance_nei.csv")
geo_csv_legacy <- file.path(RUN_OUT, "matrix_geographic_distance_km.csv")

nei_mat <- compute_nei_matrix(gi)
save_matrix_outputs(nei_mat, nei_csv, nei_rds, "Nei genetic distance matrix")
if (!is.null(nei_mat)) {
  write.csv(nei_mat, nei_csv_legacy, row.names = TRUE)
  message("Created legacy Nei matrix CSV: ", nei_csv_legacy)
}

meta_path <- file.path(INPUTS_DIR, "site_metadata.csv")
if (!file.exists(meta_path)) {
  warning("Missing site metadata file: ", meta_path, " (cannot compute geographic matrix).")
  geo_mat <- NULL
} else {
  meta_obj <- read_site_metadata(meta_path)
  geo_mat <- compute_geo_matrix(meta_obj, gi_sites = if (!is.null(nei_mat)) rownames(nei_mat) else NULL)
}
save_matrix_outputs(geo_mat, geo_csv, geo_rds, "geographic distance matrix")
if (!is.null(geo_mat)) {
  write.csv(geo_mat, geo_csv_legacy, row.names = TRUE)
  message("Created legacy geographic matrix CSV: ", geo_csv_legacy)
}

message("Distance-matrix step complete. Output folder: ", DIST_DIR)
