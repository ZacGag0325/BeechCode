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

normalize_name <- function(x) {
  x <- iconv(x, from = "", to = "ASCII//TRANSLIT")
  x[is.na(x)] <- ""
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x
}

build_name_map <- function(df) {
  raw <- names(df)
  norm <- normalize_name(raw)
  names(raw) <- norm
  raw
}

first_col_by_patterns <- function(df, patterns) {
  nm <- build_name_map(df)
  keys <- names(nm)
  for (pat in patterns) {
    hit <- which(grepl(pat, keys, perl = TRUE))
    if (length(hit) > 0) return(nm[keys[hit[1]]])
  }
  NA_character_
}

extract_site_coordinates <- function(df, source_label = "unknown") {
  if (!is.data.frame(df) || nrow(df) == 0) {
    return(list(ok = FALSE, reason = paste0(source_label, " is empty")))
  }
  
  site_col <- first_col_by_patterns(
    df,
    c(
      "^site$", "^sites$", "^site_id$", "^population$", "^pop$",
      "^numero_population$", "^numero_de_population$", "^num_population$", "^numeropopulation$", "^numerodepopulation$", "^numeropop$"
    )
  )
  
  if (is.na(site_col)) {
    return(list(ok = FALSE, reason = paste0("Could not detect site column in ", source_label)))
  }
  
  lat_col <- first_col_by_patterns(df, c("^latitude$", "^lat$", "^y$", "^northing$", "^coord_y$", "^utm_y$", "^y_coord$"))
  lon_col <- first_col_by_patterns(df, c("^longitude$", "^lon$", "^long$", "^x$", "^easting$", "^coord_x$", "^utm_x$", "^x_coord$"))
  
  if (is.na(lat_col) || is.na(lon_col)) {
    return(list(ok = FALSE, reason = paste0("Could not detect coordinate columns in ", source_label), site_col = site_col))
  }
  
  coord_df <- df %>%
    mutate(
      Site = normalize_site(.data[[site_col]]),
      coord_a = suppressWarnings(as.numeric(.data[[lat_col]])),
      coord_b = suppressWarnings(as.numeric(.data[[lon_col]]))
    ) %>%
    filter(!is.na(Site), !is.na(coord_a), !is.na(coord_b))
  
  if (!nrow(coord_df)) {
    return(list(ok = FALSE, reason = paste0("No valid coordinates after cleaning for ", source_label), site_col = site_col, lat_col = lat_col, lon_col = lon_col))
  }
  
  mode <- {
    lat_norm <- normalize_name(lat_col)
    lon_norm <- normalize_name(lon_col)
    if (grepl("lat|latitude", lat_norm) || grepl("lon|long|longitude", lon_norm)) "latlon" else "projected"
  }
  
  site_df <- coord_df %>%
    group_by(Site) %>%
    summarise(coord_a = mean(coord_a, na.rm = TRUE), coord_b = mean(coord_b, na.rm = TRUE), .groups = "drop")
  
  list(
    ok = TRUE,
    data = site_df,
    mode = mode,
    source = source_label,
    site_col = site_col,
    lat_col = lat_col,
    lon_col = lon_col
  )
}

read_site_metadata <- function(meta_path, df_ids = NULL) {
  meta_obj <- NULL
  
  if (file.exists(meta_path)) {
    meta_try <- tryCatch(readr::read_csv(meta_path, show_col_types = FALSE), error = function(e) e)
    if (inherits(meta_try, "error")) {
      warning("Failed reading site metadata file: ", meta_path, " -> ", conditionMessage(meta_try))
    } else {
      meta_obj <- extract_site_coordinates(meta_try, source_label = "inputs/site_metadata.csv")
      if (isTRUE(meta_obj$ok)) {
        message("Using inputs/site_metadata.csv")
        message("Detected columns -> site: ", meta_obj$site_col, ", coord1: ", meta_obj$lat_col, ", coord2: ", meta_obj$lon_col)
        return(meta_obj)
      }
      warning(meta_obj$reason)
    }
  } else {
    warning("Missing site metadata file: ", meta_path)
  }
  
  if (!is.null(df_ids)) {
    fallback <- extract_site_coordinates(as.data.frame(df_ids), source_label = "df_ids")
    if (isTRUE(fallback$ok)) {
      message("Falling back to df_ids site means")
      message("Detected columns -> site: ", fallback$site_col, ", coord1: ", fallback$lat_col, ", coord2: ", fallback$lon_col)
      return(fallback)
    }
    warning("df_ids fallback failed: ", fallback$reason)
  }
  
  list(ok = FALSE, reason = "No usable site coordinate source (metadata + df_ids fallback failed).")
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

compute_geo_matrix <- function(meta_obj, site_order = NULL) {
  if (is.null(meta_obj) || !isTRUE(meta_obj$ok)) {
    warning("No usable coordinate metadata; skipping geographic matrix.")
    return(NULL)
  }
  
  dat <- meta_obj$data
  dat <- dat[!is.na(dat$coord_a) & !is.na(dat$coord_b), , drop = FALSE]
  
  if (!is.null(site_order)) {
    dat <- dat[dat$Site %in% site_order, , drop = FALSE]
  }
  
  if (nrow(dat) < 2) {
    warning("Too few sites with coordinates to compute geographic distances.")
    return(NULL)
  }
  
  if (!is.null(site_order)) {
    dat <- dat[match(site_order, dat$Site), , drop = FALSE]
    dat <- dat[!is.na(dat$Site), , drop = FALSE]
  } else {
    dat <- dat[order(dat$Site), , drop = FALSE]
  }
  
  if (nrow(dat) < 2) {
    warning("Too few aligned sites to compute geographic distances.")
    return(NULL)
  }
  
  if (meta_obj$mode == "latlon") {
    mat <- geosphere::distm(as.matrix(dat[, c("coord_b", "coord_a")]), fun = geosphere::distHaversine) / 1000
    message("Computed geographic matrix with Haversine distances (km).")
  } else if (meta_obj$mode == "projected") {
    mat <- as.matrix(stats::dist(dat[, c("coord_b", "coord_a")])) / 1000
    message("Computed geographic matrix with Euclidean distances (km) from projected coordinates.")
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
meta_obj <- read_site_metadata(meta_path, df_ids = if (exists("df_ids")) df_ids else NULL)
if (!isTRUE(meta_obj$ok)) {
  warning(meta_obj$reason)
  geo_mat <- NULL
} else {
  target_sites <- if (!is.null(nei_mat)) rownames(nei_mat) else sort(unique(meta_obj$data$Site))
  geo_mat <- tryCatch(
    compute_geo_matrix(meta_obj, site_order = target_sites),
    error = function(e) {
      warning("Geographic matrix computation failed: ", conditionMessage(e))
      NULL
    }
  )
}

save_matrix_outputs(geo_mat, geo_csv, geo_rds, "geographic distance matrix")
if (!is.null(geo_mat)) {
  write.csv(geo_mat, geo_csv_legacy, row.names = TRUE)
  message("Created legacy geographic matrix CSV: ", geo_csv_legacy)
}

message("Distance-matrix step complete. Output folder: ", DIST_DIR)
