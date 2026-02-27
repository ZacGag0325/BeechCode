############################################################
# scripts/91_build_site_coordinates_from_raw.R
# Build site-level coordinates database from raw field data
############################################################

options(repos = c(CRAN = "https://cloud.r-project.org"))

pkgs <- c("readxl", "readr", "dplyr", "stringi")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

suppressPackageStartupMessages({
  library(readxl)
  library(readr)
  library(dplyr)
  library(stringi)
})

find_project_root <- function() {
  cur <- normalizePath(getwd(), mustWork = FALSE)
  repeat {
    if (file.exists(file.path(cur, "00_master_pipeline.R"))) return(cur)
    parent <- dirname(cur)
    if (identical(parent, cur)) break
    cur <- parent
  }
  stop("MISSING FROM YOUR SIDE:\n- Cannot find project root (00_master_pipeline.R).")
}

to_utf8 <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  if (!is.character(x)) return(x)
  y <- suppressWarnings(iconv(x, from = "", to = "UTF-8", sub = ""))
  y[is.na(y)] <- ""
  y
}

clean_names <- function(x) {
  x <- to_utf8(as.character(x))
  x <- trimws(x)
  x <- stringi::stri_trans_general(x, "Latin-ASCII")
  x <- gsub("\\s+", "_", x)
  x <- gsub("[^a-zA-Z0-9_]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  tolower(x)
}

read_any_tabular <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("xlsx", "xls")) {
    x <- readxl::read_excel(path)
  } else if (ext == "csv") {
    x <- suppressWarnings(readr::read_csv(path, show_col_types = FALSE, guess_max = 50000))
  } else {
    x <- utils::read.table(path, header = TRUE, fill = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  }
  x <- as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE)
  names(x) <- clean_names(names(x))
  for (j in seq_along(x)) {
    if (is.character(x[[j]]) || is.factor(x[[j]])) x[[j]] <- to_utf8(as.character(x[[j]]))
  }
  x
}

choose_col <- function(nms, candidates) {
  nms <- clean_names(nms)
  cands <- clean_names(candidates)
  i <- which(nms %in% cands)
  if (length(i) == 0) return(NA_character_)
  nms[i[1]]
}

PROJECT_ROOT <- find_project_root()
RAW_DIR <- file.path(PROJECT_ROOT, "data", "raw")
IN_DIR <- file.path(PROJECT_ROOT, "inputs")
dir.create(IN_DIR, recursive = TRUE, showWarnings = FALSE)

field_candidates <- c(
  file.path(RAW_DIR, "donnees_modifiees_west_summer2024 copie.xlsx"),
  file.path(RAW_DIR, "donnees_modifiees_west_summer2024.xlsx"),
  file.path(RAW_DIR, "field_data.xlsx"),
  file.path(RAW_DIR, "field_data.csv")
)

source_file <- field_candidates[file.exists(field_candidates)][1]
if (is.na(source_file)) {
  stop(
    "MISSING FROM YOUR SIDE:\n",
    "- No field raw file found in data/raw/.\n",
    "- Provide one of: donnees_modifiees_west_summer2024 copie.xlsx, donnees_modifiees_west_summer2024.xlsx, field_data.xlsx, field_data.csv."
  )
}

cat("[build_site_coordinates] source:", source_file, "\n")
field <- read_any_tabular(source_file)

site_candidates <- c("site", "population", "pop", "numero_population", "num_population", "site_id")
lon_candidates <- c("lon", "long", "longitude", "x")
lat_candidates <- c("lat", "latitude", "y")

site_col <- choose_col(names(field), site_candidates)
lon_col <- choose_col(names(field), lon_candidates)
lat_col <- choose_col(names(field), lat_candidates)

# fallback to site_metadata when lat/lon missing in field file
site_meta_path <- file.path(IN_DIR, "site_metadata.csv")
site_meta <- NULL
if (file.exists(site_meta_path)) {
  site_meta <- read_any_tabular(site_meta_path)
}

if (is.na(site_col) && !is.null(site_meta)) site_col <- choose_col(names(site_meta), site_candidates)
if (is.na(lon_col) && !is.null(site_meta)) lon_col <- choose_col(names(site_meta), lon_candidates)
if (is.na(lat_col) && !is.null(site_meta)) lat_col <- choose_col(names(site_meta), lat_candidates)

if (is.na(site_col) || is.na(lon_col) || is.na(lat_col)) {
  stop(
    "MISSING FROM YOUR SIDE:\n",
    "- Cannot identify site/lon/lat columns from raw field file (or site_metadata fallback).\n",
    "- Expected minimum columns: site, lon, lat."
  )
}

coords <- field %>%
  transmute(
    site = trimws(to_utf8(as.character(.data[[site_col]]))),
    lon = suppressWarnings(as.numeric(.data[[lon_col]])),
    lat = suppressWarnings(as.numeric(.data[[lat_col]]))
  ) %>%
  filter(site != "", !is.na(lon), !is.na(lat)) %>%
  group_by(site) %>%
  summarise(lon = mean(lon, na.rm = TRUE), lat = mean(lat, na.rm = TRUE), .groups = "drop")

if (nrow(coords) == 0) {
  stop("MISSING FROM YOUR SIDE:\n- No valid coordinates extracted (all lon/lat missing or invalid).")
}

coords_path <- file.path(IN_DIR, "site_coordinates.csv")
readr::write_csv(coords, coords_path)

cat("[build_site_coordinates] written:\n")
cat(" - ", coords_path, " (n_sites=", nrow(coords), ")\n", sep = "")
cat("[build_site_coordinates] columns used: site_col=", site_col, " | lon_col=", lon_col, " | lat_col=", lat_col, "\n", sep = "")
