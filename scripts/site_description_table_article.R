#!/usr/bin/env Rscript

required_packages <- c("dplyr", "readxl", "openxlsx", "stringr", "tidyr", "purrr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(
    paste0("Missing required package(s): ", paste(missing_packages, collapse = ", "), ". Please install them before running this script."),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(readxl)
  library(openxlsx)
  library(stringr)
  library(tidyr)
  library(purrr)
})

# -----------------------------
# User-setting required for prism BA (MUST be set)
# -----------------------------
prism_baf <- 2  # prism basal area factor (m²/ha). Facteur 2 prism => prism_baf <- 2.
# Basal area per prism point = n included trees (Valeur_prisme == 2) * prism_baf.
# Site-level basal area = average of prism-point basal areas across Point_prisme values (typically 3 points/site).

if (!is.finite(prism_baf) || prism_baf <= 0) {
  stop("`prism_baf` must be a positive finite value.", call. = FALSE)
}

# -----------------------------
# Helpers
# -----------------------------
normalize_ascii <- function(x) {
  x |>
    iconv(from = "", to = "ASCII//TRANSLIT") |>
    tolower() |>
    str_replace_all("[^a-z0-9]+", "_") |>
    str_replace_all("_+", "_") |>
    str_replace_all("^_|_$", "")
}

coerce_numeric <- function(x) {
  if (is.numeric(x)) return(x)
  x_chr <- as.character(x)
  x_chr <- str_replace_all(x_chr, ",", ".")
  x_chr <- str_replace_all(x_chr, "[^0-9.\\-]+", "")
  suppressWarnings(as.numeric(x_chr))
}

pick_column <- function(df, patterns, label, required = FALSE) {
  cn <- names(df)
  cn_norm <- normalize_ascii(cn)
  idx <- which(str_detect(cn_norm, str_c(patterns, collapse = "|")))
  if (length(idx) == 0) {
    if (required) warning("Could not detect required column for ", label, call. = FALSE)
    return(NA_character_)
  }
  if (length(idx) > 1) {
    message("Multiple matches for ", label, ": ", paste(cn[idx], collapse = ", "), ". Using: ", cn[idx[1]])
  }
  cn[idx[1]]
}

find_raw_workbook <- function(raw_dir) {
  if (!dir.exists(raw_dir)) stop("Directory not found: ", raw_dir, call. = FALSE)
  files <- list.files(raw_dir, pattern = "\\.(xlsx|xls)$", full.names = TRUE, ignore.case = TRUE)
  if (length(files) == 0) stop("No Excel file found in data/raw.", call. = FALSE)
  
  norm_names <- normalize_ascii(basename(files))
  preferred_idx <- which(str_detect(norm_names, "donnee|modifie|west|summer|copie"))
  chosen <- if (length(preferred_idx) > 0) files[preferred_idx[1]] else files[1]
  
  message("Excel files found in data/raw:")
  message(paste0(" - ", basename(files), collapse = "\n"))
  message("Selected workbook: ", basename(chosen))
  chosen
}

extract_dbh_midpoint <- function(x) {
  x_chr <- as.character(x)
  x_chr <- str_replace_all(x_chr, ",", ".")
  x_chr <- str_replace_all(x_chr, "\\s+", "")
  
  ranges <- str_match(x_chr, "(-?[0-9]+\\.?[0-9]*)[^0-9]+(-?[0-9]+\\.?[0-9]*)")
  low <- suppressWarnings(as.numeric(ranges[, 2]))
  high <- suppressWarnings(as.numeric(ranges[, 3]))
  midpoint <- ifelse(!is.na(low) & !is.na(high), (low + high) / 2, NA_real_)
  
  plus_single <- str_match(x_chr, "([0-9]+\\.?[0-9]*)\\+")
  plus_val <- suppressWarnings(as.numeric(plus_single[, 2]))
  midpoint <- ifelse(is.na(midpoint) & !is.na(plus_val), plus_val + 0.5, midpoint)
  
  ge_single <- str_match(x_chr, ">=?([0-9]+\\.?[0-9]*)")
  ge_val <- suppressWarnings(as.numeric(ge_single[, 2]))
  midpoint <- ifelse(is.na(midpoint) & !is.na(ge_val), ge_val + 0.5, midpoint)
  
  midpoint
}

resolve_project_root <- function() {
  script_path <- tryCatch(normalizePath(sys.frames()[[1]]$ofile, winslash = "/", mustWork = FALSE), error = function(e) NA_character_)
  if (is.na(script_path) || script_path == "") {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- args[grepl("^--file=", args)]
    if (length(file_arg) > 0) script_path <- normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = FALSE)
  }
  candidates <- unique(normalizePath(c(dirname(script_path), getwd()), winslash = "/", mustWork = FALSE))
  for (cand in candidates) {
    probe <- cand
    for (i in seq_len(8)) {
      if (dir.exists(file.path(probe, "scripts")) && dir.exists(file.path(probe, "data", "raw"))) return(probe)
      parent <- dirname(probe)
      if (identical(parent, probe)) break
      probe <- parent
    }
  }
  getwd()
}

beech_pattern <- "fagus\\s*grandifolia|^heg$|hetre|h[êe]tre|american\\s*beech"

# -----------------------------
# Paths and workbook loading
# -----------------------------
project_root <- resolve_project_root()
raw_dir <- file.path(project_root, "data", "raw")
out_dir <- file.path(project_root, "outputs", "tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

excel_path <- find_raw_workbook(raw_dir)
sheets <- excel_sheets(excel_path)
message("Workbook sheets: ", paste(sheets, collapse = ", "))

read_sheet_safe <- function(sheet_name) {
  if (!(sheet_name %in% sheets)) {
    warning("Sheet '", sheet_name, "' not found.", call. = FALSE)
    return(NULL)
  }
  read_excel(excel_path, sheet = sheet_name)
}

arbre_raw <- read_sheet_safe("arbre")
gaule_raw <- read_sheet_safe("gaule")
genetique_raw <- read_sheet_safe("genetique")

if (is.null(arbre_raw)) stop("Sheet 'arbre' is required.", call. = FALSE)
if (is.null(gaule_raw)) stop("Sheet 'gaule' is required for this article table workflow.", call. = FALSE)
if (is.null(genetique_raw)) stop("Sheet 'genetique' is required for south-to-north ordering and elevation.", call. = FALSE)

message("Detected columns in 'arbre': ", paste(names(arbre_raw), collapse = ", "))
message("Detected columns in 'gaule': ", paste(names(gaule_raw), collapse = ", "))
message("Detected columns in 'genetique': ", paste(names(genetique_raw), collapse = ", "))

# -----------------------------
# Column detection
# -----------------------------
site_patterns <- c("^site$", "site_?id", "id_?site", "placette", "parcelle", "plot", "station", "localite", "location")
species_patterns <- c("espece", "species", "taxon", "essence", "code_?ess", "^sp$")
dbh_patterns <- c("^dbh$", "dhp", "dhp_cm", "diam", "diametre", "diameter", "d130", "dhp_tige")
point_patterns <- c("point_?prisme", "prism_?point", "point")
val_prism_patterns <- c("valeur_?prisme", "inclusion", "inclu", "prism_?value")
lat_patterns <- c("^lat", "latitude", "y_coord", "coord_?y")
lon_patterns <- c("^lon", "long", "longitude", "x_coord", "coord_?x")
elev_patterns <- c("elev", "alt", "altitude", "hauteur")
count_patterns <- c("^n$", "count", "nombre", "effectif", "nb", "abond")

arbre_site_col <- pick_column(arbre_raw, site_patterns, "arbre site", required = TRUE)
arbre_species_col <- pick_column(arbre_raw, species_patterns, "arbre species", required = TRUE)
arbre_dbh_col <- pick_column(arbre_raw, dbh_patterns, "arbre dbh", required = TRUE)
arbre_point_col <- pick_column(arbre_raw, point_patterns, "arbre point_prisme", required = TRUE)
arbre_val_prism_col <- pick_column(arbre_raw, val_prism_patterns, "arbre valeur_prisme", required = TRUE)

if (any(is.na(c(arbre_site_col, arbre_species_col, arbre_dbh_col, arbre_point_col, arbre_val_prism_col)))) {
  stop("Could not detect required arbre columns (site/species/dbh/point_prisme/valeur_prisme).", call. = FALSE)
}

# -----------------------------
# Arbre: included prism trees only (Valeur_prisme == 2)
# -----------------------------
arbre_all <- arbre_raw |>
  mutate(
    Site = str_trim(as.character(.data[[arbre_site_col]])),
    Point_prisme = str_trim(as.character(.data[[arbre_point_col]])),
    species_raw = str_trim(as.character(.data[[arbre_species_col]])),
    species_norm = normalize_ascii(species_raw),
    dbh_cm = coerce_numeric(.data[[arbre_dbh_col]]),
    valeur_prisme = coerce_numeric(.data[[arbre_val_prism_col]])
  ) |>
  filter(!is.na(Site) & Site != "", !is.na(Point_prisme) & Point_prisme != "")

arbre_inc <- arbre_all |>
  filter(!is.na(valeur_prisme) & valeur_prisme == 2)

if (nrow(arbre_inc) == 0) {
  stop("No included prism trees found with Valeur_prisme == 2 in 'arbre'.", call. = FALSE)
}

message("Beech detection pattern used: ", beech_pattern)
message("Mean DBH metrics are calculated using included prism trees only (Valeur_prisme == 2).")

# per point prism BA
point_totals <- arbre_inc |>
  group_by(Site, Point_prisme) |>
  summarise(n_included_total = n(), total_basal_area_by_point = n_included_total * prism_baf, .groups = "drop")

point_beech <- arbre_inc |>
  filter(str_detect(species_norm, beech_pattern)) |>
  group_by(Site, Point_prisme) |>
  summarise(n_included_beech = n(), beech_basal_area_by_point = n_included_beech * prism_baf, .groups = "drop")

point_species <- arbre_inc |>
  group_by(Site, Point_prisme, species_raw) |>
  summarise(n_included_species = n(), species_ba_point = n_included_species * prism_baf, .groups = "drop")

# full point grid from included trees
site_points <- arbre_inc |>
  distinct(Site, Point_prisme)

point_diag <- site_points |>
  left_join(point_totals, by = c("Site", "Point_prisme")) |>
  left_join(point_beech, by = c("Site", "Point_prisme")) |>
  mutate(
    n_included_total = coalesce(n_included_total, 0L),
    total_basal_area_by_point = coalesce(total_basal_area_by_point, 0),
    n_included_beech = coalesce(n_included_beech, 0L),
    beech_basal_area_by_point = coalesce(beech_basal_area_by_point, 0)
  )

site_ba <- point_diag |>
  group_by(Site) |>
  summarise(
    n_prism_points_detected = n_distinct(Point_prisme),
    n_included_arbre_total = sum(n_included_total, na.rm = TRUE),
    n_included_arbre_beech = sum(n_included_beech, na.rm = TRUE),
    Basal_Area = mean(total_basal_area_by_point, na.rm = TRUE),
    Beech_Basal_Area = mean(beech_basal_area_by_point, na.rm = TRUE),
    n_included_by_point = paste0(Point_prisme, ":", n_included_total, collapse = " | "),
    total_basal_area_by_point = paste0(Point_prisme, ":", round(total_basal_area_by_point, 4), collapse = " | "),
    beech_basal_area_by_point = paste0(Point_prisme, ":", round(beech_basal_area_by_point, 4), collapse = " | "),
    .groups = "drop"
  )

# mean DBH on included only
dbh_site <- arbre_inc |>
  group_by(Site) |>
  summarise(
    Mean_DBH = ifelse(all(is.na(dbh_cm)), NA_real_, mean(dbh_cm, na.rm = TRUE)),
    Mean_Beech_DBH = {
      b <- dbh_cm[str_detect(species_norm, beech_pattern)]
      if (length(b) == 0 || all(is.na(b))) NA_real_ else mean(b, na.rm = TRUE)
    },
    n_beech_for_mean_dbh = sum(str_detect(species_norm, beech_pattern) & !is.na(dbh_cm), na.rm = TRUE),
    .groups = "drop"
  )

# dominant species by average prism BA across points
species_site <- point_species |>
  group_by(Site, species_raw) |>
  summarise(avg_species_ba = mean(species_ba_point, na.rm = TRUE), .groups = "drop")

dom_ties <- species_site |>
  group_by(Site) |>
  filter(avg_species_ba == max(avg_species_ba, na.rm = TRUE)) |>
  summarise(n_top = n(), top_species = paste(species_raw, collapse = ", "), .groups = "drop") |>
  filter(n_top > 1)
if (nrow(dom_ties) > 0) {
  warning(
    "Dominant species ties detected; first species after sorting kept. Sites: ",
    paste0(dom_ties$Site, " [", dom_ties$top_species, "]", collapse = "; "),
    call. = FALSE
  )
}

dominant_species <- species_site |>
  arrange(Site, desc(avg_species_ba), species_raw) |>
  group_by(Site) |>
  slice(1) |>
  ungroup() |>
  transmute(Site, Dominant_Species = species_raw)

# -----------------------------
# Gaule: method check + computation
# -----------------------------
gaule_site_col <- pick_column(gaule_raw, site_patterns, "gaule site", required = TRUE)
gaule_species_col <- pick_column(gaule_raw, species_patterns, "gaule species", required = TRUE)
gaule_dbh_col <- pick_column(gaule_raw, dbh_patterns, "gaule dbh")
gaule_count_col <- pick_column(gaule_raw, count_patterns, "gaule count")
gaule_point_col <- pick_column(gaule_raw, point_patterns, "gaule point_prisme")
gaule_r_utilise_col <- pick_column(gaule_raw, c("r_?utilise", "rayon", "radius"), "gaule r_utilise")

if (any(is.na(c(gaule_site_col, gaule_species_col)))) {
  stop("Cannot compute Beech Sapling Basal Area: missing site/species columns in gaule.", call. = FALSE)
}

if (is.na(gaule_dbh_col)) {
  stop("Cannot compute Beech Sapling Basal Area: no DBH or DBH class column detected in gaule.", call. = FALSE)
}

message("Gaule method: treating gaule as non-prism sapling observations; computing individual basal area from DBH (or DBH midpoint) and summing by site (no prism BAF applied).")
message("If this assumption is not correct for your field protocol, STOP and adjust this section before using outputs.")

gaule <- gaule_raw |>
  mutate(
    Site = str_trim(as.character(.data[[gaule_site_col]])),
    species_raw = str_trim(as.character(.data[[gaule_species_col]])),
    species_norm = normalize_ascii(species_raw),
    dbh_exact = coerce_numeric(.data[[gaule_dbh_col]]),
    dbh_mid = extract_dbh_midpoint(.data[[gaule_dbh_col]]),
    dbh_cm = ifelse(!is.na(dbh_exact), dbh_exact, dbh_mid),
    n_stems = if (!is.na(gaule_count_col)) coerce_numeric(.data[[gaule_count_col]]) else 1
  ) |>
  mutate(n_stems = ifelse(is.na(n_stems) | n_stems <= 0, 1, n_stems)) |>
  filter(!is.na(Site) & Site != "", str_detect(species_norm, beech_pattern))

if (nrow(gaule) == 0) {
  warning("No beech sapling rows detected in gaule for the selected beech code/pattern.", call. = FALSE)
}

if (all(is.na(gaule$dbh_cm)) && nrow(gaule) > 0) {
  stop("Cannot compute Beech Sapling Basal Area: gaule has beech rows but DBH values/classes could not be interpreted.", call. = FALSE)
}

gaule_summary <- gaule |>
  mutate(ba_m2 = pi * (dbh_cm / 200)^2 * n_stems) |>
  group_by(Site) |>
  summarise(
    Beech_Sapling_Basal_Area = ifelse(all(is.na(ba_m2)), NA_real_, sum(ba_m2, na.rm = TRUE)),
    n_gaule_beech = n(),
    .groups = "drop"
  )

# -----------------------------
# Genetique: elevation + latitude for ordering
# -----------------------------
gen_site_col <- pick_column(genetique_raw, site_patterns, "genetique site", required = TRUE)
gen_lat_col <- pick_column(genetique_raw, lat_patterns, "genetique latitude", required = TRUE)
gen_lon_col <- pick_column(genetique_raw, lon_patterns, "genetique longitude")
gen_elev_col <- pick_column(genetique_raw, elev_patterns, "genetique elevation", required = TRUE)

if (any(is.na(c(gen_site_col, gen_lat_col, gen_elev_col)))) {
  stop("genetique must provide site, latitude, and elevation columns for this workflow.", call. = FALSE)
}

gen <- genetique_raw |>
  mutate(
    Site = str_trim(as.character(.data[[gen_site_col]])),
    latitude = coerce_numeric(.data[[gen_lat_col]]),
    longitude = if (!is.na(gen_lon_col)) coerce_numeric(.data[[gen_lon_col]]) else NA_real_,
    elevation = coerce_numeric(.data[[gen_elev_col]])
  ) |>
  filter(!is.na(Site) & Site != "")

latlon_elev <- gen |>
  group_by(Site) |>
  summarise(
    latitude_used_for_sorting = ifelse(all(is.na(latitude)), NA_real_, mean(latitude, na.rm = TRUE)),
    mean_longitude = ifelse(all(is.na(longitude)), NA_real_, mean(longitude, na.rm = TRUE)),
    Elevation = ifelse(all(is.na(elevation)), NA_real_, mean(elevation, na.rm = TRUE)),
    n_elevation_records = sum(!is.na(elevation)),
    .groups = "drop"
  )

if (any(!is.na(latlon_elev$mean_longitude) & latlon_elev$mean_longitude > 0)) {
  bad <- latlon_elev |> filter(!is.na(mean_longitude) & mean_longitude > 0)
  warning("Longitude check: positive longitudes found (unexpected for Québec): ", paste(bad$Site, collapse = ", "), call. = FALSE)
}

# -----------------------------
# Assemble final table
# -----------------------------
final_tbl <- site_ba |>
  left_join(dbh_site, by = "Site") |>
  left_join(gaule_summary |> select(Site, Beech_Sapling_Basal_Area), by = "Site") |>
  left_join(dominant_species, by = "Site") |>
  left_join(latlon_elev |> select(Site, latitude_used_for_sorting, Elevation), by = "Site")

if (any(is.na(final_tbl$latitude_used_for_sorting))) {
  stop("Cannot sort south-to-north: missing latitude for one or more sites in genetique.", call. = FALSE)
}

final_tbl <- final_tbl |>
  arrange(latitude_used_for_sorting) |>
  transmute(
    Site,
    `Basal Area` = Basal_Area,
    `Beech Basal Area` = Beech_Basal_Area,
    `Mean DBH` = Mean_DBH,
    `Mean Beech DBH` = Mean_Beech_DBH,
    `Beech Sapling Basal Area` = Beech_Sapling_Basal_Area,
    `Dominant Species` = Dominant_Species,
    Elevation
  )

if (nrow(final_tbl) != dplyr::n_distinct(final_tbl$Site)) {
  stop("Final table does not have exactly one row per site.", call. = FALSE)
}
if (any(final_tbl$`Basal Area` < 0, na.rm = TRUE) || any(final_tbl$`Beech Basal Area` < 0, na.rm = TRUE)) {
  stop("Negative prism basal area detected in final table.", call. = FALSE)
}

# -----------------------------
# Diagnostics table
# -----------------------------
diagnostics_tbl <- site_ba |>
  left_join(dbh_site |> select(Site, n_beech_for_mean_dbh), by = "Site") |>
  left_join(gaule_summary |> select(Site, n_gaule_beech), by = "Site") |>
  left_join(latlon_elev |> select(Site, n_elevation_records, latitude_used_for_sorting), by = "Site") |>
  mutate(
    n_gaule_beech = coalesce(n_gaule_beech, 0L),
    basal_area_method_used = "Prism BA = mean across Point_prisme of (n included trees where Valeur_prisme==2 * prism_baf)",
    prism_baf_used = prism_baf
  ) |>
  select(
    Site,
    n_prism_points_detected,
    n_included_arbre_total,
    n_included_arbre_beech,
    n_included_by_point,
    total_basal_area_by_point,
    beech_basal_area_by_point,
    n_beech_for_mean_dbh,
    n_gaule_beech,
    n_elevation_records,
    latitude_used_for_sorting,
    basal_area_method_used,
    prism_baf_used
  ) |>
  arrange(latitude_used_for_sorting)

# -----------------------------
# Save outputs
# -----------------------------
csv_path <- file.path(out_dir, "site_description_table_article.csv")
xlsx_path <- file.path(out_dir, "site_description_table_article.xlsx")
diag_path <- file.path(out_dir, "site_description_table_article_diagnostics.csv")

write.csv(final_tbl, csv_path, row.names = FALSE, na = "")
openxlsx::write.xlsx(final_tbl, xlsx_path, overwrite = TRUE)
write.csv(diagnostics_tbl, diag_path, row.names = FALSE, na = "")

# -----------------------------
# Console output
# -----------------------------
message("\nFinal table:")
print(final_tbl, n = Inf, width = Inf)

message("\nDiagnostics table:")
print(diagnostics_tbl, n = Inf, width = Inf)

sites_lt3 <- diagnostics_tbl |> filter(n_prism_points_detected < 3)

range_or_na <- function(x) {
  if (all(is.na(x))) return("NA to NA")
  paste0(round(min(x, na.rm = TRUE), 4), " to ", round(max(x, na.rm = TRUE), 4))
}

message("\nSummary:")
message(" - prism_baf used: ", prism_baf)
message(" - Prism points per site: ", paste0(diagnostics_tbl$Site, "=", diagnostics_tbl$n_prism_points_detected, collapse = "; "))
message(" - Any sites with fewer than 3 prism points: ", ifelse(nrow(sites_lt3) > 0, paste(sites_lt3$Site, collapse = ", "), "No"))
message(" - Sites with Beech Basal Area = 0: ", {s <- final_tbl$Site[!is.na(final_tbl$`Beech Basal Area`) & final_tbl$`Beech Basal Area` == 0]; if (length(s)==0) "None" else paste(s, collapse=", ")})
message(" - Sites with Mean Beech DBH = NA: ", {s <- final_tbl$Site[is.na(final_tbl$`Mean Beech DBH`)]; if (length(s)==0) "None" else paste(s, collapse=", ")})
message(" - Basal Area range: ", range_or_na(final_tbl$`Basal Area`))
message(" - Beech Basal Area range: ", range_or_na(final_tbl$`Beech Basal Area`))

message("\nSaved:")
message(" - ", csv_path)
message(" - ", xlsx_path)
message(" - ", diag_path)