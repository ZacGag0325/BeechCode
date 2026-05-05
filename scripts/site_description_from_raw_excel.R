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
# Constants
# -----------------------------
fallback_plot_area_m2 <- 576 # 24m x 24m

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
  
  # examples: 0-1, 1_2, ]2-4], 4+, >=4
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

beech_pattern <- "fagus\\s*grandifolia|^heg$|hetre|h[êe]tre|american\\s*beech"

# -----------------------------
# Paths and workbook loading
# -----------------------------
project_root <- getwd()
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

message("Detected columns in 'arbre': ", paste(names(arbre_raw), collapse = ", "))
if (!is.null(gaule_raw)) message("Detected columns in 'gaule': ", paste(names(gaule_raw), collapse = ", "))
if (!is.null(genetique_raw)) message("Detected columns in 'genetique': ", paste(names(genetique_raw), collapse = ", "))

# -----------------------------
# Column detection
# -----------------------------
site_patterns <- c("^site$", "site_?id", "id_?site", "placette", "parcelle", "plot", "station", "localite", "location")
species_patterns <- c("espece", "species", "taxon", "essence", "code_?ess", "^sp$")
dbh_patterns <- c("^dbh$", "dhp", "dhp_cm", "diam", "diametre", "diameter", "d130")
plot_area_patterns <- c("plot_?area", "surface_?plac", "superficie", "area_?m2")
lat_patterns <- c("^lat", "latitude", "y_coord", "coord_?y")
lon_patterns <- c("^lon", "long", "longitude", "x_coord", "coord_?x")
elev_patterns <- c("elev", "alt", "altitude", "hauteur")
count_patterns <- c("^n$", "count", "nombre", "effectif", "nb", "abond")

arbre_site_col <- pick_column(arbre_raw, site_patterns, "arbre site", required = TRUE)
arbre_species_col <- pick_column(arbre_raw, species_patterns, "arbre species", required = TRUE)
arbre_dbh_col <- pick_column(arbre_raw, dbh_patterns, "arbre dbh", required = TRUE)
arbre_plot_area_col <- pick_column(arbre_raw, plot_area_patterns, "arbre plot area")

if (any(is.na(c(arbre_site_col, arbre_species_col, arbre_dbh_col)))) {
  stop("Could not detect required columns in 'arbre' (site/species/dbh).", call. = FALSE)
}

# -----------------------------
# Arbre metrics
# -----------------------------
arbre <- arbre_raw |>
  mutate(
    site = str_trim(as.character(.data[[arbre_site_col]])),
    species = str_trim(as.character(.data[[arbre_species_col]])),
    species_norm = normalize_ascii(species),
    dbh_cm = coerce_numeric(.data[[arbre_dbh_col]]),
    plot_area_m2 = if (!is.na(arbre_plot_area_col)) coerce_numeric(.data[[arbre_plot_area_col]]) else NA_real_
  ) |>
  mutate(
    plot_area_m2 = ifelse(is.na(plot_area_m2) | plot_area_m2 <= 0, fallback_plot_area_m2, plot_area_m2),
    ba_m2_tree = pi * (dbh_cm / 200)^2
  ) |>
  filter(!is.na(site) & site != "")

if (is.na(arbre_plot_area_col)) {
  warning("No plot area column found in 'arbre'; fallback area of 576 m² used.", call. = FALSE)
}

arbre_site_summary <- arbre |>
  group_by(site) |>
  summarise(
    Site = first(site),
    Basal_Area = sum(ba_m2_tree * (10000 / plot_area_m2), na.rm = TRUE),
    Beech_Basal_Area = sum(ifelse(str_detect(species_norm, beech_pattern), ba_m2_tree * (10000 / plot_area_m2), 0), na.rm = TRUE),
    Mean_DBH = ifelse(all(is.na(dbh_cm)), NA_real_, mean(dbh_cm, na.rm = TRUE)),
    Mean_Beech_DBH = {
      beech_dbh <- dbh_cm[str_detect(species_norm, beech_pattern)]
      if (length(beech_dbh) == 0 || all(is.na(beech_dbh))) NA_real_ else mean(beech_dbh, na.rm = TRUE)
    },
    .groups = "drop"
  )

dominant_species <- arbre |>
  mutate(ba_ha = ba_m2_tree * (10000 / plot_area_m2)) |>
  group_by(site, species) |>
  summarise(ba_ha = sum(ba_ha, na.rm = TRUE), .groups = "drop") |>
  group_by(site) |>
  slice_max(order_by = ba_ha, n = 1, with_ties = FALSE) |>
  ungroup() |>
  transmute(Site = site, Dominant_Species = species)

arbre_site_summary <- arbre_site_summary |>
  left_join(dominant_species, by = "Site")

# -----------------------------
# Gaule metrics (Beech Sapling Basal Area)
# -----------------------------
if (!is.null(gaule_raw)) {
  gaule_site_col <- pick_column(gaule_raw, site_patterns, "gaule site")
  gaule_species_col <- pick_column(gaule_raw, species_patterns, "gaule species")
  gaule_dbh_col <- pick_column(gaule_raw, dbh_patterns, "gaule dbh")
  gaule_count_col <- pick_column(gaule_raw, count_patterns, "gaule count")
  gaule_plot_area_col <- pick_column(gaule_raw, plot_area_patterns, "gaule plot area")
  
  if (is.na(gaule_site_col) || is.na(gaule_species_col) || (is.na(gaule_dbh_col) && is.na(gaule_count_col))) {
    warning("Could not detect enough columns in 'gaule' to compute Beech Sapling Basal Area; keeping NA.", call. = FALSE)
    gaule_summary <- tibble(Site = unique(arbre_site_summary$Site), Beech_Sapling_Basal_Area = NA_real_)
  } else {
    gaule <- gaule_raw |>
      mutate(
        site = str_trim(as.character(.data[[gaule_site_col]])),
        species = str_trim(as.character(.data[[gaule_species_col]])),
        species_norm = normalize_ascii(species),
        dbh_exact = if (!is.na(gaule_dbh_col)) coerce_numeric(.data[[gaule_dbh_col]]) else NA_real_,
        dbh_mid = if (!is.na(gaule_dbh_col)) extract_dbh_midpoint(.data[[gaule_dbh_col]]) else NA_real_,
        dbh_cm = ifelse(!is.na(dbh_exact), dbh_exact, dbh_mid),
        n_stems = if (!is.na(gaule_count_col)) coerce_numeric(.data[[gaule_count_col]]) else 1,
        plot_area_m2 = if (!is.na(gaule_plot_area_col)) coerce_numeric(.data[[gaule_plot_area_col]]) else fallback_plot_area_m2
      ) |>
      mutate(
        n_stems = ifelse(is.na(n_stems) | n_stems <= 0, 1, n_stems),
        plot_area_m2 = ifelse(is.na(plot_area_m2) | plot_area_m2 <= 0, fallback_plot_area_m2, plot_area_m2)
      ) |>
      filter(!is.na(site) & site != "", str_detect(species_norm, beech_pattern))
    
    if (all(is.na(gaule$dbh_cm))) {
      warning("Could not derive DBH (exact or class midpoint) in 'gaule'; keeping Beech Sapling Basal Area as NA.", call. = FALSE)
      gaule_summary <- tibble(Site = unique(arbre_site_summary$Site), Beech_Sapling_Basal_Area = NA_real_)
    } else {
      gaule_summary <- gaule |>
        mutate(ba_m2 = pi * (dbh_cm / 200)^2,
               ba_ha = ba_m2 * n_stems * (10000 / plot_area_m2)) |>
        group_by(site) |>
        summarise(Beech_Sapling_Basal_Area = sum(ba_ha, na.rm = TRUE), .groups = "drop") |>
        transmute(Site = site, Beech_Sapling_Basal_Area)
    }
  }
} else {
  warning("Sheet 'gaule' unavailable; Beech Sapling Basal Area set to NA.", call. = FALSE)
  gaule_summary <- tibble(Site = unique(arbre_site_summary$Site), Beech_Sapling_Basal_Area = NA_real_)
}

# -----------------------------
# Genetique metrics (Elevation, Latitude ordering)
# -----------------------------
latlon_elev_summary <- NULL
if (!is.null(genetique_raw)) {
  gen_site_col <- pick_column(genetique_raw, site_patterns, "genetique site")
  gen_lat_col <- pick_column(genetique_raw, lat_patterns, "genetique latitude")
  gen_lon_col <- pick_column(genetique_raw, lon_patterns, "genetique longitude")
  gen_elev_col <- pick_column(genetique_raw, elev_patterns, "genetique elevation")
  
  if (!is.na(gen_site_col)) {
    gen <- genetique_raw |>
      mutate(
        site = str_trim(as.character(.data[[gen_site_col]])),
        latitude = if (!is.na(gen_lat_col)) coerce_numeric(.data[[gen_lat_col]]) else NA_real_,
        longitude = if (!is.na(gen_lon_col)) coerce_numeric(.data[[gen_lon_col]]) else NA_real_,
        elevation = if (!is.na(gen_elev_col)) coerce_numeric(.data[[gen_elev_col]]) else NA_real_
      ) |>
      filter(!is.na(site) & site != "")
    
    latlon_elev_summary <- gen |>
      group_by(site) |>
      summarise(
        mean_latitude = ifelse(all(is.na(latitude)), NA_real_, mean(latitude, na.rm = TRUE)),
        mean_longitude = ifelse(all(is.na(longitude)), NA_real_, mean(longitude, na.rm = TRUE)),
        Elevation = ifelse(all(is.na(elevation)), NA_real_, mean(elevation, na.rm = TRUE)),
        .groups = "drop"
      ) |>
      transmute(Site = site, mean_latitude, mean_longitude, Elevation)
  } else {
    warning("No site column detected in 'genetique'; cannot compute elevation/latitude from this sheet.", call. = FALSE)
  }
}

# latitude fallback from arbre/gaule if genetique latitude unavailable
if (is.null(latlon_elev_summary) || all(is.na(latlon_elev_summary$mean_latitude))) {
  lat_source <- NULL
  lat_from_arbre <- pick_column(arbre_raw, lat_patterns, "arbre latitude")
  lon_from_arbre <- pick_column(arbre_raw, lon_patterns, "arbre longitude")
  if (!is.na(lat_from_arbre)) {
    lat_source <- arbre_raw |>
      mutate(
        Site = str_trim(as.character(.data[[arbre_site_col]])),
        mean_latitude = coerce_numeric(.data[[lat_from_arbre]]),
        mean_longitude = if (!is.na(lon_from_arbre)) coerce_numeric(.data[[lon_from_arbre]]) else NA_real_
      ) |>
      filter(!is.na(Site) & Site != "") |>
      group_by(Site) |>
      summarise(
        mean_latitude = ifelse(all(is.na(mean_latitude)), NA_real_, mean(mean_latitude, na.rm = TRUE)),
        mean_longitude = ifelse(all(is.na(mean_longitude)), NA_real_, mean(mean_longitude, na.rm = TRUE)),
        .groups = "drop"
      )
  }
  
  if (!is.null(lat_source)) {
    message("Using latitude fallback from 'arbre'.")
    if (is.null(latlon_elev_summary)) {
      latlon_elev_summary <- lat_source |> mutate(Elevation = NA_real_)
    } else {
      latlon_elev_summary <- latlon_elev_summary |>
        select(-mean_latitude, -mean_longitude) |>
        left_join(lat_source, by = "Site")
    }
  } else {
    warning("No latitude column found in genetique or other raw sheets; ordering will fall back to Site name.", call. = FALSE)
    if (is.null(latlon_elev_summary)) {
      latlon_elev_summary <- tibble(Site = unique(arbre_site_summary$Site), mean_latitude = NA_real_, mean_longitude = NA_real_, Elevation = NA_real_)
    }
  }
}

# -----------------------------
# Assemble final table
# -----------------------------
final_tbl <- arbre_site_summary |>
  left_join(gaule_summary, by = "Site") |>
  left_join(latlon_elev_summary, by = "Site") |>
  mutate(
    Site = as.character(Site)
  ) |>
  distinct(Site, .keep_all = TRUE)

if (any(!is.na(final_tbl$mean_latitude))) {
  final_tbl <- final_tbl |> arrange(mean_latitude)
} else {
  final_tbl <- final_tbl |> arrange(Site)
}

final_tbl <- final_tbl |>
  transmute(
    Site = Site,
    `Basal Area` = Basal_Area,
    `Beech Basal Area` = Beech_Basal_Area,
    `Mean DBH` = Mean_DBH,
    `Mean Beech DBH` = Mean_Beech_DBH,
    `Beech Sapling Basal Area` = Beech_Sapling_Basal_Area,
    `Dominant Species` = Dominant_Species,
    Elevation = Elevation
  )

# -----------------------------
# Quality checks
# -----------------------------
if (!is.null(latlon_elev_summary) && any(!is.na(latlon_elev_summary$mean_longitude))) {
  pos_long <- latlon_elev_summary |> filter(!is.na(mean_longitude) & mean_longitude > 0)
  if (nrow(pos_long) > 0) {
    warning("Longitude check: some longitudes are positive, unexpected for Québec. Sites: ", paste(pos_long$Site, collapse = ", "), call. = FALSE)
  }
}

if (nrow(final_tbl) != dplyr::n_distinct(final_tbl$Site)) {
  warning("Table does not have exactly one row per site.", call. = FALSE)
}

neg_ba_sites <- final_tbl |> filter(`Basal Area` < 0 | `Beech Basal Area` < 0 | `Beech Sapling Basal Area` < 0)
if (nrow(neg_ba_sites) > 0) {
  warning("Negative basal area values detected for site(s): ", paste(neg_ba_sites$Site, collapse = ", "), call. = FALSE)
}

beech_presence <- arbre |>
  group_by(site) |>
  summarise(has_beech = any(str_detect(species_norm, beech_pattern), na.rm = TRUE), .groups = "drop") |>
  transmute(Site = site, has_beech)

check_beech_dbh <- final_tbl |>
  left_join(beech_presence, by = "Site") |>
  filter(is.na(`Mean Beech DBH`) & isTRUE(has_beech))
if (nrow(check_beech_dbh) > 0) {
  warning("Mean Beech DBH is NA for site(s) where beech is present in arbre: ", paste(check_beech_dbh$Site, collapse = ", "), call. = FALSE)
}

if (!is.null(gaule_raw)) {
  if (any(is.na(final_tbl$`Beech Sapling Basal Area`))) {
    warning("Beech Sapling Basal Area has NA values. This should occur only when gaule columns were unavailable or unusable.", call. = FALSE)
  }
}

# -----------------------------
# Save outputs
# -----------------------------
csv_path <- file.path(out_dir, "site_description_table_article.csv")
xlsx_path <- file.path(out_dir, "site_description_table_article.xlsx")

write.csv(final_tbl, csv_path, row.names = FALSE, na = "")
openxlsx::write.xlsx(final_tbl, xlsx_path, overwrite = TRUE)

message("\nFinal table:")
print(final_tbl, n = Inf, width = Inf)

# Summary block
range_or_na <- function(x) {
  if (all(is.na(x))) return("NA to NA")
  paste0(round(min(x, na.rm = TRUE), 4), " to ", round(max(x, na.rm = TRUE), 4))
}

message("\nSummary:")
message(" - Number of sites included: ", nrow(final_tbl))
message(" - Basal area range: ", range_or_na(final_tbl$`Basal Area`))
message(" - Beech basal area range: ", range_or_na(final_tbl$`Beech Basal Area`))
message(" - Mean DBH range: ", range_or_na(final_tbl$`Mean DBH`))
message(" - Mean beech DBH range: ", range_or_na(final_tbl$`Mean Beech DBH`))
message(" - Beech sapling basal area range: ", range_or_na(final_tbl$`Beech Sapling Basal Area`))
message(" - Number of sites with missing elevation: ", sum(is.na(final_tbl$Elevation)))

message("\nSaved:")
message(" - ", csv_path)
message(" - ", xlsx_path)