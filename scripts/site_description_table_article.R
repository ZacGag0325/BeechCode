#!/usr/bin/env Rscript

required_packages <- c("dplyr", "readxl", "openxlsx", "stringr", "tidyr", "purrr")
optional_word_packages <- c("flextable", "officer")

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(
    paste0("Missing required package(s): ", paste(missing_packages, collapse = ", "), ". Please install them before running this script."),
    call. = FALSE
  )
}

word_packages_available <- all(vapply(optional_word_packages, requireNamespace, logical(1), quietly = TRUE))
missing_word_packages <- optional_word_packages[!vapply(optional_word_packages, requireNamespace, logical(1), quietly = TRUE)]

suppressPackageStartupMessages({
  library(dplyr)
  library(readxl)
  library(openxlsx)
  library(stringr)
  library(tidyr)
  library(purrr)
})

# -----------------------------
# Required sampling parameter
# -----------------------------
prism_baf <- 2  # mature-tree prism factor (m²/ha), facteur 2

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

pick_column <- function(df, patterns, label, required = FALSE, exclude = character()) {
  cn <- names(df)
  cn_norm <- normalize_ascii(cn)
  exclude_norm <- normalize_ascii(exclude)
  idx <- which(str_detect(cn_norm, str_c(patterns, collapse = "|")) & !(cn_norm %in% exclude_norm))
  if (length(idx) == 0) {
    if (required) warning("Could not detect required column for ", label, call. = FALSE)
    return(NA_character_)
  }
  if (length(idx) > 1) {
    message("Multiple matches for ", label, ": ", paste(cn[idx], collapse = ", "), ". Using: ", cn[idx[1]])
  }
  cn[idx[1]]
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
      if (dir.exists(file.path(probe, "scripts")) && dir.exists(file.path(probe, "data"))) return(probe)
      parent <- dirname(probe)
      if (identical(parent, probe)) break
      probe <- parent
    }
  }
  getwd()
}

find_raw_workbook <- function(raw_dir) {
  if (!dir.exists(raw_dir)) stop("Directory not found: ", raw_dir, call. = FALSE)
  files <- list.files(raw_dir, pattern = "\\.(xlsx|xls)$", full.names = TRUE, ignore.case = TRUE)
  if (length(files) == 0) stop("No Excel file found in data/raw.", call. = FALSE)
  preferred_idx <- which(str_detect(normalize_ascii(basename(files)), "donnee|modifie|west|summer|copie"))
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
  midpoint
}

range_or_na <- function(x, digits = 4) {
  if (length(x) == 0 || all(is.na(x))) return("NA to NA")
  paste0(round(min(x, na.rm = TRUE), digits), " to ", round(max(x, na.rm = TRUE), digits))
}

collapse_or_none <- function(x) {
  x <- unique(x[!is.na(x) & x != ""])
  if (length(x) == 0) "None" else paste(sort(x), collapse = ", ")
}

match_sheet <- function(sheets, target) {
  idx <- which(normalize_ascii(sheets) == normalize_ascii(target))
  if (length(idx) == 0) NA_character_ else sheets[idx[1]]
}

read_sheet_safe <- function(sheet_name) {
  matched_sheet <- match_sheet(sheets, sheet_name)
  if (is.na(matched_sheet)) stop("Sheet '", sheet_name, "' not found.", call. = FALSE)
  read_excel(excel_path, sheet = matched_sheet)
}

find_preferred_column <- function(df, candidates) {
  cn <- names(df)
  cn_norm <- normalize_ascii(cn)
  candidate_norm <- normalize_ascii(candidates)
  idx <- match(candidate_norm, cn_norm, nomatch = 0)
  idx <- idx[idx > 0]
  if (length(idx) == 0) NA_character_ else cn[idx[1]]
}

make_lookup <- function(df, from_patterns, to_patterns, lookup_name) {
  if (identical(lookup_name, "species_lookup")) {
    from_col <- find_preferred_column(df, c("species_code_fr", "code_fr", "espece_fr", "species_fr"))
    to_col <- find_preferred_column(df, c("species_code_en", "code_en", "espece_en", "species_en"))
  } else if (identical(lookup_name, "site_lookup")) {
    from_col <- find_preferred_column(df, c("site_code", "original_site_code", "site_code_original", "original_site", "old_site_code"))
    to_col <- find_preferred_column(df, c("site_label", "final_site_label", "site_final", "final_site", "new_site_label"))
  } else {
    from_col <- NA_character_
    to_col <- NA_character_
  }
  
  if (is.na(from_col)) {
    from_col <- pick_column(df, from_patterns, paste0(lookup_name, " original/code"), required = TRUE)
  }
  if (is.na(to_col)) {
    to_col <- pick_column(df, to_patterns, paste0(lookup_name, " final/label"), required = TRUE, exclude = from_col)
  }
  if (is.na(from_col) || is.na(to_col)) {
    stop("Could not detect required columns in ", lookup_name, ". Detected columns: ", paste(names(df), collapse = ", "), call. = FALSE)
  }
  
  lookup <- df |>
    transmute(
      original_code = str_trim(as.character(.data[[from_col]])),
      converted_code = str_trim(as.character(.data[[to_col]])),
      lookup_key = normalize_ascii(original_code)
    ) |>
    filter(!is.na(original_code) & original_code != "", !is.na(converted_code) & converted_code != "") |>
    distinct(lookup_key, .keep_all = TRUE)
  
  if (nrow(lookup) == 0) stop(lookup_name, " did not contain any usable lookup rows.", call. = FALSE)
  lookup
}

convert_codes <- function(codes, lookup) {
  code_chr <- str_trim(as.character(codes))
  idx <- match(normalize_ascii(code_chr), lookup$lookup_key)
  converted <- lookup$converted_code[idx]
  ifelse(is.na(converted), code_chr, converted)
}

beech_pattern <- "fagus\\s*grandifolia|^heg$|hetre|h[êe]tre|american\\s*beech"

# -----------------------------
# Load workbook and lookup sheets
# -----------------------------
project_root <- resolve_project_root()
raw_dir <- file.path(project_root, "data", "raw")
out_dir <- file.path(project_root, "outputs", "tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

excel_path <- find_raw_workbook(raw_dir)
sheets <- excel_sheets(excel_path)
message("Workbook sheets: ", paste(sheets, collapse = ", "))

site_lookup_sheet <- match_sheet(sheets, "site_lookup")
species_lookup_sheet <- match_sheet(sheets, "species_lookup")
site_lookup_used <- !is.na(site_lookup_sheet)
species_lookup_used <- !is.na(species_lookup_sheet)

if (!site_lookup_used) stop("Required sheet 'site_lookup' was not found in the workbook.", call. = FALSE)
if (!species_lookup_used) stop("Required sheet 'species_lookup' was not found in the workbook.", call. = FALSE)

message("site_lookup sheet detected: ", site_lookup_sheet)
message("species_lookup sheet detected: ", species_lookup_sheet)

site_lookup_raw <- read_excel(excel_path, sheet = site_lookup_sheet)
species_lookup_raw <- read_excel(excel_path, sheet = species_lookup_sheet)

site_lookup <- make_lookup(
  site_lookup_raw,
  from_patterns = c("original.*site", "old.*site", "original", "old", "site.*code", "code.*site", "^site$", "^code$"),
  to_patterns = c("final.*site", "site.*label", "label", "new.*site", "article.*site", "converted.*site", "^site_final$"),
  lookup_name = "site_lookup"
)

species_lookup <- make_lookup(
  species_lookup_raw,
  from_patterns = c("french", "francais", "français", "original", "old", "espece", "species.*code", "code.*species", "^species$", "^code$"),
  to_patterns = c("english", "anglais", "final", "new", "converted", "article", "^en$", "code_en"),
  lookup_name = "species_lookup"
)

arbre_raw <- read_sheet_safe("arbre")
gaule_raw <- read_sheet_safe("gaule")
genetique_raw <- read_sheet_safe("genetique")

message("Detected columns in 'arbre': ", paste(names(arbre_raw), collapse = ", "))
message("Detected columns in 'gaule': ", paste(names(gaule_raw), collapse = ", "))
message("Detected columns in 'genetique': ", paste(names(genetique_raw), collapse = ", "))

site_patterns <- c("^site$", "site_?id", "id_?site", "placette", "parcelle", "plot", "station", "localite", "location")
species_patterns <- c("espece", "species", "taxon", "essence", "code_?ess", "^sp$")
dbh_patterns <- c("^dbh$", "dhp", "dhp_cm", "diam", "diametre", "diameter", "d130", "dhp_tige")
point_patterns <- c("point_?prisme", "prism_?point", "point")
val_prism_patterns <- c("valeur_?prisme", "inclusion", "inclu", "prism_?value")
lat_patterns <- c("^lat", "latitude", "y_coord", "coord_?y")
lon_patterns <- c("^lon", "long", "longitude", "x_coord", "coord_?x")
elev_patterns <- c("elev", "alt", "altitude", "hauteur")
radius_patterns <- c("r_?utilise", "rayon", "radius")
count_patterns <- c("^n$", "count", "nombre", "effectif", "nb", "abond")

# -----------------------------
# Arbre prism calculations
# -----------------------------
message("Mature-tree basal area uses prism sampling with prism_baf <- 2 and Valeur_prisme == 2 inclusion.")

arbre_site_col <- pick_column(arbre_raw, site_patterns, "arbre site", required = TRUE)
arbre_species_col <- pick_column(arbre_raw, species_patterns, "arbre species", required = TRUE)
arbre_dbh_col <- pick_column(arbre_raw, dbh_patterns, "arbre dbh", required = TRUE)
arbre_point_col <- pick_column(arbre_raw, point_patterns, "arbre point_prisme", required = TRUE)
arbre_val_col <- pick_column(arbre_raw, val_prism_patterns, "arbre valeur_prisme", required = TRUE)

arbre_inc <- arbre_raw |>
  mutate(
    Site = str_trim(as.character(.data[[arbre_site_col]])),
    Point_prisme = str_trim(as.character(.data[[arbre_point_col]])),
    species_raw = str_trim(as.character(.data[[arbre_species_col]])),
    species_norm = normalize_ascii(species_raw),
    dbh_cm = coerce_numeric(.data[[arbre_dbh_col]]),
    valeur_prisme = coerce_numeric(.data[[arbre_val_col]])
  ) |>
  filter(!is.na(Site) & Site != "", !is.na(Point_prisme) & Point_prisme != "", valeur_prisme == 2)

if (nrow(arbre_inc) == 0) stop("No included trees found in arbre where Valeur_prisme == 2.", call. = FALSE)

point_totals <- arbre_inc |>
  group_by(Site, Point_prisme) |>
  summarise(n_included_total = n(), total_basal_area_by_point = n_included_total * prism_baf, .groups = "drop")

point_beech <- arbre_inc |>
  filter(str_detect(species_norm, beech_pattern)) |>
  group_by(Site, Point_prisme) |>
  summarise(n_included_beech = n(), beech_basal_area_by_point = n_included_beech * prism_baf, .groups = "drop")

site_points <- arbre_inc |> distinct(Site, Point_prisme)
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
    .groups = "drop"
  )

dbh_site <- arbre_inc |>
  group_by(Site) |>
  summarise(
    Mean_DBH = ifelse(all(is.na(dbh_cm)), NA_real_, mean(dbh_cm, na.rm = TRUE)),
    Mean_Beech_DBH = {
      b <- dbh_cm[str_detect(species_norm, beech_pattern)]
      if (length(b) == 0 || all(is.na(b))) NA_real_ else mean(b, na.rm = TRUE)
    },
    .groups = "drop"
  )

species_present_by_site <- arbre_inc |>
  distinct(Site, species_raw)

species_point <- site_points |>
  inner_join(species_present_by_site, by = "Site", relationship = "many-to-many") |>
  left_join(
    arbre_inc |>
      count(Site, Point_prisme, species_raw, name = "n_included_species"),
    by = c("Site", "Point_prisme", "species_raw")
  ) |>
  mutate(
    n_included_species = coalesce(n_included_species, 0L),
    ba_point = n_included_species * prism_baf
  )

species_site <- species_point |>
  group_by(Site, species_raw) |>
  summarise(avg_species_ba = mean(ba_point, na.rm = TRUE), .groups = "drop")

top_species_ranked <- species_site |>
  group_by(Site) |>
  arrange(desc(avg_species_ba), species_raw, .by_group = TRUE) |>
  mutate(
    species_order = row_number(),
    cutoff_ba = if (n() <= 3) min(avg_species_ba, na.rm = TRUE) else avg_species_ba[3],
    tied_at_cutoff = n() > 3 & avg_species_ba == cutoff_ba & sum(avg_species_ba == cutoff_ba, na.rm = TRUE) > 1
  ) |>
  ungroup() |>
  mutate(species_converted = convert_codes(species_raw, species_lookup))

top_species_by_site <- top_species_ranked |>
  filter(species_order <= 3)

tie_sites <- top_species_ranked |>
  filter(tied_at_cutoff) |>
  distinct(Site) |>
  pull(Site)
if (length(tie_sites) > 0) {
  message(
    "Tie detected at the top-3 dominant-species cutoff. The final table keeps exactly 3 species; tied species were ordered alphabetically for site(s): ",
    paste(sort(tie_sites), collapse = ", ")
  )
}

dominant_species <- top_species_by_site |>
  group_by(Site) |>
  summarise(
    top_3_species_original_codes = paste(species_raw, collapse = ", "),
    top_3_species_converted_codes = paste(species_converted, collapse = ", "),
    top_3_species_basal_area_values = paste(round(avg_species_ba, 6), collapse = ", "),
    `Top 3 Dominant Species` = top_3_species_converted_codes,
    .groups = "drop"
  )

# -----------------------------
# Gaule fixed-radius calculations from R_utilise
# -----------------------------
message("Gaule basal area uses fixed-radius subplot area from R_utilise (not prism).")
message("In gaule, Point_prisme is station ID only. Beech sapling basal area is expressed in m²/ha.")

gaule_site_col <- pick_column(gaule_raw, site_patterns, "gaule site", required = TRUE)
gaule_species_col <- pick_column(gaule_raw, species_patterns, "gaule species", required = TRUE)
gaule_point_col <- pick_column(gaule_raw, point_patterns, "gaule point_prisme", required = TRUE)
gaule_dbh_col <- pick_column(gaule_raw, dbh_patterns, "gaule dbh", required = TRUE)
gaule_r_col <- pick_column(gaule_raw, radius_patterns, "gaule R_utilise", required = TRUE)
gaule_count_col <- pick_column(gaule_raw, count_patterns, "gaule count")

gaule_beech <- gaule_raw |>
  mutate(
    Site = str_trim(as.character(.data[[gaule_site_col]])),
    Point_prisme = str_trim(as.character(.data[[gaule_point_col]])),
    species_norm = normalize_ascii(str_trim(as.character(.data[[gaule_species_col]]))),
    dbh_exact = coerce_numeric(.data[[gaule_dbh_col]]),
    dbh_mid = extract_dbh_midpoint(.data[[gaule_dbh_col]]),
    dbh_cm = ifelse(!is.na(dbh_exact), dbh_exact, dbh_mid),
    r_utilise_m = coerce_numeric(.data[[gaule_r_col]]),
    n_stems = if (!is.na(gaule_count_col)) coerce_numeric(.data[[gaule_count_col]]) else 1
  ) |>
  mutate(n_stems = ifelse(is.na(n_stems) | n_stems <= 0, 1, n_stems)) |>
  filter(!is.na(Site) & Site != "", !is.na(Point_prisme) & Point_prisme != "", str_detect(species_norm, beech_pattern))

if (nrow(gaule_beech) == 0) warning("No beech gaule rows detected.", call. = FALSE)

bad_r <- gaule_beech |>
  filter(is.na(r_utilise_m) | !is.finite(r_utilise_m) | r_utilise_m <= 0) |>
  distinct(Site, Point_prisme)
if (nrow(bad_r) > 0) {
  warning("Missing/invalid R_utilise detected for these Site×Point_prisme: ", paste0(bad_r$Site, "-", bad_r$Point_prisme, collapse = ", "), call. = FALSE)
  stop("Cannot continue: beech gaule rows have missing/invalid R_utilise.", call. = FALSE)
}

r_values <- sort(unique(gaule_beech$r_utilise_m))
message("R_utilise values detected in gaule (m): ", paste(round(r_values, 4), collapse = ", "))

gaule_subplot <- gaule_beech |>
  mutate(
    ba_individual_m2 = pi * (dbh_cm / 200)^2 * n_stems,
    subplot_area_m2 = pi * r_utilise_m^2,
    subplot_area_ha = subplot_area_m2 / 10000
  ) |>
  group_by(Site, Point_prisme) |>
  summarise(
    n_gaule_beech = sum(n_stems, na.rm = TRUE),
    summed_ba_m2 = sum(ba_individual_m2, na.rm = TRUE),
    gaule_radius_values_detected = paste(sort(unique(round(r_utilise_m, 4))), collapse = " / "),
    gaule_subplot_area_ha = first(subplot_area_ha),
    beech_sapling_basal_area_by_subplot = summed_ba_m2 / gaule_subplot_area_ha,
    .groups = "drop"
  )

gaule_site <- gaule_subplot |>
  group_by(Site) |>
  summarise(
    n_gaule_beech = sum(n_gaule_beech, na.rm = TRUE),
    `Beech Sapling Basal Area` = mean(beech_sapling_basal_area_by_subplot, na.rm = TRUE),
    gaule_radius_values_detected = paste(sort(unique(gaule_radius_values_detected)), collapse = " / "),
    gaule_subplot_area_ha = paste(round(sort(unique(gaule_subplot_area_ha)), 8), collapse = " / "),
    beech_sapling_basal_area_by_subplot = paste0(Point_prisme, ":", round(beech_sapling_basal_area_by_subplot, 6), collapse = " | "),
    .groups = "drop"
  )

# -----------------------------
# Genetique for sorting/elevation/latitude
# -----------------------------
gen_site_col <- pick_column(genetique_raw, site_patterns, "genetique site", required = TRUE)
gen_lat_col <- pick_column(genetique_raw, lat_patterns, "genetique latitude", required = TRUE)
gen_lon_col <- pick_column(genetique_raw, lon_patterns, "genetique longitude")
gen_elev_col <- pick_column(genetique_raw, elev_patterns, "genetique elevation", required = TRUE)

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
    Latitude = ifelse(all(is.na(latitude)), NA_real_, round(mean(latitude, na.rm = TRUE), 5)),
    Elevation = ifelse(all(is.na(elevation)), NA_real_, mean(elevation, na.rm = TRUE)),
    n_elevation_records = sum(!is.na(elevation)),
    mean_longitude = ifelse(all(is.na(longitude)), NA_real_, mean(longitude, na.rm = TRUE)),
    .groups = "drop"
  )

if (any(!is.na(latlon_elev$mean_longitude) & latlon_elev$mean_longitude > 0)) {
  warning("Positive longitudes detected for some sites; check coordinate sign.", call. = FALSE)
}

# -----------------------------
# Lookup diagnostics
# -----------------------------
all_site_codes <- unique(c(site_ba$Site, dbh_site$Site, gaule_site$Site, latlon_elev$Site))
site_conversion_tbl <- tibble(Site = all_site_codes) |>
  mutate(
    Site_label = convert_codes(Site, site_lookup),
    site_missing_from_lookup = !(normalize_ascii(Site) %in% site_lookup$lookup_key)
  )

missing_site_codes <- site_conversion_tbl |>
  filter(site_missing_from_lookup) |>
  pull(Site)

species_codes_in_arbre <- unique(arbre_inc$species_raw)
missing_species_codes <- species_codes_in_arbre[!(normalize_ascii(species_codes_in_arbre) %in% species_lookup$lookup_key)]

converted_site_count <- sum(!site_conversion_tbl$site_missing_from_lookup)
converted_species_count <- sum(normalize_ascii(species_codes_in_arbre) %in% species_lookup$lookup_key)

# -----------------------------
# Final table
# -----------------------------
final_tbl <- site_ba |>
  left_join(dbh_site, by = "Site") |>
  left_join(gaule_site |> select(Site, `Beech Sapling Basal Area`), by = "Site") |>
  left_join(dominant_species |> select(Site, `Top 3 Dominant Species`), by = "Site") |>
  left_join(latlon_elev |> select(Site, latitude_used_for_sorting, Latitude, Elevation), by = "Site") |>
  left_join(site_conversion_tbl |> select(Site, Site_label), by = "Site") |>
  arrange(latitude_used_for_sorting) |>
  transmute(
    Site = Site_label,
    Latitude,
    Elevation,
    `Basal Area` = Basal_Area,
    `Beech Basal Area` = Beech_Basal_Area,
    `Mean DBH` = Mean_DBH,
    `Mean Beech DBH` = Mean_Beech_DBH,
    `Beech Sapling Basal Area`,
    `Top 3 Dominant Species`
  )

# -----------------------------
# Diagnostics
# -----------------------------
diagnostics_tbl <- site_ba |>
  left_join(gaule_site |>
              select(Site, n_gaule_beech, gaule_radius_values_detected, gaule_subplot_area_ha, beech_sapling_basal_area_by_subplot),
            by = "Site") |>
  left_join(latlon_elev |> select(Site, n_elevation_records, latitude_used_for_sorting, Latitude, Elevation), by = "Site") |>
  left_join(dominant_species |>
              select(Site, top_3_species_original_codes, top_3_species_converted_codes, top_3_species_basal_area_values),
            by = "Site") |>
  left_join(site_conversion_tbl, by = "Site") |>
  mutate(
    species_codes_missing_from_species_lookup = collapse_or_none(missing_species_codes),
    site_lookup_used = site_lookup_used,
    species_lookup_used = species_lookup_used,
    prism_baf_used = prism_baf,
    R_utilise_values_detected_for_gaule = ifelse(length(r_values) == 0, "None", paste(round(r_values, 4), collapse = ", ")),
    basal_area_method_used = "Prism BA from included trees (Valeur_prisme==2): per-point n*prism_baf then averaged across points",
    gaule_method_used = "Fixed-radius subplot BA from R_utilise and DBH: subplot BA m²/ha then averaged across points"
  ) |>
  transmute(
    original_Site_code = Site,
    converted_Site_label = Site_label,
    Latitude,
    Elevation,
    top_3_species_original_codes,
    top_3_species_converted_codes,
    top_3_species_basal_area_values,
    species_codes_missing_from_species_lookup,
    site_lookup_used,
    species_lookup_used,
    prism_baf_used,
    R_utilise_values_detected_for_gaule,
    n_prism_points_detected,
    n_included_arbre_total,
    n_included_arbre_beech,
    n_gaule_beech,
    gaule_radius_values_detected,
    gaule_subplot_area_ha,
    beech_sapling_basal_area_by_subplot,
    n_elevation_records,
    latitude_used_for_sorting,
    basal_area_method_used,
    gaule_method_used
  ) |>
  arrange(latitude_used_for_sorting)

# -----------------------------
# Save outputs
# -----------------------------
csv_path <- file.path(out_dir, "site_description_table_article.csv")
xlsx_path <- file.path(out_dir, "site_description_table_article.xlsx")
diag_path <- file.path(out_dir, "site_description_table_article_diagnostics.csv")
docx_path <- file.path(out_dir, "site_description_table_article_formatted.docx")

write.csv(final_tbl, csv_path, row.names = FALSE, na = "")
openxlsx::write.xlsx(final_tbl, xlsx_path, overwrite = TRUE)
write.csv(diagnostics_tbl, diag_path, row.names = FALSE, na = "")

if (word_packages_available) {
  formatted_table <- flextable::flextable(final_tbl) |>
    flextable::theme_booktabs() |>
    flextable::autofit() |>
    flextable::align(align = "center", part = "all") |>
    flextable::align(j = "Top 3 Dominant Species", align = "left", part = "body") |>
    flextable::bold(part = "header")
  
  flextable::save_as_docx(
    `Site description table` = formatted_table,
    path = docx_path
  )
} else {
  message(
    "Skipping formatted Word table because optional package(s) are not installed: ",
    paste(missing_word_packages, collapse = ", "),
    ". Install flextable and officer to create: ",
    docx_path
  )
}

# -----------------------------
# Console summary
# -----------------------------
message("\nFinal table:")
print(final_tbl, n = Inf, width = Inf)

message("\nTop 3 dominant species per site:")
top_species_summary <- dominant_species |>
  left_join(site_conversion_tbl |> select(Site, Site_label), by = "Site") |>
  arrange(match(Site_label, final_tbl$Site)) |>
  transmute(Site = Site_label, top_3_species_converted_codes, top_3_species_basal_area_values)
print(top_species_summary, n = Inf, width = Inf)

message("\nSummary:")
message(" - site_lookup sheet detected: ", ifelse(site_lookup_used, site_lookup_sheet, "No"))
message(" - species_lookup sheet detected: ", ifelse(species_lookup_used, species_lookup_sheet, "No"))
message(" - number of site codes converted: ", converted_site_count)
message(" - missing site codes from site_lookup: ", collapse_or_none(missing_site_codes))
message(" - number of species codes converted: ", converted_species_count)
message(" - missing species codes from species_lookup: ", collapse_or_none(missing_species_codes))
message(" - top 3 dominant species per site: ", paste0(top_species_summary$Site, " = ", top_species_summary$top_3_species_converted_codes, collapse = " | "))
message(" - Latitude range: ", range_or_na(final_tbl$Latitude, digits = 5))
message(" - Elevation range: ", range_or_na(final_tbl$Elevation))
message(" - Basal Area range: ", range_or_na(final_tbl$`Basal Area`))
message(" - Beech Basal Area range: ", range_or_na(final_tbl$`Beech Basal Area`))
message(" - Beech Sapling Basal Area range: ", range_or_na(final_tbl$`Beech Sapling Basal Area`))
message(" - prism_baf used = ", prism_baf)
message(" - R_utilise values detected in gaule = ", ifelse(length(r_values) == 0, "None", paste(round(r_values, 4), collapse = ", ")))
message(" - Any missing/invalid R_utilise in beech gaule rows = No (script would stop otherwise)")
message(" - Sites with Beech Basal Area = 0: ", {s <- final_tbl$Site[!is.na(final_tbl$`Beech Basal Area`) & final_tbl$`Beech Basal Area` == 0]; if (length(s) == 0) "None" else paste(s, collapse = ", ")})
message(" - Sites with Mean Beech DBH = NA: ", {s <- final_tbl$Site[is.na(final_tbl$`Mean Beech DBH`)]; if (length(s) == 0) "None" else paste(s, collapse = ", ")})

message("\nSaved:")
message(" - ", csv_path)
message(" - ", xlsx_path)
message(" - ", diag_path)
message(
  " - ",
  docx_path,
  ifelse(
    word_packages_available,
    "",
    paste0(" (skipped: missing ", paste(missing_word_packages, collapse = ", "), ")")
  )
)