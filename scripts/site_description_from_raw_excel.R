#!/usr/bin/env Rscript

required_packages <- c("tidyverse", "readxl", "openxlsx", "janitor", "stringr", "here")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(paste0("Missing required package(s): ", paste(missing_packages, collapse = ", "), ". Please install them before running this script."), call. = FALSE)
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(openxlsx)
  library(janitor)
  library(stringr)
  library(here)
})

# -----------------------------
# User-editable constants
# -----------------------------
plot_area_m2 <- 24 * 24
beech_code <- "HEG"
manual_site_column <- NA_character_
min_expected_sites <- 8

# -----------------------------
# Helpers
# -----------------------------
resolve_project_root <- function() {
  candidate_roots <- c(here::here(), getwd())
  candidate_roots <- unique(normalizePath(candidate_roots, winslash = "/", mustWork = FALSE))
  fallback <- NULL
  
  for (cand in candidate_roots) {
    probe <- cand
    for (i in seq_len(6)) {
      has_data <- dir.exists(file.path(probe, "data"))
      has_scripts <- dir.exists(file.path(probe, "scripts"))
      has_raw <- dir.exists(file.path(probe, "data", "raw"))
      if (has_data && has_scripts) {
        if (has_raw) return(probe)
        fallback <- probe
      }
      parent <- dirname(probe)
      if (identical(parent, probe)) break
      probe <- parent
    }
  }
  
  if (!is.null(fallback)) return(fallback)
  normalizePath(here::here(), winslash = "/", mustWork = FALSE)
}

normalize_ascii <- function(x) {
  x %>%
    iconv(from = "", to = "ASCII//TRANSLIT") %>%
    tolower() %>%
    str_replace_all("[^a-z0-9]+", "_") %>%
    str_replace_all("_+", "_") %>%
    str_replace_all("^_|_$", "")
}

find_excel_file <- function(search_dir) {
  if (!dir.exists(search_dir)) {
    stop(paste0("Expected raw-data directory not found: ", search_dir), call. = FALSE)
  }
  excel_files <- list.files(search_dir, pattern = "\\.(xlsx|xls)$", full.names = TRUE, ignore.case = TRUE)
  if (length(excel_files) == 0) stop(paste0("No Excel files found in: ", search_dir), call. = FALSE)
  
  preferred <- excel_files[str_detect(normalize_ascii(basename(excel_files)), "donnee|donn?ee|west")]
  selected <- if (length(preferred) >= 1) preferred[1] else excel_files[1]
  
  message("Candidate Excel file(s):")
  message(paste0(" - ", basename(excel_files), collapse = "\n"))
  message("Selected workbook: ", basename(selected))
  selected
}

pick_column <- function(df, patterns, label) {
  cn <- names(df)
  cn_norm <- normalize_ascii(cn)
  idx <- which(str_detect(cn_norm, str_c(patterns, collapse = "|")))
  if (length(idx) == 0) return(NA_character_)
  if (length(idx) > 1) message("Multiple matches for ", label, ": ", paste(cn[idx], collapse = ", "), ". Using first: ", cn[idx[1]])
  cn[idx[1]]
}

coerce_numeric <- function(x) {
  if (is.numeric(x)) return(x)
  x_chr <- as.character(x)
  x_chr <- str_replace_all(x_chr, ",", ".")
  x_chr <- str_replace_all(x_chr, "[^0-9.\\-]+", "")
  suppressWarnings(as.numeric(x_chr))
}

non_empty_unique <- function(x) {
  x_chr <- str_trim(as.character(x))
  x_chr <- x_chr[!is.na(x_chr) & x_chr != "" & x_chr != "NA"]
  unique(x_chr)
}

infer_site_columns <- function(df) {
  cn <- names(df)
  cn_norm <- normalize_ascii(cn)
  site_like_patterns <- c("^site$", "site_?id", "id_?site", "code_?site", "num_?site", "no_?site", "placette", "parcelle", "plot", "station", "localite", "location", "bloc", "block")
  idx <- which(str_detect(cn_norm, str_c(site_like_patterns, collapse = "|")))
  
  distinct_counts <- vapply(df, function(col) length(non_empty_unique(col)), integer(1))
  broad_idx <- which(distinct_counts >= 4 & distinct_counts <= 100)
  unique(c(cn[idx], cn[broad_idx]))
}

choose_best_site_column <- function(df, candidate_cols) {
  if (length(candidate_cols) == 0) stop("No possible site column candidates were found in 'arbre'.", call. = FALSE)
  
  diagnostics <- purrr::map_dfr(candidate_cols, function(colname) {
    vals_raw <- non_empty_unique(df[[colname]])
    vals_num <- suppressWarnings(as.numeric(vals_raw))
    numeric_share <- ifelse(length(vals_raw) == 0, 0, mean(!is.na(vals_num)))
    
    tibble(
      column_name = colname,
      n_unique = length(vals_raw),
      numeric_share = numeric_share,
      unique_values_preview = paste(utils::head(vals_raw, 25), collapse = ", ")
    )
  }) %>%
    mutate(
      name_bonus = if_else(str_detect(normalize_ascii(column_name), "^site$|site_?id|id_?site|placette|parcelle|plot|station"), 1000, 0),
      numeric_penalty = if_else(numeric_share > 0.95, 200, 0),
      score = name_bonus + n_unique - numeric_penalty
    ) %>%
    arrange(desc(score), desc(n_unique), column_name)
  
  message("\nPossible site/location columns in 'arbre' (diagnostics):")
  print(diagnostics, n = Inf, width = Inf)
  
  list(best_column = diagnostics$column_name[1], diagnostics = diagnostics)
}

# -----------------------------
# Load workbook
# -----------------------------
project_root <- resolve_project_root()
message("Project root resolved to: ", project_root)

raw_dir <- file.path(project_root, "data", "raw")
excel_path <- find_excel_file(raw_dir)
message("Using Excel file: ", excel_path)

all_sheets <- readxl::excel_sheets(excel_path)
if (!("arbre" %in% all_sheets)) {
  stop(paste0("Sheet 'arbre' not found. Available sheets: ", paste(all_sheets, collapse = ", ")), call. = FALSE)
}

arbre_raw <- readxl::read_excel(excel_path, sheet = "arbre") %>% janitor::clean_names()
message("All columns in 'arbre':")
message(paste0(" - ", names(arbre_raw), collapse = "\n"))

# -----------------------------
# Detect columns
# -----------------------------
species_col <- pick_column(arbre_raw, c("espece", "species", "specie", "taxon", "code_?ess", "essence", "sp$", "^sp_"), "species")
dbh_col <- pick_column(arbre_raw, c("dbh", "dhp", "diam", "diametre", "diameter", "d130", "d_?bh"), "DBH")
ba_col <- pick_column(arbre_raw, c("basal", "surface_?terr", "g_?ind", "ba_?m2", "barea", "st_tige"), "basal area")
plot_area_col <- pick_column(arbre_raw, c("plot_?area", "surface_?plac", "area_?m2", "superficie", "quadrat_?area"), "plot area")

site_candidates <- infer_site_columns(arbre_raw)
site_choice <- choose_best_site_column(arbre_raw, site_candidates)
site_col <- site_choice$best_column

if (!is.na(manual_site_column)) {
  if (!(manual_site_column %in% names(arbre_raw))) {
    stop(paste0("manual_site_column '", manual_site_column, "' does not exist. Available columns: ", paste(names(arbre_raw), collapse = ", ")), call. = FALSE)
  }
  site_col <- manual_site_column
  message("Manual site column override applied: ", site_col)
}

message("\nDetected columns in 'arbre':")
message(" - Site (selected): ", site_col)
message(" - Species: ", species_col)
message(" - DBH: ", dbh_col)
message(" - Basal area (existing): ", ba_col)
message(" - Plot area in arbre: ", plot_area_col)

required_for_core <- c(site_col, species_col, dbh_col)
if (any(is.na(required_for_core))) {
  stop(paste0("Failed to detect required columns. Detected: site=", site_col, ", species=", species_col, ", dbh=", dbh_col), call. = FALSE)
}

# required diagnostics
message("\nSite-column diagnostics for all possible site/location columns:")
for (colname in site_candidates) {
  vals <- non_empty_unique(arbre_raw[[colname]])
  message("\nColumn: ", colname)
  message(" - Number of unique values: ", length(vals))
  message(" - Unique values: ", paste(vals, collapse = ", "))
}

arbre <- arbre_raw %>%
  mutate(
    site_id_raw = str_trim(as.character(.data[[site_col]])),
    species_clean = str_to_upper(str_trim(as.character(.data[[species_col]]))),
    dbh_cm = coerce_numeric(.data[[dbh_col]])
  ) %>%
  mutate(site_id_raw = if_else(site_id_raw == "", NA_character_, site_id_raw))

unique_sites <- length(unique(na.omit(arbre$site_id_raw)))
message("\nDiagnostic - unique sites detected in selected column '", site_col, "': ", unique_sites)

if (unique_sites < min_expected_sites && is.na(manual_site_column)) {
  stop(paste0("Automatic site-column detection found only ", unique_sites, " sites (< ", min_expected_sites, "). Set manual_site_column. Candidates: ", paste(site_candidates, collapse = ", ")), call. = FALSE)
}

site_frequency <- arbre %>% filter(!is.na(site_id_raw)) %>% count(site_id_raw, name = "n_individuals") %>% arrange(site_id_raw)
message("Diagnostic - individuals per site:")
print(site_frequency, n = Inf)

diag_preview <- arbre %>% transmute(site = site_id_raw, espece = .data[[species_col]], dhp_cm = dbh_cm)
message("Diagnostic - first 20 rows of site/espece/dhp_cm:")
print(utils::head(diag_preview, 20), n = 20)

if (all(is.na(arbre$site_id_raw))) stop("All site IDs are missing after processing selected site column.", call. = FALSE)

plot_area_m2_used <- plot_area_m2
if (!is.na(plot_area_col)) {
  pa <- coerce_numeric(arbre[[plot_area_col]])
  pa <- pa[!is.na(pa) & pa > 0]
  if (length(pa) > 0) {
    plot_area_m2_used <- unique(pa)[1]
    message("Using detected plot area from data: ", plot_area_m2_used, " m2")
  } else {
    message("Plot area column detected but invalid; using fallback: ", plot_area_m2)
  }
} else {
  message("No plot area column detected; using fallback plot_area_m2 = ", plot_area_m2)
}

if (!is.na(ba_col)) {
  arbre <- arbre %>% mutate(individual_basal_area_m2 = coerce_numeric(.data[[ba_col]]))
  idx <- is.na(arbre$individual_basal_area_m2)
  arbre$individual_basal_area_m2[idx] <- pi * (arbre$dbh_cm[idx] / 200)^2
} else {
  arbre <- arbre %>% mutate(individual_basal_area_m2 = pi * (dbh_cm / 200)^2)
}

all_sites <- arbre %>% filter(!is.na(site_id_raw), site_id_raw != "") %>% distinct(site_id_raw)

site_summary <- arbre %>%
  filter(!is.na(site_id_raw), site_id_raw != "") %>%
  group_by(site_id_raw) %>%
  summarise(
    basal_area_m2_ha = ifelse(all(is.na(individual_basal_area_m2)), NA_real_, sum(individual_basal_area_m2, na.rm = TRUE) / plot_area_m2_used * 10000),
    mean_dbh_cm = ifelse(all(is.na(dbh_cm)), NA_real_, mean(dbh_cm, na.rm = TRUE)),
    sd_dbh_cm = ifelse(all(is.na(dbh_cm)), NA_real_, sd(dbh_cm, na.rm = TRUE)),
    .groups = "drop"
  )

beech_summary <- arbre %>%
  filter(!is.na(site_id_raw), site_id_raw != "", species_clean == beech_code) %>%
  group_by(site_id_raw) %>%
  summarise(
    beech_basal_area_m2_ha = ifelse(all(is.na(individual_basal_area_m2)), NA_real_, sum(individual_basal_area_m2, na.rm = TRUE) / plot_area_m2_used * 10000),
    mean_beech_dbh_cm = ifelse(all(is.na(dbh_cm)), NA_real_, mean(dbh_cm, na.rm = TRUE)),
    sd_beech_dbh_cm = ifelse(all(is.na(dbh_cm)), NA_real_, sd(dbh_cm, na.rm = TRUE)),
    .groups = "drop"
  )

site_description <- all_sites %>%
  left_join(site_summary, by = "site_id_raw") %>%
  left_join(beech_summary, by = "site_id_raw") %>%
  arrange(site_id_raw) %>%
  transmute(
    site_label = site_id_raw,
    basal_area_m2_ha,
    mean_dbh_cm,
    sd_dbh_cm,
    beech_basal_area_m2_ha,
    mean_beech_dbh_cm,
    sd_beech_dbh_cm
  )

formatted_table <- site_description %>%
  mutate(
    `Mean DBH (cm)` = if_else(is.na(mean_dbh_cm), NA_character_, sprintf("%.1f (%.2f)", mean_dbh_cm, if_else(is.na(sd_dbh_cm), 0, sd_dbh_cm))),
    `Mean beech DBH (cm)` = if_else(is.na(mean_beech_dbh_cm), NA_character_, sprintf("%.1f (%.2f)", mean_beech_dbh_cm, if_else(is.na(sd_beech_dbh_cm), 0, sd_beech_dbh_cm)))
  ) %>%
  transmute(
    `Site #` = site_label,
    `Basal area (m²·ha⁻¹)` = round(basal_area_m2_ha, 3),
    `Mean DBH (cm)`,
    `Beech basal area (m²·ha⁻¹)` = round(beech_basal_area_m2_ha, 3),
    `Mean beech DBH (cm)`
  )

out_dir <- file.path(project_root, "outputs", "tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

readr::write_csv(site_description, file.path(out_dir, "site_description_table.csv"), na = "")
openxlsx::write.xlsx(site_description, file.path(out_dir, "site_description_table.xlsx"), overwrite = TRUE)
readr::write_csv(formatted_table, file.path(out_dir, "site_description_table_formatted.csv"), na = "")

message("\nFinal formatted table:")
print(formatted_table, n = Inf)

message("\nOutputs saved to:")
message(" - ", file.path(out_dir, "site_description_table.csv"))
message(" - ", file.path(out_dir, "site_description_table.xlsx"))
message(" - ", file.path(out_dir, "site_description_table_formatted.csv"))