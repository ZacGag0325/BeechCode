#!/usr/bin/env Rscript

# Standalone script: build a site-description table from raw Excel data.
# This script intentionally does NOT modify or hook into any master pipeline.

required_packages <- c("tidyverse", "readxl", "openxlsx", "janitor", "stringr", "here")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(
    paste0(
      "Missing required package(s): ", paste(missing_packages, collapse = ", "),
      ". Please install them before running this script."
    ),
    call. = FALSE
  )
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
plot_area_m2 <- 24 * 24  # Fallback if no valid plot area can be detected.
beech_code <- "HEG"
manual_site_column <- NA_character_  # Set to exact arbre column name to force site grouping.
min_expected_sites <- 8              # Auto-detection must find at least this many sites.

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
    stop(
      paste0("Expected raw-data directory not found: ", search_dir,
             ". Please create it or update the path in this script."),
      call. = FALSE
    )
  }
  
  excel_files <- list.files(search_dir, pattern = "\\.(xlsx|xls)$", full.names = TRUE, ignore.case = TRUE)
  if (length(excel_files) == 0) {
    stop(paste0("No Excel files found in: ", search_dir), call. = FALSE)
  }
  
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
  if (length(idx) > 1) {
    message("Multiple matches for ", label, ": ", paste(cn[idx], collapse = ", "), ". Using first: ", cn[idx[1]])
  }
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
  
  site_like_patterns <- c(
    "^site$", "site_?id", "id_?site", "code_?site", "num_?site", "no_?site",
    "placette", "parcelle", "plot", "station", "localite", "location", "bloc", "block"
  )
  
  idx <- which(str_detect(cn_norm, str_c(site_like_patterns, collapse = "|")))
  # also include any column with moderate distinctness that might be site ID
  distinct_counts <- vapply(df, function(col) length(non_empty_unique(col)), integer(1))
  broad_idx <- which(distinct_counts >= 4 & distinct_counts <= 100)
  
  unique(c(cn[idx], cn[broad_idx]))
}

choose_best_site_column <- function(df, candidate_cols, min_expected = 8) {
  if (length(candidate_cols) == 0) {
    stop("No possible site column candidates were found in 'arbre'.", call. = FALSE)
  }
  
  diagnostics <- purrr::map_dfr(candidate_cols, function(colname) {
    vals_raw <- non_empty_unique(df[[colname]])
    vals_num <- suppressWarnings(as.numeric(vals_raw))
    vals_num_non_na <- vals_num[!is.na(vals_num)]
    
    tibble(
      column_name = colname,
      n_unique_raw = length(vals_raw),
      n_unique_numeric = length(unique(vals_num_non_na)),
      numeric_share = ifelse(length(vals_raw) == 0, 0, mean(!is.na(vals_num))),
      unique_values_preview = paste(utils::head(vals_raw, 25), collapse = ", ")
    )
  }) %>%
    mutate(
      n_unique_effective = pmax(n_unique_raw, n_unique_numeric),
      score = n_unique_effective + ifelse(str_detect(normalize_ascii(column_name), "site|placette|parcelle|plot"), 0.5, 0)
    ) %>%
    arrange(desc(score), desc(n_unique_effective), column_name)
  
  message("\nPossible site/location columns in 'arbre' (diagnostics):")
  print(diagnostics, n = Inf, width = Inf)
  
  best <- diagnostics$column_name[1]
  best_n <- diagnostics$n_unique_effective[1]
  
  list(best_column = best, best_n = best_n, diagnostics = diagnostics)
}

# -----------------------------
# Locate project and input
# -----------------------------
project_root <- resolve_project_root()
message("Project root resolved to: ", project_root)

raw_dir <- file.path(project_root, "data", "raw")
excel_path <- find_excel_file(raw_dir)
message("Using Excel file: ", excel_path)

all_sheets <- readxl::excel_sheets(excel_path)
if (!("arbre" %in% all_sheets)) {
  stop(
    paste0("Sheet 'arbre' not found in workbook. Available sheets: ", paste(all_sheets, collapse = ", ")),
    call. = FALSE
  )
}

arbre_raw <- readxl::read_excel(excel_path, sheet = "arbre") %>% janitor::clean_names()

message("\nAll columns in 'arbre':")
message(paste0(" - ", names(arbre_raw), collapse = "\n"))

# -----------------------------
# Detect core columns
# -----------------------------
species_col <- pick_column(arbre_raw, c("espece", "species", "specie", "taxon", "code_?ess", "essence", "sp$", "^sp_"), "species")
dbh_col <- pick_column(arbre_raw, c("dbh", "dhp", "diam", "diametre", "diameter", "d130", "d_?bh"), "DBH")
ba_col <- pick_column(arbre_raw, c("basal", "surface_?terr", "g_?ind", "ba_?m2", "barea"), "basal area")
elev_col_arbre <- pick_column(arbre_raw, c("elev", "alt", "altitude", "z_?m", "height_?asl"), "elevation")
plot_area_col <- pick_column(arbre_raw, c("plot_?area", "surface_?plac", "area_?m2", "superficie", "quadrat_?area"), "plot area")

site_candidates <- infer_site_columns(arbre_raw)
site_choice <- choose_best_site_column(arbre_raw, site_candidates, min_expected_sites)

site_col <- site_choice$best_column
if (!is.na(manual_site_column)) {
  if (!(manual_site_column %in% names(arbre_raw))) {
    stop(
      paste0(
        "manual_site_column was provided ('", manual_site_column,
        "') but does not exist in arbre sheet. Available columns: ",
        paste(names(arbre_raw), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  site_col <- manual_site_column
  message("Manual site column override applied: ", site_col)
}

message("\nDetected columns in 'arbre':")
message(" - Site (selected): ", site_col)
message(" - Species: ", species_col)
message(" - DBH: ", dbh_col)
message(" - Basal area (existing): ", ba_col)
message(" - Elevation in arbre: ", elev_col_arbre)
message(" - Plot area in arbre: ", plot_area_col)

required_for_core <- c(site_col, species_col, dbh_col)
if (any(is.na(required_for_core))) {
  stop(
    paste0(
      "Failed to detect required columns in 'arbre'.\n",
      "Detected: site=", site_col, ", species=", species_col, ", dbh=", dbh_col, ".\n",
      "Available columns in arbre: ", paste(names(arbre_raw), collapse = ", ")
    ),
    call. = FALSE
  )
}

arbre <- arbre_raw %>%
  mutate(
    site_id_raw = str_trim(as.character(.data[[site_col]])),
    species_clean = str_to_upper(str_trim(as.character(.data[[species_col]]))),
    dbh_cm = coerce_numeric(.data[[dbh_col]])
  ) %>%
  mutate(site_id_raw = if_else(site_id_raw == "", NA_character_, site_id_raw))

# -----------------------------
# Site diagnostics (required)
# -----------------------------
message("\nSite-column diagnostics for all possible site/location columns:")
for (colname in site_candidates) {
  vals <- non_empty_unique(arbre_raw[[colname]])
  message("\nColumn: ", colname)
  message(" - Number of unique values: ", length(vals))
  message(" - Unique values: ", paste(vals, collapse = ", "))
}

unique_sites <- length(unique(na.omit(arbre$site_id_raw)))
message("\nDiagnostic - unique sites detected in selected column '", site_col, "': ", unique_sites)

site_frequency <- arbre %>%
  filter(!is.na(site_id_raw)) %>%
  count(site_id_raw, name = "n_individuals") %>%
  arrange(site_id_raw)
message("Diagnostic - individuals per site:")
print(site_frequency, n = Inf)

diag_preview <- arbre %>%
  transmute(
    site = site_id_raw,
    espece = .data[[species_col]],
    dhp_cm = dbh_cm
  )
message("Diagnostic - first 20 rows of site/espece/dhp_cm:")
print(utils::head(diag_preview, 20), n = 20)

if (unique_sites < min_expected_sites && is.na(manual_site_column)) {
  stop(
    paste0(
      "Automatic site-column detection found only ", unique_sites, " unique sites (< ", min_expected_sites, ").\n",
      "Please set manual_site_column to the correct column name.\n",
      "Possible site columns found: ", paste(site_candidates, collapse = ", ")
    ),
    call. = FALSE
  )
}

if (all(is.na(arbre$site_id_raw))) {
  stop("All site IDs are missing after processing selected site column; cannot continue.", call. = FALSE)
}

# Determine plot area
plot_area_m2_used <- plot_area_m2
if (!is.na(plot_area_col)) {
  candidate_plot_area <- arbre %>% pull(.data[[plot_area_col]]) %>% coerce_numeric()
  candidate_plot_area <- candidate_plot_area[!is.na(candidate_plot_area) & candidate_plot_area > 0]
  if (length(candidate_plot_area) > 0) {
    plot_area_m2_used <- unique(candidate_plot_area)[1]
    message("Using detected plot area from data: ", plot_area_m2_used, " m2")
  } else {
    message("Plot area column detected but no valid positive values; using fallback: ", plot_area_m2, " m2")
  }
} else {
  message("No plot area column detected; using fallback plot_area_m2 = ", plot_area_m2)
}

# Individual basal area
if (!is.na(ba_col)) {
  arbre <- arbre %>% mutate(individual_basal_area_m2 = coerce_numeric(.data[[ba_col]]))
  missing_ba <- is.na(arbre$individual_basal_area_m2)
  if (any(missing_ba)) {
    arbre$individual_basal_area_m2[missing_ba] <- pi * (arbre$dbh_cm[missing_ba] / 200)^2
  }
} else {
  arbre <- arbre %>% mutate(individual_basal_area_m2 = pi * (dbh_cm / 200)^2)
}

# Elevation extraction and matching using corrected site ID key
get_elevation_from_other_sheets <- function(path, sheets, site_values) {
  checked <- list()
  target_norm <- normalize_ascii(site_values)
  
  for (sh in sheets) {
    if (tolower(sh) == "arbre") next
    
    dat <- tryCatch(readxl::read_excel(path, sheet = sh) %>% janitor::clean_names(), error = function(e) NULL)
    if (is.null(dat) || ncol(dat) == 0) next
    
    elev_col <- pick_column(dat, c("elev", "alt", "altitude", "z_?m", "height_?asl"), "elevation")
    if (is.na(elev_col)) {
      checked[[length(checked) + 1]] <- list(sheet = sh, site_col = NA_character_, elev_col = elev_col, columns = names(dat))
      next
    }
    
    site_cols <- infer_site_columns(dat)
    if (length(site_cols) == 0) {
      checked[[length(checked) + 1]] <- list(sheet = sh, site_col = NA_character_, elev_col = elev_col, columns = names(dat))
      next
    }
    
    best_match <- NULL
    best_overlap <- -1
    for (sc in site_cols) {
      sc_vals <- non_empty_unique(dat[[sc]])
      overlap <- length(intersect(normalize_ascii(sc_vals), target_norm))
      if (overlap > best_overlap) {
        best_overlap <- overlap
        best_match <- sc
      }
    }
    
    checked[[length(checked) + 1]] <- list(sheet = sh, site_col = best_match, elev_col = elev_col, columns = names(dat), overlap = best_overlap)
    
    if (!is.null(best_match) && best_overlap > 0) {
      elev_df <- dat %>%
        transmute(
          site_id_raw = str_trim(as.character(.data[[best_match]])),
          elevation_m = coerce_numeric(.data[[elev_col]])
        ) %>%
        filter(!is.na(site_id_raw), site_id_raw != "", !is.na(elevation_m)) %>%
        group_by(site_id_raw) %>%
        summarise(elevation_m = mean(elevation_m, na.rm = TRUE), .groups = "drop")
      
      if (nrow(elev_df) > 0) {
        message("Using elevation from sheet '", sh, "' (site column: ", best_match, ", elevation column: ", elev_col, ", overlap with arbre sites: ", best_overlap, ")")
        return(list(elevation = elev_df, checked = checked))
      }
    }
  }
  
  list(elevation = NULL, checked = checked)
}

if (!is.na(elev_col_arbre)) {
  elevation_by_site <- arbre %>%
    transmute(site_id_raw, elevation_m = coerce_numeric(.data[[elev_col_arbre]])) %>%
    filter(!is.na(site_id_raw), site_id_raw != "", !is.na(elevation_m)) %>%
    group_by(site_id_raw) %>%
    summarise(elevation_m = mean(elevation_m, na.rm = TRUE), .groups = "drop")
  
  if (nrow(elevation_by_site) == 0) {
    elev_lookup <- get_elevation_from_other_sheets(excel_path, all_sheets, unique(na.omit(arbre$site_id_raw)))
    elevation_by_site <- elev_lookup$elevation
  }
} else {
  elev_lookup <- get_elevation_from_other_sheets(excel_path, all_sheets, unique(na.omit(arbre$site_id_raw)))
  elevation_by_site <- elev_lookup$elevation
}

if (is.null(elevation_by_site) || nrow(elevation_by_site) == 0) {
  checked_msg <- "Checked elevation candidates in 'arbre' and other sheets."
  if (exists("elev_lookup") && length(elev_lookup$checked) > 0) {
    checked_detail <- vapply(
      elev_lookup$checked,
      function(x) paste0("Sheet '", x$sheet, "' | site_col=", x$site_col, " | elev_col=", x$elev_col,
                         ifelse(!is.null(x$overlap), paste0(" | overlap=", x$overlap), ""),
                         " | columns=[", paste(x$columns, collapse = ", "), "]"),
      character(1)
    )
    checked_msg <- paste(c(checked_msg, checked_detail), collapse = "\n")
  }
  stop(paste0("No elevation data found. ", checked_msg), call. = FALSE)
}

# -----------------------------
# Build outputs (one row per real site)
# -----------------------------
all_sites <- arbre %>%
  filter(!is.na(site_id_raw), site_id_raw != "") %>%
  distinct(site_id_raw)

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
  left_join(elevation_by_site, by = "site_id_raw") %>%
  arrange(site_id_raw) %>%
  mutate(
    site_number = suppressWarnings(as.integer(as.numeric(site_id_raw))),
    site_number = if_else(is.na(site_number), NA_integer_, site_number),
    site_label = if_else(is.na(site_number), site_id_raw, as.character(site_number))
  ) %>%
  select(
    site_label,
    basal_area_m2_ha,
    mean_dbh_cm,
    sd_dbh_cm,
    beech_basal_area_m2_ha,
    mean_beech_dbh_cm,
    sd_beech_dbh_cm,
    elevation_m
  )

formatted_table <- site_description %>%
  mutate(
    `Mean DBH (cm)` = if_else(
      is.na(mean_dbh_cm),
      NA_character_,
      sprintf("%.1f (%.2f)", mean_dbh_cm, if_else(is.na(sd_dbh_cm), 0, sd_dbh_cm))
    ),
    `Mean beech DBH (cm)` = if_else(
      is.na(mean_beech_dbh_cm),
      NA_character_,
      sprintf("%.1f (%.2f)", mean_beech_dbh_cm, if_else(is.na(sd_beech_dbh_cm), 0, sd_beech_dbh_cm))
    )
  ) %>%
  transmute(
    `Site #` = site_label,
    `Basal area (m²·ha⁻¹)` = round(basal_area_m2_ha, 3),
    `Mean DBH (cm)`,
    `Beech basal area (m²·ha⁻¹)` = round(beech_basal_area_m2_ha, 3),
    `Mean beech DBH (cm)`,
    `Elevation (m)` = round(elevation_m, 1)
  )

out_dir <- file.path(project_root, "outputs", "tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Raw output with explicit stats
site_description_export <- site_description %>%
  rename(
    site_number = site_label
  )

readr::write_csv(site_description_export, file.path(out_dir, "site_description_table.csv"), na = "")
openxlsx::write.xlsx(site_description_export, file.path(out_dir, "site_description_table.xlsx"), overwrite = TRUE)
readr::write_csv(formatted_table, file.path(out_dir, "site_description_table_formatted.csv"), na = "")

message("\nFinal formatted table:")
print(formatted_table, n = Inf)

message("\nOutputs saved to:")
message(" - ", file.path(out_dir, "site_description_table.csv"))
message(" - ", file.path(out_dir, "site_description_table.xlsx"))
message(" - ", file.path(out_dir, "site_description_table_formatted.csv"))