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
plot_area_m2 <- 24 * 24  # Used if no plot area can be detected automatically.
beech_code <- "HEG"
manual_site_column <- NA_character_  # Set to exact arbre column name to override auto-detection.

# -----------------------------
# Helpers
# -----------------------------
resolve_project_root <- function() {
  candidate_roots <- c(here::here(), getwd())
  candidate_roots <- unique(normalizePath(candidate_roots, winslash = "/", mustWork = FALSE))
  
  for (cand in candidate_roots) {
    probe <- cand
    for (i in seq_len(5)) {
      has_data <- dir.exists(file.path(probe, "data"))
      has_scripts <- dir.exists(file.path(probe, "scripts"))
      has_raw <- dir.exists(file.path(probe, "data", "raw"))
      if (has_data && has_scripts) {
        if (has_raw) {
          return(probe)
        }
        # still keep as fallback if data/raw is created later
        fallback <- probe
      }
      parent <- dirname(probe)
      if (identical(parent, probe)) break
      probe <- parent
    }
  }
  
  if (exists("fallback", inherits = FALSE)) return(fallback)
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
  
  if (length(preferred) >= 1) {
    message("Candidate Excel file(s):")
    message(paste0(" - ", basename(preferred), collapse = "\n"))
    return(preferred[1])
  }
  
  message("No filename match for 'donnee/donnée/west'; using first Excel file found.")
  message("Selected: ", basename(excel_files[1]))
  excel_files[1]
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
  # keep digits, decimal separators, sign; replace comma by dot
  x_chr <- str_replace_all(x_chr, ",", ".")
  x_chr <- str_replace_all(x_chr, "[^0-9.\\-]+", "")
  suppressWarnings(as.numeric(x_chr))
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

# -----------------------------
# Detect relevant columns
# -----------------------------
site_col <- pick_column(arbre_raw, c("^site$", "site_?id", "id_?site", "placette", "plot", "parcelle"), "site")
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
species_col <- pick_column(arbre_raw, c("espece", "species", "specie", "taxon", "code_?ess", "essence", "sp$", "^sp_"), "species")
dbh_col <- pick_column(arbre_raw, c("dbh", "dhp", "diam", "diametre", "diameter", "d130", "d_?bh"), "DBH")
ba_col <- pick_column(arbre_raw, c("basal", "surface_?terr", "g_?ind", "ba_?m2", "barea"), "basal area")
elev_col_arbre <- pick_column(arbre_raw, c("elev", "alt", "altitude", "z_?m", "height_?asl"), "elevation")
plot_area_col <- pick_column(arbre_raw, c("plot_?area", "surface_?plac", "area_?m2", "superficie", "quadrat_?area"), "plot area")

message("Detected columns in 'arbre':")
message(" - Site: ", site_col)
message(" - Species: ", species_col)
message(" - DBH: ", dbh_col)
message(" - Basal area (existing): ", ba_col)
message(" - Elevation in arbre: ", elev_col_arbre)
message(" - Plot area in arbre: ", plot_area_col)
message("All columns in 'arbre': ", paste(names(arbre_raw), collapse = ", "))

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
    site_clean = coerce_numeric(.data[[site_col]]),
    species_clean = str_to_upper(str_trim(as.character(.data[[species_col]]))),
    dbh_cm = coerce_numeric(.data[[dbh_col]])
  )

# -----------------------------
# Diagnostics
# -----------------------------
unique_sites <- arbre %>%
  filter(!is.na(site_clean)) %>%
  distinct(site_clean) %>%
  nrow()
message("Diagnostic - unique sites detected: ", unique_sites)

site_frequency <- arbre %>%
  filter(!is.na(site_clean)) %>%
  count(site_clean, name = "n_individuals") %>%
  arrange(site_clean)
message("Diagnostic - individuals per site:")
print(site_frequency, n = Inf)

diag_preview <- arbre %>%
  transmute(
    site = .data[[site_col]],
    espece = .data[[species_col]],
    dhp_cm = dbh_cm
  )
message("Diagnostic - first 20 rows of site/espece/dhp_cm:")
print(utils::head(diag_preview, 20), n = 20)

if (all(is.na(arbre$site_clean))) {
  stop("All site values are NA after numeric conversion; cannot continue.", call. = FALSE)
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

# Elevation extraction
get_elevation_from_other_sheets <- function(path, sheets, site_values) {
  checked <- list()
  for (sh in sheets) {
    if (tolower(sh) == "arbre") next
    
    dat <- tryCatch(readxl::read_excel(path, sheet = sh) %>% janitor::clean_names(), error = function(e) NULL)
    if (is.null(dat) || ncol(dat) == 0) next
    
    sc <- pick_column(dat, c("^site$", "site_?id", "id_?site", "placette", "plot", "parcelle"), "site")
    ec <- pick_column(dat, c("elev", "alt", "altitude", "z_?m", "height_?asl"), "elevation")
    
    checked[[length(checked) + 1]] <- list(sheet = sh, site_col = sc, elev_col = ec, columns = names(dat))
    
    if (!is.na(sc) && !is.na(ec)) {
      elev_df <- dat %>%
        transmute(
          site_clean = coerce_numeric(.data[[sc]]),
          elevation_m = coerce_numeric(.data[[ec]])
        ) %>%
        filter(!is.na(site_clean), !is.na(elevation_m)) %>%
        group_by(site_clean) %>%
        summarise(elevation_m = mean(elevation_m, na.rm = TRUE), .groups = "drop")
      
      if (nrow(elev_df) > 0 && any(elev_df$site_clean %in% site_values)) {
        message("Using elevation from sheet '", sh, "' (site column: ", sc, ", elevation column: ", ec, ")")
        return(list(elevation = elev_df, checked = checked))
      }
    }
  }
  list(elevation = NULL, checked = checked)
}

if (!is.na(elev_col_arbre)) {
  elevation_by_site <- arbre %>%
    transmute(site_clean, elevation_m = coerce_numeric(.data[[elev_col_arbre]])) %>%
    filter(!is.na(site_clean), !is.na(elevation_m)) %>%
    group_by(site_clean) %>%
    summarise(elevation_m = mean(elevation_m, na.rm = TRUE), .groups = "drop")
  
  if (nrow(elevation_by_site) == 0) {
    elev_lookup <- get_elevation_from_other_sheets(excel_path, all_sheets, unique(arbre$site_clean))
    elevation_by_site <- elev_lookup$elevation
  }
} else {
  elev_lookup <- get_elevation_from_other_sheets(excel_path, all_sheets, unique(arbre$site_clean))
  elevation_by_site <- elev_lookup$elevation
}

if (is.null(elevation_by_site) || nrow(elevation_by_site) == 0) {
  checked_msg <- "Checked elevation candidates in 'arbre' and other sheets."
  if (exists("elev_lookup") && length(elev_lookup$checked) > 0) {
    checked_detail <- vapply(
      elev_lookup$checked,
      function(x) paste0("Sheet '", x$sheet, "' | site_col=", x$site_col, " | elev_col=", x$elev_col,
                         " | columns=[", paste(x$columns, collapse = ", "), "]"),
      character(1)
    )
    checked_msg <- paste(c(checked_msg, checked_detail), collapse = "\n")
  } else {
    checked_msg <- paste0(
      checked_msg,
      "\nAvailable sheets: ", paste(all_sheets, collapse = ", "),
      "\nArbre columns: ", paste(names(arbre_raw), collapse = ", ")
    )
  }
  stop(
    paste0("No elevation data found. ", checked_msg),
    call. = FALSE
  )
}

# -----------------------------
# Build outputs
# -----------------------------
all_sites <- arbre %>%
  filter(!is.na(site_clean)) %>%
  distinct(site_clean)

site_summary <- arbre %>%
  filter(!is.na(site_clean)) %>%
  group_by(site_clean) %>%
  summarise(
    basal_area_m2_ha = ifelse(
      all(is.na(individual_basal_area_m2)),
      NA_real_,
      sum(individual_basal_area_m2, na.rm = TRUE) / plot_area_m2_used * 10000
    ),
    mean_dbh_cm = ifelse(all(is.na(dbh_cm)), NA_real_, mean(dbh_cm, na.rm = TRUE)),
    sd_dbh_cm = ifelse(all(is.na(dbh_cm)), NA_real_, sd(dbh_cm, na.rm = TRUE)),
    .groups = "drop"
  )

beech_summary <- arbre %>%
  filter(!is.na(site_clean), species_clean == beech_code, !is.na(dbh_cm), !is.na(individual_basal_area_m2)) %>%
  group_by(site_clean) %>%
  summarise(
    beech_basal_area_m2_ha = sum(individual_basal_area_m2, na.rm = TRUE) / plot_area_m2_used * 10000,
    mean_beech_dbh_cm = mean(dbh_cm, na.rm = TRUE),
    sd_beech_dbh_cm = sd(dbh_cm, na.rm = TRUE),
    .groups = "drop"
  )

site_description <- all_sites %>%
  left_join(site_summary, by = "site_clean") %>%
  left_join(beech_summary, by = "site_clean") %>%
  left_join(elevation_by_site, by = "site_clean") %>%
  arrange(site_clean) %>%
  mutate(site_number = as.integer(round(site_clean))) %>%
  select(
    site_number,
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
    `Site #` = site_number,
    `Basal area (m²·ha⁻¹)` = round(basal_area_m2_ha, 3),
    `Mean DBH (cm)`,
    `Beech basal area (m²·ha⁻¹)` = round(beech_basal_area_m2_ha, 3),
    `Mean beech DBH (cm)`,
    `Elevation (m)` = round(elevation_m, 1)
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