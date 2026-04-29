#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(stringr)
  library(stringi)
  library(purrr)
  library(tidyr)
  library(janitor)
  library(readr)
  library(openxlsx)
})

# Standalone script to build site-level basal-area summary tables from raw Excel data.
# The script automatically detects arbre/arbre feuille and gaule/gaule feuille sheets,
# standardizes likely column names, computes basal area per stem, summarizes by site,
# and exports CSV + XLSX outputs.

detect_project_root <- function() {
  fallback_root <- path.expand("~/Desktop/BeechCode")
  
  normalize_root_candidate <- function(path_value) {
    if (is.null(path_value) || length(path_value) == 0 || is.na(path_value)) return(NA_character_)
    p <- normalizePath(as.character(path_value), winslash = "/", mustWork = FALSE)
    
    # If a .Rproj file path is returned, use its parent folder as project root.
    if (grepl("\\.rproj$", tolower(p))) {
      p <- dirname(p)
    }
    
    if (!dir.exists(p)) return(NA_character_)
    p
  }
  
  # 1) Prefer here::here() if available.
  if (requireNamespace("here", quietly = TRUE)) {
    candidate <- tryCatch(here::here(), error = function(e) NA_character_)
    candidate <- normalize_root_candidate(candidate)
    if (!is.na(candidate)) return(candidate)
  }
  
  # 2) Try rprojroot::find_rstudio_root_file() if available.
  if (requireNamespace("rprojroot", quietly = TRUE)) {
    candidate <- tryCatch(rprojroot::find_rstudio_root_file(), error = function(e) NA_character_)
    candidate <- normalize_root_candidate(candidate)
    if (!is.na(candidate)) return(candidate)
  }
  
  # 3) Use the script path if available (works when sourced from file).
  script_path <- tryCatch(normalizePath(sys.frames()[[1]]$ofile, winslash = "/", mustWork = FALSE),
                          error = function(e) NA_character_)
  script_path <- normalize_root_candidate(script_path)
  if (!is.na(script_path) && script_path != "") {
    script_dir <- if (file.exists(script_path)) dirname(script_path) else script_path
    cur <- script_dir
    for (i in seq_len(12)) {
      if (length(list.files(cur, pattern = "\\.Rproj$", ignore.case = TRUE)) > 0) {
        return(cur)
      }
      parent <- dirname(cur)
      if (identical(parent, cur)) break
      cur <- parent
    }
  }
  
  # 4) Manual cwd handling and upward search.
  wd <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  wd <- normalize_root_candidate(wd)
  if (!is.na(wd)) {
    cur <- wd
    for (i in seq_len(12)) {
      if (length(list.files(cur, pattern = "\\.Rproj$", ignore.case = TRUE)) > 0) {
        return(cur)
      }
      parent <- dirname(cur)
      if (identical(parent, cur)) break
      cur <- parent
    }
  }
  
  # 5) Explicit fallback requested by user.
  fallback_candidate <- normalize_root_candidate(fallback_root)
  if (!is.na(fallback_candidate)) return(fallback_candidate)
  
  # Final fallback to current working directory.
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

project_root <- detect_project_root()
raw_dir <- file.path(project_root, "data", "raw")
out_dir <- file.path(project_root, "outputs", "tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Detected project root: ", project_root)
message("Raw data directory: ", raw_dir)

if (!dir.exists(raw_dir)) {
  stop("Expected raw data directory not found: ", raw_dir,
       "\nCreate data/raw and place your Excel file(s) there.", call. = FALSE)
}

excel_files <- list.files(raw_dir, pattern = "\\.xls[xm]?$", full.names = TRUE, ignore.case = TRUE)
if (length(excel_files) == 0) {
  stop("No Excel files found in data/raw. Add at least one .xls/.xlsx/.xlsm file.", call. = FALSE)
}

clean_text <- function(x) {
  x %>%
    as.character() %>%
    stri_trans_general("Latin-ASCII") %>%
    str_to_lower() %>%
    str_replace_all("[^a-z0-9]+", " ") %>%
    str_squish()
}

matches_any <- function(x, patterns) {
  p <- paste(patterns, collapse = "|")
  str_detect(clean_text(x), p)
}

is_beech <- function(x) {
  txt <- clean_text(x)
  str_detect(txt, "fagus|grandifolia|beech|hetre|fag gr|f gra|f grandifolia|fagus grandifolia")
}

find_best_column <- function(df, patterns, required = TRUE, label = "column") {
  cn <- names(df)
  hits <- cn[matches_any(cn, patterns)]
  if (length(hits) == 0) {
    if (required) {
      stop("Could not detect a ", label, " column. Available columns: ",
           paste(cn, collapse = ", "), call. = FALSE)
    }
    return(NA_character_)
  }
  hits[[1]]
}

infer_species_column_from_values <- function(df) {
  char_cols <- names(df)[map_lgl(df, ~ is.character(.x) || is.factor(.x))]
  if (length(char_cols) == 0) return(NA_character_)
  
  scores <- map_dbl(char_cols, function(cl) {
    vals <- df[[cl]]
    vals <- vals[!is.na(vals)]
    if (length(vals) == 0) return(-Inf)
    vals <- as.character(vals)
    diversity <- dplyr::n_distinct(vals)
    beech_hits <- sum(is_beech(vals), na.rm = TRUE)
    diversity + beech_hits * 5
  })
  
  best <- char_cols[[which.max(scores)]]
  if (!is.finite(max(scores)) || max(scores) <= 1) return(NA_character_)
  best
}

as_numeric_safe <- function(x) {
  if (is.numeric(x)) return(x)
  readr::parse_number(as.character(x),
                      locale = locale(decimal_mark = ".", grouping_mark = ","),
                      na = c("", "NA", "N/A", "na", "null"))
}

read_sheet <- function(file, sheet) {
  read_excel(file, sheet = sheet) %>%
    janitor::clean_names(case = "snake")
}

score_sheet <- function(sheet_name, type = c("arbre", "gaule")) {
  type <- match.arg(type)
  nm <- clean_text(sheet_name)
  if (type == "arbre") {
    score <- 0
    if (str_detect(nm, "arbre")) score <- score + 3
    if (str_detect(nm, "feuille")) score <- score + 1
    if (str_detect(nm, "gaule")) score <- score - 3
    return(score)
  }
  score <- 0
  if (str_detect(nm, "gaule")) score <- score + 3
  if (str_detect(nm, "feuille")) score <- score + 1
  if (str_detect(nm, "arbre")) score <- score - 2
  score
}

find_target_sheet <- function(file, type = c("arbre", "gaule")) {
  type <- match.arg(type)
  sheets <- excel_sheets(file)
  scores <- map_dbl(sheets, score_sheet, type = type)
  if (all(scores <= 0)) return(NULL)
  tibble(file = file, sheet = sheets, score = scores) %>%
    arrange(desc(score), sheet) %>%
    slice(1)
}

all_sheets <- map_dfr(excel_files, function(f) {
  tibble(file = f, sheet = excel_sheets(f))
})

arbre_match <- map(excel_files, find_target_sheet, type = "arbre") %>% compact() %>% bind_rows() %>% arrange(desc(score)) %>% slice(1)
gaule_match <- map(excel_files, find_target_sheet, type = "gaule") %>% compact() %>% bind_rows() %>% arrange(desc(score)) %>% slice(1)

if (nrow(arbre_match) == 0) {
  stop("Could not identify an 'arbre'/'arbre feuille' sheet across data/raw Excel files.", call. = FALSE)
}
if (nrow(gaule_match) == 0) {
  stop("Could not identify a 'gaule'/'gaule feuille' sheet across data/raw Excel files.", call. = FALSE)
}

message("Detected arbre source: ", basename(arbre_match$file), " :: ", arbre_match$sheet)
message("Detected gaule source: ", basename(gaule_match$file), " :: ", gaule_match$sheet)

arbre <- read_sheet(arbre_match$file, arbre_match$sheet)
gaule <- read_sheet(gaule_match$file, gaule_match$sheet)

site_patterns <- c("^site$", "site", "station", "plot", "placette", "parcelle", "id_site", "site_id")
species_patterns <- c("species", "specie", "taxon", "espece", "essence", "ess", "latin", "genre", "genus", "nom_commun", "nom_latin")
dbh_patterns <- c("^dbh$", "dhp", "d hp", "diam", "diametre", "diameter", "d130", "circonference")
area_patterns <- c("area", "surface", "quadrat", "plot_area", "sampled_area", "superficie")
stem_count_patterns <- c("count", "n_stem", "stem_n", "stems", "nb", "nombre", "effectif")

summarize_basal_area <- function(df, dataset_name, has_species = TRUE) {
  site_col <- find_best_column(df, site_patterns, required = TRUE, label = paste(dataset_name, "site"))
  dbh_col <- find_best_column(df, dbh_patterns, required = TRUE, label = paste(dataset_name, "DBH/DHP/diameter"))
  species_col <- if (has_species) find_best_column(df, species_patterns, required = FALSE, label = paste(dataset_name, "species")) else NA_character_
  if (has_species && is.na(species_col)) {
    species_col <- infer_species_column_from_values(df)
  }
  area_col <- find_best_column(df, area_patterns, required = FALSE, label = paste(dataset_name, "area"))
  stem_count_col <- find_best_column(df, stem_count_patterns, required = FALSE, label = paste(dataset_name, "stem count"))
  
  if (is.na(species_col) && has_species) {
    warning("No species/taxon column detected in ", dataset_name,
            ". Beech-specific summaries from this table will be NA/0.")
  }
  
  dd <- df %>%
    mutate(
      site = as.character(.data[[site_col]]),
      dbh_cm = as_numeric_safe(.data[[dbh_col]]),
      stem_count = if (!is.na(stem_count_col)) as_numeric_safe(.data[[stem_count_col]]) else 1,
      stem_count = if_else(is.na(stem_count) | stem_count <= 0, 1, stem_count),
      sampled_area_m2_row = if (!is.na(area_col)) as_numeric_safe(.data[[area_col]]) else NA_real_,
      species_raw = if (!is.na(species_col)) as.character(.data[[species_col]]) else NA_character_,
      beech_flag = if (!is.na(species_col)) is_beech(species_raw) else TRUE
    ) %>%
    filter(!is.na(site), site != "", !is.na(dbh_cm), dbh_cm > 0) %>%
    mutate(ba_m2 = pi * (dbh_cm / 200)^2 * stem_count)
  
  area_by_site <- dd %>%
    group_by(site) %>%
    summarise(total_sampled_area_m2 = sum(sampled_area_m2_row, na.rm = TRUE),
              area_values_n = sum(!is.na(sampled_area_m2_row)), .groups = "drop") %>%
    mutate(has_valid_area = area_values_n > 0 & total_sampled_area_m2 > 0)
  
  if (any(area_by_site$has_valid_area)) {
    message("[", dataset_name, "] Some/all sites include sampling area. Basal area will be standardized to m2/ha for those sites.")
  } else {
    warning("[", dataset_name, "] No sampling area detected. Basal area is returned as raw summed m2.")
  }
  
  list(data = dd, area = area_by_site, site_col = site_col, dbh_col = dbh_col,
       species_col = species_col, area_col = area_col, stem_count_col = stem_count_col)
}

arbre_prep <- summarize_basal_area(arbre, "arbre", has_species = TRUE)
gaule_prep <- summarize_basal_area(gaule, "gaule", has_species = FALSE)

arbre_sum <- arbre_prep$data %>%
  group_by(site) %>%
  summarise(
    basal_area_all_trees_raw_m2 = sum(ba_m2, na.rm = TRUE),
    mean_dbh_all_trees_cm = mean(dbh_cm, na.rm = TRUE),
    beech_basal_area_arbre_raw_m2 = sum(ba_m2[beech_flag %in% TRUE], na.rm = TRUE),
    mean_beech_dbh_arbre_cm = ifelse(any(beech_flag %in% TRUE), mean(dbh_cm[beech_flag %in% TRUE], na.rm = TRUE), NA_real_),
    n_trees_used = n(),
    n_beech_trees_used = sum(beech_flag %in% TRUE, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(arbre_prep$area, by = "site") %>%
  mutate(
    basal_area_all_trees = if_else(has_valid_area,
                                   basal_area_all_trees_raw_m2 / total_sampled_area_m2 * 10000,
                                   basal_area_all_trees_raw_m2),
    beech_basal_area_arbre = if_else(has_valid_area,
                                     beech_basal_area_arbre_raw_m2 / total_sampled_area_m2 * 10000,
                                     beech_basal_area_arbre_raw_m2),
    basal_area_units = if_else(has_valid_area, "m2_ha", "raw_m2")
  ) %>%
  select(site,
         basal_area_all_trees,
         mean_dbh_all_trees_cm,
         beech_basal_area_arbre,
         mean_beech_dbh_arbre_cm,
         n_trees_used,
         n_beech_trees_used,
         basal_area_units,
         total_sampled_area_m2,
         basal_area_all_trees_raw_m2,
         beech_basal_area_arbre_raw_m2)

gaule_sum <- gaule_prep$data %>%
  group_by(site) %>%
  summarise(
    beech_gaule_basal_area_raw_m2 = sum(ba_m2, na.rm = TRUE),
    mean_gaule_dbh_cm = mean(dbh_cm, na.rm = TRUE),
    n_gaules_used = n(),
    .groups = "drop"
  ) %>%
  left_join(gaule_prep$area, by = "site") %>%
  mutate(
    beech_gaule_basal_area = if_else(has_valid_area,
                                     beech_gaule_basal_area_raw_m2 / total_sampled_area_m2 * 10000,
                                     beech_gaule_basal_area_raw_m2),
    gaule_basal_area_units = if_else(has_valid_area, "m2_ha", "raw_m2")
  ) %>%
  select(site,
         beech_gaule_basal_area,
         mean_gaule_dbh_cm,
         n_gaules_used,
         gaule_basal_area_units,
         total_sampled_area_m2,
         beech_gaule_basal_area_raw_m2)

combined <- full_join(arbre_sum, gaule_sum, by = "site", suffix = c("_arbre", "_gaule")) %>%
  mutate(
    basal_area_standardization = case_when(
      basal_area_units == "m2_ha" & gaule_basal_area_units == "m2_ha" ~ "both_standardized_m2_ha",
      basal_area_units == "m2_ha" & (is.na(gaule_basal_area_units) | gaule_basal_area_units == "raw_m2") ~ "arbre_standardized_gaule_raw",
      (is.na(basal_area_units) | basal_area_units == "raw_m2") & gaule_basal_area_units == "m2_ha" ~ "arbre_raw_gaule_standardized",
      TRUE ~ "raw_m2"
    )
  ) %>%
  arrange(site)

arbre_out_csv <- file.path(out_dir, "site_basal_area_arbre_summary.csv")
gaule_out_csv <- file.path(out_dir, "site_basal_area_gaule_summary.csv")
combined_out_csv <- file.path(out_dir, "site_basal_area_combined_summary.csv")
excel_out <- file.path(out_dir, "site_basal_area_summary_tables.xlsx")

write_csv(arbre_sum, arbre_out_csv, na = "")
write_csv(gaule_sum, gaule_out_csv, na = "")
write_csv(combined, combined_out_csv, na = "")

wb <- createWorkbook()
addWorksheet(wb, "arbre_summary")
addWorksheet(wb, "gaule_summary")
addWorksheet(wb, "combined_summary")
writeData(wb, "arbre_summary", arbre_sum)
writeData(wb, "gaule_summary", gaule_sum)
writeData(wb, "combined_summary", combined)
saveWorkbook(wb, excel_out, overwrite = TRUE)

message("Wrote CSV tables:")
message(" - ", arbre_out_csv)
message(" - ", gaule_out_csv)
message(" - ", combined_out_csv)
message("Wrote Excel workbook: ", excel_out)

print(combined)