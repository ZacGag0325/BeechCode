#!/usr/bin/env Rscript

# -----------------------------
# Purpose
# -----------------------------
# Format outputs/tables/site_description_table_article.csv as an article-ready
# Word table. This script is meant to be run after:
#   source("scripts/site_description_table_article.R")
#
# Important: flextable is optional. If flextable cannot be loaded on your
# machine, this script still creates a DOCX using officer::body_add_table().

# -----------------------------
# Package handling
# -----------------------------
required_packages <- c("readr", "dplyr", "officer")
optional_packages <- c("flextable")
auto_install_missing_packages <- TRUE
cran_repo <- "https://cran.rstudio.com"

install_if_missing <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0 && auto_install_missing_packages) {
    message("Missing package(s) detected: ", paste(missing, collapse = ", "))
    message("Attempting to install missing package(s) into the active R library: ", .libPaths()[1])
    utils::install.packages(missing, repos = cran_repo, dependencies = TRUE)
    missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  }
  missing
}

missing_required <- install_if_missing(required_packages)
if (length(missing_required) > 0) {
  stop(
    paste0(
      "Missing required package(s): ", paste(missing_required, collapse = ", "), "\n",
      "Please install them in the same R session/library used to run this script with:\n",
      "install.packages(c(\"", paste(missing_required, collapse = "\", \""), "\"), dependencies = TRUE)\n\n",
      "Current .libPaths():\n - ", paste(.libPaths(), collapse = "\n - ")
    ),
    call. = FALSE
  )
}

missing_optional <- install_if_missing(optional_packages)
has_flextable <- !("flextable" %in% missing_optional) && requireNamespace("flextable", quietly = TRUE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

message("DOCX package status: officer=", requireNamespace("officer", quietly = TRUE),
        ", flextable=", has_flextable)
if (!has_flextable) {
  message("flextable is not available in this R session, so an officer-only Word table will be created instead.")
}

# -----------------------------
# Helpers
# -----------------------------
resolve_project_root <- function() {
  script_path <- NA_character_
  
  frame_files <- vapply(sys.frames(), function(frame) {
    if (!is.null(frame$ofile)) return(as.character(frame$ofile))
    NA_character_
  }, character(1))
  frame_files <- frame_files[!is.na(frame_files) & frame_files != ""]
  if (length(frame_files) > 0) {
    script_path <- normalizePath(frame_files[length(frame_files)], winslash = "/", mustWork = FALSE)
  }
  
  if (is.na(script_path) || script_path == "") {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- args[grepl("^--file=", args)]
    if (length(file_arg) > 0) {
      script_path <- normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = FALSE)
    }
  }
  
  candidates <- unique(normalizePath(c(dirname(script_path), getwd()), winslash = "/", mustWork = FALSE))
  
  for (cand in candidates) {
    probe <- cand
    for (i in seq_len(8)) {
      if (dir.exists(file.path(probe, "scripts")) && dir.exists(file.path(probe, "outputs", "tables"))) {
        return(probe)
      }
      parent <- dirname(probe)
      if (identical(parent, probe)) break
      probe <- parent
    }
  }
  
  getwd()
}

format_number <- function(x, digits, na_text = "NA") {
  ifelse(is.na(x), na_text, sprintf(paste0("%.", digits, "f"), x))
}

normalize_col_name <- function(x) {
  x_norm <- iconv(x, from = "", to = "ASCII//TRANSLIT")
  x_norm <- tolower(x_norm)
  x_norm <- gsub("[^a-z0-9]+", "_", x_norm)
  x_norm <- gsub("_+", "_", x_norm)
  gsub("^_|_$", "", x_norm)
}

find_table_col <- function(df, candidates) {
  exact_idx <- match(candidates, names(df), nomatch = 0)
  exact_idx <- exact_idx[exact_idx > 0]
  if (length(exact_idx) > 0) return(names(df)[exact_idx[1]])
  
  names_norm <- normalize_col_name(names(df))
  candidate_norm <- normalize_col_name(candidates)
  norm_idx <- match(candidate_norm, names_norm, nomatch = 0)
  norm_idx <- norm_idx[norm_idx > 0]
  if (length(norm_idx) > 0) names(df)[norm_idx[1]] else NA_character_
}

require_table_col <- function(df, label, candidates) {
  col <- find_table_col(df, candidates)
  if (!is.na(col)) return(col)
  
  stop(
    "Missing required column for ", label, ". Accepted column names include: ",
    paste(candidates, collapse = ", "),
    "\nColumns found in CSV: ", paste(names(df), collapse = ", "),
    call. = FALSE
  )
}

format_table_column <- function(x, digits = NULL, na_text = "NA") {
  if (is.numeric(x)) {
    if (is.null(digits)) return(ifelse(is.na(x), na_text, as.character(x)))
    return(format_number(x, digits, na_text = na_text))
  }
  
  x_chr <- as.character(x)
  ifelse(is.na(x_chr) | x_chr == "", na_text, x_chr)
}

# -----------------------------
# Paths
# -----------------------------
project_root <- resolve_project_root()
table_dir <- file.path(project_root, "outputs", "tables")
csv_path <- file.path(table_dir, "site_description_table_article.csv")
docx_path <- file.path(table_dir, "site_description_table_article_formatted.docx")

if (!file.exists(csv_path)) {
  stop(
    "Input table not found: ", csv_path,
    "\nRun scripts/site_description_table_article.R first.",
    call. = FALSE
  )
}

# -----------------------------
# Read and validate input
# -----------------------------
site_table <- readr::read_csv(csv_path, show_col_types = FALSE)

# The upstream site_description_table_article.R script now writes an
# article-ready CSV with formatted column names such as
# "Basal area (m²·ha⁻¹)" and mean (SD) character values. Older CSVs used
# raw numeric names such as "Basal Area". Accept both shapes so this
# formatter can be run safely after either version of the table generator.
site_col <- require_table_col(site_table, "site", c("Site"))
latitude_col <- find_table_col(site_table, c("Latitude", "Lat"))
elevation_col <- require_table_col(site_table, "elevation", c("Elevation", "Elevation (m)", "elevation_mean"))
basal_area_col <- require_table_col(site_table, "basal area", c("Basal Area", "Basal_Area", "Basal area (m²·ha⁻¹)", "basal_area_mean"))
mean_dbh_col <- require_table_col(site_table, "mean DBH", c("Mean DBH", "Mean_DBH", "Mean DBH (cm)", "mean_dbh"))
beech_basal_area_col <- require_table_col(site_table, "beech basal area", c("Beech Basal Area", "Beech_Basal_Area", "Beech basal area (m²·ha⁻¹)", "beech_basal_area_mean"))
mean_beech_dbh_col <- require_table_col(site_table, "mean beech DBH", c("Mean Beech DBH", "Mean_Beech_DBH", "Mean beech DBH (cm)", "mean_beech_dbh"))
beech_sapling_basal_area_col <- require_table_col(site_table, "beech sapling basal area", c("Beech Sapling Basal Area", "Beech sapling basal area (m²·ha⁻¹)", "beech_sapling_basal_area_mean"))
dominant_species_col <- require_table_col(site_table, "dominant species", c("Dominant Species", "Top 3 dominant species", "Top 3 Dominant Species", "top_3_species_converted_codes"))

# -----------------------------
# Build article-ready table
# -----------------------------
site_table_pretty <- data.frame(
  Site = format_table_column(site_table[[site_col]]),
  check.names = FALSE
)

if (!is.na(latitude_col)) {
  site_table_pretty[["Latitude"]] <- format_table_column(site_table[[latitude_col]], digits = 5)
}

site_table_pretty[["Elevation\n(m)"]] <- format_table_column(site_table[[elevation_col]], digits = 0)
site_table_pretty[["Basal area\n(m²·ha⁻¹)"]] <- format_table_column(site_table[[basal_area_col]], digits = 1)
site_table_pretty[["Mean DBH\n(cm)"]] <- format_table_column(site_table[[mean_dbh_col]], digits = 1)
site_table_pretty[["Beech basal area\n(m²·ha⁻¹)"]] <- format_table_column(site_table[[beech_basal_area_col]], digits = 2)
site_table_pretty[["Mean beech DBH\n(cm)"]] <- format_table_column(site_table[[mean_beech_dbh_col]], digits = 1)
site_table_pretty[["Beech saplings basal area\n(m²·ha⁻¹)"]] <- format_table_column(site_table[[beech_sapling_basal_area_col]], digits = 2)
site_table_pretty[["Dominant species"]] <- format_table_column(site_table[[dominant_species_col]])

# -----------------------------
# Save Word table
# -----------------------------
if (has_flextable) {
  ft <- flextable::flextable(site_table_pretty) %>%
    flextable::theme_vanilla() %>%
    flextable::bold(part = "header") %>%
    flextable::bg(part = "header", bg = "#BFBFBF") %>%
    flextable::align(align = "center", part = "all") %>%
    flextable::valign(valign = "center", part = "all") %>%
    flextable::fontsize(size = 10, part = "all") %>%
    flextable::fontsize(size = 11, part = "header") %>%
    flextable::border_outer(border = officer::fp_border(color = "black", width = 1)) %>%
    flextable::border_inner_h(border = officer::fp_border(color = "black", width = 0.75)) %>%
    flextable::border_inner_v(border = officer::fp_border(color = "black", width = 0.75)) %>%
    flextable::autofit()
  
  print(ft)
  
  flextable::save_as_docx(
    "Table 1. Site characteristics" = ft,
    path = docx_path
  )
  
  message("Saved formatted Word table with flextable: ", docx_path)
} else {
  doc <- officer::read_docx()
  doc <- officer::body_add_par(doc, "Table 1. Site characteristics", style = "heading 1")
  doc <- officer::body_add_table(doc, value = site_table_pretty, style = "table_template")
  print(doc, target = docx_path)
  
  message("Saved formatted Word table with officer fallback: ", docx_path)
  message("Install/load flextable later if you want the grey-header flextable styling.")
}