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

required_cols <- c(
  "Site",
  "Basal Area",
  "Mean DBH",
  "Beech Basal Area",
  "Mean Beech DBH",
  "Beech Sapling Basal Area",
  "Dominant Species",
  "Elevation"
)
missing_cols <- setdiff(required_cols, names(site_table))
if (length(missing_cols) > 0) {
  stop("Missing required columns in CSV: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

# -----------------------------
# Build article-ready table
# -----------------------------
site_table_pretty <- site_table %>%
  mutate(
    `Basal area\n(m²·ha⁻¹)` = format_number(`Basal Area`, 1),
    `Mean DBH\n(cm)` = format_number(`Mean DBH`, 1),
    `Beech basal area\n(m²·ha⁻¹)` = format_number(`Beech Basal Area`, 2),
    `Mean beech DBH\n(cm)` = format_number(`Mean Beech DBH`, 1),
    `Beech saplings basal area\n(m²·ha⁻¹)` = format_number(`Beech Sapling Basal Area`, 2),
    `Dominant species` = `Dominant Species`,
    `Elevation\n(m)` = format_number(Elevation, 0)
  ) %>%
  select(
    Site,
    `Basal area\n(m²·ha⁻¹)`,
    `Mean DBH\n(cm)`,
    `Beech basal area\n(m²·ha⁻¹)`,
    `Mean beech DBH\n(cm)`,
    `Beech saplings basal area\n(m²·ha⁻¹)`,
    `Dominant species`,
    `Elevation\n(m)`
  )

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