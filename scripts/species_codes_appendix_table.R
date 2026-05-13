#!/usr/bin/env Rscript

# -----------------------------
# Purpose
# -----------------------------
# Build an appendix table defining every tree species code used in the
# "Top 3 dominant species" column of the site characteristics table.
#
# Outputs:
#   - outputs/tables/species_codes_appendix_table.csv
#   - outputs/tables/species_codes_appendix_table.docx

# -----------------------------
# Package handling
# -----------------------------
required_packages <- c("dplyr", "readr", "readxl", "stringr", "tidyr", "purrr", "officer")
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
  library(dplyr)
  library(readr)
  library(readxl)
  library(stringr)
  library(tidyr)
  library(purrr)
})

message("DOCX package status: officer=", requireNamespace("officer", quietly = TRUE),
        ", flextable=", has_flextable)
if (!has_flextable) {
  message("flextable is not available in this R session, so an officer-only Word table will be created instead.")
}

# -----------------------------
# Constants
# -----------------------------
caption <- "Appendix Table S1. Species codes used to describe dominant tree species in sampled stands. Scientific names are provided for each species code used in the site characteristics table."

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

normalize_code <- function(x) {
  str_trim(as.character(x)) |>
    str_to_upper()
}

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
      if (dir.exists(file.path(probe, "scripts")) && dir.exists(file.path(probe, "data"))) {
        return(probe)
      }
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
  preferred_idx <- which(str_detect(normalize_ascii(basename(files)), "donnee|modifie|west|summer|copie|raw|data"))
  chosen <- if (length(preferred_idx) > 0) files[preferred_idx[1]] else files[1]
  message("Excel files found in data/raw:")
  message(paste0(" - ", basename(files), collapse = "\n"))
  message("Selected workbook: ", basename(chosen))
  chosen
}

find_preferred_column <- function(df, candidates) {
  cn <- names(df)
  cn_norm <- normalize_ascii(cn)
  candidate_norm <- normalize_ascii(candidates)
  idx <- match(candidate_norm, cn_norm, nomatch = 0)
  idx <- idx[idx > 0]
  if (length(idx) == 0) NA_character_ else cn[idx[1]]
}

pick_column_regex <- function(df, patterns, label, required = FALSE, exclude = character()) {
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

match_sheet <- function(sheets, target) {
  idx <- which(normalize_ascii(sheets) == normalize_ascii(target))
  if (length(idx) == 0) NA_character_ else sheets[idx[1]]
}

read_sheet_quiet <- function(path, sheet) {
  suppressMessages(readxl::read_excel(path, sheet = sheet))
}

looks_like_species_lookup <- function(df) {
  cn_norm <- normalize_ascii(names(df))
  has_code <- any(str_detect(cn_norm, "species.*code|code.*species|species_code|code_fr|code_en|espece.*code|code.*espece"))
  has_name <- any(str_detect(cn_norm, "species.*name|name.*species|nom.*espece|espece.*nom|latin|scientific|scientifique"))
  has_code && has_name
}

find_species_lookup_sheet <- function(excel_path) {
  sheets <- readxl::excel_sheets(excel_path)
  exact <- match_sheet(sheets, "species_lookup")
  if (!is.na(exact)) return(exact)
  
  name_candidates <- sheets[str_detect(normalize_ascii(sheets), "species|espece|essence|lookup|code")]
  for (sheet in name_candidates) {
    preview <- tryCatch(read_sheet_quiet(excel_path, sheet), error = function(e) NULL)
    if (!is.null(preview) && looks_like_species_lookup(preview)) return(sheet)
  }
  
  for (sheet in sheets) {
    preview <- tryCatch(read_sheet_quiet(excel_path, sheet), error = function(e) NULL)
    if (!is.null(preview) && looks_like_species_lookup(preview)) return(sheet)
  }
  
  stop("Could not find a species lookup sheet in the raw Excel workbook. Available sheets: ", paste(sheets, collapse = ", "), call. = FALSE)
}

build_species_lookup <- function(species_lookup_raw) {
  code_en_col <- find_preferred_column(species_lookup_raw, c(
    "species_code_en", "code_en", "species_en", "espece_en", "english_code", "code_anglais"
  ))
  code_fr_col <- find_preferred_column(species_lookup_raw, c(
    "species_code_fr", "code_fr", "species_fr", "espece_fr", "french_code", "code_francais"
  ))
  name_en_col <- find_preferred_column(species_lookup_raw, c(
    "species_name_en", "name_en", "english_name", "nom_en", "nom_anglais", "common_name_en"
  ))
  name_fr_col <- find_preferred_column(species_lookup_raw, c(
    "species_name_fr", "name_fr", "french_name", "nom_fr", "nom_francais", "common_name_fr"
  ))
  latin_col <- find_preferred_column(species_lookup_raw, c(
    "species_name_latin", "latin_name", "scientific_name", "scientific_name_latin", "nom_latin", "nom_scientifique"
  ))
  
  if (is.na(code_en_col)) {
    code_en_col <- pick_column_regex(species_lookup_raw, c("species.*code.*en", "code.*en", "english.*code", "anglais.*code"), "English species code")
  }
  if (is.na(code_fr_col)) {
    code_fr_col <- pick_column_regex(species_lookup_raw, c("species.*code.*fr", "code.*fr", "french.*code", "francais.*code"), "French species code")
  }
  if (is.na(name_en_col)) {
    name_en_col <- pick_column_regex(species_lookup_raw, c("species.*name.*en", "name.*en", "english.*name", "nom.*anglais"), "English species name")
  }
  if (is.na(name_fr_col)) {
    name_fr_col <- pick_column_regex(species_lookup_raw, c("species.*name.*fr", "name.*fr", "french.*name", "nom.*francais"), "French species name")
  }
  if (is.na(latin_col)) {
    latin_col <- pick_column_regex(species_lookup_raw, c("latin", "scientific", "scientifique"), "scientific species name")
  }
  
  code_cols <- c(code_en_col, code_fr_col)
  code_cols <- code_cols[!is.na(code_cols)]
  if (length(code_cols) == 0) {
    stop(
      "Could not detect species code columns in the species lookup sheet. Detected columns: ",
      paste(names(species_lookup_raw), collapse = ", "),
      call. = FALSE
    )
  }
  
  lookup_base <- species_lookup_raw |>
    mutate(
      english_name = if (!is.na(name_en_col)) str_squish(as.character(.data[[name_en_col]])) else NA_character_,
      french_name = if (!is.na(name_fr_col)) str_squish(as.character(.data[[name_fr_col]])) else NA_character_,
      scientific_name = if (!is.na(latin_col)) str_squish(as.character(.data[[latin_col]])) else NA_character_
    )
  
  lookup_long <- map_dfr(code_cols, function(code_col) {
    lookup_base |>
      transmute(
        species_code = normalize_code(.data[[code_col]]),
        english_name,
        french_name,
        scientific_name,
        source_code_column = code_col
      )
  }) |>
    filter(!is.na(species_code), species_code != "") |>
    distinct(species_code, .keep_all = TRUE)
  
  if (nrow(lookup_long) == 0) stop("The species lookup sheet did not contain any usable species codes.", call. = FALSE)
  
  message("Species lookup columns used:")
  message(" - English code: ", ifelse(is.na(code_en_col), "not detected", code_en_col))
  message(" - French code: ", ifelse(is.na(code_fr_col), "not detected", code_fr_col))
  message(" - English name: ", ifelse(is.na(name_en_col), "not detected", name_en_col))
  message(" - French name: ", ifelse(is.na(name_fr_col), "not detected", name_fr_col))
  message(" - Scientific name: ", ifelse(is.na(latin_col), "not detected", latin_col))
  
  lookup_long
}

find_site_characteristics_table <- function(table_dir) {
  preferred <- c(
    file.path(table_dir, "site_description_table_article.csv"),
    file.path(table_dir, "site_characteristics_table.csv"),
    file.path(table_dir, "site_characteristics.csv")
  )
  existing_preferred <- preferred[file.exists(preferred)]
  if (length(existing_preferred) > 0) return(existing_preferred[1])
  
  candidates <- list.files(table_dir, pattern = "\\.csv$", full.names = TRUE, ignore.case = TRUE)
  candidates <- candidates[str_detect(normalize_ascii(basename(candidates)), "site.*(description|characteristic|char)")]
  if (length(candidates) > 0) return(candidates[1])
  
  stop(
    "Could not find a site characteristics CSV in ", table_dir, ".\n",
    "Run scripts/site_description_table_article.R first, then rerun this script.",
    call. = FALSE
  )
}

pick_dominant_species_column <- function(site_tbl) {
  exact_candidates <- c(
    "Top 3 dominant species", "Top 3 Dominant Species", "Top 3 dominant species", "Dominant Species", "Dominant species"
  )
  exact <- exact_candidates[exact_candidates %in% names(site_tbl)]
  if (length(exact) > 0) return(exact[1])
  
  cn <- names(site_tbl)
  cn_norm <- normalize_ascii(cn)
  idx <- which(str_detect(cn_norm, "top.*3.*dominant.*species|dominant.*species|top.*species"))
  if (length(idx) == 0) {
    stop(
      "Could not detect the 'Top 3 dominant species' column in the site characteristics table. Detected columns: ",
      paste(names(site_tbl), collapse = ", "),
      call. = FALSE
    )
  }
  cn[idx[1]]
}

split_species_codes <- function(x) {
  x |>
    as.character() |>
    str_split("\\s*,\\s*|\\s*;\\s*|\\s*/\\s*|\\s*\\|\\s*") |>
    unlist(use.names = FALSE) |>
    str_squish() |>
    str_replace_all("^[[:punct:]]+|[[:punct:]]+$", "") |>
    normalize_code()
}

xml_escape <- function(x) {
  x <- ifelse(is.na(x), "", as.character(x))
  x <- str_replace_all(x, "&", "&amp;")
  x <- str_replace_all(x, "<", "&lt;")
  x <- str_replace_all(x, ">", "&gt;")
  x <- str_replace_all(x, '"', "&quot;")
  x <- str_replace_all(x, "'", "&apos;")
  x
}

word_cell_xml <- function(value, bold = FALSE, italic = FALSE) {
  text <- xml_escape(value)
  bold_xml <- if (bold) "<w:b/>" else ""
  italic_xml <- if (italic) "<w:i/>" else ""
  run_pr <- if (bold || italic) paste0("<w:rPr>", bold_xml, italic_xml, "</w:rPr>") else ""
  paste0(
    "<w:tc>",
    "<w:tcPr><w:tcW w:w=\"2400\" w:type=\"dxa\"/></w:tcPr>",
    "<w:p><w:r>", run_pr, "<w:t>", text, "</w:t></w:r></w:p>",
    "</w:tc>"
  )
}

word_paragraph_xml <- function(value, bold = FALSE) {
  bold_xml <- if (bold) "<w:rPr><w:b/></w:rPr>" else ""
  paste0("<w:p><w:r>", bold_xml, "<w:t>", xml_escape(value), "</w:t></w:r></w:p>")
}

write_basic_docx_table <- function(df, path, title) {
  tmp_dir <- tempfile("species_codes_appendix_docx_")
  dir.create(file.path(tmp_dir, "_rels"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(tmp_dir, "word", "_rels"), recursive = TRUE, showWarnings = FALSE)
  
  scientific_col <- match("Scientific name", names(df))
  rows <- apply(df, 1, function(row) {
    cells <- vapply(seq_along(row), function(i) {
      word_cell_xml(row[[i]], italic = !is.na(scientific_col) && identical(i, scientific_col))
    }, character(1))
    paste0("<w:tr>", paste(cells, collapse = ""), "</w:tr>")
  })
  header <- paste0("<w:tr>", paste(vapply(names(df), word_cell_xml, character(1), bold = TRUE), collapse = ""), "</w:tr>")
  
  document_xml <- paste0(
    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
    '<w:document xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main" ',
    'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">',
    '<w:body>',
    word_paragraph_xml(title, bold = TRUE),
    '<w:tbl>',
    '<w:tblPr><w:tblStyle w:val="TableGrid"/><w:tblW w:w="0" w:type="auto"/>',
    '<w:tblBorders>',
    '<w:top w:val="single" w:sz="4" w:space="0" w:color="auto"/>',
    '<w:left w:val="single" w:sz="4" w:space="0" w:color="auto"/>',
    '<w:bottom w:val="single" w:sz="4" w:space="0" w:color="auto"/>',
    '<w:right w:val="single" w:sz="4" w:space="0" w:color="auto"/>',
    '<w:insideH w:val="single" w:sz="4" w:space="0" w:color="auto"/>',
    '<w:insideV w:val="single" w:sz="4" w:space="0" w:color="auto"/>',
    '</w:tblBorders></w:tblPr>',
    header,
    paste(rows, collapse = ""),
    '</w:tbl>',
    '<w:sectPr><w:pgSz w:w="15840" w:h="12240" w:orient="landscape"/>',
    '<w:pgMar w:top="720" w:right="720" w:bottom="720" w:left="720" w:header="360" w:footer="360" w:gutter="0"/></w:sectPr>',
    '</w:body></w:document>'
  )
  
  writeLines(c(
    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
    '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">',
    '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>',
    '<Default Extension="xml" ContentType="application/xml"/>',
    '<Override PartName="/word/document.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.document.main+xml"/>',
    '</Types>'
  ), file.path(tmp_dir, "[Content_Types].xml"), useBytes = TRUE)
  writeLines(c(
    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
    '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">',
    '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="word/document.xml"/>',
    '</Relationships>'
  ), file.path(tmp_dir, "_rels", ".rels"), useBytes = TRUE)
  writeLines(c(
    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
    '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"></Relationships>'
  ), file.path(tmp_dir, "word", "_rels", "document.xml.rels"), useBytes = TRUE)
  writeLines(document_xml, file.path(tmp_dir, "word", "document.xml"), useBytes = TRUE)
  
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(tmp_dir)
  if (file.exists(path)) unlink(path)
  utils::zip(zipfile = path, files = list.files(tmp_dir, recursive = TRUE, all.files = TRUE, no.. = TRUE), flags = "-q")
  if (!file.exists(path)) stop("Failed to create Word document at: ", path, call. = FALSE)
  invisible(path)
}

write_species_appendix_docx <- function(df, path, caption) {
  docx_status <- tryCatch(
    {
      if (has_flextable) {
        ft <- flextable::flextable(df) |>
          flextable::theme_booktabs() |>
          flextable::autofit() |>
          flextable::align(align = "center", part = "all") |>
          flextable::align(j = "Species code", align = "left", part = "body") |>
          flextable::align(j = c("English name", "French name", "Scientific name"), align = "left", part = "body") |>
          flextable::bold(part = "header") |>
          flextable::compose(
            j = "Scientific name",
            value = flextable::as_paragraph(flextable::as_i(`Scientific name`)),
            part = "body"
          )
        
        doc <- officer::read_docx()
        doc <- officer::body_add_par(doc, caption, style = "heading 1")
        doc <- officer::body_add_flextable(doc, ft)
        print(doc, target = path)
        "created with flextable"
      } else {
        doc <- officer::read_docx()
        doc <- officer::body_add_par(doc, caption, style = "heading 1")
        doc <- officer::body_add_table(doc, value = as.data.frame(df), style = "table_template")
        print(doc, target = path)
        "created with officer fallback"
      }
    },
    error = function(e) {
      message(
        "Package-based Word export failed: ",
        conditionMessage(e),
        ". Creating a basic Word-compatible .docx file without extra packages instead."
      )
      write_basic_docx_table(df, path, title = caption)
      "created with built-in docx fallback after package export failed"
    }
  )
  
  if (!file.exists(path)) {
    stop("Word appendix table was not created: ", path, call. = FALSE)
  }
  
  docx_status
}

# -----------------------------
# Paths
# -----------------------------
project_root <- resolve_project_root()
raw_dir <- file.path(project_root, "data", "raw")
out_dir <- file.path(project_root, "outputs", "tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

site_csv_path <- find_site_characteristics_table(out_dir)
csv_path <- file.path(out_dir, "species_codes_appendix_table.csv")
docx_path <- file.path(out_dir, "species_codes_appendix_table.docx")

# -----------------------------
# Read species lookup from raw Excel file
# -----------------------------
excel_path <- find_raw_workbook(raw_dir)
species_lookup_sheet <- find_species_lookup_sheet(excel_path)
message("Species lookup sheet detected: ", species_lookup_sheet)
species_lookup_raw <- read_sheet_quiet(excel_path, species_lookup_sheet)
species_lookup <- build_species_lookup(species_lookup_raw)

# -----------------------------
# Read site characteristics table and split Top 3 dominant species codes
# -----------------------------
message("Site characteristics table detected: ", site_csv_path)
site_tbl <- readr::read_csv(site_csv_path, show_col_types = FALSE)
dominant_species_col <- pick_dominant_species_column(site_tbl)
message("Dominant species column used: ", dominant_species_col)

species_codes_used <- split_species_codes(site_tbl[[dominant_species_col]])
species_codes_used <- sort(unique(species_codes_used[!is.na(species_codes_used) & species_codes_used != "" & species_codes_used != "NA"]))

if (length(species_codes_used) == 0) {
  stop("No species codes were found in column '", dominant_species_col, "'.", call. = FALSE)
}

message("Unique species codes found in site characteristics table: ", paste(species_codes_used, collapse = ", "))

# -----------------------------
# Join used species codes to lookup table
# -----------------------------
appendix_table_joined <- tibble(`Species code` = species_codes_used) |>
  left_join(
    species_lookup |>
      select(species_code, english_name, french_name, scientific_name),
    by = c("Species code" = "species_code")
  ) |>
  mutate(found_in_lookup = !is.na(english_name) | !is.na(french_name) | !is.na(scientific_name))

missing_species_codes <- appendix_table_joined |>
  filter(!found_in_lookup) |>
  pull(`Species code`)

if (length(missing_species_codes) > 0) {
  warning(
    "Species code(s) used in the site characteristics table were missing from the lookup table: ",
    paste(missing_species_codes, collapse = ", "),
    if (any(c("BL", "TA") %in% missing_species_codes)) " (includes BL and/or TA)" else "",
    call. = FALSE
  )
}

incomplete_species_codes <- appendix_table_joined |>
  filter(found_in_lookup, if_any(c(english_name, french_name, scientific_name), is.na)) |>
  pull(`Species code`)

if (length(incomplete_species_codes) > 0) {
  warning(
    "Species code(s) were found in the lookup table but are missing one or more name fields: ",
    paste(incomplete_species_codes, collapse = ", "),
    call. = FALSE
  )
}

appendix_table <- appendix_table_joined |>
  transmute(
    `Species code`,
    `English name` = english_name,
    `French name` = french_name,
    `Scientific name` = scientific_name
  )

# -----------------------------
# Save outputs
# -----------------------------
readr::write_csv(appendix_table, csv_path, na = "")
docx_status <- write_species_appendix_docx(appendix_table, docx_path, caption)

# -----------------------------
# Console output
# -----------------------------
message("\nFinal species codes appendix table:")
print(appendix_table, n = Inf, width = Inf)

message("\nSaved:")
message(" - ", csv_path)
message(" - ", docx_path, " (", docx_status, ")")
message("Species codes appendix table was saved.")