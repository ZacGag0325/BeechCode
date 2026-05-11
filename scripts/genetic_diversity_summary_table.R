#!/usr/bin/env Rscript

############################################################
# Genetic diversity summary table for thesis/article use
#
# Standalone script. Loads the final cleaned, clone-corrected
# microsatellite object from outputs/v1/objects through the shared
# BeechCode loader, calculates site-level diversity statistics, and
# exports article-ready CSV/Excel/Word tables.
#
# Main outputs:
# - outputs/tables/genetic_diversity_summary_table.csv
# - outputs/tables/genetic_diversity_summary_table.xlsx (if openxlsx is available)
# - outputs/tables/genetic_diversity_summary_table.docx
############################################################

required_packages <- c("adegenet", "hierfstat", "dplyr", "tidyr", "stringr", "readxl")
missing_required <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_required) > 0) {
  stop(
    "[genetic_diversity_summary_table] Missing required package(s): ",
    paste(missing_required, collapse = ", "),
    ". Please install them before running this script.",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(adegenet)
  library(hierfstat)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readxl)
})

# Reuse the project root, object-loading, validation, and output-directory
# strategy used by the main population-genetic scripts.
source("scripts/_load_objects.R")

SCRIPT_TAG <- "[genetic_diversity_summary_table]"
GI_MLL_FILE <- file.path(OBJ_DIR, "gi_mll.rds")
DF_IDS_FILE <- file.path(OBJ_DIR, "df_ids.rds")
META_FILE <- file.path(OBJ_DIR, "meta.rds")

message(SCRIPT_TAG, " Using final cleaned clone-corrected object: ", GI_MLL_FILE)
message(SCRIPT_TAG, " Metadata loaded from: ", DF_IDS_FILE)
message(SCRIPT_TAG, " Pipeline metadata loaded from: ", META_FILE)

# -----------------------------
# Optional-package status
# -----------------------------
has_officer <- requireNamespace("officer", quietly = TRUE)
has_flextable <- requireNamespace("flextable", quietly = TRUE)
has_openxlsx <- requireNamespace("openxlsx", quietly = TRUE)

# -----------------------------
# Helpers
# -----------------------------
normalize_lookup_key <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("\uFEFF", "", x, fixed = TRUE)
  x <- gsub("[[:cntrl:]]", "", x)
  x[is.na(x)] <- ""
  x
}

normalize_ascii <- function(x) {
  x %>%
    iconv(from = "", to = "ASCII//TRANSLIT") %>%
    tolower() %>%
    stringr::str_replace_all("[^a-z0-9]+", "_") %>%
    stringr::str_replace_all("_+", "_") %>%
    stringr::str_replace_all("^_|_$", "")
}

pick_column <- function(df, choices, label = "column", required = FALSE) {
  nms <- names(df)
  nms_norm <- normalize_ascii(nms)
  choices_norm <- normalize_ascii(choices)
  idx <- match(TRUE, nms_norm %in% choices_norm, nomatch = 0)
  if (idx == 0) {
    if (required) {
      stop(SCRIPT_TAG, " Could not find required ", label, ". Accepted names: ", paste(choices, collapse = ", "), call. = FALSE)
    }
    return(NA_character_)
  }
  nms[idx]
}

safe_mean <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (length(x) == 0 || all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

round3 <- function(x) round(as.numeric(x), 3)

format3 <- function(x) {
  x_num <- suppressWarnings(as.numeric(x))
  ifelse(is.na(x_num), NA_character_, formatC(round(x_num, 3), format = "f", digits = 3))
}

natural_site_rank <- function(x) {
  x_chr <- as.character(x)
  prefix <- stringr::str_extract(x_chr, "^[A-Za-z]+")
  num <- suppressWarnings(as.integer(stringr::str_extract(x_chr, "[0-9]+")))
  prefix_rank <- dplyr::case_when(
    toupper(prefix) == "S" ~ 1L,
    toupper(prefix) == "N" ~ 2L,
    TRUE ~ 99L
  )
  ifelse(is.na(num), 9999L, prefix_rank * 100L + num)
}

# -----------------------------
# site_lookup loading and site labels/order
# -----------------------------
SITE_LOOKUP_SHEET <- "site_lookup"

find_site_lookup_workbook <- function(raw_dir = file.path(PROJECT_ROOT, "data", "raw"), sheet = SITE_LOOKUP_SHEET) {
  if (!dir.exists(raw_dir)) {
    message(SCRIPT_TAG, " site_lookup unavailable because data/raw does not exist: ", raw_dir)
    return(NULL)
  }
  excel_files <- list.files(raw_dir, pattern = "\\.(xlsx|xls)$", full.names = TRUE, ignore.case = TRUE)
  if (length(excel_files) == 0) {
    message(SCRIPT_TAG, " site_lookup unavailable: no Excel workbook was found in data/raw.")
    return(NULL)
  }
  
  has_lookup <- vapply(
    excel_files,
    function(path) {
      sheets <- tryCatch(readxl::excel_sheets(path), error = function(e) character(0))
      any(normalize_ascii(sheets) == normalize_ascii(sheet))
    },
    logical(1)
  )
  lookup_files <- excel_files[has_lookup]
  
  if (length(lookup_files) == 0) {
    message(
      SCRIPT_TAG, " site_lookup unavailable: no workbook in data/raw contains a '", sheet,
      "' sheet. Checked: ", paste(basename(excel_files), collapse = ", ")
    )
    return(NULL)
  }
  if (length(lookup_files) > 1) {
    message(
      SCRIPT_TAG, " Multiple workbooks contain site_lookup; using ", basename(lookup_files[1]),
      ". Other matches: ", paste(basename(lookup_files[-1]), collapse = ", ")
    )
  } else {
    message(SCRIPT_TAG, " Reading site_lookup from: ", basename(lookup_files[1]))
  }
  
  lookup_files[1]
}

load_site_lookup <- function() {
  workbook <- find_site_lookup_workbook()
  if (is.null(workbook)) return(NULL)
  
  sheets <- readxl::excel_sheets(workbook)
  sheet <- sheets[match(normalize_ascii(SITE_LOOKUP_SHEET), normalize_ascii(sheets))]
  lookup <- suppressMessages(readxl::read_excel(workbook, sheet = sheet)) %>%
    as.data.frame(stringsAsFactors = FALSE)
  
  old_site_col <- pick_column(lookup, c("Site", "site", "site_code", "old_site", "code_site"), "old site code", required = TRUE)
  label_col <- pick_column(lookup, c("Site_label", "site_label", "new_site", "site_new", "label", "site_id"), "new site label")
  region_col <- pick_column(lookup, c("Region", "region"), "region")
  order_col <- pick_column(lookup, c("Site_order", "site_order", "order", "ordre", "sort", "south_north_order"), "south-to-north order")
  latitude_col <- pick_column(lookup, c("Latitude", "latitude", "lat", "LAT", "Lat", "y", "Y"), "latitude")
  
  if (is.na(label_col)) label_col <- old_site_col
  
  out <- lookup %>%
    mutate(
      old_site = normalize_lookup_key(.data[[old_site_col]]),
      site_label = normalize_lookup_key(.data[[label_col]]),
      region = if (!is.na(region_col)) normalize_lookup_key(.data[[region_col]]) else NA_character_,
      site_order = if (!is.na(order_col)) suppressWarnings(as.numeric(.data[[order_col]])) else NA_real_,
      latitude = if (!is.na(latitude_col)) suppressWarnings(as.numeric(stringr::str_replace_all(as.character(.data[[latitude_col]]), ",", "."))) else NA_real_
    ) %>%
    filter(nzchar(old_site), nzchar(site_label)) %>%
    distinct(old_site, .keep_all = TRUE)
  
  if (nrow(out) == 0) {
    message(SCRIPT_TAG, " site_lookup was found but no usable site-code/site-label rows were detected; old site codes will be retained.")
    return(NULL)
  }
  if (anyDuplicated(out$old_site)) {
    stop(SCRIPT_TAG, " site_lookup contains duplicated old site codes after cleaning.", call. = FALSE)
  }
  
  message(SCRIPT_TAG, " Loaded site_lookup rows: ", nrow(out))
  out
}

find_latitude_col <- function(df) {
  pick_column(df, c("Latitude", "latitude", "lat", "LAT", "Lat", "y", "Y"), "latitude")
}

site_lookup <- load_site_lookup()

map_site_labels <- function(site_values) {
  site_values <- normalize_lookup_key(site_values)
  if (is.null(site_lookup)) return(site_values)
  idx <- match(site_values, site_lookup$old_site)
  labels <- site_lookup$site_label[idx]
  ifelse(is.na(labels) | !nzchar(labels), site_values, labels)
}

build_site_order_table <- function(original_site_values, display_site_values, df_ids_tbl = NULL) {
  out <- data.frame(
    original_site = as.character(original_site_values),
    Site = as.character(display_site_values),
    stringsAsFactors = FALSE
  ) %>%
    distinct(original_site, Site)
  
  if (!is.null(site_lookup)) {
    lookup_order <- site_lookup %>%
      transmute(
        original_site = old_site,
        Site = site_label,
        lookup_order = site_order,
        latitude = latitude
      ) %>%
      distinct(original_site, .keep_all = TRUE)
    out <- out %>%
      left_join(lookup_order, by = c("original_site", "Site"))
  } else {
    out$lookup_order <- NA_real_
    out$latitude <- NA_real_
  }
  
  if (!any(!is.na(out$latitude)) && !is.null(df_ids_tbl)) {
    lat_col <- find_latitude_col(df_ids_tbl)
    if (!is.na(lat_col)) {
      df_ids_cols_local <- resolve_df_ids_columns(df_ids_tbl, context = SCRIPT_TAG, require = TRUE)
      lat_by_site <- df_ids_tbl %>%
        transmute(
          original_site = normalize_lookup_key(.data[[df_ids_cols_local$site_col]]),
          latitude_from_df_ids = suppressWarnings(as.numeric(stringr::str_replace_all(as.character(.data[[lat_col]]), ",", ".")))
        ) %>%
        filter(nzchar(original_site)) %>%
        group_by(original_site) %>%
        summarise(latitude_from_df_ids = safe_mean(latitude_from_df_ids), .groups = "drop")
      out <- out %>%
        left_join(lat_by_site, by = "original_site") %>%
        mutate(latitude = ifelse(is.na(latitude), latitude_from_df_ids, latitude)) %>%
        select(-latitude_from_df_ids)
    }
  }
  
  if (any(!is.na(out$lookup_order))) {
    message(SCRIPT_TAG, " Ordering sites south-to-north using site_lookup Site_order.")
    out <- out %>% arrange(is.na(lookup_order), lookup_order, natural_site_rank(Site), Site)
  } else if (any(!is.na(out$latitude))) {
    message(SCRIPT_TAG, " Ordering sites south-to-north using latitude.")
    out <- out %>% arrange(is.na(latitude), latitude, natural_site_rank(Site), Site)
  } else {
    message(SCRIPT_TAG, " Ordering sites using natural S1-S6/N1-N6 fallback order.")
    out <- out %>% arrange(natural_site_rank(Site), Site)
  }
  
  out %>% mutate(display_order = row_number())
}

# -----------------------------
# Final cleaned genetic dataset
# -----------------------------
if (!exists("gi_mll") || !inherits(gi_mll, "genind")) {
  stop(SCRIPT_TAG, " gi_mll was not loaded as a genind object. Run scripts/00_master_pipeline.R first.", call. = FALSE)
}
if (!exists("df_ids_mll") || !is.data.frame(df_ids_mll)) {
  stop(SCRIPT_TAG, " df_ids_mll metadata was not created by scripts/_load_objects.R.", call. = FALSE)
}
if (adegenet::nInd(gi_mll) != nrow(df_ids_mll)) {
  stop(SCRIPT_TAG, " nInd(gi_mll) does not match nrow(df_ids_mll).", call. = FALSE)
}

original_sites <- normalize_lookup_key(as.character(adegenet::pop(gi_mll)))
if (length(original_sites) != adegenet::nInd(gi_mll) || any(!nzchar(original_sites))) {
  stop(SCRIPT_TAG, " gi_mll does not contain a complete site/population assignment.", call. = FALSE)
}

display_sites <- map_site_labels(original_sites)
missing_lookup_sites <- if (!is.null(site_lookup)) setdiff(unique(original_sites), site_lookup$old_site) else character(0)
if (length(missing_lookup_sites) > 0) {
  message(SCRIPT_TAG, " Site code(s) not found in site_lookup and retained as-is: ", paste(missing_lookup_sites, collapse = ", "))
}

# Work on a copy so the sourced gi_mll object remains intact for any interactive session.
gen_obj <- gi_mll
adegenet::pop(gen_obj) <- as.factor(display_sites)
site_order <- build_site_order_table(unique(original_sites), unique(display_sites), df_ids_tbl = df_ids_mll)

message(SCRIPT_TAG, " Individuals retained after filtering: ", adegenet::nInd(gen_obj))
message(SCRIPT_TAG, " Loci retained after filtering: ", adegenet::nLoc(gen_obj))
message(SCRIPT_TAG, " Final site labels: ", paste(site_order$Site, collapse = ", "))

# -----------------------------
# Diversity calculations
# -----------------------------
pop_sizes <- table(as.character(adegenet::pop(gen_obj)))
if (any(pop_sizes < 2)) {
  stop(
    SCRIPT_TAG, " At least one site has fewer than two retained individuals (",
    paste(names(pop_sizes)[pop_sizes < 2], collapse = ", "),
    "); Ho/He/FIS and rarefied Ar cannot be calculated reliably.",
    call. = FALSE
  )
}

hf <- hierfstat::genind2hierfstat(gen_obj)
bs <- hierfstat::basic.stats(hf)

if (is.null(bs$Ho) || is.null(bs$Hs) || is.null(bs$Fis)) {
  stop(SCRIPT_TAG, " hierfstat::basic.stats did not return Ho, Hs, and Fis matrices.", call. = FALSE)
}

site_levels <- rownames(bs$Ho)
site_n_tbl <- table(as.character(adegenet::pop(gen_obj)))

allele_tab <- adegenet::tab(gen_obj, NA.method = "asis")
loc_fac <- adegenet::locFac(gen_obj)
loci <- adegenet::locNames(gen_obj)
site_row_indices <- split(seq_len(nrow(allele_tab)), as.character(adegenet::pop(gen_obj)))

count_alleles_for_site_locus <- function(rows, locus_name) {
  cols <- which(loc_fac == locus_name)
  if (length(cols) == 0 || length(rows) == 0) return(NA_integer_)
  mat <- allele_tab[rows, cols, drop = FALSE]
  if (all(is.na(mat))) return(NA_integer_)
  as.integer(sum(colSums(mat, na.rm = TRUE) > 0))
}

mean_alleles_by_site <- tidyr::expand_grid(
  Site = names(site_row_indices),
  Locus = loci
) %>%
  rowwise() %>%
  mutate(Allele_count = count_alleles_for_site_locus(site_row_indices[[Site]], Locus)) %>%
  ungroup() %>%
  group_by(Site) %>%
  summarise(Na = safe_mean(Allele_count), .groups = "drop")

min_n <- min(as.integer(pop_sizes), na.rm = TRUE)
message(SCRIPT_TAG, " Rarefied allelic richness sample size (min.n): ", min_n)

ar_obj <- hierfstat::allelic.richness(hf, min.n = min_n)
if (is.null(ar_obj$Ar)) {
  stop(SCRIPT_TAG, " hierfstat::allelic.richness did not return an Ar matrix.", call. = FALSE)
}

ar_mat <- ar_obj$Ar
if (all(site_levels %in% colnames(ar_mat))) {
  allelic_richness_by_site <- data.frame(
    Site = colnames(ar_mat),
    Ar = as.numeric(colMeans(ar_mat, na.rm = TRUE)),
    stringsAsFactors = FALSE
  )
} else if (all(site_levels %in% rownames(ar_mat))) {
  allelic_richness_by_site <- data.frame(
    Site = rownames(ar_mat),
    Ar = as.numeric(rowMeans(ar_mat, na.rm = TRUE)),
    stringsAsFactors = FALSE
  )
} else {
  stop(SCRIPT_TAG, " Could not match allelic.richness Ar matrix dimensions to site labels.", call. = FALSE)
}

heterozygosity_by_site <- data.frame(
  Site = site_levels,
  N = as.integer(site_n_tbl[site_levels]),
  Ho = apply(bs$Ho, 1, safe_mean),
  He = apply(bs$Hs, 1, safe_mean),
  FIS = apply(bs$Fis, 1, safe_mean),
  stringsAsFactors = FALSE
)

genetic_diversity_summary_table <- heterozygosity_by_site %>%
  left_join(mean_alleles_by_site, by = "Site") %>%
  left_join(allelic_richness_by_site, by = "Site") %>%
  left_join(site_order %>% select(Site, display_order), by = "Site") %>%
  arrange(display_order, Site) %>%
  transmute(
    Site = Site,
    N = as.integer(N),
    Na = round3(Na),
    Ar = round3(Ar),
    Ho = round3(Ho),
    He = round3(He),
    FIS = round3(FIS)
  )

genetic_diversity_summary_table_word <- genetic_diversity_summary_table %>%
  mutate(
    Na = format3(Na),
    Ar = format3(Ar),
    Ho = format3(Ho),
    He = format3(He),
    FIS = format3(FIS)
  )

# -----------------------------
# Export helpers
# -----------------------------
write_genetic_diversity_xlsx <- function(summary_tbl, path) {
  if (!has_openxlsx) {
    message(SCRIPT_TAG, " openxlsx is not available; skipping optional Excel export.")
    return(NA_character_)
  }
  
  wb <- openxlsx::createWorkbook()
  sheet <- "Genetic diversity"
  openxlsx::addWorksheet(wb, sheet)
  openxlsx::writeData(wb, sheet, summary_tbl)
  header_style <- openxlsx::createStyle(
    textDecoration = "bold",
    halign = "center",
    valign = "center",
    border = "bottom",
    fgFill = "#D9EAD3"
  )
  body_style <- openxlsx::createStyle(halign = "center", valign = "center")
  openxlsx::addStyle(wb, sheet, header_style, rows = 1, cols = seq_len(ncol(summary_tbl)), gridExpand = TRUE)
  if (nrow(summary_tbl) > 0) {
    openxlsx::addStyle(wb, sheet, body_style, rows = 2:(nrow(summary_tbl) + 1), cols = seq_len(ncol(summary_tbl)), gridExpand = TRUE)
  }
  openxlsx::freezePane(wb, sheet, firstRow = TRUE)
  openxlsx::setColWidths(wb, sheet, cols = seq_len(ncol(summary_tbl)), widths = "auto")
  openxlsx::saveWorkbook(wb, path, overwrite = TRUE)
  path
}

xml_escape_word <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub('"', "&quot;", x, fixed = TRUE)
  x <- gsub("'", "&apos;", x, fixed = TRUE)
  x
}

word_cell_xml <- function(value, bold = FALSE) {
  bold_xml <- if (bold) "<w:rPr><w:b/></w:rPr>" else ""
  paste0(
    "<w:tc>",
    "<w:tcPr><w:tcW w:w=\"1800\" w:type=\"dxa\"/></w:tcPr>",
    "<w:p><w:r>", bold_xml, "<w:t>", xml_escape_word(value), "</w:t></w:r></w:p>",
    "</w:tc>"
  )
}

word_paragraph_xml <- function(value, bold = FALSE) {
  bold_xml <- if (bold) "<w:rPr><w:b/></w:rPr>" else ""
  paste0("<w:p><w:r>", bold_xml, "<w:t>", xml_escape_word(value), "</w:t></w:r></w:p>")
}

write_basic_docx_table <- function(df, path, title, note = NULL) {
  tmp_dir <- tempfile("genetic_diversity_summary_docx_")
  dir.create(file.path(tmp_dir, "_rels"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(tmp_dir, "word", "_rels"), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  
  df <- as.data.frame(lapply(df, as.character), stringsAsFactors = FALSE, check.names = FALSE)
  rows <- apply(df, 1, function(row) {
    paste0("<w:tr>", paste(vapply(row, word_cell_xml, character(1)), collapse = ""), "</w:tr>")
  })
  header <- paste0("<w:tr>", paste(vapply(names(df), word_cell_xml, character(1), bold = TRUE), collapse = ""), "</w:tr>")
  note_xml <- if (!is.null(note) && nzchar(note)) word_paragraph_xml(note) else ""
  
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
    note_xml,
    '<w:sectPr><w:pgSz w:w="12240" w:h="15840"/>',
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
  utils::zip(
    zipfile = path,
    files = list.files(tmp_dir, recursive = TRUE, all.files = TRUE, no.. = TRUE),
    flags = "-q"
  )
  if (!file.exists(path)) stop(SCRIPT_TAG, " Failed to create Word document at: ", path, call. = FALSE)
  invisible(path)
}

write_genetic_diversity_docx <- function(summary_tbl, path, title, note) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  
  docx_status <- tryCatch(
    {
      if (has_officer && has_flextable) {
        ft <- flextable::flextable(summary_tbl) %>%
          flextable::theme_booktabs() %>%
          flextable::autofit() %>%
          flextable::align(align = "center", part = "all") %>%
          flextable::align(j = "Site", align = "left", part = "body") %>%
          flextable::bold(part = "header") %>%
          flextable::fontsize(size = 10, part = "all") %>%
          flextable::fontsize(size = 11, part = "header")
        
        doc <- officer::read_docx()
        doc <- officer::body_add_par(doc, title, style = "heading 1")
        doc <- officer::body_add_flextable(doc, ft)
        doc <- officer::body_add_par(doc, note, style = "Normal")
        print(doc, target = path)
        "created with flextable"
      } else if (has_officer) {
        message(SCRIPT_TAG, " flextable is not available/loadable; creating the Word table with officer::body_add_table instead.")
        doc <- officer::read_docx()
        doc <- officer::body_add_par(doc, title, style = "heading 1")
        doc <- officer::body_add_table(doc, value = as.data.frame(summary_tbl), style = "table_template")
        doc <- officer::body_add_par(doc, note, style = "Normal")
        print(doc, target = path)
        "created with officer fallback"
      } else {
        message(SCRIPT_TAG, " Neither officer nor flextable is loadable; creating a basic Word-compatible .docx file instead.")
        write_basic_docx_table(summary_tbl, path, title = title, note = note)
        "created with built-in docx fallback"
      }
    },
    error = function(e) {
      message(SCRIPT_TAG, " Package-based Word export failed: ", conditionMessage(e), ". Creating a basic Word-compatible .docx file instead.")
      write_basic_docx_table(summary_tbl, path, title = title, note = note)
      "created with built-in docx fallback after package export failed"
    }
  )
  
  if (!file.exists(path)) {
    stop(SCRIPT_TAG, " Word table was not created: ", path, call. = FALSE)
  }
  
  message(SCRIPT_TAG, " Saved Word table: ", path, " (", docx_status, ")")
  invisible(path)
}

# -----------------------------
# Save outputs
# -----------------------------
main_csv <- file.path(TABLES_DIR, "genetic_diversity_summary_table.csv")
main_xlsx <- file.path(TABLES_DIR, "genetic_diversity_summary_table.xlsx")
main_docx <- file.path(TABLES_DIR, "genetic_diversity_summary_table.docx")

table_title <- "Table X. Genetic diversity indices for Fagus grandifolia sites ordered from south to north."
table_note <- "N = number of individuals retained after filtering; Na = mean number of alleles per locus; Ar = allelic richness; Ho = observed heterozygosity; He = expected heterozygosity; FIS = inbreeding coefficient."

write.csv(genetic_diversity_summary_table, main_csv, row.names = FALSE, na = "")
message(SCRIPT_TAG, " Saved CSV table: ", main_csv)

xlsx_path <- write_genetic_diversity_xlsx(genetic_diversity_summary_table, main_xlsx)
if (!is.na(xlsx_path)) message(SCRIPT_TAG, " Saved Excel table: ", xlsx_path)

write_genetic_diversity_docx(genetic_diversity_summary_table_word, main_docx, title = table_title, note = table_note)

# -----------------------------
# Final console diagnostics
# -----------------------------
message("\n", SCRIPT_TAG, " Final diagnostic table:")
print(genetic_diversity_summary_table %>% select(Site, N), row.names = FALSE)

message("\n", SCRIPT_TAG, " Genetic diversity summary table:")
print(genetic_diversity_summary_table_word, row.names = FALSE)

message("\n", SCRIPT_TAG, " Files saved:")
message(" - ", main_csv)
if (!is.na(xlsx_path)) message(" - ", xlsx_path)
message(" - ", main_docx)

message("\n", SCRIPT_TAG, " Dataset used:")
message(" - ", GI_MLL_FILE)
message(" - Object: gi_mll")
message(" - Individuals retained after filtering: ", adegenet::nInd(gen_obj))
message(" - Loci retained after filtering: ", adegenet::nLoc(gen_obj))
message(" - Rarefaction min.n for Ar: ", min_n)
message("\n", SCRIPT_TAG, " Done.")