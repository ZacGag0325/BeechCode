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
# - outputs/tables/genetic_diversity_summary_table.xlsx
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

# Reuse the project root, object-loading, metadata-alignment, validation, and
# output-directory strategy used by the main BeechCode genetic scripts.
source("scripts/_load_objects.R")

SCRIPT_TAG <- "[genetic_diversity_summary_table]"
EXPECTED_SITE_LABELS <- c("S1", "S2", "S3", "S4", "S5", "S6", "N1", "N2", "N3", "N4", "N5", "N6")
TABLE_TITLE <- "Table X. Genetic diversity indices for Fagus grandifolia sites ordered from south to north."
TABLE_NOTE <- "N = number of individuals retained after filtering; Na = mean number of alleles per locus; Ar = rarefied allelic richness; Ho = observed heterozygosity; He = expected heterozygosity; FIS = inbreeding coefficient."

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

contains_x_site_labels <- function(x) {
  any(grepl("^X[0-9]+$", as.character(x)))
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
  order_col <- pick_column(lookup, c("Site_order", "site_order", "order", "ordre", "sort", "south_north_order"), "south-to-north order")
  latitude_col <- pick_column(lookup, c("Latitude", "latitude", "lat", "LAT", "Lat", "y", "Y"), "latitude")
  
  if (is.na(label_col)) label_col <- old_site_col
  
  out <- lookup %>%
    mutate(
      old_site = normalize_lookup_key(.data[[old_site_col]]),
      site_label = normalize_lookup_key(.data[[label_col]]),
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
  if (anyDuplicated(out$site_label)) {
    stop(SCRIPT_TAG, " site_lookup contains duplicated display site labels after cleaning.", call. = FALSE)
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

build_site_order_table <- function(site_map, df_ids_tbl = NULL) {
  out <- site_map %>%
    transmute(original_site = Raw_site, Site = Site) %>%
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
    out <- out %>% left_join(lookup_order, by = c("original_site", "Site"))
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
    message(SCRIPT_TAG, " Ordering sites using expected S1-S6/N1-N6 order.")
    out <- out %>% mutate(expected_order = match(Site, EXPECTED_SITE_LABELS)) %>% arrange(is.na(expected_order), expected_order, natural_site_rank(Site), Site) %>% select(-expected_order)
  }
  
  out %>% mutate(display_order = row_number())
}

# -----------------------------
# Final cleaned genetic dataset and validated site assignment
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

# Always rebuild population/site labels from df_ids_mll by matching individual ID.
# This avoids using locus-like labels (X1, X2, ...) if any upstream object has an
# incorrect pop slot and explicitly confirms every retained individual has a site.
df_ids_cols_mll <- resolve_df_ids_columns(df_ids_mll, context = SCRIPT_TAG, require = TRUE)
metadata_ids <- normalize_id(df_ids_mll[[df_ids_cols_mll$id_col]])
metadata_sites <- normalize_lookup_key(df_ids_mll[[df_ids_cols_mll$site_col]])
if (anyDuplicated(metadata_ids)) {
  dup_ids <- unique(df_ids_mll[[df_ids_cols_mll$id_col]][duplicated(metadata_ids)])
  stop(SCRIPT_TAG, " df_ids_mll contains duplicated individual IDs: ", paste(head(dup_ids, 10), collapse = ", "), call. = FALSE)
}

match_idx <- match(normalize_id(adegenet::indNames(gi_mll)), metadata_ids)
if (any(is.na(match_idx))) {
  missing_ids <- adegenet::indNames(gi_mll)[is.na(match_idx)]
  stop(SCRIPT_TAG, " Could not match every gi_mll individual to df_ids_mll. Missing examples: ", paste(head(missing_ids, 10), collapse = ", "), call. = FALSE)
}

raw_sites <- metadata_sites[match_idx]
if (length(raw_sites) != adegenet::nInd(gi_mll) || any(is.na(raw_sites)) || any(!nzchar(raw_sites))) {
  bad_ids <- adegenet::indNames(gi_mll)[is.na(raw_sites) | !nzchar(raw_sites)]
  stop(SCRIPT_TAG, " Every gi_mll individual must have a valid site in df_ids_mll. Problem examples: ", paste(head(bad_ids, 10), collapse = ", "), call. = FALSE)
}
if (contains_x_site_labels(raw_sites)) {
  stop(SCRIPT_TAG, " Raw metadata site labels contain locus-like X labels; refusing to continue.", call. = FALSE)
}

site_map <- data.frame(
  Raw_site = raw_sites,
  Site = map_site_labels(raw_sites),
  stringsAsFactors = FALSE
)
if (contains_x_site_labels(site_map$Site)) {
  stop(SCRIPT_TAG, " Final site labels contain locus-like X labels after site_lookup mapping; refusing to continue.", call. = FALSE)
}

missing_lookup_sites <- if (!is.null(site_lookup)) setdiff(unique(raw_sites), site_lookup$old_site) else character(0)
if (length(missing_lookup_sites) > 0) {
  message(SCRIPT_TAG, " Site code(s) not found in site_lookup and retained as-is: ", paste(missing_lookup_sites, collapse = ", "))
}

final_site_labels <- sort(unique(site_map$Site))
if (!setequal(final_site_labels, EXPECTED_SITE_LABELS)) {
  stop(
    SCRIPT_TAG, " Final site labels do not match the expected 12 site labels.",
    "\nExpected: ", paste(EXPECTED_SITE_LABELS, collapse = ", "),
    "\nObserved: ", paste(final_site_labels, collapse = ", "),
    call. = FALSE
  )
}

# Work on a copy so sourced gi_mll remains intact for interactive sessions.
gen_obj <- gi_mll
adegenet::pop(gen_obj) <- factor(site_map$Site, levels = EXPECTED_SITE_LABELS)

if (any(is.na(adegenet::pop(gen_obj)))) {
  stop(SCRIPT_TAG, " Some individuals have NA pop labels after assigning validated site labels.", call. = FALSE)
}

site_order <- build_site_order_table(site_map, df_ids_tbl = df_ids_mll)
if (!setequal(site_order$Site, EXPECTED_SITE_LABELS) || nrow(site_order) != length(EXPECTED_SITE_LABELS)) {
  stop(SCRIPT_TAG, " Site-order table must contain exactly the 12 expected site labels.", call. = FALSE)
}

message(SCRIPT_TAG, " Individuals retained after filtering: ", adegenet::nInd(gen_obj))
message(SCRIPT_TAG, " Loci retained after filtering: ", adegenet::nLoc(gen_obj))
message(SCRIPT_TAG, " Final site labels: ", paste(site_order$Site, collapse = ", "))

# -----------------------------
# Diversity calculations
# -----------------------------
pop_sizes <- table(as.character(adegenet::pop(gen_obj)))
pop_sizes <- pop_sizes[site_order$Site]
if (any(is.na(pop_sizes)) || any(pop_sizes <= 0)) {
  stop(SCRIPT_TAG, " N per site contains missing or zero values after population assignment.", call. = FALSE)
}
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

# hierfstat::basic.stats stores loci as rows and populations as columns in this
# project. The previous version accidentally treated rows as sites, producing
# X1, X2, ..., X15. These helpers explicitly aggregate by population/site.
extract_basic_stat_by_site <- function(stat_matrix, stat_name, site_labels) {
  mat <- as.matrix(stat_matrix)
  storage.mode(mat) <- "numeric"
  
  if (all(site_labels %in% colnames(mat))) {
    out <- vapply(site_labels, function(site) safe_mean(mat[, site]), numeric(1))
  } else if (all(site_labels %in% rownames(mat))) {
    out <- vapply(site_labels, function(site) safe_mean(mat[site, ]), numeric(1))
  } else {
    stop(
      SCRIPT_TAG, " Could not find expected site labels in hierfstat basic.stats matrix for ", stat_name,
      ". Row names: ", paste(head(rownames(mat), 20), collapse = ", "),
      ". Column names: ", paste(head(colnames(mat), 20), collapse = ", "),
      call. = FALSE
    )
  }
  
  data.frame(Site = names(out), value = as.numeric(out), stringsAsFactors = FALSE)
}

ho_by_site <- extract_basic_stat_by_site(bs$Ho, "Ho", site_order$Site) %>% rename(Ho = value)
he_by_site <- extract_basic_stat_by_site(bs$Hs, "He", site_order$Site) %>% rename(He = value)
fis_by_site <- extract_basic_stat_by_site(bs$Fis, "FIS", site_order$Site) %>% rename(FIS = value)

allele_tab <- adegenet::tab(gen_obj, NA.method = "asis")
loc_fac <- adegenet::locFac(gen_obj)
loci <- adegenet::locNames(gen_obj)
site_row_indices <- split(seq_len(nrow(allele_tab)), as.character(adegenet::pop(gen_obj)))
site_row_indices <- site_row_indices[site_order$Site]

count_alleles_for_site_locus <- function(rows, locus_name) {
  cols <- which(loc_fac == locus_name)
  if (length(cols) == 0 || length(rows) == 0) return(NA_integer_)
  mat <- allele_tab[rows, cols, drop = FALSE]
  if (all(is.na(mat))) return(NA_integer_)
  as.integer(sum(colSums(mat, na.rm = TRUE) > 0))
}

allele_counts_by_site_locus <- tidyr::expand_grid(
  Site = site_order$Site,
  Locus = loci
) %>%
  rowwise() %>%
  mutate(Allele_count = count_alleles_for_site_locus(site_row_indices[[Site]], Locus)) %>%
  ungroup()

mean_alleles_by_site <- allele_counts_by_site_locus %>%
  group_by(Site) %>%
  summarise(Na = safe_mean(Allele_count), .groups = "drop")

min_n <- min(as.integer(pop_sizes), na.rm = TRUE)
message(SCRIPT_TAG, " Rarefied allelic richness sample size (min.n): ", min_n)

ar_obj <- hierfstat::allelic.richness(hf, min.n = min_n)
if (is.null(ar_obj$Ar)) {
  stop(SCRIPT_TAG, " hierfstat::allelic.richness did not return an Ar matrix.", call. = FALSE)
}

extract_ar_by_site <- function(ar_matrix, site_labels) {
  mat <- as.matrix(ar_matrix)
  storage.mode(mat) <- "numeric"
  
  if (all(site_labels %in% colnames(mat))) {
    out <- vapply(site_labels, function(site) safe_mean(mat[, site]), numeric(1))
  } else if (all(site_labels %in% rownames(mat))) {
    out <- vapply(site_labels, function(site) safe_mean(mat[site, ]), numeric(1))
  } else {
    stop(
      SCRIPT_TAG, " Could not find expected site labels in hierfstat allelic.richness Ar matrix.",
      " Row names: ", paste(head(rownames(mat), 20), collapse = ", "),
      ". Column names: ", paste(head(colnames(mat), 20), collapse = ", "),
      call. = FALSE
    )
  }
  
  data.frame(Site = names(out), Ar = as.numeric(out), stringsAsFactors = FALSE)
}

allelic_richness_by_site <- extract_ar_by_site(ar_obj$Ar, site_order$Site)

n_by_site <- data.frame(
  Site = names(pop_sizes),
  N = as.integer(pop_sizes),
  stringsAsFactors = FALSE
)

genetic_diversity_summary_table <- n_by_site %>%
  left_join(mean_alleles_by_site, by = "Site") %>%
  left_join(allelic_richness_by_site, by = "Site") %>%
  left_join(ho_by_site, by = "Site") %>%
  left_join(he_by_site, by = "Site") %>%
  left_join(fis_by_site, by = "Site") %>%
  left_join(site_order %>% select(Site, display_order), by = "Site") %>%
  arrange(display_order) %>%
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
# Strong diagnostic checks requested for this standalone table
# -----------------------------
if (nrow(genetic_diversity_summary_table) != 12) {
  stop(SCRIPT_TAG, " Final table must have exactly 12 rows, but it has ", nrow(genetic_diversity_summary_table), ".", call. = FALSE)
}
if (!identical(genetic_diversity_summary_table$Site, site_order$Site)) {
  stop(SCRIPT_TAG, " Final table site order does not match the south-to-north site-order table.", call. = FALSE)
}
if (!setequal(genetic_diversity_summary_table$Site, EXPECTED_SITE_LABELS)) {
  stop(SCRIPT_TAG, " Final table does not contain exactly S1-S6 and N1-N6.", call. = FALSE)
}
if (contains_x_site_labels(genetic_diversity_summary_table$Site)) {
  stop(SCRIPT_TAG, " Final table contains locus-like X site labels (e.g., X1, X2), indicating a transposed/basic.stats orientation problem.", call. = FALSE)
}
if (any(is.na(genetic_diversity_summary_table$N))) {
  stop(SCRIPT_TAG, " Final table N contains NA values.", call. = FALSE)
}
if (any(is.na(genetic_diversity_summary_table$Na))) {
  stop(SCRIPT_TAG, " Final table Na contains NA values.", call. = FALSE)
}

# -----------------------------
# Export helpers
# -----------------------------
write_genetic_diversity_xlsx <- function(summary_tbl, path) {
  if (!has_openxlsx) {
    stop(SCRIPT_TAG, " openxlsx is required to save ", path, " but is not available in this R session.", call. = FALSE)
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
  invisible(path)
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
      if (has_officer && has_flextable && "body_add_flextable" %in% getNamespaceExports("flextable")) {
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
        doc <- flextable::body_add_flextable(doc, value = ft)
        doc <- officer::body_add_par(doc, note, style = "Normal")
        print(doc, target = path)
        "created with flextable::body_add_flextable"
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

write.csv(genetic_diversity_summary_table, main_csv, row.names = FALSE, na = "")
message(SCRIPT_TAG, " Saved CSV table: ", main_csv)

write_genetic_diversity_xlsx(genetic_diversity_summary_table, main_xlsx)
message(SCRIPT_TAG, " Saved Excel table: ", main_xlsx)

write_genetic_diversity_docx(genetic_diversity_summary_table_word, main_docx, title = TABLE_TITLE, note = TABLE_NOTE)

# -----------------------------
# Final console diagnostics
# -----------------------------
message("\n", SCRIPT_TAG, " Final N per site:")
print(genetic_diversity_summary_table %>% select(Site, N), row.names = FALSE)

message("\n", SCRIPT_TAG, " Final genetic diversity table:")
print(genetic_diversity_summary_table_word, row.names = FALSE)

message("\n", SCRIPT_TAG, " Saved output paths:")
message(" - ", main_csv)
message(" - ", main_xlsx)
message(" - ", main_docx)

message("\n", SCRIPT_TAG, " Dataset used:")
message(" - ", GI_MLL_FILE)
message(" - Object: gi_mll")
message(" - Individuals retained after filtering: ", adegenet::nInd(gen_obj))
message(" - Loci retained after filtering: ", adegenet::nLoc(gen_obj))
message(" - Rarefaction min.n for Ar: ", min_n)
message("\n", SCRIPT_TAG, " Done.")