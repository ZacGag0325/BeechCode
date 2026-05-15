# scripts/99_make_supplementary_tables.R
############################################################
# Final supplementary table builder for thesis/article Results
#
# This script gathers already-created analysis CSV outputs and writes
# clean supplementary CSV/XLSX/DOCX tables under:
#   outputs/tables/supplementary/
#
# It does not change primary analysis logic. If an input table is missing,
# the script reports which upstream script should be run first. The only
# exception is Table S1, where a fill-in primer/tail template is generated
# when no source primer table can be found.
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readxl)
  library(openxlsx)
})

# ------------------------------------------------------------------
# Project paths
# ------------------------------------------------------------------
find_project_root <- function() {
  candidates <- c(getwd(), normalizePath(file.path(getwd(), ".."), mustWork = FALSE))
  cmd_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
  if (length(cmd_file) > 0 && nzchar(cmd_file[1])) {
    candidates <- c(candidates, dirname(normalizePath(cmd_file[1], mustWork = FALSE)))
  }
  
  for (start in unique(candidates)) {
    cur <- normalizePath(start, mustWork = FALSE)
    repeat {
      if (file.exists(file.path(cur, "scripts", "99_make_supplementary_tables.R"))) return(cur)
      parent <- dirname(cur)
      if (identical(parent, cur)) break
      cur <- parent
    }
  }
  
  stop("Cannot find project root containing scripts/99_make_supplementary_tables.R", call. = FALSE)
}

PROJECT_ROOT <- find_project_root()
setwd(PROJECT_ROOT)

OUTPUT_DIR <- file.path(PROJECT_ROOT, "outputs")
TABLES_DIR <- file.path(OUTPUT_DIR, "tables")
FIGURES_DIR <- file.path(OUTPUT_DIR, "figures")
TABLES_SUPP_DIR <- file.path(TABLES_DIR, "supplementary")
FIGURES_SUPP_DIR <- file.path(FIGURES_DIR, "supplementary")
COMPARISON_DIR <- file.path(TABLES_DIR, "comparisons")
NO_SUSPECT_DIR <- file.path(TABLES_DIR, "noSuspectLoci")

for (d in c(OUTPUT_DIR, TABLES_DIR, TABLES_SUPP_DIR, FIGURES_DIR, FIGURES_SUPP_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

EXPECTED_SITE_LABELS <- c(paste0("S", 1:6), paste0("N", 1:6))
S1_COLUMNS <- c(
  "Locus", "Forward_primer", "Reverse_primer", "Tail_sequence", "Fluorescent_label",
  "Multiplex", "Repeat_motif", "Size_range", "Reference", "Included_in_final_dataset", "Notes"
)

created_outputs <- character(0)
missing_sources <- character(0)
scripts_to_run <- character(0)
manual_work <- character(0)

# ------------------------------------------------------------------
# General helpers
# ------------------------------------------------------------------
normalize_name <- function(x) {
  x <- iconv(as.character(x), from = "", to = "ASCII//TRANSLIT")
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  gsub("^_|_$", "", x)
}

normalize_key <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("\uFEFF", "", x, fixed = TRUE)
  x <- gsub("[[:cntrl:]]", "", x)
  x[is.na(x)] <- ""
  x
}

record_created <- function(path) {
  created_outputs <<- unique(c(created_outputs, normalizePath(path, mustWork = FALSE)))
  invisible(path)
}

record_missing <- function(path, script = NULL) {
  missing_sources <<- unique(c(missing_sources, path))
  if (!is.null(script)) scripts_to_run <<- unique(c(scripts_to_run, script))
  invisible(path)
}

read_csv_if_exists <- function(path, script = NULL) {
  if (!file.exists(path)) {
    record_missing(path, script)
    return(NULL)
  }
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

write_csv_out <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.csv(df, path, row.names = FALSE, na = "")
  message("[99_supplementary] Saved CSV: ", path)
  record_created(path)
}

write_xlsx_out <- function(sheets, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (is.data.frame(sheets)) sheets <- list(Table = sheets)
  sheets <- sheets[!vapply(sheets, is.null, logical(1))]
  if (length(sheets) == 0) return(invisible(NULL))
  openxlsx::write.xlsx(sheets, path, overwrite = TRUE)
  message("[99_supplementary] Saved XLSX: ", path)
  record_created(path)
}

round_for_docx <- function(df, digits = 3) {
  out <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
  for (nm in names(out)) {
    if (is.numeric(out[[nm]])) out[[nm]] <- ifelse(is.na(out[[nm]]), "", format(round(out[[nm]], digits), nsmall = 0, trim = TRUE, scientific = FALSE))
    if (is.logical(out[[nm]])) out[[nm]] <- ifelse(is.na(out[[nm]]), "", as.character(out[[nm]]))
  }
  out
}

xml_escape <- function(x) {
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
    "<w:tc><w:tcPr><w:tcW w:w=\"2200\" w:type=\"dxa\"/></w:tcPr><w:p><w:r>",
    bold_xml, "<w:t>", xml_escape(value), "</w:t></w:r></w:p></w:tc>"
  )
}

word_paragraph_xml <- function(value, bold = FALSE) {
  bold_xml <- if (bold) "<w:rPr><w:b/></w:rPr>" else ""
  paste0("<w:p><w:r>", bold_xml, "<w:t>", xml_escape(value), "</w:t></w:r></w:p>")
}

write_basic_docx <- function(tables, path, title, notes = NULL) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (is.data.frame(tables)) tables <- list(Table = tables)
  tables <- tables[!vapply(tables, is.null, logical(1))]
  if (length(tables) == 0) return(invisible(NULL))
  
  tmp_dir <- tempfile("supp_docx_")
  dir.create(file.path(tmp_dir, "_rels"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(tmp_dir, "word", "_rels"), recursive = TRUE, showWarnings = FALSE)
  
  table_xml <- vapply(names(tables), function(nm) {
    df <- round_for_docx(tables[[nm]])
    df <- as.data.frame(lapply(df, as.character), stringsAsFactors = FALSE, check.names = FALSE)
    header <- paste0("<w:tr>", paste(vapply(names(df), word_cell_xml, character(1), bold = TRUE), collapse = ""), "</w:tr>")
    rows <- apply(df, 1, function(row) paste0("<w:tr>", paste(vapply(row, word_cell_xml, character(1)), collapse = ""), "</w:tr>"))
    paste0(
      word_paragraph_xml(nm, bold = TRUE),
      "<w:tbl><w:tblPr><w:tblStyle w:val=\"TableGrid\"/><w:tblW w:w=\"0\" w:type=\"auto\"/><w:tblBorders>",
      "<w:top w:val=\"single\" w:sz=\"4\" w:space=\"0\" w:color=\"auto\"/><w:left w:val=\"single\" w:sz=\"4\" w:space=\"0\" w:color=\"auto\"/>",
      "<w:bottom w:val=\"single\" w:sz=\"4\" w:space=\"0\" w:color=\"auto\"/><w:right w:val=\"single\" w:sz=\"4\" w:space=\"0\" w:color=\"auto\"/>",
      "<w:insideH w:val=\"single\" w:sz=\"4\" w:space=\"0\" w:color=\"auto\"/><w:insideV w:val=\"single\" w:sz=\"4\" w:space=\"0\" w:color=\"auto\"/>",
      "</w:tblBorders></w:tblPr>", header, paste(rows, collapse = ""), "</w:tbl>"
    )
  }, character(1))
  
  note_xml <- if (!is.null(notes) && length(notes) > 0) paste(vapply(notes, word_paragraph_xml, character(1)), collapse = "") else ""
  document_xml <- paste0(
    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
    '<w:document xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><w:body>',
    word_paragraph_xml(title, bold = TRUE), note_xml, paste(table_xml, collapse = ""),
    '<w:sectPr><w:pgSz w:w="15840" w:h="12240" w:orient="landscape"/><w:pgMar w:top="720" w:right="720" w:bottom="720" w:left="720" w:header="360" w:footer="360" w:gutter="0"/></w:sectPr>',
    '</w:body></w:document>'
  )
  
  writeLines(c('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>', '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">', '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>', '<Default Extension="xml" ContentType="application/xml"/>', '<Override PartName="/word/document.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.document.main+xml"/>', '</Types>'), file.path(tmp_dir, "[Content_Types].xml"), useBytes = TRUE)
  writeLines(c('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>', '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">', '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="word/document.xml"/>', '</Relationships>'), file.path(tmp_dir, "_rels", ".rels"), useBytes = TRUE)
  writeLines(c('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>', '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"></Relationships>'), file.path(tmp_dir, "word", "_rels", "document.xml.rels"), useBytes = TRUE)
  writeLines(document_xml, file.path(tmp_dir, "word", "document.xml"), useBytes = TRUE)
  
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(tmp_dir)
  if (file.exists(path)) unlink(path)
  utils::zip(zipfile = path, files = list.files(tmp_dir, recursive = TRUE, all.files = TRUE, no.. = TRUE), flags = "-q")
  if (!file.exists(path)) stop("Failed to create DOCX: ", path, call. = FALSE)
  message("[99_supplementary] Saved DOCX: ", path)
  record_created(path)
}

order_sites <- function(df, site_col = "Site") {
  if (is.null(df) || !site_col %in% names(df)) return(df)
  site_values <- as.character(df[[site_col]])
  if (all(site_values %in% EXPECTED_SITE_LABELS)) {
    return(
      df %>%
        mutate("{site_col}" := factor(.data[[site_col]], levels = EXPECTED_SITE_LABELS, ordered = TRUE)) %>%
        arrange(.data[[site_col]]) %>%
        mutate("{site_col}" := as.character(.data[[site_col]]))
    )
  }
  df$.supp_site_order_tmp <- match(site_values, EXPECTED_SITE_LABELS)
  df$.supp_site_order_tmp[is.na(df$.supp_site_order_tmp)] <- Inf
  df <- df[order(df$.supp_site_order_tmp, seq_len(nrow(df))), , drop = FALSE]
  df$.supp_site_order_tmp <- NULL
  rownames(df) <- NULL
  df
}

load_site_lookup <- function() {
  raw_dir <- file.path(PROJECT_ROOT, "data", "raw")
  if (!dir.exists(raw_dir)) return(NULL)
  excel_files <- list.files(raw_dir, pattern = "\\.(xlsx|xls)$", full.names = TRUE, ignore.case = TRUE)
  for (path in excel_files) {
    sheets <- tryCatch(readxl::excel_sheets(path), error = function(e) character(0))
    hit <- which(normalize_name(sheets) == "site_lookup")
    if (length(hit) == 0) next
    lookup <- suppressMessages(readxl::read_excel(path, sheet = sheets[hit[1]])) %>% as.data.frame(stringsAsFactors = FALSE)
    old_col <- names(lookup)[match("site", normalize_name(names(lookup)), nomatch = 0)]
    label_col <- names(lookup)[match("site_label", normalize_name(names(lookup)), nomatch = 0)]
    order_col <- names(lookup)[match("site_order", normalize_name(names(lookup)), nomatch = 0)]
    if (length(old_col) == 0 || length(label_col) == 0) return(NULL)
    out <- data.frame(
      old_site = normalize_key(lookup[[old_col]]),
      site_label = normalize_key(lookup[[label_col]]),
      site_order = if (length(order_col) > 0) suppressWarnings(as.numeric(lookup[[order_col]])) else NA_real_,
      stringsAsFactors = FALSE
    ) %>% filter(nzchar(old_site), nzchar(site_label)) %>% distinct(old_site, .keep_all = TRUE)
    return(out)
  }
  NULL
}

site_lookup <- load_site_lookup()
convert_site_column <- function(df, site_col = "Site") {
  if (is.null(df) || !site_col %in% names(df) || is.null(site_lookup)) return(order_sites(df, site_col))
  idx <- match(normalize_key(df[[site_col]]), site_lookup$old_site)
  lab <- site_lookup$site_label[idx]
  df[[site_col]] <- ifelse(is.na(lab) | !nzchar(lab), as.character(df[[site_col]]), lab)
  order_sites(df, site_col)
}

select_existing_cols <- function(df, cols) {
  if (is.null(df)) return(NULL)
  keep <- cols[cols %in% names(df)]
  df[, keep, drop = FALSE]
}

# ------------------------------------------------------------------
# Table S1: microsatellite loci / primer / tail information
# ------------------------------------------------------------------
find_primer_source <- function() {
  search_dirs <- c(file.path(PROJECT_ROOT, "data"), file.path(PROJECT_ROOT, "sources"), PROJECT_ROOT)
  search_dirs <- search_dirs[dir.exists(search_dirs)]
  files <- unlist(lapply(search_dirs, function(d) {
    list.files(d, pattern = "\\.(xlsx|xls|csv|tsv)$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  }), use.names = FALSE)
  files <- files[!grepl("^outputs/", normalizePath(files, mustWork = FALSE))]
  if (length(files) == 0) return(NULL)
  
  score_file <- function(path) {
    nm <- normalize_name(basename(path))
    sum(vapply(c("primer", "microsatellite", "locus", "loci", "tail", "fluor", "multiplex", "appendix"), grepl, logical(1), x = nm))
  }
  candidates <- files[vapply(files, score_file, integer(1)) > 0]
  if (length(candidates) == 0) return(NULL)
  
  for (path in candidates[order(vapply(candidates, score_file, integer(1)), decreasing = TRUE)]) {
    ext <- tolower(tools::file_ext(path))
    if (ext %in% c("xlsx", "xls")) {
      sheets <- tryCatch(readxl::excel_sheets(path), error = function(e) character(0))
      for (sh in sheets) {
        dat <- tryCatch(suppressMessages(readxl::read_excel(path, sheet = sh, n_max = 5)) %>% as.data.frame(), error = function(e) NULL)
        if (is.null(dat)) next
        cn <- normalize_name(names(dat))
        if (any(grepl("primer|forward|reverse|tail|fluor|dye|multiplex|motif|reference", cn))) {
          return(list(path = path, sheet = sh, ext = ext))
        }
      }
    } else {
      dat <- tryCatch(read.csv(path, nrows = 5, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
      if (!is.null(dat) && any(grepl("primer|forward|reverse|tail|fluor|dye|multiplex|motif|reference", normalize_name(names(dat))))) {
        return(list(path = path, sheet = NA_character_, ext = ext))
      }
    }
  }
  NULL
}

standardise_s1 <- function(df) {
  nms <- normalize_name(names(df))
  pick <- function(patterns) {
    idx <- which(vapply(patterns, function(p) any(grepl(p, nms)), logical(1)))[1]
    if (is.na(idx)) return(NA_character_)
    hit <- which(grepl(patterns[idx], nms))[1]
    names(df)[hit]
  }
  mapping <- c(
    Locus = pick(c("^locus$", "loci", "marker")),
    Forward_primer = pick(c("forward.*primer", "primer.*forward", "^fwd", "^forward$")),
    Reverse_primer = pick(c("reverse.*primer", "primer.*reverse", "^rev", "^reverse$")),
    Tail_sequence = pick(c("tail.*sequence", "^tail$", "m13")),
    Fluorescent_label = pick(c("fluor", "dye", "label")),
    Multiplex = pick(c("multiplex", "plex")),
    Repeat_motif = pick(c("repeat.*motif", "motif")),
    Size_range = pick(c("size.*range", "range", "allele.*size")),
    Reference = pick(c("reference", "citation", "source")),
    Included_in_final_dataset = pick(c("included", "retained", "final")),
    Notes = pick(c("note", "comment"))
  )
  out <- setNames(data.frame(matrix(NA_character_, nrow = nrow(df), ncol = length(S1_COLUMNS))), S1_COLUMNS)
  for (col in names(mapping)) if (!is.na(mapping[[col]])) out[[col]] <- as.character(df[[mapping[[col]]]])
  out
}

primer_source <- find_primer_source()
if (is.null(primer_source)) {
  s1 <- setNames(data.frame(matrix(NA_character_, nrow = 1, ncol = length(S1_COLUMNS))), S1_COLUMNS)
  s1$Notes <- "TEMPLATE ONLY: no primer/tail source table was found automatically; fill values manually before final submission."
  warning("[99_supplementary] No microsatellite primer/tail source table found. Creating Table S1 template requiring manual completion.", call. = FALSE)
  manual_work <- unique(c(manual_work, "Table S1 primer/tail values must be filled manually unless a source primer table is added under data/ or sources/."))
} else {
  if (primer_source$ext %in% c("xlsx", "xls")) {
    s1_raw <- suppressMessages(readxl::read_excel(primer_source$path, sheet = primer_source$sheet)) %>% as.data.frame(stringsAsFactors = FALSE)
  } else if (primer_source$ext == "tsv") {
    s1_raw <- read.delim(primer_source$path, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    s1_raw <- read.csv(primer_source$path, stringsAsFactors = FALSE, check.names = FALSE)
  }
  s1 <- standardise_s1(s1_raw)
  s1$Notes <- ifelse(is.na(s1$Notes) | !nzchar(s1$Notes), paste0("Source: ", basename(primer_source$path), ifelse(is.na(primer_source$sheet), "", paste0(" / ", primer_source$sheet))), s1$Notes)
}
write_csv_out(s1, file.path(TABLES_SUPP_DIR, "Table_S1_microsatellite_loci.csv"))
write_xlsx_out(s1, file.path(TABLES_SUPP_DIR, "Table_S1_microsatellite_loci.xlsx"))
write_basic_docx(s1, file.path(TABLES_SUPP_DIR, "Table_S1_microsatellite_loci.docx"), "Table S1. Microsatellite loci, primers, fluorescent tails, and references")

# ------------------------------------------------------------------
# Table S2: HWE
# ------------------------------------------------------------------
hwe_site_locus <- read_csv_if_exists(file.path(TABLES_DIR, "hwe_by_site_by_locus.csv"), "scripts/02_hwe.R")
hwe_by_locus <- read_csv_if_exists(file.path(TABLES_DIR, "hwe_by_locus.csv"), "scripts/02_hwe.R")
hwe_by_site <- read_csv_if_exists(file.path(TABLES_DIR, "hwe_by_site.csv"), "scripts/02_hwe.R")
if (!is.null(hwe_site_locus)) {
  hwe_site_locus <- convert_site_column(hwe_site_locus, "Site")
  s2_cols <- c("Site", "Locus", "N", "missing_fraction", "n_unique_genotypes", "p_value_raw", "p_value_bonferroni", "p_value_fdr", "significant_raw", "significant_bonferroni", "significant_fdr", "status", "note")
  s2 <- select_existing_cols(hwe_site_locus, s2_cols)
  write_csv_out(s2, file.path(TABLES_SUPP_DIR, "Table_S2_HWE_by_site_by_locus.csv"))
  sheets <- list(
    by_site_by_locus = s2,
    by_locus = hwe_by_locus,
    by_site = convert_site_column(hwe_by_site, "Site")
  )
  write_xlsx_out(sheets, file.path(TABLES_SUPP_DIR, "Table_S2_HWE_summary.xlsx"))
  if (nrow(s2) <= 250) {
    write_basic_docx(s2, file.path(TABLES_SUPP_DIR, "Table_S2_HWE_by_site_by_locus.docx"), "Table S2. Hardy-Weinberg equilibrium tests by site and locus")
  } else {
    message("[99_supplementary] Skipped Table S2 DOCX because the site-by-locus table has ", nrow(s2), " rows; CSV/XLSX were written.")
  }
}

# ------------------------------------------------------------------
# Table S3: Linkage disequilibrium
# ------------------------------------------------------------------
ld_global <- read_csv_if_exists(file.path(TABLES_DIR, "linkage_disequilibrium_global.csv"), "scripts/02_hwe.R")
ld_by_site <- read_csv_if_exists(file.path(TABLES_DIR, "linkage_disequilibrium_by_site.csv"), "scripts/02_hwe.R")
ld_pairwise <- read_csv_if_exists(file.path(TABLES_DIR, "linkage_disequilibrium_pairwise.csv"), "scripts/02_hwe.R")
if (!is.null(ld_global)) {
  ld_global <- if ("Scope" %in% names(ld_global)) convert_site_column(ld_global, "Scope") else ld_global
  write_csv_out(ld_global, file.path(TABLES_SUPP_DIR, "Table_S3_linkage_disequilibrium_global.csv"))
}
if (!is.null(ld_by_site)) {
  ld_by_site <- convert_site_column(ld_by_site, "Site")
  write_csv_out(ld_by_site, file.path(TABLES_SUPP_DIR, "Table_S3_linkage_disequilibrium_by_site.csv"))
}
if (!is.null(ld_global) || !is.null(ld_by_site) || !is.null(ld_pairwise)) {
  write_xlsx_out(list(global = ld_global, by_site = ld_by_site, pairwise = ld_pairwise), file.path(TABLES_SUPP_DIR, "Table_S3_linkage_disequilibrium_summary.xlsx"))
}

# ------------------------------------------------------------------
# Table S4: Sensitivity analysis after removing suspect loci
# ------------------------------------------------------------------
locus_manifest <- read_csv_if_exists(file.path(NO_SUSPECT_DIR, "locus_manifest_noSuspectLoci.csv"), "scripts/hwe_sensitivity_analysis.R")
clonality_cmp <- read_csv_if_exists(file.path(TABLES_DIR, "hwe_sensitivity_clonality_comparison.csv"), "scripts/hwe_sensitivity_analysis.R")
clonality_overall <- read_csv_if_exists(file.path(TABLES_DIR, "hwe_sensitivity_clonality_overall_summary.csv"), "scripts/hwe_sensitivity_analysis.R")
het_overall <- read_csv_if_exists(file.path(COMPARISON_DIR, "full_vs_noSuspectLoci_heterozygosity_fis_overall_comparison.csv"), "scripts/hwe_sensitivity_analysis.R")
het_site <- read_csv_if_exists(file.path(COMPARISON_DIR, "full_vs_noSuspectLoci_heterozygosity_fis_by_site_comparison.csv"), "scripts/hwe_sensitivity_analysis.R")
ar_cmp <- read_csv_if_exists(file.path(COMPARISON_DIR, "full_vs_noSuspectLoci_allelic_richness_comparison.csv"), "scripts/hwe_sensitivity_analysis.R")
diff_cmp <- read_csv_if_exists(file.path(COMPARISON_DIR, "full_vs_noSuspectLoci_differentiation_comparison.csv"), "scripts/hwe_sensitivity_analysis.R")
if (!is.null(clonality_cmp)) write_csv_out(clonality_cmp, file.path(TABLES_SUPP_DIR, "Table_S4_sensitivity_clonality_comparison.csv"))
if (!is.null(het_site)) write_csv_out(convert_site_column(het_site, "Site"), file.path(TABLES_SUPP_DIR, "Table_S4_sensitivity_heterozygosity_fis_by_site.csv"))
if (!is.null(ar_cmp)) write_csv_out(convert_site_column(ar_cmp, "Site"), file.path(TABLES_SUPP_DIR, "Table_S4_sensitivity_allelic_richness.csv"))
if (!is.null(clonality_cmp) || !is.null(het_site) || !is.null(ar_cmp) || !is.null(diff_cmp)) {
  write_xlsx_out(
    list(
      locus_manifest = locus_manifest,
      clonality = clonality_cmp,
      clonality_overall = clonality_overall,
      heterozygosity_fis_overall = het_overall,
      heterozygosity_fis_by_site = convert_site_column(het_site, "Site"),
      allelic_richness = convert_site_column(ar_cmp, "Site"),
      differentiation = diff_cmp
    ),
    file.path(TABLES_SUPP_DIR, "Table_S4_sensitivity_analysis_summary.xlsx")
  )
  write_basic_docx(
    list(
      "S4a. Clonality sensitivity" = clonality_cmp,
      "S4b. Heterozygosity/FIS by site" = convert_site_column(het_site, "Site"),
      "S4c. Allelic richness" = convert_site_column(ar_cmp, "Site"),
      "S4d. Differentiation summary" = diff_cmp
    ),
    file.path(TABLES_SUPP_DIR, "Table_S4_sensitivity_analysis.docx"),
    "Table S4. Sensitivity analysis after removing suspect loci",
    notes = "FULL = all retained loci; REDUCED = EJV8T_A_0, ERHBI_A_0, FCM5, and FG5 removed. Clonality is calculated on gi before clone correction."
  )
}

# ------------------------------------------------------------------
# Table S5: Sampled-stem nearest-neighbour distances
# ------------------------------------------------------------------
dist_summary <- read_csv_if_exists(file.path(TABLES_DIR, "sampled_stem_distance_summary.csv"), "scripts/sampled_stem_descriptive_results.R")
dist_by_site <- read_csv_if_exists(file.path(TABLES_DIR, "sampled_stem_distance_by_site.csv"), "scripts/sampled_stem_descriptive_results.R")
nearest_pairs <- read_csv_if_exists(file.path(TABLES_DIR, "sampled_stem_nearest_neighbor_pairs.csv"), "scripts/sampled_stem_descriptive_results.R")
intended_pairs <- read_csv_if_exists(file.path(TABLES_DIR, "sampled_stem_intended_pair_distances.csv"), "scripts/sampled_stem_descriptive_results.R")
if (!is.null(dist_summary)) write_csv_out(dist_summary, file.path(TABLES_SUPP_DIR, "Table_S5_sampled_stem_distance_summary.csv"))
if (!is.null(dist_by_site)) write_csv_out(convert_site_column(dist_by_site, "Site"), file.path(TABLES_SUPP_DIR, "Table_S5_sampled_stem_distance_by_site.csv"))
if (!is.null(dist_summary) || !is.null(dist_by_site) || !is.null(nearest_pairs) || !is.null(intended_pairs)) {
  write_xlsx_out(list(summary = dist_summary, by_site = convert_site_column(dist_by_site, "Site"), nearest_neighbor_pairs = nearest_pairs, intended_pair_distances = intended_pairs), file.path(TABLES_SUPP_DIR, "Table_S5_sampled_stem_distances.xlsx"))
  write_basic_docx(list("S5a. Overall summary" = dist_summary, "S5b. By-site summary" = convert_site_column(dist_by_site, "Site")), file.path(TABLES_SUPP_DIR, "Table_S5_sampled_stem_distances.docx"), "Table S5. Sampled-stem nearest-neighbour distance summaries")
}

# ------------------------------------------------------------------
# Table S6: Sampled-stem DBH / developmental-stage summaries
# ------------------------------------------------------------------
dev_summary <- read_csv_if_exists(file.path(TABLES_DIR, "sampled_stem_developmental_stage_summary.csv"), "scripts/sampled_stem_descriptive_results.R")
dev_by_site <- read_csv_if_exists(file.path(TABLES_DIR, "sampled_stem_developmental_stage_by_site.csv"), "scripts/sampled_stem_descriptive_results.R")
dbh_distribution <- read_csv_if_exists(file.path(TABLES_DIR, "sampled_stem_dbh_distribution.csv"), "scripts/sampled_stem_descriptive_results.R or scripts/sampled_stem_dbh_distribution.R")
diameter_values <- read_csv_if_exists(file.path(TABLES_DIR, "sampled_stem_diameter_values.csv"), "scripts/sampled_stem_descriptive_results.R")
dbh_by_site <- read_csv_if_exists(file.path(TABLES_DIR, "sampled_stem_dbh_distribution_by_site.csv"), "scripts/sampled_stem_dbh_distribution.R")
if (is.null(dbh_by_site)) dbh_by_site <- dev_by_site
if (!is.null(dbh_distribution)) write_csv_out(dbh_distribution, file.path(TABLES_SUPP_DIR, "Table_S6_sampled_stem_dbh_summary.csv"))
if (!is.null(dbh_by_site)) write_csv_out(convert_site_column(dbh_by_site, "Site"), file.path(TABLES_SUPP_DIR, "Table_S6_sampled_stem_dbh_by_site.csv"))
if (!is.null(dev_summary) || !is.null(dev_by_site) || !is.null(dbh_distribution) || !is.null(diameter_values)) {
  write_xlsx_out(list(developmental_stage_summary = dev_summary, developmental_stage_by_site = convert_site_column(dev_by_site, "Site"), dbh_distribution = dbh_distribution, dbh_by_site = convert_site_column(dbh_by_site, "Site"), diameter_values = diameter_values), file.path(TABLES_SUPP_DIR, "Table_S6_sampled_stem_dbh_distribution.xlsx"))
  write_basic_docx(list("S6a. DBH distribution" = dbh_distribution, "S6b. DBH/developmental-stage by site" = convert_site_column(dbh_by_site, "Site"), "S6c. Developmental-stage summary" = dev_summary), file.path(TABLES_SUPP_DIR, "Table_S6_sampled_stem_dbh_distribution.docx"), "Table S6. Sampled-stem DBH and developmental-stage summaries")
}

# Explicit DBH figure label check from source code.
dbh_script <- file.path(PROJECT_ROOT, "scripts", "sampled_stem_dbh_distribution.R")
if (file.exists(dbh_script)) {
  dbh_code <- paste(readLines(dbh_script, warn = FALSE), collapse = "\n")
  if (grepl("Percent of stems", dbh_code, fixed = TRUE) && grepl("sprintf(\"%.1f\", `Percent of stems`)", dbh_code, fixed = TRUE)) {
    message("[99_supplementary] DBH figure label check: sampled_stem_dbh_distribution.R labels each bar using its row-specific `Percent of stems` value.")
  } else {
    warning("[99_supplementary] Could not confirm row-specific DBH percent labels in sampled_stem_dbh_distribution.R; inspect Figure 2 labels before final submission.", call. = FALSE)
    manual_work <- unique(c(manual_work, "Inspect DBH figure labels before final submission."))
  }
}

# ------------------------------------------------------------------
# Table S7: Bruvo-distance threshold diagnostic
# ------------------------------------------------------------------
bruvo_summary_path <- file.path(TABLES_SUPP_DIR, "bruvo_distance_summary.csv")
if (!file.exists(bruvo_summary_path) && file.exists(file.path(PROJECT_ROOT, "scripts", "bruvo_distance_histogram.R"))) {
  scripts_to_run <- unique(c(scripts_to_run, "scripts/bruvo_distance_histogram.R"))
}
bruvo_summary <- read_csv_if_exists(bruvo_summary_path, "scripts/bruvo_distance_histogram.R")
if (!is.null(bruvo_summary)) {
  write_csv_out(bruvo_summary, file.path(TABLES_SUPP_DIR, "Table_S7_bruvo_distance_diagnostic.csv"))
  write_basic_docx(bruvo_summary, file.path(TABLES_SUPP_DIR, "Table_S7_bruvo_distance_diagnostic.docx"), "Table S7. Bruvo-distance threshold diagnostic summary")
}

# ------------------------------------------------------------------
# Final reproducibility checklist
# ------------------------------------------------------------------
study_map_candidates <- c(
  file.path(FIGURES_DIR, "study_area_map.png"),
  file.path(FIGURES_DIR, "study_area_map.pdf"),
  file.path(FIGURES_DIR, "Figure_1_study_area_map.png"),
  file.path(FIGURES_DIR, "Figure_1_study_area_map.pdf")
)
if (!any(file.exists(study_map_candidates))) {
  manual_work <- unique(c(manual_work, "Study map / Figure 1 was not found under outputs/figures; create or copy the final map for the Results section."))
}

cat("\n============================================================\n")
cat("Supplementary table build checklist\n")
cat("============================================================\n")
expected_outputs <- c(
  "Table_S1_microsatellite_loci.csv", "Table_S1_microsatellite_loci.xlsx", "Table_S1_microsatellite_loci.docx",
  "Table_S2_HWE_by_site_by_locus.csv", "Table_S2_HWE_summary.xlsx", "Table_S2_HWE_by_site_by_locus.docx",
  "Table_S3_linkage_disequilibrium_summary.xlsx", "Table_S3_linkage_disequilibrium_global.csv", "Table_S3_linkage_disequilibrium_by_site.csv",
  "Table_S4_sensitivity_analysis_summary.xlsx", "Table_S4_sensitivity_clonality_comparison.csv", "Table_S4_sensitivity_heterozygosity_fis_by_site.csv", "Table_S4_sensitivity_allelic_richness.csv", "Table_S4_sensitivity_analysis.docx",
  "Table_S5_sampled_stem_distance_summary.csv", "Table_S5_sampled_stem_distance_by_site.csv", "Table_S5_sampled_stem_distances.xlsx", "Table_S5_sampled_stem_distances.docx",
  "Table_S6_sampled_stem_dbh_summary.csv", "Table_S6_sampled_stem_dbh_by_site.csv", "Table_S6_sampled_stem_dbh_distribution.xlsx", "Table_S6_sampled_stem_dbh_distribution.docx",
  "Table_S7_bruvo_distance_diagnostic.csv", "Table_S7_bruvo_distance_diagnostic.docx"
)
for (nm in expected_outputs) {
  path <- file.path(TABLES_SUPP_DIR, nm)
  cat(if (file.exists(path)) "[created] " else "[missing ] ", path, "\n", sep = "")
}

cat("\nMissing source files:\n")
if (length(missing_sources) == 0) cat("- None detected.\n") else cat(paste0("- ", sort(unique(missing_sources)), collapse = "\n"), "\n")

cat("\nScripts to run before scripts/99_make_supplementary_tables.R:\n")
if (length(scripts_to_run) == 0) cat("- None detected.\n") else cat(paste0("- ", sort(unique(scripts_to_run)), collapse = "\n"), "\n")

cat("\nOutputs that still need manual work:\n")
if (length(manual_work) == 0) cat("- None detected.\n") else cat(paste0("- ", sort(unique(manual_work)), collapse = "\n"), "\n")

cat("\nCreated output files in this run:\n")
if (length(created_outputs) == 0) cat("- None.\n") else cat(paste0("- ", sort(unique(created_outputs)), collapse = "\n"), "\n")
cat("============================================================\n")