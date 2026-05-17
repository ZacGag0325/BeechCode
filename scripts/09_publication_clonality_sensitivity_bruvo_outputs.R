# scripts/09_publication_clonality_sensitivity_bruvo_outputs.R
############################################################
# Publication-ready clonality, sensitivity, and Bruvo diagnostic outputs
#
# Produces:
# 1) Clonality repeated-stems percentage figure and source table
# 2) Sensitivity table comparing full 15-locus and reduced 11-locus datasets
# 3) Bruvo distance histogram/diagnostic figure and source table
#
# This script reuses the canonical BeechCode object loader and output-folder
# conventions in scripts/_load_objects.R. It does not change filtering,
# thresholds, site ordering, or final clone-assignment logic.
############################################################

suppressPackageStartupMessages({
  library(adegenet)
  library(poppr)
  library(hierfstat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readxl)
})

source("scripts/_load_objects.R")

SCRIPT_TAG <- "[09_publication_outputs]"
BRUVO_MLL_THRESHOLD <- 0.09
BRUVO_ALGORITHM <- "farthest_neighbor"
BRUVO_REPLEN_VALUE <- 2
SUSPECT_LOCI <- c("EJV8T_A_0", "ERHBI_A_0", "FCM5", "FG5")
EXPECTED_SITE_LABELS <- c("S1", "S2", "S3", "S4", "S5", "S6", "N1", "N2", "N3", "N4", "N5", "N6")
SITE_LOOKUP_SHEET <- "site_lookup"

CLONALITY_FIG_PNG <- file.path(FIGURES_DIR, "clonality_repeated_stems_percent_publication.png")
CLONALITY_FIG_PDF <- file.path(FIGURES_DIR, "clonality_repeated_stems_percent_publication.pdf")
CLONALITY_TABLE_CSV <- file.path(TABLES_DIR, "clonality_repeated_stems_percent_table.csv")
CLONALITY_TABLE_DOCX <- file.path(TABLES_DIR, "clonality_repeated_stems_percent_table.docx")

SENSITIVITY_TABLE_CSV <- file.path(TABLES_SUPP_DIR, "sensitivity_full_vs_reduced_summary.csv")
SENSITIVITY_TABLE_DOCX <- file.path(TABLES_SUPP_DIR, "sensitivity_full_vs_reduced_summary.docx")

BRUVO_FIG_PNG <- file.path(FIGURES_SUPP_DIR, "bruvo_distance_histogram_publication.png")
BRUVO_FIG_PDF <- file.path(FIGURES_SUPP_DIR, "bruvo_distance_histogram_publication.pdf")
BRUVO_SUMMARY_CSV <- file.path(TABLES_SUPP_DIR, "bruvo_distance_threshold_summary.csv")

# -----------------------------------------------------------------------------
# General helpers
# -----------------------------------------------------------------------------
normalize_lookup_key <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("\uFEFF", "", x, fixed = TRUE)
  x <- gsub("[[:cntrl:]]", "", x)
  x[is.na(x)] <- ""
  x
}

normalize_ascii_key <- function(x) {
  x <- iconv(as.character(x), from = "", to = "ASCII//TRANSLIT")
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x[is.na(x)] <- ""
  x
}

pick_column <- function(df, choices, label = "column", required = FALSE) {
  idx <- match(TRUE, normalize_ascii_key(names(df)) %in% normalize_ascii_key(choices), nomatch = 0)
  if (idx == 0) {
    if (isTRUE(required)) {
      stop(SCRIPT_TAG, " Could not find required ", label, ". Accepted names: ", paste(choices, collapse = ", "), call. = FALSE)
    }
    return(NA_character_)
  }
  names(df)[idx]
}

safe_mean <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (length(x) == 0 || all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

round3 <- function(x) round(as.numeric(x), 3)
round1 <- function(x) round(as.numeric(x), 1)
format3 <- function(x) ifelse(is.na(x), NA_character_, sprintf("%.3f", round3(x)))
format1 <- function(x) ifelse(is.na(x), NA_character_, sprintf("%.1f", round1(x)))

calc_richness <- function(n, g) {
  ifelse(is.na(n) | is.na(g) | n <= 1, NA_real_, (g - 1) / (n - 1))
}

validate_expected_site_order <- function(df, site_col = "Site", context = "site table") {
  if (!site_col %in% names(df)) {
    stop(SCRIPT_TAG, " ", context, " is missing column '", site_col, "'.", call. = FALSE)
  }
  observed <- as.character(df[[site_col]])
  if (!setequal(observed, EXPECTED_SITE_LABELS)) {
    stop(
      SCRIPT_TAG, " ", context, " must contain exactly S1-S6 and N1-N6.",
      "\nExpected: ", paste(EXPECTED_SITE_LABELS, collapse = ", "),
      "\nObserved: ", paste(sort(unique(observed)), collapse = ", "),
      call. = FALSE
    )
  }
  df %>%
    mutate("{site_col}" := factor(.data[[site_col]], levels = EXPECTED_SITE_LABELS, ordered = TRUE)) %>%
    arrange(.data[[site_col]]) %>%
    mutate("{site_col}" := as.character(.data[[site_col]]))
}

# -----------------------------------------------------------------------------
# Minimal Word export helper, matching the project strategy:
# officer/flextable when available, with a built-in docx fallback.
# -----------------------------------------------------------------------------
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
    "<w:tc><w:tcPr><w:tcW w:w=\"2400\" w:type=\"dxa\"/></w:tcPr>",
    "<w:p><w:r>", bold_xml, "<w:t>", xml_escape_word(value), "</w:t></w:r></w:p></w:tc>"
  )
}

word_paragraph_xml <- function(value, bold = FALSE) {
  bold_xml <- if (bold) "<w:rPr><w:b/></w:rPr>" else ""
  paste0("<w:p><w:r>", bold_xml, "<w:t>", xml_escape_word(value), "</w:t></w:r></w:p>")
}

write_basic_docx_table <- function(df, path, title, note = NULL) {
  tmp_dir <- tempfile("publication_outputs_docx_")
  dir.create(file.path(tmp_dir, "_rels"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(tmp_dir, "word", "_rels"), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  
  df <- as.data.frame(lapply(df, as.character), stringsAsFactors = FALSE, check.names = FALSE)
  rows <- apply(df, 1, function(row) paste0("<w:tr>", paste(vapply(row, word_cell_xml, character(1)), collapse = ""), "</w:tr>"))
  header <- paste0("<w:tr>", paste(vapply(names(df), word_cell_xml, character(1), bold = TRUE), collapse = ""), "</w:tr>")
  note_xml <- if (!is.null(note) && nzchar(note)) word_paragraph_xml(note) else ""
  
  document_xml <- paste0(
    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
    '<w:document xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main" ',
    'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><w:body>',
    word_paragraph_xml(title, bold = TRUE),
    '<w:tbl><w:tblPr><w:tblStyle w:val="TableGrid"/><w:tblW w:w="0" w:type="auto"/>',
    '<w:tblBorders><w:top w:val="single" w:sz="4" w:space="0" w:color="auto"/>',
    '<w:left w:val="single" w:sz="4" w:space="0" w:color="auto"/>',
    '<w:bottom w:val="single" w:sz="4" w:space="0" w:color="auto"/>',
    '<w:right w:val="single" w:sz="4" w:space="0" w:color="auto"/>',
    '<w:insideH w:val="single" w:sz="4" w:space="0" w:color="auto"/>',
    '<w:insideV w:val="single" w:sz="4" w:space="0" w:color="auto"/></w:tblBorders></w:tblPr>',
    header, paste(rows, collapse = ""), '</w:tbl>', note_xml,
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
  utils::zip(zipfile = path, files = list.files(tmp_dir, recursive = TRUE, all.files = TRUE, no.. = TRUE), flags = "-q")
  if (!file.exists(path)) stop(SCRIPT_TAG, " Failed to create Word document at: ", path, call. = FALSE)
  invisible(path)
}

add_flextable_to_docx <- function(doc, ft) {
  if (requireNamespace("flextable", quietly = TRUE) && "body_add_flextable" %in% getNamespaceExports("flextable")) {
    return(flextable::body_add_flextable(doc, value = ft))
  }
  if (requireNamespace("officer", quietly = TRUE) && "body_add_flextable" %in% getNamespaceExports("officer")) {
    return(officer::body_add_flextable(doc, value = ft))
  }
  stop("No compatible body_add_flextable() function is exported by officer/flextable.", call. = FALSE)
}

write_word_table <- function(df, path, title, note = NULL, left_align_cols = character()) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  status <- if (requireNamespace("officer", quietly = TRUE) && requireNamespace("flextable", quietly = TRUE)) {
    ft <- flextable::flextable(df) %>%
      flextable::theme_booktabs() %>%
      flextable::autofit() %>%
      flextable::align(align = "center", part = "all") %>%
      flextable::bold(part = "header") %>%
      flextable::fontsize(size = 10, part = "all") %>%
      flextable::fontsize(size = 11, part = "header")
    for (col in intersect(left_align_cols, names(df))) {
      ft <- flextable::align(ft, j = col, align = "left", part = "body")
    }
    doc <- officer::read_docx()
    doc <- officer::body_add_par(doc, title, style = "heading 1")
    doc <- add_flextable_to_docx(doc, ft)
    if (!is.null(note) && nzchar(note)) doc <- officer::body_add_par(doc, note, style = "Normal")
    print(doc, target = path)
    "created with officer/flextable"
  } else if (requireNamespace("officer", quietly = TRUE)) {
    doc <- officer::read_docx()
    doc <- officer::body_add_par(doc, title, style = "heading 1")
    doc <- officer::body_add_table(doc, value = as.data.frame(df), style = "table_template")
    if (!is.null(note) && nzchar(note)) doc <- officer::body_add_par(doc, note, style = "Normal")
    print(doc, target = path)
    "created with officer fallback"
  } else {
    write_basic_docx_table(df, path, title = title, note = note)
    "created with built-in docx fallback"
  }
  if (!file.exists(path)) stop(SCRIPT_TAG, " Word export was not created: ", path, call. = FALSE)
  message(SCRIPT_TAG, " Saved: ", path, " (", status, ")")
  invisible(path)
}

# -----------------------------------------------------------------------------
# Site lookup and label mapping
# -----------------------------------------------------------------------------
find_site_lookup_workbook <- function(raw_dir = file.path(PROJECT_ROOT, "data", "raw"), sheet = SITE_LOOKUP_SHEET) {
  if (!dir.exists(raw_dir)) return(NULL)
  excel_files <- list.files(raw_dir, pattern = "\\.(xlsx|xls)$", full.names = TRUE, ignore.case = TRUE)
  if (length(excel_files) == 0) return(NULL)
  has_lookup <- vapply(excel_files, function(path) {
    sheets <- tryCatch(readxl::excel_sheets(path), error = function(e) character(0))
    any(normalize_ascii_key(sheets) == normalize_ascii_key(sheet))
  }, logical(1))
  lookup_files <- excel_files[has_lookup]
  if (length(lookup_files) == 0) return(NULL)
  lookup_files[1]
}

load_site_lookup <- function() {
  workbook <- find_site_lookup_workbook()
  if (is.null(workbook)) {
    message(SCRIPT_TAG, " site_lookup workbook not found; using existing S1-S6/N1-N6 labels if already present in loaded objects.")
    return(NULL)
  }
  sheets <- readxl::excel_sheets(workbook)
  sheet <- sheets[match(normalize_ascii_key(SITE_LOOKUP_SHEET), normalize_ascii_key(sheets))]
  lookup <- suppressMessages(readxl::read_excel(workbook, sheet = sheet)) %>% as.data.frame(stringsAsFactors = FALSE)
  site_col <- pick_column(lookup, c("Site", "site", "site_code", "old_site", "code_site"), "site code", required = TRUE)
  label_col <- pick_column(lookup, c("Site_label", "site_label", "new_site", "site_new", "label", "site_id"), "site display label", required = FALSE)
  order_col <- pick_column(lookup, c("Site_order", "site_order", "order", "ordre", "sort", "south_north_order"), "site order", required = FALSE)
  if (is.na(label_col)) label_col <- site_col
  
  out <- lookup %>%
    transmute(
      old_site = normalize_lookup_key(.data[[site_col]]),
      Site = normalize_lookup_key(.data[[label_col]]),
      Site_order = if (!is.na(order_col)) suppressWarnings(as.numeric(.data[[order_col]])) else match(normalize_lookup_key(.data[[label_col]]), EXPECTED_SITE_LABELS)
    ) %>%
    filter(nzchar(old_site), nzchar(Site)) %>%
    distinct(old_site, .keep_all = TRUE)
  
  if (anyDuplicated(out$old_site)) stop(SCRIPT_TAG, " site_lookup contains duplicated site codes.", call. = FALSE)
  if (!all(EXPECTED_SITE_LABELS %in% out$Site)) {
    stop(SCRIPT_TAG, " site_lookup does not provide all expected labels S1-S6 and N1-N6.", call. = FALSE)
  }
  out
}

site_lookup <- load_site_lookup()

map_site_labels <- function(raw_site) {
  raw_site <- normalize_lookup_key(raw_site)
  if (all(raw_site %in% EXPECTED_SITE_LABELS)) return(raw_site)
  if (is.null(site_lookup)) {
    stop(
      SCRIPT_TAG, " Site labels are not already S1-S6/N1-N6 and no site_lookup workbook was found. ",
      "Run the project pipeline with data/raw/site_lookup available or provide mapped site labels in df_ids.",
      call. = FALSE
    )
  }
  mapped <- site_lookup$Site[match(raw_site, site_lookup$old_site)]
  if (any(is.na(mapped) | !nzchar(mapped))) {
    missing <- sort(unique(raw_site[is.na(mapped) | !nzchar(mapped)]))
    stop(SCRIPT_TAG, " Could not map these site codes with site_lookup: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  mapped
}

get_site_labels_for_genind <- function(gobj, context) {
  aligned <- align_df_ids_to_genind(gobj, df_ids, context = context)
  mapped <- map_site_labels(aligned$Site)
  if (!setequal(unique(mapped), EXPECTED_SITE_LABELS)) {
    stop(SCRIPT_TAG, " ", context, " did not resolve to exactly S1-S6 and N1-N6.", call. = FALSE)
  }
  mapped
}

set_expected_pop_labels <- function(gobj, context) {
  out <- gobj
  adegenet::pop(out) <- factor(get_site_labels_for_genind(out, context), levels = EXPECTED_SITE_LABELS)
  out
}

# -----------------------------------------------------------------------------
# Clonality assignment helpers
# -----------------------------------------------------------------------------
make_exact_mlg_labels <- function(gobj) {
  gc <- poppr::as.genclone(gobj)
  raw <- tryCatch(poppr::mlg.vector(gc), error = function(e) as.integer(factor(poppr::mlg(gc))))
  paste0("MLG_", as.integer(factor(raw)))
}

run_bruvo_call <- function(expr) {
  captured_warnings <- character(0)
  value <- withCallingHandlers(
    expr,
    warning = function(w) {
      msg <- conditionMessage(w)
      if (grepl("NAs introduced by coercion", msg, fixed = TRUE)) {
        captured_warnings <<- c(captured_warnings, msg)
        invokeRestart("muffleWarning")
      }
    }
  )
  list(value = value, warning_count = length(captured_warnings), warning_messages = unique(captured_warnings))
}

make_bruvo_mll_labels <- function(gobj, threshold = BRUVO_MLL_THRESHOLD, algorithm = BRUVO_ALGORITHM) {
  gc <- poppr::as.genclone(gobj)
  replen <- rep(BRUVO_REPLEN_VALUE, adegenet::nLoc(gobj))
  names(replen) <- adegenet::locNames(gobj)
  result <- run_bruvo_call({
    poppr::mlg.filter(gc, distance = poppr::bruvo.dist, replen = replen, algorithm = algorithm) <- threshold
  })
  raw <- poppr::mll(gc)
  list(
    labels = paste0("MLL_", as.integer(factor(raw))),
    warning_count = result$warning_count,
    warning_messages = result$warning_messages
  )
}

get_final_assignments <- function(gobj, dataset_label = "Full dataset", prefer_saved_assignments = TRUE) {
  gobj <- set_expected_pop_labels(gobj, paste0(SCRIPT_TAG, " ", dataset_label))
  aligned <- align_df_ids_to_genind(gobj, df_ids, context = paste0(SCRIPT_TAG, " ", dataset_label, " assignments"))
  has_saved <- isTRUE(prefer_saved_assignments) && all(c("MLG", "MLL") %in% names(aligned))
  
  if (has_saved) {
    mlg_labels <- as.character(aligned$MLG)
    mll_labels <- as.character(aligned$MLL)
    source <- "final saved MLG/MLL columns in df_ids"
  } else {
    mlg_labels <- make_exact_mlg_labels(gobj)
    mll_labels <- make_bruvo_mll_labels(gobj)$labels
    source <- "computed from genind with poppr::mlg.vector and Bruvo mlg.filter"
  }
  
  if (any(is.na(mlg_labels) | !nzchar(mlg_labels)) || any(is.na(mll_labels) | !nzchar(mll_labels))) {
    stop(SCRIPT_TAG, " Missing/blank MLG or MLL labels for ", dataset_label, ".", call. = FALSE)
  }
  
  data.frame(
    Individual = adegenet::indNames(gobj),
    Site = as.character(adegenet::pop(gobj)),
    MLG = mlg_labels,
    MLL = mll_labels,
    Assignment_source = source,
    stringsAsFactors = FALSE
  )
}

summarise_clonality_by_site <- function(assignments) {
  assignments %>%
    group_by(Site) %>%
    summarise(
      N = dplyr::n(),
      MLG = dplyr::n_distinct(MLG),
      MLL = dplyr::n_distinct(MLL),
      .groups = "drop"
    ) %>%
    mutate(
      R_MLL = calc_richness(N, MLL),
      Repeated_stems = N - MLL,
      Unique_percent = 100 * MLL / N,
      Repeated_percent = 100 * Repeated_stems / N
    ) %>%
    validate_expected_site_order("Site", "clonality repeated-stems table")
}

make_reduced_genind <- function(gobj, loci_to_remove) {
  loci_available <- adegenet::locNames(gobj)
  requested <- unique(trimws(as.character(loci_to_remove)))
  requested <- requested[nzchar(requested)]
  present <- intersect(requested, loci_available)
  missing <- setdiff(requested, loci_available)
  if (length(missing) > 0) {
    stop(SCRIPT_TAG, " Requested suspect loci are absent from the loaded dataset: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  kept <- setdiff(loci_available, present)
  if (length(kept) < 2) stop(SCRIPT_TAG, " Fewer than two loci would remain after suspect-locus removal.", call. = FALSE)
  
  loci_df <- adegenet::genind2df(gobj, sep = "/", usepop = FALSE)
  keep_df <- intersect(kept, names(loci_df))
  if (length(keep_df) != length(kept)) {
    missing_keep <- setdiff(kept, names(loci_df))
    stop(SCRIPT_TAG, " Failed to retain all expected loci while rebuilding reduced genind: ", paste(missing_keep, collapse = ", "), call. = FALSE)
  }
  
  reduced <- adegenet::df2genind(
    X = loci_df[, keep_df, drop = FALSE],
    sep = "/",
    ploidy = adegenet::ploidy(gobj),
    ind.names = adegenet::indNames(gobj),
    type = if (!is.null(gobj@type) && length(gobj@type) > 0) gobj@type else "codom",
    NA.char = "NA"
  )
  adegenet::pop(reduced) <- adegenet::pop(gobj)
  if (adegenet::nLoc(reduced) != length(kept)) {
    stop(SCRIPT_TAG, " Locus subsetting mismatch: expected ", length(kept), " retained loci, got ", adegenet::nLoc(reduced), ".", call. = FALSE)
  }
  list(reduced = reduced, kept = kept, removed = present)
}

compute_clonality_overall <- function(gobj, dataset_label, prefer_saved_assignments = FALSE) {
  assignments <- get_final_assignments(gobj, dataset_label = dataset_label, prefer_saved_assignments = prefer_saved_assignments)
  n <- nrow(assignments)
  mlg <- dplyr::n_distinct(assignments$MLG)
  mll <- dplyr::n_distinct(assignments$MLL)
  data.frame(
    Dataset = dataset_label,
    N = n,
    MLG = mlg,
    MLL = mll,
    R_MLL = calc_richness(n, mll),
    Repeated_stems = n - mll,
    stringsAsFactors = FALSE
  )
}

# -----------------------------------------------------------------------------
# Diversity helpers for full/reduced sensitivity table
# -----------------------------------------------------------------------------
extract_overall_stat <- function(overall_obj, stat_name) {
  if (is.null(overall_obj)) return(NA_real_)
  if (is.atomic(overall_obj) && !is.null(names(overall_obj))) {
    val <- overall_obj[[stat_name]]
    return(if (length(val) == 0) NA_real_ else as.numeric(val[1]))
  }
  if (is.matrix(overall_obj) || is.data.frame(overall_obj)) {
    rn <- rownames(overall_obj)
    cn <- colnames(overall_obj)
    if (!is.null(rn) && stat_name %in% rn) return(safe_mean(as.numeric(overall_obj[stat_name, , drop = TRUE])))
    if (!is.null(cn) && stat_name %in% cn) return(safe_mean(as.numeric(overall_obj[, stat_name, drop = TRUE])))
  }
  NA_real_
}

compute_diversity_overall <- function(gobj, dataset_label) {
  gobj <- set_expected_pop_labels(gobj, paste0(SCRIPT_TAG, " ", dataset_label, " diversity"))
  hf <- hierfstat::genind2hierfstat(gobj)
  bs <- hierfstat::basic.stats(hf)
  if (is.null(bs$Ho) || is.null(bs$Hs) || is.null(bs$Fis)) {
    stop(SCRIPT_TAG, " hierfstat::basic.stats did not return Ho, Hs, and Fis for ", dataset_label, ".", call. = FALSE)
  }
  ho <- extract_overall_stat(bs$overall, "Ho")
  he <- extract_overall_stat(bs$overall, "Hs")
  fis <- extract_overall_stat(bs$overall, "Fis")
  if (is.na(ho)) ho <- safe_mean(bs$Ho)
  if (is.na(he)) he <- safe_mean(bs$Hs)
  if (is.na(fis)) fis <- safe_mean(bs$Fis)
  data.frame(
    Dataset = dataset_label,
    N_diversity = adegenet::nInd(gobj),
    Loci_retained = adegenet::nLoc(gobj),
    Ho = ho,
    He = he,
    FIS = fis,
    stringsAsFactors = FALSE
  )
}

# -----------------------------------------------------------------------------
# 1) Publication-ready clonality repeated-stems percentage figure
# -----------------------------------------------------------------------------
message(SCRIPT_TAG, " Building clonality repeated-stems percentage figure from final MLL assignments...")
final_assignments <- get_final_assignments(gi, dataset_label = "Full final clonality dataset", prefer_saved_assignments = TRUE)
clonality_site_table <- summarise_clonality_by_site(final_assignments) %>%
  transmute(
    Site,
    N = as.integer(N),
    MLG = as.integer(MLG),
    MLL = as.integer(MLL),
    R_MLL = round3(R_MLL),
    Repeated_stems = as.integer(Repeated_stems),
    Repeated_percent = round1(Repeated_percent)
  )

write.csv(clonality_site_table, CLONALITY_TABLE_CSV, row.names = FALSE, na = "")
message(SCRIPT_TAG, " Saved: ", CLONALITY_TABLE_CSV)

clonality_word_table <- clonality_site_table %>%
  mutate(R_MLL = format3(R_MLL), Repeated_percent = format1(Repeated_percent))
write_word_table(
  clonality_word_table,
  CLONALITY_TABLE_DOCX,
  title = "Clonality repeated-stems percentage table",
  note = "Note. Repeated stems are calculated as N - MLL using final multilocus lineage assignments. R_MLL = (MLL - 1) / (N - 1).",
  left_align_cols = "Site"
)

clonality_plot_df <- clonality_site_table %>%
  mutate(
    Site = factor(Site, levels = EXPECTED_SITE_LABELS, ordered = TRUE),
    label = ifelse(Repeated_percent > 0, paste0(format1(Repeated_percent), "%"), NA_character_)
  )

clonality_plot <- ggplot(clonality_plot_df, aes(x = Site, y = Repeated_percent)) +
  geom_col(width = 0.72, fill = "grey35") +
  geom_text(aes(label = label), vjust = -0.35, size = 3.6, na.rm = TRUE) +
  scale_y_continuous(
    limits = c(0, max(15, max(clonality_plot_df$Repeated_percent, na.rm = TRUE) * 1.25)),
    expand = expansion(mult = c(0, 0.03))
  ) +
  labs(x = "Site", y = "Repeated stems (%)") +
  theme_classic(base_size = 12) +
  theme(
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11, color = "black"),
    axis.line = element_line(linewidth = 0.35, color = "black"),
    axis.ticks = element_line(linewidth = 0.35, color = "black"),
    panel.border = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

ggsave(CLONALITY_FIG_PNG, clonality_plot, width = 6.5, height = 4, dpi = 320, bg = "white")
ggsave(CLONALITY_FIG_PDF, clonality_plot, width = 6.5, height = 4, bg = "white")
message(SCRIPT_TAG, " Saved: ", CLONALITY_FIG_PNG)
message(SCRIPT_TAG, " Saved: ", CLONALITY_FIG_PDF)

# -----------------------------------------------------------------------------
# 2) Clean sensitivity-analysis summary table
# -----------------------------------------------------------------------------
message(SCRIPT_TAG, " Building full-vs-reduced sensitivity summary table...")
full_clonality <- compute_clonality_overall(gi, "Full dataset", prefer_saved_assignments = TRUE)
full_gi_for_reduction <- set_expected_pop_labels(gi, paste0(SCRIPT_TAG, " full clonality reduction"))
reduced_gi_info <- make_reduced_genind(full_gi_for_reduction, SUSPECT_LOCI)
reduced_clonality <- compute_clonality_overall(reduced_gi_info$reduced, "Reduced dataset", prefer_saved_assignments = FALSE)

full_diversity_gi <- set_expected_pop_labels(gi_mll, paste0(SCRIPT_TAG, " full clone-corrected diversity"))
reduced_diversity_info <- make_reduced_genind(full_diversity_gi, SUSPECT_LOCI)
full_diversity <- compute_diversity_overall(full_diversity_gi, "Full dataset")
reduced_diversity <- compute_diversity_overall(reduced_diversity_info$reduced, "Reduced dataset")

diversity_all <- bind_rows(full_diversity, reduced_diversity)
clonality_all <- bind_rows(full_clonality, reduced_clonality)

sensitivity_table <- clonality_all %>%
  left_join(diversity_all, by = "Dataset") %>%
  mutate(
    `Loci retained` = as.integer(Loci_retained),
    `Loci removed` = ifelse(Dataset == "Full dataset", "None", paste(SUSPECT_LOCI, collapse = ", ")),
    `Individuals/stems used for clonality` = as.integer(N),
    MLG = as.integer(MLG),
    MLL = as.integer(MLL),
    `Overall MLL richness` = round3(R_MLL),
    `Repeated MLL stems` = as.integer(Repeated_stems),
    Ho = round3(Ho),
    He = round3(He),
    FIS = round3(FIS)
  ) %>%
  select(
    Dataset,
    `Loci retained`,
    `Loci removed`,
    `Individuals/stems used for clonality`,
    MLG,
    MLL,
    `Overall MLL richness`,
    `Repeated MLL stems`,
    Ho,
    He,
    FIS
  )

write.csv(sensitivity_table, SENSITIVITY_TABLE_CSV, row.names = FALSE, na = "")
message(SCRIPT_TAG, " Saved: ", SENSITIVITY_TABLE_CSV)

sensitivity_word_table <- sensitivity_table %>%
  mutate(
    `Overall MLL richness` = format3(`Overall MLL richness`),
    Ho = format3(Ho),
    He = format3(He),
    FIS = format3(FIS)
  )
write_word_table(
  sensitivity_word_table,
  SENSITIVITY_TABLE_DOCX,
  title = "Table S6. Sensitivity analysis comparing clonality and genetic diversity estimates between the full 15-locus dataset and the reduced 11-locus dataset excluding four suspect loci.",
  note = "Note. The reduced dataset excluded EJV8T_A_0, ERHBI_A_0, FCM5, and FG5. MLG = multilocus genotype; MLL = multilocus lineage; Ho = observed heterozygosity; He = expected heterozygosity; FIS = inbreeding coefficient.",
  left_align_cols = c("Dataset", "Loci removed")
)

# -----------------------------------------------------------------------------
# 3) Publication-ready Bruvo distance histogram / diagnostic figure
# -----------------------------------------------------------------------------
message(SCRIPT_TAG, " Building Bruvo distance histogram and threshold diagnostic...")
bruvo_gi <- set_expected_pop_labels(gi, paste0(SCRIPT_TAG, " Bruvo diagnostic"))
bruvo_replen <- rep(BRUVO_REPLEN_VALUE, adegenet::nLoc(bruvo_gi))
names(bruvo_replen) <- adegenet::locNames(bruvo_gi)

if (any(is.na(names(bruvo_replen))) || any(!nzchar(names(bruvo_replen)))) {
  stop(SCRIPT_TAG, " Bruvo diagnostic requires named loci in gi.", call. = FALSE)
}

bruvo_dist_result <- run_bruvo_call(poppr::bruvo.dist(bruvo_gi, replen = bruvo_replen))
bruvo_values <- as.numeric(bruvo_dist_result$value)
bruvo_values <- bruvo_values[is.finite(bruvo_values)]
if (length(bruvo_values) == 0) stop(SCRIPT_TAG, " Bruvo distance calculation returned no finite pairwise distances.", call. = FALSE)

diagnostic_mll_result <- make_bruvo_mll_labels(bruvo_gi, threshold = BRUVO_MLL_THRESHOLD, algorithm = BRUVO_ALGORITHM)
final_mlg_total <- dplyr::n_distinct(final_assignments$MLG)
final_mll_total <- dplyr::n_distinct(final_assignments$MLL)
diagnostic_mll_total <- dplyr::n_distinct(diagnostic_mll_result$labels)
threshold_n <- sum(bruvo_values <= BRUVO_MLL_THRESHOLD, na.rm = TRUE)
threshold_pct <- 100 * threshold_n / length(bruvo_values)
bruvo_match <- identical(as.integer(final_mll_total), as.integer(diagnostic_mll_total))
bruvo_warning_count <- bruvo_dist_result$warning_count + diagnostic_mll_result$warning_count
bruvo_warning_messages <- unique(c(bruvo_dist_result$warning_messages, diagnostic_mll_result$warning_messages))

bruvo_summary <- data.frame(
  total_pairwise_comparisons = length(bruvo_values),
  number_pairwise_distances_le_threshold = threshold_n,
  percentage_pairwise_distances_le_threshold = threshold_pct,
  threshold_used = BRUVO_MLL_THRESHOLD,
  final_MLG_total = final_mlg_total,
  final_MLL_total = final_mll_total,
  diagnostic_MLL_total = diagnostic_mll_total,
  diagnostic_MLL_matches_final_MLL = bruvo_match,
  bruvo_algorithm = BRUVO_ALGORITHM,
  bruvo_replen = paste(paste(names(bruvo_replen), bruvo_replen, sep = "="), collapse = ";"),
  captured_bruvo_na_coercion_warning_count = bruvo_warning_count,
  captured_bruvo_na_coercion_warning_messages = ifelse(length(bruvo_warning_messages) == 0, NA_character_, paste(bruvo_warning_messages, collapse = " | ")),
  stringsAsFactors = FALSE
)
write.csv(bruvo_summary, BRUVO_SUMMARY_CSV, row.names = FALSE, na = "")
message(SCRIPT_TAG, " Saved: ", BRUVO_SUMMARY_CSV)

hist_df <- data.frame(Bruvo_distance = bruvo_values)
y_max <- max(ggplot_build(
  ggplot(hist_df, aes(x = Bruvo_distance)) + geom_histogram(binwidth = 0.01, boundary = 0, closed = "left")
)$data[[1]]$count, na.rm = TRUE)

bruvo_plot <- ggplot(hist_df, aes(x = Bruvo_distance)) +
  geom_histogram(binwidth = 0.01, boundary = 0, closed = "left", fill = "grey40", color = "white", linewidth = 0.15) +
  geom_vline(xintercept = BRUVO_MLL_THRESHOLD, linetype = "dashed", linewidth = 0.55, color = "black") +
  annotate(
    "text",
    x = BRUVO_MLL_THRESHOLD + 0.035,
    y = y_max * 0.92,
    label = paste0("Threshold = ", BRUVO_MLL_THRESHOLD),
    hjust = 0,
    vjust = 1,
    size = 3.6
  ) +
  scale_x_continuous(limits = c(0, max(1, max(bruvo_values, na.rm = TRUE))), expand = expansion(mult = c(0, 0.01))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Pairwise Bruvo distance", y = "Number of pairwise comparisons") +
  theme_classic(base_size = 12) +
  theme(
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11, color = "black"),
    axis.line = element_line(linewidth = 0.35, color = "black"),
    axis.ticks = element_line(linewidth = 0.35, color = "black"),
    panel.border = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

ggsave(BRUVO_FIG_PNG, bruvo_plot, width = 6.5, height = 4, dpi = 320, bg = "white")
ggsave(BRUVO_FIG_PDF, bruvo_plot, width = 6.5, height = 4, bg = "white")
message(SCRIPT_TAG, " Saved: ", BRUVO_FIG_PNG)
message(SCRIPT_TAG, " Saved: ", BRUVO_FIG_PDF)

# -----------------------------------------------------------------------------
# Console checks requested by user
# -----------------------------------------------------------------------------
clonality_total_stems <- nrow(final_assignments)
clonality_total_mlg <- dplyr::n_distinct(final_assignments$MLG)
clonality_total_mll <- dplyr::n_distinct(final_assignments$MLL)
clonality_repeated_stems <- clonality_total_stems - clonality_total_mll
clonality_repeated_pct <- 100 * clonality_repeated_stems / clonality_total_stems
sites_with_repeated <- clonality_site_table %>% filter(Repeated_stems > 0) %>% pull(Site)

cat("\n", SCRIPT_TAG, " Clonality console check\n", sep = "")
cat(SCRIPT_TAG, " total stems = ", clonality_total_stems, "\n", sep = "")
cat(SCRIPT_TAG, " total MLGs = ", clonality_total_mlg, "\n", sep = "")
cat(SCRIPT_TAG, " total MLLs = ", clonality_total_mll, "\n", sep = "")
cat(SCRIPT_TAG, " overall repeated stems = ", clonality_repeated_stems, "\n", sep = "")
cat(SCRIPT_TAG, " overall repeated percentage = ", sprintf("%.1f%%", clonality_repeated_pct), "\n", sep = "")
cat(SCRIPT_TAG, " sites with repeated MLLs = ", ifelse(length(sites_with_repeated) == 0, "none", paste(sites_with_repeated, collapse = ", ")), "\n", sep = "")

cat("\n", SCRIPT_TAG, " Sensitivity console check\n", sep = "")
cat(SCRIPT_TAG, " full vs reduced MLG count = ", full_clonality$MLG, " vs ", reduced_clonality$MLG, "\n", sep = "")
cat(SCRIPT_TAG, " full vs reduced MLL count = ", full_clonality$MLL, " vs ", reduced_clonality$MLL, "\n", sep = "")
cat(SCRIPT_TAG, " full Ho, He, FIS = ", sprintf("%.3f, %.3f, %.3f", full_diversity$Ho, full_diversity$He, full_diversity$FIS), "\n", sep = "")
cat(SCRIPT_TAG, " reduced Ho, He, FIS = ", sprintf("%.3f, %.3f, %.3f", reduced_diversity$Ho, reduced_diversity$He, reduced_diversity$FIS), "\n", sep = "")

cat("\n", SCRIPT_TAG, " Bruvo console check\n", sep = "")
cat(SCRIPT_TAG, " threshold used = ", BRUVO_MLL_THRESHOLD, "\n", sep = "")
cat(SCRIPT_TAG, " total pairwise comparisons = ", length(bruvo_values), "\n", sep = "")
cat(SCRIPT_TAG, " pairwise distances <= threshold = ", threshold_n, " (", sprintf("%.3f%%", threshold_pct), ")\n", sep = "")
cat(SCRIPT_TAG, " diagnostic MLL count matches final MLL count = ", bruvo_match, "\n", sep = "")
cat(SCRIPT_TAG, " All publication outputs complete.\n", sep = "")

invisible(list(
  clonality_site_table = clonality_site_table,
  sensitivity_table = sensitivity_table,
  bruvo_summary = bruvo_summary,
  clonality_plot = clonality_plot,
  bruvo_plot = bruvo_plot
))