############################################################
# Appendix H Table H.1: Pairwise genetic differentiation
# among 12 American beech (Fagus grandifolia Ehrh.) sites.
#
# This script formats existing clone-corrected microsatellite
# pairwise differentiation outputs only. It does not recalculate
# Jost's D or Weir and Cockerham's FST.
############################################################

suppressWarnings(suppressPackageStartupMessages({
  required_packages <- c("dplyr", "flextable", "officer", "purrr", "readr", "stringr", "tibble")
  missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_packages) > 0) {
    stop(
      "[appendix_H] Missing required package(s): ", paste(missing_packages, collapse = ", "),
      ". Install them before running this formatting script.",
      call. = FALSE
    )
  }
  library(dplyr)
  library(flextable)
  library(officer)
  library(purrr)
  library(readr)
  library(stringr)
  library(tibble)
}))

# ----------------------------
# Project-root and path helpers
# ----------------------------
find_project_root <- function() {
  candidates <- c(getwd(), normalizePath(file.path(getwd(), ".."), mustWork = FALSE))
  cmd_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
  if (length(cmd_file) > 0 && nzchar(cmd_file[1])) {
    candidates <- c(candidates, dirname(normalizePath(cmd_file[1], mustWork = FALSE)))
  }
  
  for (start in unique(candidates)) {
    cur <- normalizePath(start, mustWork = FALSE)
    repeat {
      if (file.exists(file.path(cur, "scripts", "make_appendix_H_pairwise_differentiation.R"))) return(cur)
      parent <- dirname(cur)
      if (identical(parent, cur)) break
      cur <- parent
    }
  }
  
  stop(
    "[appendix_H] Cannot find project root containing scripts/make_appendix_H_pairwise_differentiation.R.",
    call. = FALSE
  )
}

PROJECT_ROOT <- find_project_root()
setwd(PROJECT_ROOT)

# Clear file-path variables are intentionally kept near the top.
OUTPUTS_DIR <- file.path(PROJECT_ROOT, "outputs")
TABLES_DIR <- file.path(OUTPUTS_DIR, "tables")
SUPPLEMENTARY_DIR <- file.path(TABLES_DIR, "supplementary")
POP_STRUCTURE_DIR <- file.path(OUTPUTS_DIR, "population_structure")
SEARCH_DIRS <- c(TABLES_DIR, SUPPLEMENTARY_DIR, POP_STRUCTURE_DIR, OUTPUTS_DIR)
OUTPUT_CSV <- file.path(TABLES_DIR, "appendix_H_pairwise_differentiation.csv")
OUTPUT_DOCX <- file.path(TABLES_DIR, "appendix_H_pairwise_differentiation.docx")
JOST_MATRIX_CSV <- file.path(SUPPLEMENTARY_DIR, "jost_d_matrix_ordered.csv")
FST_MATRIX_CSV <- file.path(SUPPLEMENTARY_DIR, "fst_matrix_ordered.csv")

EXPECTED_SITES <- c("S1", "S2", "S3", "S4", "S5", "S6", "N1", "N2", "N3", "N4", "N5", "N6")
TOLERANCE <- 1e-6
TITLE_TEXT <- "Table H.1. Pairwise genetic differentiation among the 12 American beech (Fagus grandifolia Ehrh.) study sites based on clone-corrected microsatellite data."
NOTE_TEXT <- "Note. Pairwise Jost's D values are presented below the diagonal and Weir and Cockerham's FST values above the diagonal. Sites are ordered from south to north. Diagonal cells represent within-site comparisons and are therefore omitted."

message_h <- function(...) message("[appendix_H] ", ...)

# ----------------------------
# Input discovery
# ----------------------------
list_candidate_files <- function() {
  dirs <- unique(normalizePath(SEARCH_DIRS[file.exists(SEARCH_DIRS)], mustWork = TRUE))
  if (length(dirs) == 0) {
    stop("[appendix_H] None of the expected output search directories exist: ", paste(SEARCH_DIRS, collapse = ", "), call. = FALSE)
  }
  files <- unlist(lapply(dirs, list.files, pattern = "\\.(csv|rds)$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE), use.names = FALSE)
  unique(normalizePath(files, mustWork = FALSE))
}

read_candidate <- function(path) {
  ext <- tolower(tools::file_ext(path))
  obj <- switch(
    ext,
    csv = suppressMessages(readr::read_csv(path, show_col_types = FALSE, progress = FALSE)),
    rds = readRDS(path),
    stop("Unsupported file extension: ", path, call. = FALSE)
  )
  if (is.matrix(obj)) obj <- as.data.frame(obj, check.names = FALSE)
  obj
}

# Score candidate files using filenames and column names so existing outputs are reused.
score_candidate <- function(path, statistic = c("jost", "fst")) {
  statistic <- match.arg(statistic)
  score <- 0
  haystack <- tolower(basename(path))
  path_text <- tolower(path)
  obj <- tryCatch(read_candidate(path), error = function(e) NULL)
  if (is.data.frame(obj)) haystack <- paste(haystack, tolower(paste(names(obj), collapse = " ")))
  
  score <- score + 2 * str_count(haystack, "pairwise|population differentiation|differentiation")
  # Down-rank sensitivity/comparison summaries: their filenames may contain
  # "differentiation", "Jost", or "FST", but they are not the pairwise 12 x 12
  # result tables needed for Appendix H.
  if (str_detect(path_text, "comparison|comparisons|full_vs|nosuspect|no_suspect|sensitivity")) {
    score <- score - 20
  }
  if (statistic == "jost") {
    score <- score + 5 * str_count(haystack, "jost|jost_d|dest|d_est")
    score <- score - 3 * str_count(haystack, "fst|weir|cockerham")
  } else {
    score <- score + 5 * str_count(haystack, "fst|f_st|weir|cockerham")
    score <- score - 3 * str_count(haystack, "jost|jost_d|dest|d_est")
  }
  score
}

select_input_file <- function(statistic) {
  candidates <- list_candidate_files()
  if (length(candidates) == 0) {
    stop("[appendix_H] No CSV or RDS files were found under: ", paste(SEARCH_DIRS, collapse = ", "), call. = FALSE)
  }
  scored <- tibble(path = candidates, score = vapply(candidates, score_candidate, numeric(1), statistic = statistic)) |>
    arrange(desc(score), path)
  
  # Do not select a file based on keywords alone. Try candidates in score order
  # and accept only the first file that can actually be parsed and validated as
  # a complete pairwise matrix for the requested statistic. This prevents
  # high-scoring comparison summaries from being misidentified as pairwise
  # Jost's D or FST inputs.
  attempted <- character(0)
  positive_candidates <- scored |> filter(score > 0)
  if (nrow(positive_candidates) > 0) {
    for (candidate_path in positive_candidates$path) {
      parsed <- tryCatch({
        mat <- read_pairwise_matrix(candidate_path, statistic)
        validate_pairwise_matrix(mat, toupper(statistic))
        TRUE
      }, error = function(e) {
        attempted <<- c(attempted, paste0(candidate_path, " -- ", conditionMessage(e)))
        FALSE
      })
      if (isTRUE(parsed)) return(candidate_path)
    }
  }
  
  if (length(attempted) == 0) {
    stop(
      "[appendix_H] Could not identify a suitable ", statistic, " file. Searched CSV/RDS files under: ",
      paste(SEARCH_DIRS, collapse = ", "),
      ". Expected filename or column terms include jost, jost_d, dest, pairwise, fst, weir, cockerham, or population differentiation.",
      call. = FALSE
    )
  }
  
  stop(
    "[appendix_H] Candidate files were found for ", statistic,
    ", but none could be parsed and validated as a complete 12-site pairwise matrix. Attempted:\n- ",
    paste(attempted, collapse = "\n- "),
    call. = FALSE
  )
}

# ----------------------------
# Site-label standardization and matrix conversion
# ----------------------------
standardize_site <- function(x) {
  raw <- toupper(trimws(as.character(x)))
  raw[is.na(raw)] <- ""
  
  # Prefer an explicit site token anywhere in the label (e.g., "Site S1",
  # "S1 clone-corrected", "pop_N03"). If no token is found, fall back to
  # cleaning the whole string below.
  token <- str_match(raw, "(^|[^A-Z0-9])([SN]0?[1-6])([^A-Z0-9]|$)")[, 3]
  cleaned <- str_replace_all(raw, "[^A-Z0-9]", "")
  cleaned <- str_replace(cleaned, "^SITE", "")
  cleaned <- str_replace(cleaned, "^POPULATION", "")
  cleaned <- str_replace(cleaned, "^POP", "")
  cleaned <- str_replace(cleaned, "^S0([1-6])$", "S\\1")
  cleaned <- str_replace(cleaned, "^N0([1-6])$", "N\\1")
  
  out <- ifelse(!is.na(token), token, cleaned)
  out <- str_replace(out, "^S0([1-6])$", "S\\1")
  out <- str_replace(out, "^N0([1-6])$", "N\\1")
  out
}

find_long_columns <- function(df, statistic) {
  nm <- names(df)
  clean <- tolower(str_replace_all(nm, "[^A-Za-z0-9]+", "_"))
  clean <- str_replace_all(clean, "^_|_$", "")
  compact <- str_replace_all(clean, "_", "")
  
  pop1_candidates <- c(
    "population1", "population_1", "pop1", "pop_1", "site1", "site_1",
    "site_a", "site_i", "pop_a", "pop_i", "population_a", "population_i",
    "from", "row", "var1", "x"
  )
  pop2_candidates <- c(
    "population2", "population_2", "pop2", "pop_2", "site2", "site_2",
    "site_b", "site_j", "pop_b", "pop_j", "population_b", "population_j",
    "to", "column", "col", "var2", "y"
  )
  
  pop1 <- which(clean %in% pop1_candidates | compact %in% pop1_candidates)[1]
  pop2 <- which(clean %in% pop2_candidates | compact %in% pop2_candidates)[1]
  
  # If column names are non-standard, identify the two columns whose values look
  # most like the expected site labels after standardization.
  if (is.na(pop1) || is.na(pop2)) {
    site_like_counts <- vapply(df, function(col) sum(standardize_site(col) %in% EXPECTED_SITES, na.rm = TRUE), integer(1))
    site_like <- which(site_like_counts > 0)
    if (length(site_like) >= 2) {
      ordered_site_like <- site_like[order(site_like_counts[site_like], decreasing = TRUE)]
      pop1 <- ordered_site_like[1]
      pop2 <- ordered_site_like[2]
    }
  }
  
  value_patterns <- if (statistic == "jost") {
    c("jost", "jostd", "jost_d", "josts_d", "josts", "dest", "d_est")
  } else {
    c("fst", "f_st", "wc_fst", "weir", "cockerham", "theta")
  }
  value <- which(map_lgl(seq_along(clean), function(i) {
    any(str_detect(clean[i], value_patterns)) ||
      any(str_detect(compact[i], value_patterns)) ||
      (statistic == "jost" && clean[i] %in% c("d", "value")) ||
      (statistic == "fst" && clean[i] %in% c("value"))
  }))[1]
  
  # If there is only one numeric value column after the two population columns,
  # use it as a final fallback. This covers simple long tables such as
  # population_1, population_2, value.
  if (is.na(value) && !any(is.na(c(pop1, pop2)))) {
    excluded <- c(pop1, pop2)
    numeric_candidates <- which(map_lgl(df, ~ !all(is.na(suppressWarnings(as.numeric(.x))))))
    numeric_candidates <- setdiff(numeric_candidates, excluded)
    if (length(numeric_candidates) == 1) value <- numeric_candidates[1]
  }
  
  if (any(is.na(c(pop1, pop2, value)))) return(NULL)
  list(pop1 = nm[pop1], pop2 = nm[pop2], value = nm[value])
}

matrix_from_square <- function(df) {
  rn <- rownames(df)
  # If the first column contains site labels, use it as row names.
  first_col <- df[[1]]
  first_sites <- standardize_site(first_col)
  if (sum(first_sites %in% EXPECTED_SITES, na.rm = TRUE) >= length(EXPECTED_SITES) / 2) {
    rn <- first_sites
    df <- df[-1]
  }
  cn <- standardize_site(names(df))
  if (length(rn) != nrow(df) || all(rn %in% as.character(seq_len(nrow(df))))) rn <- cn[seq_len(nrow(df))]
  rn <- standardize_site(rn)
  if (nrow(df) != ncol(df) || !all(EXPECTED_SITES %in% rn) || !all(EXPECTED_SITES %in% cn)) return(NULL)
  mat <- as.matrix(df)
  suppressWarnings(storage.mode(mat) <- "numeric")
  rownames(mat) <- rn
  colnames(mat) <- cn
  mat[EXPECTED_SITES, EXPECTED_SITES, drop = FALSE]
}

matrix_from_long <- function(df, statistic) {
  cols <- find_long_columns(df, statistic)
  if (is.null(cols)) return(NULL)
  long <- df |>
    transmute(pop1 = standardize_site(.data[[cols$pop1]]), pop2 = standardize_site(.data[[cols$pop2]]), value = suppressWarnings(as.numeric(.data[[cols$value]]))) |>
    filter(pop1 %in% EXPECTED_SITES, pop2 %in% EXPECTED_SITES)
  mat <- matrix(NA_real_, length(EXPECTED_SITES), length(EXPECTED_SITES), dimnames = list(EXPECTED_SITES, EXPECTED_SITES))
  if (nrow(long) == 0) return(NULL)
  for (i in seq_len(nrow(long))) {
    mat[long$pop1[i], long$pop2[i]] <- long$value[i]
    mat[long$pop2[i], long$pop1[i]] <- long$value[i]
  }
  diag(mat) <- ifelse(is.na(diag(mat)), 0, diag(mat))
  mat
}

read_pairwise_matrix <- function(path, statistic) {
  obj <- read_candidate(path)
  if (!is.data.frame(obj)) stop("[appendix_H] Input is not a data frame or matrix: ", path, call. = FALSE)
  mat <- matrix_from_square(obj)
  if (is.null(mat)) mat <- matrix_from_long(obj, statistic)
  if (is.null(mat)) {
    stop(
      "[appendix_H] Could not parse ", statistic,
      " input as either a square matrix or long pairwise table: ", path,
      ". Available columns: ", paste(names(obj), collapse = ", "),
      call. = FALSE
    )
  }
  mat[EXPECTED_SITES, EXPECTED_SITES, drop = FALSE]
}

validate_pairwise_matrix <- function(mat, label) {
  missing_sites <- setdiff(EXPECTED_SITES, rownames(mat))
  if (length(missing_sites) > 0) stop("[appendix_H] ", label, " matrix is missing site(s): ", paste(missing_sites, collapse = ", "), call. = FALSE)
  if (!isTRUE(all.equal(mat, t(mat), tolerance = TOLERANCE, check.attributes = FALSE))) {
    stop("[appendix_H] ", label, " matrix is not symmetric within tolerance ", TOLERANCE, ".", call. = FALSE)
  }
  diag_values <- diag(mat)
  if (!all(is.na(diag_values) | abs(diag_values) <= TOLERANCE)) {
    stop("[appendix_H] ", label, " diagonal must be zero or NA before formatting.", call. = FALSE)
  }
  off_diag <- mat[row(mat) != col(mat)]
  if (any(is.na(off_diag))) stop("[appendix_H] ", label, " matrix has missing off-diagonal pairwise values.", call. = FALSE)
  invisible(TRUE)
}

format_matrix_for_csv <- function(mat) {
  out <- as.data.frame(mat, check.names = FALSE) |>
    mutate(across(everything(), ~ ifelse(is.na(.x), NA_character_, sprintf("%.3f", round(.x, 3))))) |>
    tibble::rownames_to_column("Site")
  out
}

make_combined_table <- function(jost_mat, fst_mat) {
  combined <- matrix("—", length(EXPECTED_SITES), length(EXPECTED_SITES), dimnames = list(EXPECTED_SITES, EXPECTED_SITES))
  combined[lower.tri(combined)] <- sprintf("%.3f", round(jost_mat[lower.tri(jost_mat)], 3))
  combined[upper.tri(combined)] <- sprintf("%.3f", round(fst_mat[upper.tri(fst_mat)], 3))
  as.data.frame(combined, check.names = FALSE) |>
    tibble::rownames_to_column("Site") |>
    tibble::as_tibble()
}

write_appendix_docx <- function(combined_tbl, path) {
  title_par <- officer::fpar(
    officer::ftext("Table H.1. Pairwise genetic differentiation among the 12 American beech (", fp_text(font.size = 9, font.family = "Arial")),
    officer::ftext("Fagus grandifolia", fp_text(font.size = 9, italic = TRUE, font.family = "Arial")),
    officer::ftext(" Ehrh.) study sites based on clone-corrected microsatellite data.", fp_text(font.size = 9, font.family = "Arial"))
  )
  
  ft <- flextable(combined_tbl) |>
    theme_booktabs() |>
    fontsize(size = 8, part = "all") |>
    font(fontname = "Arial", part = "all") |>
    bold(part = "header") |>
    bold(j = 1, part = "body") |>
    align(align = "center", part = "all") |>
    align(j = 1, align = "center", part = "all") |>
    bg(part = "header", bg = "#F2F2F2") |>
    border_remove() |>
    hline_top(part = "header", border = fp_border(color = "#666666", width = 1)) |>
    hline_bottom(part = "header", border = fp_border(color = "#666666", width = 0.75)) |>
    hline_bottom(part = "body", border = fp_border(color = "#666666", width = 1)) |>
    padding(padding = 1.5, part = "all") |>
    width(j = 1, width = 0.45) |>
    width(j = 2:ncol(combined_tbl), width = 0.47) |>
    set_table_properties(layout = "fixed", width = 0.98)
  
  doc <- read_docx() |>
    body_set_default_section(prop_section(page_size = page_size(orient = "landscape"), page_margins = page_mar(top = 0.35, bottom = 0.35, left = 0.35, right = 0.35))) |>
    body_add_fpar(title_par, style = "Normal") |>
    body_add_par("", style = "Normal") |>
    body_add_flextable(ft, split = FALSE) |>
    body_add_par("", style = "Normal") |>
    body_add_par(NOTE_TEXT, style = "Normal")
  
  print(doc, target = path)
}

# ----------------------------
# Main workflow
# ----------------------------
dir.create(TABLES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SUPPLEMENTARY_DIR, recursive = TRUE, showWarnings = FALSE)

jost_file <- select_input_file("jost")
fst_file <- select_input_file("fst")
message_h("Jost's D input file: ", jost_file)
message_h("FST input file: ", fst_file)

jost_mat <- read_pairwise_matrix(jost_file, "jost")
fst_mat <- read_pairwise_matrix(fst_file, "fst")
validate_pairwise_matrix(jost_mat, "Jost's D")
validate_pairwise_matrix(fst_mat, "FST")

combined_tbl <- make_combined_table(jost_mat, fst_mat)
readr::write_csv(combined_tbl, OUTPUT_CSV, na = "")
readr::write_csv(format_matrix_for_csv(jost_mat), JOST_MATRIX_CSV, na = "")
readr::write_csv(format_matrix_for_csv(fst_mat), FST_MATRIX_CSV, na = "")
write_appendix_docx(combined_tbl, OUTPUT_DOCX)

jost_off <- jost_mat[row(jost_mat) != col(jost_mat)]
fst_off <- fst_mat[row(fst_mat) != col(fst_mat)]

message_h("Jost's D matrix dimensions: ", paste(dim(jost_mat), collapse = " x "))
message_h("FST matrix dimensions: ", paste(dim(fst_mat), collapse = " x "))
message_h("Site order: ", paste(EXPECTED_SITES, collapse = ", "))
message_h(sprintf("Jost's D off-diagonal summary: min = %.3f, max = %.3f, mean = %.3f", min(jost_off), max(jost_off), mean(jost_off)))
message_h(sprintf("FST off-diagonal summary: min = %.3f, max = %.3f, mean = %.3f", min(fst_off), max(fst_off), mean(fst_off)))
message_h("Saved combined CSV: ", OUTPUT_CSV)
message_h("Saved combined DOCX: ", OUTPUT_DOCX)
message_h("Saved ordered Jost's D matrix: ", JOST_MATRIX_CSV)
message_h("Saved ordered FST matrix: ", FST_MATRIX_CSV)

message_h("Final combined table:")
print(combined_tbl, n = nrow(combined_tbl), width = Inf)