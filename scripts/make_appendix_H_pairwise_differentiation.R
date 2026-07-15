############################################################
# Appendix H Table H.1: Pairwise genetic differentiation
# among 12 American beech (Fagus grandifolia Ehrh.) sites.
#
# This script recomputes site-level pairwise Jost's D and
# Weir & Cockerham's FST from the final clone-corrected
# microsatellite dataset used by the thesis workflow.
#
# IMPORTANT:
# - Do not use outputs/tables/pairwise_jostD_long.csv for Appendix H.
#   That file can contain repeated within-site comparisons and is not
#   guaranteed to be a complete site-level 12 x 12 table.
# - The canonical clone-corrected object is gi_mll, loaded via
#   scripts/_load_objects.R from outputs/v1/objects/gi_mll.rds.
############################################################

suppressWarnings(suppressPackageStartupMessages({
  required_packages <- c(
    "adegenet", "dplyr", "flextable", "hierfstat", "mmod",
    "officer", "poppr", "readr", "tibble"
  )
  missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_packages) > 0) {
    stop(
      "[appendix_H] Missing required package(s): ", paste(missing_packages, collapse = ", "),
      ". Install them before running this script.",
      call. = FALSE
    )
  }
  
  library(adegenet)
  library(dplyr)
  library(flextable)
  library(hierfstat)
  library(mmod)
  library(officer)
  library(poppr)
  library(readr)
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

OUTPUTS_DIR <- file.path(PROJECT_ROOT, "outputs")
TABLES_DIR <- file.path(OUTPUTS_DIR, "tables")
SUPPLEMENTARY_DIR <- file.path(TABLES_DIR, "supplementary")
OBJECTS_DIR <- file.path(OUTPUTS_DIR, "v1", "objects")
LOADER_SCRIPT <- file.path(PROJECT_ROOT, "scripts", "_load_objects.R")

OUTPUT_CSV <- file.path(TABLES_DIR, "appendix_H_pairwise_differentiation.csv")
OUTPUT_DOCX <- file.path(TABLES_DIR, "appendix_H_pairwise_differentiation.docx")
JOST_MATRIX_CSV <- file.path(SUPPLEMENTARY_DIR, "jost_d_matrix_ordered.csv")
FST_MATRIX_CSV <- file.path(SUPPLEMENTARY_DIR, "fst_matrix_ordered.csv")

EXPECTED_SITES <- c("AMC", "ALB", "IKJ", "IKO", "LGG", "LGR", "ML1", "ML2", "ML3", "CPF", "PLI", "LDF")
TOLERANCE <- 1e-6
TITLE_TEXT <- "Table H.1. Pairwise genetic differentiation among the 12 American beech (Fagus grandifolia Ehrh.) study sites based on clone-corrected microsatellite data."
NOTE_TEXT <- "Note. Pairwise Jost's D values are presented below the diagonal and Weir and Cockerham's FST values above the diagonal. Sites are ordered as AMC, ALB, IKJ, IKO, LGG, LGR, ML1, ML2, ML3, CPF, PLI, and LDF. Diagonal cells represent within-site comparisons and are therefore omitted."

message_h <- function(...) message("[appendix_H] ", ...)

for (d in c(TABLES_DIR, SUPPLEMENTARY_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ----------------------------
# Load final clone-corrected data
# ----------------------------
if (!file.exists(LOADER_SCRIPT)) {
  stop("[appendix_H] Missing loader script: ", LOADER_SCRIPT, call. = FALSE)
}

required_object_files <- c("gi_mll.rds", "df_ids.rds", "meta.rds")
missing_object_files <- required_object_files[!file.exists(file.path(OBJECTS_DIR, required_object_files))]
if (length(missing_object_files) > 0) {
  stop(
    "[appendix_H] Missing final thesis object file(s) in ", OBJECTS_DIR, ": ",
    paste(missing_object_files, collapse = ", "),
    ". Run scripts/00_master_pipeline.R first to regenerate the final clone-corrected objects.",
    call. = FALSE
  )
}

# _load_objects.R loads gi_mll and df_ids_mll and aligns pop(gi_mll) to Site.
source(LOADER_SCRIPT)

if (!exists("gi_mll")) {
  stop("[appendix_H] _load_objects.R did not create gi_mll.", call. = FALSE)
}
if (!inherits(gi_mll, "genind")) {
  stop("[appendix_H] gi_mll is not a genind/genclone object.", call. = FALSE)
}

message_h("Using final clone-corrected object: ", file.path(OBJECTS_DIR, "gi_mll.rds"))
message_h("nInd(gi_mll) = ", adegenet::nInd(gi_mll), "; nLoc(gi_mll) = ", adegenet::nLoc(gi_mll))

# ----------------------------
# Validation helpers
# ----------------------------
standardize_site <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- gsub("[^A-Z0-9]+", "", x)
  x
}

validate_sites <- function(gobj) {
  site_values <- standardize_site(adegenet::pop(gobj))
  if (length(site_values) != adegenet::nInd(gobj)) {
    stop("[appendix_H] pop(gi_mll) length does not match nInd(gi_mll).", call. = FALSE)
  }
  
  observed_sites <- unique(site_values)
  missing_sites <- setdiff(EXPECTED_SITES, observed_sites)
  extra_sites <- setdiff(observed_sites, EXPECTED_SITES)
  
  if (length(missing_sites) > 0 || length(extra_sites) > 0) {
    stop(
      "[appendix_H] Site labels in gi_mll do not match the 12 expected thesis sites.",
      "\nMissing expected site(s): ", ifelse(length(missing_sites) == 0, "none", paste(missing_sites, collapse = ", ")),
      "\nUnexpected site(s): ", ifelse(length(extra_sites) == 0, "none", paste(extra_sites, collapse = ", ")),
      "\nObserved site labels: ", paste(sort(observed_sites), collapse = ", "),
      call. = FALSE
    )
  }
  
  site_counts <- table(factor(site_values, levels = EXPECTED_SITES))
  if (any(site_counts == 0)) {
    stop("[appendix_H] At least one expected site has zero individuals after clone correction.", call. = FALSE)
  }
  
  adegenet::pop(gobj) <- factor(site_values, levels = EXPECTED_SITES)
  gobj
}

coerce_numeric_matrix <- function(x, label) {
  if (inherits(x, "dist")) {
    mat <- as.matrix(x)
  } else if (is.matrix(x)) {
    mat <- x
  } else if (is.data.frame(x)) {
    mat <- as.matrix(x)
  } else {
    mat <- tryCatch(as.matrix(x), error = function(e) NULL)
  }
  
  if (is.null(mat) || !is.matrix(mat)) {
    stop("[appendix_H] Could not coerce ", label, " result to a matrix.", call. = FALSE)
  }
  
  suppressWarnings(storage.mode(mat) <- "numeric")
  if (!is.numeric(mat)) {
    stop("[appendix_H] ", label, " matrix is not numeric after coercion.", call. = FALSE)
  }
  mat
}

order_matrix <- function(mat, label) {
  if (nrow(mat) != ncol(mat)) {
    stop("[appendix_H] ", label, " result is not a square matrix.", call. = FALSE)
  }
  
  rn <- standardize_site(rownames(mat))
  cn <- standardize_site(colnames(mat))
  
  if (length(rn) == nrow(mat) && length(cn) == ncol(mat) &&
      all(EXPECTED_SITES %in% rn) && all(EXPECTED_SITES %in% cn)) {
    rownames(mat) <- rn
    colnames(mat) <- cn
    mat <- mat[EXPECTED_SITES, EXPECTED_SITES, drop = FALSE]
  } else if (nrow(mat) == length(EXPECTED_SITES)) {
    # Some genetics functions return matrices in population-level order without
    # useful labels. Because pop(gi_mll) is explicitly factored with EXPECTED_SITES,
    # the function output order is expected to follow that population order.
    rownames(mat) <- EXPECTED_SITES
    colnames(mat) <- EXPECTED_SITES
  } else {
    stop(
      "[appendix_H] Could not align ", label, " matrix to the 12 expected sites.",
      "\nRow names: ", paste(rownames(mat), collapse = ", "),
      "\nColumn names: ", paste(colnames(mat), collapse = ", "),
      call. = FALSE
    )
  }
  
  mat
}

validate_pairwise_matrix <- function(mat, label) {
  if (!identical(rownames(mat), EXPECTED_SITES) || !identical(colnames(mat), EXPECTED_SITES)) {
    stop("[appendix_H] ", label, " matrix is not ordered with the expected site labels.", call. = FALSE)
  }
  if (!isTRUE(all.equal(mat, t(mat), tolerance = TOLERANCE, check.attributes = FALSE))) {
    stop("[appendix_H] ", label, " matrix is not symmetric within tolerance ", TOLERANCE, ".", call. = FALSE)
  }
  diag_values <- diag(mat)
  if (!all(is.na(diag_values) | abs(diag_values) <= TOLERANCE)) {
    stop("[appendix_H] ", label, " diagonal must be zero or NA before formatting.", call. = FALSE)
  }
  off_diag <- mat[row(mat) != col(mat)]
  if (any(is.na(off_diag))) {
    stop("[appendix_H] ", label, " matrix has missing off-diagonal pairwise values.", call. = FALSE)
  }
  invisible(TRUE)
}

format_matrix_for_csv <- function(mat) {
  as.data.frame(mat, check.names = FALSE) |>
    mutate(across(everything(), ~ ifelse(is.na(.x), NA_character_, sprintf("%.3f", round(.x, 3))))) |>
    tibble::rownames_to_column("Site")
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
    bg(part = "header", bg = "#F2F2F2") |>
    border_remove() |>
    hline_top(part = "header", border = fp_border(color = "#666666", width = 1)) |>
    hline_bottom(part = "header", border = fp_border(color = "#666666", width = 0.75)) |>
    hline_bottom(part = "body", border = fp_border(color = "#666666", width = 1)) |>
    padding(padding = 1.5, part = "all") |>
    width(j = 1, width = 0.55) |>
    width(j = 2:ncol(combined_tbl), width = 0.48) |>
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
# Recompute pairwise statistics from gi_mll
# ----------------------------
gi_mll <- validate_sites(gi_mll)
site_counts <- table(adegenet::pop(gi_mll))
message_h("Validated expected sites exactly: ", paste(EXPECTED_SITES, collapse = ", "))
message_h("Clone-corrected individuals per site: ", paste(paste(names(site_counts), site_counts, sep = "="), collapse = ", "))

message_h("Computing pairwise Jost's D with mmod::pairwise_D(linearized = FALSE).")
jost_raw <- mmod::pairwise_D(gi_mll, linearized = FALSE)
jost_mat <- coerce_numeric_matrix(jost_raw, "Jost's D") |>
  order_matrix("Jost's D")
diag(jost_mat) <- 0
validate_pairwise_matrix(jost_mat, "Jost's D")

message_h("Computing pairwise Weir and Cockerham's FST with hierfstat::pairwise.WCfst().")
hf <- hierfstat::genind2hierfstat(gi_mll)
if (!is.data.frame(hf) || ncol(hf) < 2) {
  stop("[appendix_H] hierfstat::genind2hierfstat returned an invalid data frame.", call. = FALSE)
}
if (!is.numeric(hf[[1]])) hf[[1]] <- as.integer(factor(hf[[1]], levels = EXPECTED_SITES))

fst_raw <- hierfstat::pairwise.WCfst(hf)
fst_mat <- coerce_numeric_matrix(fst_raw, "FST") |>
  order_matrix("FST")
diag(fst_mat) <- 0
validate_pairwise_matrix(fst_mat, "FST")

# ----------------------------
# Save outputs
# ----------------------------
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