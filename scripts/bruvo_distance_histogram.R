# scripts/bruvo_distance_histogram.R
############################################################
# Diagnostic verification of Bruvo-distance MLL clustering
#
# Purpose:
# - Verify whether the final clonality code defines MLLs using Bruvo distance
#   and a 0.09 threshold.
# - Recompute pairwise Bruvo distances from the same full genetic object used
#   for final clonality (gi) without modifying the main pipeline or any RDS
#   objects.
# - Recompute diagnostic Bruvo-based MLLs with the same settings used by the
#   final code and compare them with the current final MLL labels.
#
# Outputs:
# - outputs/figures/supplementary/bruvo_distance_histogram.png
# - outputs/figures/supplementary/bruvo_distance_histogram.pdf
# - outputs/tables/supplementary/bruvo_distance_summary.csv
############################################################

suppressPackageStartupMessages({
  library(adegenet)
  library(poppr)
  library(dplyr)
  library(ggplot2)
})

source("scripts/_load_objects.R")

SCRIPT_TAG <- "[bruvo_distance_histogram]"
ORIGINAL_CODE_FILES <- c(
  file.path(PROJECT_ROOT, "scripts", "00_master_pipeline.R"),
  file.path(PROJECT_ROOT, "scripts", "01_clonality.R")
)
BRUVO_MLL_THRESHOLD <- 0.09
BRUVO_ALGORITHM <- "farthest_neighbor"
OBJECT_USED <- "gi"

PNG_OUT <- file.path(FIGURES_SUPP_DIR, "bruvo_distance_histogram.png")
PDF_OUT <- file.path(FIGURES_SUPP_DIR, "bruvo_distance_histogram.pdf")
SUMMARY_OUT <- file.path(TABLES_SUPP_DIR, "bruvo_distance_summary.csv")

read_code_file <- function(path) {
  if (!file.exists(path)) {
    return(character(0))
  }
  readLines(path, warn = FALSE)
}

code_lines <- unlist(lapply(ORIGINAL_CODE_FILES, read_code_file), use.names = FALSE)
code_text <- paste(code_lines, collapse = "\n")

original_uses_bruvo <- grepl("bruvo\\.dist", code_text) && grepl("mlg\\.filter", code_text)
original_calls_explicit_bruvo_dist <- grepl("poppr::bruvo\\.dist|[^[:alnum:]_:]bruvo\\.dist", code_text)
original_threshold_found <- grepl("0\\.09", code_text) || grepl("BRUVO_MLL_THRESHOLD\\s*<-\\s*0\\.09", code_text)
original_uses_gi <- grepl("build_mll_clone_corrected_object\\s*<-\\s*function\\(gi\\)", code_text) ||
  grepl("compute_mlg_mll_from_gi\\s*<-\\s*function\\(gi", code_text)

if (!isTRUE(original_uses_bruvo)) {
  stop(
    SCRIPT_TAG,
    " The inspected final clonality code did not contain both mlg.filter and bruvo.dist. ",
    "Use a separate diagnostic script instead of treating Bruvo distance as the final method.",
    call. = FALSE
  )
}

if (!isTRUE(original_threshold_found)) {
  stop(
    SCRIPT_TAG,
    " The inspected final clonality code did not contain threshold 0.09. ",
    "Refusing to document a 0.09 Bruvo threshold without code evidence.",
    call. = FALSE
  )
}

if (!exists(OBJECT_USED, inherits = FALSE)) {
  stop(SCRIPT_TAG, " Required object 'gi' was not loaded by scripts/_load_objects.R.", call. = FALSE)
}

gobj <- gi

if (!inherits(gobj, "genind")) {
  stop(SCRIPT_TAG, " Object 'gi' is not a genind object.", call. = FALSE)
}

n_ind <- adegenet::nInd(gobj)
n_loci <- adegenet::nLoc(gobj)

if (n_ind < 2) {
  stop(SCRIPT_TAG, " Bruvo distance requires at least two individuals.", call. = FALSE)
}

if (n_loci < 2) {
  stop(SCRIPT_TAG, " Bruvo MLL clustering requires at least two loci.", call. = FALSE)
}

if (is.null(adegenet::locNames(gobj)) || any(!nzchar(adegenet::locNames(gobj)))) {
  stop(SCRIPT_TAG, " Bruvo repeat-length settings require named loci in gi.", call. = FALSE)
}

# The final code uses a constant repeat length of 2 for every retained locus.
# This intentionally mirrors scripts/00_master_pipeline.R and scripts/01_clonality.R.
bruvo_replen <- rep(2, n_loci)
names(bruvo_replen) <- adegenet::locNames(gobj)

if (length(bruvo_replen) != n_loci || any(is.na(bruvo_replen))) {
  stop(SCRIPT_TAG, " Bruvo repeat lengths are missing or incomplete.", call. = FALSE)
}

make_exact_mlg_labels <- function(x) {
  gc <- poppr::as.genclone(x)
  raw <- poppr::mlg.vector(gc)
  paste0("MLG_", as.integer(factor(raw)))
}

make_bruvo_mll_labels <- function(x, threshold, algorithm, replen) {
  gc_mll <- poppr::as.genclone(x)
  poppr::mlg.filter(
    gc_mll,
    distance = poppr::bruvo.dist,
    replen = replen,
    algorithm = algorithm
  ) <- threshold
  raw <- poppr::mll(gc_mll)
  paste0("MLL_", as.integer(factor(raw)))
}

extract_final_mll_labels <- function(df_ids_tbl, x) {
  if (!"MLL" %in% names(df_ids_tbl)) {
    return(rep(NA_character_, adegenet::nInd(x)))
  }
  cols <- resolve_df_ids_columns(df_ids_tbl, context = SCRIPT_TAG, require = TRUE)
  id_map <- setNames(as.character(df_ids_tbl$MLL), normalize_id(df_ids_tbl[[cols$id_col]]))
  unname(id_map[normalize_id(adegenet::indNames(x))])
}

extract_final_mlg_labels <- function(df_ids_tbl, x) {
  if (!"MLG" %in% names(df_ids_tbl)) {
    return(rep(NA_character_, adegenet::nInd(x)))
  }
  cols <- resolve_df_ids_columns(df_ids_tbl, context = SCRIPT_TAG, require = TRUE)
  id_map <- setNames(as.character(df_ids_tbl$MLG), normalize_id(df_ids_tbl[[cols$id_col]]))
  unname(id_map[normalize_id(adegenet::indNames(x))])
}

exact_mlg_labels <- make_exact_mlg_labels(gobj)
diagnostic_mll_labels <- make_bruvo_mll_labels(
  gobj,
  threshold = BRUVO_MLL_THRESHOLD,
  algorithm = BRUVO_ALGORITHM,
  replen = bruvo_replen
)
final_mlg_labels <- extract_final_mlg_labels(df_ids, gobj)
final_mll_labels <- extract_final_mll_labels(df_ids, gobj)

current_final_mlg_total <- if (all(is.na(final_mlg_labels))) NA_integer_ else dplyr::n_distinct(final_mlg_labels, na.rm = TRUE)
current_final_mll_total <- if (all(is.na(final_mll_labels))) {
  if (exists("gi_mll", inherits = FALSE)) adegenet::nInd(gi_mll) else NA_integer_
} else {
  dplyr::n_distinct(final_mll_labels, na.rm = TRUE)
}

exact_mlg_total <- dplyr::n_distinct(exact_mlg_labels, na.rm = TRUE)
diagnostic_bruvo_mll_total <- dplyr::n_distinct(diagnostic_mll_labels, na.rm = TRUE)
bruvo_matches_current_final_mll <- if (is.na(current_final_mll_total)) {
  NA
} else {
  identical(as.integer(diagnostic_bruvo_mll_total), as.integer(current_final_mll_total))
}

bruvo_dist <- poppr::bruvo.dist(gobj, replen = bruvo_replen)
bruvo_values <- as.numeric(bruvo_dist)
bruvo_values <- bruvo_values[is.finite(bruvo_values)]

if (length(bruvo_values) == 0) {
  stop(SCRIPT_TAG, " Bruvo distance calculation returned no finite pairwise distances.", call. = FALSE)
}

pairwise_le_threshold <- sum(bruvo_values <= BRUVO_MLL_THRESHOLD, na.rm = TRUE)
pairwise_gt_threshold <- sum(bruvo_values > BRUVO_MLL_THRESHOLD, na.rm = TRUE)
percent_le_threshold <- 100 * pairwise_le_threshold / length(bruvo_values)
percent_gt_threshold <- 100 * pairwise_gt_threshold / length(bruvo_values)

summary_tbl <- data.frame(
  object_used = OBJECT_USED,
  n_individuals = n_ind,
  n_loci = n_loci,
  exact_MLG_total = exact_mlg_total,
  current_final_MLG_total = current_final_mlg_total,
  current_final_MLL_total = current_final_mll_total,
  bruvo_based_MLL_total_at_threshold_0_09 = diagnostic_bruvo_mll_total,
  bruvo_threshold = BRUVO_MLL_THRESHOLD,
  bruvo_algorithm = BRUVO_ALGORITHM,
  bruvo_replen = paste(paste(names(bruvo_replen), bruvo_replen, sep = "="), collapse = ";"),
  n_pairwise_distances = length(bruvo_values),
  min_bruvo_distance = min(bruvo_values, na.rm = TRUE),
  first_quartile_bruvo_distance = as.numeric(stats::quantile(bruvo_values, 0.25, na.rm = TRUE, names = FALSE)),
  median_bruvo_distance = stats::median(bruvo_values, na.rm = TRUE),
  mean_bruvo_distance = mean(bruvo_values, na.rm = TRUE),
  third_quartile_bruvo_distance = as.numeric(stats::quantile(bruvo_values, 0.75, na.rm = TRUE, names = FALSE)),
  max_bruvo_distance = max(bruvo_values, na.rm = TRUE),
  n_pairwise_distances_le_0_09 = pairwise_le_threshold,
  percent_pairwise_distances_le_0_09 = percent_le_threshold,
  n_pairwise_distances_gt_0_09 = pairwise_gt_threshold,
  percent_pairwise_distances_gt_0_09 = percent_gt_threshold,
  original_final_code_uses_bruvo_distance = original_uses_bruvo,
  original_final_code_calls_poppr_bruvo_dist = original_calls_explicit_bruvo_dist,
  original_final_code_contains_threshold_0_09 = original_threshold_found,
  original_final_code_uses_gi = original_uses_gi,
  bruvo_based_MLL_result_matches_current_final_MLL_result = bruvo_matches_current_final_mll,
  stringsAsFactors = FALSE
)

hist_df <- data.frame(Bruvo_distance = bruvo_values)

p <- ggplot(hist_df, aes(x = Bruvo_distance)) +
  geom_histogram(bins = 60, color = "white", fill = "#2C7FB8") +
  geom_vline(xintercept = BRUVO_MLL_THRESHOLD, linetype = "dashed", linewidth = 0.8, color = "#D95F0E") +
  labs(
    title = "Pairwise Bruvo distance distribution",
    subtitle = paste0("Full final genetic object (gi); dashed line = MLL threshold ", BRUVO_MLL_THRESHOLD),
    x = "Bruvo distance",
    y = "Number of pairwise comparisons"
  ) +
  theme_bw(base_size = 12)

utils::write.csv(summary_tbl, SUMMARY_OUT, row.names = FALSE)
ggsave(PNG_OUT, p, width = 7.5, height = 5, dpi = 300)
ggsave(PDF_OUT, p, width = 7.5, height = 5)

cat("\n")
cat(SCRIPT_TAG, " Verification complete\n", sep = "")
cat(SCRIPT_TAG, " Object used: ", OBJECT_USED, "\n", sep = "")
cat(SCRIPT_TAG, " Number of individuals: ", n_ind, "\n", sep = "")
cat(SCRIPT_TAG, " Number of loci: ", n_loci, "\n", sep = "")
cat(SCRIPT_TAG, " Original final clonality code uses Bruvo distance: ", original_uses_bruvo, "\n", sep = "")
cat(SCRIPT_TAG, " Original final clonality code calls poppr::bruvo.dist/bruvo.dist: ", original_calls_explicit_bruvo_dist, "\n", sep = "")
cat(SCRIPT_TAG, " Original final clonality code contains threshold 0.09: ", original_threshold_found, "\n", sep = "")
cat(SCRIPT_TAG, " Repeat length information used: constant replen = 2 for each retained locus\n", sep = "")
cat(SCRIPT_TAG, " Current final MLG total: ", current_final_mlg_total, "\n", sep = "")
cat(SCRIPT_TAG, " Current final MLL total: ", current_final_mll_total, "\n", sep = "")
cat(SCRIPT_TAG, " Diagnostic exact MLG total: ", exact_mlg_total, "\n", sep = "")
cat(SCRIPT_TAG, " Diagnostic Bruvo MLL total at 0.09: ", diagnostic_bruvo_mll_total, "\n", sep = "")
cat(SCRIPT_TAG, " Diagnostic Bruvo MLL total matches current final MLL total: ", bruvo_matches_current_final_mll, "\n", sep = "")
cat(SCRIPT_TAG, " Output files:\n", sep = "")
cat("  - ", PNG_OUT, "\n", sep = "")
cat("  - ", PDF_OUT, "\n", sep = "")
cat("  - ", SUMMARY_OUT, "\n", sep = "")
cat("\n")

invisible(list(
  summary = summary_tbl,
  distances = bruvo_values,
  plot = p,
  diagnostic_MLL = diagnostic_mll_labels,
  exact_MLG = exact_mlg_labels
))