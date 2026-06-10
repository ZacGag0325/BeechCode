# scripts/check_updated_main_poppr.R
############################################################
# QC check for the updated main BeechCode genotype workbook.
#
# Purpose:
# - Read data/raw/poppr.xlsx, sheet Genotypes-Indiv.
# - Calculate per-individual missingness across paired allele loci.
# - Apply the same retained/removed rule used by the main workflow:
#   remove individuals with >35% missing loci.
# - Report whether 278 individuals are retained without manually excluding
#   any redone LDF/N6 samples.
############################################################

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(janitor)
  library(stringr)
  library(readr)
})

INPUT_FILE <- file.path("data", "raw", "poppr.xlsx")
INPUT_SHEET <- "Genotypes-Indiv."
OUTPUT_DIR <- file.path("outputs", "qc")
OUTPUT_FILE <- file.path(OUTPUT_DIR, "updated_poppr_missingness_check.csv")
MISSINGNESS_THRESHOLD <- 0.35
EXPECTED_RETAINED_N <- 278L
EXPECTED_SAMPLE_ID_COL <- "Nom_Labo_Échantillons"

missing_tokens <- c("", "na", "n/a", "null", ".", "-", "?", "0", "-9")

clean_value <- function(x) {
  str_trim(as.character(x))
}

is_missing_allele <- function(x) {
  x_chr <- clean_value(x)
  is.na(x) | str_to_lower(x_chr) %in% missing_tokens
}

stop_if_missing_file_or_sheet <- function(path, sheet) {
  if (!file.exists(path)) {
    stop("Input file not found: ", path, call. = FALSE)
  }
  sheets <- readxl::excel_sheets(path)
  if (!(sheet %in% sheets)) {
    stop(
      "Sheet '", sheet, "' was not found in ", path, ". Available sheets: ",
      paste(sheets, collapse = ", "),
      call. = FALSE
    )
  }
}

find_column_exact_or_clean <- function(raw_names, cleaned_names, raw_target) {
  exact <- which(raw_names == raw_target)
  if (length(exact) > 0) return(raw_names[exact[1]])
  
  target_clean <- janitor::make_clean_names(raw_target)
  cleaned_match <- which(cleaned_names == target_clean)
  if (length(cleaned_match) > 0) return(raw_names[cleaned_match[1]])
  
  NA_character_
}

find_site_column <- function(raw_names, cleaned_names) {
  choices_raw <- c("Site", "site", "Population", "population", "Numéro_Population", "Numero_Population")
  raw_match <- which(str_to_lower(raw_names) %in% str_to_lower(choices_raw))
  if (length(raw_match) > 0) return(raw_names[raw_match[1]])
  
  choices_clean <- janitor::make_clean_names(choices_raw)
  clean_match <- which(cleaned_names %in% choices_clean)
  if (length(clean_match) > 0) return(raw_names[clean_match[1]])
  
  NA_character_
}

paired_locus_columns <- function(df) {
  nms <- names(df)
  allele_1 <- nms[str_detect(nms, "(_|\\.)1$")]
  allele_2 <- nms[str_detect(nms, "(_|\\.)2$")]
  loci_1 <- str_replace(allele_1, "(_|\\.)1$", "")
  loci_2 <- str_replace(allele_2, "(_|\\.)2$", "")
  loci <- sort(intersect(loci_1, loci_2))
  
  if (length(loci) == 0) {
    stop(
      "No paired allele columns ending in _1/_2 or .1/.2 were found. Available columns: ",
      paste(nms, collapse = ", "),
      call. = FALSE
    )
  }
  
  resolve_allele_col <- function(loc, suffix) {
    candidates <- c(paste0(loc, "_", suffix), paste0(loc, ".", suffix))
    hit <- candidates[candidates %in% nms]
    if (length(hit) == 0) return(NA_character_)
    hit[1]
  }
  
  out <- tibble(
    locus = loci,
    allele_1_col = vapply(loci, resolve_allele_col, character(1), suffix = "1"),
    allele_2_col = vapply(loci, resolve_allele_col, character(1), suffix = "2")
  )
  
  if (any(is.na(out$allele_1_col)) || any(is.na(out$allele_2_col))) {
    stop("Could not resolve all paired allele columns after locus detection.", call. = FALSE)
  }
  
  out
}

calculate_missingness <- function(df, sample_id_col, site_col, locus_pairs) {
  n_ind <- nrow(df)
  n_missing_loci <- integer(n_ind)
  
  for (i in seq_len(nrow(locus_pairs))) {
    m1 <- is_missing_allele(df[[locus_pairs$allele_1_col[i]]])
    m2 <- is_missing_allele(df[[locus_pairs$allele_2_col[i]]])
    # A locus is missing for this QC if both alleles are missing.
    n_missing_loci <- n_missing_loci + as.integer(m1 & m2)
  }
  
  n_loci_total <- nrow(locus_pairs)
  site_values <- if (!is.na(site_col)) clean_value(df[[site_col]]) else NA_character_
  
  tibble(
    sample_id = clean_value(df[[sample_id_col]]),
    site = site_values,
    n_loci_total = n_loci_total,
    n_loci_missing = n_missing_loci,
    missingness_prop = n_loci_missing / n_loci_total,
    missingness_percent = missingness_prop * 100,
    retained_at_35_percent = missingness_prop <= MISSINGNESS_THRESHOLD
  ) %>%
    filter(!is.na(sample_id), sample_id != "")
}

is_ldf_n6_sample <- function(sample_id, site) {
  site_chr <- str_to_upper(coalesce(as.character(site), ""))
  id_chr <- str_to_upper(coalesce(as.character(sample_id), ""))
  site_chr == "N6" | str_detect(site_chr, "LDF") | str_detect(id_chr, "LDF")
}

print_ldf_n6_results <- function(qc_tbl) {
  ldf_n6_tbl <- qc_tbl %>% filter(is_ldf_n6_sample(sample_id, site))
  
  cat("\n[check_updated_main_poppr] LDF/N6 sample missingness results:\n")
  if (nrow(ldf_n6_tbl) == 0) {
    cat("No LDF/N6 samples were identified from site labels or sample IDs.\n")
    return(invisible(ldf_n6_tbl))
  }
  
  cat("Could not automatically identify the two redone LDF/N6 samples by name; showing all LDF/N6 samples instead.\n")
  print(ldf_n6_tbl %>% arrange(site, sample_id), n = Inf)
  invisible(ldf_n6_tbl)
}

main <- function() {
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
  stop_if_missing_file_or_sheet(INPUT_FILE, INPUT_SHEET)
  
  raw <- readxl::read_excel(INPUT_FILE, sheet = INPUT_SHEET, .name_repair = "minimal") %>%
    as.data.frame(check.names = FALSE, stringsAsFactors = FALSE)
  
  raw_names <- names(raw)
  cleaned_names <- janitor::make_clean_names(raw_names)
  sample_id_col <- find_column_exact_or_clean(raw_names, cleaned_names, EXPECTED_SAMPLE_ID_COL)
  if (is.na(sample_id_col)) {
    stop(
      "Required sample ID column '", EXPECTED_SAMPLE_ID_COL, "' was not found in ",
      INPUT_FILE, " sheet '", INPUT_SHEET, "'. Available columns: ",
      paste(raw_names, collapse = ", "),
      call. = FALSE
    )
  }
  
  site_col <- find_site_column(raw_names, cleaned_names)
  if (is.na(site_col)) {
    warning("No site column was identified; site will be NA in the QC output.", call. = FALSE)
  }
  
  locus_pairs <- paired_locus_columns(raw)
  qc_tbl <- calculate_missingness(raw, sample_id_col, site_col, locus_pairs)
  
  readr::write_csv(qc_tbl, OUTPUT_FILE)
  ldf_n6_tbl <- print_ldf_n6_results(qc_tbl)
  
  retained_n <- sum(qc_tbl$retained_at_35_percent, na.rm = TRUE)
  removed_tbl <- qc_tbl %>% filter(!retained_at_35_percent) %>% arrange(desc(missingness_prop), site, sample_id)
  
  cat("\n[check_updated_main_poppr] Summary\n")
  cat("Total individuals in updated poppr.xlsx: ", nrow(qc_tbl), "\n", sep = "")
  cat("Number retained at 35% missingness threshold: ", retained_n, "\n", sep = "")
  cat("Number removed at 35% missingness threshold: ", nrow(removed_tbl), "\n", sep = "")
  cat("Retained number equals 278: ", ifelse(retained_n == EXPECTED_RETAINED_N, "yes", "no"), "\n", sep = "")
  
  if (retained_n != EXPECTED_RETAINED_N) {
    warning(
      "Updated poppr.xlsx did not produce 278 retained individuals at the 35% missingness threshold. Check outputs/qc/updated_poppr_missingness_check.csv.",
      call. = FALSE
    )
  }
  
  cat("\nAll LDF/N6 samples and whether they are retained:\n")
  if (nrow(ldf_n6_tbl) > 0) {
    print(ldf_n6_tbl %>% select(sample_id, site, missingness_percent, retained_at_35_percent), n = Inf)
  } else {
    cat("No LDF/N6 samples identified.\n")
  }
  
  cat("\nRemoved samples at >35% missingness:\n")
  if (nrow(removed_tbl) > 0) {
    print(removed_tbl %>% select(sample_id, site, missingness_percent), n = Inf)
  } else {
    cat("none\n")
  }
  
  cat("\nQC table written to: ", OUTPUT_FILE, "\n", sep = "")
}

main()