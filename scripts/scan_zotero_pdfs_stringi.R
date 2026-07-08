# Scan exported Zotero PDFs for beech regeneration, root suckering,
# clonality, seed recruitment, validation, disease, and management terms.
#
# Run from the BeechCode repository with:
# source("scripts/lit_review/scan_zotero_pdfs_stringi.R")
#
# If you copy this file to scripts/scan_zotero_pdfs_stringi.R, it will also
# try to find the repository root from the script location so relative paths
# still resolve correctly.
#
# IMPORTANT: All scores and yes/no fields created here are automatic
# suggestions for triage only. They must be manually verified against the
# full article text before being used for coding or analysis.

required_packages <- c(
  "pdftools",
  "stringi",
  "dplyr",
  "purrr",
  "tibble",
  "readr",
  "stringr"
)

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(
    "Missing required package(s): ",
    paste(missing_packages, collapse = ", "),
    ". Please install them before running this script.",
    call. = FALSE
  )
}

library(pdftools)
library(stringi)
library(dplyr)
library(purrr)
library(tibble)
library(readr)
library(stringr)

script_path <- tryCatch(normalizePath(sys.frame(1)$ofile, mustWork = TRUE), error = function(e) NA_character_)

ancestor_dirs <- function(path) {
  path <- normalizePath(path, mustWork = FALSE)
  dirs <- character(0)
  
  repeat {
    dirs <- c(dirs, path)
    parent <- dirname(path)
    if (identical(parent, path)) {
      break
    }
    path <- parent
  }
  
  unique(dirs)
}

find_repo_root <- function() {
  candidate_starts <- c(getwd())
  if (!is.na(script_path)) {
    candidate_starts <- c(candidate_starts, dirname(script_path))
  }
  
  candidate_roots <- unique(unlist(purrr::map(candidate_starts, ancestor_dirs)))
  matches <- candidate_roots[file.exists(file.path(candidate_roots, "data/lit_review"))]
  
  if (length(matches) > 0) {
    return(matches[[1]])
  }
  
  normalizePath(getwd(), mustWork = FALSE)
}

repo_root <- find_repo_root()
path_from_root <- function(...) file.path(repo_root, ...)

zotero_pdf_dir_options <- c(
  path_from_root("data/lit_review/Beech_zotero_export"),
  path_from_root("data/lit_review/Beech_Zotero_export")
)
zotero_pdf_dir <- zotero_pdf_dir_options[dir.exists(zotero_pdf_dir_options)][1]
if (is.na(zotero_pdf_dir)) {
  zotero_pdf_dir <- zotero_pdf_dir_options[[1]]
}

output_dir <- path_from_root("outputs/lit_review")
pdf_scan_dir <- file.path(output_dir, "pdf_scans")
hits_output_path <- file.path(pdf_scan_dir, "beech_pdf_stringi_hits.csv")
priority_output_path <- file.path(pdf_scan_dir, "beech_pdf_coding_priorities.csv")

# Keep the codebook path explicit so users can compare these suggested fields
# with the manual coding definitions while reviewing results.
codebook_path <- path_from_root("data/lit_review/coding/coding_codebook_revised.csv")

if (!file.exists(codebook_path)) {
  warning("Codebook not found at: ", codebook_path, call. = FALSE)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pdf_scan_dir, recursive = TRUE, showWarnings = FALSE)

term_groups <- list(
  root_suckering = c(
    "root sucker", "root suckers", "root suckering", "root sprout",
    "root sprouts", "root sprouting", "suckers", "suckering",
    "sprout origin", "sucker origin", "vegetative reproduction",
    "vegetative propagation", "asexual reproduction", "adventitious shoots"
  ),
  clonality = c(
    "clonal", "clone", "clones", "clonality", "genet", "genets",
    "ramet", "ramets", "multilocus genotype", "MLG",
    "multilocus lineage", "MLL"
  ),
  seed_recruitment = c(
    "seedling", "seedlings", "seed origin", "seed-origin", "seed recruitment",
    "sexual reproduction", "sexual recruitment", "germination", "seed rain",
    "mast", "beechnut", "beechnuts", "parentage", "pollen", "seed dispersal"
  ),
  regeneration = c(
    "regeneration", "advance regeneration", "sapling", "saplings",
    "understory", "understorey", "seedling bank", "recruitment", "thicket",
    "beech thicket", "beech regeneration"
  ),
  validation = c(
    "microsatellite", "microsatellites", "allozyme", "ISSR", "SNP",
    "genotype", "genotyping", "genetic marker", "genetic markers",
    "parentage analysis", "excavation", "excavated", "root tracing",
    "physical connection", "connected roots", "root connection", "root graft"
  ),
  bbd = c(
    "beech bark disease", "BBD", "beech scale", "Cryptococcus fagisuga",
    "Neonectria", "Nectria"
  ),
  management_disturbance = c(
    "herbicide", "glyphosate", "harvest", "shelterwood", "selection cutting",
    "canopy gap", "gap", "deer", "browsing", "fire", "ice storm",
    "hurricane", "disturbance"
  )
)

escape_regex <- function(x) {
  # stringi does not export a regex-escape helper in all installed versions,
  # so escape regex metacharacters locally before building the term pattern.
  stringr::str_replace_all(x, "([][{}()+*^$|\\.?])", "\\\\\\1")
}

make_regex <- function(terms) {
  paste0("\\b(", paste(escape_regex(terms), collapse = "|"), ")\\b")
}

term_patterns <- purrr::map_chr(term_groups, make_regex)

clean_pdf_text <- function(pdf_pages) {
  pdf_pages |>
    paste(collapse = " ") |>
    stringi::stri_replace_all_regex("\\s+", " ") |>
    stringi::stri_trim_both()
}

safe_extract_pdf_text <- function(pdf_path) {
  tryCatch(
    {
      text <- pdftools::pdf_text(pdf_path) |>
        clean_pdf_text()
      
      list(
        text = text,
        text_extracted = TRUE,
        extraction_error = NA_character_
      )
    },
    error = function(e) {
      list(
        text = NA_character_,
        text_extracted = FALSE,
        extraction_error = conditionMessage(e)
      )
    }
  )
}

count_hits <- function(text, pattern) {
  if (is.na(text) || !nzchar(text)) {
    return(0L)
  }
  
  matches <- stringi::stri_count_regex(text, pattern, opts_regex = stringi::stri_opts_regex(case_insensitive = TRUE))
  as.integer(matches)
}

first_snippet <- function(text, pattern, context_chars = 100L) {
  if (is.na(text) || !nzchar(text)) {
    return(NA_character_)
  }
  
  location <- stringi::stri_locate_first_regex(
    text,
    pattern,
    opts_regex = stringi::stri_opts_regex(case_insensitive = TRUE)
  )
  
  if (is.na(location[, "start"])) {
    return(NA_character_)
  }
  
  snippet_start <- max(1L, location[, "start"] - context_chars)
  snippet_end <- min(stringi::stri_length(text), location[, "end"] + context_chars)
  
  stringi::stri_sub(text, snippet_start, snippet_end) |>
    stringi::stri_replace_all_regex("\\s+", " ") |>
    stringi::stri_trim_both()
}

suggest_yes_no <- function(condition) {
  ifelse(condition, "yes", "no")
}

pdf_files <- character(0)
if (dir.exists(zotero_pdf_dir)) {
  pdf_files <- list.files(
    zotero_pdf_dir,
    pattern = "\\.pdf$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )
} else {
  warning("Zotero PDF directory not found at: ", zotero_pdf_dir, call. = FALSE)
}

if (length(pdf_files) == 0) {
  warning("No PDFs found under: ", zotero_pdf_dir, call. = FALSE)
}

scan_one_pdf <- function(pdf_path) {
  extracted <- safe_extract_pdf_text(pdf_path)
  text <- extracted$text
  
  hit_counts <- purrr::imap(term_patterns, ~ count_hits(text, .x))
  snippets <- purrr::imap(term_patterns, ~ first_snippet(text, .x))
  
  root_or_clonal_hits <- hit_counts$root_suckering > 0 || hit_counts$clonality > 0
  seed_or_regeneration_hits <- hit_counts$seed_recruitment > 0 || hit_counts$regeneration > 0
  validation_hits_present <- hit_counts$validation > 0
  
  genetic_validation_terms <- c(
    "microsatellite", "microsatellites", "allozyme", "ISSR", "SNP",
    "genotype", "genotyping", "genetic marker", "genetic markers",
    "parentage analysis", "parentage"
  )
  physical_validation_terms <- c(
    "excavation", "excavated", "root tracing", "physical connection",
    "connected roots", "root connection", "root graft"
  )
  
  genetic_validation_hits <- count_hits(text, make_regex(genetic_validation_terms))
  physical_validation_hits <- count_hits(text, make_regex(physical_validation_terms))
  
  suggested_priority <- dplyr::case_when(
    hit_counts$root_suckering > 0 || hit_counts$clonality > 1 || hit_counts$validation > 0 ~ "1_core",
    hit_counts$regeneration > 0 || hit_counts$bbd > 0 || hit_counts$management_disturbance > 0 ~ "2_high",
    hit_counts$seed_recruitment > 0 ~ "3_background",
    TRUE ~ "4_low"
  )
  
  suggested_root_treatment_score <- dplyr::case_when(
    validation_hits_present && root_or_clonal_hits ~ 5L,
    hit_counts$root_suckering > 2 || hit_counts$clonality > 2 ~ 4L,
    root_or_clonal_hits ~ 3L,
    TRUE ~ 0L
  )
  
  # Automatic suggestion only: parentage/genetic validation terms co-occurring
  # anywhere in the extracted text with seed terms do not prove validated sexual
  # recruitment without manual review.
  suggested_seed_treatment_score <- dplyr::case_when(
    genetic_validation_hits > 0 && hit_counts$seed_recruitment > 0 ~ 5L,
    hit_counts$seed_recruitment > 5 ~ 4L,
    hit_counts$seed_recruitment > 0 ~ 3L,
    TRUE ~ 0L
  )
  
  snippet_values <- unlist(snippets, use.names = TRUE)
  best_snippet <- snippet_values[!is.na(snippet_values) & nzchar(snippet_values)][1]
  if (length(best_snippet) == 0) {
    best_snippet <- NA_character_
  }
  
  tibble(
    pdf_file = basename(pdf_path),
    pdf_path = pdf_path,
    text_extracted = extracted$text_extracted,
    extraction_error = extracted$extraction_error,
    n_characters = ifelse(is.na(text), 0L, stringi::stri_length(text)),
    root_suckering_hits = hit_counts$root_suckering,
    clonality_hits = hit_counts$clonality,
    seed_recruitment_hits = hit_counts$seed_recruitment,
    regeneration_hits = hit_counts$regeneration,
    validation_hits = hit_counts$validation,
    bbd_hits = hit_counts$bbd,
    management_disturbance_hits = hit_counts$management_disturbance,
    root_or_clonal_terms_present = root_or_clonal_hits,
    seed_or_regeneration_terms_present = seed_or_regeneration_hits,
    validation_terms_present = validation_hits_present,
    suggested_priority = suggested_priority,
    suggested_root_treatment_score = suggested_root_treatment_score,
    suggested_seed_treatment_score = suggested_seed_treatment_score,
    root_terms_present_suggested = suggest_yes_no(root_or_clonal_hits),
    seed_terms_present_suggested = suggest_yes_no(seed_or_regeneration_hits),
    genetic_validation_suggested = suggest_yes_no(genetic_validation_hits > 0),
    physical_validation_suggested = suggest_yes_no(physical_validation_hits > 0),
    root_suckering_snippet = snippets$root_suckering,
    clonality_snippet = snippets$clonality,
    seed_recruitment_snippet = snippets$seed_recruitment,
    regeneration_snippet = snippets$regeneration,
    validation_snippet = snippets$validation,
    bbd_snippet = snippets$bbd,
    management_disturbance_snippet = snippets$management_disturbance,
    best_snippet = best_snippet
  )
}

empty_scan_results <- tibble(
  pdf_file = character(),
  pdf_path = character(),
  text_extracted = logical(),
  extraction_error = character(),
  n_characters = integer(),
  root_suckering_hits = integer(),
  clonality_hits = integer(),
  seed_recruitment_hits = integer(),
  regeneration_hits = integer(),
  validation_hits = integer(),
  bbd_hits = integer(),
  management_disturbance_hits = integer(),
  root_or_clonal_terms_present = logical(),
  seed_or_regeneration_terms_present = logical(),
  validation_terms_present = logical(),
  suggested_priority = character(),
  suggested_root_treatment_score = integer(),
  suggested_seed_treatment_score = integer(),
  root_terms_present_suggested = character(),
  seed_terms_present_suggested = character(),
  genetic_validation_suggested = character(),
  physical_validation_suggested = character(),
  root_suckering_snippet = character(),
  clonality_snippet = character(),
  seed_recruitment_snippet = character(),
  regeneration_snippet = character(),
  validation_snippet = character(),
  bbd_snippet = character(),
  management_disturbance_snippet = character(),
  best_snippet = character()
)

pdf_scan_results <- if (length(pdf_files) == 0) {
  empty_scan_results
} else {
  purrr::map_dfr(pdf_files, scan_one_pdf)
}

manual_coding_priorities <- pdf_scan_results |>
  dplyr::select(
    pdf_file,
    pdf_path,
    text_extracted,
    n_characters,
    suggested_priority,
    suggested_root_treatment_score,
    suggested_seed_treatment_score,
    root_terms_present_suggested,
    seed_terms_present_suggested,
    genetic_validation_suggested,
    physical_validation_suggested,
    root_suckering_hits,
    clonality_hits,
    seed_recruitment_hits,
    regeneration_hits,
    validation_hits,
    bbd_hits,
    management_disturbance_hits,
    best_snippet
  )

readr::write_csv(pdf_scan_results, hits_output_path, na = "")
readr::write_csv(manual_coding_priorities, priority_output_path, na = "")

priority_counts <- pdf_scan_results |>
  dplyr::count(suggested_priority, name = "n")

message("\nBeech Zotero PDF stringi scan complete")
message("--------------------------------------")
message("PDFs found: ", length(pdf_files))
message("PDFs successfully extracted: ", sum(pdf_scan_results$text_extracted, na.rm = TRUE))
message(
  "PDFs with root/clonal hits: ",
  sum(pdf_scan_results$root_or_clonal_terms_present, na.rm = TRUE)
)
message(
  "PDFs with validation hits: ",
  sum(pdf_scan_results$validation_terms_present, na.rm = TRUE)
)
message("Suggested priority counts:")
if (nrow(priority_counts) == 0) {
  message("  none")
} else {
  purrr::pwalk(priority_counts, ~ message("  ", ..1, ": ", ..2))
}
message("Full hit output: ", hits_output_path)
message("Manual coding priority output: ", priority_output_path)
message("\nReminder: suggested scores and yes/no fields are automatic triage aids only and require manual verification.")