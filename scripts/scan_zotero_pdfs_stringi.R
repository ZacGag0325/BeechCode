# Scan exported Zotero PDFs for beech regeneration, root suckering,
# clonality, seed-origin regeneration, validation, disease, and management terms.
#
# Run from the BeechCode repository root with either:
# source("scripts/scan_zotero_pdfs_stringi.R")
# or:
# source("scripts/lit_review/scan_zotero_pdfs_stringi.R")
#
# IMPORTANT: All scores and yes/no/unclear fields created here are automatic
# text-mining triage suggestions only. They must be manually verified against
# the full article text before being used as final evidence or coded data.

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

normalize_start_dir <- function(path) {
  path <- normalizePath(path, mustWork = FALSE)
  if (file.exists(path) && !dir.exists(path)) {
    path <- dirname(path)
  }
  path
}

ancestor_dirs <- function(path) {
  path <- normalize_start_dir(path)
  dirs <- character(0)
  repeat {
    dirs <- c(dirs, path)
    parent <- dirname(path)
    if (identical(parent, path)) break
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
  candidate_roots <- candidate_roots[!grepl("\\.Rproj$", candidate_roots)]
  matches <- candidate_roots[file.exists(file.path(candidate_roots, "data/lit_review"))]
  if (length(matches) > 0) return(matches[[1]])
  normalize_start_dir(getwd())
}

repo_root <- find_repo_root()
path_from_root <- function(...) file.path(repo_root, ...)

zotero_pdf_dir_options <- c(
  path_from_root("data/lit_review/Beech_zotero_export"),
  path_from_root("data/lit_review/Beech_Zotero_export")
)
zotero_pdf_dir <- zotero_pdf_dir_options[dir.exists(zotero_pdf_dir_options)][1]
if (is.na(zotero_pdf_dir)) zotero_pdf_dir <- zotero_pdf_dir_options[[1]]

output_dir <- path_from_root("outputs/lit_review")
pdf_scan_dir <- file.path(output_dir, "pdf_scans")
hits_output_path <- file.path(pdf_scan_dir, "beech_pdf_stringi_hits.csv")
priority_output_path <- file.path(pdf_scan_dir, "beech_pdf_coding_priorities.csv")
summary_output_path <- file.path(pdf_scan_dir, "beech_pdf_text_mining_summary.csv")

# pdftools/Poppler can print non-fatal PDF repair/font diagnostics while still
# extracting text successfully. Set to FALSE if you need to debug one PDF.
suppress_pdf_diagnostics <- TRUE

codebook_path <- path_from_root("data/lit_review/coding/coding_codebook_revised.csv")
if (!file.exists(codebook_path)) {
  warning("Codebook not found at: ", codebook_path, call. = FALSE)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pdf_scan_dir, recursive = TRUE, showWarnings = FALSE)

term_groups <- list(
  root_suckering = c(
    "root sucker", "root suckers", "root suckering", "root sprout",
    "root sprouts", "root sprouting", "suckering", "sucker origin",
    "sprout origin", "originated from root", "root-origin", "root origin",
    "root-derived", "root derived", "vegetative reproduction",
    "vegetative propagation", "asexual reproduction", "adventitious shoot",
    "adventitious shoots"
  ),
  clonality = c(
    "clonal", "clonality", "clone", "clones", "genet", "genets",
    "ramet", "ramets", "multilocus genotype", "MLG",
    "multilocus lineage", "MLL"
  ),
  broad_seedling_regeneration = c(
    "seedling", "seedlings", "sapling", "saplings", "regeneration",
    "advance regeneration", "recruitment", "seedling bank", "understory",
    "understorey"
  ),
  explicit_seed_origin = c(
    "seed origin", "seed-origin", "seed-origin regeneration", "seed-derived",
    "seed derived", "from seed", "originated from seed", "sexual reproduction",
    "sexual recruitment", "sexual regeneration", "seedling origin",
    "germination", "seed rain", "seed dispersal", "beechnut", "beechnuts",
    "mast", "parentage", "pollen dispersal"
  ),
  distinguish_origin = c(
    "distinguish seed", "distinguished seed", "seedling versus sucker",
    "seedlings versus suckers", "seedling vs sucker", "seedlings vs root suckers",
    "seed origin and root sucker", "seed-origin and root-sucker", "seedling or sucker",
    "sucker or seedling", "root sucker origin", "seedling origin",
    "origin of regeneration", "mode of regeneration", "regeneration origin",
    "sexual and vegetative", "vegetative and sexual", "seed and root sucker",
    "root sucker and seed"
  ),
  genetic_validation = c(
    "microsatellite", "microsatellites", "allozyme", "allozymes", "ISSR",
    "SNP", "genotyped", "genotyping", "genotype", "genotypes",
    "genetic marker", "genetic markers", "multilocus genotype", "MLG",
    "multilocus lineage", "MLL", "parentage analysis", "genetic validation",
    "clonal assignment"
  ),
  physical_validation = c(
    "excavation", "excavated", "root tracing", "traced roots",
    "physical connection", "physically connected", "connected roots",
    "root connection", "root graft", "root grafting", "attached roots"
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

measurement_action_terms <- c(
  "density", "abundance", "counted", "measured", "compared", "sampled",
  "frequency", "proportion", "percent", "percentage", "survival", "growth",
  "morphology", "origin", "classified"
)

escape_regex <- function(x) {
  # Escape regex metacharacters before combining literal search terms into one pattern.
  metacharacters <- c("\\", ".", "|", "(", ")", "[", "]", "{", "}", "+", "*", "^", "$", "?")
  
  vapply(
    strsplit(x, split = "", fixed = TRUE),
    function(chars) {
      escaped_chars <- ifelse(chars %in% metacharacters, paste0("\\", chars), chars)
      paste0(escaped_chars, collapse = "")
    },
    character(1)
  )
}

make_regex <- function(terms) {
  paste0("\\b(", paste(escape_regex(terms), collapse = "|"), ")\\b")
}

term_patterns <- purrr::map_chr(term_groups, make_regex)
measurement_action_pattern <- make_regex(measurement_action_terms)

near_regex <- function(left_terms, right_terms, window_words = 12L) {
  left <- paste(escape_regex(left_terms), collapse = "|")
  right <- paste(escape_regex(right_terms), collapse = "|")
  between <- paste0("(?:\\W+\\w+){0,", window_words, "}\\W+")
  paste0("\\b(?:", left, ")\\b", between, "\\b(?:", right, ")\\b|\\b(?:", right, ")\\b", between, "\\b(?:", left, ")\\b")
}

root_measurement_pattern <- near_regex(term_groups$root_suckering, measurement_action_terms, window_words = 15L)
clonality_context_pattern <- near_regex(
  term_groups$clonality,
  c("Fagus", "beech", "regeneration", "root", "sucker", "suckering", "sprout", "seedling"),
  window_words = 12L
)

clean_pdf_text <- function(pdf_pages) {
  pdf_pages |>
    paste(collapse = " ") |>
    stringi::stri_replace_all_regex("\\s+", " ") |>
    stringi::stri_trim_both()
}

safe_extract_pdf_text <- function(pdf_path) {
  tryCatch(
    {
      pdf_pages <- NULL
      if (isTRUE(suppress_pdf_diagnostics)) {
        invisible(utils::capture.output(
          pdf_pages <- suppressWarnings(pdftools::pdf_text(pdf_path)),
          type = "message"
        ))
      } else {
        pdf_pages <- pdftools::pdf_text(pdf_path)
      }
      text <- clean_pdf_text(pdf_pages)
      list(text = text, text_extracted = TRUE, extraction_error = NA_character_)
    },
    error = function(e) {
      list(text = NA_character_, text_extracted = FALSE, extraction_error = conditionMessage(e))
    }
  )
}

count_hits <- function(text, pattern) {
  if (is.na(text) || !nzchar(text)) return(0L)
  as.integer(stringi::stri_count_regex(text, pattern, opts_regex = stringi::stri_opts_regex(case_insensitive = TRUE)))
}

first_snippet <- function(text, pattern, context_chars = 100L) {
  if (is.na(text) || !nzchar(text)) return(NA_character_)
  location <- stringi::stri_locate_first_regex(text, pattern, opts_regex = stringi::stri_opts_regex(case_insensitive = TRUE))
  if (is.na(location[, "start"])) return(NA_character_)
  snippet_start <- max(1L, location[, "start"] - context_chars)
  snippet_end <- min(stringi::stri_length(text), location[, "end"] + context_chars)
  stringi::stri_sub(text, snippet_start, snippet_end) |>
    stringi::stri_replace_all_regex("\\s+", " ") |>
    stringi::stri_trim_both()
}

suggest_yes_no <- function(condition) ifelse(condition, "yes", "no")

suggest_root_measured <- function(root_hits, root_measurement_hits) {
  dplyr::case_when(
    root_hits > 1 && root_measurement_hits > 0 ~ "yes",
    root_hits > 0 ~ "unclear",
    TRUE ~ "no"
  )
}

suggest_distinguishes <- function(distinguish_hits, root_hits, explicit_seed_hits) {
  dplyr::case_when(
    distinguish_hits > 0 ~ "yes",
    root_hits > 0 && explicit_seed_hits > 0 ~ "unclear",
    TRUE ~ "no"
  )
}

suggest_main_framing <- function(root_hits, clonality_hits, explicit_seed_hits, broad_seedling_hits) {
  root_clonal_hits <- root_hits + clonality_hits
  dplyr::case_when(
    root_hits > 0 && explicit_seed_hits > 0 ~ "both_root_and_seed",
    root_clonal_hits >= 3 && root_clonal_hits > explicit_seed_hits ~ "root_suckering_or_clonal",
    explicit_seed_hits >= 3 && root_clonal_hits <= 1 ~ "seed_origin_or_sexual",
    broad_seedling_hits > 0 && root_hits == 0 && explicit_seed_hits == 0 ~ "broad_seedling_regeneration",
    TRUE ~ "background_or_other"
  )
}

pdf_files <- character(0)
if (dir.exists(zotero_pdf_dir)) {
  pdf_files <- list.files(zotero_pdf_dir, pattern = "\\.pdf$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  pdf_files <- unique(normalizePath(pdf_files, mustWork = FALSE))
} else {
  warning("Zotero PDF directory not found at: ", zotero_pdf_dir, call. = FALSE)
}
if (length(pdf_files) == 0) warning("No PDFs found under: ", zotero_pdf_dir, call. = FALSE)

scan_one_pdf <- function(pdf_path) {
  extracted <- safe_extract_pdf_text(pdf_path)
  text <- extracted$text
  
  hit_counts <- purrr::imap(term_patterns, ~ count_hits(text, .x))
  snippets <- purrr::imap(term_patterns, ~ first_snippet(text, .x))
  root_measurement_hits <- count_hits(text, root_measurement_pattern)
  clonality_context_hits <- count_hits(text, clonality_context_pattern)
  
  root_suckering_discussed_suggested <- suggest_yes_no(hit_counts$root_suckering > 0)
  root_suckering_measured_suggested <- suggest_root_measured(hit_counts$root_suckering, root_measurement_hits)
  clonality_discussed_suggested <- suggest_yes_no(hit_counts$clonality > 0)
  seedling_regeneration_discussed_suggested <- suggest_yes_no(hit_counts$broad_seedling_regeneration > 0)
  explicit_seed_origin_discussed_suggested <- suggest_yes_no(hit_counts$explicit_seed_origin > 0)
  distinguishes_seed_vs_sucker_suggested <- suggest_distinguishes(hit_counts$distinguish_origin, hit_counts$root_suckering, hit_counts$explicit_seed_origin)
  genetic_validation_suggested <- suggest_yes_no(hit_counts$genetic_validation > 0)
  physical_validation_suggested <- suggest_yes_no(hit_counts$physical_validation > 0)
  main_framing_suggested <- suggest_main_framing(hit_counts$root_suckering, hit_counts$clonality, hit_counts$explicit_seed_origin, hit_counts$broad_seedling_regeneration)
  
  suggested_priority <- dplyr::case_when(
    hit_counts$root_suckering > 0 || hit_counts$distinguish_origin > 0 || hit_counts$genetic_validation > 0 || hit_counts$physical_validation > 0 ~ "1_core",
    hit_counts$broad_seedling_regeneration > 0 || hit_counts$explicit_seed_origin > 0 || hit_counts$bbd > 0 || hit_counts$management_disturbance > 0 ~ "2_high",
    hit_counts$bbd > 0 || hit_counts$management_disturbance > 0 ~ "3_background",
    TRUE ~ "4_low"
  )
  
  suggested_root_treatment_score <- dplyr::case_when(
    (hit_counts$genetic_validation > 0 || hit_counts$physical_validation > 0) && (hit_counts$root_suckering > 0 || hit_counts$clonality > 0) ~ 5L,
    root_suckering_measured_suggested == "yes" || hit_counts$root_suckering > 2 || clonality_context_hits > 2 ~ 4L,
    hit_counts$root_suckering > 0 || clonality_context_hits > 0 ~ 3L,
    TRUE ~ 0L
  )
  
  suggested_seed_treatment_score <- dplyr::case_when(
    hit_counts$genetic_validation > 0 && hit_counts$explicit_seed_origin > 0 ~ 5L,
    hit_counts$explicit_seed_origin > 5 ~ 4L,
    hit_counts$explicit_seed_origin > 0 ~ 3L,
    TRUE ~ 0L
  )
  
  snippet_priority <- c(
    snippets$distinguish_origin,
    snippets$root_suckering,
    snippets$genetic_validation,
    snippets$physical_validation,
    snippets$explicit_seed_origin,
    snippets$broad_seedling_regeneration
  )
  best_snippet <- snippet_priority[!is.na(snippet_priority) & nzchar(snippet_priority)][1]
  if (length(best_snippet) == 0) best_snippet <- NA_character_
  
  tibble(
    pdf_file = basename(pdf_path),
    pdf_path = pdf_path,
    text_extracted = extracted$text_extracted,
    extraction_error = extracted$extraction_error,
    n_characters = ifelse(is.na(text), 0L, stringi::stri_length(text)),
    root_suckering_hits = hit_counts$root_suckering,
    clonality_hits = hit_counts$clonality,
    clonality_context_hits = clonality_context_hits,
    broad_seedling_regeneration_hits = hit_counts$broad_seedling_regeneration,
    explicit_seed_origin_hits = hit_counts$explicit_seed_origin,
    distinguish_origin_hits = hit_counts$distinguish_origin,
    genetic_validation_hits = hit_counts$genetic_validation,
    physical_validation_hits = hit_counts$physical_validation,
    root_measurement_context_hits = root_measurement_hits,
    bbd_hits = hit_counts$bbd,
    management_disturbance_hits = hit_counts$management_disturbance,
    root_suckering_discussed_suggested = root_suckering_discussed_suggested,
    root_suckering_measured_suggested = root_suckering_measured_suggested,
    clonality_discussed_suggested = clonality_discussed_suggested,
    seedling_regeneration_discussed_suggested = seedling_regeneration_discussed_suggested,
    explicit_seed_origin_discussed_suggested = explicit_seed_origin_discussed_suggested,
    distinguishes_seed_vs_sucker_suggested = distinguishes_seed_vs_sucker_suggested,
    genetic_validation_suggested = genetic_validation_suggested,
    physical_validation_suggested = physical_validation_suggested,
    main_framing_suggested = main_framing_suggested,
    root_or_clonal_terms_present = hit_counts$root_suckering > 0 || hit_counts$clonality > 0,
    seed_or_regeneration_terms_present = hit_counts$broad_seedling_regeneration > 0 || hit_counts$explicit_seed_origin > 0,
    validation_terms_present = hit_counts$genetic_validation > 0 || hit_counts$physical_validation > 0,
    suggested_priority = suggested_priority,
    suggested_root_treatment_score = suggested_root_treatment_score,
    suggested_seed_treatment_score = suggested_seed_treatment_score,
    root_terms_present_suggested = root_suckering_discussed_suggested,
    seed_terms_present_suggested = suggest_yes_no(hit_counts$broad_seedling_regeneration > 0 || hit_counts$explicit_seed_origin > 0),
    root_suckering_snippet = snippets$root_suckering,
    clonality_snippet = snippets$clonality,
    broad_seedling_regeneration_snippet = snippets$broad_seedling_regeneration,
    explicit_seed_origin_snippet = snippets$explicit_seed_origin,
    distinguish_origin_snippet = snippets$distinguish_origin,
    genetic_validation_snippet = snippets$genetic_validation,
    physical_validation_snippet = snippets$physical_validation,
    bbd_snippet = snippets$bbd,
    management_disturbance_snippet = snippets$management_disturbance,
    best_snippet = best_snippet
  )
}

empty_scan_results <- tibble(
  pdf_file = character(), pdf_path = character(), text_extracted = logical(), extraction_error = character(), n_characters = integer(),
  root_suckering_hits = integer(), clonality_hits = integer(), clonality_context_hits = integer(), broad_seedling_regeneration_hits = integer(),
  explicit_seed_origin_hits = integer(), distinguish_origin_hits = integer(), genetic_validation_hits = integer(), physical_validation_hits = integer(),
  root_measurement_context_hits = integer(), bbd_hits = integer(), management_disturbance_hits = integer(),
  root_suckering_discussed_suggested = character(), root_suckering_measured_suggested = character(), clonality_discussed_suggested = character(),
  seedling_regeneration_discussed_suggested = character(), explicit_seed_origin_discussed_suggested = character(), distinguishes_seed_vs_sucker_suggested = character(),
  genetic_validation_suggested = character(), physical_validation_suggested = character(), main_framing_suggested = character(),
  root_or_clonal_terms_present = logical(), seed_or_regeneration_terms_present = logical(), validation_terms_present = logical(),
  suggested_priority = character(), suggested_root_treatment_score = integer(), suggested_seed_treatment_score = integer(),
  root_terms_present_suggested = character(), seed_terms_present_suggested = character(),
  root_suckering_snippet = character(), clonality_snippet = character(), broad_seedling_regeneration_snippet = character(), explicit_seed_origin_snippet = character(),
  distinguish_origin_snippet = character(), genetic_validation_snippet = character(), physical_validation_snippet = character(), bbd_snippet = character(),
  management_disturbance_snippet = character(), best_snippet = character()
)

pdf_scan_results <- if (length(pdf_files) == 0) empty_scan_results else purrr::map_dfr(pdf_files, scan_one_pdf)

manual_coding_priorities <- pdf_scan_results |>
  dplyr::select(
    pdf_file, pdf_path, text_extracted, n_characters, suggested_priority,
    suggested_root_treatment_score, suggested_seed_treatment_score,
    root_suckering_discussed_suggested, root_suckering_measured_suggested,
    clonality_discussed_suggested, seedling_regeneration_discussed_suggested,
    explicit_seed_origin_discussed_suggested, distinguishes_seed_vs_sucker_suggested,
    genetic_validation_suggested, physical_validation_suggested, main_framing_suggested,
    root_suckering_hits, clonality_hits, clonality_context_hits,
    broad_seedling_regeneration_hits, explicit_seed_origin_hits, distinguish_origin_hits,
    genetic_validation_hits, physical_validation_hits, root_measurement_context_hits,
    bbd_hits, management_disturbance_hits, best_snippet
  )

count_yes <- function(column) sum(column == "yes", na.rm = TRUE)
summary_metrics <- tibble(
  metric = c(
    "total_pdfs", "pdfs_successfully_extracted", "root_suckering_discussed_yes",
    "clonality_discussed_yes", "root_or_clonal_discussed_yes",
    "explicit_seed_origin_discussed_yes", "broad_seedling_regeneration_discussed_yes",
    "distinguishes_seed_vs_sucker_yes", "genetic_validation_yes", "physical_validation_yes"
  ),
  value = c(
    length(pdf_files), sum(pdf_scan_results$text_extracted, na.rm = TRUE),
    count_yes(pdf_scan_results$root_suckering_discussed_suggested),
    count_yes(pdf_scan_results$clonality_discussed_suggested),
    sum(pdf_scan_results$root_suckering_discussed_suggested == "yes" | pdf_scan_results$clonality_discussed_suggested == "yes", na.rm = TRUE),
    count_yes(pdf_scan_results$explicit_seed_origin_discussed_suggested),
    count_yes(pdf_scan_results$seedling_regeneration_discussed_suggested),
    count_yes(pdf_scan_results$distinguishes_seed_vs_sucker_suggested),
    count_yes(pdf_scan_results$genetic_validation_suggested),
    count_yes(pdf_scan_results$physical_validation_suggested)
  ),
  interpretation = c(
    "Unique PDF files discovered recursively in the Zotero export folder.",
    "PDFs where pdftools returned text without stopping the workflow.",
    "Automatic yes count for root suckering / vegetative reproduction terms.",
    "Automatic yes count for clonality terms; generic clone/genotype hits require manual review.",
    "Automatic yes count for either root suckering or clonality discussion.",
    "Automatic yes count for explicit seed-origin / sexual regeneration terms.",
    "Automatic yes count for broad seedling/sapling/regeneration terms.",
    "Automatic yes count for language distinguishing seed-origin stems from root suckers.",
    "Automatic yes count for genetic validation terms.",
    "Automatic yes count for physical validation terms."
  )
)

main_framing_counts <- pdf_scan_results |>
  dplyr::count(main_framing_suggested, name = "value") |>
  dplyr::mutate(metric = paste0("main_framing_", main_framing_suggested), interpretation = "Number of PDFs assigned to this automatic main framing category.") |>
  dplyr::select(metric, value, interpretation)

text_mining_summary <- dplyr::bind_rows(summary_metrics, main_framing_counts)

readr::write_csv(pdf_scan_results, hits_output_path, na = "")
readr::write_csv(manual_coding_priorities, priority_output_path, na = "")
readr::write_csv(text_mining_summary, summary_output_path, na = "")

priority_counts <- pdf_scan_results |> dplyr::count(suggested_priority, name = "n")
main_framing_summary <- pdf_scan_results |> dplyr::count(main_framing_suggested, name = "n")

message("\nBeech Zotero PDF stringi scan complete")
message("--------------------------------------")
message("PDFs found: ", length(pdf_files))
message("PDF diagnostics suppressed: ", suppress_pdf_diagnostics)
message("PDFs successfully extracted: ", sum(pdf_scan_results$text_extracted, na.rm = TRUE))
message("PDFs with root suckering discussed: ", count_yes(pdf_scan_results$root_suckering_discussed_suggested))
message("PDFs with clonality discussed: ", count_yes(pdf_scan_results$clonality_discussed_suggested))
message("PDFs with explicit seed-origin discussed: ", count_yes(pdf_scan_results$explicit_seed_origin_discussed_suggested))
message("PDFs with broad seedling regeneration discussed: ", count_yes(pdf_scan_results$seedling_regeneration_discussed_suggested))
message("PDFs that distinguish seed vs sucker: ", count_yes(pdf_scan_results$distinguishes_seed_vs_sucker_suggested))
message("PDFs with genetic validation: ", count_yes(pdf_scan_results$genetic_validation_suggested))
message("PDFs with physical validation: ", count_yes(pdf_scan_results$physical_validation_suggested))
message("Suggested priority counts:")
if (nrow(priority_counts) == 0) message("  none") else purrr::pwalk(priority_counts, ~ message("  ", ..1, ": ", ..2))
message("Main framing counts:")
if (nrow(main_framing_summary) == 0) message("  none") else purrr::pwalk(main_framing_summary, ~ message("  ", ..1, ": ", ..2))
message("Full hit output: ", hits_output_path)
message("Manual coding priority output: ", priority_output_path)
message("Text-mining summary output: ", summary_output_path)
message("\nWarning: this is text-mining triage, not final evidence. Manually verify all suggested fields against the articles.")