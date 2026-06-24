# 01_merge_scopus_exports_revised.R
# Merge revised Scopus searches S1-S5 for American beech systematic review.

library(readr)
library(dplyr)
library(stringr)
library(purrr)
library(tibble)

# -----------------------------
# Paths
# -----------------------------

get_script_path <- function() {
  frame_files <- map_chr(sys.frames(), ~ {
    ofile <- .x$ofile
    if (is.null(ofile)) {
      NA_character_
    } else {
      as.character(ofile)[1]
    }
  })
  frame_files <- frame_files[!is.na(frame_files) & frame_files != ""]
  
  if (length(frame_files) > 0) {
    return(normalizePath(frame_files[length(frame_files)], winslash = "/", mustWork = FALSE))
  }
  
  command_file <- commandArgs(trailingOnly = FALSE)
  file_arg <- command_file[str_detect(command_file, "^--file=")]
  if (length(file_arg) > 0) {
    return(normalizePath(str_remove(file_arg[1], "^--file="), winslash = "/", mustWork = FALSE))
  }
  
  NA_character_
}

script_path <- get_script_path()
project_root <- if (!is.na(script_path)) {
  normalizePath(file.path(dirname(script_path), "..", ".."), winslash = "/", mustWork = FALSE)
} else {
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

scopus_dir <- file.path(project_root, "data", "lit_review", "scopus_exports")
screening_dir <- file.path(project_root, "data", "lit_review", "screening")
coding_dir <- file.path(project_root, "data", "lit_review", "coding")
output_table_dir <- file.path(project_root, "outputs", "lit_review", "tables")

walk(
  c(scopus_dir, screening_dir, coding_dir, output_table_dir),
  ~ dir.create(.x, recursive = TRUE, showWarnings = FALSE)
)

# -----------------------------
# Expected revised Scopus files
# -----------------------------

file_specs <- tibble(
  search_id = c("S1", "S2", "S3", "S4", "S5"),
  search_name = c(
    "Main regeneration corpus",
    "Root suckering clonality",
    "Seed recruitment",
    "BBD management proliferation",
    "Exact root sucker terms"
  ),
  expected_file = c(
    "scopus_S1_main_regeneration_corpus.csv",
    "scopus_S2_root_suckering_clonality.csv",
    "scopus_S3_seed_recruitment.csv",
    "scopus_S4_bbd_management_proliferation.csv",
    "scopus_S5_exact_root_sucker_terms.csv"
  )
) %>%
  mutate(file_prefix = str_remove(expected_file, "\\.csv$"))

print_scopus_diagnostics <- function() {
  found_files <- list.files(scopus_dir, full.names = FALSE, all.files = FALSE)
  
  message("Current working directory: ", normalizePath(getwd(), winslash = "/", mustWork = FALSE))
  if (!is.na(script_path)) {
    message("Script path: ", script_path)
  }
  message("Project root used by script: ", project_root)
  message(
    "Expected Scopus folder path: ",
    normalizePath(scopus_dir, winslash = "/", mustWork = FALSE)
  )
  message("Files currently found in ", scopus_dir, ":")
  
  if (length(found_files) == 0) {
    message("  <none>")
  } else {
    walk(sort(found_files), ~ message("  - ", .x))
  }
}

find_file <- function(prefix) {
  hits <- list.files(
    scopus_dir,
    pattern = paste0("^", prefix, "\\.csv(\\.csv)?$"),
    full.names = TRUE,
    ignore.case = TRUE
  )
  
  if (length(hits) == 0) {
    return(NA_character_)
  }
  
  if (length(hits) > 1) {
    message("Multiple files found for ", prefix, ". Using first one:")
    message(paste(hits, collapse = "\n"))
  }
  
  hits[1]
}

file_specs <- file_specs %>%
  mutate(file_path = map_chr(file_prefix, find_file))

missing_files <- file_specs %>%
  filter(is.na(file_path))

if (nrow(missing_files) > 0) {
  print_scopus_diagnostics()
  stop(
    "Missing required Scopus export(s):\n",
    paste0("  - ", missing_files$expected_file, collapse = "\n"),
    "\nAccidental '.csv.csv' endings are accepted, but each S1-S5 export must be present.",
    call. = FALSE
  )
}

# -----------------------------
# Helper functions
# -----------------------------

norm_name <- function(x) {
  x %>%
    str_to_lower() %>%
    str_replace_all("[^a-z0-9]", "")
}

pick_col <- function(df, candidates) {
  nms <- names(df)
  nms_norm <- norm_name(nms)
  cand_norm <- norm_name(candidates)
  
  idx <- match(cand_norm, nms_norm)
  idx <- idx[!is.na(idx)]
  
  if (length(idx) == 0) {
    return(rep(NA_character_, nrow(df)))
  }
  
  as.character(df[[idx[1]]])
}

first_nonempty <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & str_squish(x) != ""]
  if (length(x) == 0) return(NA_character_)
  x[1]
}

clean_title <- function(x) {
  x %>%
    str_to_lower() %>%
    str_replace_all("[^a-z0-9]+", " ") %>%
    str_squish() %>%
    na_if("")
}

clean_doi <- function(x) {
  x %>%
    str_to_lower() %>%
    str_trim() %>%
    str_remove("^https?://(dx\\.)?doi\\.org/") %>%
    str_remove("^doi:\\s*") %>%
    na_if("")
}

make_citation <- function(authors, year, title) {
  case_when(
    !is.na(authors) & authors != "" & !is.na(year) & year != "" ~ paste0(authors, " (", year, "). ", title),
    !is.na(authors) & authors != "" ~ paste0(authors, ". ", title),
    TRUE ~ title
  )
}

# -----------------------------
# Read one Scopus file
# -----------------------------

read_scopus_file <- function(file_path, search_id, search_name) {
  message("Reading ", search_id, ": ", basename(file_path))
  
  raw <- read_csv(file_path, show_col_types = FALSE, guess_max = 10000)
  
  tibble(
    search_id = search_id,
    search_name = search_name,
    search_file = basename(file_path),
    authors = pick_col(raw, c("Authors", "Author(s)", "Author full names", "Author Names")),
    year = pick_col(raw, c("Year", "Publication Year", "PubYear", "Date")),
    title = pick_col(raw, c("Title", "Document Title", "Article Title", "Article title")),
    source_title = pick_col(raw, c("Source title", "Source Title", "Journal", "Publication Name", "Source")),
    doi = pick_col(raw, c("DOI", "Digital Object Identifier", "DOI Link")),
    abstract = pick_col(raw, c("Abstract", "Abstract Note", "Description")),
    author_keywords = pick_col(raw, c("Author Keywords", "Author keywords", "Authors Keywords", "Keywords")),
    index_keywords = pick_col(raw, c("Index Keywords", "Indexed Keywords", "Index keywords", "Index Terms")),
    document_type = pick_col(raw, c("Document Type", "Document type", "Publication Type", "Type")),
    cited_by = pick_col(raw, c("Cited by", "Cited By", "Cited by count", "Times Cited"))
  )
}

# -----------------------------
# Merge all searches
# -----------------------------

all_records <- pmap_dfr(
  list(file_specs$file_path, file_specs$search_id, file_specs$search_name),
  read_scopus_file
)

raw_n <- nrow(all_records)

all_records <- all_records %>%
  mutate(
    doi_clean = clean_doi(doi),
    title_clean = clean_title(title),
    dedup_key = case_when(
      !is.na(doi_clean) ~ paste0("doi:", doi_clean),
      !is.na(title_clean) ~ paste0("title:", title_clean),
      TRUE ~ NA_character_
    )
  )

undedupeable_n <- sum(is.na(all_records$dedup_key))
if (undedupeable_n > 0) {
  message("Dropping ", undedupeable_n, " record(s) with no DOI or title to deduplicate.")
}

all_records_for_dedup <- all_records %>%
  filter(!is.na(dedup_key))

master <- all_records_for_dedup %>%
  group_by(dedup_key) %>%
  summarise(
    search_source = paste(sort(unique(search_id)), collapse = "; "),
    search_names = paste(sort(unique(search_name)), collapse = "; "),
    search_files = paste(sort(unique(search_file)), collapse = "; "),
    authors = first_nonempty(authors),
    year = first_nonempty(year),
    title = first_nonempty(title),
    source_title = first_nonempty(source_title),
    doi = first_nonempty(doi),
    abstract = first_nonempty(abstract),
    author_keywords = first_nonempty(author_keywords),
    index_keywords = first_nonempty(index_keywords),
    document_type = first_nonempty(document_type),
    cited_by = first_nonempty(cited_by),
    .groups = "drop"
  ) %>%
  arrange(str_to_lower(title), year) %>%
  mutate(
    article_id = sprintf("AB_%04d", row_number()),
    include_title_abstract = "",
    exclusion_reason = "",
    screening_notes = ""
  ) %>%
  select(
    article_id,
    search_source,
    search_names,
    authors,
    year,
    title,
    source_title,
    doi,
    abstract,
    author_keywords,
    index_keywords,
    document_type,
    cited_by,
    include_title_abstract,
    exclusion_reason,
    screening_notes,
    search_files
  )

unique_n <- nrow(master)
duplicates_removed <- raw_n - unique_n

# -----------------------------
# Summary tables and coding files
# -----------------------------

search_counts <- all_records %>%
  count(search_id, search_name, search_file, name = "raw_records") %>%
  arrange(search_id)

summary_table <- tibble(
  metric = c(
    "Raw Scopus records total",
    "Records without DOI or title dropped before deduplication",
    "Unique records after deduplication",
    "Duplicates removed"
  ),
  count = c(raw_n, undedupeable_n, unique_n, duplicates_removed)
)

coding_template <- master %>%
  transmute(
    article_id,
    citation = make_citation(authors, year, title),
    year,
    title,
    region = "",
    study_context = "",
    study_type = "",
    main_regeneration_frame = "",
    root_terms_present = "",
    seed_terms_present = "",
    root_term_count_rough = "",
    seed_term_count_rough = "",
    root_terms_location = "",
    seed_terms_location = "",
    root_treatment_score = "",
    seed_treatment_score = "",
    root_framing = "",
    seed_framing = "",
    root_evidence_type = "",
    seed_evidence_type = "",
    does_article_assume_suckering = "",
    does_article_distinguish_seed_vs_sucker = "",
    compares_seed_and_sucker = "",
    genetic_validation = "",
    physical_validation = "",
    main_claim_about_regeneration = "",
    important_quote_root = "",
    important_quote_seed = "",
    page_root = "",
    page_seed = "",
    relevance_to_thesis = "",
    notes = ""
  )

codebook <- tribble(
  ~field, ~allowed_value_or_score, ~description,
  "root_treatment_score", "0", "Root suckering / clonality not mentioned.",
  "root_treatment_score", "1", "Brief mention only.",
  "root_treatment_score", "2", "Discussed as background/context.",
  "root_treatment_score", "3", "Framed as important, common, dominant, or problematic.",
  "root_treatment_score", "4", "Empirically measured.",
  "root_treatment_score", "5", "Genetically or physically validated.",
  "seed_treatment_score", "0", "Seed recruitment not mentioned.",
  "seed_treatment_score", "1", "Brief mention only.",
  "seed_treatment_score", "2", "Discussed as background/context.",
  "seed_treatment_score", "3", "Framed as important, common, dominant, or limiting.",
  "seed_treatment_score", "4", "Empirically measured.",
  "seed_treatment_score", "5", "Genetically/parentage-validated sexual recruitment.",
  "main_regeneration_frame", "free text", "Coder summary of whether the article frames regeneration mainly as root suckering/clonality, seed recruitment, both, neither, or uncertain.",
  "root_terms_present", "yes/no/unclear", "Whether root suckering, sprouting, clonal, proliferation, ramet, genet, or equivalent root-origin terms appear.",
  "seed_terms_present", "yes/no/unclear", "Whether seed, seedling, sexual recruitment, parentage, germination, mast, or equivalent seed-origin terms appear.",
  "does_article_assume_suckering", "yes/no/unclear", "Whether the article assumes stems/regeneration are suckers without validation.",
  "does_article_distinguish_seed_vs_sucker", "yes/no/unclear", "Whether the article explicitly separates seed-origin and sucker-origin regeneration.",
  "compares_seed_and_sucker", "yes/no/unclear", "Whether the article compares seed recruitment with root suckering/clonal regeneration.",
  "genetic_validation", "yes/no/unclear", "Whether genetic, parentage, marker, or clonality analyses validate origin.",
  "physical_validation", "yes/no/unclear", "Whether excavation, root tracing, or direct physical evidence validates origin."
)

# -----------------------------
# Write outputs
# -----------------------------

master_out <- file.path(screening_dir, "01_master_scopus_deduplicated_screening_revised.csv")
summary_out <- file.path(output_table_dir, "scopus_dedup_summary_revised.csv")
search_counts_out <- file.path(output_table_dir, "scopus_search_counts_revised.csv")
coding_template_out <- file.path(coding_dir, "02_full_text_coding_template_revised.csv")
codebook_out <- file.path(coding_dir, "coding_codebook_revised.csv")

write_csv(master, master_out)
write_csv(summary_table, summary_out)
write_csv(search_counts, search_counts_out)
write_csv(coding_template, coding_template_out)
write_csv(codebook, codebook_out)

# -----------------------------
# Clean terminal summary
# -----------------------------

message("\nDone.")
message("Raw records total: ", raw_n)
message("Records per search:")
pwalk(search_counts, ~ message("  ", ..1, " (", ..2, "): ", ..4, " [", ..3, "]"))
message("Unique records after deduplication: ", unique_n)
message("Duplicates removed: ", duplicates_removed)
message("\nWrote files:")
walk(
  c(master_out, summary_out, search_counts_out, coding_template_out, codebook_out),
  ~ message("  - ", .x)
)