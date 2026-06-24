# 01_merge_scopus_exports_revised.R
# Merge revised Scopus searches S1-S5 for American beech literature review

library(readr)
library(dplyr)
library(stringr)
library(purrr)
library(tibble)

# -----------------------------
# Paths
# -----------------------------

scopus_dir <- file.path("data", "lit_review", "scopus_exports")
screening_dir <- file.path("data", "lit_review", "screening")
output_table_dir <- file.path("outputs", "lit_review", "tables")

dir.create(screening_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_table_dir, recursive = TRUE, showWarnings = FALSE)

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
  file_prefix = c(
    "scopus_S1_main_regeneration_corpus",
    "scopus_S2_root_suckering_clonality",
    "scopus_S3_seed_recruitment",
    "scopus_S4_bbd_management_proliferation",
    "scopus_S5_exact_root_sucker_terms"
  )
)

find_file <- function(prefix) {
  hits <- list.files(
    scopus_dir,
    pattern = paste0("^", prefix, "\\.csv(\\.csv)?$"),
    full.names = TRUE
  )
  
  if (length(hits) == 0) {
    stop("Missing file for prefix: ", prefix)
  }
  
  if (length(hits) > 1) {
    message("Multiple files found for ", prefix, ". Using first one:")
    message(paste(hits, collapse = "\n"))
  }
  
  hits[1]
}

file_specs <- file_specs %>%
  mutate(file_path = map_chr(file_prefix, find_file))

# -----------------------------
# Helper functions
# -----------------------------

norm_name <- function(x) {
  x %>%
    tolower() %>%
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
    str_squish()
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
    
    authors = pick_col(raw, c("Authors", "Author(s)", "Author full names")),
    year = pick_col(raw, c("Year", "Publication Year")),
    title = pick_col(raw, c("Title", "Document Title", "Article Title")),
    source_title = pick_col(raw, c("Source title", "Journal", "Publication Name")),
    doi = pick_col(raw, c("DOI", "Digital Object Identifier")),
    abstract = pick_col(raw, c("Abstract", "Abstract Note")),
    author_keywords = pick_col(raw, c("Author Keywords", "Author keywords")),
    index_keywords = pick_col(raw, c("Index Keywords", "Indexed Keywords", "Index keywords")),
    document_type = pick_col(raw, c("Document Type", "Document type")),
    cited_by = pick_col(raw, c("Cited by", "Cited by count"))
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

# Clean DOI and title for deduplication
all_records <- all_records %>%
  mutate(
    doi_clean = doi %>%
      str_to_lower() %>%
      str_trim() %>%
      na_if(""),
    title_clean = clean_title(title),
    dedup_key = if_else(
      !is.na(doi_clean),
      paste0("doi:", doi_clean),
      paste0("title:", title_clean)
    )
  ) %>%
  filter(!is.na(title_clean), title_clean != "")

# Deduplicate, keeping search source information
master <- all_records %>%
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
  arrange(title) %>%
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
# Summary table
# -----------------------------

search_counts <- all_records %>%
  count(search_id, search_name, search_file, name = "raw_records")

summary_table <- tibble(
  metric = c(
    "Raw Scopus records total",
    "Unique records after deduplication",
    "Duplicates removed"
  ),
  count = c(raw_n, unique_n, duplicates_removed)
)

# -----------------------------
# Write outputs
# -----------------------------

master_out <- file.path(
  screening_dir,
  "01_master_scopus_deduplicated_screening_revised.csv"
)

summary_out <- file.path(
  output_table_dir,
  "scopus_dedup_summary_revised.csv"
)

search_counts_out <- file.path(
  output_table_dir,
  "scopus_search_counts_revised.csv"
)

write_csv(master, master_out)
write_csv(summary_table, summary_out)
write_csv(search_counts, search_counts_out)

message("\nDone.")
message("Raw records: ", raw_n)
message("Unique records: ", unique_n)
message("Duplicates removed: ", duplicates_removed)
message("\nWrote:")
message(master_out)
message(summary_out)
message(search_counts_out)