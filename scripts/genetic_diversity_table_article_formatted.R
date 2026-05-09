#!/usr/bin/env Rscript

############################################################
# Genetic diversity and inbreeding table for article/thesis use
#
# Standalone script. It searches for an MLL-/clone-corrected
# genetic object, computes site-level diversity metrics, writes
# CSV/RDS appendix tables, and exports a Word-formatted main table.
#
# Main outputs:
# - outputs/tables/genetic_diversity_table_article_formatted.csv
# - outputs/tables/genetic_diversity_table_article_formatted.rds
# - outputs/word/genetic_diversity_table_article_formatted.docx
#
# Appendix outputs:
# - outputs/tables/allele_counts_by_locus_and_site.csv
# - outputs/tables/allele_counts_by_locus_and_site.rds
# - outputs/tables/genetic_diversity_by_locus_and_site.csv
# - outputs/tables/genetic_diversity_by_locus_and_site.rds
############################################################

required_packages <- c("adegenet", "hierfstat", "dplyr", "tidyr", "stringr", "readxl")
optional_packages <- c("officer", "flextable")

missing_required <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_required) > 0) {
  stop(
    "Missing required package(s): ", paste(missing_required, collapse = ", "),
    ". Please install them before running this script.",
    call. = FALSE
  )
}

has_officer <- requireNamespace("officer", quietly = TRUE)
has_flextable <- requireNamespace("flextable", quietly = TRUE)

suppressPackageStartupMessages({
  library(adegenet)
  library(hierfstat)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readxl)
})

# -----------------------------
# Console warning collector
# -----------------------------
script_warnings <- character(0)

add_warning <- function(...) {
  msg <- paste0(...)
  script_warnings <<- unique(c(script_warnings, msg))
  warning(msg, call. = FALSE, immediate. = TRUE)
}

add_note_warning <- function(...) {
  msg <- paste0(...)
  script_warnings <<- unique(c(script_warnings, msg))
  message("WARNING: ", msg)
}

# -----------------------------
# Generic helpers
# -----------------------------
normalize_ascii <- function(x) {
  x |>
    iconv(from = "", to = "ASCII//TRANSLIT") |>
    tolower() |>
    str_replace_all("[^a-z0-9]+", "_") |>
    str_replace_all("_+", "_") |>
    str_replace_all("^_|_$", "")
}

normalize_key <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("\uFEFF", "", x, fixed = TRUE)
  x <- gsub("[[:cntrl:]]", "", x)
  x[is.na(x)] <- ""
  x
}

normalize_id <- function(x) toupper(normalize_key(x))

safe_mean <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (length(x) == 0 || all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

round3 <- function(x) round(as.numeric(x), 3)

format3 <- function(x) {
  x_num <- suppressWarnings(as.numeric(x))
  ifelse(is.na(x_num), NA_character_, formatC(round(x_num, 3), format = "f", digits = 3))
}

collapse_or_none <- function(x) {
  x <- unique(x[!is.na(x) & nzchar(as.character(x))])
  if (length(x) == 0) "None" else paste(x, collapse = "; ")
}

pick_column <- function(df, patterns, label = "column", required = FALSE, exclude = character()) {
  cn <- names(df)
  cn_norm <- normalize_ascii(cn)
  exclude_norm <- normalize_ascii(exclude)
  idx <- which(str_detect(cn_norm, str_c(patterns, collapse = "|")) & !(cn_norm %in% exclude_norm))
  if (length(idx) == 0) {
    if (required) add_warning("Could not detect required ", label, ".")
    return(NA_character_)
  }
  if (length(idx) > 1) {
    message("Multiple matches for ", label, ": ", paste(cn[idx], collapse = ", "), ". Using: ", cn[idx[1]])
  }
  cn[idx[1]]
}

# -----------------------------
# Project root and directories
# -----------------------------
resolve_script_path <- function() {
  cmd_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
  if (length(cmd_file) > 0 && nzchar(cmd_file[1])) {
    return(normalizePath(cmd_file[1], winslash = "/", mustWork = FALSE))
  }
  frame_files <- vapply(sys.frames(), function(x) {
    val <- x$ofile
    if (is.null(val)) NA_character_ else as.character(val)
  }, character(1))
  frame_files <- frame_files[!is.na(frame_files) & nzchar(frame_files)]
  if (length(frame_files) > 0) {
    return(normalizePath(frame_files[length(frame_files)], winslash = "/", mustWork = FALSE))
  }
  NA_character_
}

find_project_root <- function() {
  script_path <- resolve_script_path()
  starts <- c(getwd())
  if (!is.na(script_path)) starts <- c(dirname(script_path), starts)
  starts <- unique(normalizePath(starts, winslash = "/", mustWork = FALSE))
  
  root_markers <- function(path) {
    dir.exists(file.path(path, "scripts")) &&
      (dir.exists(file.path(path, "data")) || file.exists(file.path(path, "scripts", "_load_objects.R")) ||
         length(list.files(path, pattern = "\\.Rproj$", full.names = TRUE)) > 0)
  }
  
  for (start in starts) {
    cur <- start
    repeat {
      if (root_markers(cur)) return(cur)
      parent <- dirname(cur)
      if (identical(parent, cur)) break
      cur <- parent
    }
  }
  
  stop("Cannot detect BeechCode project root. Run this script from inside the project or keep it under scripts/.", call. = FALSE)
}

PROJECT_ROOT <- find_project_root()
setwd(PROJECT_ROOT)

TABLES_DIR <- file.path(PROJECT_ROOT, "outputs", "tables")
WORD_DIR <- file.path(PROJECT_ROOT, "outputs", "word")
dir.create(TABLES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(WORD_DIR, recursive = TRUE, showWarnings = FALSE)

message("[genetic_diversity] Project root: ", PROJECT_ROOT)
message("[genetic_diversity] Tables directory: ", TABLES_DIR)
message("[genetic_diversity] Word directory: ", WORD_DIR)

# -----------------------------
# site_lookup loading and site labels/order
# -----------------------------
find_site_lookup_workbook <- function(raw_dir = file.path(PROJECT_ROOT, "data", "raw"), sheet = "site_lookup") {
  if (!dir.exists(raw_dir)) {
    add_note_warning("site_lookup unavailable because data/raw does not exist: ", raw_dir)
    return(NULL)
  }
  excel_files <- list.files(raw_dir, pattern = "\\.(xlsx|xls)$", full.names = TRUE, ignore.case = TRUE)
  if (length(excel_files) == 0) {
    add_note_warning("site_lookup unavailable because no Excel workbook was found in data/raw.")
    return(NULL)
  }
  
  has_lookup <- vapply(excel_files, function(path) {
    sheets <- tryCatch(readxl::excel_sheets(path), error = function(e) character(0))
    any(normalize_ascii(sheets) == normalize_ascii(sheet))
  }, logical(1))
  
  lookup_files <- excel_files[has_lookup]
  if (length(lookup_files) == 0) {
    add_note_warning(
      "site_lookup unavailable because no workbook in data/raw contains a site_lookup sheet. Checked: ",
      paste(basename(excel_files), collapse = ", ")
    )
    return(NULL)
  }
  if (length(lookup_files) > 1) {
    add_note_warning(
      "Multiple workbooks contain site_lookup; using ", basename(lookup_files[1]),
      ". Other matches: ", paste(basename(lookup_files[-1]), collapse = ", ")
    )
  } else {
    message("[genetic_diversity] Reading site_lookup from: ", basename(lookup_files[1]))
  }
  lookup_files[1]
}

load_site_lookup <- function() {
  workbook <- find_site_lookup_workbook()
  if (is.null(workbook)) return(NULL)
  sheets <- readxl::excel_sheets(workbook)
  sheet <- sheets[match(normalize_ascii("site_lookup"), normalize_ascii(sheets))]
  lookup <- suppressMessages(readxl::read_excel(workbook, sheet = sheet)) |>
    as.data.frame(stringsAsFactors = FALSE)
  
  site_col <- pick_column(lookup, c("^site$", "site_code", "old_site", "code_site"), "site_lookup old Site", required = TRUE)
  label_col <- pick_column(lookup, c("site_label", "new_site", "label", "site_id", "^site_new$"), "site_lookup new site label")
  if (is.na(label_col)) label_col <- site_col
  order_col <- pick_column(lookup, c("site_order", "^order$", "ordre", "sort", "south.*north", "north.*south"), "site_lookup order")
  lat_col <- pick_column(lookup, c("^lat", "latitude", "y_coord", "gps_lat"), "site_lookup latitude")
  
  out <- lookup |>
    mutate(
      old_site = normalize_key(.data[[site_col]]),
      site_label = normalize_key(.data[[label_col]]),
      site_order = if (!is.na(order_col)) suppressWarnings(as.numeric(.data[[order_col]])) else NA_real_,
      latitude = if (!is.na(lat_col)) suppressWarnings(as.numeric(str_replace_all(as.character(.data[[lat_col]]), ",", "."))) else NA_real_
    ) |>
    filter(nzchar(old_site), nzchar(site_label)) |>
    distinct(old_site, .keep_all = TRUE)
  
  if (nrow(out) == 0) {
    add_note_warning("site_lookup was found but no usable Site/Site_label rows were detected.")
    return(NULL)
  }
  
  message("[genetic_diversity] Loaded site_lookup rows: ", nrow(out))
  out
}

site_lookup <- load_site_lookup()

natural_site_rank <- function(x) {
  x_chr <- as.character(x)
  prefix <- str_extract(x_chr, "^[A-Za-z]+")
  num <- suppressWarnings(as.integer(str_extract(x_chr, "[0-9]+")))
  prefix_rank <- case_when(
    toupper(prefix) == "S" ~ 1L,
    toupper(prefix) == "N" ~ 2L,
    TRUE ~ 99L
  )
  ifelse(is.na(num), 999L, prefix_rank * 100L + num)
}

map_site_labels <- function(site_values) {
  site_values <- normalize_key(site_values)
  if (is.null(site_lookup)) return(site_values)
  idx <- match(site_values, site_lookup$old_site)
  labels <- site_lookup$site_label[idx]
  ifelse(is.na(labels) | !nzchar(labels), site_values, labels)
}

site_order_table <- function(site_values) {
  labels <- unique(map_site_labels(site_values))
  out <- data.frame(Site = labels, stringsAsFactors = FALSE)
  
  if (!is.null(site_lookup)) {
    lookup_order <- site_lookup |>
      transmute(Site = site_label, lookup_order = site_order, latitude = latitude) |>
      distinct(Site, .keep_all = TRUE)
    out <- out |>
      left_join(lookup_order, by = "Site")
  } else {
    out$lookup_order <- NA_real_
    out$latitude <- NA_real_
  }
  
  if (any(!is.na(out$lookup_order))) {
    out <- out |> arrange(is.na(lookup_order), lookup_order, natural_site_rank(Site), Site)
    message("[genetic_diversity] Sorting sites using site_lookup order column.")
  } else if (any(!is.na(out$latitude))) {
    out <- out |> arrange(is.na(latitude), latitude, natural_site_rank(Site), Site)
    message("[genetic_diversity] Sorting sites south-to-north using site_lookup latitude.")
  } else {
    out <- out |> arrange(natural_site_rank(Site), Site)
    message("[genetic_diversity] Sorting sites by natural S1-S6, N1-N6 fallback order.")
  }
  
  out |> mutate(display_order = row_number())
}

# -----------------------------
# Genetic object discovery
# -----------------------------
likely_dirs <- c(
  file.path(PROJECT_ROOT, "outputs", "v1", "objects"),
  file.path(PROJECT_ROOT, "outputs", "objects"),
  file.path(PROJECT_ROOT, "outputs", "tables"),
  file.path(PROJECT_ROOT, "data", "processed")
)

existing_likely_dirs <- likely_dirs[dir.exists(likely_dirs)]
files_searched <- unlist(lapply(existing_likely_dirs, function(d) {
  list.files(d, pattern = "\\.(rds|RDS|rda|RData)$", full.names = TRUE, recursive = FALSE)
}), use.names = FALSE)

candidate_object_log <- character(0)

is_genetic_object <- function(x) inherits(x, c("genind", "genclone"))

score_candidate <- function(name, path, object, source_type) {
  text <- tolower(paste(name, basename(path), dirname(path), source_type, collapse = " "))
  score <- 0
  if (str_detect(text, "gi_mll|mll|clone_correct|clone.correct|clonecorrect|clonal|clone_corrected")) score <- score + 100
  if (str_detect(text, "genind|genclone|gi")) score <- score + 20
  if (str_detect(text, "full|raw|ramet|non_clone|nonclone")) score <- score - 30
  if (inherits(object, "genclone")) score <- score + 10
  score
}

collect_candidates_from_file <- function(path) {
  out <- list()
  ext <- tolower(tools::file_ext(path))
  
  if (ext == "rds") {
    obj <- tryCatch(readRDS(path), error = function(e) e)
    if (is_genetic_object(obj)) {
      nm <- tools::file_path_sans_ext(basename(path))
      out[[length(out) + 1]] <- list(
        object = obj,
        object_name = nm,
        path = path,
        source_type = "rds",
        score = score_candidate(nm, path, obj, "rds")
      )
      candidate_object_log <<- c(candidate_object_log, paste0(path, " :: ", nm, " [genetic object]"))
    } else {
      candidate_object_log <<- c(candidate_object_log, paste0(path, " [not a genind/genclone RDS]"))
    }
  } else if (ext %in% c("rda", "rdata")) {
    env <- new.env(parent = emptyenv())
    loaded_names <- tryCatch(load(path, envir = env), error = function(e) character(0))
    if (length(loaded_names) == 0) {
      candidate_object_log <<- c(candidate_object_log, paste0(path, " [could not load or empty]"))
    }
    for (nm in loaded_names) {
      obj <- get(nm, envir = env)
      if (is_genetic_object(obj)) {
        out[[length(out) + 1]] <- list(
          object = obj,
          object_name = nm,
          path = path,
          source_type = "RData",
          score = score_candidate(nm, path, obj, "RData")
        )
        candidate_object_log <<- c(candidate_object_log, paste0(path, " :: ", nm, " [genetic object]"))
      } else {
        candidate_object_log <<- c(candidate_object_log, paste0(path, " :: ", nm, " [not genind/genclone]"))
      }
    }
  }
  
  out
}

candidates <- unlist(lapply(files_searched, collect_candidates_from_file), recursive = FALSE)

if (length(candidates) == 0) {
  stop(
    "No genind/genclone genetic object could be found.\n",
    "Directories searched:\n - ", paste(likely_dirs, collapse = "\n - "), "\n",
    "Files/objects checked:\n - ", ifelse(length(candidate_object_log) == 0, "None", paste(candidate_object_log, collapse = "\n - ")),
    call. = FALSE
  )
}

candidate_scores <- vapply(candidates, function(x) x$score, numeric(1))
selected <- candidates[[which.max(candidate_scores)]]
gen_obj <- selected$object
used_dataset_label <- paste0(selected$object_name, " from ", selected$path)
clone_corrected <- selected$score >= 100

if (!clone_corrected) {
  add_warning(
    "No clearly clone-corrected/MLL-corrected genetic dataset was found. Using the best available full dataset instead: ",
    used_dataset_label
  )
} else {
  message("[genetic_diversity] Selected clone-corrected/MLL-corrected dataset: ", used_dataset_label)
}

# -----------------------------
# Metadata discovery for Site assignment if pop(gen_obj) is missing
# -----------------------------
collect_metadata_candidates <- function() {
  metadata <- list()
  
  for (path in files_searched) {
    ext <- tolower(tools::file_ext(path))
    if (ext == "rds") {
      obj <- tryCatch(readRDS(path), error = function(e) NULL)
      if (is.data.frame(obj)) {
        metadata[[length(metadata) + 1]] <- list(df = obj, name = tools::file_path_sans_ext(basename(path)), path = path)
      }
    } else if (ext %in% c("rda", "rdata")) {
      env <- new.env(parent = emptyenv())
      loaded_names <- tryCatch(load(path, envir = env), error = function(e) character(0))
      for (nm in loaded_names) {
        obj <- get(nm, envir = env)
        if (is.data.frame(obj)) metadata[[length(metadata) + 1]] <- list(df = obj, name = nm, path = path)
      }
    }
  }
  
  metadata
}

assign_pop_from_metadata <- function(gobj) {
  p <- adegenet::pop(gobj)
  if (!is.null(p) && length(p) == adegenet::nInd(gobj) && any(!is.na(p))) {
    adegenet::pop(gobj) <- as.factor(normalize_key(as.character(p)))
    return(gobj)
  }
  
  metadata <- collect_metadata_candidates()
  if (length(metadata) == 0) {
    stop("The selected genetic object has no usable population/site vector and no metadata data.frame was found in searched object files.", call. = FALSE)
  }
  
  for (md in metadata) {
    df <- md$df
    id_col <- pick_column(df, c("^ind$", "individual", "sample", "sample_id", "sampleid", "^id$", "ind_id", "tree_id"), "metadata individual ID")
    site_col <- pick_column(df, c("^site$", "site_code", "population", "^pop$", "numero_population", "num_ro_population"), "metadata Site")
    if (is.na(id_col) || is.na(site_col)) next
    
    mapper <- setNames(normalize_key(df[[site_col]]), normalize_id(df[[id_col]]))
    sites <- mapper[normalize_id(adegenet::indNames(gobj))]
    if (length(sites) == adegenet::nInd(gobj) && all(!is.na(sites)) && all(nzchar(sites))) {
      adegenet::pop(gobj) <- as.factor(sites)
      message("[genetic_diversity] Assigned Site/pop from metadata object ", md$name, " in ", md$path)
      return(gobj)
    }
  }
  
  stop("Could not assign Site/pop to the selected genetic object from available metadata. Checked metadata in searched RDS/RData files.", call. = FALSE)
}

gen_obj <- assign_pop_from_metadata(gen_obj)

original_sites <- normalize_key(as.character(adegenet::pop(gen_obj)))
new_sites <- map_site_labels(original_sites)
adegenet::pop(gen_obj) <- as.factor(new_sites)
site_order <- site_order_table(new_sites)

missing_lookup_sites <- if (!is.null(site_lookup)) setdiff(unique(original_sites), site_lookup$old_site) else character(0)
if (length(missing_lookup_sites) > 0) {
  add_note_warning("These genetic object site codes were not found in site_lookup and were kept as-is: ", paste(missing_lookup_sites, collapse = ", "))
}

# -----------------------------
# Diversity calculations
# -----------------------------
message("[genetic_diversity] Number of individuals used: ", adegenet::nInd(gen_obj))
message("[genetic_diversity] Number of loci used: ", adegenet::nLoc(gen_obj))
message("[genetic_diversity] Sites detected: ", paste(levels(adegenet::pop(gen_obj)), collapse = ", "))

pop_sizes <- table(as.character(adegenet::pop(gen_obj)))
low_sample_sites <- names(pop_sizes)[pop_sizes < 5]
if (length(low_sample_sites) > 0) {
  add_note_warning("Low sample size (N < 5) for site(s): ", paste(low_sample_sites, collapse = ", "))
}
if (any(pop_sizes < 2)) {
  add_note_warning("At least one site has N < 2; some Ho/He/FIS/Ar estimates may be NA.")
}

hf <- hierfstat::genind2hierfstat(gen_obj)
bs <- hierfstat::basic.stats(hf)

site_levels <- rownames(bs$Ho)
site_n_tbl <- table(as.character(adegenet::pop(gen_obj)))

allele_tab <- adegenet::tab(gen_obj, NA.method = "asis")
loc_fac <- adegenet::locFac(gen_obj)
loci <- adegenet::locNames(gen_obj)

allele_count_one <- function(rows, locus_name) {
  cols <- which(loc_fac == locus_name)
  if (length(cols) == 0 || length(rows) == 0) return(NA_integer_)
  mat <- allele_tab[rows, cols, drop = FALSE]
  if (all(is.na(mat))) return(NA_integer_)
  as.integer(sum(colSums(mat, na.rm = TRUE) > 0))
}

site_row_indices <- split(seq_len(nrow(allele_tab)), as.character(adegenet::pop(gen_obj)))

allele_counts_long <- tidyr::expand_grid(
  Locus = loci,
  Site = names(site_row_indices)
) |>
  rowwise() |>
  mutate(Allele_Count = allele_count_one(site_row_indices[[Site]], Locus)) |>
  ungroup()

overall_counts <- data.frame(
  Locus = loci,
  Overall = vapply(loci, function(loc) allele_count_one(seq_len(nrow(allele_tab)), loc), integer(1)),
  stringsAsFactors = FALSE
)

allele_counts_by_locus_and_site <- allele_counts_long |>
  left_join(site_order, by = "Site") |>
  arrange(display_order, Locus) |>
  select(Locus, Site, Allele_Count) |>
  tidyr::pivot_wider(names_from = Site, values_from = Allele_Count) |>
  left_join(overall_counts, by = "Locus") |>
  arrange(Locus)

mean_alleles_by_site <- allele_counts_long |>
  group_by(Site) |>
  summarise(A = safe_mean(Allele_Count), .groups = "drop")

min_n <- suppressWarnings(min(as.integer(pop_sizes), na.rm = TRUE))
allelic_richness_by_site <- data.frame(Site = names(pop_sizes), Ar = NA_real_, stringsAsFactors = FALSE)
if (is.finite(min_n) && min_n >= 2) {
  ar_obj <- tryCatch(hierfstat::allelic.richness(hf, min.n = min_n), error = function(e) e)
  if (inherits(ar_obj, "error")) {
    add_note_warning("Rarefied allelic richness calculation failed: ", conditionMessage(ar_obj))
  } else if (!is.null(ar_obj$Ar)) {
    ar_mat <- ar_obj$Ar
    allelic_richness_by_site <- data.frame(
      Site = colnames(ar_mat),
      Ar = as.numeric(colMeans(ar_mat, na.rm = TRUE)),
      stringsAsFactors = FALSE
    )
  } else {
    add_note_warning("hierfstat::allelic.richness did not return an Ar matrix; Ar set to NA.")
  }
} else {
  add_note_warning("Rarefied allelic richness was not calculated because the smallest site sample size is < 2.")
}

heterozygosity_by_site <- data.frame(
  Site = site_levels,
  N = as.integer(site_n_tbl[site_levels]),
  Ho = apply(bs$Ho, 1, safe_mean),
  He = apply(bs$Hs, 1, safe_mean),
  FIS = apply(bs$Fis, 1, safe_mean),
  stringsAsFactors = FALSE
)

main_table_numeric <- heterozygosity_by_site |>
  left_join(mean_alleles_by_site, by = "Site") |>
  left_join(allelic_richness_by_site, by = "Site") |>
  left_join(site_order, by = "Site") |>
  arrange(display_order, Site) |>
  transmute(
    Site = Site,
    N = as.integer(N),
    A = round3(A),
    Ar = round3(Ar),
    Ho = round3(Ho),
    He = round3(He),
    FIS = round3(FIS)
  )

# Keep numeric CSV/RDS values rounded, but format Word/console values with trailing zeros.
main_table_word <- main_table_numeric |>
  mutate(
    A = format3(A),
    Ar = format3(Ar),
    Ho = format3(Ho),
    He = format3(He),
    FIS = format3(FIS)
  )

per_locus_diversity <- NULL
if (!is.null(bs$Ho) && !is.null(bs$Hs) && !is.null(bs$Fis)) {
  per_locus_diversity <- tidyr::expand_grid(
    Site = rownames(bs$Ho),
    Locus = colnames(bs$Ho)
  ) |>
    rowwise() |>
    mutate(
      N = as.integer(site_n_tbl[Site]),
      A = allele_count_one(site_row_indices[[Site]], Locus),
      Ho = suppressWarnings(as.numeric(bs$Ho[Site, Locus])),
      He = suppressWarnings(as.numeric(bs$Hs[Site, Locus])),
      FIS = suppressWarnings(as.numeric(bs$Fis[Site, Locus]))
    ) |>
    ungroup() |>
    left_join(site_order, by = "Site") |>
    arrange(display_order, Locus) |>
    select(Site, Locus, N, A, Ho, He, FIS)
} else {
  add_note_warning("Per-locus Ho/He/FIS table could not be created because basic.stats output did not include expected matrices.")
}

# -----------------------------
# Word export helpers
# -----------------------------
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
  text <- xml_escape_word(value)
  bold_xml <- if (bold) "<w:rPr><w:b/></w:rPr>" else ""
  paste0(
    "<w:tc>",
    "<w:tcPr><w:tcW w:w=\"1800\" w:type=\"dxa\"/></w:tcPr>",
    "<w:p><w:r>", bold_xml, "<w:t>", text, "</w:t></w:r></w:p>",
    "</w:tc>"
  )
}

word_paragraph_xml <- function(value, bold = FALSE) {
  bold_xml <- if (bold) "<w:rPr><w:b/></w:rPr>" else ""
  paste0("<w:p><w:r>", bold_xml, "<w:t>", xml_escape_word(value), "</w:t></w:r></w:p>")
}

write_basic_docx_table <- function(df, path, title, note = NULL) {
  tmp_dir <- tempfile("genetic_diversity_docx_")
  dir.create(file.path(tmp_dir, "_rels"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(tmp_dir, "word", "_rels"), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  
  df <- as.data.frame(lapply(df, as.character), stringsAsFactors = FALSE, check.names = FALSE)
  rows <- apply(df, 1, function(row) {
    paste0("<w:tr>", paste(vapply(row, word_cell_xml, character(1)), collapse = ""), "</w:tr>")
  })
  header <- paste0("<w:tr>", paste(vapply(names(df), word_cell_xml, character(1), bold = TRUE), collapse = ""), "</w:tr>")
  note_xml <- if (!is.null(note) && nzchar(note)) word_paragraph_xml(note) else ""
  
  document_xml <- paste0(
    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
    '<w:document xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main" ',
    'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">',
    '<w:body>',
    word_paragraph_xml(title, bold = TRUE),
    '<w:tbl>',
    '<w:tblPr><w:tblStyle w:val="TableGrid"/><w:tblW w:w="0" w:type="auto"/>',
    '<w:tblBorders>',
    '<w:top w:val="single" w:sz="4" w:space="0" w:color="auto"/>',
    '<w:left w:val="single" w:sz="4" w:space="0" w:color="auto"/>',
    '<w:bottom w:val="single" w:sz="4" w:space="0" w:color="auto"/>',
    '<w:right w:val="single" w:sz="4" w:space="0" w:color="auto"/>',
    '<w:insideH w:val="single" w:sz="4" w:space="0" w:color="auto"/>',
    '<w:insideV w:val="single" w:sz="4" w:space="0" w:color="auto"/>',
    '</w:tblBorders></w:tblPr>',
    header,
    paste(rows, collapse = ""),
    '</w:tbl>',
    note_xml,
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
  utils::zip(
    zipfile = path,
    files = list.files(tmp_dir, recursive = TRUE, all.files = TRUE, no.. = TRUE),
    flags = "-q"
  )
  if (!file.exists(path)) stop("Failed to create Word document at: ", path, call. = FALSE)
  invisible(path)
}

write_word_table <- function(df, path, title, note) {
  note_for_doc <- note
  if (!clone_corrected) {
    note_for_doc <- paste0(
      note_for_doc,
      " Warning: a clearly clone-corrected/MLL-corrected genetic dataset was not found, so the best available full dataset was used."
    )
  }
  if (length(script_warnings) > 0) {
    note_for_doc <- paste0(note_for_doc, " Additional warnings: ", paste(script_warnings, collapse = " | "))
  }
  
  status <- tryCatch({
    if (has_officer && has_flextable) {
      ft <- flextable::flextable(df) |>
        flextable::theme_booktabs() |>
        flextable::autofit() |>
        flextable::align(align = "center", part = "all") |>
        flextable::align(j = "Site", align = "left", part = "body") |>
        flextable::bold(part = "header") |>
        flextable::fit_to_width(max_width = 6.5)
      
      doc <- officer::read_docx()
      doc <- officer::body_add_par(doc, title, style = "heading 1")
      doc <- officer::body_add_flextable(doc, ft)
      doc <- officer::body_add_par(doc, note_for_doc, style = "Normal")
      print(doc, target = path)
      "created with flextable"
    } else if (has_officer) {
      add_note_warning("flextable is not available/loadable; creating the Word table with officer::body_add_table instead.")
      doc <- officer::read_docx()
      doc <- officer::body_add_par(doc, title, style = "heading 1")
      doc <- officer::body_add_table(doc, value = as.data.frame(df), style = "table_template")
      doc <- officer::body_add_par(doc, note_for_doc, style = "Normal")
      print(doc, target = path)
      "created with officer fallback"
    } else {
      add_note_warning("Neither officer nor flextable is loadable; creating a basic Word-compatible .docx file instead.")
      write_basic_docx_table(df, path, title = title, note = note_for_doc)
      "created with built-in docx fallback"
    }
  }, error = function(e) {
    add_note_warning("Package-based Word export failed: ", conditionMessage(e), ". Creating a basic Word-compatible .docx file instead.")
    write_basic_docx_table(df, path, title = title, note = note_for_doc)
    "created with built-in docx fallback after package export failed"
  })
  
  if (!file.exists(path)) stop("Word table was not created: ", path, call. = FALSE)
  status
}

# -----------------------------
# Save outputs
# -----------------------------
main_csv <- file.path(TABLES_DIR, "genetic_diversity_table_article_formatted.csv")
main_rds <- file.path(TABLES_DIR, "genetic_diversity_table_article_formatted.rds")
main_docx <- file.path(WORD_DIR, "genetic_diversity_table_article_formatted.docx")
allele_csv <- file.path(TABLES_DIR, "allele_counts_by_locus_and_site.csv")
allele_rds <- file.path(TABLES_DIR, "allele_counts_by_locus_and_site.rds")
per_locus_csv <- file.path(TABLES_DIR, "genetic_diversity_by_locus_and_site.csv")
per_locus_rds <- file.path(TABLES_DIR, "genetic_diversity_by_locus_and_site.rds")

table_caption <- "Table X. Genetic diversity and inbreeding coefficients by site."
table_note <- "N = number of retained individuals or multilocus lineages used in the analysis; A = mean number of alleles per locus; Ar = rarefied allelic richness; Ho = observed heterozygosity; He = expected heterozygosity; FIS = inbreeding coefficient."

write.csv(main_table_numeric, main_csv, row.names = FALSE, na = "")
saveRDS(main_table_numeric, main_rds)
write.csv(allele_counts_by_locus_and_site, allele_csv, row.names = FALSE, na = "")
saveRDS(allele_counts_by_locus_and_site, allele_rds)

created_files <- c(main_csv, main_rds, allele_csv, allele_rds)
if (!is.null(per_locus_diversity)) {
  write.csv(per_locus_diversity, per_locus_csv, row.names = FALSE, na = "")
  saveRDS(per_locus_diversity, per_locus_rds)
  created_files <- c(created_files, per_locus_csv, per_locus_rds)
}

docx_status <- write_word_table(main_table_word, main_docx, title = table_caption, note = table_note)
created_files <- c(created_files, main_docx)

# -----------------------------
# Final console output
# -----------------------------
message("\nFinal formatted main table:")
print(main_table_word, row.names = FALSE)

message("\nFiles created:")
for (path in created_files) {
  extra <- if (identical(path, main_docx)) paste0(" (", docx_status, ")") else ""
  message(" - ", path, extra)
}

message("\nGenetic dataset used:")
message(" - ", used_dataset_label)
message(" - clone-corrected/MLL-corrected detected: ", ifelse(clone_corrected, "Yes", "No"))
message(" - individuals used: ", adegenet::nInd(gen_obj))
message(" - loci used: ", adegenet::nLoc(gen_obj))
message(" - rarefaction sample size for Ar: ", ifelse(is.finite(min_n) && min_n >= 2, min_n, "not calculated"))

message("\nWarnings and calculation notes:")
if (length(script_warnings) == 0) {
  message(" - None")
} else {
  for (msg in script_warnings) message(" - ", msg)
}

message("\nDone.")