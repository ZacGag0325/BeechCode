# scripts/00_master_pipeline.R
############################################################
# Build canonical objects for downstream analyses:
# - outputs/v1/objects/gi.rds
# - outputs/v1/objects/gi_mll.rds
# - outputs/v1/objects/df_ids.rds
# - outputs/v1/objects/meta.rds
#
# Source policy:
# - Prefer real workbook sources in data/raw/ and inputs/.
# - Do not require inputs/meta_ind.csv if workbook metadata exist.
# - Exclude genotype IDs outside the canonical metadata universe.
#
# Microsatellite clone-correction policy:
# - gi = full dataset for clonality and individual-level diagnostics.
#   Individuals with >35% missing allele data are removed from gi immediately
#   after import and before any clone detection or downstream analysis.
# - gi_mll = clone-corrected dataset based on MLLs defined from Bruvo distance.
# - Exact MLG assignments are retained alongside MLL assignments in df_ids.
############################################################

suppressPackageStartupMessages({
  library(adegenet)
  library(poppr)
})

BRUVO_MLL_THRESHOLD <- 0.09
BRUVO_ALGORITHM <- "farthest_neighbor"
MISSING_ALLELE_FILTER_THRESHOLD <- 0.35

find_project_root <- function() {
  candidates <- c(getwd(), normalizePath(file.path(getwd(), ".."), mustWork = FALSE))
  cmd_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
  if (length(cmd_file) > 0 && nzchar(cmd_file[1])) {
    candidates <- c(candidates, dirname(normalizePath(cmd_file[1], mustWork = FALSE)))
  }
  for (start in unique(candidates)) {
    cur <- normalizePath(start, mustWork = FALSE)
    repeat {
      if (file.exists(file.path(cur, "scripts", "00_master_pipeline.R"))) return(cur)
      parent <- dirname(cur)
      if (identical(parent, cur)) break
      cur <- parent
    }
  }
  stop("[00_master_pipeline] Cannot locate project root.")
}

PROJECT_ROOT <- find_project_root()
setwd(PROJECT_ROOT)

OBJ_DIR <- file.path(PROJECT_ROOT, "outputs", "v1", "objects")
TABLES_SUPP_DIR <- file.path(PROJECT_ROOT, "outputs", "tables", "supplementary")
dir.create(OBJ_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLES_SUPP_DIR, recursive = TRUE, showWarnings = FALSE)

resolve_col <- function(df, choices) {
  nms <- names(df)
  idx <- match(TRUE, tolower(nms) %in% tolower(choices), nomatch = 0)
  if (idx == 0) return(NA_character_)
  nms[idx]
}

DF_IDS_ID_CHOICES <- c("ind", "individual", "sample", "sampleid", "id", "ind_id")
DF_IDS_SITE_CHOICES <- c("Site", "site", "pop", "population")

normalize_id <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("\uFEFF", "", x, fixed = TRUE)
  x <- gsub("[[:cntrl:]]", "", x)
  toupper(x)
}

clean_column_names <- function(nms) {
  nms <- trimws(as.character(nms))
  nms[is.na(nms) | !nzchar(nms)] <- "unnamed_col"
  nms <- gsub("\r|\n|\t", "_", nms)
  nms <- gsub("\\s+", "_", nms)
  make.unique(nms, sep = "__dup")
}

read_csv_with_comments <- function(path) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  lines <- lines[!grepl("^\\s*#", lines)]
  lines <- lines[nzchar(trimws(lines))]
  if (length(lines) == 0) stop("[00_master_pipeline] No tabular content in ", path)
  con <- textConnection(lines)
  on.exit(close(con), add = TRUE)
  read.csv(con, stringsAsFactors = FALSE, check.names = FALSE)
}

source_dirs <- function() {
  c(
    file.path(PROJECT_ROOT, "data", "raw"),
    file.path(PROJECT_ROOT, "inputs")
  )
}

list_candidate_files <- function() {
  roots <- source_dirs()
  files <- character(0)
  for (d in roots) {
    if (!dir.exists(d)) next
    f <- list.files(d, recursive = TRUE, full.names = TRUE)
    if (length(f) == 0) next
    f <- f[file.info(f)$isdir %in% FALSE]
    f <- f[grepl("\\.(xlsx|xls|csv|tsv|txt)$", basename(f), ignore.case = TRUE)]
    files <- c(files, f)
  }
  unique(normalizePath(files, winslash = "/", mustWork = FALSE))
}

needs_readxl <- function(path) {
  grepl("\\.(xlsx|xls)$", path, ignore.case = TRUE)
}

read_table_file <- function(path, sheet = NULL) {
  ext <- tolower(sub("^.*\\.", "", path))
  if (ext %in% c("xlsx", "xls")) {
    if (!requireNamespace("readxl", quietly = TRUE)) {
      stop("[00_master_pipeline] readxl package is required to read Excel sources: ", path)
    }
    df <- as.data.frame(readxl::read_excel(path, sheet = sheet, .name_repair = "minimal"), stringsAsFactors = FALSE, check.names = FALSE)
  } else if (ext == "csv") {
    df <- read_csv_with_comments(path)
  } else if (ext %in% c("tsv", "txt")) {
    df <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    stop("Unsupported source file extension: ", ext)
  }
  if (ncol(df) == 0) return(df)
  names(df) <- clean_column_names(names(df))
  df
}

get_sheets <- function(path) {
  if (!needs_readxl(path)) return(NA_character_)
  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop("[00_master_pipeline] readxl package is required to inspect workbook sheets: ", path)
  }
  readxl::excel_sheets(path)
}

summarize_table <- function(df) {
  id_col <- resolve_col(df, c(
    DF_IDS_ID_CHOICES,
    "Nom_Labo_Échantillons", "Nom_Labo_Echantillons", "Nom_Labo_Echantillon"
  ))
  site_col <- resolve_col(df, c(
    DF_IDS_SITE_CHOICES,
    "Numéro_Population", "Numero_Population"
  ))
  lat_col <- resolve_col(df, c("latitude", "lat"))
  lon_col <- resolve_col(df, c("longitude", "lon", "long"))
  
  nms <- names(df)
  allele_1 <- grep("(_|\\.)1$", nms, perl = TRUE)
  allele_2 <- grep("(_|\\.)2$", nms, perl = TRUE)
  loci1 <- sub("(_|\\.)1$", "", nms[allele_1], perl = TRUE)
  loci2 <- sub("(_|\\.)2$", "", nms[allele_2], perl = TRUE)
  paired_loci <- sort(intersect(loci1, loci2))
  
  list(
    nrow = nrow(df),
    ncol = ncol(df),
    id_col = id_col,
    site_col = site_col,
    lat_col = lat_col,
    lon_col = lon_col,
    paired_loci = paired_loci,
    n_paired_loci = length(paired_loci)
  )
}

scan_sources <- function() {
  files <- list_candidate_files()
  if (length(files) == 0) {
    stop("[00_master_pipeline] No candidate source files found in data/raw or inputs.")
  }
  
  out <- list()
  for (path in files) {
    sheets <- get_sheets(path)
    if (length(sheets) == 1 && is.na(sheets[1])) {
      tab <- tryCatch(read_table_file(path), error = function(e) NULL)
      if (is.null(tab) || nrow(tab) == 0 || ncol(tab) == 0) next
      sm <- summarize_table(tab)
      out[[length(out) + 1]] <- list(path = path, sheet = NA_character_, df = tab, summary = sm)
    } else {
      for (sh in sheets) {
        tab <- tryCatch(read_table_file(path, sheet = sh), error = function(e) NULL)
        if (is.null(tab) || nrow(tab) == 0 || ncol(tab) == 0) next
        sm <- summarize_table(tab)
        out[[length(out) + 1]] <- list(path = path, sheet = sh, df = tab, summary = sm)
      }
    }
  }
  out
}

select_meta_ind <- function(scanned) {
  cands <- Filter(function(x) {
    !is.na(x$summary$id_col) && !is.na(x$summary$site_col) && x$summary$n_paired_loci < 3
  }, scanned)
  
  if (length(cands) == 0) {
    stop("[00_master_pipeline] Could not find an individual metadata table (ID + Site) in workbook sources.")
  }
  
  score <- function(x) {
    b <- tolower(basename(x$path))
    s <- 0
    if (needs_readxl(x$path)) s <- s + 50
    if (grepl("meta|field|terrain|sample|ind", b)) s <- s + 20
    if (grepl("geno|poppr|microsat", b)) s <- s - 10
    if (grepl("/inputs/", x$path)) s <- s + 5
    s + x$summary$nrow / 10000
  }
  
  best <- cands[[which.max(vapply(cands, score, numeric(1)))]]
  df <- best$df
  
  id_col <- best$summary$id_col
  site_col <- best$summary$site_col
  out <- data.frame(
    ind_id = trimws(as.character(df[[id_col]])),
    Site = trimws(as.character(df[[site_col]])),
    stringsAsFactors = FALSE
  )
  out <- out[nzchar(out$ind_id), , drop = FALSE]
  
  dup <- out$ind_id[duplicated(out$ind_id)]
  if (length(dup) > 0) {
    stop("[00_master_pipeline] Duplicate IDs in selected metadata table: ",
         paste(head(unique(dup), 30), collapse = ", "))
  }
  
  if (any(!nzchar(out$Site) | is.na(out$Site))) {
    stop("[00_master_pipeline] Selected metadata table has missing Site labels.")
  }
  
  message("[00_master_pipeline] Chosen individual metadata source: ", best$path)
  message("[00_master_pipeline] Chosen individual metadata sheet: ", ifelse(is.na(best$sheet), "<file>", best$sheet))
  message("[00_master_pipeline] Chosen individual metadata ID column: ", id_col)
  message("[00_master_pipeline] Chosen individual metadata Site column: ", site_col)
  message("[00_master_pipeline] Imported metadata dimensions: ", nrow(df), " x ", ncol(df))
  
  list(data = out, source = best)
}

select_site_meta <- function(scanned, meta_ind_tbl) {
  cands <- Filter(function(x) {
    !is.na(x$summary$site_col) &&
      (x$summary$n_paired_loci < 2) &&
      (x$summary$nrow <= 2000)
  }, scanned)
  
  if (length(cands) == 0) {
    fallback <- unique(meta_ind_tbl[, "Site", drop = FALSE])
    fallback <- data.frame(Site = fallback$Site, stringsAsFactors = FALSE)
    message("[00_master_pipeline] No explicit site metadata table found; deriving site table from individual metadata Site labels.")
    return(fallback)
  }
  
  score <- function(x) {
    b <- tolower(basename(x$path))
    s <- 0
    if (needs_readxl(x$path)) s <- s + 10
    if (grepl("site|metadata|field|coord|location", b)) s <- s + 15
    if (!is.na(x$summary$lat_col)) s <- s + 20
    if (!is.na(x$summary$lon_col)) s <- s + 20
    if (!is.na(x$summary$id_col)) s <- s - 5
    s - x$summary$nrow / 100000
  }
  
  best <- cands[[which.max(vapply(cands, score, numeric(1)))]]
  df <- best$df
  site_col <- best$summary$site_col
  df[[site_col]] <- trimws(as.character(df[[site_col]]))
  df <- df[nzchar(df[[site_col]]), , drop = FALSE]
  names(df)[names(df) == site_col] <- "Site"
  
  message("[00_master_pipeline] Chosen site metadata source: ", best$path)
  message("[00_master_pipeline] Chosen site metadata sheet: ", ifelse(is.na(best$sheet), "<file>", best$sheet))
  message("[00_master_pipeline] Chosen site metadata Site column: ", site_col)
  message("[00_master_pipeline] Imported site metadata dimensions: ", nrow(df), " x ", ncol(df))
  
  df
}

select_genotype_source <- function(scanned, canonical_ids_norm) {
  cands <- Filter(function(x) {
    !is.na(x$summary$id_col) && x$summary$n_paired_loci >= 2
  }, scanned)
  
  if (length(cands) == 0) {
    stop("[00_master_pipeline] Could not find genotype table with paired allele columns in workbook sources.")
  }
  
  score <- function(x) {
    ids <- normalize_id(x$df[[x$summary$id_col]])
    ids <- ids[nzchar(ids)]
    match_n <- sum(ids %in% canonical_ids_norm)
    b <- tolower(basename(x$path))
    s <- 0
    if (needs_readxl(x$path)) s <- s + 20
    if (grepl("poppr|geno|genotype|microsat", b)) s <- s + 10
    match_n * 1000 + x$summary$n_paired_loci * 10 + s
  }
  
  best <- cands[[which.max(vapply(cands, score, numeric(1)))]]
  
  message("[00_master_pipeline] Chosen genotype source file: ", best$path)
  message("[00_master_pipeline] Chosen genotype source sheet: ", ifelse(is.na(best$sheet), "<file>", best$sheet))
  message("[00_master_pipeline] Chosen genotype ID column: ", best$summary$id_col)
  message("[00_master_pipeline] Chosen genotype Site column: ", ifelse(is.na(best$summary$site_col), "<not found>", best$summary$site_col))
  message("[00_master_pipeline] Imported genotype dimensions: ", best$summary$nrow, " x ", best$summary$ncol)
  
  best
}

build_genind_from_table <- function(tbl, table_summary, canonical_ids_norm, id_to_site) {
  id_col <- table_summary$id_col
  pop_col <- table_summary$site_col
  paired_loci <- table_summary$paired_loci
  
  ids_raw <- trimws(as.character(tbl[[id_col]]))
  ids_norm <- normalize_id(ids_raw)
  keep_nonempty <- nzchar(ids_norm)
  
  tbl <- tbl[keep_nonempty, , drop = FALSE]
  ids_raw <- ids_raw[keep_nonempty]
  ids_norm <- ids_norm[keep_nonempty]
  
  in_canonical <- ids_norm %in% canonical_ids_norm
  
  n_raw_ids <- length(ids_norm)
  n_unique_raw_ids <- length(unique(ids_norm))
  n_matched <- sum(in_canonical)
  n_excluded <- sum(!in_canonical)
  excluded_examples <- unique(ids_raw[!in_canonical])
  
  tbl <- tbl[in_canonical, , drop = FALSE]
  ids_raw <- ids_raw[in_canonical]
  ids_norm <- ids_norm[in_canonical]
  
  if (nrow(tbl) == 0) {
    stop("[00_master_pipeline] No genotype rows matched canonical metadata IDs.")
  }
  
  paired_loci_n <- length(paired_loci)
  if (paired_loci_n == 0) {
    stop("[00_master_pipeline] No paired microsatellite loci were available to build the genind object.")
  }
  
  dup <- duplicated(ids_norm)
  n_dups <- sum(dup)
  if (n_dups > 0) {
    message("[00_master_pipeline] Duplicate genotype rows for canonical IDs: ", n_dups,
            " (keeping first row per normalized ID).")
    keep_first <- !duplicated(ids_norm)
    tbl <- tbl[keep_first, , drop = FALSE]
    ids_raw <- ids_raw[keep_first]
    ids_norm <- ids_norm[keep_first]
  }
  
  missingness_tbl <- calculate_individual_missingness(
    tbl = tbl,
    table_summary = table_summary,
    ids_raw = ids_raw
  )
  
  locus_pattern <- function(loc, suffix) {
    paste0("^", gsub("([\\W])", "\\\\\\1", loc, perl = TRUE), "(_|\\\\.)", suffix, "$")
  }
  
  geno <- data.frame(row.names = ids_raw, stringsAsFactors = FALSE)
  nms <- names(tbl)
  for (loc in paired_loci) {
    c1 <- nms[grep(locus_pattern(loc, "1"), nms, perl = TRUE)][1]
    c2 <- nms[grep(locus_pattern(loc, "2"), nms, perl = TRUE)][1]
    
    a1 <- trimws(as.character(tbl[[c1]]))
    a2 <- trimws(as.character(tbl[[c2]]))
    miss <- (!nzchar(a1) | toupper(a1) == "NA") | (!nzchar(a2) | toupper(a2) == "NA")
    g <- paste(a1, a2, sep = "/")
    g[miss] <- NA_character_
    geno[[loc]] <- g
  }
  
  site_vec <- id_to_site[ids_norm]
  if (any(is.na(site_vec) | !nzchar(site_vec))) {
    stop("[00_master_pipeline] Missing Site mapping for one or more canonical genotype IDs.")
  }
  
  gi <- adegenet::df2genind(
    geno,
    ploidy = 2,
    ind.names = rownames(geno),
    pop = as.factor(site_vec),
    type = "codom",
    sep = "/",
    ncode = 3
  )
  
  list(
    gi = gi,
    missingness = missingness_tbl,
    id_col = id_col,
    pop_col = pop_col,
    diagnostics = list(
      n_raw_ids = n_raw_ids,
      n_unique_raw_ids = n_unique_raw_ids,
      n_matched = n_matched,
      n_excluded = n_excluded,
      excluded_examples = head(excluded_examples, 20)
    )
  )
}

calculate_individual_missingness <- function(tbl, table_summary, ids_raw) {
  paired_loci <- table_summary$paired_loci
  if (length(paired_loci) == 0) {
    stop("[00_master_pipeline] Cannot calculate individual missingness: no paired loci were detected.")
  }
  
  n_ind <- nrow(tbl)
  total_alleles_expected <- 2L * length(paired_loci)
  missing_allele_count <- integer(n_ind)
  missing_locus_count <- integer(n_ind)
  nms <- names(tbl)
  
  locus_pattern <- function(loc, suffix) {
    paste0("^", gsub("([\\W])", "\\\\\\1", loc, perl = TRUE), "(_|\\\\.)", suffix, "$")
  }
  
  is_missing_allele <- function(x) {
    x_chr <- trimws(as.character(x))
    is.na(x) | !nzchar(x_chr) | tolower(x_chr) %in% c("na", "n/a", "null", ".", "-", "?")
  }
  
  for (loc in paired_loci) {
    c1 <- nms[grep(locus_pattern(loc, "1"), nms, perl = TRUE)][1]
    c2 <- nms[grep(locus_pattern(loc, "2"), nms, perl = TRUE)][1]
    if (is.na(c1) || is.na(c2)) {
      stop("[00_master_pipeline] Could not resolve allele columns for locus: ", loc)
    }
    m1 <- is_missing_allele(tbl[[c1]])
    m2 <- is_missing_allele(tbl[[c2]])
    missing_allele_count <- missing_allele_count + as.integer(m1) + as.integer(m2)
    missing_locus_count <- missing_locus_count + as.integer(m1 | m2)
  }
  
  data.frame(
    ind_id = ids_raw,
    missing_allele_count = missing_allele_count,
    total_alleles_expected = total_alleles_expected,
    percent_missing = missing_allele_count / total_alleles_expected,
    missing_locus_count = missing_locus_count,
    total_loci = length(paired_loci),
    stringsAsFactors = FALSE
  )
}

apply_missing_data_filter <- function(gi, missingness_tbl, threshold = MISSING_ALLELE_FILTER_THRESHOLD, output_dir = TABLES_SUPP_DIR) {
  if (!inherits(gi, "genind")) {
    stop("[00_master_pipeline] Missing-data filter expects a genind object.")
  }
  if (!is.data.frame(missingness_tbl) || nrow(missingness_tbl) != adegenet::nInd(gi)) {
    stop("[00_master_pipeline] Missingness table must be a data.frame aligned to all individuals in gi.")
  }
  required_cols <- c("ind_id", "percent_missing", "missing_allele_count", "total_alleles_expected")
  if (!all(required_cols %in% names(missingness_tbl))) {
    stop("[00_master_pipeline] Missingness table is missing required columns: ",
         paste(setdiff(required_cols, names(missingness_tbl)), collapse = ", "))
  }
  
  gi_ids <- adegenet::indNames(gi)
  idx <- match(normalize_id(gi_ids), normalize_id(missingness_tbl$ind_id))
  if (any(is.na(idx))) {
    stop("[00_master_pipeline] Could not align missingness table to gi. Example missing IDs: ",
         paste(head(gi_ids[is.na(idx)], 10), collapse = ", "))
  }
  missingness_tbl <- missingness_tbl[idx, , drop = FALSE]
  if (!all(normalize_id(gi_ids) == normalize_id(missingness_tbl$ind_id))) {
    stop("[00_master_pipeline] Missingness table order does not match gi individual order after alignment.")
  }
  if (any(!is.finite(missingness_tbl$percent_missing))) {
    stop("[00_master_pipeline] Non-finite per-individual missing-data proportions detected.")
  }
  
  remove_idx <- missingness_tbl$percent_missing > threshold
  retained_n <- sum(!remove_idx)
  removed_n <- sum(remove_idx)
  
  if (retained_n == 0) {
    stop(
      "[00_master_pipeline] Missing-data filter would remove all individuals at threshold > ",
      sprintf("%.1f%%", threshold * 100),
      ". Check genotype import and missing-data coding before continuing."
    )
  }
  
  removed_tbl <- data.frame(
    individual_id = missingness_tbl$ind_id[remove_idx],
    site = as.character(adegenet::pop(gi))[remove_idx],
    percent_missing = round(missingness_tbl$percent_missing[remove_idx] * 100, 2),
    removed_threshold = threshold * 100,
    stringsAsFactors = FALSE
  )
  removed_tbl <- removed_tbl[order(-removed_tbl$percent_missing, removed_tbl$site, removed_tbl$individual_id), , drop = FALSE]
  
  removed_file <- file.path(output_dir, "individuals_removed_missing_gt35pct_allele_data.csv")
  write.csv(removed_tbl, removed_file, row.names = FALSE)
  
  cat("[00_master_pipeline] Missing-data filter (>35% missing allele data; applied before clonality):\n")
  cat("[00_master_pipeline] Individuals before filtering: ", adegenet::nInd(gi), "\n", sep = "")
  cat("[00_master_pipeline] Individuals removed: ", removed_n, "\n", sep = "")
  cat("[00_master_pipeline] Individuals retained: ", retained_n, "\n", sep = "")
  if (removed_n > 0) {
    removed_lines <- paste0(
      removed_tbl$individual_id,
      " (Site=", removed_tbl$site,
      ", Missing=", sprintf("%.2f", removed_tbl$percent_missing), "%)"
    )
    cat("[00_master_pipeline] IDs removed: ", paste(removed_tbl$individual_id, collapse = ", "), "\n", sep = "")
    cat("[00_master_pipeline] Removed individual missingness: ", paste(removed_lines, collapse = "; "), "\n", sep = "")
  } else {
    cat("[00_master_pipeline] IDs removed: none\n")
  }
  cat("[00_master_pipeline] Removed-individual table saved to: ", removed_file, "\n", sep = "")
  
  gi_filtered <- gi[!remove_idx, , drop = TRUE]
  if (adegenet::nInd(gi_filtered) != retained_n) {
    stop("[00_master_pipeline] Filtered gi size does not match retained-individual count.")
  }
  if (!all(adegenet::indNames(gi_filtered) == missingness_tbl$ind_id[!remove_idx])) {
    stop("[00_master_pipeline] Filtered gi IDs are inconsistent with the retained missingness table.")
  }
  
  retained_tbl <- missingness_tbl[!remove_idx, , drop = FALSE]
  rownames(retained_tbl) <- NULL
  
  list(
    gi = gi_filtered,
    retained_missingness = retained_tbl,
    removed_missingness = removed_tbl,
    removed_ids = removed_tbl$individual_id,
    threshold = threshold,
    output_file = removed_file,
    n_before = adegenet::nInd(gi),
    n_removed = removed_n,
    n_retained = retained_n
  )
}

build_mll_clone_corrected_object <- function(gi) {
  gc_mlg <- poppr::as.genclone(gi)
  mlg_raw <- tryCatch(poppr::mlg.vector(gc_mlg), error = function(e) as.integer(factor(poppr::mlg(gc_mlg))))
  mlg_labels <- paste0("MLG_", as.integer(factor(mlg_raw)))
  
  replen <- rep(2, adegenet::nLoc(gi))
  names(replen) <- adegenet::locNames(gi)
  
  gc_mll <- gc_mlg
  poppr::mlg.filter(
    gc_mll,
    distance = poppr::bruvo.dist,
    replen = replen,
    algorithm = BRUVO_ALGORITHM
  ) <- BRUVO_MLL_THRESHOLD
  
  mll_raw <- poppr::mll(gc_mll)
  mll_labels <- paste0("MLL_", as.integer(factor(mll_raw)))
  
  if (length(mll_labels) != adegenet::nInd(gi)) {
    stop("[00_master_pipeline] MLL assignment length does not match nInd(gi).")
  }
  
  keep_mll <- !duplicated(mll_labels)
  gi_mll <- gi[keep_mll, ]
  
  if (adegenet::nInd(gi_mll) != length(unique(mll_labels))) {
    stop("[00_master_pipeline] gi_mll size is inconsistent with unique MLL count.")
  }
  
  list(
    gi_mll = gi_mll,
    mlg_labels = mlg_labels,
    mll_labels = mll_labels,
    threshold = BRUVO_MLL_THRESHOLD,
    algorithm = BRUVO_ALGORITHM,
    n_mlg = length(unique(mlg_labels)),
    n_mll = length(unique(mll_labels)),
    n_clonal_repeats = adegenet::nInd(gi) - length(unique(mll_labels))
  )
}

is_missing_allele_value <- function(x) {
  x_chr <- trimws(as.character(x))
  is.na(x) | !nzchar(x_chr) | tolower(x_chr) %in% c("na", "n/a", "null", ".", "-", "?", "0", "-9")
}

format_allele_to_3digits <- function(x) {
  x_chr <- trimws(as.character(x))
  x_chr <- gsub(",", "", x_chr, fixed = TRUE)
  if (!nzchar(x_chr)) {
    return(NA_character_)
  }
  if (!grepl("^-?[0-9]+(?:\\.0+)?$", x_chr)) {
    stop("[00_master_pipeline] Non-numeric allele value encountered during Micro-Checker Genepop export: '", x_chr, "'.")
  }
  x_num <- suppressWarnings(as.numeric(x_chr))
  if (!is.finite(x_num) || x_num < 0) {
    stop("[00_master_pipeline] Invalid allele value encountered during Micro-Checker Genepop export: '", x_chr, "'.")
  }
  sprintf("%03d", as.integer(round(x_num)))
}

detect_allele_pairs <- function(tbl, id_col, pop_col) {
  protected_cols <- c(id_col, pop_col)
  protected_cols <- protected_cols[!is.na(protected_cols)]
  
  candidate_cols <- setdiff(names(tbl), protected_cols)
  allele_cols <- candidate_cols[grepl("(_|\\.)[12]$", candidate_cols, perl = TRUE)]
  
  if (length(allele_cols) == 0) {
    stop("[00_master_pipeline] Could not detect paired allele columns ending in '_1/_2' or '.1/.2'.")
  }
  if ((length(allele_cols) %% 2) != 0) {
    stop("[00_master_pipeline] Detected an odd number of allele columns (", length(allele_cols), "). Each locus must have exactly two allele columns.")
  }
  
  loci <- sub("(_|\\.)[12]$", "", allele_cols, perl = TRUE)
  loci_tbl <- split(allele_cols, loci)
  
  missing_pairs <- names(loci_tbl)[vapply(loci_tbl, function(cols) {
    has1 <- any(grepl("(_|\\.)1$", cols, perl = TRUE))
    has2 <- any(grepl("(_|\\.)2$", cols, perl = TRUE))
    !(has1 && has2)
  }, logical(1))]
  
  if (length(missing_pairs) > 0) {
    stop("[00_master_pipeline] Some loci are missing one of the required allele columns (_1/_2 or .1/.2): ",
         paste(head(missing_pairs, 20), collapse = ", "))
  }
  
  ordered_loci <- sort(unique(loci))
  pair_tbl <- data.frame(
    locus = ordered_loci,
    a1_col = vapply(ordered_loci, function(loc) {
      cols <- loci_tbl[[loc]]
      cols[grepl("(_|\\.)1$", cols, perl = TRUE)][1]
    }, character(1)),
    a2_col = vapply(ordered_loci, function(loc) {
      cols <- loci_tbl[[loc]]
      cols[grepl("(_|\\.)2$", cols, perl = TRUE)][1]
    }, character(1)),
    stringsAsFactors = FALSE
  )
  
  pair_tbl
}

write_microchecker_genepop_export <- function(tbl,
                                              id_col = NULL,
                                              pop_col = NULL,
                                              allowed_ids_norm = NULL,
                                              id_to_pop = NULL,
                                              output_path = file.path(PROJECT_ROOT, "data", "derived", "microchecker_genepop.txt"),
                                              title_line = "BeechCode MicroChecker export") {
  if (!is.data.frame(tbl) || nrow(tbl) == 0 || ncol(tbl) == 0) {
    stop("[00_master_pipeline] Micro-Checker Genepop export expects a non-empty data.frame.")
  }
  
  if (is.null(id_col) || is.na(id_col) || !id_col %in% names(tbl)) {
    id_col <- resolve_col(tbl, c(
      DF_IDS_ID_CHOICES,
      "Nom_Labo_Échantillons", "Nom_Labo_Echantillons", "Nom_Labo_Echantillon"
    ))
  }
  if (is.na(id_col) || !nzchar(id_col)) {
    stop("[00_master_pipeline] Could not detect sample ID column for Micro-Checker Genepop export.")
  }
  
  if (is.null(pop_col) || is.na(pop_col) || !pop_col %in% names(tbl)) {
    pop_col <- resolve_col(tbl, c(
      DF_IDS_SITE_CHOICES,
      "Numéro_Population", "Numero_Population"
    ))
  }
  
  allele_pairs <- detect_allele_pairs(tbl, id_col = id_col, pop_col = pop_col)
  
  work <- tbl
  work$.__id_raw <- trimws(as.character(work[[id_col]]))
  work$.__id_norm <- normalize_id(work$.__id_raw)
  work <- work[nzchar(work$.__id_norm), , drop = FALSE]
  
  if (!is.null(allowed_ids_norm)) {
    allowed_ids_norm <- unique(normalize_id(allowed_ids_norm))
    work <- work[work$.__id_norm %in% allowed_ids_norm, , drop = FALSE]
  }
  
  if (nrow(work) == 0) {
    stop("[00_master_pipeline] No rows available for Micro-Checker Genepop export after filtering IDs.")
  }
  
  if (!is.null(id_to_pop)) {
    id_to_pop <- as.character(id_to_pop)
    names(id_to_pop) <- normalize_id(names(id_to_pop))
  }
  
  if (!is.null(pop_col) && pop_col %in% names(work)) {
    pop_vals <- trimws(as.character(work[[pop_col]]))
  } else {
    pop_vals <- rep(NA_character_, nrow(work))
  }
  
  if (!is.null(id_to_pop)) {
    mapped <- unname(id_to_pop[work$.__id_norm])
    need_fill <- is.na(pop_vals) | !nzchar(pop_vals)
    pop_vals[need_fill] <- mapped[need_fill]
  }
  
  if (any(is.na(pop_vals) | !nzchar(pop_vals))) {
    bad_ids <- unique(work$.__id_raw[is.na(pop_vals) | !nzchar(pop_vals)])
    stop("[00_master_pipeline] Missing population/site assignment for one or more individuals in Genepop export. Examples: ",
         paste(head(bad_ids, 20), collapse = ", "))
  }
  work$.__pop <- pop_vals
  
  dup_rows <- duplicated(work$.__id_norm)
  if (any(dup_rows)) {
    message("[00_master_pipeline] Duplicate IDs detected in genotype source for Genepop export; keeping first occurrence per normalized ID.")
    work <- work[!dup_rows, , drop = FALSE]
  }
  
  locus_names <- allele_pairs$locus
  if (length(locus_names) == 0) {
    stop("[00_master_pipeline] No loci detected for Micro-Checker Genepop export.")
  }
  
  partial_missing_counter <- integer(length(locus_names))
  genotype_by_locus <- vector("list", length(locus_names))
  
  for (i in seq_along(locus_names)) {
    a1_raw <- work[[allele_pairs$a1_col[i]]]
    a2_raw <- work[[allele_pairs$a2_col[i]]]
    
    m1 <- is_missing_allele_value(a1_raw)
    m2 <- is_missing_allele_value(a2_raw)
    partial_missing <- xor(m1, m2)
    partial_missing_counter[i] <- sum(partial_missing)
    
    a1_fmt <- rep("000", nrow(work))
    a2_fmt <- rep("000", nrow(work))
    
    idx1 <- which(!m1)
    idx2 <- which(!m2)
    
    if (length(idx1) > 0) {
      a1_fmt[idx1] <- vapply(a1_raw[idx1], format_allele_to_3digits, character(1))
    }
    if (length(idx2) > 0) {
      a2_fmt[idx2] <- vapply(a2_raw[idx2], format_allele_to_3digits, character(1))
    }
    
    code <- ifelse(m1 & m2, "000000", paste0(a1_fmt, a2_fmt))
    genotype_by_locus[[i]] <- code
  }
  
  total_partial_missing <- sum(partial_missing_counter)
  if (total_partial_missing > 0) {
    bad_loci <- locus_names[partial_missing_counter > 0]
    warning(
      "[00_master_pipeline] Found ", total_partial_missing,
      " locus calls with one allele missing and the other present. ",
      "These were exported with the missing allele coded as 000. Loci affected: ",
      paste(head(bad_loci, 20), collapse = ", "),
      call. = FALSE
    )
  }
  
  genotype_matrix <- do.call(cbind, genotype_by_locus)
  colnames(genotype_matrix) <- locus_names
  
  out_tbl <- data.frame(
    ind_id = work$.__id_raw,
    pop = work$.__pop,
    stringsAsFactors = FALSE
  )
  out_tbl <- cbind(out_tbl, as.data.frame(genotype_matrix, stringsAsFactors = FALSE))
  out_tbl <- out_tbl[order(out_tbl$pop, out_tbl$ind_id), , drop = FALSE]
  
  genepop_lines <- c(title_line, locus_names)
  for (p in unique(out_tbl$pop)) {
    genepop_lines <- c(genepop_lines, "Pop")
    sub <- out_tbl[out_tbl$pop == p, , drop = FALSE]
    for (r in seq_len(nrow(sub))) {
      geno_str <- paste(sub[r, locus_names, drop = TRUE], collapse = " ")
      genepop_lines <- c(genepop_lines, paste0(sub$ind_id[r], " , ", geno_str))
    }
  }
  
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  writeLines(genepop_lines, con = output_path, useBytes = TRUE)
  
  if (!file.exists(output_path)) {
    stop("[00_master_pipeline] Failed to create Genepop output file: ", output_path)
  }
  
  message("[00_master_pipeline] Micro-Checker Genepop export complete.")
  message("[00_master_pipeline] Samples: ", nrow(out_tbl))
  message("[00_master_pipeline] Populations: ", length(unique(out_tbl$pop)))
  message("[00_master_pipeline] Loci: ", length(locus_names))
  message("[00_master_pipeline] Output: ", output_path)
  
  invisible(list(
    n_samples = nrow(out_tbl),
    n_populations = length(unique(out_tbl$pop)),
    n_loci = length(locus_names),
    output_path = output_path,
    locus_names = locus_names
  ))
}

build_objects <- function() {
  scanned <- scan_sources()
  
  meta_ind_info <- select_meta_ind(scanned)
  meta_ind <- meta_ind_info$data
  meta <- select_site_meta(scanned, meta_ind)
  
  canonical_ids <- meta_ind$ind_id
  canonical_ids_norm <- normalize_id(canonical_ids)
  id_to_site <- setNames(meta_ind$Site, canonical_ids_norm)
  
  message("[00_master_pipeline] Canonical sample-universe policy: IDs from selected metadata workbook table")
  message("[00_master_pipeline] Intended canonical dataset size = ", length(canonical_ids))
  
  geno_source <- select_genotype_source(scanned, canonical_ids_norm)
  built <- build_genind_from_table(
    tbl = geno_source$df,
    table_summary = geno_source$summary,
    canonical_ids_norm = canonical_ids_norm,
    id_to_site = id_to_site
  )
  gi_unfiltered <- built$gi
  gi_unfiltered_ids <- adegenet::indNames(gi_unfiltered)
  gi_unfiltered_ids_norm <- normalize_id(gi_unfiltered_ids)
  missing_in_gi_before_filter <- canonical_ids[!(canonical_ids_norm %in% gi_unfiltered_ids_norm)]
  
  # Publication-defensible QC step:
  # calculate per-individual missingness from the original paired allele
  # columns (2 alleles per locus), then remove only those individuals whose
  # missing allele proportion is strictly greater than 35% before clonality.
  missing_filter <- apply_missing_data_filter(
    gi = gi_unfiltered,
    missingness_tbl = built$missingness,
    threshold = MISSING_ALLELE_FILTER_THRESHOLD,
    output_dir = TABLES_SUPP_DIR
  )
  gi <- missing_filter$gi
  
  mll_build <- build_mll_clone_corrected_object(gi)
  gi_mll <- mll_build$gi_mll
  
  gi_ids <- adegenet::indNames(gi)
  gi_ids_norm <- normalize_id(gi_ids)
  missing_in_gi_after_filter <- canonical_ids[!(canonical_ids_norm %in% gi_ids_norm)]
  retained_missingness <- missing_filter$retained_missingness
  match_idx <- match(gi_ids_norm, normalize_id(retained_missingness$ind_id))
  if (any(is.na(match_idx))) {
    stop("[00_master_pipeline] Failed to align retained missingness metrics to filtered gi IDs.")
  }
  retained_missingness <- retained_missingness[match_idx, , drop = FALSE]
  if (!all(normalize_id(retained_missingness$ind_id) == gi_ids_norm)) {
    stop("[00_master_pipeline] Retained missingness metrics are not in the same order as filtered gi IDs.")
  }
  
  df_ids <- data.frame(
    ind_id = gi_ids,
    Site = as.character(adegenet::pop(gi)),
    missing_allele_count = retained_missingness$missing_allele_count,
    total_alleles_expected = retained_missingness$total_alleles_expected,
    percent_missing_allele_data = retained_missingness$percent_missing,
    missing_data_filter_threshold = missing_filter$threshold,
    MLG = mll_build$mlg_labels,
    MLL = mll_build$mll_labels,
    Bruvo_MLL_threshold = mll_build$threshold,
    Bruvo_algorithm = mll_build$algorithm,
    stringsAsFactors = FALSE
  )
  
  # Export filtered full dataset as a true Genepop plain-text file for
  # Micro-Checker null-allele diagnostics.
  write_microchecker_genepop_export(
    tbl = geno_source$df,
    id_col = geno_source$summary$id_col,
    pop_col = geno_source$summary$site_col,
    allowed_ids_norm = normalize_id(adegenet::indNames(gi)),
    id_to_pop = id_to_site,
    output_path = file.path(PROJECT_ROOT, "data", "derived", "microchecker_genepop.txt"),
    title_line = "BeechCode MicroChecker export"
  )
  
  message("[00_master_pipeline] Number of raw genotype IDs: ", built$diagnostics$n_raw_ids,
          " (unique normalized: ", built$diagnostics$n_unique_raw_ids, ")")
  message("[00_master_pipeline] Number of metadata IDs: ", length(canonical_ids))
  message("[00_master_pipeline] Number matched before missing-data filtering: ", built$diagnostics$n_matched)
  message("[00_master_pipeline] Number excluded (outside canonical metadata universe): ", built$diagnostics$n_excluded)
  if (built$diagnostics$n_excluded > 0) {
    message("[00_master_pipeline] Excluded genotype ID examples (up to 20): ",
            paste(built$diagnostics$excluded_examples, collapse = ", "))
  }
  message("[00_master_pipeline] Canonical metadata IDs missing genotype rows before missing-data filtering: ", length(missing_in_gi_before_filter))
  if (length(missing_in_gi_before_filter) > 0) {
    message("[00_master_pipeline] Missing genotype ID examples before filtering (up to 20): ",
            paste(head(missing_in_gi_before_filter, 20), collapse = ", "))
  }
  message("[00_master_pipeline] Individuals removed for >", MISSING_ALLELE_FILTER_THRESHOLD * 100,
          "% missing allele data: ", missing_filter$n_removed)
  message("[00_master_pipeline] Canonical metadata IDs absent from final filtered gi: ", length(missing_in_gi_after_filter))
  if (length(missing_in_gi_after_filter) > 0) {
    message("[00_master_pipeline] Final absent ID examples (up to 20): ",
            paste(head(missing_in_gi_after_filter, 20), collapse = ", "))
  }
  cat("[00_master_pipeline] Bruvo MLL threshold: ", BRUVO_MLL_THRESHOLD, "\n", sep = "")
  cat("[00_master_pipeline] Bruvo clustering algorithm: ", BRUVO_ALGORITHM, "\n", sep = "")
  cat("[00_master_pipeline] Number of unique MLGs: ", mll_build$n_mlg, "\n", sep = "")
  cat("[00_master_pipeline] Number of unique MLLs: ", mll_build$n_mll, "\n", sep = "")
  cat("[00_master_pipeline] Number of clones (nInd - unique MLL): ", mll_build$n_clonal_repeats, "\n", sep = "")
  cat("[00_master_pipeline] Confirmed df_ids contains columns: ind_id, Site, missing_allele_count, total_alleles_expected, percent_missing_allele_data, MLG, MLL\n", sep = "")
  
  saveRDS(gi, file.path(OBJ_DIR, "gi.rds"))
  saveRDS(gi_mll, file.path(OBJ_DIR, "gi_mll.rds"))
  saveRDS(df_ids, file.path(OBJ_DIR, "df_ids.rds"))
  saveRDS(meta, file.path(OBJ_DIR, "meta.rds"))
  
  message("[00_master_pipeline] Final nInd(gi) = ", adegenet::nInd(gi))
  message("[00_master_pipeline] Final nInd(gi_mll) = ", adegenet::nInd(gi_mll))
  message("[00_master_pipeline] Final nrow(df_ids) = ", nrow(df_ids))
  message("[00_master_pipeline] Final nrow(meta) = ", nrow(meta))
  message("[00_master_pipeline] Saved objects:")
  message("  - ", file.path(OBJ_DIR, "gi.rds"))
  message("  - ", file.path(OBJ_DIR, "gi_mll.rds"))
  message("  - ", file.path(OBJ_DIR, "df_ids.rds"))
  message("  - ", file.path(OBJ_DIR, "meta.rds"))
}

build_objects()