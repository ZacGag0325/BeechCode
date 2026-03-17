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
############################################################

suppressPackageStartupMessages({
  library(adegenet)
  library(poppr)
})

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
dir.create(OBJ_DIR, recursive = TRUE, showWarnings = FALSE)

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
    df <- as.data.frame(readxl::read_excel(path, sheet = sheet), stringsAsFactors = FALSE, check.names = FALSE)
  } else if (ext == "csv") {
    df <- read_csv_with_comments(path)
  } else if (ext %in% c("tsv", "txt")) {
    df <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    stop("Unsupported source file extension: ", ext)
  }
  if (ncol(df) == 0) return(df)
  names(df) <- trimws(names(df))
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
  gi <- built$gi
  
  mlg_raw <- tryCatch(poppr::mlg.vector(gi), error = function(e) as.integer(factor(poppr::mlg(gi))))
  mll_id <- as.integer(factor(mlg_raw))
  mll_label <- paste0("MLL", mll_id)
  keep_mll <- !duplicated(mll_label)
  gi_mll <- gi[keep_mll, ]
  
  gi_ids <- adegenet::indNames(gi)
  gi_ids_norm <- normalize_id(gi_ids)
  missing_in_gi <- canonical_ids[!(canonical_ids_norm %in% gi_ids_norm)]
  
  df_ids <- data.frame(
    ind_id = gi_ids,
    Site = as.character(adegenet::pop(gi)),
    mll = mll_label,
    stringsAsFactors = FALSE
  )
  
  message("[00_master_pipeline] Number of raw genotype IDs: ", built$diagnostics$n_raw_ids,
          " (unique normalized: ", built$diagnostics$n_unique_raw_ids, ")")
  message("[00_master_pipeline] Number of metadata IDs: ", length(canonical_ids))
  message("[00_master_pipeline] Number matched: ", built$diagnostics$n_matched)
  message("[00_master_pipeline] Number excluded (outside canonical metadata universe): ", built$diagnostics$n_excluded)
  if (built$diagnostics$n_excluded > 0) {
    message("[00_master_pipeline] Excluded genotype ID examples (up to 20): ",
            paste(built$diagnostics$excluded_examples, collapse = ", "))
  }
  message("[00_master_pipeline] Canonical metadata IDs missing genotype rows: ", length(missing_in_gi))
  if (length(missing_in_gi) > 0) {
    message("[00_master_pipeline] Missing genotype ID examples (up to 20): ",
            paste(head(missing_in_gi, 20), collapse = ", "))
  }
  
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