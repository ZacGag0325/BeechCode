# scripts/00_master_pipeline.R
############################################################
# Build canonical objects for downstream analyses:
# - outputs/v1/objects/gi.rds
# - outputs/v1/objects/gi_mll.rds
# - outputs/v1/objects/df_ids.rds
# - outputs/v1/objects/meta.rds
#
# Canonical individual universe policy:
# - Use curated metadata universe from inputs/meta_ind.csv as canonical.
# - In this project that curated universe is expected to be 278 individuals.
# - Raw genotype rows not in curated metadata are reported and excluded.
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

normalize_id <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("^\\ufeff", "", x, perl = TRUE)
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

load_meta_ind <- function() {
  path <- file.path(PROJECT_ROOT, "inputs", "meta_ind.csv")
  if (!file.exists(path)) stop("[00_master_pipeline] Missing required file: ", path)
  x <- read_csv_with_comments(path)
  id_col <- resolve_col(x, c("ind", "individual", "sample", "sampleid", "id", "ind_id"))
  site_col <- resolve_col(x, c("Site", "site", "pop", "population"))
  if (is.na(id_col) || is.na(site_col)) {
    stop("[00_master_pipeline] inputs/meta_ind.csv must contain ID and Site columns.")
  }
  
  out <- data.frame(
    ind_id = trimws(as.character(x[[id_col]])),
    Site = trimws(as.character(x[[site_col]])),
    stringsAsFactors = FALSE
  )
  out <- out[nzchar(out$ind_id), , drop = FALSE]
  
  dup <- out$ind_id[duplicated(out$ind_id)]
  if (length(dup) > 0) {
    stop("[00_master_pipeline] Duplicate IDs in inputs/meta_ind.csv: ",
         paste(head(unique(dup), 30), collapse = ", "))
  }
  if (any(!nzchar(out$Site) | is.na(out$Site))) {
    stop("[00_master_pipeline] Missing Site labels in inputs/meta_ind.csv for one or more individuals.")
  }
  
  out
}

load_site_meta <- function() {
  path <- file.path(PROJECT_ROOT, "inputs", "site_metadata.csv")
  if (!file.exists(path)) stop("[00_master_pipeline] Missing required file: ", path)
  x <- read_csv_with_comments(path)
  site_col <- resolve_col(x, c("Site", "site", "pop", "population"))
  if (is.na(site_col)) stop("[00_master_pipeline] inputs/site_metadata.csv must contain Site column.")
  x[[site_col]] <- trimws(as.character(x[[site_col]]))
  x <- x[nzchar(x[[site_col]]), , drop = FALSE]
  names(x)[names(x) == site_col] <- "Site"
  x
}

choose_genotype_source <- function() {
  preferred <- c(
    file.path(PROJECT_ROOT, "data", "raw", "poppr.xlsx"),
    file.path(PROJECT_ROOT, "data", "raw", "poppr_avec_dup_E1&2.xlsx")
  )
  
  existing_pref <- preferred[file.exists(preferred)]
  if (length(existing_pref) > 0) return(existing_pref[1])
  
  # fallback discovery
  raw_dir <- file.path(PROJECT_ROOT, "data", "raw")
  if (!dir.exists(raw_dir)) {
    stop("[00_master_pipeline] data/raw not found; cannot rebuild genotype object.")
  }
  files <- list.files(raw_dir, full.names = TRUE, recursive = TRUE)
  files <- files[file.info(files)$isdir %in% FALSE]
  files <- files[grepl("\\.(xlsx|xls|csv|tsv|txt)$", basename(files), ignore.case = TRUE)]
  if (length(files) == 0) {
    stop("[00_master_pipeline] No genotype table found in data/raw.")
  }
  
  # prefer file names containing poppr then shorter names
  score <- function(p) {
    b <- tolower(basename(p))
    s <- 0
    if (grepl("poppr", b)) s <- s + 10
    if (grepl("dup", b)) s <- s - 1
    s - nchar(b) / 10000
  }
  files[order(vapply(files, score, numeric(1)), decreasing = TRUE)][1]
}

read_genotype_table <- function(path) {
  ext <- tolower(sub("^.*\\.", "", path))
  if (ext %in% c("xlsx", "xls")) {
    if (!requireNamespace("readxl", quietly = TRUE)) {
      stop("[00_master_pipeline] readxl package is required to read Excel genotype file: ", path)
    }
    as.data.frame(readxl::read_xlsx(path), stringsAsFactors = FALSE, check.names = FALSE)
  } else if (ext == "csv") {
    read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  } else if (ext %in% c("tsv", "txt")) {
    read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    stop("[00_master_pipeline] Unsupported genotype file extension: ", ext)
  }
}

build_genind_from_table <- function(tbl, canonical_ids_norm, id_to_site) {
  id_col <- resolve_col(tbl, c(
    "Nom_Labo_Échantillons", "Nom_Labo_Echantillons", "Nom_Labo_Echantillon",
    "ind_id", "ind", "individual", "sample", "sampleid", "id"
  ))
  pop_col <- resolve_col(tbl, c("Numéro_Population", "Numero_Population", "Site", "site", "pop", "population"))
  if (is.na(id_col)) stop("[00_master_pipeline] Could not find ID column in genotype table.")
  
  allele_1 <- grep("(_|\\.)1$", names(tbl), perl = TRUE)
  allele_2 <- grep("(_|\\.)2$", names(tbl), perl = TRUE)
  if (length(allele_1) == 0 || length(allele_2) == 0) {
    stop("[00_master_pipeline] Genotype table must contain paired allele columns ending with _1/_2 or .1/.2")
  }
  
  locus1 <- sub("(_|\\.)1$", "", names(tbl)[allele_1], perl = TRUE)
  locus2 <- sub("(_|\\.)2$", "", names(tbl)[allele_2], perl = TRUE)
  shared_loci <- sort(intersect(locus1, locus2))
  if (length(shared_loci) == 0) stop("[00_master_pipeline] No paired loci detected.")
  
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
    stop("[00_master_pipeline] No genotype rows matched curated metadata IDs.")
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
  
  # Build codominant genotype strings per locus
  geno <- data.frame(row.names = ids_raw, stringsAsFactors = FALSE)
  for (loc in shared_loci) {
    c1 <- names(tbl)[grep(paste0("^", gsub("([\\W])", "\\\\\\1", loc, perl = TRUE), "(_|\\\\.)1$"), names(tbl), perl = TRUE)][1]
    c2 <- names(tbl)[grep(paste0("^", gsub("([\\W])", "\\\\\\1", loc, perl = TRUE), "(_|\\\\.)2$"), names(tbl), perl = TRUE)][1]
    
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
  meta_ind <- load_meta_ind()
  meta <- load_site_meta()
  
  canonical_ids <- meta_ind$ind_id
  canonical_ids_norm <- normalize_id(canonical_ids)
  id_to_site <- setNames(meta_ind$Site, canonical_ids_norm)
  
  message("[00_master_pipeline] Canonical sample-universe policy: curated metadata IDs from inputs/meta_ind.csv")
  message("[00_master_pipeline] Intended canonical dataset size = ", length(canonical_ids),
          " (project-consistent target: 278)")
  
  source_file <- choose_genotype_source()
  message("[00_master_pipeline] Chosen genotype source file: ", source_file)
  
  tbl <- read_genotype_table(source_file)
  built <- build_genind_from_table(tbl, canonical_ids_norm, id_to_site)
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
  
  message("[00_master_pipeline] Chosen genotype ID column: ", built$id_col)
  message("[00_master_pipeline] Chosen genotype Site/Pop column: ", ifelse(is.na(built$pop_col), "<not used>", built$pop_col))
  message("[00_master_pipeline] Number of raw genotype IDs: ", built$diagnostics$n_raw_ids,
          " (unique normalized: ", built$diagnostics$n_unique_raw_ids, ")")
  message("[00_master_pipeline] Number of metadata IDs: ", length(canonical_ids))
  message("[00_master_pipeline] Number matched to canonical IDs: ", built$diagnostics$n_matched)
  message("[00_master_pipeline] Number excluded (not in curated metadata universe): ", built$diagnostics$n_excluded)
  if (built$diagnostics$n_excluded > 0) {
    message("[00_master_pipeline] Excluded ID examples (up to 20): ",
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