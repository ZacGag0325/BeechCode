# scripts/00_master_pipeline.R
############################################################
# Build canonical objects required by downstream scripts:
# - gi.rds
# - gi_mll.rds
# - df_ids.rds
# - meta.rds
# Saved to: outputs/v1/objects
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
  
  stop("Cannot locate project root containing scripts/00_master_pipeline.R")
}

PROJECT_ROOT <- find_project_root()
setwd(PROJECT_ROOT)

OUTPUT_OBJ_DIR <- file.path(PROJECT_ROOT, "outputs", "v1", "objects")
dir.create(OUTPUT_OBJ_DIR, recursive = TRUE, showWarnings = FALSE)

message("[00_master_pipeline] Project root: ", PROJECT_ROOT)
message("[00_master_pipeline] Canonical object output dir: ", OUTPUT_OBJ_DIR)

resolve_col <- function(df, choices) {
  nms <- names(df)
  idx <- match(TRUE, tolower(nms) %in% tolower(choices), nomatch = 0)
  if (idx == 0) return(NA_character_)
  nms[idx]
}

sanitize_names <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x)] <- ""
  x <- gsub("^\ufeff", "", x, perl = TRUE)
  x <- gsub("[[:cntrl:]]", "", x)
  x[x == ""] <- "unnamed"
  make.unique(x, sep = "_")
}

normalize_id <- function(x) {
  x <- as.character(x)
  x <- gsub("^\ufeff", "", x, perl = TRUE)
  x <- gsub("[[:cntrl:]]", "", x)
  x <- trimws(x)
  x <- iconv(x, from = "", to = "ASCII//TRANSLIT", sub = "")
  x <- gsub("\u00A0", " ", x, perl = TRUE)
  x <- gsub("\\s+", "", x, perl = TRUE)
  x <- toupper(x)
  x
}

trim_empty_columns <- function(df) {
  if (ncol(df) == 0) return(df)
  keep <- vapply(df, function(col) {
    x <- trimws(as.character(col))
    any(nzchar(x))
  }, logical(1))
  if (!any(keep)) return(df[, 0, drop = FALSE])
  df[, keep, drop = FALSE]
}

read_tabular_robust <- function(path, required_choices = list(), label = "table") {
  if (!file.exists(path)) stop("[00_master_pipeline] Missing required file: ", path)
  
  raw_lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  if (length(raw_lines) == 0) {
    stop("[00_master_pipeline] Input file is empty: ", path)
  }
  
  raw_lines[1] <- sub("^\ufeff", "", raw_lines[1])
  
  lines <- raw_lines
  while (length(lines) > 0) {
    top <- trimws(lines[1])
    if (!nzchar(top) || grepl("^#", top)) {
      lines <- lines[-1]
    } else {
      break
    }
  }
  
  if (length(lines) == 0) {
    stop("[00_master_pipeline] Input file has no tabular content after removing leading comments: ", path)
  }
  
  seps <- c(",", "\t", ";", "|")
  names(seps) <- c("comma", "tab", "semicolon", "pipe")
  
  header_score <- function(df, choices) {
    found <- integer(length(choices))
    if (length(choices) > 0) {
      for (i in seq_along(choices)) {
        found[i] <- !is.na(resolve_col(df, choices[[i]]))
      }
    }
    sum(found)
  }
  
  eval_candidate <- function(sep, sep_label) {
    raw_df <- tryCatch(
      read.table(
        text = lines,
        sep = sep,
        header = FALSE,
        fill = TRUE,
        quote = "\"'",
        comment.char = "",
        stringsAsFactors = FALSE,
        check.names = FALSE,
        na.strings = c("", "NA")
      ),
      error = function(e) NULL
    )
    if (is.null(raw_df) || nrow(raw_df) == 0 || ncol(raw_df) == 0) return(NULL)
    
    header_candidates <- list()
    max_header_row <- min(10L, nrow(raw_df) - 1L)
    if (max_header_row >= 1) {
      for (h in seq_len(max_header_row)) {
        data_h <- raw_df[(h + 1):nrow(raw_df), , drop = FALSE]
        if (nrow(data_h) == 0) next
        
        names(data_h) <- sanitize_names(raw_df[h, , drop = TRUE])
        data_h <- trim_empty_columns(data_h)
        
        req <- header_score(data_h, required_choices)
        score <- req * 1000 + ncol(data_h) + nrow(data_h) / 10000
        header_candidates[[length(header_candidates) + 1]] <- list(
          df = data_h,
          score = score,
          header_row = h,
          header_mode = "header_from_detected_row"
        )
      }
    }
    
    data_b <- raw_df
    names(data_b) <- paste0("V", seq_len(ncol(data_b)))
    data_b <- trim_empty_columns(data_b)
    req_b <- header_score(data_b, required_choices)
    score_b <- req_b * 1000 + ncol(data_b) + nrow(data_b) / 10000
    no_header <- list(
      df = data_b,
      score = score_b,
      header_row = NA_integer_,
      header_mode = "no_header_detected"
    )
    
    all_opts <- c(header_candidates, list(no_header))
    best_i <- which.max(vapply(all_opts, function(x) x$score, numeric(1)))
    best <- all_opts[[best_i]]
    
    list(
      df = best$df,
      sep_label = sep_label,
      header_mode = best$header_mode,
      header_row = best$header_row,
      score = best$score
    )
  }
  
  candidates <- lapply(seq_along(seps), function(i) eval_candidate(seps[[i]], names(seps)[i]))
  candidates <- Filter(Negate(is.null), candidates)
  if (length(candidates) == 0) {
    stop(
      "[00_master_pipeline] Could not parse tabular input file: ", path,
      "\nTried separators: comma, tab, semicolon, pipe."
    )
  }
  
  best_idx <- which.max(vapply(candidates, function(x) x$score, numeric(1)))
  best <- candidates[[best_idx]]
  df <- best$df
  df[] <- lapply(df, function(col) as.character(col))
  
  required_ok <- logical(length(required_choices))
  required_names <- character(length(required_choices))
  if (length(required_choices) > 0) {
    for (i in seq_along(required_choices)) {
      col_hit <- resolve_col(df, required_choices[[i]])
      required_ok[i] <- !is.na(col_hit)
      required_names[i] <- paste(required_choices[[i]], collapse = "/")
    }
  }
  
  message("[00_master_pipeline] Imported ", label, ": ", path)
  message("  parse mode: sep=", best$sep_label,
          ", header_mode=", best$header_mode,
          ifelse(is.na(best$header_row), "", paste0(", header_row=", best$header_row)))
  message("  dimensions: ", nrow(df), " x ", ncol(df))
  message("  first columns: ", paste(head(names(df), 8), collapse = ", "))
  if (length(required_ok) > 0) {
    message("  required columns found: ",
            paste0(required_names, "=", ifelse(required_ok, "YES", "NO"), collapse = " | "))
  }
  
  if (length(required_ok) > 0 && !all(required_ok)) {
    missing_req <- required_names[!required_ok]
    stop(
      "[00_master_pipeline] Malformed input file: ", path,
      "\nMissing required column groups: ", paste(missing_req, collapse = ", "),
      "\nParsed columns were: ", paste(names(df), collapse = ", ")
    )
  }
  
  df
}

print_id_mismatch_diagnostics <- function(gi_ids, df_ids) {
  gi_ids <- as.character(gi_ids)
  df_ids <- as.character(df_ids)
  
  gi_norm <- normalize_id(gi_ids)
  df_norm <- normalize_id(df_ids)
  
  exact_match <- sum(gi_ids %in% df_ids)
  norm_match <- sum(gi_norm %in% df_norm)
  
  gi_only_norm <- setdiff(gi_norm, df_norm)
  df_only_norm <- setdiff(df_norm, gi_norm)
  
  gi_example <- unique(gi_ids[gi_norm %in% gi_only_norm])
  df_example <- unique(df_ids[df_norm %in% df_only_norm])
  
  message("[00_master_pipeline] ID diagnostic summary:")
  message("  gi count: ", length(gi_ids), " | df_ids count: ", length(df_ids))
  message("  exact ID matches: ", exact_match)
  message("  normalized ID matches: ", norm_match)
  message("  gi IDs missing in df_ids (normalized): ", length(gi_only_norm))
  message("  df_ids missing in gi (normalized): ", length(df_only_norm))
  
  if (length(gi_example) > 0) {
    message("  examples gi->missing (up to 25): ", paste(head(gi_example, 25), collapse = ", "))
  }
  if (length(df_example) > 0) {
    message("  examples df_ids->missing (up to 25): ", paste(head(df_example, 25), collapse = ", "))
  }
  
  # quick heuristics
  trim_match <- sum(trimws(gi_ids) %in% trimws(df_ids))
  upper_match <- sum(toupper(trimws(gi_ids)) %in% toupper(trimws(df_ids)))
  de_punct <- function(x) gsub("[^A-Z0-9]", "", toupper(trimws(x)))
  punct_match <- sum(de_punct(gi_ids) %in% de_punct(df_ids))
  message("  heuristic matches after trim only: ", trim_match)
  message("  heuristic matches after trim+case: ", upper_match)
  message("  heuristic matches after trim+case+punctuation-strip: ", punct_match)
}

relabel_mll_ids <- function(raw_id, sample_ids) {
  if (length(raw_id) != length(sample_ids)) {
    stop("MLL relabeling error: raw_id length does not match sample_ids length.")
  }
  ord <- order(sample_ids)
  first_ids <- tapply(seq_along(raw_id), raw_id, function(i) min(ord[match(i, seq_along(raw_id))]))
  stable_levels <- names(sort(unlist(first_ids), decreasing = FALSE))
  map <- setNames(seq_along(stable_levels), stable_levels)
  as.integer(map[as.character(raw_id)])
}

make_mll_from_bruvo <- function(genind_obj, threshold) {
  replen_vec <- rep(2, adegenet::nLoc(genind_obj))
  names(replen_vec) <- adegenet::locNames(genind_obj)
  
  bruvo <- poppr::bruvo.dist(genind_obj, replen = replen_vec)
  bruvo_mat <- as.matrix(bruvo)
  
  if (!is.matrix(bruvo_mat) || nrow(bruvo_mat) != ncol(bruvo_mat)) {
    stop("Invalid Bruvo distance matrix: expected square matrix.")
  }
  if (nrow(bruvo_mat) != adegenet::nInd(genind_obj)) {
    stop("Invalid Bruvo distance matrix: nrow does not match nInd(genind).")
  }
  if (any(is.na(bruvo_mat))) {
    stop("Invalid Bruvo distance matrix: contains NA values.")
  }
  if (any(!is.finite(bruvo_mat))) {
    stop("Invalid Bruvo distance matrix: contains non-finite values.")
  }
  
  rn <- rownames(bruvo_mat)
  cn <- colnames(bruvo_mat)
  ids <- adegenet::indNames(genind_obj)
  
  if (is.null(rn) || is.null(cn) || !identical(rn, cn)) {
    stop("Bruvo matrix row/column names are missing or inconsistent.")
  }
  if (!setequal(rn, ids)) {
    stop("Bruvo matrix sample names do not match indNames(genind).")
  }
  
  bruvo_mat <- bruvo_mat[ids, ids, drop = FALSE]
  bruvo_mat <- (bruvo_mat + t(bruvo_mat)) / 2
  diag(bruvo_mat) <- 0
  
  hc <- stats::hclust(stats::as.dist(bruvo_mat), method = "average")
  raw_id <- stats::cutree(hc, h = threshold)
  stable_id <- relabel_mll_ids(raw_id, sample_ids = ids)
  
  list(
    bruvo_mat = bruvo_mat,
    mll_id = stable_id,
    mll_label = paste0("MLL", stable_id)
  )
}

load_df_ids <- function() {
  path <- file.path(PROJECT_ROOT, "inputs", "meta_ind.csv")
  x <- read_tabular_robust(
    path = path,
    required_choices = list(
      c("ind", "individual", "sample", "sampleid", "id", "ind_id"),
      c("site", "population", "pop")
    ),
    label = "individual metadata"
  )
  
  id_col <- resolve_col(x, c("ind", "individual", "sample", "sampleid", "id", "ind_id"))
  site_col <- resolve_col(x, c("site", "population", "pop"))
  
  x[[id_col]] <- trimws(as.character(x[[id_col]]))
  x[[site_col]] <- trimws(as.character(x[[site_col]]))
  x <- x[nzchar(x[[id_col]]), , drop = FALSE]
  rownames(x) <- NULL
  
  x
}

load_meta <- function() {
  path <- file.path(PROJECT_ROOT, "inputs", "site_metadata.csv")
  x <- read_tabular_robust(
    path = path,
    required_choices = list(
      c("site", "population", "pop"),
      c("latitude", "lat")
    ),
    label = "site metadata"
  )
  
  x
}

find_object_file <- function(fname) {
  search_dirs <- c(
    file.path(PROJECT_ROOT, "outputs", "v1", "objects"),
    file.path(PROJECT_ROOT, "outputs", "objects"),
    file.path(PROJECT_ROOT, "outputs", "v1"),
    file.path(PROJECT_ROOT, "outputs"),
    file.path(PROJECT_ROOT, "data", "objects"),
    file.path(PROJECT_ROOT, "data"),
    file.path(PROJECT_ROOT, "inputs")
  )
  
  for (d in unique(search_dirs)) {
    p <- file.path(d, fname)
    if (file.exists(p)) return(p)
  }
  NA_character_
}

GENOTYPE_SEARCH_DIRS <- c(
  file.path(PROJECT_ROOT, "outputs", "v1", "objects"),
  file.path(PROJECT_ROOT, "outputs", "objects"),
  file.path(PROJECT_ROOT, "outputs", "v1"),
  file.path(PROJECT_ROOT, "outputs"),
  file.path(PROJECT_ROOT, "inputs"),
  file.path(PROJECT_ROOT, "data")
)

discover_genotype_candidates <- function() {
  exts <- "\\.(rds|csv|tsv|txt|xlsx|xls)$"
  all_files <- character(0)
  for (d in unique(GENOTYPE_SEARCH_DIRS)) {
    if (!dir.exists(d)) next
    f <- list.files(d, recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
    if (length(f) == 0) next
    f <- f[file.info(f)$isdir %in% FALSE]
    f <- f[grepl(exts, basename(f), ignore.case = TRUE)]
    all_files <- c(all_files, f)
  }
  all_files <- unique(normalizePath(all_files, winslash = "/", mustWork = FALSE))
  
  if (length(all_files) == 0) {
    return(data.frame(path = character(0), ext = character(0), score = numeric(0), stringsAsFactors = FALSE))
  }
  
  score_file <- function(path) {
    b <- tolower(basename(path))
    s <- 0
    if (grepl("gi\\.rds$", b)) s <- s + 1000
    if (grepl("genind|geno|genotype|microsat|poppr|allele", b)) s <- s + 200
    if (grepl("structure|evanno|likelihood|delta|medk", b)) s <- s - 400
    if (grepl("meta|site_metadata|ids_order", b)) s <- s - 200
    if (grepl("\\.rds$", b)) s <- s + 80
    if (grepl("\\.xlsx$|\\.xls$", b)) s <- s + 40
    s
  }
  
  out <- data.frame(
    path = all_files,
    ext = tolower(sub("^.*\\.", "", all_files)),
    score = vapply(all_files, score_file, numeric(1)),
    stringsAsFactors = FALSE
  )
  
  out <- out[order(-out$score, out$path), , drop = FALSE]
  rownames(out) <- NULL
  out
}

extract_genind_from_object <- function(obj) {
  if (inherits(obj, "genind")) return(obj)
  
  if (is.list(obj)) {
    likely_names <- c("gi", "gen", "genind", "g", "data")
    hit <- intersect(likely_names, names(obj))
    for (nm in hit) {
      if (inherits(obj[[nm]], "genind")) return(obj[[nm]])
    }
    idx <- which(vapply(obj, function(x) inherits(x, "genind"), logical(1)))
    if (length(idx) > 0) return(obj[[idx[1]]])
  }
  
  NULL
}

read_delimited_guess <- function(path) {
  readers <- list(
    function() read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
    function() read.delim(path, stringsAsFactors = FALSE, check.names = FALSE),
    function() read.table(path, sep = ";", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, fill = TRUE, quote = "\"'", comment.char = ""),
    function() read.table(path, sep = "",  header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, fill = TRUE, quote = "\"'", comment.char = "")
  )
  
  for (rd in readers) {
    dat <- tryCatch(rd(), error = function(e) NULL)
    if (is.null(dat) || !is.data.frame(dat) || nrow(dat) == 0 || ncol(dat) == 0) next
    names(dat) <- sanitize_names(names(dat))
    return(dat)
  }
  NULL
}

build_genind_from_table <- function(tbl) {
  nms <- names(tbl)
  
  id_priority <- c(
    "Nom_Labo_Échantillons", "Nom_Labo_Echantillons", "Nom_Labo_Echantillon",
    "ind_id", "sampleid", "sample", "individual", "id", "ind"
  )
  pop_priority <- c("Numéro_Population", "Numero_Population", "site", "population", "pop", "site_id")
  
  pick_col <- function(priority) {
    nms_low <- tolower(nms)
    for (p in priority) {
      i <- match(tolower(p), nms_low, nomatch = 0)
      if (i > 0) return(nms[i])
    }
    NA_character_
  }
  
  id_col <- pick_col(id_priority)
  pop_col <- pick_col(pop_priority)
  
  if (is.na(id_col)) id_col <- resolve_col(tbl, c("ind", "individual", "sample", "sampleid", "id", "ind_id", "nom_labo_échantillons", "nom_labo_echantillons", "nom_labo_echantillon"))
  if (is.na(pop_col)) pop_col <- resolve_col(tbl, c("site", "population", "pop", "numero_population", "numéro_population", "site_id"))
  
  if (is.na(id_col) || is.na(pop_col)) {
    return(NULL)
  }
  
  message("[00_master_pipeline] Genotype table selected columns: id=", id_col, " | pop=", pop_col)
  
  allele_idx <- grep("(_|\\.)[12]$", nms, perl = TRUE)
  if (length(allele_idx) < 2 || length(allele_idx) %% 2 != 0) return(NULL)
  
  locus_base <- sub("(_|\\.)[12]$", "", nms[allele_idx], perl = TRUE)
  locus_tab <- table(locus_base)
  if (any(locus_tab != 2)) return(NULL)
  
  col_1 <- allele_idx[grepl("(_|\\.)1$", nms[allele_idx], perl = TRUE)]
  col_1 <- col_1[order(sub("(_|\\.)1$", "", nms[col_1], perl = TRUE))]
  
  ids <- trimws(as.character(tbl[[id_col]]))
  pops <- trimws(as.character(tbl[[pop_col]]))
  keep <- nzchar(ids)
  ids <- ids[keep]
  pops <- pops[keep]
  sub_tbl <- tbl[keep, , drop = FALSE]
  
  geno_locus <- data.frame(row.names = ids)
  for (idx in col_1) {
    nm1 <- nms[idx]
    nm2 <- sub("1$", "2", nm1)
    j <- match(nm2, nms)
    if (is.na(j)) return(NULL)
    
    locus <- sub("(_|\\.)1$", "", nm1, perl = TRUE)
    a1 <- as.character(sub_tbl[[nm1]])
    a2 <- as.character(sub_tbl[[nms[j]]])
    g <- paste(a1, a2, sep = "/")
    g[grepl("NA", g, fixed = TRUE)] <- NA_character_
    geno_locus[[locus]] <- g
  }
  
  if (ncol(geno_locus) == 0 || nrow(geno_locus) == 0) return(NULL)
  names(geno_locus) <- sanitize_names(names(geno_locus))
  
  gi <- tryCatch(
    adegenet::df2genind(
      geno_locus,
      ploidy = 2,
      ind.names = rownames(geno_locus),
      pop = as.factor(pops),
      type = "codom",
      sep = "/",
      ncode = 3
    ),
    error = function(e) NULL
  )
  
  if (is.null(gi)) return(NULL)
  attr(gi, "source_id_col") <- id_col
  attr(gi, "source_pop_col") <- pop_col
  gi
}

build_gi_from_sources <- function() {
  candidates <- discover_genotype_candidates()
  
  message("[00_master_pipeline] Genotype source search dirs:")
  for (d in unique(GENOTYPE_SEARCH_DIRS)) message("  - ", d)
  
  if (nrow(candidates) == 0) {
    stop("[00_master_pipeline] No genotype candidate files found in configured search directories.")
  }
  
  message("[00_master_pipeline] Genotype candidate files found (top 25 by priority):")
  top_n <- min(25, nrow(candidates))
  for (i in seq_len(top_n)) {
    message("  [", i, "] score=", candidates$score[i], " | ", candidates$path[i])
  }
  
  for (i in seq_len(nrow(candidates))) {
    f <- candidates$path[i]
    ext <- candidates$ext[i]
    
    if (ext == "rds") {
      obj <- tryCatch(readRDS(f), error = function(e) NULL)
      gi <- extract_genind_from_object(obj)
      if (!is.null(gi)) {
        message("[00_master_pipeline] Chosen genotype source: ", f)
        message("[00_master_pipeline] Source mode: loaded genind from RDS")
        attr(gi, "source_path") <- f
        attr(gi, "source_mode") <- "loaded_genind_rds"
        return(gi)
      }
      next
    }
    
    if (ext %in% c("xlsx", "xls")) {
      if (!requireNamespace("readxl", quietly = TRUE)) next
      tbl <- tryCatch(readxl::read_xlsx(f), error = function(e) NULL)
      if (is.null(tbl)) next
      tbl <- as.data.frame(tbl, stringsAsFactors = FALSE, check.names = FALSE)
      names(tbl) <- sanitize_names(names(tbl))
      gi <- build_genind_from_table(tbl)
      if (!is.null(gi)) {
        message("[00_master_pipeline] Chosen genotype source: ", f)
        message("[00_master_pipeline] Source mode: rebuilt genind from spreadsheet table")
        attr(gi, "source_path") <- f
        attr(gi, "source_mode") <- "rebuilt_from_spreadsheet"
        return(gi)
      }
      next
    }
    
    if (ext %in% c("csv", "tsv", "txt")) {
      tbl <- read_delimited_guess(f)
      if (is.null(tbl)) next
      gi <- build_genind_from_table(tbl)
      if (!is.null(gi)) {
        message("[00_master_pipeline] Chosen genotype source: ", f)
        message("[00_master_pipeline] Source mode: rebuilt genind from delimited table")
        attr(gi, "source_path") <- f
        attr(gi, "source_mode") <- "rebuilt_from_delimited"
        return(gi)
      }
      next
    }
  }
  
  stop(
    "[00_master_pipeline] Could not build 'gi' from available files.\n",
    "Checked ", nrow(candidates), " candidate files across search directories.\n",
    "Expected either:\n",
    "  - an RDS containing a genind object (or list element gi/genind),\n",
    "  - or a genotype table with ID + population + paired allele columns (*_1/*_2 or *.1/*.2)."
  )
}

try_load_legacy_objects <- function() {
  targets <- c("gi.rds", "gi_mll.rds", "df_ids.rds", "meta.rds")
  paths <- setNames(vapply(targets, find_object_file, character(1)), targets)
  if (any(is.na(paths))) return(NULL)
  
  message("[00_master_pipeline] Found existing object files; normalizing to outputs/v1/objects")
  objs <- lapply(paths, readRDS)
  names(objs) <- names(paths)
  objs
}

assign_pop_from_df_ids <- function(gi, df_ids) {
  id_col <- resolve_col(df_ids, c("ind", "individual", "sample", "sampleid", "id", "ind_id"))
  site_col <- resolve_col(df_ids, c("site", "population", "pop"))
  
  gi_ids <- as.character(adegenet::indNames(gi))
  df_ids_raw <- as.character(df_ids[[id_col]])
  df_site_raw <- as.character(df_ids[[site_col]])
  
  print_id_mismatch_diagnostics(gi_ids, df_ids_raw)
  
  gi_norm <- normalize_id(gi_ids)
  df_norm <- normalize_id(df_ids_raw)
  
  df_map <- data.frame(df_id = df_ids_raw, df_norm = df_norm, site = df_site_raw, stringsAsFactors = FALSE)
  df_map <- df_map[nzchar(df_map$df_norm), , drop = FALSE]
  
  dup_norm <- df_map$df_norm[duplicated(df_map$df_norm)]
  if (length(dup_norm) > 0) {
    dup_sites <- split(df_map$site[df_map$df_norm %in% dup_norm], df_map$df_norm[df_map$df_norm %in% dup_norm])
    conflicting <- vapply(dup_sites, function(v) length(unique(na.omit(v))) > 1, logical(1))
    if (any(conflicting)) {
      bad <- names(dup_sites)[conflicting]
      stop("[00_master_pipeline] Conflicting Site labels after ID normalization for IDs: ",
           paste(head(bad, 25), collapse = ", "))
    }
    df_map <- df_map[!duplicated(df_map$df_norm), , drop = FALSE]
  }
  
  site_by_norm <- setNames(df_map$site, df_map$df_norm)
  mapped_site <- site_by_norm[gi_norm]
  
  missing_idx <- which(is.na(mapped_site) | !nzchar(mapped_site))
  recovered_n <- 0L
  if (length(missing_idx) > 0) {
    gi_pop <- tryCatch(as.character(adegenet::pop(gi)), error = function(e) rep(NA_character_, adegenet::nInd(gi)))
    can_recover <- missing_idx[!is.na(gi_pop[missing_idx]) & nzchar(gi_pop[missing_idx])]
    if (length(can_recover) > 0) {
      mapped_site[can_recover] <- gi_pop[can_recover]
      recovered_n <- length(can_recover)
    }
  }
  
  final_missing <- which(is.na(mapped_site) | !nzchar(mapped_site))
  message("[00_master_pipeline] Site assignment diagnostics:")
  message("  matched via df_ids: ", length(gi_ids) - length(missing_idx))
  message("  recovered from genotype source pop column: ", recovered_n)
  message("  still unmatched after harmonization: ", length(final_missing))
  
  if (length(final_missing) > 0) {
    miss <- gi_ids[final_missing]
    stop(
      "[00_master_pipeline] Could not map Site for all gi individuals after normalization/recovery.\n",
      "Unmatched IDs (up to 50): ", paste(head(miss, 50), collapse = ", "),
      "\nAttempted normalization: trim, control-char removal, transliteration, whitespace removal, case normalization."
    )
  }
  
  adegenet::pop(gi) <- as.factor(mapped_site)
  gi
}

harmonize_df_ids_to_gi <- function(df_ids, gi) {
  id_col <- resolve_col(df_ids, c("ind", "individual", "sample", "sampleid", "id", "ind_id"))
  site_col <- resolve_col(df_ids, c("site", "population", "pop"))
  
  gi_ids <- as.character(adegenet::indNames(gi))
  gi_site <- as.character(adegenet::pop(gi))
  
  df_ids[[id_col]] <- as.character(df_ids[[id_col]])
  df_ids[[site_col]] <- as.character(df_ids[[site_col]])
  
  df_norm <- normalize_id(df_ids[[id_col]])
  gi_norm <- normalize_id(gi_ids)
  
  missing_in_df <- !(gi_norm %in% df_norm)
  if (any(missing_in_df)) {
    n_add <- sum(missing_in_df)
    add <- as.data.frame(
      matrix(NA_character_, nrow = n_add, ncol = ncol(df_ids)),
      stringsAsFactors = FALSE
    )
    names(add) <- names(df_ids)
    
    add[[id_col]] <- gi_ids[missing_in_df]
    add[[site_col]] <- gi_site[missing_in_df]
    
    # Optional metadata from source pop labels if a numeric site column exists.
    if (site_col %in% names(add)) {
      add[[site_col]] <- as.character(add[[site_col]])
    }
    
    df_ids <- rbind(df_ids, add)
    message("[00_master_pipeline] Added ", n_add, " individuals to df_ids from gi universe.")
  }
  
  rownames(df_ids) <- NULL
  df_ids
}

build_objects <- function() {
  df_ids <- load_df_ids()
  meta <- load_meta()
  
  gi <- build_gi_from_sources()
  message("[00_master_pipeline] Source ID column used (if rebuilt): ", ifelse(is.null(attr(gi, "source_id_col")), "NA", attr(gi, "source_id_col")))
  message("[00_master_pipeline] Source POP column used (if rebuilt): ", ifelse(is.null(attr(gi, "source_pop_col")), "NA", attr(gi, "source_pop_col")))
  
  gi <- assign_pop_from_df_ids(gi, df_ids)
  df_ids <- harmonize_df_ids_to_gi(df_ids, gi)
  
  mll <- make_mll_from_bruvo(gi, threshold = 0)
  mll_labels <- mll$mll_label
  
  id_col <- resolve_col(df_ids, c("ind", "individual", "sample", "sampleid", "id", "ind_id"))
  mll_col <- resolve_col(df_ids, c("mll"))
  if (is.na(mll_col)) {
    df_ids$mll <- NA_character_
    mll_col <- "mll"
  }
  
  mll_map <- setNames(mll_labels, adegenet::indNames(gi))
  df_ids[[mll_col]] <- mll_map[as.character(df_ids[[id_col]])]
  
  keep_mll <- !duplicated(mll_labels)
  gi_mll <- gi[keep_mll, ]
  
  list(
    gi.rds = gi,
    gi_mll.rds = gi_mll,
    df_ids.rds = df_ids,
    meta.rds = meta
  )
}

objs <- try_load_legacy_objects()
if (is.null(objs)) {
  message("[00_master_pipeline] No complete legacy object set found; building objects now.")
  objs <- build_objects()
}

if (!inherits(objs[["gi.rds"]], "genind")) stop("[00_master_pipeline] gi.rds is not a genind object.")
if (!inherits(objs[["gi_mll.rds"]], "genind")) stop("[00_master_pipeline] gi_mll.rds is not a genind object.")
if (!is.data.frame(objs[["df_ids.rds"]])) stop("[00_master_pipeline] df_ids.rds is not a data.frame.")
if (!is.data.frame(objs[["meta.rds"]])) stop("[00_master_pipeline] meta.rds is not a data.frame.")

saveRDS(objs[["gi.rds"]], file.path(OUTPUT_OBJ_DIR, "gi.rds"))
saveRDS(objs[["gi_mll.rds"]], file.path(OUTPUT_OBJ_DIR, "gi_mll.rds"))
saveRDS(objs[["df_ids.rds"]], file.path(OUTPUT_OBJ_DIR, "df_ids.rds"))
saveRDS(objs[["meta.rds"]], file.path(OUTPUT_OBJ_DIR, "meta.rds"))

missing_after <- c("gi.rds", "gi_mll.rds", "df_ids.rds", "meta.rds")
missing_after <- missing_after[!file.exists(file.path(OUTPUT_OBJ_DIR, missing_after))]
if (length(missing_after) > 0) {
  stop("[00_master_pipeline] Failed to write required object files: ", paste(missing_after, collapse = ", "))
}

message("[00_master_pipeline] Saved required objects:")
message("  - ", file.path(OUTPUT_OBJ_DIR, "gi.rds"))
message("  - ", file.path(OUTPUT_OBJ_DIR, "gi_mll.rds"))
message("  - ", file.path(OUTPUT_OBJ_DIR, "df_ids.rds"))
message("  - ", file.path(OUTPUT_OBJ_DIR, "meta.rds"))
message("[00_master_pipeline] nInd(gi)=", adegenet::nInd(objs[["gi.rds"]]),
        " | nInd(gi_mll)=", adegenet::nInd(objs[["gi_mll.rds"]]),
        " | nrow(df_ids)=", nrow(objs[["df_ids.rds"]]),
        " | nrow(meta)=", nrow(objs[["meta.rds"]]))