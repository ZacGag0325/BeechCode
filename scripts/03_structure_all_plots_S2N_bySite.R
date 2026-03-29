# scripts/03_structure_all_plots_S2N_bySite.R
############################################################
# STRUCTURE helper outputs (final individual barplots only)
#
# Biological rationale:
# - STRUCTURE was run externally on the full microsatellite dataset (gi),
#   not on the clone-corrected object, so this script is strictly a reader /
#   plotting / interpretation-support layer for those external results.
# - Ordering sites from south to north is used only to make geographic
#   gradients easier to compare visually across K values.
# - STRUCTURE cluster labels are arbitrary between analyses and across K
#   (for example, Q1 is not inherently the "northern" cluster and Q2 is not
#   inherently the "southern" cluster). Cluster interpretation must always be
#   compared explicitly with geography and the other population-genetic results.
# - For mapping in QGIS, the biologically relevant unit is the site rather than
#   the individual tree, so this script now also collapses individual ancestry
#   coefficients to site-level mean ancestry and computes one geographic
#   centroid per site from per-individual GPS coordinates.
#
# This script IMPROVES readability while keeping the existing analysis intact:
# - Reads externally generated STRUCTURE Q files
# - Builds ONE final individual plot per K (no per-run figures)
# - Builds ONE combined all-K figure
# - Computes mean latitude / longitude centroids per site from per-individual GPS
# - Orders sites strictly SOUTH -> NORTH by mean latitude
# - Reorders individuals so they are grouped by site in that same order
# - Adds separator lines and site labels for publication-ready readability
# - Exports site-level mean Q summaries to support biological interpretation
# - Exports QGIS-ready CSV files for pie-chart ancestry maps
# - Optionally exports a cross-method site diagnostic summary
#
# Outputs:
# - outputs/figures/structure_individual_barplot_K{K}.jpeg
# - outputs/figures/structure_individual_barplot_allK.jpeg
# - outputs/tables/supplementary/structure_run_inventory.csv
# - outputs/tables/supplementary/structure_selected_runs.csv
# - outputs/tables/structure_meanQ_by_site.csv
# - outputs/tables/qgis/structure_meanQ_by_site_K2.csv
# - outputs/tables/qgis/structure_meanQ_by_site_K8.csv
# - outputs/tables/site_genetic_diagnostic_summary.csv (optional, best-effort)
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

source("scripts/_load_objects.R")


cluster_colors <- c(
  Q1 = "#1b9e77",
  Q2 = "#d95f02",
  Q3 = "#7570b3",
  Q4 = "#e7298a",
  Q5 = "#66a61e",
  Q6 = "#e6ab02",
  Q7 = "#a6761d",
  Q8 = "#666666",
  Q9 = "#1f78b4",
  Q10 = "#b2df8a",
  Q11 = "#fb9a99",
  Q12 = "#cab2d6"
)

message("[03_structure] Preparing final STRUCTURE individual barplots and interpretation-support outputs...")
message("[03_structure] Species context: Fagus grandifolia")

STRUCTURE_INTERPRETATION_NOTE <- paste(
  "STRUCTURE cluster labels are arbitrary (for example, Q1 is not inherently north and Q2 is not inherently south).",
  "Interpretation must be compared with geography, AMOVA grouping, DAPC, and differentiation metrics."
)

QGIS_EXPORT_K_VALUES <- c(2L, 8L)
QGIS_DIR <- file.path(TABLES_DIR, "qgis")
dir.create(QGIS_DIR, recursive = TRUE, showWarnings = FALSE)

clamp01 <- function(x, tol = 1e-6) {
  x[x < 0 & x >= -tol] <- 0
  x[x > 1 & x <= 1 + tol] <- 1
  x
}

fail_validation <- function(k, file_base, run_id, msg) {
  stop(
    paste0(
      "[03_structure] Validation failed for K=", k,
      " | file=", file_base,
      " | run=", ifelse(is.na(run_id), "NA", run_id),
      "\n", msg
    ),
    call. = FALSE
  )
}

normalize_site <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("\\uFEFF", "", x, fixed = TRUE)
  x <- gsub("[[:cntrl:]]", "", x)
  x <- gsub("\\s+", " ", x)
  toupper(x)
}

normalize_name <- function(x) {
  x <- iconv(as.character(x), from = "", to = "ASCII//TRANSLIT")
  x[is.na(x)] <- as.character(x[is.na(x)])
  x <- tolower(trimws(x))
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

coerce_numeric_warn <- function(x, column_name, context) {
  raw <- trimws(as.character(x))
  raw[raw == ""] <- NA_character_
  raw <- gsub(",", ".", raw, fixed = TRUE)
  out <- suppressWarnings(as.numeric(raw))
  bad <- !is.na(raw) & is.na(out)
  if (any(bad)) {
    warning(
      context,
      " Non-numeric values detected in ", column_name,
      "; affected rows will be excluded. Example values: ",
      paste(head(unique(raw[bad]), 5), collapse = ", ")
    )
  }
  out
}

coerce_numeric_strict <- function(x, column_name, context) {
  raw <- trimws(as.character(x))
  raw[raw == ""] <- NA_character_
  raw <- gsub(",", ".", raw, fixed = TRUE)
  out <- suppressWarnings(as.numeric(raw))
  bad <- !is.na(raw) & is.na(out)
  if (any(bad)) {
    stop(
      context,
      " Non-numeric values detected in column '", column_name,
      "'. Example values: ",
      paste(head(unique(raw[bad]), 5), collapse = ", "),
      call. = FALSE
    )
  }
  out
}

assert_site_alignment <- function(reference_sites, candidate_sites, context) {
  ref_norm <- unique(normalize_site(reference_sites))
  cand_norm <- unique(normalize_site(candidate_sites))
  missing_in_candidate <- setdiff(ref_norm, cand_norm)
  extra_in_candidate <- setdiff(cand_norm, ref_norm)
  
  if (length(missing_in_candidate) > 0) {
    warning(
      context,
      " Site names present in STRUCTURE-linked metadata but absent from the comparison object: ",
      paste(head(missing_in_candidate, 10), collapse = ", ")
    )
  }
  if (length(extra_in_candidate) > 0) {
    message(
      context,
      " Extra site names present in comparison object but not required for plotting: ",
      paste(head(extra_in_candidate, 10), collapse = ", ")
    )
  }
}

assert_exact_site_strings <- function(reference_tbl, candidate_tbl, context) {
  ref <- reference_tbl %>% distinct(Site_norm, Site)
  cand <- candidate_tbl %>% distinct(Site_norm, Site)
  merged <- inner_join(ref, cand, by = "Site_norm", suffix = c("_reference", "_candidate"))
  
  mismatch <- merged %>%
    filter(Site_reference != Site_candidate)
  
  if (nrow(mismatch) > 0) {
    stop(
      context,
      " Site labels must match exactly between STRUCTURE-linked metadata and coordinate metadata. Mismatch examples: ",
      paste(
        paste0(mismatch$Site_reference[seq_len(min(10, nrow(mismatch)))], " != ", mismatch$Site_candidate[seq_len(min(10, nrow(mismatch)))]),
        collapse = "; "
      ),
      call. = FALSE
    )
  }
}

needs_readxl <- function(path) {
  tolower(tools::file_ext(path)) %in% c("xlsx", "xls")
}

clean_column_names <- function(x) {
  x <- gsub("^\\ufeff", "", x)
  x <- trimws(x)
  x
}

read_csv_with_comments <- function(path) {
  raw_lines <- readLines(path, warn = FALSE)
  if (length(raw_lines) == 0) {
    return(data.frame())
  }
  keep <- !(grepl("^\\s*#", raw_lines) | grepl("^\\s*$", raw_lines))
  cleaned <- raw_lines[keep]
  if (length(cleaned) == 0) {
    return(data.frame())
  }
  utils::read.csv(
    text = paste(cleaned, collapse = "\n"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

read_table_file <- function(path, sheet = NA_character_) {
  ext <- tolower(tools::file_ext(path))
  
  df <- if (ext == "csv") {
    read_csv_with_comments(path)
  } else if (ext %in% c("tsv", "txt")) {
    read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  } else if (ext %in% c("xlsx", "xls")) {
    if (!requireNamespace("readxl", quietly = TRUE)) {
      stop("[03_structure] readxl is required to import Excel metadata files: ", path, call. = FALSE)
    }
    readxl::read_excel(path, sheet = if (is.na(sheet)) 1 else sheet, .name_repair = "minimal") %>%
      as.data.frame(stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    stop("Unsupported metadata file extension: ", ext, call. = FALSE)
  }
  
  if (ncol(df) > 0) {
    names(df) <- clean_column_names(names(df))
  }
  df
}

get_sheets <- function(path) {
  if (!needs_readxl(path)) return(NA_character_)
  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop("[03_structure] readxl is required to inspect workbook sheets: ", path, call. = FALSE)
  }
  readxl::excel_sheets(path)
}

list_candidate_metadata_files <- function() {
  roots <- c(
    file.path(PROJECT_ROOT, "inputs"),
    file.path(PROJECT_ROOT, "data", "raw"),
    file.path(PROJECT_ROOT, "data")
  )
  roots <- unique(roots[dir.exists(roots)])
  if (length(roots) == 0) {
    return(character(0))
  }
  
  files <- unlist(lapply(roots, function(root) {
    list.files(
      root,
      pattern = "\\.(csv|tsv|txt|xlsx|xls)$",
      recursive = TRUE,
      full.names = TRUE,
      ignore.case = TRUE
    )
  }), use.names = FALSE)
  
  files <- unique(files[file.info(files)$isdir %in% FALSE])
  files[!grepl("structure|q_files|qgis|deltak|medk|lnpk", basename(files), ignore.case = TRUE)]
}

summarize_metadata_table <- function(df) {
  id_col <- resolve_col_ci(df, c(
    DF_IDS_ID_CHOICES,
    "Nom_Labo_Échantillons", "Nom_Labo_Echantillons", "Nom_Labo_Echantillon"
  ))
  site_col <- resolve_col_ci(df, c(
    DF_IDS_SITE_CHOICES,
    "Numéro_Population", "Numero_Population"
  ))
  lat_col <- resolve_col_ci(df, c("latitude", "lat"))
  lon_col <- resolve_col_ci(df, c("longitude", "lon", "long"))
  
  list(
    nrow = nrow(df),
    ncol = ncol(df),
    id_col = id_col,
    site_col = site_col,
    lat_col = lat_col,
    lon_col = lon_col
  )
}

scan_coordinate_sources <- function() {
  files <- list_candidate_metadata_files()
  out <- list()
  
  for (path in files) {
    sheets <- tryCatch(get_sheets(path), error = function(e) NULL)
    if (is.null(sheets)) next
    
    if (length(sheets) == 1 && is.na(sheets[1])) {
      df <- tryCatch(read_table_file(path), error = function(e) NULL)
      if (is.null(df) || nrow(df) == 0 || ncol(df) == 0) next
      out[[length(out) + 1]] <- list(
        path = path,
        sheet = NA_character_,
        df = df,
        summary = summarize_metadata_table(df)
      )
    } else {
      for (sheet_name in sheets) {
        df <- tryCatch(read_table_file(path, sheet = sheet_name), error = function(e) NULL)
        if (is.null(df) || nrow(df) == 0 || ncol(df) == 0) next
        out[[length(out) + 1]] <- list(
          path = path,
          sheet = sheet_name,
          df = df,
          summary = summarize_metadata_table(df)
        )
      }
    }
  }
  
  out
}

select_individual_coordinate_source <- function(scanned, n_reference_ids) {
  candidates <- Filter(function(x) {
    !is.na(x$summary$id_col) &&
      !is.na(x$summary$site_col) &&
      !is.na(x$summary$lat_col) &&
      !is.na(x$summary$lon_col)
  }, scanned)
  
  if (length(candidates) == 0) {
    stop(
      "[03_structure] Could not find an individual metadata table containing ID, Site, Latitude, and Longitude columns. ",
      "Add or expose one in inputs/ or data/raw/ so site centroids can be computed from per-individual GPS.",
      call. = FALSE
    )
  }
  
  score_candidate <- function(x) {
    df <- x$df
    sm <- x$summary
    id_values <- trimws(as.character(df[[sm$id_col]]))
    id_values <- id_values[nzchar(id_values) & !is.na(id_values)]
    site_values <- trimws(as.character(df[[sm$site_col]]))
    site_values <- site_values[nzchar(site_values) & !is.na(site_values)]
    
    uniq_id_prop <- if (length(id_values) == 0) 0 else length(unique(id_values)) / length(id_values)
    uniq_site_n <- length(unique(site_values))
    b <- tolower(basename(x$path))
    
    score <- 0
    if (grepl("field|meta|metadata|sample|ind|coord|location|gps", b)) score <- score + 30
    if (grepl("site", b)) score <- score + 5
    if (needs_readxl(x$path)) score <- score + 5
    if (x$summary$nrow >= n_reference_ids) score <- score + 40
    score <- score + min(x$summary$nrow / 10, 50)
    score <- score + uniq_id_prop * 100
    score <- score + min(uniq_site_n, 100)
    score
  }
  
  candidates[[which.max(vapply(candidates, score_candidate, numeric(1)))]]
}

resolve_individual_coordinate_table <- function(reference_ids, canonical_site_map) {
  scanned <- scan_coordinate_sources()
  selected <- select_individual_coordinate_source(scanned, n_reference_ids = length(reference_ids))
  sm <- selected$summary
  df <- selected$df
  
  message("[03_structure] Coordinate metadata source: ", selected$path)
  message("[03_structure] Coordinate metadata sheet: ", ifelse(is.na(selected$sheet), "<file>", selected$sheet))
  message("[03_structure] Coordinate metadata columns: ID=", sm$id_col, ", Site=", sm$site_col, ", Latitude=", sm$lat_col, ", Longitude=", sm$lon_col)
  
  coord_tbl <- data.frame(
    ind_id = trimws(as.character(df[[sm$id_col]])),
    Site = trimws(as.character(df[[sm$site_col]])),
    Latitude = coerce_numeric_strict(df[[sm$lat_col]], sm$lat_col, "[03_structure]"),
    Longitude = coerce_numeric_strict(df[[sm$lon_col]], sm$lon_col, "[03_structure]"),
    stringsAsFactors = FALSE
  ) %>%
    filter(nzchar(ind_id))
  
  if (nrow(coord_tbl) == 0) {
    stop("[03_structure] The selected coordinate metadata table is empty after removing blank individual IDs.", call. = FALSE)
  }
  
  if (anyDuplicated(normalize_id(coord_tbl$ind_id))) {
    dup_ids <- unique(coord_tbl$ind_id[duplicated(normalize_id(coord_tbl$ind_id))])
    stop(
      "[03_structure] Coordinate metadata contains duplicated individual IDs. Examples: ",
      paste(head(dup_ids, 10), collapse = ", "),
      call. = FALSE
    )
  }
  
  if (any(!is.na(coord_tbl$Latitude) & (coord_tbl$Latitude < -90 | coord_tbl$Latitude > 90))) {
    bad_ids <- coord_tbl$ind_id[!is.na(coord_tbl$Latitude) & (coord_tbl$Latitude < -90 | coord_tbl$Latitude > 90)]
    stop(
      "[03_structure] Latitude values outside [-90, 90] detected. Example IDs: ",
      paste(head(unique(bad_ids), 10), collapse = ", "),
      call. = FALSE
    )
  }
  
  if (any(!is.na(coord_tbl$Longitude) & (coord_tbl$Longitude < -180 | coord_tbl$Longitude > 180))) {
    bad_ids <- coord_tbl$ind_id[!is.na(coord_tbl$Longitude) & (coord_tbl$Longitude < -180 | coord_tbl$Longitude > 180)]
    stop(
      "[03_structure] Longitude values outside [-180, 180] detected. Example IDs: ",
      paste(head(unique(bad_ids), 10), collapse = ", "),
      call. = FALSE
    )
  }
  
  reference_tbl <- data.frame(
    ind_id = reference_ids,
    Site_structure = unname(canonical_site_map[normalize_id(reference_ids)]),
    stringsAsFactors = FALSE
  )
  
  coord_tbl$ind_norm <- normalize_id(coord_tbl$ind_id)
  reference_tbl$ind_norm <- normalize_id(reference_tbl$ind_id)
  joined <- reference_tbl %>%
    left_join(coord_tbl, by = "ind_norm", suffix = c("_structure", "_metadata"))
  
  missing_coord_ids <- joined$ind_id_structure[
    is.na(joined$Latitude) | is.na(joined$Longitude) | is.na(joined$Site)
  ]
  if (length(missing_coord_ids) > 0) {
    stop(
      "[03_structure] Missing coordinate metadata for one or more STRUCTURE individuals. Example IDs: ",
      paste(head(unique(missing_coord_ids), 10), collapse = ", "),
      call. = FALSE
    )
  }
  
  site_norm_structure <- normalize_site(joined$Site_structure)
  site_norm_metadata <- normalize_site(joined$Site)
  mismatch_norm <- which(site_norm_structure != site_norm_metadata)
  if (length(mismatch_norm) > 0) {
    bad_rows <- joined[mismatch_norm[seq_len(min(10, length(mismatch_norm)))], , drop = FALSE]
    stop(
      "[03_structure] Site names do not match between STRUCTURE-linked metadata and coordinate metadata for some individuals. Examples: ",
      paste(
        paste0(bad_rows$ind_id_structure, " [", bad_rows$Site_structure, " vs ", bad_rows$Site, "]"),
        collapse = "; "
      ),
      call. = FALSE
    )
  }
  
  exact_mismatch <- which(joined$Site_structure != joined$Site)
  if (length(exact_mismatch) > 0) {
    bad_rows <- joined[exact_mismatch[seq_len(min(10, length(exact_mismatch)))], , drop = FALSE]
    stop(
      "[03_structure] Site names must match exactly between STRUCTURE-linked metadata and coordinate metadata. Examples: ",
      paste(
        paste0(bad_rows$ind_id_structure, " [", bad_rows$Site_structure, " vs ", bad_rows$Site, "]"),
        collapse = "; "
      ),
      call. = FALSE
    )
  }
  
  joined %>%
    transmute(
      ind_id = ind_id_structure,
      Site = Site_structure,
      Site_norm = normalize_site(Site_structure),
      Latitude = Latitude,
      Longitude = Longitude
    )
}

build_site_coordinate_table <- function(individual_coord_tbl) {
  if (nrow(individual_coord_tbl) == 0) {
    stop("[03_structure] Individual coordinate table is empty; site centroids cannot be computed.", call. = FALSE)
  }
  
  site_tbl <- individual_coord_tbl %>%
    group_by(Site_norm) %>%
    summarise(
      Site = dplyr::first(Site),
      mean_latitude = mean(Latitude, na.rm = TRUE),
      mean_longitude = mean(Longitude, na.rm = TRUE),
      n_individuals = n(),
      .groups = "drop"
    ) %>%
    arrange(mean_latitude, Site) %>%
    mutate(site_order_south_to_north = seq_len(n()))
  
  if (any(!is.finite(site_tbl$mean_latitude) | !is.finite(site_tbl$mean_longitude))) {
    bad_sites <- site_tbl$Site[!is.finite(site_tbl$mean_latitude) | !is.finite(site_tbl$mean_longitude)]
    stop(
      "[03_structure] Missing site centroid coordinates after averaging individuals. Site(s): ",
      paste(head(unique(bad_sites), 20), collapse = ", "),
      call. = FALSE
    )
  }
  
  if (any(site_tbl$mean_latitude < -90 | site_tbl$mean_latitude > 90)) {
    bad_sites <- site_tbl$Site[site_tbl$mean_latitude < -90 | site_tbl$mean_latitude > 90]
    stop(
      "[03_structure] Invalid site centroid latitude values detected for site(s): ",
      paste(head(unique(bad_sites), 20), collapse = ", "),
      call. = FALSE
    )
  }
  
  if (any(site_tbl$mean_longitude < -180 | site_tbl$mean_longitude > 180)) {
    bad_sites <- site_tbl$Site[site_tbl$mean_longitude < -180 | site_tbl$mean_longitude > 180]
    stop(
      "[03_structure] Invalid site centroid longitude values detected for site(s): ",
      paste(head(unique(bad_sites), 20), collapse = ", "),
      call. = FALSE
    )
  }
  
  site_tbl
}

build_base_order <- function(ids_vec, site_map_final, site_coord_tbl) {
  out <- data.frame(
    Individual = ids_vec,
    Site = unname(site_map_final[normalize_id(ids_vec)]),
    stringsAsFactors = FALSE
  )
  
  if (any(is.na(out$Site) | !nzchar(out$Site))) {
    missing_ids <- out$Individual[is.na(out$Site) | !nzchar(out$Site)]
    stop(
      "[03_structure] Could not assign Site labels to all STRUCTURE individuals. Example IDs: ",
      paste(head(unique(missing_ids), 10), collapse = ", "),
      call. = FALSE
    )
  }
  
  out$Site_norm <- normalize_site(out$Site)
  out$SiteLat <- site_coord_tbl$mean_latitude[match(out$Site_norm, site_coord_tbl$Site_norm)]
  out$SiteLon <- site_coord_tbl$mean_longitude[match(out$Site_norm, site_coord_tbl$Site_norm)]
  out$SiteOrder <- site_coord_tbl$site_order_south_to_north[match(out$Site_norm, site_coord_tbl$Site_norm)]
  
  if (any(is.na(out$SiteLat) | is.na(out$SiteLon) | is.na(out$SiteOrder))) {
    bad_sites <- unique(out$Site[is.na(out$SiteLat) | is.na(out$SiteLon) | is.na(out$SiteOrder)])
    stop(
      "[03_structure] Missing site centroid coordinates for site(s): ",
      paste(head(bad_sites, 20), collapse = ", "),
      call. = FALSE
    )
  }
  out <- out %>%
    arrange(SiteOrder, SiteLat, Site, Individual)
  
  out$PlotIndex <- seq_len(nrow(out))
  out
}

numeric_columns <- function(df) {
  as.data.frame(lapply(df, function(col) {
    x <- trimws(as.character(col))
    x[x == ""] <- NA_character_
    suppressWarnings(as.numeric(x))
  }), check.names = FALSE, stringsAsFactors = FALSE)
}

window_stats <- function(mat, tol = 1e-6) {
  na_count <- sum(is.na(mat))
  lo <- sum(mat < -tol, na.rm = TRUE)
  hi <- sum(mat > 1 + tol, na.rm = TRUE)
  valid_rows <- complete.cases(mat)
  if (any(valid_rows)) {
    rs <- rowSums(mat[valid_rows, , drop = FALSE])
    row_min <- min(rs)
    row_max <- max(rs)
    row_bad <- sum(abs(rs - 1) > 1e-3)
    row_dev <- mean(abs(rs - 1))
  } else {
    row_min <- NA_real_
    row_max <- NA_real_
    row_bad <- Inf
    row_dev <- Inf
  }
  list(na = na_count, lo = lo, hi = hi, row_min = row_min, row_max = row_max, row_bad = row_bad, row_dev = row_dev)
}

choose_q_columns <- function(num_df, k_hint = NA_integer_, tol = 1e-6) {
  m <- ncol(num_df)
  if (m < 2) return(NULL)
  
  get_candidates <- function(k_values) {
    out <- list()
    for (k in k_values) {
      if (k < 2 || k > m) next
      for (s in seq_len(m - k + 1)) {
        idx <- s:(s + k - 1)
        mat <- as.matrix(num_df[, idx, drop = FALSE])
        storage.mode(mat) <- "numeric"
        st <- window_stats(mat, tol = tol)
        score <- st$na * 1e7 + (st$lo + st$hi) * 1e6 + st$row_bad * 1e3 + st$row_dev * 100
        out[[length(out) + 1]] <- list(idx = idx, k = k, q = mat, stats = st, score = score)
      }
    }
    out
  }
  
  cands <- list()
  if (!is.na(k_hint)) cands <- get_candidates(k_hint)
  if (length(cands) == 0) cands <- get_candidates(2:m)
  if (length(cands) == 0) return(NULL)
  
  best_i <- which.min(vapply(cands, function(z) z$score, numeric(1)))
  cands[[best_i]]
}

inspect_q_file_lines <- function(path, expected_k = NA_integer_) {
  txt <- tryCatch(readLines(path, warn = FALSE), error = function(e) character(0))
  if (length(txt) == 0) {
    return(list(total = 0L, blank = 0L, numeric_like = 0L, wrong_token_count = 0L, duplicate_lines = 0L))
  }
  
  raw <- trimws(txt)
  nonblank <- raw[nzchar(raw)]
  blank_n <- sum(!nzchar(raw))
  
  token_n <- vapply(strsplit(nonblank, "\\s+"), length, integer(1))
  numeric_like <- vapply(nonblank, function(line) {
    toks <- strsplit(line, "\\s+")[[1]]
    all(!is.na(suppressWarnings(as.numeric(toks))))
  }, logical(1))
  
  if (is.na(expected_k)) {
    numeric_rows <- sum(numeric_like)
    wrong_tok <- NA_integer_
  } else {
    numeric_rows <- sum(numeric_like & token_n == expected_k)
    wrong_tok <- sum(numeric_like & token_n != expected_k)
  }
  
  list(
    total = length(txt),
    blank = blank_n,
    numeric_like = numeric_rows,
    wrong_token_count = wrong_tok,
    duplicate_lines = sum(duplicated(nonblank))
  )
}

read_q_matrix <- function(path, k_hint = NA_integer_, tol = 1e-6) {
  readers <- list(
    list(label = "whitespace_no_header", fun = function() read.table(path, header = FALSE, sep = "", stringsAsFactors = FALSE, check.names = FALSE, fill = TRUE, comment.char = "", quote = "")),
    list(label = "whitespace_header",    fun = function() read.table(path, header = TRUE,  sep = "", stringsAsFactors = FALSE, check.names = FALSE, fill = TRUE, comment.char = "", quote = "")),
    list(label = "tab_no_header",        fun = function() read.table(path, header = FALSE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, fill = TRUE, comment.char = "", quote = "")),
    list(label = "tab_header",           fun = function() read.table(path, header = TRUE,  sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, fill = TRUE, comment.char = "", quote = "")),
    list(label = "csv_header",           fun = function() read.csv(path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)),
    list(label = "csv_no_header",        fun = function() read.csv(path, header = FALSE, stringsAsFactors = FALSE, check.names = FALSE))
  )
  
  for (r in readers) {
    dat <- tryCatch(r$fun(), error = function(e) NULL)
    if (is.null(dat) || nrow(dat) == 0 || ncol(dat) == 0) next
    
    num_df <- numeric_columns(dat)
    if (nrow(num_df) == 0 || ncol(num_df) < 2) next
    
    keep_non_empty <- vapply(num_df, function(x) sum(!is.na(x)) > 0, logical(1))
    num_df <- num_df[, keep_non_empty, drop = FALSE]
    if (ncol(num_df) < 2) next
    
    first_col <- num_df[[1]]
    if (!anyNA(first_col) && length(first_col) > 1 && all(abs(first_col - seq_len(length(first_col))) < 1e-8)) {
      num_df <- num_df[, -1, drop = FALSE]
      if (ncol(num_df) < 2) next
    }
    
    sel <- choose_q_columns(num_df, k_hint = k_hint, tol = tol)
    if (is.null(sel)) next
    
    q <- sel$q
    q <- q[rowSums(!is.na(q)) > 0, , drop = FALSE]
    colnames(q) <- paste0("Q", seq_len(ncol(q)))
    
    line_diag <- inspect_q_file_lines(path, expected_k = ncol(q))
    
    return(list(ok = TRUE, matrix = q, reader = r$label, q_indices = sel$idx, stats = sel$stats, line_diag = line_diag))
  }
  
  list(ok = FALSE, matrix = NULL, reader = NA_character_, q_indices = integer(0), stats = NULL, line_diag = inspect_q_file_lines(path, expected_k = k_hint))
}

summarize_candidate <- function(file_base, k_infer, parsed, reason) {
  if (!parsed$ok) {
    reason_print <- if (grepl("^REJECTED:", reason)) reason else paste0("REJECTED: ", reason)
    message("[03_structure] Candidate: ", file_base,
            " | inferred K=", ifelse(is.na(k_infer), "NA", k_infer),
            " | ", reason_print)
    return(invisible(NULL))
  }
  
  q <- parsed$matrix
  st <- parsed$stats
  ld <- parsed$line_diag
  message("[03_structure] Candidate: ", file_base,
          " | inferred K=", ifelse(is.na(k_infer), "NA", k_infer),
          " | rows=", nrow(q),
          " | cols=", ncol(q),
          " | q_col_indices=", paste(parsed$q_indices, collapse = ","),
          " | numeric=TRUE",
          " | min=", signif(min(q, na.rm = TRUE), 6),
          " | max=", signif(max(q, na.rm = TRUE), 6),
          " | rowSumRange=", signif(st$row_min, 6), "-", signif(st$row_max, 6),
          " | lines(total/blank/numeric_like)=", ld$total, "/", ld$blank, "/", ld$numeric_like,
          " | duplicate_nonblank_rows=", ld$duplicate_lines,
          " | ", reason)
}

validate_selected_run <- function(k_num, run_info, q_df, q_cols, id_reference,
                                  tol = 1e-6, row_tol = 5e-3,
                                  hard_row_tol = 5e-2,
                                  renormalize_rows = TRUE) {
  q_mat <- as.matrix(q_df[, q_cols, drop = FALSE])
  storage.mode(q_mat) <- "numeric"
  
  na_count <- sum(is.na(q_mat))
  lt_count <- sum(q_mat < -tol, na.rm = TRUE)
  gt_count <- sum(q_mat > 1 + tol, na.rm = TRUE)
  
  q_mat <- clamp01(q_mat, tol = tol)
  row_sums <- rowSums(q_mat)
  abs_dev <- abs(row_sums - 1)
  max_abs_dev <- max(abs_dev)
  bad_row_sums <- sum(abs_dev > row_tol)
  hard_bad_row_sums <- sum(abs_dev > hard_row_tol)
  
  duplicated_ids <- sum(duplicated(q_df$Individual))
  missing_ids <- sum(!(id_reference %in% q_df$Individual))
  
  message(sprintf("[03_structure] K=%d validation details:", k_num))
  message("  selected file/run: ", run_info$file_base, " / ", ifelse(is.na(run_info$run), "NA", run_info$run))
  message("  number of individuals: ", nrow(q_df))
  message("  number of Q columns: ", length(q_cols))
  message("  Q column names: ", paste(q_cols, collapse = ", "))
  message("  all Q columns numeric: TRUE")
  message("  total NA in Q columns: ", na_count)
  message("  Q min/max: ", signif(min(q_mat, na.rm = TRUE), 6), " / ", signif(max(q_mat, na.rm = TRUE), 6))
  message("  row-sum range: ", signif(min(row_sums), 6), " - ", signif(max(row_sums), 6))
  message("  max |rowSum-1|: ", signif(max_abs_dev, 6))
  message("  rows with row sums not close to 1: ", bad_row_sums)
  message("  duplicated individual IDs after merging: ", duplicated_ids)
  message("  missing IDs after matching to reference order: ", missing_ids)
  
  if (na_count > 0) fail_validation(k_num, run_info$file_base, run_info$run, paste0("NA values in Q columns: ", na_count))
  if (lt_count > 0 || gt_count > 0) {
    fail_validation(k_num, run_info$file_base, run_info$run,
                    paste0("Out-of-range Q values detected. <0 count=", lt_count, ", >1 count=", gt_count))
  }
  if (duplicated_ids > 0 || missing_ids > 0) {
    fail_validation(k_num, run_info$file_base, run_info$run,
                    paste0("ID matching issue: duplicated_ids=", duplicated_ids, ", missing_ids=", missing_ids))
  }
  if (hard_bad_row_sums > 0) {
    fail_validation(k_num, run_info$file_base, run_info$run,
                    paste0("Row sums show true invalidity for ", hard_bad_row_sums,
                           " rows (hard tol=", hard_row_tol,
                           "; max |rowSum-1|=", signif(max_abs_dev, 6), ")."))
  }
  
  if (renormalize_rows) {
    renorm_idx <- which(row_sums > 0 & abs_dev > 1e-12)
    renorm_rows <- length(renorm_idx)
    if (renorm_rows > 0) {
      message("  renormalizing ", renorm_rows,
              " rows to sum exactly 1 (prevents stacked-bar clipping at y=1).")
      q_mat[renorm_idx, ] <- q_mat[renorm_idx, , drop = FALSE] / row_sums[renorm_idx]
      q_mat <- clamp01(q_mat, tol = tol)
      row_sums_after <- rowSums(q_mat)
      message("  post-renormalization row-sum range: ",
              signif(min(row_sums_after), 6), " - ", signif(max(row_sums_after), 6))
    }
  }
  
  non_finite_after <- sum(!is.finite(q_mat))
  if (non_finite_after > 0) {
    fail_validation(k_num, run_info$file_base, run_info$run,
                    paste0("Non-finite values introduced after renormalization/clamping: ", non_finite_after))
  }
  
  q_df[, q_cols] <- as.data.frame(q_mat, stringsAsFactors = FALSE)
  q_df
}

validate_plotting_df <- function(plot_df, y_col = "Q", context = "", tol = 1e-6,
                                 check_k_col = NULL, check_cluster = TRUE) {
  y <- suppressWarnings(as.numeric(plot_df[[y_col]]))
  y_na <- sum(is.na(y))
  y_non_finite <- sum(!is.finite(y))
  y_lt0 <- sum(y < 0, na.rm = TRUE)
  y_gt1 <- sum(y > 1, na.rm = TRUE)
  y_gt1_tiny <- sum(y > 1 & y <= 1 + tol, na.rm = TRUE)
  y_lt0_tiny <- sum(y < 0 & y >= -tol, na.rm = TRUE)
  y_min <- if (all(is.na(y))) NA_real_ else min(y, na.rm = TRUE)
  y_max <- if (all(is.na(y))) NA_real_ else max(y, na.rm = TRUE)
  
  missing_cluster <- if (check_cluster && "Cluster" %in% names(plot_df)) {
    sum(is.na(plot_df$Cluster) | !nzchar(as.character(plot_df$Cluster)))
  } else 0L
  missing_x <- if ("PlotIndex" %in% names(plot_df)) sum(is.na(plot_df$PlotIndex)) else NA_integer_
  missing_k <- if (!is.null(check_k_col) && check_k_col %in% names(plot_df)) sum(is.na(plot_df[[check_k_col]])) else 0L
  
  message("[03_structure] Plot data diagnostics (", context, "):")
  message("  rows: ", nrow(plot_df))
  message("  NA in ", y_col, ": ", y_na)
  message("  non-finite in ", y_col, ": ", y_non_finite)
  message("  ", y_col, " < 0: ", y_lt0, " (tiny within tol: ", y_lt0_tiny, ")")
  message("  ", y_col, " > 1: ", y_gt1, " (tiny within tol: ", y_gt1_tiny, ")")
  message("  ", y_col, " min/max: ", signif(y_min, 8), " / ", signif(y_max, 8))
  message("  missing Cluster assignment rows: ", missing_cluster)
  message("  missing PlotIndex rows: ", missing_x)
  if (!is.null(check_k_col)) message("  missing ", check_k_col, " rows: ", missing_k)
  
  if (y_non_finite > 0 || y_na > 0) {
    stop("[03_structure] Invalid plotting data in ", context,
         ": non-finite or NA ancestry values detected (NA=", y_na,
         ", non-finite=", y_non_finite, ").", call. = FALSE)
  }
  
  hard_lt0 <- sum(y < -tol, na.rm = TRUE)
  hard_gt1 <- sum(y > 1 + tol, na.rm = TRUE)
  if (hard_lt0 > 0 || hard_gt1 > 0) {
    stop("[03_structure] Invalid plotting data in ", context,
         ": ancestry values outside [0,1] beyond tolerance (lt0=", hard_lt0,         ", gt1=", hard_gt1, ").", call. = FALSE)
  }
  
  y <- clamp01(y, tol = tol)
  if (any(!is.finite(y))) {
    stop("[03_structure] Invalid plotting data in ", context,
         ": non-finite values remain after clamping.", call. = FALSE)
  }
  
  if (missing_cluster > 0 || (!is.na(missing_x) && missing_x > 0) || missing_k > 0) {
    stop("[03_structure] Invalid plotting data in ", context,
         ": missing plotting keys (Cluster/PlotIndex/KLabel).", call. = FALSE)
  }
  
  plot_df[[y_col]] <- y
  plot_df
}

build_structure_mean_q_table <- function(per_k_q_list, site_summary_tbl) {
  if (length(per_k_q_list) == 0) {
    return(data.frame())
  }
  
  out <- bind_rows(lapply(per_k_q_list, function(x) {
    q_df <- x$q_df
    q_cols <- x$q_cols
    q_df %>%
      select(Site, Site_norm, all_of(q_cols)) %>%
      group_by(Site, Site_norm) %>%
      summarise(across(all_of(q_cols), ~mean(.x, na.rm = TRUE), .names = "mean_{.col}"), .groups = "drop") %>%
      mutate(K = x$K)
  })) %>%
    left_join(site_summary_tbl, by = c("Site", "Site_norm")) %>%
    relocate(K, site_order_south_to_north, Site, mean_latitude, mean_longitude, n_individuals)
  
  out$interpretation_note <- STRUCTURE_INTERPRETATION_NOTE
  out %>% arrange(K, site_order_south_to_north, Site)
}

validate_qgis_site_table <- function(site_q_tbl, k_num, tol = 1e-4) {
  required_cols <- c("Site", "Latitude", "Longitude", "n_individuals", paste0("Q", seq_len(k_num)))
  missing_cols <- setdiff(required_cols, names(site_q_tbl))
  if (length(missing_cols) > 0) {
    stop(
      "[03_structure] QGIS export table for K=", k_num,
      " is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
  
  if (any(!nzchar(site_q_tbl$Site) | is.na(site_q_tbl$Site))) {
    stop("[03_structure] QGIS export table for K=", k_num, " has blank Site values.", call. = FALSE)
  }
  if (any(is.na(site_q_tbl$Latitude) | is.na(site_q_tbl$Longitude))) {
    bad_sites <- site_q_tbl$Site[is.na(site_q_tbl$Latitude) | is.na(site_q_tbl$Longitude)]
    stop(
      "[03_structure] Missing site centroid coordinates for K=", k_num,
      " export. Site(s): ", paste(head(unique(bad_sites), 20), collapse = ", "),
      call. = FALSE
    )
  }
  if (any(site_q_tbl$Latitude < -90 | site_q_tbl$Latitude > 90)) {
    bad_sites <- site_q_tbl$Site[site_q_tbl$Latitude < -90 | site_q_tbl$Latitude > 90]
    stop(
      "[03_structure] Invalid centroid latitude for K=", k_num,
      " export. Site(s): ", paste(head(unique(bad_sites), 20), collapse = ", "),
      call. = FALSE
    )
  }
  if (any(site_q_tbl$Longitude < -180 | site_q_tbl$Longitude > 180)) {
    bad_sites <- site_q_tbl$Site[site_q_tbl$Longitude < -180 | site_q_tbl$Longitude > 180]
    stop(
      "[03_structure] Invalid centroid longitude for K=", k_num,
      " export. Site(s): ", paste(head(unique(bad_sites), 20), collapse = ", "),
      call. = FALSE
    )
  }
  if (any(is.na(site_q_tbl$n_individuals) | site_q_tbl$n_individuals < 1)) {
    bad_sites <- site_q_tbl$Site[is.na(site_q_tbl$n_individuals) | site_q_tbl$n_individuals < 1]
    stop(
      "[03_structure] Invalid n_individuals for K=", k_num,
      " export. Site(s): ", paste(head(unique(bad_sites), 20), collapse = ", "),
      call. = FALSE
    )
  }
  
  q_cols <- paste0("Q", seq_len(k_num))
  q_mat <- as.matrix(site_q_tbl[, q_cols, drop = FALSE])
  storage.mode(q_mat) <- "numeric"
  
  if (any(!is.finite(q_mat))) {
    bad_sites <- site_q_tbl$Site[rowSums(!is.finite(q_mat)) > 0]
    stop(
      "[03_structure] Missing or non-finite mean ancestry values in QGIS export for K=", k_num,
      ". Site(s): ", paste(head(unique(bad_sites), 20), collapse = ", "),
      call. = FALSE
    )
  }
  
  row_sums <- rowSums(q_mat)
  bad_sum <- abs(row_sums - 1) > tol
  if (any(bad_sum)) {
    stop(
      "[03_structure] Site mean ancestry proportions do not sum approximately to 1 for K=", k_num,
      ". Site(s): ", paste(head(site_q_tbl$Site[bad_sum], 20), collapse = ", "),
      "; row sums: ", paste(signif(row_sums[bad_sum], 6), collapse = ", "),
      call. = FALSE
    )
  }
  
  message("[03_structure] K=", k_num, " site-level mean Q row-sum check passed. Range = ", signif(min(row_sums), 6), " - ", signif(max(row_sums), 6))
  
  invisible(site_q_tbl)
}

export_qgis_site_mean_q <- function(per_k_q_list, site_summary_tbl, export_k_values, out_dir) {
  available_k <- sort(vapply(per_k_q_list, function(x) x$K, integer(1)))
  missing_k <- setdiff(export_k_values, available_k)
  if (length(missing_k) > 0) {
    stop(
      "[03_structure] Requested QGIS export K value(s) are not available among selected STRUCTURE runs: ",
      paste(missing_k, collapse = ", "),
      ". Available selected K values: ",
      paste(available_k, collapse = ", "),
      call. = FALSE
    )
  }
  
  exported_files <- character(0)
  
  for (k_num in export_k_values) {
    per_k <- per_k_q_list[[match(k_num, vapply(per_k_q_list, function(x) x$K, integer(1)))]]
    q_cols <- per_k$q_cols
    
    # QGIS pie-chart outputs use one row per site, with a site centroid derived
    # from the mean latitude/longitude of all individuals assigned to that site.
    site_q_tbl <- per_k$q_df %>%
      select(Site, Site_norm, all_of(q_cols)) %>%
      group_by(Site, Site_norm) %>%
      summarise(
        n_individuals_from_q = n(),
        across(all_of(q_cols), ~mean(.x, na.rm = TRUE)),
        .groups = "drop"
      ) %>%
      left_join(site_summary_tbl, by = c("Site", "Site_norm")) %>%
      arrange(site_order_south_to_north, Site) %>%
      transmute(
        Site = Site,
        Latitude = mean_latitude,
        Longitude = mean_longitude,
        n_individuals = n_individuals,
        n_individuals_from_q = n_individuals_from_q,
        !!!rlang::syms(q_cols)
      )
    
    n_mismatch <- site_q_tbl$n_individuals != site_q_tbl$n_individuals_from_q
    if (any(is.na(n_mismatch) | n_mismatch)) {
      bad_sites <- site_q_tbl$Site[is.na(n_mismatch) | n_mismatch]
      stop(
        "[03_structure] Site individual counts disagreed between coordinate metadata and STRUCTURE rows for K=", k_num,
        ". Site(s): ", paste(head(unique(bad_sites), 20), collapse = ", "),
        call. = FALSE
      )
    }
    
    site_q_tbl <- site_q_tbl %>%
      select(-n_individuals_from_q)
    
    validate_qgis_site_table(site_q_tbl, k_num = k_num, tol = 1e-4)
    
    out_file <- file.path(out_dir, sprintf("structure_meanQ_by_site_K%d.csv", k_num))
    write.csv(site_q_tbl, out_file, row.names = FALSE)
    exported_files <- c(exported_files, out_file)
    message("[03_structure] Saved QGIS pie-chart table for K=", k_num, ": ", out_file)
  }
  
  message("[03_structure] Exported QGIS K values: ", paste(export_k_values, collapse = ", "))
  message("[03_structure] Ordered site list used for QGIS exports: ", paste(site_summary_tbl$Site, collapse = " -> "))
  
  invisible(exported_files)
}

build_optional_site_diagnostic <- function(site_summary_tbl, structure_mean_q) {
  amova_path <- file.path(TABLES_DIR, "amova_site_region_groups.csv")
  dapc_centroids_path <- file.path(TABLES_DIR, "dapc_group_centroids.csv")
  jost_path <- file.path(MATRICES_DIR, "pairwise_jostD.csv")
  
  missing_inputs <- c(
    if (!file.exists(amova_path)) "amova_site_region_groups.csv" else character(0),
    if (!file.exists(dapc_centroids_path)) "dapc_group_centroids.csv" else character(0),
    if (!file.exists(jost_path)) "pairwise_jostD.csv" else character(0)
  )
  
  if (length(missing_inputs) > 0) {
    warning(
      "[03_structure] Optional site genetic diagnostic summary will be skipped because required input(s) are missing: ",
      paste(missing_inputs, collapse = ", ")
    )
    return(NULL)
  }
  
  amova_tbl <- read.csv(amova_path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c("Site", "Region") %in% names(amova_tbl))) {
    warning("[03_structure] Optional diagnostic summary skipped: amova_site_region_groups.csv lacks required Site/Region columns.")
    return(NULL)
  }
  
  dapc_tbl <- read.csv(dapc_centroids_path, stringsAsFactors = FALSE, check.names = FALSE)
  required_dapc_cols <- c("Site", "LD1_centroid", "LD2_centroid")
  if (!all(required_dapc_cols %in% names(dapc_tbl))) {
    warning("[03_structure] Optional diagnostic summary skipped: dapc_group_centroids.csv lacks required centroid columns.")
    return(NULL)
  }
  
  jost_df <- read.csv(jost_path, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  jost_mat <- as.matrix(jost_df)
  suppressWarnings(storage.mode(jost_mat) <- "numeric")
  if (!is.matrix(jost_mat) || nrow(jost_mat) != ncol(jost_mat)) {
    warning("[03_structure] Optional diagnostic summary skipped: pairwise_jostD.csv is not a square numeric matrix.")
    return(NULL)
  }
  
  mean_jost <- data.frame(
    Site = rownames(jost_mat),
    mean_pairwise_JostD_to_other_sites = apply(jost_mat, 1, function(x) {
      x <- as.numeric(x)
      mean(x[is.finite(x) & x != 0], na.rm = TRUE)
    }),
    stringsAsFactors = FALSE
  )
  
  assert_site_alignment(site_summary_tbl$Site, amova_tbl$Site, "[03_structure diagnostic]")
  assert_site_alignment(site_summary_tbl$Site, dapc_tbl$Site, "[03_structure diagnostic]")
  assert_site_alignment(site_summary_tbl$Site, rownames(jost_mat), "[03_structure diagnostic]")
  
  structure_wide <- structure_mean_q %>%
    select(-interpretation_note) %>%
    pivot_longer(
      cols = starts_with("mean_Q"),
      names_to = "cluster",
      values_to = "mean_q"
    ) %>%
    mutate(cluster = paste0("K", K, "_", cluster)) %>%
    select(-K) %>%
    pivot_wider(names_from = cluster, values_from = mean_q)
  
  out <- site_summary_tbl %>%
    transmute(
      site_order_south_to_north = site_order_south_to_north,
      Site = Site,
      mean_latitude = mean_latitude,
      mean_longitude = mean_longitude,
      n_individuals = n_individuals
    ) %>%
    left_join(amova_tbl %>% select(Site, Region), by = "Site") %>%
    left_join(structure_wide, by = c("Site", "site_order_south_to_north")) %>%
    left_join(dapc_tbl %>% select(Site, LD1_centroid, LD2_centroid), by = "Site") %>%
    left_join(mean_jost, by = "Site") %>%
    arrange(site_order_south_to_north, Site)
  
  out$interpretation_note <- STRUCTURE_INTERPRETATION_NOTE
  out
}

# ----------------------------
# 1) Individual IDs and Site map
# ----------------------------
df_ids_cols <- resolve_df_ids_columns(df_ids, context = "[03_structure]", require = TRUE)
id_col_dfids <- df_ids_cols$id_col
site_col_dfids <- df_ids_cols$site_col

ids_dfids <- trimws(as.character(df_ids[[id_col_dfids]]))
ids_dfids <- ids_dfids[nzchar(ids_dfids)]
site_map_dfids <- setNames(as.character(df_ids[[site_col_dfids]]), normalize_id(ids_dfids))

load_ids_order_from_raw <- function() {
  ids_path <- file.path(PROJECT_ROOT, "data", "structure", "ids_order_from_raw.csv")
  if (!file.exists(ids_path)) {
    return(list(ids = character(0), source = "ids_order_from_raw.csv (missing)"))
  }
  
  ids_df <- tryCatch(
    read.csv(ids_path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL
  )
  if (is.null(ids_df) || nrow(ids_df) == 0) {
    return(list(ids = character(0), source = "ids_order_from_raw.csv (unreadable)"))
  }
  
  id_col <- resolve_col_ci(ids_df, c("ind", "individual", "sample", "sampleid", "id", "ind_id"))
  if (is.na(id_col)) id_col <- names(ids_df)[1]
  
  ids <- trimws(as.character(ids_df[[id_col]]))
  ids <- ids[nzchar(ids)]
  ids <- unique(ids)
  
  list(ids = ids, source = "data/structure/ids_order_from_raw.csv")
}

extract_structure_order_from_results <- function() {
  res_dir <- file.path(PROJECT_ROOT, "data", "structure", "STRUCTURE_ZG_HEG-results")
  if (!dir.exists(res_dir)) {
    return(list(ids = character(0), pop = integer(0), source = "STRUCTURE results (missing)"))
  }
  
  files <- list.files(res_dir, full.names = TRUE)
  files <- files[file.info(files)$isdir %in% FALSE]
  if (length(files) == 0) {
    return(list(ids = character(0), pop = integer(0), source = "STRUCTURE results (empty)"))
  }
  
  b <- basename(files)
  k_vals <- suppressWarnings(as.integer(gsub(".*K([0-9]+).*", "\\1", b, perl = TRUE)))
  ord <- order(ifelse(is.na(k_vals) | k_vals < 2, Inf, k_vals), b)
  
  for (i in ord) {
    f <- files[i]
    txt <- tryCatch(readLines(f, warn = FALSE), error = function(e) character(0))
    if (length(txt) == 0) next
    
    anc_start <- grep("^Inferred ancestry of individuals:", txt)
    if (length(anc_start) == 0) next
    block <- txt[(anc_start[1] + 1):length(txt)]
    
    keep <- grepl("^\\s*[0-9]+\\s+\\S+\\s+\\([^)]*\\)\\s+[0-9]+\\s*:\\s+", block)
    lines <- block[keep]
    if (length(lines) == 0) next
    
    labels <- sub("^\\s*[0-9]+\\s+(\\S+)\\s+\\([^)]*\\)\\s+[0-9]+\\s*:.*$", "\\1", lines, perl = TRUE)
    pops <- suppressWarnings(as.integer(sub("^\\s*[0-9]+\\s+\\S+\\s+\\([^)]*\\)\\s+([0-9]+)\\s*:.*$", "\\1", lines, perl = TRUE)))
    labels <- trimws(labels)
    labels <- labels[nzchar(labels)]
    
    if (length(labels) > 0) {
      return(list(ids = labels, pop = pops[seq_along(labels)], source = paste0("", f)))
    }
  }
  
  list(ids = character(0), pop = integer(0), source = "STRUCTURE results (no ancestry block parsed)")
}

# ----------------------------
# 2) Q file discovery
# ----------------------------
find_q_files <- function() {
  q_root <- file.path(PROJECT_ROOT, "data", "structure", "Q_Files")
  if (!dir.exists(q_root)) {
    q_root_alt <- file.path(PROJECT_ROOT, "data", "structure", "Q_files")
    if (dir.exists(q_root_alt)) q_root <- q_root_alt
  }
  
  message("[03_structure] Searching folder: ", q_root)
  
  if (!dir.exists(q_root)) {
    return(list(root = q_root, files = character(0), ignored = character(0)))
  }
  
  files <- list.files(q_root, recursive = TRUE, full.names = TRUE, pattern = "(?i)\\.(q|txt|csv)$")
  files <- files[file.info(files)$isdir %in% FALSE]
  b <- basename(files)
  
  ignore_idx <- grepl("(?i)summary|extracted|evanno|likelihood|readme", b)
  ignored <- unique(files[ignore_idx])
  kept <- unique(files[!ignore_idx])
  
  list(root = q_root, files = kept, ignored = ignored)
}

extract_k_run <- function(path) {
  b <- basename(path)
  k_hit <- regmatches(b, regexpr("(?i)K[0-9]+", b, perl = TRUE))
  K <- if (length(k_hit) == 0 || !nzchar(k_hit)) NA_integer_ else as.integer(gsub("[^0-9]", "", k_hit))
  
  run <- suppressWarnings(as.integer(sub("(?i).*K[0-9]+[-_]?([0-9]+).*$", "\\1", b, perl = TRUE)))
  if (is.na(run)) run <- suppressWarnings(as.integer(sub("(?i).*[_-]run[_-]?([0-9]+).*$", "\\1", b, perl = TRUE)))
  
  list(K = K, run = run)
}

# ----------------------------
# 3) Parse all runs and choose ONE final run per K
# ----------------------------
ids_order_raw <- load_ids_order_from_raw()
ids_results <- extract_structure_order_from_results()

scan <- find_q_files()
if (length(scan$ignored) > 0) message("[03_structure] Ignored helper files: ", length(scan$ignored))

if (length(scan$files) == 0) {
  message("[03_structure] No external STRUCTURE Q files found; checked: ", scan$root)
  message("[03_structure] Done.")
} else {
  parsed_candidates <- list()
  run_inventory <- list()
  
  for (f in scan$files) {
    meta_name <- extract_k_run(f)
    parsed <- read_q_matrix(f, k_hint = meta_name$K, tol = 1e-6)
    
    if (!parsed$ok) {
      summarize_candidate(basename(f), meta_name$K, parsed, "could not parse a numeric Q matrix")
      run_inventory[[length(run_inventory) + 1]] <- data.frame(
        file = basename(f), full_path = f, K = meta_name$K, run = meta_name$run,
        reader = parsed$reader, n_rows = NA_integer_, n_cols = NA_integer_,
        q_col_indices = NA_character_, mean_abs_row_sum_deviation = NA_real_,
        status = "unreadable", reason = "could_not_parse_numeric_q_matrix",
        stringsAsFactors = FALSE
      )
      next
    }
    
    q_mat <- parsed$matrix
    K_detected <- ncol(q_mat)
    K_final <- if (!is.na(meta_name$K)) meta_name$K else K_detected
    if (!is.na(meta_name$K) && meta_name$K != K_detected) K_final <- K_detected
    
    row_dev <- mean(abs(rowSums(clamp01(q_mat, 1e-6), na.rm = TRUE) - 1), na.rm = TRUE)
    
    reject_reason <- NULL
    status <- "parsed_candidate"
    if (K_final == 1) {
      reject_reason <- "K1 skipped"
      status <- "skipped_k1"
    }
    
    ld <- parsed$line_diag
    cause_msg <- "ACCEPTED as parsed candidate"
    if (ld$blank > 0 && ld$numeric_like == nrow(q_mat)) {
      cause_msg <- paste0(cause_msg, " (blank lines ignored; no extra non-data rows kept)")
    }
    if (!is.na(ld$wrong_token_count) && ld$wrong_token_count > 0) {
      cause_msg <- paste0(cause_msg, " (non-Q token-count rows excluded)")
    }
    
    summarize_candidate(
      basename(f),
      meta_name$K,
      parsed,
      ifelse(is.null(reject_reason), cause_msg, paste0("REJECTED: ", reject_reason))
    )
    
    run_inventory[[length(run_inventory) + 1]] <- data.frame(
      file = basename(f), full_path = f, K = K_final, run = meta_name$run,
      reader = parsed$reader, n_rows = nrow(q_mat), n_cols = ncol(q_mat),
      q_col_indices = paste(parsed$q_indices, collapse = ";"),
      mean_abs_row_sum_deviation = row_dev,
      status = status,
      reason = ifelse(is.null(reject_reason), "candidate", reject_reason),
      stringsAsFactors = FALSE
    )
    
    if (is.null(reject_reason)) {
      parsed_candidates[[length(parsed_candidates) + 1]] <- list(
        file = f,
        file_base = basename(f),
        K = K_final,
        run = meta_name$run,
        reader = parsed$reader,
        q = q_mat,
        mad = row_dev
      )
    }
  }
  
  if (length(run_inventory) > 0) {
    inv_df <- bind_rows(run_inventory) %>% arrange(K, run, file)
    inv_file <- file.path(TABLES_SUPP_DIR, "structure_run_inventory.csv")
    write.csv(inv_df, inv_file, row.names = FALSE)
    message("[03_structure] Saved run inventory: ", inv_file)
  }
  
  if (length(parsed_candidates) == 0) {
    message("[03_structure] No usable K>=2 Q runs after parsing.")
    message("[03_structure] Done.")
  } else {
    parsed_runs <- Filter(function(x) x$K >= 2, parsed_candidates)
    row_counts <- sort(unique(vapply(parsed_runs, function(x) nrow(x$q), integer(1))))
    if (length(row_counts) != 1) {
      stop("[03_structure] Parsed Q files have inconsistent row counts across runs: ", paste(row_counts, collapse = ", "), call. = FALSE)
    }
    q_n <- row_counts[1]
    
    ref_candidates <- list(
      list(name = ids_order_raw$source, ids = ids_order_raw$ids),
      list(name = ids_results$source, ids = unique(ids_results$ids)),
      list(name = "df_ids object order", ids = unique(ids_dfids))
    )
    
    ref_lengths <- vapply(ref_candidates, function(x) length(x$ids), integer(1))
    message("[03_structure] Reference lengths (for matching parsed n=", q_n, "): ",
            paste0(vapply(ref_candidates, function(x) x$name, character(1)), "=", ref_lengths, collapse = " | "))
    
    idx_match <- which(ref_lengths == q_n)
    if (length(idx_match) == 0) {
      stop("[03_structure] No reference ID source matches parsed Q row count n=", q_n,
           ". Check ids_order_from_raw.csv and STRUCTURE input export order.", call. = FALSE)
    }
    
    ref_pick <- idx_match[1]
    id_reference <- ref_candidates[[ref_pick]]$ids
    ref_source_used <- ref_candidates[[ref_pick]]$name
    message("[03_structure] Using reference ID source: ", ref_source_used, " (n=", length(id_reference), ")")
    
    missing_in_df <- setdiff(normalize_id(id_reference), normalize_id(ids_dfids))
    extra_in_df <- setdiff(normalize_id(ids_dfids), normalize_id(id_reference))
    if (length(missing_in_df) > 0 || length(extra_in_df) > 0) {
      message("[03_structure] ID mismatch diagnostic:")
      message("  IDs in reference but not df_ids: ", length(missing_in_df),
              ifelse(length(missing_in_df) > 0, paste0(" (e.g., ", paste(head(missing_in_df, 5), collapse = ", "), ")"), ""))
      message("  IDs in df_ids but not reference: ", length(extra_in_df),
              ifelse(length(extra_in_df) > 0, paste0(" (e.g., ", paste(head(extra_in_df, 5), collapse = ", "), ")"), ""))
    }
    
    site_map_final <- site_map_dfids
    missing_site_before <- sum(is.na(site_map_final[normalize_id(id_reference)]))
    
    if (length(ids_results$ids) > 0 && length(ids_results$pop) == length(ids_results$ids)) {
      lab_pop <- data.frame(Individual = ids_results$ids, PopIdx = ids_results$pop, stringsAsFactors = FALSE)
      known <- lab_pop %>%
        mutate(Site = site_map_dfids[normalize_id(Individual)]) %>%
        filter(!is.na(Site), !is.na(PopIdx))
      
      if (nrow(known) > 0) {
        pop_to_site <- known %>%
          count(PopIdx, Site, name = "n") %>%
          group_by(PopIdx) %>%
          arrange(desc(n), Site, .by_group = TRUE) %>%
          slice(1) %>%
          ungroup() %>%
          select(PopIdx, Site)
        
        inferred <- lab_pop %>%
          left_join(pop_to_site, by = "PopIdx")
        inferred_map <- setNames(inferred$Site, normalize_id(inferred$Individual))
        
        need_fill <- id_reference[is.na(site_map_final[normalize_id(id_reference)])]
        if (length(need_fill) > 0) {
          fill_vals <- inferred_map[normalize_id(need_fill)]
          site_map_final[normalize_id(need_fill)] <- fill_vals
        }
      }
    }
    
    missing_site_after <- sum(is.na(site_map_final[normalize_id(id_reference)]))
    recovered_site_n <- missing_site_before - missing_site_after
    if (missing_site_before > 0) {
      message("[03_structure] Site metadata recovery diagnostic:")
      message("  missing Site labels before recovery: ", missing_site_before)
      message("  recovered via STRUCTURE results metadata: ", recovered_site_n)
      message("  still unmatched after recovery: ", missing_site_after)
    }
    
    if (missing_site_after > 0) {
      unresolved_ids <- id_reference[is.na(site_map_final[normalize_id(id_reference)])]
      stop(
        "[03_structure] Some STRUCTURE individuals still lack Site labels after metadata recovery. Example IDs: ",
        paste(head(unique(unresolved_ids), 10), collapse = ", "),
        call. = FALSE
      )
    }
    
    individual_coord_tbl <- resolve_individual_coordinate_table(
      reference_ids = id_reference,
      canonical_site_map = site_map_final
    )
    site_coord_tbl <- build_site_coordinate_table(individual_coord_tbl)
    
    assert_exact_site_strings(
      reference_tbl = site_coord_tbl %>% select(Site_norm, Site),
      candidate_tbl = individual_coord_tbl %>% select(Site_norm, Site),
      context = "[03_structure]"
    )
    
    base_order_df <- build_base_order(id_reference, site_map_final, site_coord_tbl)
    if (nrow(base_order_df) != q_n) {
      stop("[03_structure] Internal alignment error: reference order n=", nrow(base_order_df), " but parsed Q n=", q_n, call. = FALSE)
    }
    
    site_blocks <- base_order_df %>%
      group_by(Site, Site_norm) %>%
      summarise(
        n = n(),
        SiteLat = dplyr::first(SiteLat),
        SiteLon = dplyr::first(SiteLon),
        SiteOrder = dplyr::first(SiteOrder),
        .groups = "drop"
      ) %>%
      arrange(SiteOrder, SiteLat, Site) %>%
      mutate(
        xmin = cumsum(dplyr::lag(n, default = 0)) + 0.5,
        xmax = cumsum(n) + 0.5,
        xmid = (xmin + xmax) / 2
      )
    
    final_site_order <- site_blocks$Site
    message("[03_structure] Final ordered site list (south -> north):")
    message("[03_structure] ", paste(final_site_order, collapse = " -> "))
    
    site_summary_tbl <- site_blocks %>%
      transmute(
        site_order_south_to_north = SiteOrder,
        Site = Site,
        Site_norm = Site_norm,
        mean_latitude = SiteLat,
        mean_longitude = SiteLon,
        n_individuals = n
      )
    
    assert_site_alignment(site_summary_tbl$Site, unique(as.character(df_ids_mll$Site)), "[03_structure]")
    
    separators <- site_blocks$xmax[-nrow(site_blocks)]
    k_levels <- sort(unique(as.integer(names(split(parsed_runs, sapply(parsed_runs, function(x) x$K))))))
    
    runs_by_k <- split(parsed_runs, sapply(parsed_runs, function(x) x$K))
    selected_rows <- list()
    allk_plot_data <- list()
    validations_passed <- 0L
    per_k_q_for_summary <- list()
    
    for (k_name in names(runs_by_k)) {
      k_runs <- runs_by_k[[k_name]]
      k_num <- as.integer(k_name)
      
      score_df <- bind_rows(lapply(k_runs, function(x) {
        data.frame(
          file = x$file,
          file_base = x$file_base,
          K = x$K,
          run = ifelse(is.na(x$run), Inf, x$run),
          run_report = x$run,
          reader = x$reader,
          mad = x$mad,
          stringsAsFactors = FALSE
        )
      })) %>% arrange(mad, run, file_base)
      
      best <- k_runs[[match(score_df$file[1], sapply(k_runs, function(x) x$file))]]
      
      message("[03_structure] Found ", length(k_runs), " usable run(s) for K=", k_num)
      message("[03_structure] Using best run for K=", k_num,
              " (run=", ifelse(is.na(best$run), "NA", best$run),
              ", reader=", best$reader,
              ", mean|rowSum-1|=", signif(best$mad, 4), ")")
      message("[03_structure] Reminder for K=", k_num, ": ", STRUCTURE_INTERPRETATION_NOTE)
      
      q_df <- cbind(base_order_df, as.data.frame(best$q, stringsAsFactors = FALSE, check.names = FALSE))
      q_cols <- grep("^Q", names(q_df), value = TRUE)
      
      q_df <- validate_selected_run(
        k_num = k_num,
        run_info = best,
        q_df = q_df,
        q_cols = q_cols,
        id_reference = base_order_df$Individual,
        tol = 1e-6,
        row_tol = 5e-3,
        hard_row_tol = 5e-2,
        renormalize_rows = TRUE
      )
      
      per_k_q_for_summary[[length(per_k_q_for_summary) + 1]] <- list(K = k_num, q_df = q_df, q_cols = q_cols)
      
      plot_df <- q_df %>%
        select(PlotIndex, Site, all_of(q_cols)) %>%
        pivot_longer(cols = all_of(q_cols), names_to = "Cluster", values_to = "Q") %>%
        mutate(Cluster = factor(Cluster, levels = q_cols))
      
      plot_df <- validate_plotting_df(
        plot_df,
        y_col = "Q",
        context = paste0("single-K K=", k_num, " (", best$file_base, ")"),
        tol = 1e-6,
        check_cluster = TRUE
      )
      
      plot_df$K <- k_num
      plot_df$KLabel <- factor(paste0("K=", k_num), levels = paste0("K=", k_levels))
      allk_plot_data[[length(allk_plot_data) + 1]] <- plot_df
      
      p <- ggplot(plot_df, aes(x = PlotIndex, y = Q, fill = Cluster)) +
        geom_col(width = 1) +
        geom_vline(xintercept = separators, linewidth = 0.25, color = "grey25") +
        scale_fill_manual(values = cluster_colors[q_cols], breaks = q_cols, drop = FALSE) +
        scale_x_continuous(
          breaks = site_blocks$xmid,
          labels = site_blocks$Site,
          expand = c(0, 0)
        ) +
        scale_y_continuous(expand = c(0, 0)) +
        coord_cartesian(ylim = c(0, 1), expand = FALSE, clip = "on") +
        theme_bw(base_size = 11) +
        theme(
          panel.grid = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
          axis.ticks.x = element_blank(),
          legend.position = "right",
          plot.title = element_text(face = "bold")
        ) +
        labs(
          title = sprintf("STRUCTURE ancestry barplot (K=%d)", k_num),
          subtitle = paste(
            "Individuals grouped by site; sites ordered left-to-right from south to north by mean site latitude.",
            STRUCTURE_INTERPRETATION_NOTE
          ),
          x = "Site blocks (south -> north)",
          y = "Ancestry proportion"
        )
      
      out_fig <- file.path(FIGURES_DIR, sprintf("structure_individual_barplot_K%d.jpeg", k_num))
      ggsave(out_fig, p, width = 13.5, height = 5.8, dpi = 320)
      message("[03_structure] Saved: ", out_fig)
      
      selected_rows[[length(selected_rows) + 1]] <- data.frame(
        K = k_num,
        selected_file = best$file_base,
        selected_path = best$file,
        selected_run = best$run,
        reader = best$reader,
        mean_abs_row_sum_deviation = best$mad,
        n_runs_available = length(k_runs),
        method = "best_single_run_no_cluster_relabeling",
        interpretation_note = STRUCTURE_INTERPRETATION_NOTE,
        stringsAsFactors = FALSE
      )
      
      validations_passed <- validations_passed + 1L
    }
    
    if (length(allk_plot_data) > 0) {
      combined_plot_df <- bind_rows(allk_plot_data)
      combined_clusters <- paste0("Q", seq_len(max(combined_plot_df$K, na.rm = TRUE)))
      combined_plot_df$Cluster <- factor(combined_plot_df$Cluster, levels = combined_clusters)
      
      combined_plot_df <- validate_plotting_df(
        combined_plot_df,
        y_col = "Q",
        context = "combined all-K plot",
        tol = 1e-6,
        check_k_col = "KLabel",
        check_cluster = TRUE
      )
      
      p_all <- ggplot(combined_plot_df, aes(x = PlotIndex, y = Q, fill = Cluster)) +
        geom_col(width = 1) +
        geom_vline(xintercept = separators, linewidth = 0.2, color = "grey25") +
        scale_fill_manual(values = cluster_colors[combined_clusters], breaks = combined_clusters, drop = FALSE) +
        scale_x_continuous(
          breaks = site_blocks$xmid,
          labels = site_blocks$Site,
          expand = c(0, 0)
        ) +
        scale_y_continuous(expand = c(0, 0)) +
        coord_cartesian(ylim = c(0, 1), expand = FALSE, clip = "on") +
        facet_grid(rows = vars(KLabel)) +
        theme_bw(base_size = 11) +
        theme(
          panel.grid = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
          axis.ticks.x = element_blank(),
          legend.position = "right",
          strip.text.y = element_text(face = "bold")
        ) +
        labs(
          title = "STRUCTURE ancestry barplots across K",
          subtitle = paste(
            "The same individual and site order is used for every K so spatial comparisons are direct.",
            STRUCTURE_INTERPRETATION_NOTE
          ),
          x = "Site blocks (south -> north)",
          y = "Ancestry proportion"
        )
      
      out_all <- file.path(FIGURES_DIR, "structure_individual_barplot_allK.jpeg")
      ggsave(out_all, p_all, width = 13.5, height = max(5.8, 2.3 * length(unique(combined_plot_df$K))), dpi = 320)
      message("[03_structure] Saved: ", out_all)
    }
    
    if (length(selected_rows) > 0) {
      selected_df <- bind_rows(selected_rows) %>% arrange(K)
      selected_file <- file.path(TABLES_SUPP_DIR, "structure_selected_runs.csv")
      write.csv(selected_df, selected_file, row.names = FALSE)
      message("[03_structure] Saved selected-run summary: ", selected_file)
    }
    
    structure_mean_q <- build_structure_mean_q_table(per_k_q_for_summary, site_summary_tbl)
    if (nrow(structure_mean_q) > 0) {
      mean_q_file <- file.path(TABLES_DIR, "structure_meanQ_by_site.csv")
      write.csv(structure_mean_q, mean_q_file, row.names = FALSE)
      message("[03_structure] Saved site-level mean Q summary: ", mean_q_file)
    } else {
      warning("[03_structure] structure_meanQ_by_site.csv was not written because no per-K Q summaries were available.")
    }
    
    export_qgis_site_mean_q(
      per_k_q_list = per_k_q_for_summary,
      site_summary_tbl = site_summary_tbl,
      export_k_values = QGIS_EXPORT_K_VALUES,
      out_dir = QGIS_DIR
    )
    
    site_diag <- build_optional_site_diagnostic(site_summary_tbl, structure_mean_q)
    if (!is.null(site_diag) && nrow(site_diag) > 0) {
      diag_file <- file.path(TABLES_DIR, "site_genetic_diagnostic_summary.csv")
      write.csv(site_diag, diag_file, row.names = FALSE)
      message("[03_structure] Saved optional site genetic diagnostic summary: ", diag_file)
    }
    
    message("[03_structure] Validation summary: all K plots passed validation (", validations_passed, " K values).")
    message("[03_structure] Validation summary: plotting data verified finite and within [0,1] before each geom_col call.")
    message("[03_structure] Validation summary: QGIS exports wrote one clean centroid-based site table per requested K.")
    message("[03_structure] Interpretation reminder: ", STRUCTURE_INTERPRETATION_NOTE)
    message("[03_structure] Done.")
  }
}