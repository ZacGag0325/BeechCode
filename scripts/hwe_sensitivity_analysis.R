# scripts/hwe_sensitivity_analysis.R
############################################################
# Sensitivity analysis: rerun core genetics summaries with
# suspect loci removed.
#
# Main fix:
# - FULL vs REDUCED heterozygosity/FIS comparison is truly
#   by site, not by locus.
# - Site order is forced to:
#   S1, S2, S3, S4, S5, S6, N1, N2, N3, N4, N5, N6
# - FULL diversity = clone-corrected gi_mll with all retained loci.
# - REDUCED diversity = clone-corrected gi_mll after removing:
#   EJV8T_A_0, ERHBI_A_0, FCM5, FG5.
############################################################

suppressPackageStartupMessages({
  library(adegenet)
  library(poppr)
  library(hierfstat)
  library(mmod)
  library(pegas)
  library(dplyr)
  library(tidyr)
  library(readxl)
  library(ggplot2)
  library(rlang)
})

source("scripts/_load_objects.R")

# ------------------------------------------------------------------
# Settings
# ------------------------------------------------------------------
suspect_null_loci <- c("EJV8T_A_0", "ERHBI_A_0", "FCM5", "FG5")
analysis_suffix <- "noSuspectLoci"

EXPECTED_SITE_LABELS <- c(
  "S1", "S2", "S3", "S4", "S5", "S6",
  "N1", "N2", "N3", "N4", "N5", "N6"
)

SITE_LOOKUP_SHEET <- "site_lookup"
HWE_MONTE_CARLO_REPS <- 9999L

SENS_TABLES_DIR <- file.path(TABLES_DIR, analysis_suffix)
SENS_FIGURES_DIR <- file.path(FIGURES_DIR, analysis_suffix)
SENS_MATRICES_DIR <- file.path(MATRICES_DIR, analysis_suffix)
COMPARISON_DIR <- file.path(TABLES_DIR, "comparisons")

for (d in c(SENS_TABLES_DIR, SENS_FIGURES_DIR, SENS_MATRICES_DIR, COMPARISON_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

message("[hwe_sensitivity] Starting reduced-loci sensitivity branch.")
message("[hwe_sensitivity] Suspect loci requested for removal: ", paste(suspect_null_loci, collapse = ", "))

# ------------------------------------------------------------------
# General helpers
# ------------------------------------------------------------------
normalize_lookup_key <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("\uFEFF", "", x, fixed = TRUE)
  x <- gsub("[[:cntrl:]]", "", x)
  x[is.na(x)] <- ""
  x
}

normalize_ascii <- function(x) {
  x <- iconv(as.character(x), from = "", to = "ASCII//TRANSLIT")
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  gsub("^_|_$", "", x)
}

pick_column <- function(df, choices, label = "column", required = FALSE) {
  idx <- match(TRUE, normalize_ascii(names(df)) %in% normalize_ascii(choices), nomatch = 0)
  if (idx == 0) {
    if (required) {
      stop(
        "[hwe_sensitivity] Could not find required ", label,
        ". Accepted names: ", paste(choices, collapse = ", "),
        call. = FALSE
      )
    }
    return(NA_character_)
  }
  names(df)[idx]
}

safe_mean <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (length(x) == 0 || all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

safe_se <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[!is.na(x)]
  if (length(x) <= 1) return(NA_real_)
  stats::sd(x) / sqrt(length(x))
}

contains_x_site_labels <- function(x) {
  any(grepl("^X[0-9]+$", as.character(x)))
}

write_csv_msg <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.csv(df, path, row.names = FALSE, na = "")
  message("[hwe_sensitivity] Saved: ", path)
}

force_expected_site_order <- function(df, site_col = "Site", context = "site table") {
  if (!site_col %in% names(df)) {
    stop("[hwe_sensitivity] ", context, " is missing column: ", site_col, call. = FALSE)
  }
  
  observed <- as.character(df[[site_col]])
  if (!setequal(observed, EXPECTED_SITE_LABELS)) {
    stop(
      "[hwe_sensitivity] ", context, " must contain exactly S1-S6 and N1-N6.",
      "\nExpected: ", paste(EXPECTED_SITE_LABELS, collapse = ", "),
      "\nObserved: ", paste(sort(unique(observed)), collapse = ", "),
      call. = FALSE
    )
  }
  
  df %>%
    mutate("{site_col}" := factor(.data[[site_col]], levels = EXPECTED_SITE_LABELS)) %>%
    arrange(.data[[site_col]]) %>%
    mutate("{site_col}" := as.character(.data[[site_col]]))
}

validate_site_table <- function(df, context = "site table") {
  df <- force_expected_site_order(df, "Site", context)
  
  if (nrow(df) != length(EXPECTED_SITE_LABELS)) {
    stop(
      "[hwe_sensitivity] ", context, " must have 12 rows, but has ",
      nrow(df), ".",
      call. = FALSE
    )
  }
  
  if (!identical(as.character(df$Site), EXPECTED_SITE_LABELS)) {
    stop(
      "[hwe_sensitivity] ", context, " site labels/order are incorrect.",
      "\nExpected: ", paste(EXPECTED_SITE_LABELS, collapse = ", "),
      "\nObserved: ", paste(as.character(df$Site), collapse = ", "),
      call. = FALSE
    )
  }
  
  if (contains_x_site_labels(df$Site)) {
    stop(
      "[hwe_sensitivity] ", context,
      " contains locus-like X labels; refusing to write output.",
      call. = FALSE
    )
  }
  
  df
}

# ------------------------------------------------------------------
# site_lookup loading and label mapping
# ------------------------------------------------------------------
find_site_lookup_workbook <- function(raw_dir = file.path(PROJECT_ROOT, "data", "raw"),
                                      sheet = SITE_LOOKUP_SHEET) {
  if (!dir.exists(raw_dir)) {
    message("[hwe_sensitivity] site_lookup unavailable because data/raw does not exist: ", raw_dir)
    return(NULL)
  }
  
  excel_files <- list.files(
    raw_dir,
    pattern = "\\.(xlsx|xls)$",
    full.names = TRUE,
    ignore.case = TRUE
  )
  
  if (length(excel_files) == 0) {
    message("[hwe_sensitivity] site_lookup unavailable: no Excel workbook was found in data/raw.")
    return(NULL)
  }
  
  has_lookup <- vapply(
    excel_files,
    function(path) {
      sheets <- tryCatch(readxl::excel_sheets(path), error = function(e) character(0))
      any(normalize_ascii(sheets) == normalize_ascii(sheet))
    },
    logical(1)
  )
  
  lookup_files <- excel_files[has_lookup]
  
  if (length(lookup_files) == 0) {
    message("[hwe_sensitivity] site_lookup unavailable: no workbook in data/raw contains a '", sheet, "' sheet.")
    return(NULL)
  }
  
  if (length(lookup_files) > 1) {
    message("[hwe_sensitivity] Multiple workbooks contain site_lookup; using ", basename(lookup_files[1]), ".")
  } else {
    message("[hwe_sensitivity] Reading site_lookup from: ", basename(lookup_files[1]))
  }
  
  lookup_files[1]
}

load_site_lookup <- function() {
  workbook <- find_site_lookup_workbook()
  if (is.null(workbook)) return(NULL)
  
  sheets <- readxl::excel_sheets(workbook)
  sheet <- sheets[match(normalize_ascii(SITE_LOOKUP_SHEET), normalize_ascii(sheets))]
  
  lookup <- suppressMessages(readxl::read_excel(workbook, sheet = sheet)) %>%
    as.data.frame(stringsAsFactors = FALSE)
  
  old_site_col <- pick_column(
    lookup,
    c("Site", "site", "site_code", "old_site", "code_site"),
    "old site code",
    required = TRUE
  )
  
  label_col <- pick_column(
    lookup,
    c("Site_label", "site_label", "new_site", "site_new", "label", "site_id"),
    "display site label"
  )
  
  order_col <- pick_column(
    lookup,
    c("Site_order", "site_order", "order", "ordre", "sort", "south_north_order"),
    "site order"
  )
  
  if (is.na(label_col)) label_col <- old_site_col
  
  out <- lookup %>%
    mutate(
      old_site = normalize_lookup_key(.data[[old_site_col]]),
      site_label = normalize_lookup_key(.data[[label_col]]),
      site_order = if (!is.na(order_col)) suppressWarnings(as.numeric(.data[[order_col]])) else NA_real_
    ) %>%
    filter(nzchar(old_site), nzchar(site_label)) %>%
    distinct(old_site, .keep_all = TRUE)
  
  if (nrow(out) == 0) {
    message("[hwe_sensitivity] site_lookup was found but no usable rows were detected; raw site labels will be retained.")
    return(NULL)
  }
  
  if (anyDuplicated(out$old_site)) {
    stop("[hwe_sensitivity] site_lookup contains duplicated old site codes after cleaning.", call. = FALSE)
  }
  
  if (anyDuplicated(out$site_label)) {
    stop("[hwe_sensitivity] site_lookup contains duplicated display site labels after cleaning.", call. = FALSE)
  }
  
  message("[hwe_sensitivity] Loaded site_lookup rows: ", nrow(out))
  out
}

site_lookup <- load_site_lookup()

map_site_labels <- function(site_values) {
  site_values <- normalize_lookup_key(site_values)
  if (is.null(site_lookup)) return(site_values)
  
  idx <- match(site_values, site_lookup$old_site)
  labels <- site_lookup$site_label[idx]
  
  ifelse(is.na(labels) | !nzchar(labels), site_values, labels)
}

get_clone_corrected_site_labels <- function(gobj) {
  validate_columns(df_ids_mll, c("ind_id", "Site"), df_name = "[hwe_sensitivity] df_ids_mll")
  
  if (!all(adegenet::indNames(gobj) == df_ids_mll$ind_id)) {
    stop("[hwe_sensitivity] Clone-corrected genind object is not aligned with df_ids_mll.", call. = FALSE)
  }
  
  raw_sites <- normalize_lookup_key(df_ids_mll$Site)
  
  if (any(!nzchar(raw_sites))) {
    bad_ids <- df_ids_mll$ind_id[!nzchar(raw_sites)]
    stop(
      "[hwe_sensitivity] df_ids_mll contains missing/blank Site labels for: ",
      paste(head(bad_ids, 10), collapse = ", "),
      call. = FALSE
    )
  }
  
  if (contains_x_site_labels(raw_sites)) {
    stop("[hwe_sensitivity] Raw df_ids_mll site labels contain locus-like X labels; refusing to continue.", call. = FALSE)
  }
  
  site_labels <- map_site_labels(raw_sites)
  
  if (contains_x_site_labels(site_labels)) {
    stop("[hwe_sensitivity] Final mapped site labels contain locus-like X labels; refusing to continue.", call. = FALSE)
  }
  
  if (!setequal(unique(site_labels), EXPECTED_SITE_LABELS)) {
    stop(
      "[hwe_sensitivity] Final clone-corrected site labels must be exactly S1-S6 and N1-N6.",
      "\nExpected: ", paste(EXPECTED_SITE_LABELS, collapse = ", "),
      "\nObserved: ", paste(sort(unique(site_labels)), collapse = ", "),
      call. = FALSE
    )
  }
  
  site_labels
}

# ------------------------------------------------------------------
# Genind reduction helper
# ------------------------------------------------------------------
make_reduced_genind <- function(gobj, loci_to_remove) {
  loci_available <- adegenet::locNames(gobj)
  
  requested <- unique(trimws(as.character(loci_to_remove)))
  requested <- requested[nzchar(requested)]
  
  found <- intersect(requested, loci_available)
  missing <- setdiff(requested, loci_available)
  retained <- setdiff(loci_available, found)
  
  if (length(retained) < 2) {
    stop("[hwe_sensitivity] Too few loci retained after filtering (<2).", call. = FALSE)
  }
  
  if (length(missing) > 0) {
    warning(
      "[hwe_sensitivity] Requested locus/loci not found and therefore not removed: ",
      paste(missing, collapse = ", ")
    )
  }
  
  if (length(found) == 0) {
    warning("[hwe_sensitivity] None of the requested loci were found; reduced dataset equals full dataset.")
  }
  
  loci_df <- adegenet::genind2df(gobj, sep = "/", usepop = FALSE)
  keep_df <- intersect(retained, names(loci_df))
  
  if (length(keep_df) != length(retained)) {
    missing_keep <- setdiff(retained, names(loci_df))
    stop(
      "[hwe_sensitivity] Failed to retain all expected loci when rebuilding genind. Missing: ",
      paste(missing_keep, collapse = ", "),
      call. = FALSE
    )
  }
  
  rebuilt <- adegenet::df2genind(
    X = loci_df[, keep_df, drop = FALSE],
    sep = "/",
    ploidy = adegenet::ploidy(gobj),
    ind.names = adegenet::indNames(gobj),
    type = if (!is.null(gobj@type) && length(gobj@type) > 0) gobj@type else "codom",
    NA.char = "NA"
  )
  
  adegenet::pop(rebuilt) <- adegenet::pop(gobj)
  
  if (adegenet::nLoc(rebuilt) != length(retained)) {
    stop(
      "[hwe_sensitivity] Locus subsetting mismatch: expected ",
      length(retained), " got ", adegenet::nLoc(rebuilt), ".",
      call. = FALSE
    )
  }
  
  message("[hwe_sensitivity] Loci successfully removed: ", if (length(found) > 0) paste(found, collapse = ", ") else "none")
  message("[hwe_sensitivity] Loci not found: ", if (length(missing) > 0) paste(missing, collapse = ", ") else "none")
  message("[hwe_sensitivity] Locus count before filtering: ", length(loci_available))
  message("[hwe_sensitivity] Locus count after filtering: ", length(retained))
  
  list(
    reduced = rebuilt,
    found = found,
    missing = missing,
    retained = retained,
    all_full = loci_available
  )
}

# ------------------------------------------------------------------
# Diversity helpers
# ------------------------------------------------------------------
extract_basic_stat_by_site <- function(stat_matrix, site_labels, stat_name) {
  mat <- as.matrix(stat_matrix)
  storage.mode(mat) <- "numeric"
  
  if (all(site_labels %in% colnames(mat))) {
    values <- vapply(site_labels, function(site) safe_mean(mat[, site]), numeric(1))
  } else if (all(site_labels %in% rownames(mat))) {
    values <- vapply(site_labels, function(site) safe_mean(mat[site, ]), numeric(1))
  } else {
    stop(
      "[hwe_sensitivity] Could not locate site labels in basic.stats matrix for ", stat_name, ".",
      " Row names: ", paste(head(rownames(mat), 20), collapse = ", "),
      ". Column names: ", paste(head(colnames(mat), 20), collapse = ", "),
      call. = FALSE
    )
  }
  
  data.frame(Site = names(values), value = as.numeric(values), stringsAsFactors = FALSE)
}

extract_basic_stat_by_locus <- function(stat_matrix, stat_name, site_labels) {
  mat <- as.matrix(stat_matrix)
  storage.mode(mat) <- "numeric"
  
  if (all(site_labels %in% colnames(mat))) {
    values <- apply(mat, 1, safe_mean)
    loci <- rownames(mat)
  } else if (all(site_labels %in% rownames(mat))) {
    values <- apply(mat, 2, safe_mean)
    loci <- colnames(mat)
  } else {
    values <- numeric(0)
    loci <- character(0)
  }
  
  data.frame(Locus = loci, value = as.numeric(values), stringsAsFactors = FALSE) %>%
    rename(!!stat_name := value)
}

extract_ar_by_site <- function(ar_matrix, site_labels) {
  mat <- as.matrix(ar_matrix)
  storage.mode(mat) <- "numeric"
  
  if (all(site_labels %in% colnames(mat))) {
    out <- vapply(site_labels, function(site) safe_mean(mat[, site]), numeric(1))
    out_se <- vapply(site_labels, function(site) safe_se(mat[, site]), numeric(1))
  } else if (all(site_labels %in% rownames(mat))) {
    out <- vapply(site_labels, function(site) safe_mean(mat[site, ]), numeric(1))
    out_se <- vapply(site_labels, function(site) safe_se(mat[site, ]), numeric(1))
  } else {
    stop(
      "[hwe_sensitivity] Could not locate site labels in allelic.richness Ar matrix.",
      " Row names: ", paste(head(rownames(mat), 20), collapse = ", "),
      ". Column names: ", paste(head(colnames(mat), 20), collapse = ", "),
      call. = FALSE
    )
  }
  
  data.frame(
    Site = names(out),
    Allelic_Richness = as.numeric(out),
    Allelic_Richness_SE = as.numeric(out_se),
    stringsAsFactors = FALSE
  )
}

compute_na_by_site <- function(gobj, site_order) {
  allele_tab <- adegenet::tab(gobj, NA.method = "asis")
  loc_fac <- adegenet::locFac(gobj)
  loci <- adegenet::locNames(gobj)
  
  site_row_indices <- split(seq_len(nrow(allele_tab)), as.character(adegenet::pop(gobj)))
  site_row_indices <- site_row_indices[site_order]
  
  count_alleles_for_site_locus <- function(rows, locus_name) {
    cols <- which(loc_fac == locus_name)
    if (length(cols) == 0 || length(rows) == 0) return(NA_integer_)
    mat <- allele_tab[rows, cols, drop = FALSE]
    if (all(is.na(mat))) return(NA_integer_)
    as.integer(sum(colSums(mat, na.rm = TRUE) > 0))
  }
  
  expand.grid(Site = site_order, Locus = loci, stringsAsFactors = FALSE) %>%
    rowwise() %>%
    mutate(Allele_count = count_alleles_for_site_locus(site_row_indices[[Site]], Locus)) %>%
    ungroup() %>%
    group_by(Site) %>%
    summarise(Na = safe_mean(Allele_count), .groups = "drop")
}

extract_overall_stat <- function(overall_obj, stat_name) {
  if (is.null(overall_obj)) return(NA_real_)
  
  if (is.atomic(overall_obj) && !is.null(names(overall_obj))) {
    val <- overall_obj[[stat_name]]
    return(if (length(val) == 0) NA_real_ else as.numeric(val[1]))
  }
  
  if (is.matrix(overall_obj) || is.data.frame(overall_obj)) {
    rn <- rownames(overall_obj)
    cn <- colnames(overall_obj)
    
    if (!is.null(rn) && stat_name %in% rn) {
      return(safe_mean(as.numeric(overall_obj[stat_name, , drop = TRUE])))
    }
    
    if (!is.null(cn) && stat_name %in% cn) {
      return(safe_mean(as.numeric(overall_obj[, stat_name, drop = TRUE])))
    }
  }
  
  NA_real_
}

run_diversity_bundle <- function(gobj, dataset_label) {
  site_order <- EXPECTED_SITE_LABELS
  observed_sites <- unique(as.character(adegenet::pop(gobj)))
  
  if (!setequal(observed_sites, site_order)) {
    stop(
      "[hwe_sensitivity] ", dataset_label, " pop labels are not the expected S1-S6/N1-N6 set.",
      "\nExpected: ", paste(site_order, collapse = ", "),
      "\nObserved: ", paste(sort(observed_sites), collapse = ", "),
      call. = FALSE
    )
  }
  
  hf <- hierfstat::genind2hierfstat(gobj)
  bs <- hierfstat::basic.stats(hf)
  
  if (is.null(bs$Ho) || is.null(bs$Hs) || is.null(bs$Fis)) {
    stop(
      "[hwe_sensitivity] hierfstat::basic.stats did not return Ho, Hs, and Fis matrices for ",
      dataset_label, ".",
      call. = FALSE
    )
  }
  
  ho_by_site <- extract_basic_stat_by_site(bs$Ho, site_order, "Ho") %>%
    rename(Ho = value)
  
  he_by_site <- extract_basic_stat_by_site(bs$Hs, site_order, "He") %>%
    rename(He = value)
  
  fis_by_site <- extract_basic_stat_by_site(bs$Fis, site_order, "FIS") %>%
    rename(FIS = value)
  
  pop_sizes <- table(as.character(adegenet::pop(gobj)))
  pop_sizes <- pop_sizes[site_order]
  
  if (any(is.na(pop_sizes)) || any(pop_sizes <= 0)) {
    stop(
      "[hwe_sensitivity] ", dataset_label,
      " N per site contains missing or zero values after site assignment.",
      call. = FALSE
    )
  }
  
  n_by_site <- data.frame(
    Site = names(pop_sizes),
    N = as.integer(pop_sizes),
    stringsAsFactors = FALSE
  )
  
  na_by_site <- compute_na_by_site(gobj, site_order)
  
  min_n <- min(as.integer(pop_sizes), na.rm = TRUE)
  ar <- if (!is.na(min_n) && min_n >= 2) {
    hierfstat::allelic.richness(hf, min.n = min_n)
  } else {
    NULL
  }
  
  if (is.null(ar) || is.null(ar$Ar)) {
    ar_by_site <- data.frame(
      Site = site_order,
      Allelic_Richness = NA_real_,
      Allelic_Richness_SE = NA_real_,
      stringsAsFactors = FALSE
    )
  } else {
    ar_by_site <- extract_ar_by_site(ar$Ar, site_order)
  }
  
  ar_by_site <- validate_site_table(ar_by_site, paste0(dataset_label, " allelic richness"))
  
  by_site <- n_by_site %>%
    left_join(na_by_site, by = "Site") %>%
    left_join(ar_by_site, by = "Site") %>%
    left_join(ho_by_site, by = "Site") %>%
    left_join(he_by_site, by = "Site") %>%
    left_join(fis_by_site, by = "Site") %>%
    mutate(Dataset = dataset_label) %>%
    select(Dataset, Site, N, Na, Ho, He, FIS, Allelic_Richness, Allelic_Richness_SE) %>%
    validate_site_table(paste0(dataset_label, " heterozygosity/FIS by site"))
  
  by_locus <- extract_basic_stat_by_locus(bs$Ho, "Ho", site_order) %>%
    full_join(extract_basic_stat_by_locus(bs$Hs, "He", site_order), by = "Locus") %>%
    full_join(extract_basic_stat_by_locus(bs$Fis, "FIS", site_order), by = "Locus") %>%
    mutate(Dataset = dataset_label) %>%
    select(Dataset, Locus, Ho, He, FIS)
  
  overall_ho <- extract_overall_stat(bs$overall, "Ho")
  overall_he <- extract_overall_stat(bs$overall, "Hs")
  overall_fis <- extract_overall_stat(bs$overall, "Fis")
  
  if (is.na(overall_ho)) overall_ho <- safe_mean(by_site$Ho)
  if (is.na(overall_he)) overall_he <- safe_mean(by_site$He)
  if (is.na(overall_fis)) overall_fis <- safe_mean(by_site$FIS)
  
  overall <- data.frame(
    Dataset = dataset_label,
    N = adegenet::nInd(gobj),
    N_loci = adegenet::nLoc(gobj),
    Ho = overall_ho,
    He = overall_he,
    FIS = overall_fis,
    stringsAsFactors = FALSE
  )
  
  ar_by_site <- ar_by_site %>%
    mutate(Dataset = dataset_label) %>%
    select(Dataset, Site, Allelic_Richness, Allelic_Richness_SE) %>%
    validate_site_table(paste0(dataset_label, " allelic richness by site"))
  
  list(
    by_site = by_site,
    by_locus = by_locus,
    overall = overall,
    allelic_richness = ar_by_site
  )
}

# ------------------------------------------------------------------
# HWE, differentiation, PCA/DAPC helper summaries
# ------------------------------------------------------------------
matrix_to_long_unique <- function(m, value_name) {
  long_all <- as.data.frame(as.table(m), stringsAsFactors = FALSE)
  names(long_all) <- c("Site1", "Site2", value_name)
  
  long_all %>%
    filter(Site1 != Site2) %>%
    mutate(site_pair = paste(pmin(Site1, Site2), pmax(Site1, Site2), sep = "__")) %>%
    group_by(site_pair) %>%
    summarise(
      Site1 = first(pmin(Site1, Site2)),
      Site2 = first(pmax(Site1, Site2)),
      value = mean(.data[[value_name]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rename(!!value_name := value) %>%
    select(Site1, Site2, all_of(value_name))
}

run_differentiation_bundle <- function(gobj, dataset_label) {
  jost_mat <- as.matrix(mmod::pairwise_D(gobj, linearized = FALSE))
  diag(jost_mat) <- 0
  
  hf <- hierfstat::genind2hierfstat(gobj)
  if (!is.numeric(hf[[1]])) hf[[1]] <- as.integer(factor(hf[[1]]))
  
  fst_raw <- tryCatch(
    hierfstat::pairwise.WCfst(hf),
    error = function(e) hierfstat::pairwise.neifst(hf)
  )
  
  fst_mat <- as.matrix(fst_raw)
  diag(fst_mat) <- 0
  
  site_lookup_tmp <- tapply(
    as.character(adegenet::pop(gobj)),
    hf[[1]],
    function(x) names(sort(table(x), decreasing = TRUE))[1]
  )
  
  site_lookup_tmp <- as.character(site_lookup_tmp)
  
  if (length(site_lookup_tmp) == nrow(fst_mat)) {
    rownames(fst_mat) <- site_lookup_tmp
    colnames(fst_mat) <- site_lookup_tmp
  }
  
  jost_long <- matrix_to_long_unique(jost_mat, "JostD") %>%
    mutate(Dataset = dataset_label)
  
  fst_long <- matrix_to_long_unique(fst_mat, "FST") %>%
    mutate(Dataset = dataset_label)
  
  summary_tbl <- data.frame(
    Dataset = dataset_label,
    mean_pairwise_JostD = mean(jost_long$JostD, na.rm = TRUE),
    median_pairwise_JostD = median(jost_long$JostD, na.rm = TRUE),
    mean_pairwise_FST = mean(fst_long$FST, na.rm = TRUE),
    median_pairwise_FST = median(fst_long$FST, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  
  list(
    jost_mat = jost_mat,
    fst_mat = fst_mat,
    jost_long = jost_long,
    fst_long = fst_long,
    summary = summary_tbl
  )
}

run_hwe_by_site_locus <- function(gobj, dataset_label) {
  site_vec <- as.character(adegenet::pop(gobj))
  loci_df <- adegenet::genind2df(gobj, sep = "/", usepop = FALSE)
  loci_names <- names(loci_df)
  
  all_rows <- list()
  idx <- 1L
  
  for (site in EXPECTED_SITE_LABELS) {
    use <- site_vec == site
    if (sum(use) < 2) next
    
    site_df <- loci_df[use, , drop = FALSE]
    
    for (loc in loci_names) {
      vals <- trimws(as.character(site_df[[loc]]))
      vals[vals %in% c("", "NA", "0", "0/0", "NA/NA", "-")] <- NA_character_
      vals <- vals[!is.na(vals)]
      if (length(vals) < 2) next
      
      geno_fac <- factor(vals)
      
      p_val <- tryCatch(
        pegas::hw.test(geno_fac, B = HWE_MONTE_CARLO_REPS)$p.value,
        error = function(e) NA_real_
      )
      
      all_rows[[idx]] <- data.frame(
        Dataset = dataset_label,
        Site = site,
        Locus = loc,
        N_non_missing = length(vals),
        N_genotype_classes = nlevels(geno_fac),
        p_value_raw = as.numeric(p_val),
        stringsAsFactors = FALSE
      )
      
      idx <- idx + 1L
    }
  }
  
  out <- bind_rows(all_rows)
  if (nrow(out) == 0) return(out)
  
  out %>%
    mutate(
      Site = factor(Site, levels = EXPECTED_SITE_LABELS),
      p_value_adj_fdr = p.adjust(p_value_raw, method = "BH"),
      hwe_reject_raw_0_05 = !is.na(p_value_raw) & p_value_raw < 0.05,
      hwe_reject_fdr_0_05 = !is.na(p_value_adj_fdr) & p_value_adj_fdr < 0.05
    ) %>%
    arrange(Site, Locus) %>%
    mutate(Site = as.character(Site))
}

run_pca_bundle <- function(gobj, dataset_label) {
  X <- adegenet::tab(gobj, freq = TRUE, NA.method = "mean")
  keep_cols <- apply(X, 2, function(v) stats::var(v, na.rm = TRUE) > 0)
  X <- X[, keep_cols, drop = FALSE]
  
  pca_fit <- stats::prcomp(X, center = TRUE, scale. = FALSE)
  var_exp <- (pca_fit$sdev^2 / sum(pca_fit$sdev^2)) * 100
  
  scores <- as.data.frame(pca_fit$x[, 1:2, drop = FALSE], stringsAsFactors = FALSE)
  scores$Individual <- rownames(scores)
  scores$Site <- as.character(adegenet::pop(gobj))
  scores$Dataset <- dataset_label
  
  variance <- data.frame(
    Dataset = dataset_label,
    PC = paste0("PC", seq_along(var_exp)),
    Percent_Variance = as.numeric(var_exp),
    stringsAsFactors = FALSE
  )
  
  list(scores = scores, variance = variance)
}

run_dapc_bundle <- function(gobj, dataset_label) {
  grp <- as.factor(adegenet::pop(gobj))
  
  dapc_fit <- adegenet::dapc(
    gobj,
    pop = grp,
    n.pca = min(50, adegenet::nInd(gobj) - 1),
    n.da = min(nlevels(grp) - 1, 10)
  )
  
  coords <- as.data.frame(dapc_fit$ind.coord[, 1:2, drop = FALSE], stringsAsFactors = FALSE)
  coords$Individual <- rownames(coords)
  coords$Site <- as.character(grp)
  coords$Dataset <- dataset_label
  
  if (ncol(coords) >= 2) names(coords)[1:2] <- c("LD1", "LD2")
  
  eig <- as.numeric(dapc_fit$eig)
  
  eig_df <- data.frame(
    Dataset = dataset_label,
    Axis = paste0("LD", seq_along(eig)),
    Eigenvalue = eig,
    stringsAsFactors = FALSE
  )
  
  list(coords = coords, eig = eig_df)
}

# ------------------------------------------------------------------
# Comparison helpers
# ------------------------------------------------------------------
compare_two_dataset_table <- function(df, by_cols, value_cols) {
  wide <- df %>%
    select(Dataset, all_of(by_cols), all_of(value_cols)) %>%
    pivot_wider(names_from = Dataset, values_from = all_of(value_cols), names_sep = "__")
  
  for (v in value_cols) {
    full_col <- paste0(v, "__FULL")
    red_col <- paste0(v, "__REDUCED")
    delta_col <- paste0("delta_", v, "_REDUCED_minus_FULL")
    
    if (all(c(full_col, red_col) %in% names(wide))) {
      wide[[delta_col]] <- suppressWarnings(as.numeric(wide[[red_col]]) - as.numeric(wide[[full_col]]))
    }
  }
  
  wide
}

compare_two_dataset_table_clean <- function(df, by_cols, value_cols) {
  wide <- df %>%
    select(Dataset, all_of(by_cols), all_of(value_cols)) %>%
    pivot_wider(names_from = Dataset, values_from = all_of(value_cols), names_sep = "__")
  
  for (v in value_cols) {
    full_col <- paste0(v, "__FULL")
    red_col <- paste0(v, "__REDUCED")
    
    if (all(c(full_col, red_col) %in% names(wide))) {
      wide[[paste0(v, "_delta")]] <- suppressWarnings(as.numeric(wide[[red_col]]) - as.numeric(wide[[full_col]]))
    }
  }
  
  names(wide) <- gsub("__FULL$", "_full", names(wide))
  names(wide) <- gsub("__REDUCED$", "_reduced", names(wide))
  
  if ("Site" %in% names(wide) && setequal(as.character(wide$Site), EXPECTED_SITE_LABELS)) {
    wide <- force_expected_site_order(wide, "Site", "comparison table")
  }
  
  ordered_cols <- c(
    by_cols,
    unlist(lapply(
      value_cols,
      function(v) c(paste0(v, "_full"), paste0(v, "_reduced"), paste0(v, "_delta"))
    ))
  )
  
  ordered_cols <- ordered_cols[ordered_cols %in% names(wide)]
  
  wide %>%
    select(all_of(ordered_cols), everything())
}

validate_by_site_comparison <- function(tbl, context = "by-site comparison") {
  required_cols <- c(
    "Site",
    "N_full", "N_reduced",
    "Ho_full", "Ho_reduced", "Ho_delta",
    "He_full", "He_reduced", "He_delta",
    "FIS_full", "FIS_reduced", "FIS_delta"
  )
  
  missing_cols <- setdiff(required_cols, names(tbl))
  
  if (length(missing_cols) > 0) {
    stop(
      "[hwe_sensitivity] ", context,
      " is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
  
  tbl <- validate_site_table(tbl, context)
  
  invisible(tbl)
}

# ------------------------------------------------------------------
# Build clone-corrected full and reduced diversity objects
# ------------------------------------------------------------------
full_gi <- gi_mll

reduced_info <- make_reduced_genind(full_gi, suspect_null_loci)
reduced_gi <- reduced_info$reduced

clone_corrected_site_labels <- get_clone_corrected_site_labels(full_gi)

adegenet::pop(full_gi) <- factor(clone_corrected_site_labels, levels = EXPECTED_SITE_LABELS)

if (!all(adegenet::indNames(reduced_gi) == adegenet::indNames(full_gi))) {
  stop("[hwe_sensitivity] Reduced clone-corrected genind object is not aligned with full_gi.", call. = FALSE)
}

adegenet::pop(reduced_gi) <- factor(clone_corrected_site_labels, levels = EXPECTED_SITE_LABELS)

message("[hwe_sensitivity] Diversity datasets are clone-corrected: FULL = gi_mll; REDUCED = gi_mll without suspect loci.")
message("[hwe_sensitivity] Site labels/order for diversity summaries: ", paste(EXPECTED_SITE_LABELS, collapse = ", "))

# ------------------------------------------------------------------
# Locus manifest
# ------------------------------------------------------------------
locus_manifest <- data.frame(
  dataset = c("FULL", "REDUCED"),
  n_loci = c(adegenet::nLoc(full_gi), adegenet::nLoc(reduced_gi)),
  loci = c(
    paste(adegenet::locNames(full_gi), collapse = ";"),
    paste(adegenet::locNames(reduced_gi), collapse = ";")
  ),
  removed_loci = c("", paste(reduced_info$found, collapse = ";")),
  requested_but_not_found = c("", paste(reduced_info$missing, collapse = ";")),
  stringsAsFactors = FALSE
)

write_csv_msg(
  locus_manifest,
  file.path(SENS_TABLES_DIR, paste0("locus_manifest_", analysis_suffix, ".csv"))
)

# ------------------------------------------------------------------
# Run full/reduced summaries
# ------------------------------------------------------------------
full_div <- run_diversity_bundle(full_gi, "FULL")
red_div <- run_diversity_bundle(reduced_gi, "REDUCED")

full_diff <- run_differentiation_bundle(full_gi, "FULL")
red_diff <- run_differentiation_bundle(reduced_gi, "REDUCED")

full_hwe <- run_hwe_by_site_locus(full_gi, "FULL")
red_hwe <- run_hwe_by_site_locus(reduced_gi, "REDUCED")

full_pca <- run_pca_bundle(full_gi, "FULL")
red_pca <- run_pca_bundle(reduced_gi, "REDUCED")

full_dapc <- run_dapc_bundle(full_gi, "FULL")
red_dapc <- run_dapc_bundle(reduced_gi, "REDUCED")

# ------------------------------------------------------------------
# Write reduced-only outputs
# ------------------------------------------------------------------
red_div$by_site <- validate_site_table(
  red_div$by_site,
  "heterozygosity_fis_by_site_noSuspectLoci"
)

red_div$allelic_richness <- validate_site_table(
  red_div$allelic_richness,
  "allelic_richness_by_site_noSuspectLoci"
)

write_csv_msg(
  red_div$by_site,
  file.path(SENS_TABLES_DIR, paste0("heterozygosity_fis_by_site_", analysis_suffix, ".csv"))
)

write_csv_msg(
  red_div$overall,
  file.path(SENS_TABLES_DIR, paste0("heterozygosity_fis_overall_", analysis_suffix, ".csv"))
)

write_csv_msg(
  red_div$allelic_richness,
  file.path(SENS_TABLES_DIR, paste0("allelic_richness_by_site_", analysis_suffix, ".csv"))
)

write_csv_msg(
  red_diff$fst_long,
  file.path(SENS_TABLES_DIR, paste0("pairwise_fst_", analysis_suffix, ".csv"))
)

write_csv_msg(
  red_diff$jost_long,
  file.path(SENS_TABLES_DIR, paste0("pairwise_jostD_", analysis_suffix, ".csv"))
)

write_csv_msg(
  red_hwe,
  file.path(SENS_TABLES_DIR, paste0("hwe_by_site_by_locus_", analysis_suffix, ".csv"))
)

write_csv_msg(
  red_pca$variance,
  file.path(SENS_TABLES_DIR, paste0("pca_variance_", analysis_suffix, ".csv"))
)

write_csv_msg(
  red_dapc$eig,
  file.path(SENS_TABLES_DIR, paste0("dapc_eigenvalues_", analysis_suffix, ".csv"))
)

write.csv(
  red_diff$fst_mat,
  file.path(SENS_MATRICES_DIR, paste0("pairwise_fst_", analysis_suffix, ".csv"))
)

write.csv(
  red_diff$jost_mat,
  file.path(SENS_MATRICES_DIR, paste0("pairwise_jostD_", analysis_suffix, ".csv"))
)

# ------------------------------------------------------------------
# FULL vs REDUCED comparisons
# ------------------------------------------------------------------
div_by_site_all <- bind_rows(full_div$by_site, red_div$by_site)
div_by_locus_all <- bind_rows(full_div$by_locus, red_div$by_locus)
div_overall_all <- bind_rows(full_div$overall, red_div$overall)
ar_all <- bind_rows(full_div$allelic_richness, red_div$allelic_richness)

diff_summary_all <- bind_rows(full_diff$summary, red_diff$summary)
fst_long_all <- bind_rows(full_diff$fst_long, red_diff$fst_long)
jost_long_all <- bind_rows(full_diff$jost_long, red_diff$jost_long)

hwe_all <- bind_rows(full_hwe, red_hwe)
pca_var_all <- bind_rows(full_pca$variance, red_pca$variance)
dapc_eig_all <- bind_rows(full_dapc$eig, red_dapc$eig)

summary_stats <- data.frame(
  Dataset = c("FULL", "REDUCED"),
  N_individuals = c(adegenet::nInd(full_gi), adegenet::nInd(reduced_gi)),
  N_loci = c(adegenet::nLoc(full_gi), adegenet::nLoc(reduced_gi)),
  stringsAsFactors = FALSE
)

summary_stats_cmp <- compare_two_dataset_table(
  summary_stats,
  by_cols = character(0),
  value_cols = c("N_individuals", "N_loci")
)

div_overall_cmp <- compare_two_dataset_table(
  div_overall_all,
  by_cols = character(0),
  value_cols = c("N", "N_loci", "Ho", "He", "FIS")
)

div_site_cmp <- compare_two_dataset_table_clean(
  div_by_site_all,
  by_cols = c("Site"),
  value_cols = c(
    "N",
    "Ho",
    "He",
    "FIS",
    "Na",
    "Allelic_Richness",
    "Allelic_Richness_SE"
  )
)

div_site_cmp <- validate_site_table(
  div_site_cmp,
  "full_vs_noSuspectLoci_heterozygosity_fis_by_site_comparison"
)

validate_by_site_comparison(
  div_site_cmp,
  "full_vs_noSuspectLoci_heterozygosity_fis_by_site_comparison"
)

div_locus_cmp <- compare_two_dataset_table_clean(
  div_by_locus_all,
  by_cols = c("Locus"),
  value_cols = c("Ho", "He", "FIS")
)

ar_cmp <- compare_two_dataset_table_clean(
  ar_all,
  by_cols = c("Site"),
  value_cols = c("Allelic_Richness", "Allelic_Richness_SE")
)

ar_cmp <- validate_site_table(
  ar_cmp,
  "full_vs_noSuspectLoci_allelic_richness_comparison"
)

diff_cmp <- compare_two_dataset_table(
  diff_summary_all,
  by_cols = character(0),
  value_cols = c(
    "mean_pairwise_JostD",
    "median_pairwise_JostD",
    "mean_pairwise_FST",
    "median_pairwise_FST"
  )
)

pca_cmp <- compare_two_dataset_table(
  pca_var_all,
  by_cols = c("PC"),
  value_cols = c("Percent_Variance")
)

dapc_cmp <- compare_two_dataset_table(
  dapc_eig_all,
  by_cols = c("Axis"),
  value_cols = c("Eigenvalue")
)

write_csv_msg(
  summary_stats_cmp,
  file.path(COMPARISON_DIR, "full_vs_noSuspectLoci_locus_count_comparison.csv")
)

write_csv_msg(
  div_overall_cmp,
  file.path(COMPARISON_DIR, "full_vs_noSuspectLoci_heterozygosity_fis_overall_comparison.csv")
)

write_csv_msg(
  div_site_cmp,
  file.path(COMPARISON_DIR, "full_vs_noSuspectLoci_heterozygosity_fis_by_site_comparison.csv")
)

write_csv_msg(
  div_locus_cmp,
  file.path(COMPARISON_DIR, "full_vs_noSuspectLoci_heterozygosity_fis_by_locus_comparison.csv")
)

write_csv_msg(
  ar_cmp,
  file.path(COMPARISON_DIR, "full_vs_noSuspectLoci_allelic_richness_comparison.csv")
)

write_csv_msg(
  diff_cmp,
  file.path(COMPARISON_DIR, "full_vs_noSuspectLoci_differentiation_comparison.csv")
)

write_csv_msg(
  pca_cmp,
  file.path(COMPARISON_DIR, "full_vs_noSuspectLoci_pca_variance_comparison.csv")
)

write_csv_msg(
  dapc_cmp,
  file.path(COMPARISON_DIR, "full_vs_noSuspectLoci_dapc_comparison.csv")
)

write_csv_msg(
  fst_long_all,
  file.path(COMPARISON_DIR, "full_vs_noSuspectLoci_pairwise_fst_long.csv")
)

write_csv_msg(
  jost_long_all,
  file.path(COMPARISON_DIR, "full_vs_noSuspectLoci_pairwise_jostD_long.csv")
)

write_csv_msg(
  hwe_all,
  file.path(COMPARISON_DIR, "full_vs_noSuspectLoci_hwe_by_site_by_locus_long.csv")
)

write_csv_msg(
  pca_var_all,
  file.path(COMPARISON_DIR, "full_vs_noSuspectLoci_pca_variance_long.csv")
)

write_csv_msg(
  dapc_eig_all,
  file.path(COMPARISON_DIR, "full_vs_noSuspectLoci_dapc_eigenvalues_long.csv")
)

# ------------------------------------------------------------------
# Figures
# ------------------------------------------------------------------
pca_scores_all <- bind_rows(full_pca$scores, red_pca$scores)
dapc_coords_all <- bind_rows(full_dapc$coords, red_dapc$coords)

pca_plot <- ggplot(pca_scores_all, aes(PC1, PC2, color = Site)) +
  geom_point(alpha = 0.8, size = 1.8) +
  facet_wrap(~ Dataset, scales = "free") +
  theme_bw(base_size = 11) +
  labs(
    title = "PCA comparison: FULL vs noSuspectLoci",
    x = "PC1",
    y = "PC2"
  )

ggsave(
  file.path(SENS_FIGURES_DIR, paste0("pca_", analysis_suffix, ".pdf")),
  pca_plot,
  width = 8.2,
  height = 4.6
)

dapc_plot <- ggplot(dapc_coords_all, aes(LD1, LD2, color = Site)) +
  geom_point(alpha = 0.8, size = 1.8) +
  facet_wrap(~ Dataset, scales = "free") +
  theme_bw(base_size = 11) +
  labs(
    title = "DAPC comparison: FULL vs noSuspectLoci",
    x = "LD1",
    y = "LD2"
  )

ggsave(
  file.path(SENS_FIGURES_DIR, paste0("dapc_", analysis_suffix, ".pdf")),
  dapc_plot,
  width = 8.2,
  height = 4.6
)

delta_tbl <- bind_rows(
  div_overall_cmp %>% transmute(metric = "Ho", delta = delta_Ho_REDUCED_minus_FULL),
  div_overall_cmp %>% transmute(metric = "He", delta = delta_He_REDUCED_minus_FULL),
  div_overall_cmp %>% transmute(metric = "FIS", delta = delta_FIS_REDUCED_minus_FULL),
  diff_cmp %>% transmute(metric = "mean_pairwise_FST", delta = delta_mean_pairwise_FST_REDUCED_minus_FULL),
  diff_cmp %>% transmute(metric = "mean_pairwise_JostD", delta = delta_mean_pairwise_JostD_REDUCED_minus_FULL)
)

delta_plot <- ggplot(delta_tbl, aes(metric, delta, fill = metric)) +
  geom_col(show.legend = FALSE, width = 0.72) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  theme_bw(base_size = 11) +
  labs(
    title = "Sensitivity deltas (REDUCED - FULL)",
    x = NULL,
    y = "Delta"
  )

ggsave(
  file.path(SENS_FIGURES_DIR, "full_vs_noSuspectLoci_metric_deltas.pdf"),
  delta_plot,
  width = 8,
  height = 4.5
)

# ------------------------------------------------------------------
# Final console verification
# ------------------------------------------------------------------
cat("\n[hwe_sensitivity] Corrected FULL vs REDUCED by-site diversity comparison (should be 12 rows: S1-S6, N1-N6):\n")
print(div_site_cmp)
cat("[hwe_sensitivity] Corrected by-site comparison row count: ", nrow(div_site_cmp), "\n", sep = "")
cat("[hwe_sensitivity] Corrected by-site comparison site order: ", paste(div_site_cmp$Site, collapse = ", "), "\n", sep = "")

message("[hwe_sensitivity] Completed reduced-loci sensitivity branch successfully.")