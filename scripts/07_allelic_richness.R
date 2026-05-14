# scripts/07_allelic_richness.R
############################################################
# Genetic diversity by site (clone-corrected genind: gi_mll)
#
# This script intentionally calculates site/population summaries from the
# clone-corrected object `gi_mll` (one representative per MLL). It rebuilds
# the population labels from df_ids_mll and the current site_lookup sheet so
# the final tables use display labels S1-S6 and N1-N6, not locus names X1-X15.
#
# Main outputs:
# - outputs/tables/heterozygosity_fis_by_site.csv
# - outputs/tables/heterozygosity_by_site.csv        (legacy mirror)
# - outputs/tables/allelic_richness_by_site.csv
# - outputs/tables/site_genetic_summary.csv
# - outputs/tables/heterozygosity_fis_overall.csv
# - outputs/tables/heterozygosity_fis_by_locus.csv   (explicit locus-level QC)
# Optional Excel mirrors (.xlsx) are written when writexl is available.
############################################################

suppressPackageStartupMessages({
  library(adegenet)
  library(hierfstat)
  library(dplyr)
  library(readxl)
})

source("scripts/_load_objects.R")

SCRIPT_TAG <- "[07_allelic_richness]"
EXPECTED_SITE_LABELS <- c("S1", "S2", "S3", "S4", "S5", "S6", "N1", "N2", "N3", "N4", "N5", "N6")
SITE_LOOKUP_SHEET <- "site_lookup"

message(SCRIPT_TAG, " Calculating by-site Ho/He/FIS and allelic richness from clone-corrected gi_mll...")

# ----------------------------
# General helpers
# ----------------------------
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
      stop(SCRIPT_TAG, " Could not find required ", label, ". Accepted names: ", paste(choices, collapse = ", "), call. = FALSE)
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

round4 <- function(x) round(as.numeric(x), 4)

round3 <- function(x) round(as.numeric(x), 3)

contains_x_site_labels <- function(x) {
  any(grepl("^X[0-9]+$", as.character(x)))
}

write_table_with_optional_excel <- function(df, csv_path) {
  dir.create(dirname(csv_path), recursive = TRUE, showWarnings = FALSE)
  write.csv(df, csv_path, row.names = FALSE, na = "")
  message(SCRIPT_TAG, " Saved: ", csv_path)
  
  if (requireNamespace("writexl", quietly = TRUE)) {
    xlsx_path <- sub("\\.csv$", ".xlsx", csv_path)
    writexl::write_xlsx(df, path = xlsx_path)
    message(SCRIPT_TAG, " Saved: ", xlsx_path)
  }
}

# ----------------------------
# site_lookup loading and label mapping
# ----------------------------
find_site_lookup_workbook <- function(raw_dir = file.path(PROJECT_ROOT, "data", "raw"), sheet = SITE_LOOKUP_SHEET) {
  if (!dir.exists(raw_dir)) {
    message(SCRIPT_TAG, " site_lookup unavailable because data/raw does not exist: ", raw_dir)
    return(NULL)
  }
  
  excel_files <- list.files(raw_dir, pattern = "\\.(xlsx|xls)$", full.names = TRUE, ignore.case = TRUE)
  if (length(excel_files) == 0) {
    message(SCRIPT_TAG, " site_lookup unavailable: no Excel workbook was found in data/raw.")
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
    message(SCRIPT_TAG, " site_lookup unavailable: no workbook in data/raw contains a '", sheet, "' sheet.")
    return(NULL)
  }
  if (length(lookup_files) > 1) {
    message(SCRIPT_TAG, " Multiple workbooks contain site_lookup; using ", basename(lookup_files[1]), ".")
  } else {
    message(SCRIPT_TAG, " Reading site_lookup from: ", basename(lookup_files[1]))
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
  
  old_site_col <- pick_column(lookup, c("Site", "site", "site_code", "old_site", "code_site"), "old site code", required = TRUE)
  label_col <- pick_column(lookup, c("Site_label", "site_label", "new_site", "site_new", "label", "site_id"), "display site label")
  order_col <- pick_column(lookup, c("Site_order", "site_order", "order", "ordre", "sort", "south_north_order"), "south-to-north site order")
  
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
    message(SCRIPT_TAG, " site_lookup was found but no usable rows were detected; raw site labels will be retained.")
    return(NULL)
  }
  if (anyDuplicated(out$old_site)) {
    stop(SCRIPT_TAG, " site_lookup contains duplicated old site codes after cleaning.", call. = FALSE)
  }
  if (anyDuplicated(out$site_label)) {
    stop(SCRIPT_TAG, " site_lookup contains duplicated display site labels after cleaning.", call. = FALSE)
  }
  
  message(SCRIPT_TAG, " Loaded site_lookup rows: ", nrow(out))
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

build_site_order <- function(site_labels) {
  sites <- unique(as.character(site_labels))
  out <- data.frame(Site = sites, stringsAsFactors = FALSE)
  
  if (!is.null(site_lookup)) {
    lookup_order <- site_lookup %>%
      transmute(Site = site_label, site_order = site_order) %>%
      distinct(Site, .keep_all = TRUE)
    out <- out %>% left_join(lookup_order, by = "Site")
  } else {
    out$site_order <- NA_real_
  }
  
  out %>%
    mutate(expected_order = match(Site, EXPECTED_SITE_LABELS)) %>%
    arrange(is.na(site_order), site_order, is.na(expected_order), expected_order, Site) %>%
    transmute(Site)
}

# ----------------------------
# Validate and rebuild the clone-corrected site assignments
# ----------------------------
validate_columns(df_ids_mll, c("ind_id", "Site"), df_name = "07_allelic_richness df_ids_mll")
if (!all(adegenet::indNames(gi_mll) == df_ids_mll$ind_id)) {
  stop(SCRIPT_TAG, " gi_mll and df_ids_mll are not aligned.", call. = FALSE)
}

raw_sites <- normalize_lookup_key(df_ids_mll$Site)
if (any(!nzchar(raw_sites))) {
  bad_ids <- df_ids_mll$ind_id[!nzchar(raw_sites)]
  stop(SCRIPT_TAG, " df_ids_mll contains missing/blank Site labels for: ", paste(head(bad_ids, 10), collapse = ", "), call. = FALSE)
}
if (contains_x_site_labels(raw_sites)) {
  stop(SCRIPT_TAG, " Raw df_ids_mll site labels contain locus-like X labels; refusing to continue.", call. = FALSE)
}

site_labels <- map_site_labels(raw_sites)
if (contains_x_site_labels(site_labels)) {
  stop(SCRIPT_TAG, " Final mapped site labels contain locus-like X labels; refusing to continue.", call. = FALSE)
}
if (!setequal(unique(site_labels), EXPECTED_SITE_LABELS)) {
  stop(
    SCRIPT_TAG, " Final site labels must be exactly S1-S6 and N1-N6.",
    "\nExpected: ", paste(EXPECTED_SITE_LABELS, collapse = ", "),
    "\nObserved: ", paste(sort(unique(site_labels)), collapse = ", "),
    call. = FALSE
  )
}

site_order <- build_site_order(site_labels)
if (!setequal(site_order$Site, EXPECTED_SITE_LABELS) || nrow(site_order) != length(EXPECTED_SITE_LABELS)) {
  stop(SCRIPT_TAG, " Site-order table must contain exactly the 12 expected sites.", call. = FALSE)
}

# Work on a local clone so sourced objects remain unchanged in interactive use.
gen_obj <- gi_mll
adegenet::pop(gen_obj) <- factor(site_labels, levels = site_order$Site)

message(SCRIPT_TAG, " Dataset used: gi_mll (clone-corrected; one representative per MLL)")
message(SCRIPT_TAG, " Individuals retained: ", adegenet::nInd(gen_obj))
message(SCRIPT_TAG, " Loci retained: ", adegenet::nLoc(gen_obj))
message(SCRIPT_TAG, " Site order: ", paste(site_order$Site, collapse = ", "))

# ----------------------------
# Ho / He / FIS by site
# ----------------------------
hf <- hierfstat::genind2hierfstat(gen_obj)
bs <- hierfstat::basic.stats(hf)
if (is.null(bs$Ho) || is.null(bs$Hs) || is.null(bs$Fis)) {
  stop(SCRIPT_TAG, " hierfstat::basic.stats did not return Ho, Hs, and Fis matrices.", call. = FALSE)
}

# hierfstat::basic.stats commonly stores loci as rows and populations as columns.
# The old script treated row names as sites, which created X1-X15 rows. This
# helper detects the orientation and always aggregates over the locus dimension
# for each site/population label.
extract_basic_stat_by_site <- function(stat_matrix, site_labels, stat_name) {
  mat <- as.matrix(stat_matrix)
  storage.mode(mat) <- "numeric"
  
  if (all(site_labels %in% colnames(mat))) {
    values <- vapply(site_labels, function(site) safe_mean(mat[, site]), numeric(1))
  } else if (all(site_labels %in% rownames(mat))) {
    values <- vapply(site_labels, function(site) safe_mean(mat[site, ]), numeric(1))
  } else {
    stop(
      SCRIPT_TAG, " Could not locate site labels in basic.stats matrix for ", stat_name, ".",
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

ho_by_site <- extract_basic_stat_by_site(bs$Ho, site_order$Site, "Ho") %>% rename(Ho = value)
he_by_site <- extract_basic_stat_by_site(bs$Hs, site_order$Site, "He") %>% rename(He = value)
fis_by_site <- extract_basic_stat_by_site(bs$Fis, site_order$Site, "FIS") %>% rename(FIS = value)

pop_sizes <- table(as.character(adegenet::pop(gen_obj)))
pop_sizes <- pop_sizes[site_order$Site]
if (any(is.na(pop_sizes)) || any(pop_sizes <= 0)) {
  stop(SCRIPT_TAG, " N per site contains missing or zero values after site assignment.", call. = FALSE)
}
if (any(pop_sizes < 2)) {
  stop(SCRIPT_TAG, " At least one site has fewer than two clone-corrected individuals: ", paste(names(pop_sizes)[pop_sizes < 2], collapse = ", "), call. = FALSE)
}

n_by_site <- data.frame(Site = names(pop_sizes), N = as.integer(pop_sizes), stringsAsFactors = FALSE)

# Explicitly save locus-level values under an unambiguous filename so they are
# not mistaken for by-site summaries.
heterozygosity_fis_by_locus <- extract_basic_stat_by_locus(bs$Ho, "Ho", site_order$Site) %>%
  full_join(extract_basic_stat_by_locus(bs$Hs, "He", site_order$Site), by = "Locus") %>%
  full_join(extract_basic_stat_by_locus(bs$Fis, "FIS", site_order$Site), by = "Locus") %>%
  mutate(Ho = round4(Ho), He = round4(He), FIS = round4(FIS))
write_table_with_optional_excel(heterozygosity_fis_by_locus, file.path(TABLES_DIR, "heterozygosity_fis_by_locus.csv"))

# ----------------------------
# Na and rarefied allelic richness by site
# ----------------------------
allele_tab <- adegenet::tab(gen_obj, NA.method = "asis")
loc_fac <- adegenet::locFac(gen_obj)
loci <- adegenet::locNames(gen_obj)
site_row_indices <- split(seq_len(nrow(allele_tab)), as.character(adegenet::pop(gen_obj)))
site_row_indices <- site_row_indices[site_order$Site]

count_alleles_for_site_locus <- function(rows, locus_name) {
  cols <- which(loc_fac == locus_name)
  if (length(cols) == 0 || length(rows) == 0) return(NA_integer_)
  mat <- allele_tab[rows, cols, drop = FALSE]
  if (all(is.na(mat))) return(NA_integer_)
  as.integer(sum(colSums(mat, na.rm = TRUE) > 0))
}

allele_counts_by_site_locus <- expand.grid(
  Site = site_order$Site,
  Locus = loci,
  stringsAsFactors = FALSE
) %>%
  rowwise() %>%
  mutate(Allele_count = count_alleles_for_site_locus(site_row_indices[[Site]], Locus)) %>%
  ungroup()

mean_alleles_by_site <- allele_counts_by_site_locus %>%
  group_by(Site) %>%
  summarise(Na = safe_mean(Allele_count), .groups = "drop")

min_n <- min(as.integer(pop_sizes), na.rm = TRUE)
message(SCRIPT_TAG, " Rarefied allelic richness sample size (min.n): ", min_n)

ar <- hierfstat::allelic.richness(hf, min.n = min_n)
if (is.null(ar$Ar)) {
  stop(SCRIPT_TAG, " hierfstat::allelic.richness did not return an Ar matrix.", call. = FALSE)
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
      SCRIPT_TAG, " Could not locate site labels in allelic.richness Ar matrix.",
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

allelic_richness_by_site <- extract_ar_by_site(ar$Ar, site_order$Site) %>%
  mutate(
    Allelic_Richness = round3(Allelic_Richness),
    Allelic_Richness_SE = round3(Allelic_Richness_SE)
  )
write_table_with_optional_excel(allelic_richness_by_site, file.path(TABLES_DIR, "allelic_richness_by_site.csv"))

# ----------------------------
# Final by-site Ho/He/FIS table
# ----------------------------
heterozygosity_fis_by_site <- n_by_site %>%
  left_join(mean_alleles_by_site, by = "Site") %>%
  left_join(allelic_richness_by_site, by = "Site") %>%
  left_join(ho_by_site, by = "Site") %>%
  left_join(he_by_site, by = "Site") %>%
  left_join(fis_by_site, by = "Site") %>%
  mutate(
    Na = round3(Na),
    Ho = round4(Ho),
    He = round4(He),
    FIS = round4(FIS),
    FIS_interpretation = case_when(
      is.na(FIS) ~ NA_character_,
      FIS > 0 ~ "heterozygote_deficit",
      FIS < 0 ~ "heterozygote_excess",
      TRUE ~ "no_deficit_or_excess"
    )
  ) %>%
  select(Site, N, Na, Allelic_Richness, Allelic_Richness_SE, Ho, He, FIS, FIS_interpretation)

if (nrow(heterozygosity_fis_by_site) != length(EXPECTED_SITE_LABELS)) {
  stop(SCRIPT_TAG, " heterozygosity_fis_by_site must have 12 rows, but has ", nrow(heterozygosity_fis_by_site), ".", call. = FALSE)
}
if (!identical(heterozygosity_fis_by_site$Site, site_order$Site)) {
  stop(SCRIPT_TAG, " heterozygosity_fis_by_site site order does not match the validated site order.", call. = FALSE)
}
if (contains_x_site_labels(heterozygosity_fis_by_site$Site)) {
  stop(SCRIPT_TAG, " heterozygosity_fis_by_site contains locus-like X labels; refusing to write output.", call. = FALSE)
}
if (any(is.na(heterozygosity_fis_by_site[, c("N", "Ho", "He", "FIS")]))) {
  stop(SCRIPT_TAG, " heterozygosity_fis_by_site contains missing N/Ho/He/FIS values.", call. = FALSE)
}

write_table_with_optional_excel(heterozygosity_fis_by_site, file.path(TABLES_DIR, "heterozygosity_fis_by_site.csv"))
write_table_with_optional_excel(heterozygosity_fis_by_site, file.path(TABLES_DIR, "heterozygosity_by_site.csv"))

# ----------------------------
# Overall Ho / He / FIS (clone-corrected, all sites pooled)
# ----------------------------
extract_overall_stat <- function(overall_obj, stat_name) {
  if (is.null(overall_obj)) return(NA_real_)
  
  if (is.atomic(overall_obj) && !is.null(names(overall_obj))) {
    val <- overall_obj[[stat_name]]
    return(if (length(val) == 0) NA_real_ else as.numeric(val[1]))
  }
  
  if (is.matrix(overall_obj) || is.data.frame(overall_obj)) {
    rn <- rownames(overall_obj)
    cn <- colnames(overall_obj)
    if (!is.null(rn) && stat_name %in% rn) return(safe_mean(overall_obj[stat_name, , drop = TRUE]))
    if (!is.null(cn) && stat_name %in% cn) return(safe_mean(overall_obj[, stat_name, drop = TRUE]))
  }
  
  NA_real_
}

overall_Ho <- extract_overall_stat(bs$overall, "Ho")
overall_He <- extract_overall_stat(bs$overall, "Hs")
overall_FIS <- extract_overall_stat(bs$overall, "Fis")
if (is.na(overall_Ho)) overall_Ho <- safe_mean(heterozygosity_fis_by_site$Ho)
if (is.na(overall_He)) overall_He <- safe_mean(heterozygosity_fis_by_site$He)
if (is.na(overall_FIS)) overall_FIS <- safe_mean(heterozygosity_fis_by_site$FIS)

overall_heterozygosity_fis <- data.frame(
  N = adegenet::nInd(gen_obj),
  Ho = round4(overall_Ho),
  He = round4(overall_He),
  FIS = round4(overall_FIS),
  FIS_interpretation = dplyr::case_when(
    is.na(overall_FIS) ~ NA_character_,
    overall_FIS > 0 ~ "heterozygote_deficit",
    overall_FIS < 0 ~ "heterozygote_excess",
    TRUE ~ "no_deficit_or_excess"
  ),
  stringsAsFactors = FALSE
)
write_table_with_optional_excel(overall_heterozygosity_fis, file.path(TABLES_DIR, "heterozygosity_fis_overall.csv"))

# ----------------------------
# Site-level summary including clonality columns when available
# ----------------------------
standardize_clonality_by_site <- function(clonality_raw) {
  if (!is.data.frame(clonality_raw) || nrow(clonality_raw) == 0) return(NULL)
  if (!"Level" %in% names(clonality_raw)) return(NULL)
  
  site_rows <- clonality_raw %>% filter(.data$Level == "site")
  if (nrow(site_rows) == 0) return(NULL)
  
  site_display <- if ("Site_label" %in% names(site_rows)) {
    as.character(site_rows$Site_label)
  } else if ("Site" %in% names(site_rows)) {
    map_site_labels(site_rows$Site)
  } else {
    rep(NA_character_, nrow(site_rows))
  }
  
  site_rows %>%
    mutate(Site = as.character(site_display)) %>%
    transmute(
      Site = Site,
      Clonality_N_individuals = if ("N_individuals" %in% names(site_rows)) as.numeric(N_individuals) else if ("N" %in% names(site_rows)) as.numeric(N) else NA_real_,
      G = if ("G" %in% names(site_rows)) as.numeric(G) else if ("N_MLL" %in% names(site_rows)) as.numeric(N_MLL) else NA_real_,
      Clonal_Richness_R = if ("Clonal_Richness_R" %in% names(site_rows)) as.numeric(Clonal_Richness_R) else if ("Clonal_Richness_MLL" %in% names(site_rows)) as.numeric(Clonal_Richness_MLL) else NA_real_,
      N_MLG = if ("N_MLG" %in% names(site_rows)) as.numeric(N_MLG) else NA_real_,
      N_MLL = if ("N_MLL" %in% names(site_rows)) as.numeric(N_MLL) else if ("G" %in% names(site_rows)) as.numeric(G) else NA_real_,
      Clonal_Richness_MLG = if ("Clonal_Richness_MLG" %in% names(site_rows)) as.numeric(Clonal_Richness_MLG) else NA_real_,
      Clonal_Richness_MLL = if ("Clonal_Richness_MLL" %in% names(site_rows)) as.numeric(Clonal_Richness_MLL) else if ("Clonal_Richness_R" %in% names(site_rows)) as.numeric(Clonal_Richness_R) else NA_real_
    ) %>%
    distinct(Site, .keep_all = TRUE)
}

clonality_file <- file.path(TABLES_DIR, "clonality_summary.csv")
clonality_by_site <- NULL
if (file.exists(clonality_file)) {
  clonality_by_site <- standardize_clonality_by_site(read.csv(clonality_file, stringsAsFactors = FALSE, check.names = FALSE))
} else {
  message(SCRIPT_TAG, " Missing clonality_summary.csv; site_genetic_summary.csv will contain genetic diversity columns only.")
}

site_genetic_summary <- heterozygosity_fis_by_site %>%
  select(Site, N, Na, Allelic_Richness, Allelic_Richness_SE, Ho, He, FIS)

if (!is.null(clonality_by_site)) {
  site_genetic_summary <- site_genetic_summary %>%
    left_join(clonality_by_site, by = "Site") %>%
    select(
      Site,
      N,
      Clonality_N_individuals,
      G,
      Clonal_Richness_R,
      N_MLG,
      N_MLL,
      Clonal_Richness_MLG,
      Clonal_Richness_MLL,
      Na,
      Allelic_Richness,
      Allelic_Richness_SE,
      Ho,
      He,
      FIS
    )
}

if (any(is.na(site_genetic_summary[, c("Ho", "He", "FIS")]))) {
  stop(SCRIPT_TAG, " site_genetic_summary contains missing Ho/He/FIS values after by-site join.", call. = FALSE)
}
if (contains_x_site_labels(site_genetic_summary$Site)) {
  stop(SCRIPT_TAG, " site_genetic_summary contains locus-like X labels; refusing to write output.", call. = FALSE)
}

write_table_with_optional_excel(site_genetic_summary, file.path(TABLES_DIR, "site_genetic_summary.csv"))

message("\n", SCRIPT_TAG, " Final by-site heterozygosity/FIS table:")
print(heterozygosity_fis_by_site, row.names = FALSE)

message("\n", SCRIPT_TAG, " Output files to check:")
message(" - ", file.path(TABLES_DIR, "heterozygosity_fis_by_site.csv"))
message(" - ", file.path(TABLES_DIR, "heterozygosity_by_site.csv"))
message(" - ", file.path(TABLES_DIR, "allelic_richness_by_site.csv"))
message(" - ", file.path(TABLES_DIR, "site_genetic_summary.csv"))
message(" - ", file.path(TABLES_DIR, "heterozygosity_fis_by_locus.csv"), " (locus-level QC, not by-site)")
message(SCRIPT_TAG, " Done.")