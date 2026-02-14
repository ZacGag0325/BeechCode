############################################################
# 00_master_pipeline.R  (BeechCode)
############################################################

### 0) PROJECT ROOT + PATHS (ROBUST) ----
suppressPackageStartupMessages(library(here))

FALLBACK_ROOT <- file.path(path.expand("~"), "Desktop", "BeechCode")
candidate_root <- here::here()

if (dir.exists(file.path(candidate_root, "data", "raw"))) {
  PROJECT_ROOT <- candidate_root
} else if (dir.exists(file.path(FALLBACK_ROOT, "data", "raw"))) {
  PROJECT_ROOT <- FALLBACK_ROOT
  setwd(PROJECT_ROOT)
} else {
  stop(
    "Can't find project root.\n",
    "Expected to find a folder with data/raw.\n\n",
    "Tried:\n  - ", candidate_root, "\n  - ", FALLBACK_ROOT, "\n\n",
    "Fix: open your BeechCode.Rproj OR change FALLBACK_ROOT."
  )
}

message("Project root detected as: ", PROJECT_ROOT)

RAW_DIR       <- file.path(PROJECT_ROOT, "data", "raw")
PROCESSED_DIR <- file.path(PROJECT_ROOT, "data", "processed")
OUTPUTS_DIR   <- file.path(PROJECT_ROOT, "outputs")

dir.create(RAW_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(PROCESSED_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(OUTPUTS_DIR, showWarnings = FALSE, recursive = TRUE)

RUN_TAG <- "v1"
RUN_OUT <- file.path(OUTPUTS_DIR, RUN_TAG)
dir.create(RUN_OUT, showWarnings = FALSE, recursive = TRUE)

message("Outputs will be saved to: ", RUN_OUT)

############################################################
# GRAPHICS SAFETY (macOS): avoid X11/XQuartz errors
############################################################
options(device = function(...) grDevices::pdf(file = file.path(RUN_OUT, "Rplots.pdf")))

x11 <- function(...) stop("This script avoids X11. Use pdf()/png() outputs instead.", call. = FALSE)
X11 <- x11

dev.new <- function(...) {
  fn <- file.path(RUN_OUT, paste0("plot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"))
  grDevices::png(filename = fn, width = 1200, height = 900, res = 150)
  message("Opened a PNG device instead of X11: ", fn)
  invisible()
}

### Pretty printing helpers ----
PRINT_N <- 12
VIEW_TABLES <- FALSE

hr <- function(txt = "") {
  cat("\n", paste(rep("=", 70), collapse = ""), "\n", sep = "")
  if (nzchar(txt)) cat(txt, "\n")
}
`%||%` <- function(a, b) if (!is.null(a)) a else b

show_tbl <- function(x, title = NULL, n = PRINT_N) {
  if (!is.null(title)) cat("\n--- ", title, " ---\n", sep = "")
  print(utils::head(x, n))
  if (isTRUE(VIEW_TABLES)) {
    if (is.data.frame(x)) utils::View(x, title = title %||% deparse(substitute(x)))
  }
}

### 1) PACKAGES ----
pkgs <- c("adegenet","poppr","readxl","dplyr","tidyr","hierfstat","pegas","ape")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)

suppressPackageStartupMessages({
  library(adegenet)
  library(poppr)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(hierfstat)
  library(pegas)
  library(ape)
})

### 2) SETTINGS ----
save_outputs <- TRUE
missing_threshold <- 0.35
chosen_thresh <- 0.09

HWE_B <- 10000
HWE_ALPHA <- 0.05

### Helper: always get per-individual MLG vector across poppr versions ----
get_mlg_vec <- function(x) {
  mv <- tryCatch(poppr::mlg.vector(x), error = function(e) NULL)
  if (!is.null(mv) && length(mv) == adegenet::nInd(x)) return(mv)
  
  m <- poppr::mlg(x)
  if (length(m) == adegenet::nInd(x)) return(m)
  
  stop("Could not extract per-individual MLG vector from poppr. poppr version: ",
       as.character(packageVersion("poppr")))
}

### Helper: deterministic relabeling of cluster IDs ----
relabel_mll_ids <- function(raw_id, sample_ids) {
  stopifnot(length(raw_id) == length(sample_ids))
  tmp <- data.frame(
    raw_id = as.integer(raw_id),
    sample_id = as.character(sample_ids),
    stringsAsFactors = FALSE
  )
  
  map <- tmp %>%
    group_by(raw_id) %>%
    summarise(
      cluster_size = n(),
      first_sample = min(sample_id),
      .groups = "drop"
    ) %>%
    arrange(desc(cluster_size), first_sample, raw_id) %>%
    mutate(stable_id = row_number())
  
  out <- map$stable_id[match(tmp$raw_id, map$raw_id)]
  if (any(is.na(out))) stop("Failed to relabel MLL cluster IDs deterministically.")
  as.integer(out)
}

### Helper: build deterministic MLL from Bruvo distances + threshold ----
make_mll_from_bruvo <- function(genind_obj, threshold) {
  bruvo <- poppr::bruvo.dist(genind_obj)
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

### Helper: make genotype df safe for pegas::as.loci + hw.test ----
safe_genind2loci_df <- function(gpop, sep = "/") {
  gdf <- adegenet::genind2df(gpop, sep = sep)
  gdf <- as.data.frame(lapply(gdf, function(col) {
    x <- trimws(as.character(col))
    x[x %in% c("", "NA", "NaN", "0", "0/0", "NA/NA")] <- NA
    ok <- is.na(x) | grepl("^[^/]+/[^/]+$", x)
    x[!ok] <- NA
    factor(x)
  }), check.names = FALSE, stringsAsFactors = FALSE)
  gdf
}

############################################################
# FIXED ALLELE COUNTS (version-proof)
############################################################
allele_count_summary <- function(genind_obj) {
  M <- adegenet::tab(genind_obj, freq = FALSE, NA.method = "zero")
  allele_counts <- colSums(M, na.rm = TRUE)
  
  lf <- adegenet::locFac(genind_obj)
  if (length(lf) != length(allele_counts)) {
    lf <- factor(rep(locNames(genind_obj), genind_obj@loc.n.all))
  }
  
  per_locus <- data.frame(
    Locus = levels(lf),
    n_alleles = as.integer(table(lf)[levels(lf)]),
    total_allele_obs = as.numeric(tapply(allele_counts, lf, sum)[levels(lf)]),
    stringsAsFactors = FALSE
  )
  
  overall <- data.frame(
    nInd = adegenet::nInd(genind_obj),
    nLoc = adegenet::nLoc(genind_obj),
    total_unique_alleles = sum(per_locus$n_alleles),
    stringsAsFactors = FALSE
  )
  
  list(overall = overall, per_locus = per_locus, raw_allele_counts = allele_counts)
}

### Helper: print clone groups ----
print_clone_groups <- function(df_ids, label_col = "MLL", title = NULL, n = 20) {
  stopifnot(label_col %in% names(df_ids))
  
  gr <- df_ids %>%
    count(.data[[label_col]], Site, name = "n_inds") %>%
    arrange(desc(n_inds))
  
  if (!is.null(title)) hr(title)
  cat("\nTop groups by ", label_col, " (largest first):\n", sep = "")
  show_tbl(gr, n = n)
  
  multi <- df_ids %>%
    group_by(.data[[label_col]]) %>%
    filter(n() > 1) %>%
    arrange(.data[[label_col]], Site, ind) %>%
    ungroup()
  
  if (nrow(multi) == 0) {
    cat("\nNo multi-individual groups for ", label_col, ".\n", sep = "")
  } else {
    cat("\nIndividuals in multi-individual groups for ", label_col, ":\n", sep = "")
    show_tbl(multi %>% select(ind, Site, MLG, MLL, Latitude, Longitude, Dhp_tige), n = n)
  }
  
  invisible(list(groups = gr, multi = multi))
}

### 3) INPUT CHECKS ----
geno_file  <- file.path(RAW_DIR, "poppr.xlsx")
field_file <- file.path(RAW_DIR, "donnees_modifiees_west_summer2024 copie.xlsx")
field_sheet <- "genetique"

if (!file.exists(geno_file)) {
  stop("Missing genotyping file at: ", geno_file,
       "\nPut poppr.xlsx into: ", file.path(PROJECT_ROOT, "data", "raw"))
}
if (!file.exists(field_file)) {
  stop("Missing field file at: ", field_file,
       "\nPut the field .xlsx into: ", file.path(PROJECT_ROOT, "data", "raw"))
}

### 4) READ GENOTYPE DATA ----
hr("STEP 4) READ GENOTYPE DATA")

geno_raw <- read_xlsx(geno_file) %>%
  dplyr::select(where(~ !all(is.na(.)))) %>%
  dplyr::select(-matches("^Unnamed")) %>%
  dplyr::select(-matches("^\\.\\.\\.[0-9]+$"))

names(geno_raw) <- trimws(names(geno_raw))

req_geno <- c("Nom_Labo_Échantillons", "Numéro_Population")
if (!all(req_geno %in% names(geno_raw))) {
  stop("poppr.xlsx must contain: ", paste(req_geno, collapse = ", "),
       "\nColumns found: ", paste(names(geno_raw), collapse = ", "))
}

geno_raw$Nom_Labo_Échantillons <- as.character(geno_raw$Nom_Labo_Échantillons)
geno <- geno_raw

cat("\nColumns in geno_raw:\n")
print(names(geno))
show_tbl(geno, "geno (first rows)", n = 8)

cat("\nPop counts (before QC):\n")
print(sort(table(geno$Numéro_Population), decreasing = TRUE))

### 5) READ FIELD DATA ----
hr("STEP 5) READ FIELD DATA")

field <- read_xlsx(field_file, sheet = field_sheet)
names(field) <- trimws(names(field))

req_field <- c("Nom_Labo_Échantillons", "Site")
if (!all(req_field %in% names(field))) {
  stop("Field sheet must contain: ", paste(req_field, collapse = ", "),
       "\nColumns found: ", paste(names(field), collapse = ", "))
}

field <- field %>%
  mutate(
    Nom_Labo_Échantillons = as.character(Nom_Labo_Échantillons),
    Site = as.character(Site)
  )

has_lat <- "Latitude"  %in% names(field)
has_lon <- "Longitude" %in% names(field)
has_dbh <- "Dhp_tige"  %in% names(field)

if (has_lat) field$Latitude  <- as.numeric(field$Latitude)
if (has_lon) field$Longitude <- as.numeric(field$Longitude)
if (has_dbh) field$Dhp_tige  <- as.numeric(field$Dhp_tige)

show_tbl(field, "field (first rows)", n = 8)

cat("\nField Site counts (raw):\n")
print(sort(table(field$Site), decreasing = TRUE))

### 6) JOIN FIELD INTO GENO ----
hr("STEP 6) JOIN FIELD INTO GENO")

field_keep <- c("Nom_Labo_Échantillons", "Site")
if (has_lat) field_keep <- c(field_keep, "Latitude")
if (has_lon) field_keep <- c(field_keep, "Longitude")
if (has_dbh) field_keep <- c(field_keep, "Dhp_tige")

field2 <- field %>%
  dplyr::select(all_of(field_keep)) %>%
  dplyr::distinct(Nom_Labo_Échantillons, .keep_all = TRUE)

geno <- geno %>% left_join(field2, by = "Nom_Labo_Échantillons")

n_missing_site <- sum(is.na(geno$Site))
cat("\nAfter join: missing Site count =", n_missing_site, "out of", nrow(geno), "\n")

if (n_missing_site > 0) {
  cat("\nExamples of missing Site IDs:\n")
  print(head(geno$Nom_Labo_Échantillons[is.na(geno$Site)], 20))
}

cat("\nJoined Site counts (including NA):\n")
print(sort(table(geno$Site, useNA = "ifany"), decreasing = TRUE))
show_tbl(geno, "geno after join (first rows)", n = 8)

### 7) QC: MISSINGNESS FILTER ----
hr("STEP 7) QC: MISSINGNESS FILTER")

meta_cols <- c("Nom_Labo_Échantillons","Numéro_Population","Site","Latitude","Longitude","Dhp_tige")
loci_cols <- setdiff(names(geno), meta_cols)
loci_cols <- loci_cols[!grepl("^\\.\\.\\.[0-9]+$", loci_cols)]

cat("\n# Loci columns detected =", length(loci_cols), "\n")
cat("Example loci columns:\n")
print(head(loci_cols, 20))

miss_prop <- apply(geno[, loci_cols], 1, function(r) mean(is.na(r) | r == ""))
geno$missing_prop <- miss_prop

cat("\nMissingness summary:\n")
print(summary(geno$missing_prop))

show_tbl(
  geno %>% dplyr::select(Nom_Labo_Échantillons, Site, missing_prop) %>% arrange(desc(missing_prop)),
  "Top missingness individuals (highest first)",
  n = 12
)

geno_filt <- geno %>% filter(missing_prop <= missing_threshold)

cat("\nRemoved ", nrow(geno) - nrow(geno_filt),
    " individuals over missing_threshold=", missing_threshold, "\n", sep = "")
cat("Remaining individuals: ", nrow(geno_filt), "\n", sep = "")

cat("\nRemaining Site counts after missingness filter:\n")
print(sort(table(geno_filt$Site, useNA = "ifany"), decreasing = TRUE))

if (save_outputs) {
  write.csv(geno_filt, file.path(RUN_OUT, "geno_after_missingness_filter.csv"), row.names = FALSE)
}

### 8) BUILD GENIND ----
hr("STEP 8) BUILD GENIND")

geno_tmp <- geno_filt[, loci_cols, drop = FALSE]
geno_tmp <- as.data.frame(lapply(geno_tmp, function(x) {
  x <- trimws(as.character(x))
  x[x == ""] <- NA
  x
}))

locus_bases <- sort(unique(sub("_[12]$", "", colnames(geno_tmp))))

missing_pairs <- locus_bases[
  !paste0(locus_bases, "_1") %in% colnames(geno_tmp) |
    !paste0(locus_bases, "_2") %in% colnames(geno_tmp)
]
if (length(missing_pairs) > 0) {
  stop("Some loci are missing _1 or _2 columns: ", paste(missing_pairs, collapse = ", "))
}

cat("\n# Locus bases =", length(locus_bases), "\n")
print(head(locus_bases, 30))

geno_mat <- as.data.frame(lapply(locus_bases, function(locus) {
  a1 <- geno_tmp[[paste0(locus, "_1")]]
  a2 <- geno_tmp[[paste0(locus, "_2")]]
  
  out <- ifelse(is.na(a1) | is.na(a2), NA_character_, paste(a1, a2, sep = "/"))
  
  out <- ifelse(
    is.na(out), NA_character_,
    vapply(strsplit(out, "/", fixed = TRUE), function(z) {
      zz <- suppressWarnings(as.numeric(z))
      if (all(!is.na(zz))) paste(sort(zz), collapse = "/") else paste(sort(z), collapse = "/")
    }, character(1))
  )
  out
}))
names(geno_mat) <- locus_bases
colnames(geno_mat) <- make.unique(gsub("\\.+", "_", colnames(geno_mat)))

cat("\nGenotype matrix dimensions:\n")
print(dim(geno_mat))
show_tbl(geno_mat, "geno_mat (first rows)", n = 6)

gi <- df2genind(
  geno_mat,
  ploidy = 2,
  sep = "/",
  ind.names = geno_filt$Nom_Labo_Échantillons,
  pop = as.factor(geno_filt$Numéro_Population),
  type = "codom"
)

strata(gi) <- data.frame(Site = as.factor(geno_filt$Site))
setPop(gi) <- ~Site

cat("\nGenind created.\n")
cat("N inds:", nInd(gi), " | N loci:", nLoc(gi), " | N pops (Site):", nPop(gi), "\n")
cat("Site sizes (nInd per Site):\n")
print(sort(table(pop(gi)), decreasing = TRUE))

# ---- ALLELE COUNTS (WHOLE DATASET) ----
hr("ALLELE COUNTS (WHOLE DATASET)")
acs <- allele_count_summary(gi)
print(acs$overall)
show_tbl(acs$per_locus %>% arrange(desc(n_alleles)),
         "Alleles per locus (sorted by n_alleles)",
         n = 25)

### 9) MLG/MLL ----
hr("STEP 9) MLG/MLL")

mlg_vec <- get_mlg_vec(gi)

mll_obj <- make_mll_from_bruvo(gi, threshold = chosen_thresh)
bruvo_mat <- mll_obj$bruvo_mat
mll_id <- mll_obj$mll_id
mll_lab <- mll_obj$mll_label

df_ids <- data.frame(
  ind  = indNames(gi),
  Nom_Labo_Échantillons = indNames(gi),
  Site = pop(gi),
  MLG  = sprintf("MLG_%03d", as.integer(factor(mlg_vec))),
  MLL  = mll_lab,
  stringsAsFactors = FALSE
) %>%
  left_join(
    geno_filt %>%
      dplyr::select(Nom_Labo_Échantillons, Latitude, Longitude, Dhp_tige) %>%
      dplyr::distinct(Nom_Labo_Échantillons, .keep_all = TRUE),
    by = c("ind" = "Nom_Labo_Échantillons")
  )

show_tbl(df_ids, "df_ids (first rows)", n = 12)

cat("\nUnique MLG:", length(unique(df_ids$MLG)),
    " | Unique MLL:", length(unique(df_ids$MLL)), "\n")


if (length(unique(df_ids$MLL)) > length(unique(df_ids$MLG))) {
  stop("Invalid MLL assignment: number of MLL exceeds number of MLG.")
}
if (length(unique(df_ids$MLL)) == length(unique(df_ids$MLG))) {
  cat("\n[WARNING] MLL count equals MLG count at threshold=", chosen_thresh,
      ". This can be biologically real but indicates no cluster collapsing.\n", sep = "")
}

print_clone_groups(df_ids, "MLL", title = "CLONE CHECK (MLL groups)")
print_clone_groups(df_ids, "MLG", title = "CLONE CHECK (MLG groups)")

if (save_outputs) {
  write.csv(df_ids, file.path(RUN_OUT, "ID_MLG_MLL_with_field.csv"), row.names = FALSE)
  write.csv(df_ids %>% dplyr::select(ind, Site, MLG, MLL),
            file.path(RUN_OUT, "MLG_MLL_outputs.csv"), row.names = FALSE)
}

### 10) CLONALITY BY SITE ----
hr("STEP 10) CLONALITY BY SITE")

clonality_by_site <- df_ids %>%
  group_by(Site) %>%
  summarise(
    n_inds = n(),
    n_MLG  = n_distinct(MLG),
    n_MLL  = n_distinct(MLL),
    genotypic_richness_MLG = n_MLG / n_inds,
    genotypic_richness_MLL = n_MLL / n_inds,
    .groups = "drop"
  ) %>%
  arrange(desc(n_inds))

print(clonality_by_site)

if (save_outputs) {
  write.csv(clonality_by_site, file.path(RUN_OUT, "clonality_by_site.csv"), row.names = FALSE)
}

############################################################
# 12) HARDY-WEINBERG (per Site, per locus)
############################################################
hr("STEP 12) HARDY-WEINBERG (pegas::hw.test)")

sites <- levels(pop(gi))
cat("\nSites in pop(gi):\n")
print(sites)

extract_hw_pvals <- function(ht) {
  if (is.list(ht) && !is.null(ht$p.value)) {
    pv <- ht$p.value
    return(list(pvals = as.numeric(pv), locus = names(pv)))
  }
  if (is.matrix(ht) || is.data.frame(ht)) {
    pmat <- as.matrix(ht)
    cn <- tolower(colnames(pmat))
    col_pick <- which(
      cn %in% c("p.value", "pvalue", "p", "pval", "pvals") |
        grepl("^pr", cn) |
        grepl("p\\s*value", cn)
    )
    if (length(col_pick) == 0) col_pick <- ncol(pmat)
    pvals <- suppressWarnings(as.numeric(pmat[, col_pick[1]]))
    locus <- rownames(pmat)
    return(list(pvals = pvals, locus = locus))
  }
  pv <- suppressWarnings(as.numeric(ht))
  return(list(pvals = pv, locus = names(ht)))
}

hwe_rows <- list()

for (s in sites) {
  gi_s <- popsub(gi, sublist = s)
  n_s <- nInd(gi_s)
  
  if (n_s < 8) {
    cat("\n[HWE] Site", s, ": nInd =", n_s, " < 8 -> skipped\n")
    next
  }
  
  cat("\n[HWE] Site", s, ": nInd =", n_s, "\n")
  
  gdf <- safe_genind2loci_df(gi_s, sep = "/")
  
  allele_counts_seen <- sapply(unlist(gdf), function(z) {
    if (is.na(z) || z == "") return(NA_integer_)
    length(strsplit(as.character(z), "/", fixed = TRUE)[[1]])
  })
  
  all_diploid_for_pegas <- all(is.na(allele_counts_seen) | allele_counts_seen == 2)
  B_use <- if (all_diploid_for_pegas) HWE_B else 0
  if (!all_diploid_for_pegas) cat("[HWE] Note: non-2 allele entries -> B=0\n")
  
  loci_obj <- pegas::as.loci(gdf)
  
  ht <- tryCatch(
    pegas::hw.test(loci_obj, B = B_use),
    error = function(e) {
      cat("[HWE] ERROR at site ", s, ": ", conditionMessage(e), "\n", sep = "")
      NULL
    }
  )
  if (is.null(ht)) next
  
  ex <- extract_hw_pvals(ht)
  pvals <- ex$pvals
  locus_names <- ex$locus
  
  if (length(pvals) == 0) next
  if (is.null(locus_names) || length(locus_names) != length(pvals)) {
    locus_names <- paste0("Locus_", seq_along(pvals))
  }
  
  tmp <- data.frame(
    Site  = rep(s, length(pvals)),
    Locus = locus_names,
    p_HWE = as.numeric(pvals),
    stringsAsFactors = FALSE
  ) %>%
    filter(!is.na(p_HWE))
  
  cat("[HWE] Site", s, ": loci tested =", nrow(tmp),
      " | n significant (p<", HWE_ALPHA, ") =",
      sum(tmp$p_HWE < HWE_ALPHA, na.rm = TRUE), "\n")
  
  hwe_rows[[s]] <- tmp
}

hwe_results <- if (length(hwe_rows) > 0) do.call(rbind, hwe_rows) else data.frame()

if (save_outputs) {
  write.csv(hwe_results, file.path(RUN_OUT, "hwe_by_site_by_locus.csv"), row.names = FALSE)
  cat("\nSaved HWE table to: ", file.path(RUN_OUT, "hwe_by_site_by_locus.csv"), "\n", sep = "")
}

############################################################
# 12B) GLOBAL HWE (MLG vs clone-corrected MLL)
############################################################
hr("STEP 12B) HWE WHOLE DATASET (MLG vs MLL)")

run_hwe_global_pegas <- function(genind_obj, label = "DATASET") {
  cat("\n[HWE-GLOBAL] Running on:", label, "\n")
  gdf <- safe_genind2loci_df(genind_obj, sep = "/")
  
  allele_counts_seen <- sapply(unlist(gdf), function(z) {
    if (is.na(z) || z == "") return(NA_integer_)
    length(strsplit(as.character(z), "/", fixed = TRUE)[[1]])
  })
  
  all_diploid_for_pegas <- all(is.na(allele_counts_seen) | allele_counts_seen == 2)
  B_use <- if (all_diploid_for_pegas) HWE_B else 0
  
  loci_obj <- pegas::as.loci(gdf)
  ht <- tryCatch(pegas::hw.test(loci_obj, B = B_use), error = function(e) NULL)
  if (is.null(ht)) return(data.frame())
  
  ex <- extract_hw_pvals(ht)
  pvals <- ex$pvals
  locus_names <- ex$locus
  if (length(pvals) == 0) return(data.frame())
  if (is.null(locus_names) || length(locus_names) != length(pvals)) {
    locus_names <- paste0("Locus_", seq_along(pvals))
  }
  
  data.frame(
    Dataset = label,
    Locus   = locus_names,
    p_HWE   = as.numeric(pvals),
    stringsAsFactors = FALSE
  ) %>% filter(!is.na(p_HWE)) %>% arrange(p_HWE)
}

subset_one_rep_per_group_per_site <- function(genind_obj, df_ids, label_col = "MLL") {
  reps <- df_ids %>%
    group_by(Site, .data[[label_col]]) %>%
    summarise(ind = dplyr::first(ind), .groups = "drop")
  keep_inds <- reps$ind
  keep_inds <- keep_inds[keep_inds %in% indNames(genind_obj)]
  genind_obj[keep_inds, , drop = FALSE]
}

hwe_global_mlg <- run_hwe_global_pegas(gi, label = "GLOBAL_MLG")
gi_mll <- subset_one_rep_per_group_per_site(gi, df_ids, label_col = "MLL")
hwe_global_mll <- run_hwe_global_pegas(gi_mll, label = "GLOBAL_MLL")

hwe_compare <- full_join(
  hwe_global_mlg %>% select(Locus, p_MLG = p_HWE),
  hwe_global_mll %>% select(Locus, p_MLL = p_HWE),
  by = "Locus"
) %>%
  mutate(
    sig_MLG = !is.na(p_MLG) & p_MLG < HWE_ALPHA,
    sig_MLL = !is.na(p_MLL) & p_MLL < HWE_ALPHA
  )

if (save_outputs) {
  write.csv(hwe_global_mlg, file.path(RUN_OUT, "hwe_global_mlg.csv"), row.names = FALSE)
  write.csv(hwe_global_mll, file.path(RUN_OUT, "hwe_global_mll.csv"), row.names = FALSE)
  write.csv(hwe_compare, file.path(RUN_OUT, "hwe_global_compare_mlg_vs_mll.csv"), row.names = FALSE)
}

############################################################
# SAVE SHARED OBJECTS (used by scripts/01_*, 02_*, 04_*, 05_*)
############################################################
hr("SAVE SHARED OBJECTS (gi + df_ids + meta + gi_mll)")

OBJ_DIR <- file.path(RUN_OUT, "objects")
dir.create(OBJ_DIR, showWarnings = FALSE, recursive = TRUE)

saveRDS(gi,     file.path(OBJ_DIR, "gi.rds"))
saveRDS(df_ids, file.path(OBJ_DIR, "df_ids.rds"))

meta_out <- field2 %>%
  dplyr::rename(ind = Nom_Labo_Échantillons) %>%
  dplyr::distinct(ind, .keep_all = TRUE)

saveRDS(meta_out, file.path(OBJ_DIR, "meta.rds"))

if (!(exists("gi_mll") && inherits(gi_mll, "genind"))) {
  gi_mll <- subset_one_rep_per_group_per_site(gi, df_ids, label_col = "MLL")
}
if (!inherits(gi_mll, "genind") || adegenet::nInd(gi_mll) < 1) {
  stop("Could not build gi_mll for object saving.")
}
saveRDS(gi_mll, file.path(OBJ_DIR, "gi_mll.rds"))

cat("Saved objects in: ", OBJ_DIR, "\n", sep = "")
cat(" - gi.rds\n - df_ids.rds\n - meta.rds\n - gi_mll.rds\n")

hr("PIPELINE COMPLETE")
cat("Key outputs written to:\n- ", RUN_OUT, "\n", sep = "")

