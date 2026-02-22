############################################################
# scripts/01_clone_correction.R
# MLL-based clone correction using Bruvo distance (site-level)
############################################################

options(repos = c(CRAN = "https://cloud.r-project.org"))
pkgs <- c("adegenet", "poppr", "dplyr", "tibble", "here")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)

suppressPackageStartupMessages({
  library(adegenet)
  library(poppr)
  library(dplyr)
  library(tibble)
  library(here)
})

# Compatibility guard: avoid unexported package internals (no ::: and no non-exported :: symbols)
if (packageVersion("adegenet") < "2.1.0") {
  warning("adegenet version is old (", as.character(packageVersion("adegenet")), "). Script uses exported APIs only, but please consider updating adegenet.")
}

set.seed(123)

find_project_root <- function() {
  candidate_root <- here::here()
  fallback_root <- file.path(path.expand("~"), "Desktop", "BeechCode")
  if (dir.exists(file.path(candidate_root, "scripts"))) return(candidate_root)
  if (dir.exists(file.path(fallback_root, "scripts"))) return(fallback_root)
  stop("Cannot locate project root containing scripts/.")
}

PROJECT_ROOT <- find_project_root()
setwd(PROJECT_ROOT)

RUN_TAG <- "v1"
RUN_OUT <- file.path(PROJECT_ROOT, "outputs", RUN_TAG)
OBJ_DIR <- file.path(RUN_OUT, "objects")
OUTDIR <- file.path(RUN_OUT, "clone_correction")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

need <- c("gi.rds", "df_ids.rds", "meta.rds")
missing <- need[!file.exists(file.path(OBJ_DIR, need))]
if (length(missing) > 0) {
  stop(
    "Missing required object files in ", OBJ_DIR, ": ", paste(missing, collapse = ", "),
    "\nRun 00_master_pipeline.R first to create these objects."
  )
}

gi <- readRDS(file.path(OBJ_DIR, "gi.rds"))
df_ids <- readRDS(file.path(OBJ_DIR, "df_ids.rds"))
meta <- readRDS(file.path(OBJ_DIR, "meta.rds"))

if (!exists("gi") || !inherits(gi, "genind")) {
  stop("Object 'gi' is missing or not a genind. Upstream producer: 00_master_pipeline.R")
}
if (adegenet::nLoc(gi) == 0) {
  stop("genind object 'gi' has zero loci. Upstream data/preprocessing must be fixed before clone correction.")
}
if (!is.data.frame(df_ids)) {
  stop("df_ids.rds does not contain a data.frame. Upstream producer: 00_master_pipeline.R")
}
if (!is.data.frame(meta)) {
  stop("meta.rds does not contain a data.frame. Upstream producer: 00_master_pipeline.R")
}

id_candidates <- c("ind", "SampleID", "Nom_Labo_Échantillons")
id_col <- id_candidates[id_candidates %in% names(df_ids)][1]
if (is.na(id_col)) {
  stop("No individual ID column found in df_ids. Expected one of: ", paste(id_candidates, collapse = ", "))
}
if (!"Site" %in% names(df_ids)) stop("df_ids must contain a Site column.")

ids_gi <- indNames(gi)
ids_df <- as.character(df_ids[[id_col]])

if (anyDuplicated(ids_df) > 0) {
  stop("Duplicate IDs detected in df_ids column ", id_col, ".")
}
if (!setequal(ids_gi, ids_df)) {
  only_gi <- setdiff(ids_gi, ids_df)
  only_df <- setdiff(ids_df, ids_gi)
  stop(
    "ID mismatch between gi and df_ids. ",
    "Only in gi: ", if (length(only_gi)) paste(only_gi, collapse = ", ") else "none", "; ",
    "Only in df_ids: ", if (length(only_df)) paste(only_df, collapse = ", ") else "none"
  )
}

if (!identical(ids_gi, ids_df)) {
  idx <- match(ids_gi, ids_df)
  df_ids <- df_ids[idx, , drop = FALSE]
}

if (!identical(indNames(gi), as.character(df_ids[[id_col]]))) {
  stop("Failed to align df_ids to gi individual order.")
}

pop(gi) <- as.factor(df_ids$Site)
if (any(is.na(pop(gi)))) stop("Some individuals have missing Site after alignment.")

# Missing data diagnostics
# Replaced non-exported adegenet::info_table() with exported adegenet::tab(..., NA.method='asis')
# and compute per-individual missingness directly from the genotype table.
geno_tab_asis <- adegenet::tab(gi, freq = FALSE, NA.method = "asis")
if (!is.matrix(geno_tab_asis)) {
  stop("Failed to build genotype table with adegenet::tab().")
}
miss_pct <- rowMeans(is.na(geno_tab_asis)) * 100
high_missing_n <- sum(miss_pct > 35, na.rm = TRUE)
if (high_missing_n > 0) {
  warning(high_missing_n, " individuals have >35% missing data.")
}

# Missingness summaries (allele-level, locus-level, individual-level)
QC_OUTDIR <- file.path(RUN_OUT, "qc")
dir.create(QC_OUTDIR, showWarnings = FALSE, recursive = TRUE)

allele_missing <- colMeans(is.na(geno_tab_asis))
loci_fac <- as.character(adegenet::locFac(gi))
if (length(allele_missing) != length(loci_fac)) {
  stop("Mismatch between allele columns and locus factor length when computing missingness summaries.")
}
locus_missing <- tapply(allele_missing, loci_fac, mean)
locus_missing_df <- data.frame(
  locus = names(locus_missing),
  missing_prop = as.numeric(locus_missing),
  stringsAsFactors = FALSE
)

individual_missing_df <- data.frame(
  individual = rownames(geno_tab_asis),
  missing_prop = as.numeric(rowMeans(is.na(geno_tab_asis))),
  stringsAsFactors = FALSE
)

write.csv(locus_missing_df, file.path(QC_OUTDIR, "locus_missingness.csv"), row.names = FALSE)
write.csv(individual_missing_df, file.path(QC_OUTDIR, "individual_missingness.csv"), row.names = FALSE)

# MLG
mlg_vec <- tryCatch(poppr::mlg.vector(gi), error = function(e) poppr::mlg(gi))
if (length(mlg_vec) != nInd(gi)) stop("Could not compute per-individual MLG vector.")

# Bruvo distance with explicit repeat lengths (future-proof poppr behavior)
locs <- adegenet::locNames(gi)
repeat_lengths <- rep(2L, length(locs))
names(repeat_lengths) <- locs
stopifnot(length(repeat_lengths) == length(adegenet::locNames(gi)))
message("Using replen length = ", length(repeat_lengths), " for ", length(locs), " loci")
bruvo <- poppr::bruvo.dist(gi, replen = repeat_lengths)
bruvo_mat <- as.matrix(bruvo)
if (any(!is.finite(bruvo_mat))) stop("Bruvo distance matrix has non-finite values.")
bruvo_mat <- (bruvo_mat + t(bruvo_mat)) / 2
diag(bruvo_mat) <- 0

if (!identical(rownames(bruvo_mat), indNames(gi))) {
  if (!setequal(rownames(bruvo_mat), indNames(gi))) {
    stop("Bruvo matrix names do not match genind individual names.")
  }
  bruvo_mat <- bruvo_mat[indNames(gi), indNames(gi), drop = FALSE]
}

# Automatic threshold heuristic + fallback
fallback_threshold <- 0.09
heuristic_threshold <- NA_real_
threshold_method <- "fallback"

if (nInd(gi) >= 3) {
  nn <- apply(bruvo_mat + diag(Inf, nInd(gi)), 1, min, na.rm = TRUE)
  nn <- nn[is.finite(nn)]
  if (length(nn) >= 5) {
    q1 <- stats::quantile(nn, probs = 0.25, na.rm = TRUE)
    med <- stats::median(nn, na.rm = TRUE)
    madv <- stats::mad(nn, center = med, constant = 1, na.rm = TRUE)
    candidate <- as.numeric(med + 2.5 * madv)
    lower <- as.numeric(q1)
    upper <- as.numeric(stats::quantile(nn, probs = 0.95, na.rm = TRUE))
    heuristic_threshold <- min(max(candidate, lower), upper)
    
    if (is.finite(heuristic_threshold) && heuristic_threshold > 0) {
      threshold_method <- "auto_nn_mad"
    }
  }
}

chosen_threshold <- if (threshold_method == "auto_nn_mad") heuristic_threshold else fallback_threshold
if (!is.finite(chosen_threshold) || chosen_threshold <= 0) {
  warning("Automatic MLL threshold detection failed; using fallback threshold = ", fallback_threshold)
  chosen_threshold <- fallback_threshold
  threshold_method <- "fallback"
}

hc <- hclust(as.dist(bruvo_mat), method = "average")
mll_raw <- cutree(hc, h = chosen_threshold)

# stable MLL labels
map_mll <- tibble(ind = indNames(gi), mll_raw = mll_raw) %>%
  group_by(mll_raw) %>%
  summarise(cluster_size = dplyr::n(), first_ind = min(ind), .groups = "drop") %>%
  arrange(desc(cluster_size), first_ind, mll_raw) %>%
  mutate(MLL = paste0("MLL", dplyr::row_number()))

mll_lab <- map_mll$MLL[match(mll_raw, map_mll$mll_raw)]
if (any(is.na(mll_lab))) stop("Failed to create stable MLL labels.")

# build id table
if (!"ind" %in% names(df_ids)) df_ids$ind <- as.character(df_ids[[id_col]])
df_ids <- df_ids %>%
  mutate(
    MLG = as.integer(mlg_vec),
    MLL = as.character(mll_lab)
  )

# clone correction at Site x MLL
# FIX: replaced mutate(.keep = row_number() == 1L) pattern with stable slice_head() + anti_join()
# to avoid dplyr .keep argument collision and n()/row_number() masking errors.
base_tbl <- df_ids %>%
  mutate(ind = as.character(ind), Site = as.character(Site))

kept_tbl <- base_tbl %>%
  group_by(Site, MLL) %>%
  arrange(ind, .by_group = TRUE) %>%
  slice_head(n = 1) %>%
  ungroup()

dropped_tbl <- base_tbl %>%
  anti_join(kept_tbl %>% select(ind, Site, MLL), by = c("ind", "Site", "MLL")) %>%
  select(Site, ind, MLG, MLL)

kept_ids <- kept_tbl$ind
if (length(kept_ids) < 2) {
  stop("Clone correction retained too few individuals (<2). Check threshold settings.")
}

gi_cc <- gi[match(kept_ids, indNames(gi)), , drop = FALSE]
pop(gi_cc) <- as.factor(kept_tbl$Site)

df_ids_cc <- kept_tbl

# clone summary by site
simpson_from_counts <- function(x) {
  x <- as.numeric(x)
  if (sum(x) <= 1) return(NA_real_)
  p <- x / sum(x)
  1 - sum(p^2)
}

clone_summary <- base_tbl %>%
  group_by(Site) %>%
  summarise(
    N = dplyr::n(),
    G_MLG = n_distinct(MLG),
    G_MLL = n_distinct(MLL),
    R = ifelse(N > 1, (G_MLL - 1) / (N - 1), NA_real_),
    Simpson = simpson_from_counts(table(MLL)),
    .groups = "drop"
  ) %>%
  arrange(Site)

threshold_info <- data.frame(
  threshold_method = threshold_method,
  threshold_auto = heuristic_threshold,
  threshold_fallback = fallback_threshold,
  threshold_used = chosen_threshold,
  n_ind_raw = nInd(gi),
  n_ind_clone_corrected = nInd(gi_cc),
  stringsAsFactors = FALSE
)

# save objects (backward compatible)
saveRDS(gi_cc, file.path(OBJ_DIR, "gi_cc.rds"))
saveRDS(gi_cc, file.path(OBJ_DIR, "gi_mll.rds"))
saveRDS(df_ids, file.path(OBJ_DIR, "df_ids.rds"))
saveRDS(df_ids_cc, file.path(OBJ_DIR, "df_ids_cc.rds"))

# save tabular outputs
write.csv(clone_summary, file.path(OUTDIR, "clone_summary_by_site.csv"), row.names = FALSE)
write.csv(dropped_tbl, file.path(OUTDIR, "clone_dropped_individuals.csv"), row.names = FALSE)
write.csv(threshold_info, file.path(OUTDIR, "mll_threshold_info.csv"), row.names = FALSE)

cat(
  "DONE clone correction.\n",
  "Threshold method: ", threshold_method, " | used: ", signif(chosen_threshold, 4), "\n",
  "Individuals raw: ", nInd(gi), " | clone-corrected: ", nInd(gi_cc), "\n",
  "Outputs in: ", OUTDIR, "\n",
  sep = ""
)
