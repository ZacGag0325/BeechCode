# filename: scripts/04_amova.R
############################################################
# scripts/04_amova.R
# AMOVA among Sites (+ permutation test)
############################################################

options(repos = c(CRAN = "https://cloud.r-project.org"))

pkgs <- c("adegenet", "poppr", "ade4", "here", "dplyr")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)

suppressPackageStartupMessages(library(adegenet))
suppressPackageStartupMessages(library(poppr))
suppressPackageStartupMessages(library(ade4))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(dplyr))

# -----------------------------------------------------------------------------
# Load objects only if not already present in environment
# -----------------------------------------------------------------------------
if (!exists("RUN_OUT", inherits = FALSE)) {
  candidate_root <- here::here()
  fallback_root <- file.path(path.expand("~"), "Desktop", "BeechCode")
  
  if (dir.exists(file.path(candidate_root, "data", "raw"))) {
    project_root <- candidate_root
  } else if (dir.exists(file.path(fallback_root, "data", "raw"))) {
    project_root <- fallback_root
  } else {
    stop("Cannot determine project root to build RUN_OUT. Set RUN_OUT or run from BeechCode project root.")
  }
  
  RUN_TAG <- "v1"
  RUN_OUT <- file.path(project_root, "outputs", RUN_TAG)
}

OBJ_DIR <- file.path(RUN_OUT, "objects")

# Load gi_mll first (preferred for clone-corrected analyses), fallback to gi
if (!exists("gi_mll", inherits = FALSE) && file.exists(file.path(OBJ_DIR, "gi_mll.rds"))) {
  gi_mll <- readRDS(file.path(OBJ_DIR, "gi_mll.rds"))
}
if (!exists("gi", inherits = FALSE) && file.exists(file.path(OBJ_DIR, "gi.rds"))) {
  gi <- readRDS(file.path(OBJ_DIR, "gi.rds"))
}
if (!exists("df_ids", inherits = FALSE)) {
  if (file.exists(file.path(OBJ_DIR, "df_ids.rds"))) {
    df_ids <- readRDS(file.path(OBJ_DIR, "df_ids.rds"))
  } else {
    stop("Missing df_ids object and file: ", file.path(OBJ_DIR, "df_ids.rds"))
  }
}

if (exists("gi_mll", inherits = FALSE)) {
  gi_use <- gi_mll
  gi_label <- "gi_mll"
} else if (exists("gi", inherits = FALSE)) {
  gi_use <- gi
  gi_label <- "gi"
} else {
  stop("Neither gi_mll nor gi is available in memory or in ", OBJ_DIR)
}

OUTDIR <- file.path(RUN_OUT, "amova_only")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# Validate columns and robustly align individuals between gi_use and df_ids
# -----------------------------------------------------------------------------
if (!("Site" %in% names(df_ids))) stop("df_ids must contain a 'Site' column.")
if (!("Latitude" %in% names(df_ids))) stop("df_ids must contain a 'Latitude' column.")

# Choose ID column that best matches gi individual names.
id_candidates <- c("ind", "Nom_Labo_Échantillons", "SampleID")
id_candidates <- id_candidates[id_candidates %in% names(df_ids)]
if (length(id_candidates) == 0) {
  stop("No individual ID column found in df_ids. Expected one of: ind, Nom_Labo_Échantillons, SampleID")
}

gi_ids <- as.character(adegenet::indNames(gi_use))

score_match <- function(colname) {
  vals <- as.character(df_ids[[colname]])
  vals <- vals[!is.na(vals)]
  length(intersect(unique(vals), unique(gi_ids)))
}

id_scores <- sapply(id_candidates, score_match)
id_col <- id_candidates[which.max(id_scores)]

cat("[AMOVA] Using genetic object:", gi_label, "\n")
cat("[AMOVA] ID column selected from df_ids:", id_col, "\n")
cat("[AMOVA] N individuals in", gi_label, ":", length(gi_ids), "\n")
cat("[AMOVA] N rows in df_ids before alignment:", nrow(df_ids), "\n")

# Build working ID vector
if (!("ind" %in% names(df_ids))) {
  df_ids$ind <- as.character(df_ids[[id_col]])
} else {
  # prefer 'ind' if it actually matches better than selected column
  ind_score <- score_match("ind")
  if (ind_score < id_scores[which.max(id_scores)]) {
    df_ids$ind <- as.character(df_ids[[id_col]])
  } else {
    df_ids$ind <- as.character(df_ids$ind)
  }
}

# De-duplicate metadata if needed (keep first occurrence and message)
if (anyDuplicated(df_ids$ind) > 0) {
  dup_vals <- unique(df_ids$ind[duplicated(df_ids$ind)])
  cat("[AMOVA] Duplicate IDs found in df_ids. Keeping first row per ID.\n")
  cat("[AMOVA] Duplicate ID(s):", paste(head(dup_vals, 20), collapse = ", "), if (length(dup_vals) > 20) " ..." else "", "\n")
  df_ids <- df_ids %>% distinct(ind, .keep_all = TRUE)
}

# ---- ALIGN IDS block ----
# This prevents mismatch errors by:
# 1) dropping metadata IDs not present in gi,
# 2) checking if any gi IDs are missing in metadata,
# 3) reordering metadata to EXACTLY match gi order.

extra_in_df <- setdiff(df_ids$ind, gi_ids)
if (length(extra_in_df) > 0) {
  cat("[AMOVA] IDs present in df_ids but not in", gi_label, ":", length(extra_in_df), "(dropping)\n")
  cat("[AMOVA] Dropped IDs:", paste(extra_in_df, collapse = ", "), "\n")
}

# Filter out extras from metadata
df_ids_aligned <- df_ids %>%
  mutate(ind = as.character(ind)) %>%
  filter(ind %in% gi_ids)

# If gi has IDs absent in df_ids, optionally try recovery via meta
missing_in_df <- setdiff(gi_ids, df_ids_aligned$ind)
if (length(missing_in_df) > 0) {
  cat("[AMOVA] IDs present in", gi_label, "but missing in df_ids:", length(missing_in_df), "\n")
  
  if (exists("meta", inherits = FALSE) && is.data.frame(meta) && "Nom_Labo_Échantillons" %in% names(meta)) {
    cat("[AMOVA] Attempting to recover missing metadata from 'meta' using Nom_Labo_Échantillons...\n")
    meta_work <- meta
    meta_work$ind <- as.character(meta_work$Nom_Labo_Échantillons)
    
    to_add <- meta_work %>%
      filter(ind %in% missing_in_df) %>%
      select(any_of(names(df_ids_aligned)))
    
    if (nrow(to_add) > 0) {
      df_ids_aligned <- bind_rows(df_ids_aligned, to_add)
      cat("[AMOVA] Recovered", nrow(to_add), "row(s) from meta.\n")
    }
  }
  
  missing_after_recovery <- setdiff(gi_ids, df_ids_aligned$ind)
  if (length(missing_after_recovery) > 0) {
    stop(
      "Cannot align IDs: some individuals in ", gi_label, " are missing in df_ids metadata. Missing IDs: ",
      paste(missing_after_recovery, collapse = ", ")
    )
  }
}

# Enforce exact order of metadata to match gi
df_ids_aligned <- df_ids_aligned %>%
  mutate(ind = as.character(ind)) %>%
  slice(match(gi_ids, ind))

if (any(is.na(df_ids_aligned$ind))) {
  stop("Alignment failed: NA IDs produced after ordering metadata to gi IDs.")
}

stopifnot(identical(as.character(df_ids_aligned$ind), as.character(gi_ids)))
cat("[AMOVA] N rows in df_ids after alignment:", nrow(df_ids_aligned), "\n")
cat("[AMOVA] ID alignment successful: IDs are identical and ordered.\n")

# replace df_ids with aligned version for all downstream analysis
df_ids <- df_ids_aligned

# Ensure pop is Site (consistent)
pop(gi_use) <- as.factor(df_ids$Site)

# -----------------------------------------------------------------------------
# Existing AMOVA by Site
# -----------------------------------------------------------------------------
am <- poppr::poppr.amova(gi_use, ~Site, within = TRUE)
print(am)

set.seed(123)
am_test <- ade4::randtest(am, nrepet = 999)
print(am_test)

capture.output(am,      file = file.path(OUTDIR, "amova_site.txt"))
capture.output(am_test, file = file.path(OUTDIR, "amova_site_randtest.txt"))

pdf(file.path(OUTDIR, "amova_site_randtest.pdf"), width = 7, height = 5)
plot(am_test)
dev.off()

############################################################
# STEP: AMOVA North vs South (Region_NS/Site)
############################################################

# Missing latitude check by site
lat_na_sites <- sort(unique(as.character(df_ids$Site[is.na(df_ids$Latitude)])))
if (length(lat_na_sites) > 0) {
  stop("Missing latitude values detected for Site(s): ", paste(lat_na_sites, collapse = ", "))
}

# Site mean latitude (works whether latitude is individual-level or constant per site)
site_lat <- stats::aggregate(Latitude ~ Site, data = df_ids, FUN = mean)
site_lat <- site_lat[order(site_lat$Latitude), , drop = FALSE]

if (nrow(site_lat) != 12) {
  stop("Expected 12 sites for North/South split, found ", nrow(site_lat), ".")
}

# Assign 6 southernmost + 6 northernmost sites
site_lat$Region_NS <- c(rep("South", 6), rep("North", 6))
site_to_region <- site_lat[, c("Site", "Region_NS")]

# Attach Region_NS to individuals
map_idx <- match(df_ids$Site, site_to_region$Site)
if (any(is.na(map_idx))) {
  missing_sites <- sort(unique(as.character(df_ids$Site[is.na(map_idx)])))
  stop("Could not map some Site values to Region_NS: ", paste(missing_sites, collapse = ", "))
}

df_ids$Region_NS <- factor(site_to_region$Region_NS[map_idx], levels = c("South", "North"))

# Build hierarchy in gi
strata(gi_use) <- data.frame(
  Region_NS = df_ids$Region_NS,
  Site = df_ids$Site,
  row.names = gi_ids,
  stringsAsFactors = FALSE
)
setPop(gi_use) <- ~Site

# Hierarchical AMOVA and permutation test
am_ns <- poppr::poppr.amova(gi_use, ~Region_NS/Site, within = TRUE)
print(am_ns)

cat("\nVariance components and percentages (Region_NS/Site):\n")
print(am_ns$componentsofcovariance)

if ("Sigma" %in% colnames(am_ns$componentsofcovariance)) {
  sigma_vals <- am_ns$componentsofcovariance[, "Sigma"]
} else {
  sigma_vals <- am_ns$componentsofcovariance[, ncol(am_ns$componentsofcovariance)]
}

var_pct <- 100 * sigma_vals / sum(sigma_vals)
print(data.frame(
  Component = rownames(am_ns$componentsofcovariance),
  Percent = var_pct,
  row.names = NULL
))

cat("\nPhi-statistics (Region_NS/Site):\n")
print(am_ns$statphi)

set.seed(123)
am_ns_test <- ade4::randtest(am_ns, nrepet = 999)
cat("\nRandtest p-values (Region_NS/Site):\n")
print(am_ns_test)

capture.output(am_ns,      file = file.path(OUTDIR, "amova_regionNS_within_site.txt"))
capture.output(am_ns_test, file = file.path(OUTDIR, "amova_regionNS_randtest.txt"))

pdf(file.path(OUTDIR, "amova_regionNS_randtest.pdf"), width = 7, height = 5)
plot(am_ns_test)
dev.off()

cat("DONE AMOVA. Outputs in: ", OUTDIR, "\n", sep = "")
