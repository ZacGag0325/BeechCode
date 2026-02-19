# filename: scripts/04_amova.R
############################################################
# scripts/04_amova.R
# AMOVA among Sites (+ permutation test) using clone-corrected dataset
############################################################

options(repos = c(CRAN = "https://cloud.r-project.org"))

pkgs <- c("adegenet", "poppr", "ade4", "here")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)

suppressPackageStartupMessages(library(adegenet))
suppressPackageStartupMessages(library(poppr))
suppressPackageStartupMessages(library(ade4))
suppressPackageStartupMessages(library(here))

# -----------------------------------------------------------------------------
# Load gi/df_ids only if not already present in environment
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

if (!exists("gi_mll", inherits = FALSE) || !exists("df_ids", inherits = FALSE)) {
  need <- c("gi_mll.rds", "df_ids.rds")
  missing <- need[!file.exists(file.path(OBJ_DIR, need))]
  if (length(missing) > 0) {
    stop(
      "Missing saved objects in ", OBJ_DIR, ". Missing file(s): ",
      paste(missing, collapse = ", "),
      ". Ensure gi_mll and df_ids are in memory or create these files first."
    )
  }
  
  if (!exists("gi_mll", inherits = FALSE)) {
    gi_mll <- readRDS(file.path(OBJ_DIR, "gi_mll.rds"))
  }
  if (!exists("df_ids", inherits = FALSE)) {
    df_ids <- readRDS(file.path(OBJ_DIR, "df_ids.rds"))
  }
}

OUTDIR <- file.path(RUN_OUT, "amova_only")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# Validate columns and align individuals between gi and df_ids
# -----------------------------------------------------------------------------
id_candidates <- c("ind", "SampleID", "Nom_Labo_Échantillons")
id_col <- id_candidates[id_candidates %in% names(df_ids)][1]
if (is.na(id_col)) {
  stop(
    "No individual ID column found in df_ids. Expected one of: ",
    paste(id_candidates, collapse = ", ")
  )
}

for (nm in c("Site", "Latitude")) {
  if (!nm %in% names(df_ids)) stop("df_ids must contain a '", nm, "' column.")
}

gi_use <- gi_mll
ids_gi <- indNames(gi_use)
ids_df <- as.character(df_ids[[id_col]])

if (anyDuplicated(ids_df) > 0) {
  dups <- unique(ids_df[duplicated(ids_df)])
  stop("Duplicate individual ID(s) found in df_ids (", id_col, "): ", paste(dups, collapse = ", "))
}

if (!setequal(ids_gi, ids_df)) {
  only_gi <- setdiff(ids_gi, ids_df)
  only_df <- setdiff(ids_df, ids_gi)
  stop(
    "Individual IDs mismatch between gi and df_ids. ",
    "Only in gi: ", if (length(only_gi)) paste(only_gi, collapse = ", ") else "none", "; ",
    "Only in df_ids: ", if (length(only_df)) paste(only_df, collapse = ", ") else "none"
  )
}

if (!identical(ids_gi, ids_df)) {
  reord <- match(ids_gi, ids_df)
  if (any(is.na(reord))) stop("Failed to reorder df_ids to match indNames(gi).")
  df_ids <- df_ids[reord, , drop = FALSE]
  ids_df <- as.character(df_ids[[id_col]])
}

if (!identical(ids_gi, ids_df)) {
  stop("Unresolvable mismatch after reordering: df_ids IDs do not match indNames(gi).")
}

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
  row.names = ids_gi,
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
