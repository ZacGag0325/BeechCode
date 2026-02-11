############################################################
# scripts/_load_objects.R
# Loads gi/df_ids/meta/gi_mll from outputs/v1/objects
############################################################

suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(adegenet))

FALLBACK_ROOT <- file.path(path.expand("~"), "Desktop", "BeechCode")
candidate_root <- here::here()

if (dir.exists(file.path(candidate_root, "data", "raw"))) {
  PROJECT_ROOT <- candidate_root
} else if (dir.exists(file.path(FALLBACK_ROOT, "data", "raw"))) {
  PROJECT_ROOT <- FALLBACK_ROOT
  setwd(PROJECT_ROOT)
} else {
  stop("Can't find project root. Open BeechCode.Rproj or set FALLBACK_ROOT.")
}

RUN_TAG <- "v1"
RUN_OUT <- file.path(PROJECT_ROOT, "outputs", RUN_TAG)
OBJ_DIR <- file.path(RUN_OUT, "objects")

need <- c("gi.rds", "df_ids.rds", "meta.rds", "gi_mll.rds")
missing <- need[!file.exists(file.path(OBJ_DIR, need))]
if (length(missing) > 0) {
  stop("Missing saved objects in: ", OBJ_DIR, "\nMissing: ", paste(missing, collapse = ", "),
       "\nFix: run 00_master_pipeline.R (with SAVE SHARED OBJECTS block) first.")
}

gi     <- readRDS(file.path(OBJ_DIR, "gi.rds"))
df_ids <- readRDS(file.path(OBJ_DIR, "df_ids.rds"))
meta   <- readRDS(file.path(OBJ_DIR, "meta.rds"))
gi_mll <- readRDS(file.path(OBJ_DIR, "gi_mll.rds"))

cat("Loaded objects:\n")
cat(" - nInd(gi)     =", nInd(gi), "\n")
cat(" - nLoc(gi)     =", nLoc(gi), "\n")
cat(" - df_ids rows  =", nrow(df_ids), "\n")
cat(" - meta rows    =", nrow(meta), "\n")
cat(" - nInd(gi_mll) =", nInd(gi_mll), "\n")
