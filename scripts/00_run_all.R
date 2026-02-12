############################################################
# scripts/00_run_all.R
# One-click runner: master pipeline + all analysis modules
############################################################

options(repos = c(CRAN = "https://cloud.r-project.org"))

# ---- Robust project root ----
suppressPackageStartupMessages(library(here))

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

cat("Project root:", PROJECT_ROOT, "\n")

# ---- Helper: safe source ----
src <- function(path) {
  f <- file.path(PROJECT_ROOT, path)
  if (!file.exists(f)) stop("Missing script: ", f)
  cat("\n========== RUNNING:", path, "==========\n")
  source(f, local = FALSE, encoding = "UTF-8")
}

# ---- 0) Run master pipeline (build gi + outputs + save objects) ----
src("00_master_pipeline.R")

# ---- 1) Run modular analyses (they load outputs/v1/objects/*.rds) ----
src("scripts/01_clonality.R")
src("scripts/02_hwe.R")
src("scripts/03_structure_all_plots_S2N_bySite.R")
src("scripts/04_amova.R")
src("scripts/05_pca_dapc.R")

cat("\n✅ ALL DONE.\nCheck outputs in: ", file.path(PROJECT_ROOT, "outputs", "v1"), "\n", sep = "")
