# scripts/00_run_all.R
############################################################
# Master runner: BeechCode genetics pipeline
# - Executes scripts in fixed order
# - Stops on first error
# - Anchors paths at project root
############################################################

find_project_root <- function() {
  candidates <- c(getwd(), normalizePath(file.path(getwd(), ".."), mustWork = FALSE))
  cmd_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
  if (length(cmd_file) > 0 && nzchar(cmd_file[1])) {
    candidates <- c(candidates, dirname(normalizePath(cmd_file[1], mustWork = FALSE)))
  }
  for (start in unique(candidates)) {
    cur <- normalizePath(start, mustWork = FALSE)
    repeat {
      if (file.exists(file.path(cur, "scripts", "00_run_all.R"))) return(cur)
      parent <- dirname(cur)
      if (identical(parent, cur)) break
      cur <- parent
    }
  }
  stop("[00_run_all] Cannot find project root containing scripts/00_run_all.R")
}

PROJECT_ROOT <- find_project_root()
setwd(PROJECT_ROOT)

cat("[00_run_all] Project root:", PROJECT_ROOT, "\n")

run_script <- function(script_name) {
  path <- file.path(PROJECT_ROOT, "scripts", script_name)
  if (!file.exists(path)) stop("[00_run_all] Missing script: ", path)
  cat("\n--- Running:", script_name, "\n")
  source(path, local = FALSE, encoding = "UTF-8")
  cat("--- Completed:", script_name, "\n")
}

for (d in c(
  file.path(PROJECT_ROOT, "outputs"),
  file.path(PROJECT_ROOT, "outputs", "tables"),
  file.path(PROJECT_ROOT, "outputs", "figures"),
  file.path(PROJECT_ROOT, "outputs", "matrices"),
  file.path(PROJECT_ROOT, "outputs", "tables", "supplementary"),
  file.path(PROJECT_ROOT, "outputs", "figures", "supplementary")
)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

run_script("00_master_pipeline.R")

obj_dir <- file.path(PROJECT_ROOT, "outputs", "v1", "objects")
required_obj <- c("gi.rds", "gi_mll.rds", "df_ids.rds", "meta.rds")
missing_obj <- required_obj[!file.exists(file.path(obj_dir, required_obj))]
if (length(missing_obj) > 0) {
  stop("[00_run_all] Object build incomplete. Missing: ", paste(missing_obj, collapse = ", "))
}
cat("[00_run_all] Verified canonical objects in", obj_dir, "\n")

run_script("01_clonality.R")
run_script("07_allelic_richness.R")
run_script("02_hwe.R")
run_script("05_pca_dapc.R")
run_script("04_amova.R")
run_script("06_distance_matrices.R")
run_script("11_isolation_by_distance.R")
run_script("03_structure_all_plots_S2N_bySite.R")

cat("\n[00_run_all] Pipeline completed successfully.\n")