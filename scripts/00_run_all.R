# scripts/00_run_all.R
############################################################
# Master runner: population genetics pipeline ONLY
# - Runs scripts in a fixed sequence
# - Stops on first error
# - Excludes non-genetic scripts
############################################################

if (!requireNamespace("here", quietly = TRUE)) {
  stop("The 'here' package is required. Install it with install.packages('here').")
}

get_runner_path <- function() {
  frame_files <- vapply(sys.frames(), function(x) x$ofile %||% NA_character_, character(1))
  frame_files <- frame_files[!is.na(frame_files)]
  
  if (length(frame_files) > 0) {
    return(normalizePath(frame_files[[length(frame_files)]], mustWork = TRUE))
  }
  
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0) {
    return(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE))
  }
  
  stop("Unable to determine path to 00_run_all.R; cannot resolve project root robustly.")
}

`%||%` <- function(x, y) if (is.null(x)) y else x

runner_path <- get_runner_path()
project_root <- normalizePath(file.path(dirname(runner_path), ".."), mustWork = TRUE)

setwd(project_root)

# Reinitialize `here` after setting the working directory so it anchors correctly
if ("here" %in% loadedNamespaces()) {
  unloadNamespace("here")
}
if (!requireNamespace("here", quietly = TRUE)) {
  stop("The 'here' package is required. Install it with install.packages('here').")
}

cat("Project root detected:", project_root, "\n")
cat("here::here() resolves to:", here::here(), "\n")
cat("Runner script path:", runner_path, "\n")
cat("Working directory set to:", getwd(), "\n")

cat("========================================\n")
cat("Starting BeechCode genetics-only pipeline\n")
cat("========================================\n")

run_script <- function(script) {
  script_path <- here::here("scripts", script)
  
  if (!file.exists(script_path)) {
    stop("Missing script: ", script_path)
  }
  
  cat("\n--- Running:", script, "\n")
  source(script_path, local = FALSE, encoding = "UTF-8")
  cat("--- Completed:", script, "\n")
}

# Ensure base output folders exist early
for (d in c(
  here::here("outputs"),
  here::here("outputs", "tables"),
  here::here("outputs", "figures"),
  here::here("outputs", "matrices"),
  here::here("outputs", "tables", "supplementary"),
  here::here("outputs", "figures", "supplementary")
)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# 0) Build required objects (gi, gi_mll, df_ids, meta)
run_script("00_master_pipeline.R")

# 1) Population genetics analyses (publication workflow)
run_script("01_clonality.R")
run_script("07_allelic_richness.R")
run_script("02_hwe.R")
run_script("05_pca_dapc.R")
run_script("04_amova.R")
run_script("06_distance_matrices.R")
run_script("11_isolation_by_distance.R")
run_script("03_structure_all_plots_S2N_bySite.R")

cat("\n========================================\n")
cat("Genetics pipeline completed successfully\n")
cat("Outputs written to:\n")
cat(" - outputs/tables\n")
cat(" - outputs/figures\n")
cat(" - outputs/matrices\n")
cat(" - outputs/tables/supplementary\n")
cat(" - outputs/figures/supplementary\n")
cat("========================================\n")