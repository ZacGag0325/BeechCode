# scripts/00_run_all.R
############################################################
# Master runner: population genetics pipeline ONLY
# - Runs scripts in a fixed sequence
# - Stops on first error
# - Excludes non-genetic scripts
############################################################

message("========================================")
message("Starting BeechCode genetics-only pipeline")
message("========================================")

run_script <- function(path) {
  message("\n--- Running: ", path)
  if (!file.exists(path)) stop("Missing script: ", path)
  source(path, local = FALSE, encoding = "UTF-8")
  message("--- Completed: ", path)
}

# Ensure base output folders exist early
for (d in c("outputs", "outputs/tables", "outputs/figures", "outputs/matrices", "outputs/tables/supplementary", "outputs/figures/supplementary")) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# 0) Build required objects (gi, gi_mll, df_ids, meta)
run_script("00_master_pipeline.R")

# 1) Population genetics analyses (publication workflow)
run_script("scripts/01_clonality.R")
run_script("scripts/07_allelic_richness.R")
run_script("scripts/02_hwe.R")
run_script("scripts/05_pca_dapc.R")
run_script("scripts/04_amova.R")
run_script("scripts/06_distance_matrices.R")
run_script("scripts/11_isolation_by_distance.R")
run_script("scripts/03_structure_all_plots_S2N_bySite.R")

message("\n========================================")
message("Genetics pipeline completed successfully")
message("Outputs written to:")
message(" - outputs/tables")
message(" - outputs/figures")
message(" - outputs/matrices")
message(" - outputs/tables/supplementary")
message(" - outputs/figures/supplementary")
message("========================================")