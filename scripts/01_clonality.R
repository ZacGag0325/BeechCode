# scripts/01_clonality.R
############################################################
# Clonality summary (full dataset: gi)
# Why gi here?
# - Clonality must be quantified on the full (non-clone-corrected) dataset.
# - Exact MLG identity is reported for threshold-free clone counts.
# - Bruvo-based MLL identity is also reported for microsatellite-aware clone
#   lineages under the configured threshold.
#
# Outputs:
# - outputs/tables/clonality_summary.csv
# - outputs/tables/clonality_individual_assignments.csv
############################################################

suppressPackageStartupMessages({
  library(poppr)
  library(adegenet)
  library(dplyr)
})

source("scripts/_load_objects.R")

message("[01_clonality] Calculating clonality summaries for both MLG and Bruvo-based MLL on gi...")

DEFAULT_BRUVO_MLL_THRESHOLD <- 0.09
DEFAULT_BRUVO_ALGORITHM <- "farthest_neighbor"

find_synonym_col <- function(df, choices) {
  resolve_col_ci(df, choices)
}

compute_mlg_mll_from_gi <- function(gi, threshold = DEFAULT_BRUVO_MLL_THRESHOLD, algorithm = DEFAULT_BRUVO_ALGORITHM) {
  gc_mlg <- poppr::as.genclone(gi)
  mlg_raw <- tryCatch(poppr::mlg.vector(gc_mlg), error = function(e) as.integer(factor(poppr::mlg(gc_mlg))))
  mlg_labels <- paste0("MLG_", as.integer(factor(mlg_raw)))
  
  replen <- rep(2, adegenet::nLoc(gi))
  names(replen) <- adegenet::locNames(gi)
  
  gc_mll <- gc_mlg
  poppr::mlg.filter(
    gc_mll,
    distance = poppr::bruvo.dist,
    replen = replen,
    algorithm = algorithm
  ) <- threshold
  
  mll_raw <- poppr::mll(gc_mll)
  mll_labels <- paste0("MLL_", as.integer(factor(mll_raw)))
  
  list(
    MLG = mlg_labels,
    MLL = mll_labels,
    Bruvo_MLL_threshold = threshold,
    Bruvo_algorithm = algorithm
  )
}

recover_or_recompute_clonality_columns <- function(df_ids_tbl, gi, gi_mll) {
  out <- df_ids_tbl
  
  if (!all(c("MLG", "MLL") %in% names(out))) {
    mlg_syn <- find_synonym_col(out, c("MLG", "mlg", "MLG_id", "mlg_id", "genotype_id"))
    mll_syn <- find_synonym_col(out, c("MLL", "mll", "clone_id", "lineage", "multilocus_lineage"))
    
    if (!is.na(mlg_syn) && mlg_syn != "MLG") names(out)[names(out) == mlg_syn] <- "MLG"
    if (!is.na(mll_syn) && mll_syn != "MLL") names(out)[names(out) == mll_syn] <- "MLL"
    
    if (all(c("MLG", "MLL") %in% names(out))) {
      message("[01_clonality] Recovered MLG/MLL from existing columns.")
    }
  }
  
  if (!all(c("MLG", "MLL") %in% names(out))) {
    message("[01_clonality] Recomputed MLG/MLL from gi because columns were missing in df_ids.")
    cols <- resolve_df_ids_columns(out, context = "[01_clonality]", require = TRUE)
    recomputed <- compute_mlg_mll_from_gi(gi)
    key <- normalize_id(out[[cols$id_col]])
    gi_key <- normalize_id(adegenet::indNames(gi))
    idx <- match(key, gi_key)
    
    if (any(is.na(idx))) {
      stop("[01_clonality] Failed to align df_ids rows to gi while recomputing MLG/MLL.")
    }
    
    out$MLG <- recomputed$MLG[idx]
    out$MLL <- recomputed$MLL[idx]
    out$Bruvo_MLL_threshold <- recomputed$Bruvo_MLL_threshold
    out$Bruvo_algorithm <- recomputed$Bruvo_algorithm
  }
  
  if (!all(c("MLG", "MLL") %in% names(out))) {
    stop("[01_clonality] Failed to recover or recompute MLG/MLL columns.")
  }
  
  n_mll_df <- length(unique(out$MLL[!is.na(out$MLL)]))
  n_mll_obj <- adegenet::nInd(gi_mll)
  if (!identical(n_mll_df, n_mll_obj)) {
    message("[01_clonality] MLL count in df_ids did not match gi_mll; recomputing MLG/MLL from gi.")
    cols <- resolve_df_ids_columns(out, context = "[01_clonality]", require = TRUE)
    recomputed <- compute_mlg_mll_from_gi(gi)
    key <- normalize_id(out[[cols$id_col]])
    gi_key <- normalize_id(adegenet::indNames(gi))
    idx <- match(key, gi_key)
    out$MLG <- recomputed$MLG[idx]
    out$MLL <- recomputed$MLL[idx]
    out$Bruvo_MLL_threshold <- recomputed$Bruvo_MLL_threshold
    out$Bruvo_algorithm <- recomputed$Bruvo_algorithm
    n_mll_df <- length(unique(out$MLL[!is.na(out$MLL)]))
  }
  
  if (!identical(n_mll_df, n_mll_obj)) {
    stop("[01_clonality] df_ids$MLL is inconsistent with gi_mll after recovery/recomputation.")
  }
  
  out
}

df_ids <- recover_or_recompute_clonality_columns(df_ids, gi, gi_mll)

df_ids_cols <- resolve_df_ids_columns(df_ids, context = "[01_clonality]", require = TRUE)
id_col <- df_ids_cols$id_col
site_col <- df_ids_cols$site_col
threshold_col <- resolve_col_ci(df_ids, c("Bruvo_MLL_threshold", "bruvo_mll_threshold"))
algorithm_col <- resolve_col_ci(df_ids, c("Bruvo_algorithm", "bruvo_algorithm"))

id_key <- normalize_id(df_ids[[id_col]])
site_map <- setNames(as.character(df_ids[[site_col]]), id_key)
mlg_map <- setNames(as.character(df_ids[["MLG"]]), id_key)
mll_map <- setNames(as.character(df_ids[["MLL"]]), id_key)

ind_key <- normalize_id(adegenet::indNames(gi))
site_labels <- site_map[ind_key]
mlg_labels <- mlg_map[ind_key]
mll_labels <- mll_map[ind_key]

if (any(is.na(site_labels))) stop("[01_clonality] Could not map all individuals to Site.")
if (any(is.na(mlg_labels))) stop("[01_clonality] Could not map all individuals to MLG.")
if (any(is.na(mll_labels))) stop("[01_clonality] Could not map all individuals to MLL.")

bruvo_threshold <- if (!is.na(threshold_col)) unique(stats::na.omit(df_ids[[threshold_col]])) else numeric(0)
bruvo_algorithm <- if (!is.na(algorithm_col)) unique(stats::na.omit(df_ids[[algorithm_col]])) else character(0)

clonality_df <- data.frame(
  Individual = adegenet::indNames(gi),
  Site = site_labels,
  MLG = mlg_labels,
  MLL = mll_labels,
  stringsAsFactors = FALSE
)

calc_R <- function(N, G) ifelse(N > 1, (G - 1) / (N - 1), NA_real_)

add_clonality_metrics <- function(dat) {
  dat %>%
    mutate(
      Clonal_Richness_MLG = calc_R(N_individuals, N_MLG),
      Clonal_Richness_MLL = calc_R(N_individuals, N_MLL),
      Genotypic_Richness_MLG = N_MLG / N_individuals,
      Genotypic_Richness_MLL = N_MLL / N_individuals
    )
}

overall <- clonality_df %>%
  summarise(
    N_individuals = dplyr::n(),
    N_MLG = dplyr::n_distinct(MLG),
    N_MLL = dplyr::n_distinct(MLL)
  ) %>%
  add_clonality_metrics() %>%
  mutate(Level = "overall", Site = "ALL") %>%
  select(Level, Site, everything())

by_site <- clonality_df %>%
  group_by(Site) %>%
  summarise(
    N_individuals = dplyr::n(),
    N_MLG = dplyr::n_distinct(MLG),
    N_MLL = dplyr::n_distinct(MLL),
    .groups = "drop"
  ) %>%
  add_clonality_metrics() %>%
  mutate(Level = "site") %>%
  select(Level, Site, everything())

clonality_summary <- bind_rows(overall, by_site) %>%
  mutate(
    N = N_individuals,
    G = N_MLL,
    Clonal_Richness_R = Clonal_Richness_MLL
  )

if (length(bruvo_threshold) > 0) {
  clonality_summary$Bruvo_MLL_threshold <- bruvo_threshold[1]
}
if (length(bruvo_algorithm) > 0) {
  clonality_summary$Bruvo_algorithm <- bruvo_algorithm[1]
}

out_file <- file.path(TABLES_DIR, "clonality_summary.csv")
write.csv(clonality_summary, out_file, row.names = FALSE)

if (length(bruvo_threshold) > 0) {
  clonality_df$Bruvo_MLL_threshold <- bruvo_threshold[1]
}
if (length(bruvo_algorithm) > 0) {
  clonality_df$Bruvo_algorithm <- bruvo_algorithm[1]
}

assign_file <- file.path(TABLES_DIR, "clonality_individual_assignments.csv")
write.csv(clonality_df, assign_file, row.names = FALSE)

message("[01_clonality] Saved: ", out_file)
message("[01_clonality] Saved: ", assign_file)