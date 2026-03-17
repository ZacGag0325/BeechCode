# scripts/01_clonality.R
############################################################
# Clonality summary (uses full genind: gi)
# Required output:
# - outputs/tables/clonality_summary.csv
# Metrics:
# - MLG
# - MLL
# - clonal richness R
############################################################

suppressPackageStartupMessages({
  library(adegenet)
  library(poppr)
  library(dplyr)
})

source("scripts/_load_objects.R")

message("[01_clonality] Calculating MLG / MLL / clonal richness...")

df_ids_cols <- resolve_df_ids_columns(df_ids, context = "[01_clonality]", require = TRUE)
id_col <- df_ids_cols$id_col
site_col <- df_ids_cols$site_col
mll_col <- resolve_col_ci(df_ids, c("mll"))


# MLG labels from full genotype data (gi)
detected_mlg <- tryCatch(poppr::mlg.vector(gi), error = function(e) as.integer(factor(poppr::mlg(gi))))
mlg_labels <- paste0("MLG_", as.integer(factor(detected_mlg)))

# MLL labels from df_ids when available; fallback to MLG labels
if (!is.na(mll_col)) {
  mll_map <- setNames(as.character(df_ids[[mll_col]]), normalize_id(df_ids[[id_col]]))
  mll_labels <- mll_map[normalize_id(adegenet::indNames(gi))]
  if (any(is.na(mll_labels))) {
    warning("Missing MLL labels for some individuals; filling with MLG labels.")
    mll_labels[is.na(mll_labels)] <- mlg_labels[is.na(mll_labels)]
  }
} else {
  warning("MLL column not found in df_ids; using MLG labels as fallback.")
  mll_labels <- mlg_labels
}

site_map <- setNames(as.character(df_ids[[site_col]]), normalize_id(df_ids[[id_col]]))
site_labels <- site_map[normalize_id(adegenet::indNames(gi))]
if (any(is.na(site_labels))) stop("Could not map all individuals to Site.")

clonality_df <- data.frame(
  Individual = adegenet::indNames(gi),
  Site = site_labels,
  MLG = mlg_labels,
  MLL = mll_labels,
  stringsAsFactors = FALSE
)

calc_R <- function(N, G) ifelse(N > 1, (G - 1) / (N - 1), NA_real_)

overall <- clonality_df %>%
  summarise(
    Level = "overall",
    Site = "ALL",
    N_individuals = n(),
    MLG = n_distinct(MLG),
    MLL = n_distinct(MLL)
  ) %>%
  mutate(
    Clonal_Richness_R = calc_R(N_individuals, MLL),
    Genotypic_Richness_MLG = MLG / N_individuals,
    Genotypic_Richness_MLL = MLL / N_individuals
  )

by_site <- clonality_df %>%
  group_by(Site) %>%
  summarise(
    Level = "site",
    N_individuals = n(),
    MLG = n_distinct(MLG),
    MLL = n_distinct(MLL),
    .groups = "drop"
  ) %>%
  mutate(
    Clonal_Richness_R = calc_R(N_individuals, MLL),
    Genotypic_Richness_MLG = MLG / N_individuals,
    Genotypic_Richness_MLL = MLL / N_individuals
  ) %>%
  select(Level, Site, N_individuals, MLG, MLL, Clonal_Richness_R, Genotypic_Richness_MLG, Genotypic_Richness_MLL)

clonality_summary <- bind_rows(overall, by_site)
out_file <- file.path(TABLES_DIR, "clonality_summary.csv")
write.csv(clonality_summary, out_file, row.names = FALSE)

message("[01_clonality] Saved: ", out_file)