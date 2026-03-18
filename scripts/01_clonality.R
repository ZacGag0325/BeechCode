# scripts/01_clonality.R
############################################################
# Clonality summary (full dataset: gi)
# Why gi here?
# - Clonality must be quantified on the full (non-clone-corrected) dataset.
# - We use exact multilocus genotype (MLG) identity for clone counting.
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

message("[01_clonality] Calculating clonality from exact MLG identity on gi...")

df_ids_cols <- resolve_df_ids_columns(df_ids, context = "[01_clonality]", require = TRUE)
id_col <- df_ids_cols$id_col
site_col <- df_ids_cols$site_col

# Exact multilocus genotype labels from full genotype data (gi)
mlg_raw <- tryCatch(poppr::mlg.vector(gi), error = function(e) as.integer(factor(poppr::mlg(gi))))
mlg_labels <- paste0("MLG_", as.integer(factor(mlg_raw)))

site_map <- setNames(as.character(df_ids[[site_col]]), normalize_id(df_ids[[id_col]]))
site_labels <- site_map[normalize_id(adegenet::indNames(gi))]
if (any(is.na(site_labels))) stop("Could not map all individuals to Site.")

clonality_df <- data.frame(
  Individual = adegenet::indNames(gi),
  Site = site_labels,
  MLG = mlg_labels,
  stringsAsFactors = FALSE
)

calc_R <- function(N, G) ifelse(N > 1, (G - 1) / (N - 1), NA_real_)

overall <- clonality_df %>%
  summarise(
    Level = "overall",
    Site = "ALL",
    N_individuals = n(),
    N_MLG = n_distinct(MLG)
  ) %>%
  mutate(
    Clonal_Richness_R = calc_R(N_individuals, N_MLG),
    Genotypic_Richness = N_MLG / N_individuals
  )

by_site <- clonality_df %>%
  group_by(Site) %>%
  summarise(
    Level = "site",
    N_individuals = n(),
    N_MLG = n_distinct(MLG),
    .groups = "drop"
  ) %>%
  mutate(
    Clonal_Richness_R = calc_R(N_individuals, N_MLG),
    Genotypic_Richness = N_MLG / N_individuals
  ) %>%
  select(Level, Site, N_individuals, N_MLG, Clonal_Richness_R, Genotypic_Richness)

clonality_summary <- bind_rows(overall, by_site)
out_file <- file.path(TABLES_DIR, "clonality_summary.csv")
write.csv(clonality_summary, out_file, row.names = FALSE)

assign_file <- file.path(TABLES_DIR, "clonality_individual_assignments.csv")
write.csv(clonality_df, assign_file, row.names = FALSE)

message("[01_clonality] Saved: ", out_file)
message("[01_clonality] Saved: ", assign_file)