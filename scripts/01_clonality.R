# filename: scripts/01_clonality.R
############################################################
# scripts/01_clonality.R
# Clonality summaries only (uses saved df_ids with MLG/MLL)
############################################################

options(repos = c(CRAN = "https://cloud.r-project.org"))
pkgs <- c("dplyr")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
suppressPackageStartupMessages(library(dplyr))

source("scripts/_load_objects.R")

OUTDIR <- file.path(RUN_OUT, "clonality_only")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

if (!all(c("Site", "MLG", "MLL") %in% names(df_ids))) {
  stop(
    "df_ids must contain Site, MLG, and MLL columns. ",
    "Run scripts/01_clone_correction.R first to populate clone IDs."
  )
}

clonality_by_site <- df_ids %>%
  group_by(Site) %>%
  summarise(
    n_inds = n(),
    n_MLG  = n_distinct(MLG),
    n_MLL  = n_distinct(MLL),
    genotypic_richness_MLG = n_MLG / n_inds,
    genotypic_richness_MLL = n_MLL / n_inds,
    clonal_richness_R = ifelse(n_inds > 1, (n_MLL - 1) / (n_inds - 1), NA_real_),
    .groups = "drop"
  ) %>%
  arrange(desc(n_inds))

overall <- data.frame(
  n_inds = nrow(df_ids),
  n_MLG  = dplyr::n_distinct(df_ids$MLG),
  n_MLL  = dplyr::n_distinct(df_ids$MLL),
  genotypic_richness_MLG = dplyr::n_distinct(df_ids$MLG) / nrow(df_ids),
  genotypic_richness_MLL = dplyr::n_distinct(df_ids$MLL) / nrow(df_ids),
  clonal_richness_R = ifelse(nrow(df_ids) > 1,
                             (dplyr::n_distinct(df_ids$MLL) - 1) / (nrow(df_ids) - 1),
                             NA_real_)
)

print(clonality_by_site)
print(overall)

write.csv(clonality_by_site, file.path(OUTDIR, "clonality_by_site.csv"), row.names = FALSE)
write.csv(overall,           file.path(OUTDIR, "clonality_overall.csv"), row.names = FALSE)

cat("DONE clonality. Outputs in: ", OUTDIR, "\n", sep = "")

