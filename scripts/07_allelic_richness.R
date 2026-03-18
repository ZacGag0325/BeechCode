# scripts/07_allelic_richness.R
############################################################
# Genetic diversity by site (clone-corrected genind: gi_mll)
# Required outputs:
# - outputs/tables/heterozygosity_by_site.csv
# - outputs/tables/allelic_richness_by_site.csv
# - outputs/tables/site_genetic_summary.csv
# Metrics:
# - Ho, He, FIS
# - Allelic richness (+ SE)
# - merged with clonality summary (N, N_MLG, N_MLL, R)
############################################################

suppressPackageStartupMessages({
  library(adegenet)
  library(hierfstat)
  library(dplyr)
})

source("scripts/_load_objects.R")

message("[07_allelic_richness] Calculating heterozygosity/FIS and allelic richness...")
validate_columns(df_ids_mll, c("ind_id", "Site"), df_name = "07_allelic_richness df_ids_mll")
if (!all(adegenet::indNames(gi_mll) == df_ids_mll$ind_id)) stop("[07_allelic_richness] gi_mll and df_ids_mll are not aligned.")

hf <- hierfstat::genind2hierfstat(gi_mll)
bs <- hierfstat::basic.stats(hf)

# ----------------------------
# 1) Ho / He / FIS by site
# ----------------------------
heterozygosity_by_site <- data.frame(
  Site = rownames(bs$Ho),
  Ho = rowMeans(bs$Ho, na.rm = TRUE),
  He = rowMeans(bs$Hs, na.rm = TRUE),
  FIS = rowMeans(bs$Fis, na.rm = TRUE),
  stringsAsFactors = FALSE
)

het_file <- file.path(TABLES_DIR, "heterozygosity_by_site.csv")
write.csv(heterozygosity_by_site, het_file, row.names = FALSE)
message("[07_allelic_richness] Saved: ", het_file)

# ----------------------------
# 2) Allelic richness + SE
# ----------------------------
pop_sizes <- table(adegenet::pop(gi_mll))
min_n <- min(pop_sizes)
if (is.na(min_n) || min_n < 2) {
  stop("At least one site has fewer than 2 individuals; cannot compute robust allelic richness.")
}

ar <- hierfstat::allelic.richness(hf, min.n = min_n)
Ar <- ar$Ar

allelic_richness_by_site <- data.frame(
  Site = colnames(Ar),
  Allelic_Richness = as.numeric(colMeans(Ar, na.rm = TRUE)),
  Allelic_Richness_SE = as.numeric(apply(Ar, 2, sd, na.rm = TRUE) / sqrt(nrow(Ar))),
  stringsAsFactors = FALSE
)

ar_file <- file.path(TABLES_DIR, "allelic_richness_by_site.csv")
write.csv(allelic_richness_by_site, ar_file, row.names = FALSE)
message("[07_allelic_richness] Saved: ", ar_file)

# ----------------------------
# 3) Site-level paper-ready summary
# ----------------------------
# We merge diversity metrics (clone-corrected gi_mll) with clonality metrics (from gi).
clonality_file <- file.path(TABLES_DIR, "clonality_summary.csv")
if (!file.exists(clonality_file)) {
  stop("Missing clonality_summary.csv. Run scripts/01_clonality.R first.")
}

clonality_raw <- read.csv(clonality_file, stringsAsFactors = FALSE)
validate_columns(clonality_raw, c("Level", "Site"), df_name = "clonality_summary.csv")

if (!"N" %in% names(clonality_raw) && "N_individuals" %in% names(clonality_raw)) {
  clonality_raw$N <- as.numeric(clonality_raw$N_individuals)
}
if (!"G" %in% names(clonality_raw) && "N_MLL" %in% names(clonality_raw)) {
  clonality_raw$G <- as.numeric(clonality_raw$N_MLL)
}
if (!"Clonal_Richness_R" %in% names(clonality_raw) && all(c("N", "G") %in% names(clonality_raw))) {
  clonality_raw$Clonal_Richness_R <- ifelse(clonality_raw$N > 1, (clonality_raw$G - 1) / (clonality_raw$N - 1), NA_real_)
}
validate_columns(clonality_raw, c("N", "G", "Clonal_Richness_R"), df_name = "clonality_summary.csv")

clonality_by_site <- clonality_raw %>%
  filter(Level == "site") %>%
  transmute(
    Site = as.character(Site),
    N = as.numeric(N),
    G = as.numeric(G),
    Clonal_Richness_R = as.numeric(Clonal_Richness_R),
    N_MLG = if ("N_MLG" %in% names(clonality_raw)) as.numeric(N_MLG) else NA_real_,
    N_MLL = if ("N_MLL" %in% names(clonality_raw)) as.numeric(N_MLL) else as.numeric(G),
    Clonal_Richness_MLG = if ("Clonal_Richness_MLG" %in% names(clonality_raw)) as.numeric(Clonal_Richness_MLG) else NA_real_,
    Clonal_Richness_MLL = if ("Clonal_Richness_MLL" %in% names(clonality_raw)) as.numeric(Clonal_Richness_MLL) else as.numeric(Clonal_Richness_R)
  )

site_genetic_summary <- clonality_by_site %>%
  left_join(heterozygosity_by_site, by = "Site") %>%
  left_join(allelic_richness_by_site, by = "Site") %>%
  arrange(Site) %>%
  select(
    Site,
    N,
    G,
    Clonal_Richness_R,
    N_MLG,
    N_MLL,
    Clonal_Richness_MLG,
    Clonal_Richness_MLL,
    Ho,
    He,
    FIS,
    Allelic_Richness,
    Allelic_Richness_SE
  )

summary_file <- file.path(TABLES_DIR, "site_genetic_summary.csv")
write.csv(site_genetic_summary, summary_file, row.names = FALSE)
message("[07_allelic_richness] Saved: ", summary_file)