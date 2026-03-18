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

clonality_by_site <- read.csv(clonality_file, stringsAsFactors = FALSE) %>%
  filter(Level == "site") %>%
  transmute(
    Site = as.character(Site),
    N_individuals = as.numeric(N_individuals),
    N_MLG = as.numeric(N_MLG),
    N_MLL = as.numeric(N_MLL),
    Clonal_Richness_MLG = as.numeric(Clonal_Richness_MLG),
    Clonal_Richness_MLL = as.numeric(Clonal_Richness_MLL)
  )

site_genetic_summary <- clonality_by_site %>%
  left_join(heterozygosity_by_site, by = "Site") %>%
  left_join(allelic_richness_by_site, by = "Site") %>%
  arrange(Site) %>%
  select(
    Site,
    N_individuals,
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