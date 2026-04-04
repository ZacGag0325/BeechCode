# scripts/07_allelic_richness.R
############################################################
# Genetic diversity by site (clone-corrected genind: gi_mll)
# Required outputs:
# - outputs/tables/heterozygosity_by_site.csv
# - outputs/tables/allelic_richness_by_site.csv
# - outputs/tables/site_genetic_summary.csv
# Added outputs (Ho/He/FIS focused):
# - outputs/tables/heterozygosity_fis_by_site.csv
# - outputs/tables/heterozygosity_fis_overall.csv
# Optional Excel mirrors (.xlsx) are written when writexl is available.
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

# ------------------------------------------------------------
# New helper section: consistent output writing (CSV + optional XLSX)
# ------------------------------------------------------------
write_table_with_optional_excel <- function(df, csv_path) {
  write.csv(df, csv_path, row.names = FALSE)
  message("[07_allelic_richness] Saved: ", csv_path)
  
  # Optional Excel output to support spreadsheet workflows without adding
  # a hard dependency. If writexl is unavailable, CSV remains authoritative.
  if (requireNamespace("writexl", quietly = TRUE)) {
    xlsx_path <- sub("\\.csv$", ".xlsx", csv_path)
    writexl::write_xlsx(df, path = xlsx_path)
    message("[07_allelic_richness] Saved: ", xlsx_path)
  }
}

# ------------------------------------------------------------
# New helper section: safe means that return NA (not NaN) when all-missing
# ------------------------------------------------------------
safe_row_mean <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

# ------------------------------------------------------------
# New helper section: robust extractor for basic.stats()$overall
# ------------------------------------------------------------
# hierfstat::basic.stats() can return `overall` either as a named numeric
# vector or as a matrix/data.frame depending on version/context. This helper
# safely extracts a requested metric (`Ho`, `Hs`, `Fis`) in both cases.
extract_overall_stat <- function(overall_obj, stat_name) {
  if (is.null(overall_obj)) return(NA_real_)
  
  if (is.atomic(overall_obj) && !is.null(names(overall_obj))) {
    val <- overall_obj[[stat_name]]
    return(if (length(val) == 0) NA_real_ else as.numeric(val[1]))
  }
  
  if (is.matrix(overall_obj) || is.data.frame(overall_obj)) {
    rn <- rownames(overall_obj)
    cn <- colnames(overall_obj)
    
    if (!is.null(rn) && stat_name %in% rn) {
      return(safe_row_mean(as.numeric(overall_obj[stat_name, , drop = TRUE])))
    }
    if (!is.null(cn) && stat_name %in% cn) {
      return(safe_row_mean(as.numeric(overall_obj[, stat_name, drop = TRUE])))
    }
  }
  
  # Fallback for rare structures (e.g., unnamed vectors/lists of length 1).
  if (length(overall_obj) == 1 && is.finite(suppressWarnings(as.numeric(overall_obj)))) {
    return(as.numeric(overall_obj))
  }
  
  NA_real_
}

hf <- hierfstat::genind2hierfstat(gi_mll)
bs <- hierfstat::basic.stats(hf)

# ----------------------------
# 1) Ho / He / FIS by site
# ----------------------------
# Standard implementation for microsatellite genotypes:
# - hierfstat::basic.stats()
# - Ho: observed heterozygosity (mean across loci)
# - He: expected heterozygosity within populations (Hs; mean across loci)
# - FIS: inbreeding coefficient (mean across loci)
#
# Interpretation notes for downstream reading:
# - Positive FIS = heterozygote deficit
# - Negative FIS = heterozygote excess
# - Used with HWE outputs to assess whether Wahlund effect is plausible.
site_levels <- rownames(bs$Ho)
site_n_tbl <- table(as.character(adegenet::pop(gi_mll)))

heterozygosity_by_site <- data.frame(
  Site = site_levels,
  N = as.integer(site_n_tbl[site_levels]),
  Ho = apply(bs$Ho, 1, safe_row_mean),
  He = apply(bs$Hs, 1, safe_row_mean),
  FIS = apply(bs$Fis, 1, safe_row_mean),
  stringsAsFactors = FALSE
) %>%
  mutate(
    FIS_interpretation = case_when(
      is.na(FIS) ~ NA_character_,
      FIS > 0 ~ "heterozygote_deficit",
      FIS < 0 ~ "heterozygote_excess",
      TRUE ~ "no_deficit_or_excess"
    )
  )

# Rounded copy for reporting tables (full precision retained in-memory objects).
heterozygosity_fis_by_site <- heterozygosity_by_site %>%
  mutate(
    Ho = round(Ho, 4),
    He = round(He, 4),
    FIS = round(FIS, 4)
  )

# Preserve legacy filename for backward compatibility with existing workflow.
het_file_legacy <- file.path(TABLES_DIR, "heterozygosity_by_site.csv")
write_table_with_optional_excel(heterozygosity_fis_by_site, het_file_legacy)

# New explicit filename requested for Ho/He/FIS reporting.
het_file_new <- file.path(TABLES_DIR, "heterozygosity_fis_by_site.csv")
write_table_with_optional_excel(heterozygosity_fis_by_site, het_file_new)

# ----------------------------
# 1b) New overall Ho / He / FIS section (all individuals pooled)
# ----------------------------
# Clone correction policy note:
# This script operates on gi_mll (clone-corrected), consistent with HWE and
# diversity summaries in this pipeline. If a raw-ramet estimate is needed,
# run analogous calculations on gi in a separate, explicitly-labeled section.
overall_Ho <- extract_overall_stat(bs$overall, "Ho")
overall_He <- extract_overall_stat(bs$overall, "Hs")
overall_FIS <- extract_overall_stat(bs$overall, "Fis")

# Defensive fallback: if overall metrics are absent in a given hierfstat build,
# compute pooled summaries as means of site-level estimates.
if (is.na(overall_Ho)) overall_Ho <- safe_row_mean(as.numeric(heterozygosity_by_site$Ho))
if (is.na(overall_He)) overall_He <- safe_row_mean(as.numeric(heterozygosity_by_site$He))
if (is.na(overall_FIS)) overall_FIS <- safe_row_mean(as.numeric(heterozygosity_by_site$FIS))

overall_heterozygosity_fis <- data.frame(
  N = nInd(gi_mll),
  Ho = overall_Ho,
  He = overall_He,
  FIS = overall_FIS,
  FIS_interpretation = dplyr::case_when(
    is.na(overall_FIS) ~ NA_character_,
    overall_FIS > 0 ~ "heterozygote_deficit",
    overall_FIS < 0 ~ "heterozygote_excess",
    TRUE ~ "no_deficit_or_excess"
  ),
  stringsAsFactors = FALSE
)

overall_heterozygosity_fis_out <- overall_heterozygosity_fis %>%
  mutate(
    Ho = round(Ho, 4),
    He = round(He, 4),
    FIS = round(FIS, 4)
  )

overall_file <- file.path(TABLES_DIR, "heterozygosity_fis_overall.csv")
write_table_with_optional_excel(overall_heterozygosity_fis_out, overall_file)

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
write_table_with_optional_excel(allelic_richness_by_site, ar_file)

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
  left_join(heterozygosity_fis_by_site %>% select(Site, Ho, He, FIS), by = "Site") %>%
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
write_table_with_optional_excel(site_genetic_summary, summary_file)