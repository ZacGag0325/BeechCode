# scripts/04_amova.R
############################################################
# AMOVA (clone-corrected: gi_mll)
# Why gi_mll here?
# - AMOVA is population-level inference and assumes independent genotypes.
# Outputs:
# - outputs/tables/amova_results.csv
# - outputs/tables/supplementary/amova_randtest_summary.csv
############################################################

suppressPackageStartupMessages({
  library(poppr)
  library(ade4)
  library(dplyr)
})

source("scripts/_load_objects.R")

message("[04_amova] Running AMOVA on gi_mll...")

amova_fit <- poppr::poppr.amova(gi_mll, ~pop)
randtest_fit <- ade4::randtest(amova_fit, nrepet = 999)

components <- as.data.frame(amova_fit$componentsofcovariance)
components$Source <- rownames(components)
rownames(components) <- NULL
names(components)[1] <- "Sigma"

phi_stats <- as.data.frame(amova_fit$statphi)
phi_stats$Source <- rownames(phi_stats)
rownames(phi_stats) <- NULL
names(phi_stats)[1] <- "Phi"

amova_results <- dplyr::full_join(components, phi_stats, by = "Source") %>%
  mutate(
    Randtest_observed = as.numeric(randtest_fit$obs),
    Randtest_p_value = as.numeric(randtest_fit$pvalue)
  ) %>%
  select(Source, everything())

amova_file <- file.path(TABLES_DIR, "amova_results.csv")
write.csv(amova_results, amova_file, row.names = FALSE)
message("[04_amova] Saved: ", amova_file)

amova_randtest_summary <- data.frame(
  test = "AMOVA_randtest",
  statistic = as.numeric(randtest_fit$obs),
  p_value = as.numeric(randtest_fit$pvalue),
  permutations = 999,
  stringsAsFactors = FALSE
)

supp_file <- file.path(TABLES_SUPP_DIR, "amova_randtest_summary.csv")
write.csv(amova_randtest_summary, supp_file, row.names = FALSE)
message("[04_amova] Saved: ", supp_file)