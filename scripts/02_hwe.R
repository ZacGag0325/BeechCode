# scripts/02_hwe.R
############################################################
# Hardy-Weinberg Equilibrium tests (clone-corrected: gi_mll)
# Required outputs:
# - outputs/tables/hwe_by_locus.csv
# - outputs/tables/hwe_by_site.csv
# - outputs/tables/hwe_by_locus_within_site.csv
# Supplementary duplicates:
# - outputs/tables/supplementary/*.csv
# Includes multiple testing correction (BH + Bonferroni)
############################################################

suppressPackageStartupMessages({
  library(adegenet)
  library(pegas)
  library(dplyr)
})

source("scripts/_load_objects.R")

message("[02_hwe] Running HWE tests on gi_mll (clone-corrected)...")

run_hwe_by_locus <- function(gobj, n_permutations = 1000) {
  gdf <- adegenet::genind2df(gobj, sep = "/")
  loci <- setdiff(names(gdf), "pop")
  
  out <- lapply(loci, function(loc) {
    loc_obj <- pegas::loci(gdf[loc])
    pval <- tryCatch(
      as.numeric(pegas::hw.test(loc_obj, B = n_permutations)[1, "Pr(Prob)"]),
      error = function(e) NA_real_
    )
    data.frame(Locus = loc, p_value = pval, stringsAsFactors = FALSE)
  })
  
  bind_rows(out)
}

apply_corrections <- function(df) {
  df %>%
    mutate(
      p_BH = p.adjust(p_value, method = "BH"),
      p_Bonferroni = p.adjust(p_value, method = "bonferroni")
    )
}

combine_site_pvalues_fisher <- function(site_locus_df) {
  site_locus_df %>%
    group_by(Site) %>%
    summarise(
      n_loci_tested = sum(!is.na(p_value)),
      fisher_chisq = -2 * sum(log(p_value), na.rm = TRUE),
      fisher_df = 2 * n_loci_tested,
      p_value = pchisq(fisher_chisq, fisher_df, lower.tail = FALSE),
      .groups = "drop"
    ) %>%
    apply_corrections()
}

save_main_and_supp <- function(df, filename) {
  main_file <- file.path(TABLES_DIR, filename)
  supp_file <- file.path(TABLES_SUPP_DIR, filename)
  write.csv(df, main_file, row.names = FALSE)
  write.csv(df, supp_file, row.names = FALSE)
  message("[02_hwe] Saved: ", main_file)
  message("[02_hwe] Saved: ", supp_file)
}

# ----------------------------
# 1) By locus (overall)
# ----------------------------
hwe_by_locus <- run_hwe_by_locus(gi_mll, n_permutations = 1000) %>%
  apply_corrections()
save_main_and_supp(hwe_by_locus, "hwe_by_locus.csv")

# ----------------------------
# 2) By locus within site
# ----------------------------
sites <- levels(adegenet::pop(gi_mll))
site_locus_list <- lapply(sites, function(s) {
  g_sub <- gi_mll[adegenet::pop(gi_mll) == s, , drop = FALSE]
  if (adegenet::nInd(g_sub) < 5) return(NULL)
  tmp <- run_hwe_by_locus(g_sub, n_permutations = 1000)
  tmp$Site <- s
  tmp
})

hwe_by_locus_within_site <- bind_rows(site_locus_list) %>%
  group_by(Site) %>%
  apply_corrections() %>%
  ungroup() %>%
  select(Site, Locus, p_value, p_BH, p_Bonferroni)
save_main_and_supp(hwe_by_locus_within_site, "hwe_by_locus_within_site.csv")

# ----------------------------
# 3) By site (overall)
# ----------------------------
hwe_by_site <- combine_site_pvalues_fisher(hwe_by_locus_within_site)
save_main_and_supp(hwe_by_site, "hwe_by_site.csv")