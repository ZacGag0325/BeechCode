# --- HWE (robust, with reasons) ------------------------------------------
suppressPackageStartupMessages({
  library(adegenet)
  library(dplyr)
  library(tidyr)
  library(pegas)
})

# gi must be a genind, with pop(gi) = your Site factor
stopifnot(inherits(gi, "genind"))
stopifnot(!is.null(pop(gi)))

hw_one_genind <- function(gi_sub, label = "POP", B = 10000, min_n = 8) {
  # Convert to pegas loci object (genotypes per individual)
  loci_obj <- tryCatch(pegas::as.loci(gi_sub), error = function(e) e)
  if (inherits(loci_obj, "error")) {
    return(tibble(
      Site = label, Locus = locNames(gi_sub),
      p_value = NA_real_, reason = paste("as.loci error:", loci_obj$message)
    ))
  }
  
  loci_names <- colnames(loci_obj)
  
  # Pre-checks per locus: n typed + number of alleles
  pre <- lapply(loci_names, function(L) {
    g <- loci_obj[, L]
    n_typed <- sum(!is.na(g))
    # count alleles from genotype strings
    alleles <- unique(unlist(strsplit(as.character(g[!is.na(g)]), "/")))
    n_alleles <- length(alleles)
    
    tibble(Locus = L, n_typed = n_typed, n_alleles = n_alleles)
  }) |> bind_rows()
  
  # Run HWE test (hw.test returns a result per locus for loci objects)
  hw <- tryCatch(pegas::hw.test(loci_obj, B = B), error = function(e) e)
  
  if (inherits(hw, "error")) {
    out <- pre %>%
      mutate(
        Site = label,
        p_value = NA_real_,
        reason = paste("hw.test error:", hw$message)
      ) %>%
      select(Site, Locus, p_value, reason, n_typed, n_alleles)
    return(out)
  }
  
  # Extract p-values safely
  # For pegas::hw.test on loci, p-values are typically in hw$p.value named by locus
  pvec <- tryCatch(hw$p.value, error = function(e) NULL)
  
  out <- pre %>%
    mutate(
      Site = label,
      p_value = if (!is.null(pvec)) as.numeric(pvec[Locus]) else NA_real_,
      reason = case_when(
        n_typed < min_n ~ paste0("too_few_typed (n_typed<", min_n, ")"),
        n_alleles < 2 ~ "monomorphic_or_no_data",
        is.na(p_value) ~ "p_value_missing_from_test",
        TRUE ~ NA_character_
      )
    ) %>%
    select(Site, Locus, p_value, reason, n_typed, n_alleles)
  
  out
}

# 1) By site
sites <- levels(pop(gi))
hwe_by_site <- lapply(sites, function(s) hw_one_genind(gi[pop(gi) == s, ], label = s)) |> bind_rows()

# 2) Global (MLG vs MLL versions if you have them; otherwise just run gi)
# If you have gi_mlg and gi_mll as separate genind objects, run hw_one_genind on each.
hwe_global <- hw_one_genind(gi, label = "GLOBAL")

# BH adjust within each Site (site-by-locus table)
hwe_by_site <- hwe_by_site %>%
  group_by(Site) %>%
  mutate(
    p_adj_bh = p.adjust(p_value, method = "BH"),
    sig_p05  = !is.na(p_value) & p_value < 0.05,
    sig_bh05 = !is.na(p_adj_bh) & p_adj_bh < 0.05
  ) %>%
  ungroup()

# Summaries
hwe_site_summary <- hwe_by_site %>%
  group_by(Site) %>%
  summarise(
    loci_tested = sum(!is.na(p_value)),
    n_sig_p05   = sum(sig_p05, na.rm = TRUE),
    n_sig_bh05  = sum(sig_bh05, na.rm = TRUE),
    n_NA        = sum(is.na(p_value)),
    .groups = "drop"
  )

hwe_global_by_locus <- hwe_global %>%
  mutate(
    p_adj_bh = p.adjust(p_value, method = "BH"),
    sig_p05  = !is.na(p_value) & p_value < 0.05,
    sig_bh05 = !is.na(p_adj_bh) & p_adj_bh < 0.05
  )

# Print the “why are these NA?” diagnostics first
cat("\n--- HWE NA reasons (by site) ---\n")
print(hwe_by_site %>% filter(is.na(p_value)) %>% count(reason, sort = TRUE), n = 50)

cat("\n--- Site summary ---\n")
print(hwe_site_summary, n = Inf)
