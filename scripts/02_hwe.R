# scripts/02_hwe.R
############################################################
# Hardy-Weinberg Equilibrium tests (clone-corrected: gi_mll)
#
# Biological/statistical rationale:
# - HWE is evaluated on the clone-corrected object (gi_mll), not on gi,
#   because repeated ramets inflate genotype counts and can create
#   artifactual departures from equilibrium.
# - adegenet is retained for handling genind objects and site labels.
# - pegas::hw.test is used for the actual HWE test because pegas provides
#   an exact/Monte Carlo HWE workflow once locus data are converted into a
#   loci-compatible object. adegenet does not provide a dedicated exact HWE
#   testing function for genind objects, and hierfstat is retained elsewhere
#   in the pipeline for F-statistics rather than HWE testing.
#
# Main required output:
# - outputs/tables/hwe_by_site_by_locus.csv
#
# Compatibility outputs:
# - outputs/tables/hwe_by_locus.csv
# - outputs/tables/hwe_by_site.csv
# - outputs/tables/hwe_by_locus_within_site.csv
# - duplicated to outputs/tables/supplementary/
############################################################

suppressPackageStartupMessages({
  library(adegenet)
  library(pegas)
  library(dplyr)
})

source("scripts/_load_objects.R")

HWE_MONTE_CARLO_REPS <- 9999L
HWE_MIN_NON_MISSING_N <- 8L
HWE_MIN_UNIQUE_GENOTYPES <- 2L
HWE_ALPHA <- 0.05

message("[02_hwe] Running HWE tests on gi_mll (clone-corrected)...")
message("[02_hwe] genind handler: adegenet")
message("[02_hwe] HWE test function: pegas::hw.test via pegas::as.loci")
message("[02_hwe] Monte Carlo replicates per exact test: ", HWE_MONTE_CARLO_REPS)

if (!inherits(gi_mll, "genind")) {
  stop("[02_hwe] gi_mll must be a genind object. Run scripts/00_master_pipeline.R first.")
}

validate_columns(df_ids_mll, c("ind_id", "Site"), df_name = "[02_hwe] df_ids_mll")
if (!all(adegenet::indNames(gi_mll) == df_ids_mll$ind_id)) {
  stop("[02_hwe] gi_mll and df_ids_mll are not aligned.")
}

site_vec <- as.character(adegenet::pop(gi_mll))
if (all(is.na(site_vec)) || !any(nzchar(site_vec))) {
  stop("[02_hwe] pop(gi_mll) is missing; Site assignments are required for by-site HWE.")
}
if (!identical(site_vec, as.character(df_ids_mll$Site))) {
  stop("[02_hwe] pop(gi_mll) does not match df_ids_mll$Site. This would invalidate by-site HWE tests.")
}

# ------------------------------------------------------------
# Helper: robust allele-pair parsing from genind2df() strings
# ------------------------------------------------------------
split_genotype <- function(x) {
  x <- trimws(as.character(x))
  if (is.na(x) || x == "" || x %in% c("NA", "0", "0/0", "NA/NA", "-")) {
    return(c(NA_character_, NA_character_))
  }
  
  parts <- strsplit(x, "/", fixed = TRUE)[[1]]
  if (length(parts) != 2) {
    return(c(NA_character_, NA_character_))
  }
  
  parts <- trimws(parts)
  if (any(!nzchar(parts))) {
    return(c(NA_character_, NA_character_))
  }
  
  sort(parts)
}

# ------------------------------------------------------------
# Helper: identify missing/invalid genotype strings
# ------------------------------------------------------------
is_missing_genotype <- function(x) {
  x <- trimws(as.character(x))
  is.na(x) | x == "" | x %in% c("NA", "0", "0/0", "NA/NA", "-")
}

# ------------------------------------------------------------
# Helper: multiple-testing correction with NA-safe handling
# ------------------------------------------------------------
add_hwe_multipletest_columns <- function(df,
                                         p_col = "p_value_raw",
                                         alpha = HWE_ALPHA,
                                         include_legacy_aliases = FALSE) {
  if (!p_col %in% names(df)) {
    stop("[02_hwe] Requested p-value column not found: ", p_col)
  }
  
  pv <- suppressWarnings(as.numeric(df[[p_col]]))
  ok <- is.finite(pv) & !is.na(pv)
  
  p_bonf <- rep(NA_real_, length(pv))
  p_fdr <- rep(NA_real_, length(pv))
  
  if (sum(ok) > 0) {
    p_bonf[ok] <- p.adjust(pv[ok], method = "bonferroni")
    p_fdr[ok] <- p.adjust(pv[ok], method = "BH")
  }
  
  df$p_value_bonferroni <- p_bonf
  df$p_value_fdr <- p_fdr
  df$significant_raw <- !is.na(pv) & pv < alpha
  df$significant_bonferroni <- !is.na(p_bonf) & p_bonf < alpha
  df$significant_fdr <- !is.na(p_fdr) & p_fdr < alpha
  
  if (include_legacy_aliases) {
    df$p_adjust_bh <- df$p_value_fdr
    df$significant_bh <- df$significant_fdr
  }
  
  df
}

# ------------------------------------------------------------
# Helper: single-locus HWE test with defensive checks
# ------------------------------------------------------------
run_single_locus_hwe <- function(geno_vec,
                                 B = HWE_MONTE_CARLO_REPS,
                                 min_n = HWE_MIN_NON_MISSING_N,
                                 min_unique_genotypes = HWE_MIN_UNIQUE_GENOTYPES) {
  parsed <- t(vapply(geno_vec, split_genotype, character(2)))
  a1 <- parsed[, 1]
  a2 <- parsed[, 2]
  
  keep <- !(is.na(a1) | is.na(a2))
  n_non_missing <- sum(keep)
  raw_missing_fraction <- 1 - (n_non_missing / length(geno_vec))
  
  base_row <- list(
    N = as.integer(n_non_missing),
    missing_fraction = as.numeric(raw_missing_fraction),
    n_alleles = NA_integer_,
    n_unique_genotypes = NA_integer_,
    p_value = NA_real_,
    status = "not_tested",
    note = ""
  )
  
  if (length(geno_vec) == 0) {
    base_row$status <- "failed"
    base_row$note <- "empty_genotype_vector"
    return(base_row)
  }
  
  if (n_non_missing < min_n) {
    base_row$status <- "skipped"
    base_row$note <- sprintf("too_few_individuals_non_missing (N=%d < %d)", n_non_missing, min_n)
    return(base_row)
  }
  
  genotype_strings <- paste(a1[keep], a2[keep], sep = "/")
  n_genotypes <- dplyr::n_distinct(genotype_strings)
  base_row$n_unique_genotypes <- as.integer(n_genotypes)
  
  if (n_genotypes < min_unique_genotypes) {
    base_row$status <- "skipped"
    base_row$note <- sprintf("insufficient_genotype_diversity (unique_genotypes=%d)", n_genotypes)
    return(base_row)
  }
  
  alleles <- c(a1[keep], a2[keep])
  n_alleles <- dplyr::n_distinct(alleles)
  base_row$n_alleles <- as.integer(n_alleles)
  
  if (n_alleles < 2) {
    base_row$status <- "skipped"
    base_row$note <- "invariant_locus"
    return(base_row)
  }
  
  # pegas::hw.test operates on a loci-style data.frame/object.
  # We convert the single locus genotype strings into a factor column and
  # then into a loci object. This keeps the workflow compatible with a genind
  # source object without forcing a full file-format conversion.
  loc_df <- data.frame(Locus = factor(genotype_strings), stringsAsFactors = TRUE)
  
  error_msg <- NULL
  pval <- tryCatch({
    loci_obj <- pegas::as.loci(loc_df)
    hw_fit <- pegas::hw.test(loci_obj, B = B)
    
    if (is.list(hw_fit) && !is.null(hw_fit$p.value)) {
      as.numeric(hw_fit$p.value[1])
    } else if (is.matrix(hw_fit) || is.data.frame(hw_fit)) {
      hw_mat <- as.matrix(hw_fit)
      cn <- tolower(colnames(hw_mat))
      pick <- which(cn %in% c("p.value", "pvalue", "p", "pr(>chi)", "pr(prob)") | grepl("p", cn))
      if (length(pick) == 0) {
        suppressWarnings(as.numeric(hw_mat[1, ncol(hw_mat)]))
      } else {
        suppressWarnings(as.numeric(hw_mat[1, pick[1]]))
      }
    } else {
      suppressWarnings(as.numeric(hw_fit[1]))
    }
  }, error = function(e) {
    error_msg <<- conditionMessage(e)
    NA_real_
  })
  
  if (!is.finite(pval) || is.na(pval)) {
    base_row$status <- "failed"
    base_row$note <- if (!is.null(error_msg) && nzchar(error_msg)) {
      paste0("hw_test_failed: ", error_msg)
    } else {
      "hw_test_failed"
    }
    return(base_row)
  }
  
  base_row$p_value <- pval
  base_row$status <- "ok"
  base_row$note <- "tested"
  base_row
}

# ------------------------------------------------------------
# Build Site x Locus table from clone-corrected genotypes
# ------------------------------------------------------------
gdf <- adegenet::genind2df(gi_mll, sep = "/")
loci <- setdiff(names(gdf), "pop")

if (length(loci) == 0) {
  stop("[02_hwe] No loci detected in gi_mll.")
}

sites <- sort(unique(site_vec[!is.na(site_vec) & nzchar(site_vec)]))
if (length(sites) == 0) {
  stop("[02_hwe] No valid Site labels found in gi_mll.")
}

site_sizes <- data.frame(
  Site = names(table(site_vec)),
  N_clone_corrected = as.integer(table(site_vec)),
  stringsAsFactors = FALSE
) %>%
  arrange(Site)

message("[02_hwe] Sites detected: ", length(sites))
message("[02_hwe] Loci detected: ", length(loci))

rows <- vector("list", length(sites) * length(loci))
k <- 1L

for (s in sites) {
  idx_site <- which(site_vec == s)
  
  for (loc in loci) {
    geno_site <- gdf[idx_site, loc]
    raw_non_missing <- sum(!is.na(geno_site) & trimws(as.character(geno_site)) != "")
    
    if (raw_non_missing == 0) {
      rows[[k]] <- data.frame(
        Site = s,
        Locus = loc,
        N = 0L,
        missing_fraction = 1,
        n_alleles = NA_integer_,
        n_unique_genotypes = NA_integer_,
        p_value = NA_real_,
        p_value_raw = NA_real_,
        significant_raw = FALSE,
        p_value_bonferroni = NA_real_,
        p_value_fdr = NA_real_,
        significant_bonferroni = FALSE,
        significant_fdr = FALSE,
        p_adjust_bh = NA_real_,
        significant_bh = FALSE,
        status = "skipped",
        note = "all_missing_locus_within_site",
        stringsAsFactors = FALSE
      )
      k <- k + 1L
      next
    }
    
    test <- run_single_locus_hwe(geno_site)
    
    rows[[k]] <- data.frame(
      Site = s,
      Locus = loc,
      N = as.integer(test$N),
      missing_fraction = as.numeric(test$missing_fraction),
      n_alleles = as.integer(test$n_alleles),
      n_unique_genotypes = as.integer(test$n_unique_genotypes),
      p_value = as.numeric(test$p_value),
      p_value_raw = as.numeric(test$p_value),
      significant_raw = isTRUE(!is.na(test$p_value) && test$p_value < HWE_ALPHA),
      p_value_bonferroni = NA_real_,
      p_value_fdr = NA_real_,
      significant_bonferroni = FALSE,
      significant_fdr = FALSE,
      p_adjust_bh = NA_real_,
      significant_bh = FALSE,
      status = as.character(test$status),
      note = as.character(test$note),
      stringsAsFactors = FALSE
    )
    k <- k + 1L
  }
}

hwe_site_locus <- dplyr::bind_rows(rows) %>%
  add_hwe_multipletest_columns(
    p_col = "p_value_raw",
    alpha = HWE_ALPHA,
    include_legacy_aliases = TRUE
  ) %>%
  arrange(Site, Locus)

if (all(hwe_site_locus$status != "ok")) {
  warning("[02_hwe] No Site x Locus HWE tests could be evaluated. Review sample sizes, missingness, and locus variability.")
}

# Main required output
main_required <- file.path(TABLES_DIR, "hwe_by_site_by_locus.csv")
write.csv(hwe_site_locus, main_required, row.names = FALSE)
message("[02_hwe] Saved: ", main_required)

# Also write a compact provenance table so methods choices remain explicit.
hwe_methods <- data.frame(
  analysis = "Hardy_Weinberg_test",
  genotype_object = "gi_mll",
  genotype_object_class = paste(class(gi_mll), collapse = ";"),
  clone_correction = "MLL_clone_corrected",
  grouping_variable = "Site",
  genind_handler_package = "adegenet",
  hwe_test_package = "pegas",
  hwe_test_function = "pegas::hw.test",
  loci_conversion_function = "pegas::as.loci",
  monte_carlo_replicates = HWE_MONTE_CARLO_REPS,
  multiple_testing_alpha = HWE_ALPHA,
  by_locus_correction_scope = "across_loci",
  by_site_correction_scope = "across_sites",
  by_site_by_locus_correction_scope = "across_all_valid_site_locus_tests",
  multiple_testing_methods = "Bonferroni and Benjamini-Hochberg_FDR",
  stringsAsFactors = FALSE
)

methods_file <- file.path(TABLES_DIR, "hwe_methods_metadata.csv")
write.csv(hwe_methods, methods_file, row.names = FALSE)
message("[02_hwe] Saved: ", methods_file)

# ------------------------------------------------------------
# Compatibility outputs for downstream scripts / supplements
# ------------------------------------------------------------
combine_fisher <- function(pv) {
  pv <- pv[is.finite(pv) & !is.na(pv) & pv > 0]
  if (length(pv) == 0) {
    return(NA_real_)
  }
  stat <- -2 * sum(log(pv))
  pchisq(stat, df = 2 * length(pv), lower.tail = FALSE)
}

hwe_by_locus <- hwe_site_locus %>%
  group_by(Locus) %>%
  summarise(
    n_sites_tested = sum(status == "ok"),
    n_sites_significant_raw = sum(significant_raw, na.rm = TRUE),
    n_sites_significant_bh = sum(significant_bh, na.rm = TRUE),
    mean_clone_corrected_N = mean(N[status == "ok"], na.rm = TRUE),
    mean_missing_fraction = mean(missing_fraction, na.rm = TRUE),
    p_value = combine_fisher(p_value_raw[status == "ok"]),
    .groups = "drop"
  ) %>%
  mutate(p_value_raw = p_value) %>%
  add_hwe_multipletest_columns(
    p_col = "p_value_raw",
    alpha = HWE_ALPHA,
    include_legacy_aliases = TRUE
  )

hwe_by_site <- hwe_site_locus %>%
  group_by(Site) %>%
  summarise(
    loci_tested = sum(status == "ok"),
    loci_significant_raw = sum(significant_raw, na.rm = TRUE),
    loci_significant_bh = sum(significant_bh, na.rm = TRUE),
    mean_missing_fraction = mean(missing_fraction, na.rm = TRUE),
    p_value = combine_fisher(p_value_raw[status == "ok"]),
    .groups = "drop"
  ) %>%
  mutate(p_value_raw = p_value) %>%
  left_join(site_sizes, by = "Site") %>%
  add_hwe_multipletest_columns(
    p_col = "p_value_raw",
    alpha = HWE_ALPHA,
    include_legacy_aliases = TRUE
  )

hwe_by_locus_within_site <- hwe_site_locus

save_main_and_supp <- function(df, filename) {
  main_file <- file.path(TABLES_DIR, filename)
  supp_file <- file.path(TABLES_SUPP_DIR, filename)
  write.csv(df, main_file, row.names = FALSE)
  write.csv(df, supp_file, row.names = FALSE)
  message("[02_hwe] Saved: ", main_file)
  message("[02_hwe] Saved: ", supp_file)
}

save_main_and_supp(hwe_by_locus, "hwe_by_locus.csv")
save_main_and_supp(hwe_by_site, "hwe_by_site.csv")
save_main_and_supp(hwe_by_locus_within_site, "hwe_by_locus_within_site.csv")

# ------------------------------------------------------------
# NEW SECTION: Multilocus linkage disequilibrium (LD) testing
# ------------------------------------------------------------
# Rationale and interpretation guardrails:
# - This LD section uses the same clone-corrected object (gi_mll) as HWE,
#   so repeated clonemates do not inflate genotypic associations.
# - For diploid, codominant microsatellite loci, we test non-random
#   association among multilocus genotypes between pairs of loci using
#   contingency-table tests (chi-squared with Monte Carlo p-values).
# - Elevated LD across ostensibly unlinked loci can be consistent with a
#   pooled-structure (Wahlund-effect-type) pattern, but LD alone is not
#   definitive proof of a Wahlund effect.
# - LD can also arise from physical linkage, drift, small sample size,
#   admixture history, or other demographic processes.
# - Therefore LD results should be treated as supportive evidence and
#   interpreted jointly with HWE/FIS, STRUCTURE/PCA/DAPC, and FST/Jost's D.

LD_MONTE_CARLO_REPS <- 9999L
LD_MIN_NON_MISSING_N <- 8L
LD_ALPHA <- 0.05

message("[02_hwe] Running multilocus LD tests on gi_mll (clone-corrected)...")
message("[02_hwe] Pairwise LD test: stats::chisq.test (Monte Carlo p-values), B = ", LD_MONTE_CARLO_REPS)
message("[02_hwe] Multiple-testing correction for pairwise LD: Benjamini-Hochberg (FDR)")

# ------------------------------------------------------------
# Helper: pairwise locus-by-locus LD within a given genind scope
# ------------------------------------------------------------
run_pairwise_ld_for_scope <- function(gobj,
                                      scope_name,
                                      min_n = LD_MIN_NON_MISSING_N,
                                      B = LD_MONTE_CARLO_REPS,
                                      alpha = LD_ALPHA) {
  gdf_scope <- adegenet::genind2df(gobj, sep = "/")
  loci_scope <- setdiff(names(gdf_scope), "pop")
  
  if (length(loci_scope) < 2) {
    return(data.frame(
      Scope = character(0),
      Locus1 = character(0),
      Locus2 = character(0),
      Test_statistic = numeric(0),
      P_value = numeric(0),
      P_value_adjusted = numeric(0),
      Significant = logical(0),
      N_non_missing_pair = integer(0),
      Missing_fraction_pair = numeric(0),
      Status = character(0),
      Notes = character(0),
      stringsAsFactors = FALSE
    ))
  }
  
  pairs <- utils::combn(loci_scope, 2, simplify = FALSE)
  out <- vector("list", length(pairs))
  
  for (i in seq_along(pairs)) {
    l1 <- pairs[[i]][1]
    l2 <- pairs[[i]][2]
    g1 <- trimws(as.character(gdf_scope[[l1]]))
    g2 <- trimws(as.character(gdf_scope[[l2]]))
    
    keep <- !(is_missing_genotype(g1) | is_missing_genotype(g2))
    n_non_missing <- sum(keep)
    miss_frac <- 1 - (n_non_missing / nrow(gdf_scope))
    
    if (n_non_missing < min_n) {
      out[[i]] <- data.frame(
        Scope = scope_name,
        Locus1 = l1,
        Locus2 = l2,
        Test_statistic = NA_real_,
        P_value = NA_real_,
        P_value_adjusted = NA_real_,
        Significant = FALSE,
        N_non_missing_pair = as.integer(n_non_missing),
        Missing_fraction_pair = as.numeric(miss_frac),
        Status = "skipped",
        Notes = sprintf("too_few_individuals_non_missing_for_pair (N=%d < %d)", n_non_missing, min_n),
        stringsAsFactors = FALSE
      )
      next
    }
    
    tbl <- table(g1[keep], g2[keep])
    
    if (nrow(tbl) < 2 || ncol(tbl) < 2) {
      out[[i]] <- data.frame(
        Scope = scope_name,
        Locus1 = l1,
        Locus2 = l2,
        Test_statistic = NA_real_,
        P_value = NA_real_,
        P_value_adjusted = NA_real_,
        Significant = FALSE,
        N_non_missing_pair = as.integer(n_non_missing),
        Missing_fraction_pair = as.numeric(miss_frac),
        Status = "skipped",
        Notes = "invariant_or_single_level_locus_pair",
        stringsAsFactors = FALSE
      )
      next
    }
    
    err <- NULL
    ld_fit <- tryCatch({
      stats::chisq.test(tbl, simulate.p.value = TRUE, B = B)
    }, error = function(e) {
      err <<- conditionMessage(e)
      NULL
    })
    
    if (is.null(ld_fit) || is.null(ld_fit$p.value) || !is.finite(ld_fit$p.value)) {
      out[[i]] <- data.frame(
        Scope = scope_name,
        Locus1 = l1,
        Locus2 = l2,
        Test_statistic = NA_real_,
        P_value = NA_real_,
        P_value_adjusted = NA_real_,
        Significant = FALSE,
        N_non_missing_pair = as.integer(n_non_missing),
        Missing_fraction_pair = as.numeric(miss_frac),
        Status = "failed",
        Notes = if (!is.null(err) && nzchar(err)) paste0("ld_test_failed: ", err) else "ld_test_failed",
        stringsAsFactors = FALSE
      )
      next
    }
    
    out[[i]] <- data.frame(
      Scope = scope_name,
      Locus1 = l1,
      Locus2 = l2,
      Test_statistic = as.numeric(ld_fit$statistic),
      P_value = as.numeric(ld_fit$p.value),
      P_value_adjusted = NA_real_,
      Significant = FALSE,
      N_non_missing_pair = as.integer(n_non_missing),
      Missing_fraction_pair = as.numeric(miss_frac),
      Status = "ok",
      Notes = "tested",
      stringsAsFactors = FALSE
    )
  }
  
  res <- dplyr::bind_rows(out)
  ok <- which(res$Status == "ok" & is.finite(res$P_value) & !is.na(res$P_value))
  if (length(ok) > 0) {
    res$P_value_adjusted[ok] <- p.adjust(res$P_value[ok], method = "BH")
    res$Significant[ok] <- res$P_value_adjusted[ok] < alpha
  }
  
  res %>%
    arrange(Scope, Locus1, Locus2)
}

# ------------------------------------------------------------
# Overall pooled and by-site pairwise LD tests
# ------------------------------------------------------------
ld_pairwise_overall <- run_pairwise_ld_for_scope(
  gobj = gi_mll,
  scope_name = "overall"
)

ld_pairwise_by_site <- dplyr::bind_rows(lapply(sites, function(s) {
  idx <- which(site_vec == s)
  if (length(idx) == 0) {
    return(NULL)
  }
  gi_site <- gi_mll[idx, ]
  run_pairwise_ld_for_scope(
    gobj = gi_site,
    scope_name = s
  )
}))

ld_pairwise_all <- dplyr::bind_rows(ld_pairwise_overall, ld_pairwise_by_site) %>%
  select(
    Scope,
    Locus1,
    Locus2,
    Test_statistic,
    P_value,
    P_value_adjusted,
    Significant,
    N_non_missing_pair,
    Missing_fraction_pair,
    Status,
    Notes
  )

# Required pairwise output
ld_pairwise_file <- file.path(TABLES_DIR, "linkage_disequilibrium_pairwise.csv")
write.csv(ld_pairwise_all, ld_pairwise_file, row.names = FALSE)
message("[02_hwe] Saved: ", ld_pairwise_file)

# Optional by-site convenience output (requested when feasible)
ld_by_site_summary <- ld_pairwise_all %>%
  filter(Scope != "overall") %>%
  group_by(Scope) %>%
  summarise(
    Number_of_loci = dplyr::n_distinct(c(Locus1, Locus2)),
    Number_of_pairwise_tests = sum(Status == "ok", na.rm = TRUE),
    Number_significant_before_correction = sum(Status == "ok" & !is.na(P_value) & P_value < LD_ALPHA, na.rm = TRUE),
    Number_significant_after_correction = sum(Status == "ok" & Significant, na.rm = TRUE),
    Notes = "Pairwise multilocus genotype-association tests by Site (Monte Carlo chi-squared, BH-FDR within Site).",
    .groups = "drop"
  ) %>%
  rename(Site = Scope) %>%
  arrange(Site)

ld_by_site_file <- file.path(TABLES_DIR, "linkage_disequilibrium_by_site.csv")
write.csv(ld_by_site_summary, ld_by_site_file, row.names = FALSE)
message("[02_hwe] Saved: ", ld_by_site_file)

# ------------------------------------------------------------
# Global LD summary table (overall + each Site)
# ------------------------------------------------------------
build_ld_global_row <- function(df_scope, scope_label) {
  if (nrow(df_scope) == 0) {
    return(data.frame(
      Scope = scope_label,
      Number_of_loci = 0L,
      Number_of_pairwise_tests = 0L,
      Number_significant_before_correction = 0L,
      Number_significant_after_correction = 0L,
      Notes = "No locus pairs available for LD testing in this scope.",
      stringsAsFactors = FALSE
    ))
  }
  
  tested <- df_scope$Status == "ok"
  data.frame(
    Scope = scope_label,
    Number_of_loci = as.integer(dplyr::n_distinct(c(df_scope$Locus1, df_scope$Locus2))),
    Number_of_pairwise_tests = as.integer(sum(tested, na.rm = TRUE)),
    Number_significant_before_correction = as.integer(sum(tested & !is.na(df_scope$P_value) & df_scope$P_value < LD_ALPHA, na.rm = TRUE)),
    Number_significant_after_correction = as.integer(sum(tested & !is.na(df_scope$P_value_adjusted) & df_scope$P_value_adjusted < LD_ALPHA, na.rm = TRUE)),
    Notes = "Elevated LD among unlinked loci can be consistent with pooled structure (Wahlund effect) but is supportive only; interpret with HWE/FIS, structure analyses, and differentiation metrics.",
    stringsAsFactors = FALSE
  )
}

ld_global <- dplyr::bind_rows(
  build_ld_global_row(ld_pairwise_all %>% filter(Scope == "overall"), "overall"),
  dplyr::bind_rows(lapply(sort(unique(ld_pairwise_all$Scope[ld_pairwise_all$Scope != "overall"])), function(s) {
    build_ld_global_row(ld_pairwise_all %>% filter(Scope == s), s)
  }))
)

ld_global_file <- file.path(TABLES_DIR, "linkage_disequilibrium_global.csv")
write.csv(ld_global, ld_global_file, row.names = FALSE)
message("[02_hwe] Saved: ", ld_global_file)

# ------------------------------------------------------------
# Optional Excel exports (only if writexl is installed)
# ------------------------------------------------------------
if (requireNamespace("writexl", quietly = TRUE)) {
  writexl::write_xlsx(ld_pairwise_all, path = file.path(TABLES_DIR, "linkage_disequilibrium_pairwise.xlsx"))
  writexl::write_xlsx(ld_global, path = file.path(TABLES_DIR, "linkage_disequilibrium_global.xlsx"))
  writexl::write_xlsx(ld_by_site_summary, path = file.path(TABLES_DIR, "linkage_disequilibrium_by_site.xlsx"))
  message("[02_hwe] Saved Excel LD outputs to: ", TABLES_DIR)
} else {
  message("[02_hwe] Package 'writexl' not installed; skipping Excel LD exports.")
}