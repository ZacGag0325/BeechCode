# scripts/02_hwe.R
############################################################
# Hardy-Weinberg Equilibrium tests (clone-corrected: gi_mll)
#
# Design goals:
# - Use clone-corrected data for population-level inference.
# - Avoid pegas::loci (not exported in current pegas versions).
# - Run by Site x Locus with robust failure handling.
# - Continue even if one locus fails.
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

message("[02_hwe] Running HWE tests on gi_mll (clone-corrected)...")
validate_columns(df_ids_mll, c("ind_id", "Site"), df_name = "02_hwe df_ids_mll")
if (!all(adegenet::indNames(gi_mll) == df_ids_mll$ind_id)) stop("[02_hwe] gi_mll and df_ids_mll are not aligned.")

# ------------------------------------------------------------
# Helper: robust allele-pair parsing
# ------------------------------------------------------------
split_genotype <- function(x) {
  x <- trimws(as.character(x))
  if (is.na(x) || x == "" || x %in% c("NA", "0", "0/0", "NA/NA", "-")) return(c(NA_character_, NA_character_))
  parts <- strsplit(x, "/", fixed = TRUE)[[1]]
  if (length(parts) != 2) return(c(NA_character_, NA_character_))
  parts <- trimws(parts)
  if (any(parts == "")) return(c(NA_character_, NA_character_))
  sort(parts)
}

# ------------------------------------------------------------
# Helper: single-locus HWE test with defensive checks
# ------------------------------------------------------------
run_single_locus_hwe <- function(geno_vec, B = 1000, min_n = 8, min_unique_genotypes = 2) {
  parsed <- t(vapply(geno_vec, split_genotype, character(2)))
  a1 <- parsed[, 1]
  a2 <- parsed[, 2]
  
  keep <- !(is.na(a1) | is.na(a2))
  n_non_missing <- sum(keep)
  
  base_row <- list(
    N = as.integer(n_non_missing),
    p_value = NA_real_,
    status = "not_tested",
    note = ""
  )
  
  if (n_non_missing < min_n) {
    base_row$status <- "skipped"
    base_row$note <- sprintf("too_few_individuals_non_missing (N=%d < %d)", n_non_missing, min_n)
    return(base_row)
  }
  
  g <- paste(a1[keep], a2[keep], sep = "/")
  n_genotypes <- dplyr::n_distinct(g)
  if (n_genotypes < min_unique_genotypes) {
    base_row$status <- "skipped"
    base_row$note <- sprintf("insufficient_genotype_diversity (unique_genotypes=%d)", n_genotypes)
    return(base_row)
  }
  
  # Invariant locus check
  alleles <- c(a1[keep], a2[keep])
  n_alleles <- dplyr::n_distinct(alleles)
  if (n_alleles < 2) {
    base_row$status <- "skipped"
    base_row$note <- "invariant_locus"
    return(base_row)
  }
  
  # Build a single-locus data.frame and convert with pegas::as.loci (NOT pegas::loci)
  loc_df <- data.frame(Locus = factor(g), stringsAsFactors = TRUE)
  
  pval <- tryCatch({
    loci_obj <- pegas::as.loci(loc_df)
    ht <- pegas::hw.test(loci_obj, B = B)
    
    # Robust p-value extraction across pegas output shapes
    if (is.list(ht) && !is.null(ht$p.value)) {
      as.numeric(ht$p.value[1])
    } else if (is.matrix(ht) || is.data.frame(ht)) {
      htm <- as.matrix(ht)
      cn <- tolower(colnames(htm))
      pick <- which(cn %in% c("p.value", "pvalue", "p", "pr(>chi)", "pr(prob)") | grepl("p", cn))
      if (length(pick) == 0) {
        suppressWarnings(as.numeric(htm[1, ncol(htm)]))
      } else {
        suppressWarnings(as.numeric(htm[1, pick[1]]))
      }
    } else {
      suppressWarnings(as.numeric(ht[1]))
    }
  }, error = function(e) {
    attr(base_row, "error_msg") <<- conditionMessage(e)
    NA_real_
  })
  
  if (is.na(pval) || !is.finite(pval)) {
    base_row$status <- "failed"
    err <- attr(base_row, "error_msg")
    base_row$note <- if (!is.null(err) && nzchar(err)) paste0("hw_test_failed: ", err) else "hw_test_failed"
    return(base_row)
  }
  
  base_row$p_value <- pval
  base_row$status <- "ok"
  base_row$note <- "tested"
  base_row
}

# ------------------------------------------------------------
# Build Site x Locus table
# ------------------------------------------------------------
gdf <- adegenet::genind2df(gi_mll, sep = "/")
loci <- setdiff(names(gdf), "pop")
site_vec <- as.character(adegenet::pop(gi_mll))

if (length(loci) == 0) stop("[02_hwe] No loci detected in gi_mll.")
if (all(is.na(site_vec))) stop("[02_hwe] pop(gi_mll) is all NA; cannot run by-site HWE.")

sites <- sort(unique(site_vec))

rows <- vector("list", length(sites) * length(loci))
k <- 1L

for (s in sites) {
  idx_site <- which(site_vec == s)
  for (loc in loci) {
    geno_site <- gdf[idx_site, loc]
    
    # all-missing check before full parsing
    raw_non_missing <- sum(!is.na(geno_site) & trimws(as.character(geno_site)) != "")
    if (raw_non_missing == 0) {
      rows[[k]] <- data.frame(
        Site = s,
        Locus = loc,
        N = 0L,
        p_value = NA_real_,
        significant_raw = FALSE,
        p_adjust_bh = NA_real_,
        significant_bh = FALSE,
        status = "skipped",
        note = "all_missing_locus_within_site",
        stringsAsFactors = FALSE
      )
      k <- k + 1L
      next
    }
    
    test <- run_single_locus_hwe(geno_site, B = 1000, min_n = 8, min_unique_genotypes = 2)
    
    rows[[k]] <- data.frame(
      Site = s,
      Locus = loc,
      N = as.integer(test$N),
      p_value = as.numeric(test$p_value),
      significant_raw = isTRUE(!is.na(test$p_value) && test$p_value < 0.05),
      p_adjust_bh = NA_real_,
      significant_bh = FALSE,
      status = as.character(test$status),
      note = as.character(test$note),
      stringsAsFactors = FALSE
    )
    k <- k + 1L
  }
}

hwe_site_locus <- dplyr::bind_rows(rows)

# BH correction within each Site across loci actually tested (status == ok)
hwe_site_locus <- hwe_site_locus %>%
  group_by(Site) %>%
  mutate(
    p_adjust_bh = {
      pv <- p_value
      ok <- !is.na(pv)
      out <- rep(NA_real_, length(pv))
      if (sum(ok) > 0) out[ok] <- p.adjust(pv[ok], method = "BH")
      out
    },
    significant_bh = !is.na(p_adjust_bh) & p_adjust_bh < 0.05
  ) %>%
  ungroup() %>%
  arrange(Site, Locus)

# Main required output
main_required <- file.path(TABLES_DIR, "hwe_by_site_by_locus.csv")
write.csv(hwe_site_locus, main_required, row.names = FALSE)
message("[02_hwe] Saved: ", main_required)

# Compatibility outputs for downstream scripts / supplements
combine_fisher <- function(pv) {
  pv <- pv[is.finite(pv) & !is.na(pv) & pv > 0]
  if (length(pv) == 0) return(NA_real_)
  stat <- -2 * sum(log(pv))
  pchisq(stat, df = 2 * length(pv), lower.tail = FALSE)
}

site_sizes <- data.frame(
  Site = as.character(levels(adegenet::pop(gi_mll))),
  N = as.integer(table(adegenet::pop(gi_mll))[levels(adegenet::pop(gi_mll))]),
  stringsAsFactors = FALSE
)

hwe_by_locus <- hwe_site_locus %>%
  group_by(Locus) %>%
  summarise(
    n_sites_tested = sum(status == "ok"),
    n_sites_significant_raw = sum(significant_raw, na.rm = TRUE),
    n_sites_significant_bh = sum(significant_bh, na.rm = TRUE),
    p_value = combine_fisher(p_value[status == "ok"]),
    .groups = "drop"
  )

hwe_by_locus$p_adjust_bh <- p.adjust(hwe_by_locus$p_value, method = "BH")
hwe_by_locus$significant_raw <- !is.na(hwe_by_locus$p_value) & hwe_by_locus$p_value < 0.05
hwe_by_locus$significant_bh <- !is.na(hwe_by_locus$p_adjust_bh) & hwe_by_locus$p_adjust_bh < 0.05

hwe_by_site <- hwe_site_locus %>%
  group_by(Site) %>%
  summarise(
    loci_tested = sum(status == "ok"),
    loci_significant_raw = sum(significant_raw, na.rm = TRUE),
    loci_significant_bh = sum(significant_bh, na.rm = TRUE),
    p_value = combine_fisher(p_value[status == "ok"]),
    .groups = "drop"
  ) %>%
  left_join(site_sizes, by = "Site")

hwe_by_site$p_adjust_bh <- p.adjust(hwe_by_site$p_value, method = "BH")
hwe_by_site$significant_raw <- !is.na(hwe_by_site$p_value) & hwe_by_site$p_value < 0.05
hwe_by_site$significant_bh <- !is.na(hwe_by_site$p_adjust_bh) & hwe_by_site$p_adjust_bh < 0.05

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