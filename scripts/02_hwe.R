# --- PATCH for scripts/02_hwe.R ---
# Replace your current run_hwe_for_genind() with this exact function.

run_hwe_for_genind <- function(genind_obj, label = "POP", min_n = MIN_N) {
  loci <- unique(adegenet::locNames(genind_obj))
  n_s <- adegenet::nInd(genind_obj)
  if (n_s < min_n) {
    return(tibble(
      Site = label,
      Locus = loci,
      p_value = NA_real_,
      n_inds = n_s,
      n_non_missing = 0L,
      reason = paste0("too_few_individuals (n<", min_n, ")")
    ))
  }
  
  gdf <- safe_genind2loci_df(genind_obj, sep = "/")
  
  missing_loci <- setdiff(loci, colnames(gdf))
  if (length(missing_loci) > 0) {
    stop("HWE input error: missing loci in genotype table for ", label, ": ", paste(missing_loci, collapse = ", "))
  }
  gdf <- gdf[, loci, drop = FALSE]
  
  n_non_missing <- vapply(loci, function(loc) sum(!is.na(gdf[[loc]])), integer(1))
  reason_vec <- vapply(loci, function(loc) locus_reason(gdf[[loc]], min_n = min_n), character(1))
  
  base_tbl <- tibble(
    Site = label,
    Locus = loci,
    n_inds = n_s,
    n_non_missing = n_non_missing,
    reason = reason_vec
  )
  
  testable <- base_tbl %>% filter(is.na(reason))
  if (nrow(testable) == 0) {
    return(base_tbl %>% mutate(p_value = NA_real_) %>% select(Site, Locus, p_value, n_inds, n_non_missing, reason))
  }
  
  loci_obj <- tryCatch(pegas::as.loci(gdf[, testable$Locus, drop = FALSE]), error = function(e) e)
  if (inherits(loci_obj, "error")) {
    return(base_tbl %>% mutate(
      p_value = NA_real_,
      reason = ifelse(is.na(reason), paste("as.loci_conversion_error:", loci_obj$message), reason)
    ) %>% select(Site, Locus, p_value, n_inds, n_non_missing, reason))
  }
  
  ht <- tryCatch(pegas::hw.test(loci_obj, B = HWE_B), error = function(e) e)
  if (inherits(ht, "error")) {
    return(base_tbl %>% mutate(
      p_value = NA_real_,
      reason = ifelse(is.na(reason), paste("hw.test_error:", ht$message), reason)
    ) %>% select(Site, Locus, p_value, n_inds, n_non_missing, reason))
  }
  
  ex <- extract_hw_pvals(ht)
  ptab <- tibble(Locus = ex$locus %||% character(0), p_value = as.numeric(ex$pvals))
  if (nrow(ptab) == 0) {
    return(base_tbl %>% mutate(
      p_value = NA_real_,
      reason = ifelse(is.na(reason), "no_pvalues_returned", reason)
    ) %>% select(Site, Locus, p_value, n_inds, n_non_missing, reason))
  }
  
  out <- base_tbl %>%
    left_join(ptab, by = "Locus") %>%
    mutate(reason = ifelse(is.na(reason) & is.na(p_value), "p_value_missing_from_test", reason)) %>%
    select(Site, Locus, p_value, n_inds, n_non_missing, reason)
  
  out
}

