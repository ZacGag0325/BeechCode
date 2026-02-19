# filename: scripts/02_hwe.R
############################################################
# scripts/02_hwe.R
# Hardy-Weinberg tests by Site + global (clone-corrected)
############################################################

options(repos = c(CRAN = "https://cloud.r-project.org"))
pkgs <- c("adegenet", "dplyr", "pegas")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)

suppressPackageStartupMessages({
  library(adegenet)
  library(dplyr)
  library(pegas)
})

find_project_root <- function() {
  candidates <- c(getwd(), normalizePath(file.path(getwd(), ".."), mustWork = FALSE))
  cmd_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
  if (length(cmd_file) > 0 && nzchar(cmd_file[1])) {
    candidates <- c(candidates, dirname(normalizePath(cmd_file[1], mustWork = FALSE)))
  }
  for (start in unique(candidates)) {
    cur <- normalizePath(start, mustWork = FALSE)
    repeat {
      if (file.exists(file.path(cur, "scripts", "_load_objects.R"))) return(cur)
      parent <- dirname(cur)
      if (identical(parent, cur)) break
      cur <- parent
    }
  }
  stop("Cannot find project root containing scripts/_load_objects.R. Open BeechCode project first.")
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

setwd(find_project_root())
source(file.path("scripts", "_load_objects.R"))

OUTDIR <- file.path(RUN_OUT, "hwe_only")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# Use clone-corrected dataset for independence-assuming analyses
gi_use <- gi_mll


HWE_ALPHA <- 0.05
HWE_B <- 10000
MIN_N <- 8

safe_genind2loci_df <- function(gpop, sep = "/") {
  gdf <- adegenet::genind2df(gpop, sep = sep)
  gdf <- as.data.frame(lapply(gdf, function(col) {
    x <- trimws(as.character(col))
    x[x %in% c("", "NA", "NaN", "0", "0/0", "NA/NA", "-")] <- NA
    ok <- is.na(x) | grepl("^[^/]+/[^/]+$", x)
    x[!ok] <- NA
    x <- ifelse(is.na(x), NA_character_, vapply(strsplit(x, "/", fixed = TRUE), function(z) {
      paste(sort(z), collapse = "/")
    }, character(1)))
    factor(x)
  }), check.names = FALSE, stringsAsFactors = FALSE)
  gdf
}

extract_hw_pvals <- function(ht) {
  if (is.list(ht) && !is.null(ht$p.value)) {
    pv <- ht$p.value
    return(list(pvals = as.numeric(pv), locus = names(pv)))
  }
  if (is.matrix(ht) || is.data.frame(ht)) {
    pmat <- as.matrix(ht)
    cn <- tolower(colnames(pmat))
    col_pick <- which(cn %in% c("p.value", "pvalue", "p", "pval", "pvals") |
                        grepl("^pr", cn) |
                        grepl("p\\s*value", cn))
    if (length(col_pick) == 0) col_pick <- ncol(pmat)
    pvals <- suppressWarnings(as.numeric(pmat[, col_pick[1]]))
    locus <- rownames(pmat)
    return(list(pvals = pvals, locus = locus))
  }
  pv <- suppressWarnings(as.numeric(ht))
  list(pvals = pv, locus = names(ht))
}

all_diploid_genotypes <- function(df_loci) {
  vals <- unlist(df_loci, use.names = FALSE)
  vals <- as.character(vals)
  vals <- vals[!is.na(vals) & nzchar(vals)]
  if (length(vals) == 0) return(FALSE)
  all(grepl("^[^/]+/[^/]+$", vals))
}

locus_reason <- function(x, min_n = MIN_N) {
  xv <- as.character(x)
  xv <- xv[!is.na(xv) & nzchar(xv)]
  if (length(xv) < min_n) return(paste0("too_few_non_missing_genotypes (n<", min_n, ")"))
  if (any(!grepl("^[^/]+/[^/]+$", xv))) return("invalid_genotype_format")
  alleles <- unique(unlist(strsplit(xv, "/", fixed = TRUE), use.names = FALSE))
  if (length(alleles) <= 1) return("monomorphic_locus")
  NA_character_
}

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
  
  base_tbl <- tibble(
    Site = label,
    Locus = loci,
    n_inds = n_s,
    n_non_missing = vapply(gdf[loci], function(x) sum(!is.na(x)), integer(1)),
    reason = vapply(gdf[loci], locus_reason, character(1), min_n = min_n)
  )
  
  testable <- base_tbl %>% filter(is.na(reason))
  if (nrow(testable) == 0) {
    return(base_tbl %>% mutate(p_value = NA_real_) %>% select(Site, Locus, p_value, n_inds, n_non_missing, reason))
  }
  
  loci_obj <- tryCatch(pegas::as.loci(gdf[, testable$Locus, drop = FALSE], col.pop = NULL), error = function(e) e)
  if (inherits(loci_obj, "error")) {
    return(base_tbl %>% mutate(
      p_value = NA_real_,
      reason = ifelse(is.na(reason), paste("as.loci_conversion_error:", loci_obj$message), reason)
    ) %>% select(Site, Locus, p_value, n_inds, n_non_missing, reason))
  }
  
  diploid_ok <- all_diploid_genotypes(gdf[, testable$Locus, drop = FALSE])
  B_use <- if (isTRUE(diploid_ok)) HWE_B else 0
  if (!isTRUE(diploid_ok)) {
    cat("[HWE] Note:", label, "contains non-diploid/ambiguous genotypes for tested loci; using exact test without Monte Carlo (B=0).\n")
  }
  
  ht <- tryCatch(pegas::hw.test(loci_obj, B = B_use), error = function(e) e)
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

sites <- levels(pop(gi_use))
cat("Sites in pop(gi_mll):\n")
print(sites)

hwe_by_site <- lapply(sites, function(s) {
  gi_s <- gi_use[pop(gi_use) == s, , drop = FALSE]
  out <- run_hwe_for_genind(gi_s, label = s)
  cat(
    "[HWE] Site", s,
    ": nInd =", nInd(gi_s),
    "| loci tested =", sum(!is.na(out$p_value)),
    "| n significant (p<", HWE_ALPHA, ") =",
    sum(out$p_value < HWE_ALPHA, na.rm = TRUE), "\n"
  )
  out
}) %>% bind_rows()

hwe_by_site <- hwe_by_site %>%
  group_by(Site) %>%
  mutate(
    p_adj_bh = p.adjust(p_value, method = "BH"),
    p_adj_bonf = p.adjust(p_value, method = "bonferroni"),
    sig_p05 = !is.na(p_value) & p_value < HWE_ALPHA,
    sig_bh05 = !is.na(p_adj_bh) & p_adj_bh < HWE_ALPHA,
    sig_bonf05 = !is.na(p_adj_bonf) & p_adj_bonf < HWE_ALPHA
  ) %>%
  ungroup()

hwe_site_summary <- hwe_by_site %>%
  group_by(Site) %>%
  summarise(
    n_inds = dplyr::first(n_inds),
    loci_total = n(),
    loci_tested = sum(!is.na(p_value)),
    n_sig_p05 = sum(sig_p05, na.rm = TRUE),
    n_sig_bh05 = sum(sig_bh05, na.rm = TRUE),
    n_sig_bonf05 = sum(sig_bonf05, na.rm = TRUE),
    n_NA = sum(is.na(p_value)),
    .groups = "drop"
  )


add_hwe_summary_cols <- function(tbl, alpha = HWE_ALPHA) {
  n_sig_p05 <- sum(!is.na(tbl$p_value) & tbl$p_value < alpha, na.rm = TRUE)
  n_sig_bh05 <- sum(!is.na(tbl$p_adj_bh) & tbl$p_adj_bh < alpha, na.rm = TRUE)
  n_sig_bonf05 <- sum(!is.na(tbl$p_adj_bonf) & tbl$p_adj_bonf < alpha, na.rm = TRUE)
  tbl %>%
    mutate(
      n_sig_p05 = n_sig_p05,
      n_sig_bh05 = n_sig_bh05,
      n_sig_bonf05 = n_sig_bonf05
    )
}

run_hwe_global <- function(genind_obj, label) {
  out <- run_hwe_for_genind(genind_obj, label = label)
  out %>%
    mutate(
      p_adj_bh = p.adjust(p_value, method = "BH"),
      p_adj_bonf = p.adjust(p_value, method = "bonferroni"),
      sig_p05 = !is.na(p_value) & p_value < HWE_ALPHA,
      sig_bh05 = !is.na(p_adj_bh) & p_adj_bh < HWE_ALPHA,
      sig_bonf05 = !is.na(p_adj_bonf) & p_adj_bonf < HWE_ALPHA
    ) %>%
    add_hwe_summary_cols(alpha = HWE_ALPHA)
}

hwe_global <- run_hwe_global(gi_use, "GLOBAL_MLL")

cat("\n--- HWE NA reasons (by site) ---\n")
print(hwe_by_site %>% filter(is.na(p_value)) %>% count(reason, sort = TRUE), n = 50)

cat("\n--- Site summary ---\n")
print(hwe_site_summary, n = Inf)

write.csv(hwe_by_site, file.path(OUTDIR, "hwe_by_site_by_locus.csv"), row.names = FALSE)
write.csv(hwe_site_summary, file.path(OUTDIR, "hwe_site_summary.csv"), row.names = FALSE)
write.csv(hwe_global, file.path(OUTDIR, "hwe_global_mll.csv"), row.names = FALSE)

cat("DONE HWE. Outputs in: ", OUTDIR, "\n", sep = "")

