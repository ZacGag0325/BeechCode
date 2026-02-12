############################################################
# scripts/02_hwe.R
# Hardy-Weinberg tests by Site + global (MLG vs MLL)
############################################################

options(repos = c(CRAN = "https://cloud.r-project.org"))
pkgs <- c("adegenet", "dplyr", "pegas")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)

suppressPackageStartupMessages({
  library(adegenet)
  library(dplyr)
  library(pegas)
})

# Make the script runnable from RStudio/Terminal regardless of working directory
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

setwd(find_project_root())
source(file.path("scripts", "_load_objects.R"))

OUTDIR <- file.path(RUN_OUT, "hwe_only")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

HWE_ALPHA <- 0.05
HWE_B <- 10000
MIN_N <- 8

safe_genind2loci_df <- function(gpop, sep = "/") {
  gdf <- adegenet::genind2df(gpop, sep = sep)
  gdf <- as.data.frame(lapply(gdf, function(col) {
    x <- trimws(as.character(col))
    x[x %in% c("", "NA", "NaN", "0", "0/0", "NA/NA")] <- NA
    ok <- is.na(x) | grepl("^[^/]+/[^/]+$", x)
    x[!ok] <- NA
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
    col_pick <- which(
      cn %in% c("p.value", "pvalue", "p", "pval", "pvals") |
        grepl("^pr", cn) |
        grepl("p\\s*value", cn)
    )
    if (length(col_pick) == 0) col_pick <- ncol(pmat)
    pvals <- suppressWarnings(as.numeric(pmat[, col_pick[1]]))
    locus <- rownames(pmat)
    return(list(pvals = pvals, locus = locus))
  }
  
  pv <- suppressWarnings(as.numeric(ht))
  list(pvals = pv, locus = names(ht))
}

run_hwe_for_genind <- function(genind_obj, label = "POP", min_n = MIN_N) {
  n_s <- nInd(genind_obj)
  if (n_s < min_n) {
    return(tibble(
      Site = label,
      Locus = locNames(genind_obj),
      p_value = NA_real_,
      n_inds = n_s,
      reason = paste0("too_few_individuals (n<", min_n, ")")
    ))
  }
  
  gdf <- safe_genind2loci_df(genind_obj, sep = "/")
  
  gdf_chr <- unlist(lapply(gdf, as.character), use.names = FALSE)
  allele_counts_seen <- sapply(gdf_chr, function(z) {
    if (is.na(z) || z == "") return(NA_integer_)
    length(strsplit(as.character(z), "/", fixed = TRUE)[[1]])
  })
  
  all_diploid_for_pegas <- all(is.na(allele_counts_seen) | allele_counts_seen == 2)
  B_use <- if (all_diploid_for_pegas) HWE_B else 0
  
  loci_obj <- pegas::as.loci(gdf)
  ht <- tryCatch(
    pegas::hw.test(loci_obj, B = B_use),
    error = function(e) e
  )
  
  if (inherits(ht, "error")) {
    return(tibble(
      Site = label,
      Locus = locNames(genind_obj),
      p_value = NA_real_,
      n_inds = n_s,
      reason = paste("hw.test error:", ht$message)
    ))
  }
  
  ex <- extract_hw_pvals(ht)
  pvals <- ex$pvals
  locus_names <- ex$locus
  
  if (length(pvals) == 0) {
    return(tibble(
      Site = label,
      Locus = locNames(genind_obj),
      p_value = NA_real_,
      n_inds = n_s,
      reason = "no_pvalues_returned"
    ))
  }
  
  if (is.null(locus_names) || length(locus_names) != length(pvals)) {
    locus_names <- paste0("Locus_", seq_along(pvals))
  }
  
  expected_loci <- locNames(genind_obj)
  ptab <- tibble(Locus = locus_names, p_value = as.numeric(pvals))
  
  tibble(
    Site = label,
    Locus = expected_loci,
    n_inds = n_s
  ) %>%
    left_join(ptab, by = "Locus") %>%
    mutate(
      reason = case_when(
        is.na(p_value) ~ "p_value_missing_from_test",
        TRUE ~ NA_character_
      )
    )
}

sites <- levels(pop(gi))
cat("Sites in pop(gi):\n")
print(sites)

hwe_by_site <- lapply(sites, function(s) {
  gi_s <- gi[pop(gi) == s, , drop = FALSE]
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
    sig_p05 = !is.na(p_value) & p_value < HWE_ALPHA,
    sig_bh05 = !is.na(p_adj_bh) & p_adj_bh < HWE_ALPHA
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
    n_NA = sum(is.na(p_value)),
    .groups = "drop"
  )

run_hwe_global <- function(genind_obj, label) {
  out <- run_hwe_for_genind(genind_obj, label = label)
  out %>%
    mutate(
      p_adj_bh = p.adjust(p_value, method = "BH"),
      sig_p05 = !is.na(p_value) & p_value < HWE_ALPHA,
      sig_bh05 = !is.na(p_adj_bh) & p_adj_bh < HWE_ALPHA
    )
}

hwe_global_mlg <- run_hwe_global(gi, "GLOBAL_MLG")
hwe_global_mll <- run_hwe_global(gi_mll, "GLOBAL_MLL")

hwe_compare <- full_join(
  hwe_global_mlg %>% select(Locus, p_MLG = p_value),
  hwe_global_mll %>% select(Locus, p_MLL = p_value),
  by = "Locus"
) %>%
  mutate(
    sig_MLG = !is.na(p_MLG) & p_MLG < HWE_ALPHA,
    sig_MLL = !is.na(p_MLL) & p_MLL < HWE_ALPHA
  )

cat("\n--- HWE NA reasons (by site) ---\n")
print(hwe_by_site %>% filter(is.na(p_value)) %>% count(reason, sort = TRUE), n = 50)

cat("\n--- Site summary ---\n")
print(hwe_site_summary, n = Inf)

write.csv(hwe_by_site, file.path(OUTDIR, "hwe_by_site_by_locus.csv"), row.names = FALSE)
write.csv(hwe_site_summary, file.path(OUTDIR, "hwe_site_summary.csv"), row.names = FALSE)
write.csv(hwe_global_mlg, file.path(OUTDIR, "hwe_global_mlg.csv"), row.names = FALSE)
write.csv(hwe_global_mll, file.path(OUTDIR, "hwe_global_mll.csv"), row.names = FALSE)
write.csv(hwe_compare, file.path(OUTDIR, "hwe_global_compare_mlg_vs_mll.csv"), row.names = FALSE)

cat("DONE HWE. Outputs in: ", OUTDIR, "\n", sep = "")
