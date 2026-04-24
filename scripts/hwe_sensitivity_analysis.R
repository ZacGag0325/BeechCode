# scripts/hwe_sensitivity_analysis.R
############################################################
# Sensitivity analysis: rerun core genetics analyses with suspect loci removed
#
# Purpose:
# - Keep the main full-data workflow intact.
# - Build a reduced clone-corrected object excluding suspect loci.
# - Re-run key summaries and differentiation analyses on reduced data.
# - Save outputs with explicit `_noSuspectLoci` naming for side-by-side comparison.
############################################################

suppressPackageStartupMessages({
  library(adegenet)
  library(poppr)
  library(ade4)
  library(hierfstat)
  library(mmod)
  library(pegas)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

source("scripts/_load_objects.R")

# ------------------------------------------------------------------
# User-editable sensitivity settings
# ------------------------------------------------------------------
suspect_null_loci <- c("EJV8T_A_0", "ERHBI_A_0", "FCM5", "FG5")
analysis_suffix <- "noSuspectLoci"
AMOVA_PERMUTATIONS <- 999
HWE_MONTE_CARLO_REPS <- 9999L

SENS_TABLES_DIR <- file.path(TABLES_DIR, analysis_suffix)
SENS_FIGURES_DIR <- file.path(FIGURES_DIR, analysis_suffix)
SENS_MATRICES_DIR <- file.path(MATRICES_DIR, analysis_suffix)
COMPARISON_DIR <- file.path(TABLES_DIR, "comparisons")

for (d in c(SENS_TABLES_DIR, SENS_FIGURES_DIR, SENS_MATRICES_DIR, COMPARISON_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

message("[hwe_sensitivity] Starting reduced-loci sensitivity branch.")
message("[hwe_sensitivity] Suspect loci requested for removal: ", paste(suspect_null_loci, collapse = ", "))

# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------
normalize_site <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("\\uFEFF", "", x, fixed = TRUE)
  x <- gsub("[[:cntrl:]]", "", x)
  x <- gsub("\\s+", " ", x)
  toupper(x)
}

safe_row_mean <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

extract_overall_stat <- function(overall_obj, stat_name) {
  if (is.null(overall_obj)) return(NA_real_)
  
  if (is.atomic(overall_obj) && !is.null(names(overall_obj))) {
    val <- overall_obj[[stat_name]]
    return(if (length(val) == 0) NA_real_ else as.numeric(val[1]))
  }
  
  if (is.matrix(overall_obj) || is.data.frame(overall_obj)) {
    rn <- rownames(overall_obj)
    cn <- colnames(overall_obj)
    
    if (!is.null(rn) && stat_name %in% rn) return(safe_row_mean(as.numeric(overall_obj[stat_name, , drop = TRUE])))
    if (!is.null(cn) && stat_name %in% cn) return(safe_row_mean(as.numeric(overall_obj[, stat_name, drop = TRUE])))
  }
  
  NA_real_
}

write_csv_msg <- function(df, path) {
  write.csv(df, path, row.names = FALSE)
  message("[hwe_sensitivity] Saved: ", path)
}

matrix_to_long_unique <- function(m, value_name) {
  long_all <- as.data.frame(as.table(m), stringsAsFactors = FALSE)
  names(long_all) <- c("Site1", "Site2", value_name)
  
  long_all %>%
    filter(Site1 != Site2) %>%
    mutate(site_pair = paste(pmin(Site1, Site2), pmax(Site1, Site2), sep = "__")) %>%
    group_by(site_pair) %>%
    summarise(
      Site1 = first(pmin(Site1, Site2)),
      Site2 = first(pmax(Site1, Site2)),
      value = mean(.data[[value_name]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rename(!!value_name := value) %>%
    select(Site1, Site2, all_of(value_name))
}

resolve_site_latitude <- function(meta_df) {
  site_col <- resolve_col_ci(meta_df, c("site", "population", "pop"))
  lat_col <- resolve_col_ci(meta_df, c("latitude", "lat"))
  
  if (is.na(site_col) || is.na(lat_col)) {
    stop("[hwe_sensitivity] Could not find Site and Latitude columns in meta for AMOVA hierarchy.")
  }
  
  data.frame(
    Site = trimws(as.character(meta_df[[site_col]])),
    Latitude = suppressWarnings(as.numeric(meta_df[[lat_col]])),
    stringsAsFactors = FALSE
  ) %>%
    filter(nzchar(Site), !is.na(Latitude)) %>%
    mutate(Site_norm = normalize_site(Site)) %>%
    group_by(Site_norm) %>%
    summarise(
      Site = dplyr::first(Site),
      Latitude = mean(Latitude, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(Latitude, Site)
}

build_latitude_regions <- function(site_levels, meta_df) {
  site_lat_tbl <- resolve_site_latitude(meta_df)
  idx <- match(normalize_site(site_levels), site_lat_tbl$Site_norm)
  
  if (any(is.na(idx))) {
    missing_sites <- site_levels[is.na(idx)]
    stop("[hwe_sensitivity] Missing latitude for AMOVA site(s): ", paste(missing_sites, collapse = ", "))
  }
  
  ranked_tbl <- site_lat_tbl[idx, c("Site_norm", "Site", "Latitude"), drop = FALSE] %>%
    mutate(Site_from_AMOVA = site_levels) %>%
    arrange(Latitude, Site_from_AMOVA) %>%
    mutate(
      Rank_south_to_north = dplyr::row_number(),
      Region = ifelse(Rank_south_to_north <= floor(n() / 2), "South", "North")
    )
  
  setNames(ranked_tbl$Region, ranked_tbl$Site_from_AMOVA)
}

make_reduced_genind <- function(gobj, loci_to_remove) {
  loci_available <- adegenet::locNames(gobj)
  requested <- unique(trimws(as.character(loci_to_remove)))
  requested <- requested[nzchar(requested)]
  
  found <- intersect(requested, loci_available)
  missing <- setdiff(requested, loci_available)
  retained <- setdiff(loci_available, found)
  
  if (length(retained) < 2) {
    stop("[hwe_sensitivity] Too few loci retained after filtering (<2).")
  }
  
  if (length(missing) > 0) {
    warning(
      "[hwe_sensitivity] Requested locus/loci not found and therefore not removed: ",
      paste(missing, collapse = ", ")
    )
  }
  
  if (length(found) == 0) {
    warning("[hwe_sensitivity] None of the requested loci were found; reduced dataset equals full dataset.")
  }
  
  reduced <- local({
    loci_df <- adegenet::genind2df(gobj, sep = "/", usepop = FALSE)
    df_loci <- names(loci_df)
    keep_df <- intersect(retained, df_loci)
    
    if (length(keep_df) != length(retained)) {
      missing_keep <- setdiff(retained, df_loci)
      stop(
        "[hwe_sensitivity] Failed to retain all expected loci when rebuilding genind. Missing: ",
        paste(missing_keep, collapse = ", ")
      )
    }
    
    rebuilt <- adegenet::df2genind(
      X = loci_df[, keep_df, drop = FALSE],
      sep = "/",
      ploidy = adegenet::ploidy(gobj),
      ind.names = adegenet::indNames(gobj),
      type = adegenet::type(gobj),
      NA.char = c("NA", "", "-", "NA/NA", "0/0")
    )
    
    adegenet::pop(rebuilt) <- adegenet::pop(gobj)
    rebuilt
  })
  
  if (adegenet::nLoc(reduced) != length(retained)) {
    stop("[hwe_sensitivity] Locus subsetting mismatch: expected ", length(retained), " got ", adegenet::nLoc(reduced), ".")
  }
  
  message("[hwe_sensitivity] Loci successfully removed: ", if (length(found) > 0) paste(found, collapse = ", ") else "none")
  message("[hwe_sensitivity] Loci not found: ", if (length(missing) > 0) paste(missing, collapse = ", ") else "none")
  message("[hwe_sensitivity] Locus count before filtering: ", length(loci_available))
  message("[hwe_sensitivity] Locus count after filtering: ", length(retained))
  
  list(
    reduced = reduced,
    found = found,
    missing = missing,
    retained = retained,
    all_full = loci_available
  )
}

run_diversity_bundle <- function(gobj, dataset_label) {
  hf <- hierfstat::genind2hierfstat(gobj)
  bs <- hierfstat::basic.stats(hf)
  
  site_levels <- rownames(bs$Ho)
  site_n_tbl <- table(as.character(adegenet::pop(gobj)))
  
  by_site <- data.frame(
    Dataset = dataset_label,
    Site = site_levels,
    N = as.integer(site_n_tbl[site_levels]),
    Ho = apply(bs$Ho, 1, safe_row_mean),
    He = apply(bs$Hs, 1, safe_row_mean),
    FIS = apply(bs$Fis, 1, safe_row_mean),
    stringsAsFactors = FALSE
  )
  
  overall_ho <- extract_overall_stat(bs$overall, "Ho")
  overall_he <- extract_overall_stat(bs$overall, "Hs")
  overall_fis <- extract_overall_stat(bs$overall, "Fis")
  
  if (is.na(overall_ho)) overall_ho <- safe_row_mean(by_site$Ho)
  if (is.na(overall_he)) overall_he <- safe_row_mean(by_site$He)
  if (is.na(overall_fis)) overall_fis <- safe_row_mean(by_site$FIS)
  
  overall <- data.frame(
    Dataset = dataset_label,
    N = adegenet::nInd(gobj),
    N_loci = adegenet::nLoc(gobj),
    Ho = overall_ho,
    He = overall_he,
    FIS = overall_fis,
    stringsAsFactors = FALSE
  )
  
  pop_sizes <- table(adegenet::pop(gobj))
  min_n <- min(pop_sizes)
  ar <- if (!is.na(min_n) && min_n >= 2) hierfstat::allelic.richness(hf, min.n = min_n) else NULL
  
  if (is.null(ar)) {
    ar_by_site <- data.frame(
      Dataset = dataset_label,
      Site = site_levels,
      Allelic_Richness = NA_real_,
      Allelic_Richness_SE = NA_real_,
      stringsAsFactors = FALSE
    )
  } else {
    Ar <- ar$Ar
    ar_by_site <- data.frame(
      Dataset = dataset_label,
      Site = colnames(Ar),
      Allelic_Richness = as.numeric(colMeans(Ar, na.rm = TRUE)),
      Allelic_Richness_SE = as.numeric(apply(Ar, 2, sd, na.rm = TRUE) / sqrt(nrow(Ar))),
      stringsAsFactors = FALSE
    )
  }
  
  list(by_site = by_site, overall = overall, allelic_richness = ar_by_site)
}

run_differentiation_bundle <- function(gobj, dataset_label) {
  jost_mat <- as.matrix(mmod::pairwise_D(gobj, linearized = FALSE))
  diag(jost_mat) <- 0
  
  hf <- hierfstat::genind2hierfstat(gobj)
  if (!is.numeric(hf[[1]])) hf[[1]] <- as.integer(factor(hf[[1]]))
  
  fst_raw <- tryCatch(
    hierfstat::pairwise.WCfst(hf),
    error = function(e) hierfstat::pairwise.neifst(hf)
  )
  fst_mat <- as.matrix(fst_raw)
  diag(fst_mat) <- 0
  
  site_lookup <- tapply(as.character(adegenet::pop(gobj)), hf[[1]], function(x) names(sort(table(x), decreasing = TRUE))[1])
  site_lookup <- as.character(site_lookup)
  if (length(site_lookup) == nrow(fst_mat)) {
    rownames(fst_mat) <- site_lookup
    colnames(fst_mat) <- site_lookup
  }
  
  jost_long <- matrix_to_long_unique(jost_mat, "JostD") %>% mutate(Dataset = dataset_label)
  fst_long <- matrix_to_long_unique(fst_mat, "FST") %>% mutate(Dataset = dataset_label)
  
  summary_tbl <- data.frame(
    Dataset = dataset_label,
    mean_pairwise_JostD = mean(jost_long$JostD, na.rm = TRUE),
    median_pairwise_JostD = median(jost_long$JostD, na.rm = TRUE),
    mean_pairwise_FST = mean(fst_long$FST, na.rm = TRUE),
    median_pairwise_FST = median(fst_long$FST, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  
  list(jost_mat = jost_mat, fst_mat = fst_mat, jost_long = jost_long, fst_long = fst_long, summary = summary_tbl)
}

run_amova_bundle <- function(gi_use, dataset_label) {
  validate_columns(df_ids_mll, c("ind_id", "Site"), df_name = "[hwe_sensitivity] df_ids_mll")
  id_to_site <- setNames(as.character(df_ids_mll$Site), normalize_id(df_ids_mll$ind_id))
  
  inds <- adegenet::indNames(gi_use)
  site_from_dfids <- id_to_site[normalize_id(inds)]
  site_labels <- if (all(!is.na(site_from_dfids))) site_from_dfids else as.character(adegenet::pop(gi_use))
  
  valid <- !is.na(site_labels) & nzchar(site_labels)
  gi_f <- gi_use[valid, , drop = FALSE]
  site_f <- droplevels(as.factor(site_labels[valid]))
  
  group_tab <- table(site_f)
  keep_groups <- names(group_tab)[group_tab >= 2]
  if (length(keep_groups) < 2) {
    stop("[hwe_sensitivity] Need >=2 sites with >=2 individuals for AMOVA.")
  }
  
  keep_idx <- site_f %in% keep_groups
  gi_f <- gi_f[keep_idx, , drop = FALSE]
  site_f <- droplevels(site_f[keep_idx])
  adegenet::pop(gi_f) <- site_f
  
  run_amova_model <- function(gi_obj, strata_df, formula_obj, model_label) {
    adegenet::strata(gi_obj) <- strata_df
    fit <- poppr::poppr.amova(gi_obj, formula_obj)
    rand <- ade4::randtest(fit, nrepet = AMOVA_PERMUTATIONS)
    
    components <- as.data.frame(fit$componentsofcovariance)
    components$Source <- rownames(components)
    rownames(components) <- NULL
    names(components)[1] <- "Sigma"
    
    phi_stats <- as.data.frame(fit$statphi)
    phi_stats$Source <- rownames(phi_stats)
    rownames(phi_stats) <- NULL
    names(phi_stats)[1] <- "Phi"
    
    results <- full_join(components, phi_stats, by = "Source") %>%
      mutate(
        Dataset = dataset_label,
        Model = model_label,
        N_individuals_used = adegenet::nInd(gi_obj),
        N_groups_used = dplyr::n_distinct(strata_df$Site),
        Permutations = AMOVA_PERMUTATIONS
      ) %>%
      select(Dataset, Model, Source, everything())
    
    rand_df <- data.frame(
      Dataset = dataset_label,
      Model = model_label,
      component = names(rand$obs),
      statistic = as.numeric(rand$obs),
      p_value = as.numeric(rand$pvalue),
      permutations = AMOVA_PERMUTATIONS,
      stringsAsFactors = FALSE
    )
    
    list(results = results, rand = rand_df)
  }
  
  strata_site <- data.frame(
    pop = site_f,
    Site = site_f,
    row.names = adegenet::indNames(gi_f),
    stringsAsFactors = TRUE
  )
  
  site_fit <- run_amova_model(gi_f, strata_site, ~pop, "Site_only")
  amova_results <- site_fit$results
  amova_rand <- site_fit$rand
  
  site_region_map <- build_latitude_regions(levels(site_f), meta)
  region_f <- factor(site_region_map[as.character(site_f)], levels = c("South", "North"))
  
  if (nlevels(droplevels(region_f)) >= 2) {
    keep_region <- !is.na(region_f)
    gi_h <- gi_f[keep_region, , drop = FALSE]
    site_h <- droplevels(site_f[keep_region])
    region_h <- droplevels(region_f[keep_region])
    adegenet::pop(gi_h) <- site_h
    
    strata_h <- data.frame(
      Region = region_h,
      Site = site_h,
      pop = site_h,
      row.names = adegenet::indNames(gi_h),
      stringsAsFactors = TRUE
    )
    
    h_fit <- run_amova_model(gi_h, strata_h, ~Region/Site, "NorthSouth_Site_hierarchical")
    amova_results <- bind_rows(amova_results, h_fit$results)
    amova_rand <- bind_rows(amova_rand, h_fit$rand)
  }
  
  list(results = amova_results, rand = amova_rand)
}

run_hwe_by_site_locus <- function(gobj, dataset_label) {
  site_vec <- as.character(adegenet::pop(gobj))
  loci_df <- adegenet::genind2df(gobj, sep = "/", usepop = FALSE)
  loci_names <- names(loci_df)
  
  all_rows <- list()
  idx <- 1L
  
  for (site in sort(unique(site_vec))) {
    use <- site_vec == site
    if (sum(use) < 2) next
    
    site_df <- loci_df[use, , drop = FALSE]
    for (loc in loci_names) {
      vals <- trimws(as.character(site_df[[loc]]))
      vals[vals %in% c("", "NA", "0", "0/0", "NA/NA", "-")] <- NA_character_
      vals <- vals[!is.na(vals)]
      if (length(vals) < 2) next
      
      geno_fac <- factor(vals)
      p_val <- tryCatch(
        pegas::hw.test(geno_fac, B = HWE_MONTE_CARLO_REPS)$p.value,
        error = function(e) NA_real_
      )
      
      all_rows[[idx]] <- data.frame(
        Dataset = dataset_label,
        Site = site,
        Locus = loc,
        N_non_missing = length(vals),
        N_genotype_classes = nlevels(geno_fac),
        p_value_raw = as.numeric(p_val),
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }
  
  out <- bind_rows(all_rows)
  if (nrow(out) == 0) return(out)
  
  out <- out %>%
    mutate(
      p_value_adj_fdr = p.adjust(p_value_raw, method = "BH"),
      hwe_reject_raw_0_05 = !is.na(p_value_raw) & p_value_raw < 0.05,
      hwe_reject_fdr_0_05 = !is.na(p_value_adj_fdr) & p_value_adj_fdr < 0.05
    )
  
  out
}

run_pca_bundle <- function(gobj, dataset_label) {
  X <- adegenet::tab(gobj, freq = TRUE, NA.method = "mean")
  keep_cols <- apply(X, 2, function(v) stats::var(v, na.rm = TRUE) > 0)
  X <- X[, keep_cols, drop = FALSE]
  
  pca_fit <- stats::prcomp(X, center = TRUE, scale. = FALSE)
  var_exp <- (pca_fit$sdev^2 / sum(pca_fit$sdev^2)) * 100
  
  scores <- as.data.frame(pca_fit$x[, 1:2, drop = FALSE], stringsAsFactors = FALSE)
  scores$Individual <- rownames(scores)
  scores$Site <- as.character(adegenet::pop(gobj))
  scores$Dataset <- dataset_label
  
  variance <- data.frame(
    Dataset = dataset_label,
    PC = paste0("PC", seq_along(var_exp)),
    Percent_Variance = as.numeric(var_exp),
    stringsAsFactors = FALSE
  )
  
  list(scores = scores, variance = variance)
}

run_dapc_bundle <- function(gobj, dataset_label) {
  grp <- as.factor(adegenet::pop(gobj))
  dapc_fit <- adegenet::dapc(gobj, pop = grp, n.pca = min(50, nInd(gobj) - 1), n.da = min(nlevels(grp) - 1, 10))
  
  coords <- as.data.frame(dapc_fit$ind.coord[, 1:2, drop = FALSE], stringsAsFactors = FALSE)
  coords$Individual <- rownames(coords)
  coords$Site <- as.character(grp)
  coords$Dataset <- dataset_label
  
  if (ncol(coords) >= 2) names(coords)[1:2] <- c("LD1", "LD2")
  
  eig <- as.numeric(dapc_fit$eig)
  eig_df <- data.frame(
    Dataset = dataset_label,
    Axis = paste0("LD", seq_along(eig)),
    Eigenvalue = eig,
    stringsAsFactors = FALSE
  )
  
  list(coords = coords, eig = eig_df)
}

compare_two_dataset_table <- function(df, by_cols, value_cols) {
  wide <- df %>%
    select(Dataset, all_of(by_cols), all_of(value_cols)) %>%
    pivot_wider(names_from = Dataset, values_from = all_of(value_cols), names_sep = "__")
  
  for (v in value_cols) {
    full_col <- paste0(v, "__FULL")
    red_col <- paste0(v, "__REDUCED")
    delta_col <- paste0("delta_", v, "_REDUCED_minus_FULL")
    if (all(c(full_col, red_col) %in% names(wide))) {
      wide[[delta_col]] <- suppressWarnings(as.numeric(wide[[red_col]]) - as.numeric(wide[[full_col]]))
    }
  }
  
  wide
}

# ------------------------------------------------------------------
# Build reduced object and manifest
# ------------------------------------------------------------------
full_gi <- gi_mll
reduced_info <- make_reduced_genind(full_gi, suspect_null_loci)
reduced_gi <- reduced_info$reduced

locus_manifest <- data.frame(
  dataset = c("FULL", "REDUCED"),
  n_loci = c(adegenet::nLoc(full_gi), adegenet::nLoc(reduced_gi)),
  loci = c(
    paste(adegenet::locNames(full_gi), collapse = ";"),
    paste(adegenet::locNames(reduced_gi), collapse = ";")
  ),
  removed_loci = c("", paste(reduced_info$found, collapse = ";")),
  requested_but_not_found = c("", paste(reduced_info$missing, collapse = ";")),
  stringsAsFactors = FALSE
)

write_csv_msg(locus_manifest, file.path(SENS_TABLES_DIR, paste0("locus_manifest_", analysis_suffix, ".csv")))

# ------------------------------------------------------------------
# Run parallel full/reduced analyses
# ------------------------------------------------------------------
full_div <- run_diversity_bundle(full_gi, "FULL")
red_div <- run_diversity_bundle(reduced_gi, "REDUCED")

full_diff <- run_differentiation_bundle(full_gi, "FULL")
red_diff <- run_differentiation_bundle(reduced_gi, "REDUCED")

full_amova <- run_amova_bundle(full_gi, "FULL")
red_amova <- run_amova_bundle(reduced_gi, "REDUCED")

full_hwe <- run_hwe_by_site_locus(full_gi, "FULL")
red_hwe <- run_hwe_by_site_locus(reduced_gi, "REDUCED")

full_pca <- run_pca_bundle(full_gi, "FULL")
red_pca <- run_pca_bundle(reduced_gi, "REDUCED")

full_dapc <- run_dapc_bundle(full_gi, "FULL")
red_dapc <- run_dapc_bundle(reduced_gi, "REDUCED")

# ------------------------------------------------------------------
# Write reduced-only outputs requested for direct comparison with main outputs
# ------------------------------------------------------------------
write_csv_msg(red_div$by_site, file.path(SENS_TABLES_DIR, paste0("heterozygosity_fis_by_site_", analysis_suffix, ".csv")))
write_csv_msg(red_div$overall, file.path(SENS_TABLES_DIR, paste0("heterozygosity_fis_overall_", analysis_suffix, ".csv")))
write_csv_msg(red_div$allelic_richness, file.path(SENS_TABLES_DIR, paste0("allelic_richness_by_site_", analysis_suffix, ".csv")))
write_csv_msg(red_diff$fst_long, file.path(SENS_TABLES_DIR, paste0("pairwise_fst_", analysis_suffix, ".csv")))
write_csv_msg(red_diff$jost_long, file.path(SENS_TABLES_DIR, paste0("pairwise_jostD_", analysis_suffix, ".csv")))
write_csv_msg(red_amova$results, file.path(SENS_TABLES_DIR, paste0("amova_", analysis_suffix, ".csv")))
write_csv_msg(red_amova$rand, file.path(SENS_TABLES_DIR, paste0("amova_randtest_", analysis_suffix, ".csv")))
write_csv_msg(red_hwe, file.path(SENS_TABLES_DIR, paste0("hwe_by_site_by_locus_", analysis_suffix, ".csv")))
write_csv_msg(red_pca$variance, file.path(SENS_TABLES_DIR, paste0("pca_variance_", analysis_suffix, ".csv")))
write_csv_msg(red_dapc$eig, file.path(SENS_TABLES_DIR, paste0("dapc_eigenvalues_", analysis_suffix, ".csv")))

write.csv(red_diff$fst_mat, file.path(SENS_MATRICES_DIR, paste0("pairwise_fst_", analysis_suffix, ".csv")))
write.csv(red_diff$jost_mat, file.path(SENS_MATRICES_DIR, paste0("pairwise_jostD_", analysis_suffix, ".csv")))

# ------------------------------------------------------------------
# Write comparison outputs (FULL vs REDUCED)
# ------------------------------------------------------------------
div_by_site_all <- bind_rows(full_div$by_site, red_div$by_site)
div_overall_all <- bind_rows(full_div$overall, red_div$overall)
ar_all <- bind_rows(full_div$allelic_richness, red_div$allelic_richness)

diff_summary_all <- bind_rows(full_diff$summary, red_diff$summary)
fst_long_all <- bind_rows(full_diff$fst_long, red_diff$fst_long)
jost_long_all <- bind_rows(full_diff$jost_long, red_diff$jost_long)

amova_results_all <- bind_rows(full_amova$results, red_amova$results)
amova_rand_all <- bind_rows(full_amova$rand, red_amova$rand)

hwe_all <- bind_rows(full_hwe, red_hwe)
pca_var_all <- bind_rows(full_pca$variance, red_pca$variance)
dapc_eig_all <- bind_rows(full_dapc$eig, red_dapc$eig)

summary_stats <- data.frame(
  Dataset = c("FULL", "REDUCED"),
  N_individuals = c(nInd(full_gi), nInd(reduced_gi)),
  N_loci = c(nLoc(full_gi), nLoc(reduced_gi)),
  stringsAsFactors = FALSE
)

summary_stats_cmp <- compare_two_dataset_table(summary_stats, by_cols = character(0), value_cols = c("N_individuals", "N_loci"))
div_overall_cmp <- compare_two_dataset_table(div_overall_all, by_cols = character(0), value_cols = c("N", "N_loci", "Ho", "He", "FIS"))
div_site_cmp <- compare_two_dataset_table(div_by_site_all, by_cols = c("Site"), value_cols = c("N", "Ho", "He", "FIS"))
ar_cmp <- compare_two_dataset_table(ar_all, by_cols = c("Site"), value_cols = c("Allelic_Richness", "Allelic_Richness_SE"))
diff_cmp <- compare_two_dataset_table(diff_summary_all, by_cols = character(0), value_cols = c("mean_pairwise_JostD", "median_pairwise_JostD", "mean_pairwise_FST", "median_pairwise_FST"))
pca_cmp <- compare_two_dataset_table(pca_var_all, by_cols = c("PC"), value_cols = c("Percent_Variance"))
dapc_cmp <- compare_two_dataset_table(dapc_eig_all, by_cols = c("Axis"), value_cols = c("Eigenvalue"))

write_csv_msg(summary_stats_cmp, file.path(COMPARISON_DIR, "full_vs_noSuspectLoci_locus_count_comparison.csv"))
write_csv_msg(div_overall_cmp, file.path(COMPARISON_DIR, "full_vs_noSuspectLoci_heterozygosity_fis_overall_comparison.csv"))
write_csv_msg(div_site_cmp, file.path(COMPARISON_DIR, "full_vs_noSuspectLoci_heterozygosity_fis_by_site_comparison.csv"))
write_csv_msg(ar_cmp, file.path(COMPARISON_DIR, "full_vs_noSuspectLoci_allelic_richness_comparison.csv"))
write_csv_msg(diff_cmp, file.path(COMPARISON_DIR, "full_vs_noSuspectLoci_differentiation_comparison.csv"))
write_csv_msg(pca_cmp, file.path(COMPARISON_DIR, "full_vs_noSuspectLoci_pca_variance_comparison.csv"))
write_csv_msg(dapc_cmp, file.path(COMPARISON_DIR, "full_vs_noSuspectLoci_dapc_comparison.csv"))

write_csv_msg(fst_long_all, file.path(COMPARISON_DIR, "full_vs_noSuspectLoci_pairwise_fst_long.csv"))
write_csv_msg(jost_long_all, file.path(COMPARISON_DIR, "full_vs_noSuspectLoci_pairwise_jostD_long.csv"))
write_csv_msg(amova_results_all, file.path(COMPARISON_DIR, "full_vs_noSuspectLoci_amova_long.csv"))
write_csv_msg(amova_rand_all, file.path(COMPARISON_DIR, "full_vs_noSuspectLoci_amova_randtest_long.csv"))
write_csv_msg(hwe_all, file.path(COMPARISON_DIR, "full_vs_noSuspectLoci_hwe_by_site_by_locus_long.csv"))
write_csv_msg(pca_var_all, file.path(COMPARISON_DIR, "full_vs_noSuspectLoci_pca_variance_long.csv"))
write_csv_msg(dapc_eig_all, file.path(COMPARISON_DIR, "full_vs_noSuspectLoci_dapc_eigenvalues_long.csv"))

# ------------------------------------------------------------------
# Figures for reduced-only and comparison views
# ------------------------------------------------------------------
pca_scores_all <- bind_rows(full_pca$scores, red_pca$scores)
dapc_coords_all <- bind_rows(full_dapc$coords, red_dapc$coords)

pca_plot <- ggplot(pca_scores_all, aes(PC1, PC2, color = Site)) +
  geom_point(alpha = 0.8, size = 1.8) +
  facet_wrap(~ Dataset, scales = "free") +
  theme_bw(base_size = 11) +
  labs(title = "PCA comparison: FULL vs noSuspectLoci", x = "PC1", y = "PC2")

ggsave(file.path(SENS_FIGURES_DIR, paste0("pca_", analysis_suffix, ".pdf")), pca_plot, width = 8.2, height = 4.6)

dapc_plot <- ggplot(dapc_coords_all, aes(LD1, LD2, color = Site)) +
  geom_point(alpha = 0.8, size = 1.8) +
  facet_wrap(~ Dataset, scales = "free") +
  theme_bw(base_size = 11) +
  labs(title = "DAPC comparison: FULL vs noSuspectLoci", x = "LD1", y = "LD2")

ggsave(file.path(SENS_FIGURES_DIR, paste0("dapc_", analysis_suffix, ".pdf")), dapc_plot, width = 8.2, height = 4.6)

delta_tbl <- bind_rows(
  div_overall_cmp %>% transmute(metric = "Ho", delta = delta_Ho_REDUCED_minus_FULL),
  div_overall_cmp %>% transmute(metric = "He", delta = delta_He_REDUCED_minus_FULL),
  div_overall_cmp %>% transmute(metric = "FIS", delta = delta_FIS_REDUCED_minus_FULL),
  diff_cmp %>% transmute(metric = "mean_pairwise_FST", delta = delta_mean_pairwise_FST_REDUCED_minus_FULL),
  diff_cmp %>% transmute(metric = "mean_pairwise_JostD", delta = delta_mean_pairwise_JostD_REDUCED_minus_FULL)
)

delta_plot <- ggplot(delta_tbl, aes(metric, delta, fill = metric)) +
  geom_col(show.legend = FALSE, width = 0.72) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  theme_bw(base_size = 11) +
  labs(title = "Sensitivity deltas (REDUCED - FULL)", x = NULL, y = "Delta")

ggsave(file.path(SENS_FIGURES_DIR, "full_vs_noSuspectLoci_metric_deltas.pdf"), delta_plot, width = 8, height = 4.5)

message("[hwe_sensitivity] Completed reduced-loci sensitivity branch successfully.")