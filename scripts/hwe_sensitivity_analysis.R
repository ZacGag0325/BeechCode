# scripts/hwe_sensitivity_analysis.R
############################################################
# HWE sensitivity analysis (FULL vs REDUCED loci sets)
#
# Purpose:
# - Re-run core population-genetic analyses with the exact same object lineage,
#   site definitions, and missing-data handling already used in BeechCode.
# - Compare results between:
#   1) FULL dataset (all loci in gi_mll)
#   2) REDUCED dataset (excluding user-specified HWE-deviating loci)
#
# Loci excluded in REDUCED dataset (provided by user):
# - EJV8T_A_0
# - ERHBI_A_0
# - FCM5
# - FG5
#
# Outputs:
# - results/hwe_sensitivity/*.csv
# - figures/hwe_sensitivity/*.jpeg
# - results/hwe_sensitivity/hwe_sensitivity_summary.md
#
# Notes:
# - This is a robustness/sensitivity analysis only.
# - It does NOT prove null alleles or Wahlund effect.
############################################################

suppressPackageStartupMessages({
  library(adegenet)
  library(poppr)
  library(ade4)
  library(hierfstat)
  library(mmod)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

source("scripts/_load_objects.R")

# ------------------------------------------------------------------
# User-configurable section
# ------------------------------------------------------------------
LOCI_TO_EXCLUDE <- c("EJV8T_A_0", "ERHBI_A_0", "FCM5", "FG5")
AMOVA_PERMUTATIONS <- 999

# ------------------------------------------------------------------
# Output folders requested by user
# ------------------------------------------------------------------
RESULTS_DIR <- file.path(PROJECT_ROOT, "results", "hwe_sensitivity")
FIGURES_HWE_DIR <- file.path(PROJECT_ROOT, "figures", "hwe_sensitivity")
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURES_HWE_DIR, recursive = TRUE, showWarnings = FALSE)

message("[hwe_sensitivity] Results directory: ", RESULTS_DIR)
message("[hwe_sensitivity] Figures directory: ", FIGURES_HWE_DIR)

# ------------------------------------------------------------------
# Defensive helpers
# ------------------------------------------------------------------
normalize_site <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("\\uFEFF", "", x, fixed = TRUE)
  x <- gsub("[[:cntrl:]]", "", x)
  x <- gsub("\\s+", " ", x)
  toupper(x)
}

validate_locus_filter <- function(gobj, loci_to_remove) {
  loci_available <- adegenet::locNames(gobj)
  missing <- setdiff(loci_to_remove, loci_available)
  placeholders <- grepl("^\\[LOCUS[0-9]+\\]$", loci_to_remove)
  
  if (any(placeholders)) {
    stop(
      "[hwe_sensitivity] LOCI_TO_EXCLUDE still contains placeholder names.",
      " Replace [LOCUS1]...[LOCUS4] with exact locus names from your dataset."
    )
  }
  
  if (length(unique(loci_to_remove)) != length(loci_to_remove)) {
    stop("[hwe_sensitivity] LOCI_TO_EXCLUDE contains duplicates. Provide 4 unique locus names.")
  }
  
  if (length(missing) > 0) {
    stop(
      "[hwe_sensitivity] Requested loci not found: ", paste(missing, collapse = ", "),
      "\nAvailable loci are: ", paste(loci_available, collapse = ", ")
    )
  }
  
  retained <- setdiff(loci_available, loci_to_remove)
  if (length(retained) < 2) {
    stop("[hwe_sensitivity] Too few loci retained after filtering (<2).")
  }
  
  list(removed = loci_to_remove, retained = retained)
}

subset_genind_loci <- function(gobj, loci_keep) {
  keep_idx <- adegenet::locNames(gobj) %in% loci_keep
  out <- gobj[, keep_idx, drop = FALSE]
  if (!inherits(out, "genind")) stop("[hwe_sensitivity] Locus subsetting failed (not a genind object).")
  out
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
  
  if (length(overall_obj) == 1 && is.finite(suppressWarnings(as.numeric(overall_obj)))) {
    return(as.numeric(overall_obj))
  }
  
  NA_real_
}

matrix_to_long_unique <- function(m, value_name) {
  long_all <- as.data.frame(as.table(m), stringsAsFactors = FALSE)
  names(long_all) <- c("Site1", "Site2", value_name)
  
  long_unique <- long_all %>%
    filter(Site1 != Site2) %>%
    mutate(site_pair = paste(pmin(Site1, Site2), pmax(Site1, Site2), sep = "__")) %>%
    group_by(site_pair) %>%
    summarise(
      Site1 = first(pmin(Site1, Site2)),
      Site2 = first(pmax(Site1, Site2)),
      value = mean(.data[[value_name]], na.rm = TRUE),
      .groups = "drop"
    )
  
  names(long_unique)[names(long_unique) == "value"] <- value_name
  long_unique %>% select(Site1, Site2, all_of(value_name))
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

run_amova_bundle <- function(gi_use, grouping_source, n_drop, dataset_label) {
  run_amova_model <- function(gi_obj, strata_df, model_formula, model_label, grouping_src, n_drop_local) {
    adegenet::strata(gi_obj) <- strata_df
    fit <- poppr::poppr.amova(gi_obj, model_formula)
    
    rand <- tryCatch(
      ade4::randtest(fit, nrepet = AMOVA_PERMUTATIONS),
      error = function(e) {
        warning("[hwe_sensitivity] randtest failed for ", dataset_label, " / ", model_label, ": ", conditionMessage(e))
        NULL
      }
    )
    
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
        Grouping_source = grouping_src,
        Permutations = AMOVA_PERMUTATIONS
      ) %>%
      select(Dataset, Model, Source, everything())
    
    if (is.null(rand)) {
      rand_df <- data.frame(
        Dataset = dataset_label,
        model = model_label,
        test = "AMOVA_randtest",
        component = NA_character_,
        statistic = NA_real_,
        p_value = NA_real_,
        permutations = AMOVA_PERMUTATIONS,
        grouping_source = grouping_src,
        n_individuals_used = adegenet::nInd(gi_obj),
        n_groups_used = dplyr::n_distinct(strata_df$Site),
        dropped_for_missing_group = as.integer(n_drop_local),
        note = "randtest_failed",
        stringsAsFactors = FALSE
      )
    } else {
      obs_vec <- as.numeric(rand$obs)
      p_vec <- as.numeric(rand$pvalue)
      n_comp <- max(length(obs_vec), length(p_vec), 1)
      if (length(obs_vec) == 0) obs_vec <- rep(NA_real_, n_comp)
      if (length(p_vec) == 0) p_vec <- rep(NA_real_, n_comp)
      if (length(obs_vec) < n_comp) obs_vec <- c(obs_vec, rep(NA_real_, n_comp - length(obs_vec)))
      if (length(p_vec) < n_comp) p_vec <- c(p_vec, rep(NA_real_, n_comp - length(p_vec)))
      
      comp_names <- names(rand$obs)
      if (is.null(comp_names) || length(comp_names) == 0) comp_names <- paste0("component_", seq_len(n_comp))
      if (length(comp_names) < n_comp) comp_names <- c(comp_names, paste0("component_", (length(comp_names) + 1):n_comp))
      
      rand_df <- data.frame(
        Dataset = dataset_label,
        model = model_label,
        test = "AMOVA_randtest",
        component = comp_names[seq_len(n_comp)],
        statistic = obs_vec[seq_len(n_comp)],
        p_value = p_vec[seq_len(n_comp)],
        permutations = AMOVA_PERMUTATIONS,
        grouping_source = grouping_src,
        n_individuals_used = adegenet::nInd(gi_obj),
        n_groups_used = dplyr::n_distinct(strata_df$Site),
        dropped_for_missing_group = as.integer(n_drop_local),
        note = "ok",
        stringsAsFactors = FALSE
      )
    }
    
    list(results = results, rand = rand_df)
  }
  
  validate_columns(df_ids_mll, c("ind_id", "Site"), df_name = "[hwe_sensitivity] df_ids_mll")
  id_to_site <- setNames(as.character(df_ids_mll$Site), normalize_id(df_ids_mll$ind_id))
  
  inds <- adegenet::indNames(gi_use)
  site_from_dfids <- id_to_site[normalize_id(inds)]
  site_from_pop <- as.character(adegenet::pop(gi_use))
  
  site_labels <- if (all(!is.na(site_from_dfids))) site_from_dfids else site_from_pop
  valid <- !is.na(site_labels) & nzchar(site_labels)
  n_drop_local <- sum(!valid)
  
  gi_f <- gi_use[valid, , drop = FALSE]
  site_f <- droplevels(as.factor(site_labels[valid]))
  
  group_tab <- table(site_f)
  keep_groups <- names(group_tab)[group_tab >= 2]
  if (length(keep_groups) < 2) stop("[hwe_sensitivity] Need >=2 sites with >=2 individuals for AMOVA.")
  
  keep_idx <- site_f %in% keep_groups
  gi_f <- gi_f[keep_idx, , drop = FALSE]
  site_f <- droplevels(site_f[keep_idx])
  
  adegenet::pop(gi_f) <- site_f
  
  strata_site <- data.frame(
    pop = site_f,
    Site = site_f,
    row.names = adegenet::indNames(gi_f),
    stringsAsFactors = TRUE
  )
  
  site_fit <- run_amova_model(gi_f, strata_site, ~pop, "Site_only", grouping_source, n_drop_local)
  amova_results <- site_fit$results
  amova_rand <- site_fit$rand
  
  site_region_map <- build_latitude_regions(levels(site_f), meta)
  region_f <- factor(site_region_map[as.character(site_f)], levels = c("South", "North"))
  
  has_region <- !all(is.na(region_f) | !nzchar(as.character(region_f)))
  if (has_region && nlevels(droplevels(region_f)) >= 2) {
    keep_region <- !is.na(region_f) & nzchar(as.character(region_f))
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
    
    h_fit <- run_amova_model(
      gi_h,
      strata_h,
      ~Region/Site,
      "NorthSouth_Site_hierarchical",
      "derived_from_site_latitude_rank",
      n_drop_local
    )
    
    amova_results <- bind_rows(amova_results, h_fit$results)
    amova_rand <- bind_rows(amova_rand, h_fit$rand)
  }
  
  list(results = amova_results, randtest = amova_rand)
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
  
  list(jost_long = jost_long, fst_long = fst_long, summary = summary_tbl)
}

run_pca_bundle <- function(gobj, dataset_label) {
  X <- adegenet::tab(gobj, freq = TRUE, NA.method = "mean")
  keep_cols <- apply(X, 2, function(v) stats::var(v, na.rm = TRUE) > 0)
  X <- X[, keep_cols, drop = FALSE]
  
  pca_fit <- stats::prcomp(X, center = TRUE, scale. = FALSE)
  var_exp <- (pca_fit$sdev^2 / sum(pca_fit$sdev^2)) * 100
  
  out <- data.frame(
    Dataset = dataset_label,
    PC = paste0("PC", seq_along(var_exp)),
    Percent_Variance = as.numeric(var_exp),
    stringsAsFactors = FALSE
  )
  
  out
}

count_alleles_by_locus <- function(gobj, dataset_label) {
  geno_df <- adegenet::genind2df(gobj, sep = "/")
  loci <- adegenet::locNames(gobj)
  
  parse_alleles <- function(x) {
    parts <- strsplit(as.character(x), split = "/", fixed = TRUE)
    unique(unlist(parts, use.names = FALSE))
  }
  
  out <- lapply(loci, function(loc) {
    vals <- geno_df[[loc]]
    vals <- vals[!is.na(vals) & nzchar(vals)]
    alleles <- parse_alleles(vals)
    alleles <- alleles[!is.na(alleles) & nzchar(alleles)]
    data.frame(Dataset = dataset_label, Locus = loc, N_alleles = dplyr::n_distinct(alleles), stringsAsFactors = FALSE)
  })
  
  bind_rows(out)
}

# ------------------------------------------------------------------
# Build FULL and REDUCED datasets (same source object + same Site labels)
# ------------------------------------------------------------------
full_gi <- gi_mll
locus_check <- validate_locus_filter(full_gi, LOCI_TO_EXCLUDE)
reduced_gi <- subset_genind_loci(full_gi, loci_keep = locus_check$retained)

message("[hwe_sensitivity] FULL loci count: ", adegenet::nLoc(full_gi))
message("[hwe_sensitivity] REDUCED loci count: ", adegenet::nLoc(reduced_gi))
message("[hwe_sensitivity] Removed loci: ", paste(locus_check$removed, collapse = ", "))
message("[hwe_sensitivity] Retained loci: ", paste(locus_check$retained, collapse = ", "))

locus_manifest <- data.frame(
  dataset = c("FULL", "REDUCED"),
  n_loci = c(adegenet::nLoc(full_gi), adegenet::nLoc(reduced_gi)),
  loci = c(
    paste(adegenet::locNames(full_gi), collapse = ";"),
    paste(adegenet::locNames(reduced_gi), collapse = ";")
  ),
  removed_loci = c("", paste(locus_check$removed, collapse = ";")),
  stringsAsFactors = FALSE
)
write.csv(locus_manifest, file.path(RESULTS_DIR, "locus_manifest.csv"), row.names = FALSE)

# ------------------------------------------------------------------
# Run analysis bundles
# ------------------------------------------------------------------
full_amova <- run_amova_bundle(full_gi, grouping_source = "df_ids$Site", n_drop = 0, dataset_label = "FULL")
red_amova <- run_amova_bundle(reduced_gi, grouping_source = "df_ids$Site", n_drop = 0, dataset_label = "REDUCED")

amova_results_all <- bind_rows(full_amova$results, red_amova$results)
amova_rand_all <- bind_rows(full_amova$randtest, red_amova$randtest)

full_div <- run_diversity_bundle(full_gi, "FULL")
red_div <- run_diversity_bundle(reduced_gi, "REDUCED")

div_by_site_all <- bind_rows(full_div$by_site, red_div$by_site)
div_overall_all <- bind_rows(full_div$overall, red_div$overall)
ar_all <- bind_rows(full_div$allelic_richness, red_div$allelic_richness)

full_diff <- run_differentiation_bundle(full_gi, "FULL")
red_diff <- run_differentiation_bundle(reduced_gi, "REDUCED")

diff_summary_all <- bind_rows(full_diff$summary, red_diff$summary)
fst_long_all <- bind_rows(full_diff$fst_long, red_diff$fst_long)
jost_long_all <- bind_rows(full_diff$jost_long, red_diff$jost_long)

pca_var_all <- bind_rows(run_pca_bundle(full_gi, "FULL"), run_pca_bundle(reduced_gi, "REDUCED"))
allele_counts_all <- bind_rows(count_alleles_by_locus(full_gi, "FULL"), count_alleles_by_locus(reduced_gi, "REDUCED"))

# ------------------------------------------------------------------
# Side-by-side comparison tables
# ------------------------------------------------------------------
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

amova_side_by_side <- compare_two_dataset_table(
  amova_results_all,
  by_cols = c("Model", "Source"),
  value_cols = setdiff(names(amova_results_all), c("Dataset", "Model", "Source"))
)

amova_rand_side_by_side <- compare_two_dataset_table(
  amova_rand_all %>% rename(Model = model, Source = component),
  by_cols = c("Model", "Source"),
  value_cols = setdiff(names(amova_rand_all), c("Dataset", "model", "component"))
)

div_overall_side_by_side <- compare_two_dataset_table(
  div_overall_all,
  by_cols = character(0),
  value_cols = c("N", "N_loci", "Ho", "He", "FIS")
)

div_site_side_by_side <- compare_two_dataset_table(
  div_by_site_all,
  by_cols = c("Site"),
  value_cols = c("N", "Ho", "He", "FIS")
)

ar_side_by_side <- compare_two_dataset_table(
  ar_all,
  by_cols = c("Site"),
  value_cols = c("Allelic_Richness", "Allelic_Richness_SE")
)

diff_summary_side_by_side <- compare_two_dataset_table(
  diff_summary_all,
  by_cols = character(0),
  value_cols = c("mean_pairwise_JostD", "median_pairwise_JostD", "mean_pairwise_FST", "median_pairwise_FST")
)

alleles_side_by_side <- compare_two_dataset_table(
  allele_counts_all,
  by_cols = c("Locus"),
  value_cols = c("N_alleles")
)

pca_side_by_side <- compare_two_dataset_table(
  pca_var_all,
  by_cols = c("PC"),
  value_cols = c("Percent_Variance")
)

# Site-level diversity summary aligned with existing output structure
site_summary_full <- full_div$by_site %>%
  left_join(full_div$allelic_richness %>% select(Site, Allelic_Richness, Allelic_Richness_SE), by = "Site") %>%
  mutate(Dataset = "FULL")

site_summary_red <- red_div$by_site %>%
  left_join(red_div$allelic_richness %>% select(Site, Allelic_Richness, Allelic_Richness_SE), by = "Site") %>%
  mutate(Dataset = "REDUCED")

site_genetic_summary_comparison <- bind_rows(site_summary_full, site_summary_red) %>%
  arrange(Site, Dataset)

site_genetic_summary_side_by_side <- compare_two_dataset_table(
  site_genetic_summary_comparison,
  by_cols = c("Site"),
  value_cols = c("N", "Ho", "He", "FIS", "Allelic_Richness", "Allelic_Richness_SE")
)

# ------------------------------------------------------------------
# Write CSV outputs
# ------------------------------------------------------------------
write.csv(amova_results_all, file.path(RESULTS_DIR, "amova_results_full_vs_reduced_long.csv"), row.names = FALSE)
write.csv(amova_rand_all, file.path(RESULTS_DIR, "amova_randtest_full_vs_reduced_long.csv"), row.names = FALSE)
write.csv(amova_side_by_side, file.path(RESULTS_DIR, "amova_comparison.csv"), row.names = FALSE)
write.csv(amova_rand_side_by_side, file.path(RESULTS_DIR, "amova_randtest_comparison.csv"), row.names = FALSE)

write.csv(div_overall_all, file.path(RESULTS_DIR, "heterozygosity_fis_overall_full_vs_reduced_long.csv"), row.names = FALSE)
write.csv(div_by_site_all, file.path(RESULTS_DIR, "heterozygosity_fis_by_site_full_vs_reduced_long.csv"), row.names = FALSE)
write.csv(div_overall_side_by_side, file.path(RESULTS_DIR, "heterozygosity_fis_overall_comparison.csv"), row.names = FALSE)
write.csv(div_site_side_by_side, file.path(RESULTS_DIR, "heterozygosity_fis_by_site_comparison.csv"), row.names = FALSE)

write.csv(ar_all, file.path(RESULTS_DIR, "allelic_richness_full_vs_reduced_long.csv"), row.names = FALSE)
write.csv(ar_side_by_side, file.path(RESULTS_DIR, "allelic_richness_comparison.csv"), row.names = FALSE)

write.csv(diff_summary_all, file.path(RESULTS_DIR, "differentiation_summary_full_vs_reduced_long.csv"), row.names = FALSE)
write.csv(diff_summary_side_by_side, file.path(RESULTS_DIR, "differentiation_summary_comparison.csv"), row.names = FALSE)
write.csv(fst_long_all, file.path(RESULTS_DIR, "pairwise_fst_full_vs_reduced_long.csv"), row.names = FALSE)
write.csv(jost_long_all, file.path(RESULTS_DIR, "pairwise_jostD_full_vs_reduced_long.csv"), row.names = FALSE)

write.csv(allele_counts_all, file.path(RESULTS_DIR, "n_alleles_by_locus_full_vs_reduced_long.csv"), row.names = FALSE)
write.csv(alleles_side_by_side, file.path(RESULTS_DIR, "n_alleles_by_locus_comparison.csv"), row.names = FALSE)

write.csv(pca_var_all, file.path(RESULTS_DIR, "pca_variance_full_vs_reduced_long.csv"), row.names = FALSE)
write.csv(pca_side_by_side, file.path(RESULTS_DIR, "pca_variance_comparison.csv"), row.names = FALSE)

write.csv(site_genetic_summary_comparison, file.path(RESULTS_DIR, "site_genetic_summary_full_vs_reduced_long.csv"), row.names = FALSE)
write.csv(site_genetic_summary_side_by_side, file.path(RESULTS_DIR, "site_genetic_summary_comparison.csv"), row.names = FALSE)

# ------------------------------------------------------------------
# Simple figure: deltas for overall Ho/He/FIS + differentiation means
# ------------------------------------------------------------------
overall_delta_for_plot <- bind_rows(
  div_overall_side_by_side %>% transmute(metric = "Ho", delta = delta_Ho_REDUCED_minus_FULL),
  div_overall_side_by_side %>% transmute(metric = "He", delta = delta_He_REDUCED_minus_FULL),
  div_overall_side_by_side %>% transmute(metric = "FIS", delta = delta_FIS_REDUCED_minus_FULL),
  diff_summary_side_by_side %>% transmute(metric = "mean_pairwise_FST", delta = delta_mean_pairwise_FST_REDUCED_minus_FULL),
  diff_summary_side_by_side %>% transmute(metric = "mean_pairwise_JostD", delta = delta_mean_pairwise_JostD_REDUCED_minus_FULL)
)

delta_plot <- ggplot(overall_delta_for_plot, aes(x = metric, y = delta, fill = metric)) +
  geom_col(width = 0.7, alpha = 0.9, show.legend = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey35") +
  theme_bw(base_size = 11) +
  labs(
    title = "Sensitivity deltas (REDUCED - FULL)",
    subtitle = "Values near zero suggest limited effect of excluding HWE-deviating loci",
    x = NULL,
    y = "Delta"
  )

delta_plot_file <- file.path(FIGURES_HWE_DIR, "overall_metric_deltas.jpeg")
ggsave(delta_plot_file, delta_plot, width = 8, height = 4.5, dpi = 320)

# ------------------------------------------------------------------
# Markdown summary with cautious interpretation language
# ------------------------------------------------------------------
safe_delta <- function(x, digits = 4) {
  ifelse(is.na(x), "NA", format(round(x, digits), nsmall = digits))
}

summary_md <- c(
  "# HWE sensitivity analysis summary",
  "",
  "This analysis compares the original clone-corrected dataset (FULL; all loci) with a reduced dataset (REDUCED; four globally/mostly-site HWE-deviating loci removed).",
  "",
  "## Caution",
  "- This is a **sensitivity analysis** for robustness of inference.",
  "- It does **not** prove null alleles.",
  "- It does **not** prove a Wahlund effect.",
  "",
  "## Key tendencies (REDUCED - FULL)",
  paste0("- Overall Ho delta: ", safe_delta(div_overall_side_by_side$delta_Ho_REDUCED_minus_FULL[1])),
  paste0("- Overall He delta: ", safe_delta(div_overall_side_by_side$delta_He_REDUCED_minus_FULL[1])),
  paste0("- Overall FIS delta: ", safe_delta(div_overall_side_by_side$delta_FIS_REDUCED_minus_FULL[1])),
  paste0("- Mean pairwise FST delta: ", safe_delta(diff_summary_side_by_side$delta_mean_pairwise_FST_REDUCED_minus_FULL[1])),
  paste0("- Mean pairwise Jost's D delta: ", safe_delta(diff_summary_side_by_side$delta_mean_pairwise_JostD_REDUCED_minus_FULL[1])),
  "",
  "## AMOVA comparison",
  "- AMOVA was rerun with the same model structure as the main workflow (Site-only + derived North/South/Site hierarchical model when eligible).",
  "- Side-by-side AMOVA table: `results/hwe_sensitivity/amova_comparison.csv`.",
  "- Side-by-side AMOVA permutation table: `results/hwe_sensitivity/amova_randtest_comparison.csv`.",
  "",
  "## Robustness interpretation",
  "If the direction and broad magnitude of key metrics (AMOVA structure, FST/Jost's D, Ho/He/FIS, and site-level summaries) remain similar between FULL and REDUCED, then your biological interpretation is likely robust to exclusion of these potentially problematic loci. If shifts are substantial, conclusions should be framed with additional caution and potential locus-driven sensitivity should be acknowledged explicitly."
)

summary_file <- file.path(RESULTS_DIR, "hwe_sensitivity_summary.md")
writeLines(summary_md, summary_file)

message("[hwe_sensitivity] Saved summary: ", summary_file)
message("[hwe_sensitivity] Completed successfully.")