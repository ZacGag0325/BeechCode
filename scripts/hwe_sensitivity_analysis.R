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
  library(rlang)
})

source("scripts/_load_objects.R")

# ------------------------------------------------------------------
# User-editable sensitivity settings
# ------------------------------------------------------------------
suspect_null_loci <- c("EJV8T_A_0", "ERHBI_A_0", "FCM5", "FG5")
analysis_suffix <- "noSuspectLoci"
AMOVA_PERMUTATIONS <- 999
HWE_MONTE_CARLO_REPS <- 9999L
BRUVO_MLL_THRESHOLD <- 0.09
BRUVO_ALGORITHM <- "farthest_neighbor"

SENS_TABLES_DIR <- file.path(TABLES_DIR, analysis_suffix)
SENS_FIGURES_DIR <- file.path(FIGURES_DIR, analysis_suffix)
SENS_MATRICES_DIR <- file.path(MATRICES_DIR, analysis_suffix)
COMPARISON_DIR <- file.path(TABLES_DIR, "comparisons")
WORD_DIR <- file.path(OUTPUT_DIR, "word")

for (d in c(SENS_TABLES_DIR, SENS_FIGURES_DIR, SENS_MATRICES_DIR, COMPARISON_DIR, WORD_DIR)) {
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

write_failure_log <- function(step_name, dataset_label, model_label = NA_character_, err_msg = NA_character_) {
  log_df <- data.frame(
    timestamp_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    step = as.character(step_name),
    dataset = as.character(dataset_label),
    model = as.character(model_label),
    error_message = as.character(err_msg),
    stringsAsFactors = FALSE
  )
  log_path <- file.path(SENS_TABLES_DIR, paste0("failure_log_", analysis_suffix, ".csv"))
  if (file.exists(log_path)) {
    existing <- tryCatch(read.csv(log_path, stringsAsFactors = FALSE), error = function(e) NULL)
    if (!is.null(existing)) log_df <- bind_rows(existing, log_df)
  }
  write.csv(log_df, log_path, row.names = FALSE)
  message("[hwe_sensitivity] Wrote failure log: ", log_path)
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
      type = if (!is.null(gobj@type) && length(gobj@type) > 0) gobj@type else "codom",
      NA.char = "NA"
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
    tryCatch({
      adegenet::strata(gi_obj) <- strata_df
      fit <- poppr::poppr.amova(gi_obj, formula_obj)
      
      components <- as.data.frame(fit$componentsofcovariance, stringsAsFactors = FALSE)
      components$Source <- rownames(components)
      rownames(components) <- NULL
      if (ncol(components) > 0) names(components)[1] <- "Sigma"
      
      phi_stats <- as.data.frame(fit$statphi, stringsAsFactors = FALSE)
      phi_stats$Source <- rownames(phi_stats)
      rownames(phi_stats) <- NULL
      if (ncol(phi_stats) > 0) names(phi_stats)[1] <- "Phi"
      
      results <- full_join(components, phi_stats, by = "Source") %>%
        mutate(
          Dataset = dataset_label,
          Model = model_label,
          N_individuals_used = adegenet::nInd(gi_obj),
          N_groups_used = dplyr::n_distinct(strata_df$Site),
          Permutations = AMOVA_PERMUTATIONS
        ) %>%
        select(Dataset, Model, Source, everything())
      
      rand <- tryCatch(
        ade4::randtest(fit, nrepet = AMOVA_PERMUTATIONS),
        error = function(e) {
          cat(
            "[hwe_sensitivity][WARN] randtest failed for dataset=", dataset_label,
            " model=", model_label, " : ", conditionMessage(e), "\n", sep = ""
          )
          write_failure_log("amova_randtest", dataset_label, model_label, conditionMessage(e))
          NULL
        }
      )
      
      cat("[hwe_sensitivity][DEBUG] dataset=", dataset_label, " model=", model_label, " class(rand)=", paste(class(rand), collapse = ","), "\n", sep = "")
      cat("[hwe_sensitivity][DEBUG] dataset=", dataset_label, " model=", model_label, " names(rand)=", if (is.null(rand)) "NULL" else paste(names(rand), collapse = ","), "\n", sep = "")
      cat("[hwe_sensitivity][DEBUG] dataset=", dataset_label, " model=", model_label, " length(rand$obs)=", if (is.null(rand) || is.null(rand$obs)) 0L else length(rand$obs), "\n", sep = "")
      cat("[hwe_sensitivity][DEBUG] dataset=", dataset_label, " model=", model_label, " length(rand$pvalue)=", if (is.null(rand) || is.null(rand$pvalue)) 0L else length(rand$pvalue), "\n", sep = "")
      cat("[hwe_sensitivity][DEBUG] dataset=", dataset_label, " model=", model_label, " names(rand$obs)=", if (is.null(rand) || is.null(rand$obs) || is.null(names(rand$obs))) "NULL" else paste(names(rand$obs), collapse = ","), "\n", sep = "")
      
      variance_sources <- if ("Source" %in% names(results)) as.character(results$Source) else character(0)
      obs_vals <- if (!is.null(rand) && !is.null(rand$obs)) as.numeric(rand$obs) else numeric(0)
      p_vals <- if (!is.null(rand) && !is.null(rand$pvalue)) as.numeric(rand$pvalue) else numeric(0)
      obs_names <- if (!is.null(rand) && !is.null(rand$obs)) names(rand$obs) else NULL
      obs_names <- if (is.null(obs_names)) character(0) else as.character(obs_names)
      
      base_components <- unique(c(obs_names, variance_sources))
      if (length(base_components) == 0) base_components <- "AMOVA_component_unavailable"
      
      rand_df <- data.frame(
        Dataset = dataset_label,
        Model = model_label,
        component = base_components,
        statistic = NA_real_,
        p_value = NA_real_,
        permutations = AMOVA_PERMUTATIONS,
        stringsAsFactors = FALSE
      )
      
      if (length(obs_names) == 0 || length(obs_vals) == 0) {
        cat(
          "[hwe_sensitivity][WARN] rand$obs missing/empty; writing NA statistics for dataset=",
          dataset_label, " model=", model_label, ".\n", sep = ""
        )
      } else {
        n_assign <- min(length(obs_names), length(obs_vals))
        rand_df$statistic[match(obs_names[seq_len(n_assign)], rand_df$component)] <- obs_vals[seq_len(n_assign)]
      }
      
      if (length(obs_names) == 0 || length(p_vals) == 0) {
        cat(
          "[hwe_sensitivity][WARN] rand$pvalue missing/empty; writing NA p-values for dataset=",
          dataset_label, " model=", model_label, ".\n", sep = ""
        )
      } else {
        n_assign_p <- min(length(obs_names), length(p_vals))
        rand_df$p_value[match(obs_names[seq_len(n_assign_p)], rand_df$component)] <- p_vals[seq_len(n_assign_p)]
      }
      
      list(results = results, rand = rand_df)
    }, error = function(e) {
      cat(
        "[hwe_sensitivity][WARN] AMOVA model failed for dataset=", dataset_label,
        " model=", model_label, " : ", conditionMessage(e), "\n", sep = ""
      )
      write_failure_log("amova_model", dataset_label, model_label, conditionMessage(e))
      fail_results <- data.frame(
        Dataset = dataset_label,
        Model = model_label,
        Source = "AMOVA_failed",
        Sigma = NA_real_,
        Phi = NA_real_,
        N_individuals_used = adegenet::nInd(gi_obj),
        N_groups_used = dplyr::n_distinct(strata_df$Site),
        Permutations = AMOVA_PERMUTATIONS,
        stringsAsFactors = FALSE
      )
      fail_rand <- data.frame(
        Dataset = dataset_label,
        Model = model_label,
        component = "AMOVA_failed",
        statistic = NA_real_,
        p_value = NA_real_,
        permutations = AMOVA_PERMUTATIONS,
        stringsAsFactors = FALSE
      )
      list(results = fail_results, rand = fail_rand)
    })
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
# Clonality sensitivity helpers (full gi vs gi with suspect loci removed)
# ------------------------------------------------------------------
calc_clonal_richness <- function(N, G) {
  ifelse(!is.na(N) & !is.na(G) & N > 1, (G - 1) / (N - 1), NA_real_)
}

compute_mlg_mll_from_gi <- function(gobj, threshold = BRUVO_MLL_THRESHOLD, algorithm = BRUVO_ALGORITHM) {
  if (!inherits(gobj, "genind")) {
    stop("[hwe_sensitivity] Clonality input must be a genind object.")
  }
  if (adegenet::nInd(gobj) < 1) {
    stop("[hwe_sensitivity] Clonality input contains no individuals.")
  }
  if (adegenet::nLoc(gobj) < 2) {
    stop("[hwe_sensitivity] Bruvo MLL clustering requires at least two loci.")
  }
  
  gc_mlg <- poppr::as.genclone(gobj)
  mlg_raw <- tryCatch(
    poppr::mlg.vector(gc_mlg),
    error = function(e) as.integer(factor(poppr::mlg(gc_mlg)))
  )
  mlg_labels <- paste0("MLG_", as.integer(factor(mlg_raw)))
  
  replen <- rep(2, adegenet::nLoc(gobj))
  names(replen) <- adegenet::locNames(gobj)
  
  gc_mll <- gc_mlg
  poppr::mlg.filter(
    gc_mll,
    distance = poppr::bruvo.dist,
    replen = replen,
    algorithm = algorithm
  ) <- threshold
  
  mll_raw <- poppr::mll(gc_mll)
  mll_labels <- paste0("MLL_", as.integer(factor(mll_raw)))
  
  if (length(mlg_labels) != adegenet::nInd(gobj) || length(mll_labels) != adegenet::nInd(gobj)) {
    stop("[hwe_sensitivity] MLG/MLL assignment length does not match the number of genind individuals.")
  }
  
  data.frame(
    Individual = adegenet::indNames(gobj),
    Site = as.character(adegenet::pop(gobj)),
    MLG = mlg_labels,
    MLL = mll_labels,
    stringsAsFactors = FALSE
  )
}

add_site_labels_for_clonality <- function(assignments_df, df_ids_tbl) {
  if (!all(c("Individual", "Site") %in% names(assignments_df))) {
    stop("[hwe_sensitivity] Internal clonality assignments are missing Individual or Site columns.")
  }
  aligned_ids <- align_df_ids_to_genind(gi, df_ids_tbl, context = "[hwe_sensitivity] clonality metadata")
  meta_tbl <- data.frame(
    Individual = aligned_ids$ind_id,
    Site = aligned_ids$Site,
    stringsAsFactors = FALSE
  )
  
  optional_cols <- intersect(c("Site_label", "Region", "Site_order"), names(aligned_ids))
  if (length(optional_cols) > 0) {
    meta_tbl <- cbind(meta_tbl, aligned_ids[, optional_cols, drop = FALSE])
  }
  
  out <- assignments_df %>%
    select(Individual, MLG, MLL) %>%
    left_join(meta_tbl, by = "Individual")
  
  if (!"Site_label" %in% names(out)) out$Site_label <- out$Site
  if (!"Site_order" %in% names(out)) out$Site_order <- seq_len(nrow(out))
  if (!"Region" %in% names(out)) out$Region <- NA_character_
  
  out$Site <- ifelse(is.na(out$Site) | !nzchar(out$Site), as.character(assignments_df$Site), as.character(out$Site))
  out$Site_label <- ifelse(is.na(out$Site_label) | !nzchar(as.character(out$Site_label)), out$Site, as.character(out$Site_label))
  out$Site_order <- suppressWarnings(as.numeric(out$Site_order))
  out
}

make_repeated_clone_table <- function(assignments_df, clone_col) {
  clone_sym <- rlang::sym(clone_col)
  assignments_df %>%
    filter(!is.na(!!clone_sym)) %>%
    group_by(!!clone_sym) %>%
    mutate(Group_Size = dplyr::n()) %>%
    ungroup() %>%
    filter(Group_Size > 1) %>%
    arrange(!!clone_sym, Individual)
}

make_repeated_group_signature <- function(repeated_df, clone_col) {
  if (nrow(repeated_df) == 0) return(character(0))
  split(repeated_df$Individual, repeated_df[[clone_col]]) %>%
    lapply(function(x) paste(sort(unique(as.character(x))), collapse = "|")) %>%
    unlist(use.names = FALSE) %>%
    sort()
}

summarise_clonality_by_site <- function(assignments_df, dataset_label) {
  assignments_df %>%
    group_by(Site, Site_label, Site_order) %>%
    summarise(
      Dataset = dataset_label,
      N = dplyr::n(),
      MLG = dplyr::n_distinct(MLG, na.rm = TRUE),
      MLL = dplyr::n_distinct(MLL, na.rm = TRUE),
      repeated_MLG_groups = sum(table(MLG) > 1),
      repeated_MLL_groups = sum(table(MLL) > 1),
      .groups = "drop"
    ) %>%
    mutate(
      R_MLG = calc_clonal_richness(N, MLG),
      R_MLL = calc_clonal_richness(N, MLL)
    ) %>%
    arrange(Site_order, Site_label)
}

summarise_clonality_overall <- function(assignments_df, dataset_label, repeated_mll_same_as_full = NA) {
  repeated_mlg <- make_repeated_clone_table(assignments_df, "MLG")
  repeated_mll <- make_repeated_clone_table(assignments_df, "MLL")
  repeated_mll_sites <- sort(unique(as.character(repeated_mll$Site_label)))
  repeated_mll_inds <- sort(unique(as.character(repeated_mll$Individual)))
  N <- nrow(assignments_df)
  total_MLG <- dplyr::n_distinct(assignments_df$MLG, na.rm = TRUE)
  total_MLL <- dplyr::n_distinct(assignments_df$MLL, na.rm = TRUE)
  
  data.frame(
    Dataset = dataset_label,
    total_N = N,
    total_MLG = total_MLG,
    total_MLL = total_MLL,
    R_MLG = calc_clonal_richness(N, total_MLG),
    R_MLL = calc_clonal_richness(N, total_MLL),
    total_repeated_MLGs = dplyr::n_distinct(repeated_mlg$MLG, na.rm = TRUE),
    total_repeated_MLLs = dplyr::n_distinct(repeated_mll$MLL, na.rm = TRUE),
    total_MLG_clone_copies = N - total_MLG,
    total_MLL_clone_copies = N - total_MLL,
    repeated_MLL_individuals = if (length(repeated_mll_inds) == 0) "none" else paste(repeated_mll_inds, collapse = ";"),
    repeated_MLL_sites = if (length(repeated_mll_sites) == 0) "none" else paste(repeated_mll_sites, collapse = ";"),
    same_individuals_sites_in_repeated_MLLs_as_full = as.character(repeated_mll_same_as_full),
    stringsAsFactors = FALSE
  )
}

build_clonality_sensitivity <- function(full_gobj, reduced_gobj, df_ids_tbl) {
  message("[hwe_sensitivity] Recomputing full-loci MLGs and Bruvo MLLs from gi for clonality sensitivity.")
  message("[hwe_sensitivity] Recomputing reduced-loci MLGs and Bruvo MLLs after removing suspect loci.")
  message("[hwe_sensitivity] Bruvo MLL threshold: ", BRUVO_MLL_THRESHOLD)
  message("[hwe_sensitivity] Bruvo clustering algorithm: ", BRUVO_ALGORITHM)
  
  full_assign <- compute_mlg_mll_from_gi(full_gobj) %>%
    add_site_labels_for_clonality(df_ids_tbl) %>%
    mutate(Dataset = "FULL")
  reduced_assign <- compute_mlg_mll_from_gi(reduced_gobj) %>%
    add_site_labels_for_clonality(df_ids_tbl) %>%
    mutate(Dataset = "REDUCED")
  
  if (!identical(full_assign$Individual, reduced_assign$Individual)) {
    stop("[hwe_sensitivity] Full and reduced clonality assignments are not in the same individual order.")
  }
  
  full_site <- summarise_clonality_by_site(full_assign, "FULL")
  reduced_site <- summarise_clonality_by_site(reduced_assign, "REDUCED")
  
  comparison <- full_site %>%
    select(Site, Site_label, Site_order, N_full = N, MLG_full = MLG, R_MLG_full = R_MLG, MLL_full = MLL, R_MLL_full = R_MLL) %>%
    full_join(
      reduced_site %>%
        select(Site, N_reduced = N, MLG_reduced = MLG, R_MLG_reduced = R_MLG, MLL_reduced = MLL, R_MLL_reduced = R_MLL),
      by = "Site"
    ) %>%
    mutate(
      delta_MLG = MLG_reduced - MLG_full,
      delta_MLL = MLL_reduced - MLL_full,
      delta_R_MLG = R_MLG_reduced - R_MLG_full,
      delta_R_MLL = R_MLL_reduced - R_MLL_full,
      interpretation = dplyr::case_when(
        is.na(delta_MLG) | is.na(delta_MLL) ~ "Could not compare this site because one dataset is missing site-level values.",
        delta_MLG == 0 & delta_MLL == 0 & abs(delta_R_MLG) < 1e-12 & abs(delta_R_MLL) < 1e-12 ~ "No meaningful clonality change after removing HWE-deviating loci.",
        delta_MLL != 0 | abs(delta_R_MLL) >= 0.01 ~ "Meaningful clonality change: Bruvo-based MLL count or R_MLL changed.",
        delta_MLG != 0 | abs(delta_R_MLG) >= 0.01 ~ "Minor clonality change: exact MLG count or R_MLG changed, but MLL conclusion is stable.",
        TRUE ~ "Very small numerical change only; clonality conclusion is stable."
      )
    ) %>%
    arrange(Site_order, Site_label) %>%
    transmute(
      Site = as.character(Site_label),
      N_full,
      MLG_full,
      R_MLG_full,
      MLL_full,
      R_MLL_full,
      N_reduced,
      MLG_reduced,
      R_MLG_reduced,
      MLL_reduced,
      R_MLL_reduced,
      delta_MLG,
      delta_MLL,
      delta_R_MLG,
      delta_R_MLL,
      interpretation
    )
  
  full_mll_repeated <- make_repeated_clone_table(full_assign, "MLL")
  reduced_mll_repeated <- make_repeated_clone_table(reduced_assign, "MLL")
  same_repeated_mll_individuals <- identical(
    sort(unique(as.character(full_mll_repeated$Individual))),
    sort(unique(as.character(reduced_mll_repeated$Individual)))
  )
  same_repeated_mll_sites <- identical(
    sort(unique(as.character(full_mll_repeated$Site_label))),
    sort(unique(as.character(reduced_mll_repeated$Site_label)))
  )
  same_repeated_mll_groups <- identical(
    make_repeated_group_signature(full_mll_repeated, "MLL"),
    make_repeated_group_signature(reduced_mll_repeated, "MLL")
  )
  same_individuals_sites <- same_repeated_mll_individuals && same_repeated_mll_sites
  
  overall <- bind_rows(
    summarise_clonality_overall(full_assign, "FULL", repeated_mll_same_as_full = TRUE),
    summarise_clonality_overall(reduced_assign, "REDUCED", repeated_mll_same_as_full = same_individuals_sites)
  )
  overall$same_repeated_MLL_group_memberships_as_full <- as.character(c(TRUE, same_repeated_mll_groups))
  
  changed_sites <- comparison %>%
    filter(
      !is.na(delta_MLG) & !is.na(delta_MLL) &
        (delta_MLG != 0 | delta_MLL != 0 | abs(delta_R_MLG) >= 0.01 | abs(delta_R_MLL) >= 0.01)
    ) %>%
    pull(Site)
  
  list(
    comparison = comparison,
    overall = overall,
    assignments = bind_rows(full_assign, reduced_assign),
    repeated_MLG = bind_rows(
      make_repeated_clone_table(full_assign, "MLG") %>% mutate(Dataset = "FULL"),
      make_repeated_clone_table(reduced_assign, "MLG") %>% mutate(Dataset = "REDUCED")
    ),
    repeated_MLL = bind_rows(
      full_mll_repeated %>% mutate(Dataset = "FULL"),
      reduced_mll_repeated %>% mutate(Dataset = "REDUCED")
    ),
    changed_sites = changed_sites,
    same_repeated_mll_individuals = same_repeated_mll_individuals,
    same_repeated_mll_sites = same_repeated_mll_sites,
    same_repeated_mll_groups = same_repeated_mll_groups
  )
}

xml_escape_word <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub('"', "&quot;", x, fixed = TRUE)
  x <- gsub("'", "&apos;", x, fixed = TRUE)
  x
}

word_cell_xml <- function(value, bold = FALSE) {
  bold_xml <- if (bold) "<w:rPr><w:b/></w:rPr>" else ""
  paste0(
    "<w:tc>",
    "<w:tcPr><w:tcW w:w=\"2400\" w:type=\"dxa\"/></w:tcPr>",
    "<w:p><w:r>", bold_xml, "<w:t>", xml_escape_word(value), "</w:t></w:r></w:p>",
    "</w:tc>"
  )
}

word_paragraph_xml <- function(value, bold = FALSE) {
  bold_xml <- if (bold) "<w:rPr><w:b/></w:rPr>" else ""
  paste0("<w:p><w:r>", bold_xml, "<w:t>", xml_escape_word(value), "</w:t></w:r></w:p>")
}

word_table_xml <- function(df) {
  df <- as.data.frame(lapply(df, as.character), stringsAsFactors = FALSE, check.names = FALSE)
  header <- paste0("<w:tr>", paste(vapply(names(df), word_cell_xml, character(1), bold = TRUE), collapse = ""), "</w:tr>")
  rows <- apply(df, 1, function(row) paste0("<w:tr>", paste(vapply(row, word_cell_xml, character(1)), collapse = ""), "</w:tr>"))
  paste0(
    "<w:tbl>",
    "<w:tblPr><w:tblStyle w:val=\"TableGrid\"/><w:tblW w:w=\"0\" w:type=\"auto\"/>",
    "<w:tblBorders>",
    "<w:top w:val=\"single\" w:sz=\"4\" w:space=\"0\" w:color=\"auto\"/>",
    "<w:left w:val=\"single\" w:sz=\"4\" w:space=\"0\" w:color=\"auto\"/>",
    "<w:bottom w:val=\"single\" w:sz=\"4\" w:space=\"0\" w:color=\"auto\"/>",
    "<w:right w:val=\"single\" w:sz=\"4\" w:space=\"0\" w:color=\"auto\"/>",
    "<w:insideH w:val=\"single\" w:sz=\"4\" w:space=\"0\" w:color=\"auto\"/>",
    "<w:insideV w:val=\"single\" w:sz=\"4\" w:space=\"0\" w:color=\"auto\"/>",
    "</w:tblBorders></w:tblPr>",
    header,
    paste(rows, collapse = ""),
    "</w:tbl>"
  )
}

write_clonality_docx <- function(comparison_tbl, overall_tbl, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp_dir <- tempfile("hwe_clonality_docx_")
  dir.create(file.path(tmp_dir, "_rels"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(tmp_dir, "word", "_rels"), recursive = TRUE, showWarnings = FALSE)
  
  document_xml <- paste0(
    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
    '<w:document xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main" ',
    'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">',
    '<w:body>',
    word_paragraph_xml("HWE sensitivity clonality comparison", bold = TRUE),
    word_paragraph_xml("Site-level FULL vs REDUCED comparison", bold = TRUE),
    word_table_xml(comparison_tbl),
    word_paragraph_xml("Overall summary", bold = TRUE),
    word_table_xml(overall_tbl),
    '<w:sectPr><w:pgSz w:w="15840" w:h="12240" w:orient="landscape"/>',
    '<w:pgMar w:top="720" w:right="720" w:bottom="720" w:left="720" w:header="360" w:footer="360" w:gutter="0"/></w:sectPr>',
    '</w:body></w:document>'
  )
  
  writeLines(c(
    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
    '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">',
    '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>',
    '<Default Extension="xml" ContentType="application/xml"/>',
    '<Override PartName="/word/document.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.document.main+xml"/>',
    '</Types>'
  ), file.path(tmp_dir, "[Content_Types].xml"), useBytes = TRUE)
  writeLines(c(
    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
    '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">',
    '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="word/document.xml"/>',
    '</Relationships>'
  ), file.path(tmp_dir, "_rels", ".rels"), useBytes = TRUE)
  writeLines(c(
    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
    '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"></Relationships>'
  ), file.path(tmp_dir, "word", "_rels", "document.xml.rels"), useBytes = TRUE)
  writeLines(document_xml, file.path(tmp_dir, "word", "document.xml"), useBytes = TRUE)
  
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(tmp_dir)
  if (file.exists(path)) unlink(path)
  utils::zip(
    zipfile = path,
    files = list.files(tmp_dir, recursive = TRUE, all.files = TRUE, no.. = TRUE),
    flags = "-q"
  )
  if (!file.exists(path)) stop("[hwe_sensitivity] Failed to create Word document at: ", path)
  message("[hwe_sensitivity] Saved: ", path)
  invisible(path)
}

# ------------------------------------------------------------------
# Build reduced object and manifest
# ------------------------------------------------------------------
full_gi <- gi_mll
reduced_info <- make_reduced_genind(full_gi, suspect_null_loci)
reduced_gi <- reduced_info$reduced

full_clonality_gi <- gi
reduced_clonality_info <- make_reduced_genind(full_clonality_gi, suspect_null_loci)
reduced_clonality_gi <- reduced_clonality_info$reduced

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

clonality_locus_manifest <- data.frame(
  dataset = c("FULL_CLONALITY_GI", "REDUCED_CLONALITY_GI"),
  n_loci = c(adegenet::nLoc(full_clonality_gi), adegenet::nLoc(reduced_clonality_gi)),
  loci = c(
    paste(adegenet::locNames(full_clonality_gi), collapse = ";"),
    paste(adegenet::locNames(reduced_clonality_gi), collapse = ";")
  ),
  removed_loci = c("", paste(reduced_clonality_info$found, collapse = ";")),
  requested_but_not_found = c("", paste(reduced_clonality_info$missing, collapse = ";")),
  stringsAsFactors = FALSE
)
write_csv_msg(clonality_locus_manifest, file.path(SENS_TABLES_DIR, paste0("locus_manifest_clonality_", analysis_suffix, ".csv")))

# ------------------------------------------------------------------
# Run parallel full/reduced analyses
# ------------------------------------------------------------------
full_div <- run_diversity_bundle(full_gi, "FULL")
red_div <- run_diversity_bundle(reduced_gi, "REDUCED")

full_diff <- run_differentiation_bundle(full_gi, "FULL")
red_diff <- run_differentiation_bundle(reduced_gi, "REDUCED")

full_amova <- tryCatch(
  run_amova_bundle(full_gi, "FULL"),
  error = function(e) {
    write_failure_log("run_amova_bundle", "FULL", "bundle", conditionMessage(e))
    list(
      results = data.frame(
        Dataset = "FULL",
        Model = "bundle_failed",
        Source = "bundle_failed",
        Sigma = NA_real_,
        Phi = NA_real_,
        N_individuals_used = nInd(full_gi),
        N_groups_used = NA_integer_,
        Permutations = AMOVA_PERMUTATIONS,
        stringsAsFactors = FALSE
      ),
      rand = data.frame(
        Dataset = "FULL",
        Model = "bundle_failed",
        component = "bundle_failed",
        statistic = NA_real_,
        p_value = NA_real_,
        permutations = AMOVA_PERMUTATIONS,
        stringsAsFactors = FALSE
      )
    )
  }
)
red_amova <- tryCatch(
  run_amova_bundle(reduced_gi, "REDUCED"),
  error = function(e) {
    write_failure_log("run_amova_bundle", "REDUCED", "bundle", conditionMessage(e))
    list(
      results = data.frame(
        Dataset = "REDUCED",
        Model = "bundle_failed",
        Source = "bundle_failed",
        Sigma = NA_real_,
        Phi = NA_real_,
        N_individuals_used = nInd(reduced_gi),
        N_groups_used = NA_integer_,
        Permutations = AMOVA_PERMUTATIONS,
        stringsAsFactors = FALSE
      ),
      rand = data.frame(
        Dataset = "REDUCED",
        Model = "bundle_failed",
        component = "bundle_failed",
        statistic = NA_real_,
        p_value = NA_real_,
        permutations = AMOVA_PERMUTATIONS,
        stringsAsFactors = FALSE
      )
    )
  }
)

full_hwe <- run_hwe_by_site_locus(full_gi, "FULL")
red_hwe <- run_hwe_by_site_locus(reduced_gi, "REDUCED")

full_pca <- run_pca_bundle(full_gi, "FULL")
red_pca <- run_pca_bundle(reduced_gi, "REDUCED")

full_dapc <- run_dapc_bundle(full_gi, "FULL")
red_dapc <- run_dapc_bundle(reduced_gi, "REDUCED")

clonality_sensitivity <- build_clonality_sensitivity(full_clonality_gi, reduced_clonality_gi, df_ids)

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

clonality_output <- list(
  site_comparison = clonality_sensitivity$comparison,
  overall_summary = clonality_sensitivity$overall,
  individual_assignments = clonality_sensitivity$assignments,
  repeated_MLGs = clonality_sensitivity$repeated_MLG,
  repeated_MLLs = clonality_sensitivity$repeated_MLL,
  settings = list(
    suspect_loci_requested = suspect_null_loci,
    suspect_loci_removed_for_clonality = reduced_clonality_info$found,
    suspect_loci_not_found_for_clonality = reduced_clonality_info$missing,
    bruvo_mll_threshold = BRUVO_MLL_THRESHOLD,
    bruvo_algorithm = BRUVO_ALGORITHM
  )
)
write_csv_msg(clonality_sensitivity$comparison, file.path(TABLES_DIR, "hwe_sensitivity_clonality_comparison.csv"))
saveRDS(clonality_output, file.path(TABLES_DIR, "hwe_sensitivity_clonality_comparison.rds"))
message("[hwe_sensitivity] Saved: ", file.path(TABLES_DIR, "hwe_sensitivity_clonality_comparison.rds"))
write_clonality_docx(
  clonality_sensitivity$comparison,
  clonality_sensitivity$overall,
  file.path(WORD_DIR, "hwe_sensitivity_clonality_comparison.docx")
)
write_csv_msg(clonality_sensitivity$overall, file.path(TABLES_DIR, "hwe_sensitivity_clonality_overall_summary.csv"))
write_csv_msg(clonality_sensitivity$assignments, file.path(TABLES_DIR, "hwe_sensitivity_clonality_individual_assignments_long.csv"))
write_csv_msg(clonality_sensitivity$repeated_MLL, file.path(TABLES_DIR, "hwe_sensitivity_clonality_repeated_MLLs_long.csv"))

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

meaningful_site_changes <- clonality_sensitivity$changed_sites
overall_full <- clonality_sensitivity$overall[clonality_sensitivity$overall$Dataset == "FULL", , drop = FALSE]
overall_reduced <- clonality_sensitivity$overall[clonality_sensitivity$overall$Dataset == "REDUCED", , drop = FALSE]
overall_changed <- !identical(overall_full$total_MLG, overall_reduced$total_MLG) ||
  !identical(overall_full$total_MLL, overall_reduced$total_MLL) ||
  !isTRUE(all.equal(overall_full$R_MLG, overall_reduced$R_MLG, tolerance = 1e-12)) ||
  !isTRUE(all.equal(overall_full$R_MLL, overall_reduced$R_MLL, tolerance = 1e-12)) ||
  !isTRUE(clonality_sensitivity$same_repeated_mll_individuals) ||
  !isTRUE(clonality_sensitivity$same_repeated_mll_sites)

cat("\n")
cat("============================================================\n")
cat("HWE SENSITIVITY CLONALITY RESULT\n")
cat("============================================================\n")
cat("Removed loci: ", if (length(reduced_clonality_info$found) > 0) paste(reduced_clonality_info$found, collapse = ", ") else "none", "\n", sep = "")
cat("Bruvo MLL threshold: ", BRUVO_MLL_THRESHOLD, "\n", sep = "")
cat("Bruvo clustering algorithm: ", BRUVO_ALGORITHM, "\n", sep = "")
cat("Did clonality change after removing HWE-deviating loci? ", if (overall_changed) "YES" else "NO", "\n", sep = "")
cat("Sites changed: ", if (length(meaningful_site_changes) > 0) paste(meaningful_site_changes, collapse = ", ") else "none", "\n", sep = "")
cat("Same individuals in repeated MLLs? ", if (clonality_sensitivity$same_repeated_mll_individuals) "YES" else "NO", "\n", sep = "")
cat("Same sites in repeated MLLs? ", if (clonality_sensitivity$same_repeated_mll_sites) "YES" else "NO", "\n", sep = "")
cat("Main clonality conclusion robust? ", if (!overall_changed && length(meaningful_site_changes) == 0) "YES" else "CHECK SITE-LEVEL CHANGES", "\n", sep = "")
cat("Primary outputs:\n")
cat("  - ", file.path(TABLES_DIR, "hwe_sensitivity_clonality_comparison.csv"), "\n", sep = "")
cat("  - ", file.path(TABLES_DIR, "hwe_sensitivity_clonality_comparison.rds"), "\n", sep = "")
cat("  - ", file.path(WORD_DIR, "hwe_sensitivity_clonality_comparison.docx"), "\n", sep = "")
cat("============================================================\n\n")

message("[hwe_sensitivity] Completed reduced-loci sensitivity branch successfully.")