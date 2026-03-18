# scripts/04_amova.R
############################################################
# AMOVA (clone-corrected: gi_mll)
# Why gi_mll here?
# - AMOVA is a population-level analysis and assumes independent genotypes.
# - Clone-corrected genotypes reduce bias from repeated ramets.
#
# Main analysis: Site-level AMOVA.
# Optional extension: hierarchical AMOVA (Region/Site) if Region metadata exist.
#
# Outputs:
# - outputs/tables/amova_results.csv
# - outputs/tables/supplementary/amova_randtest_summary.csv
############################################################

suppressPackageStartupMessages({
  library(adegenet)
  library(poppr)
  library(ade4)
  library(dplyr)
})

source("scripts/_load_objects.R")

message("[04_amova] Running AMOVA on gi_mll...")

resolve_region_col <- function(df) {
  resolve_col_ci(df, c("Region", "region", "REGION", "Area", "area", "Zone", "zone"))
}

run_amova_model <- function(gi_use, strata_df, model_formula, model_label, grouping_source, n_drop) {
  adegenet::strata(gi_use) <- strata_df
  amova_fit <- poppr::poppr.amova(gi_use, model_formula)
  
  randtest_fit <- tryCatch(
    ade4::randtest(amova_fit, nrepet = 999),
    error = function(e) {
      warning("[04_amova] randtest failed for ", model_label, ": ", conditionMessage(e))
      NULL
    }
  )
  
  components <- as.data.frame(amova_fit$componentsofcovariance)
  components$Source <- rownames(components)
  rownames(components) <- NULL
  names(components)[1] <- "Sigma"
  
  phi_stats <- as.data.frame(amova_fit$statphi)
  phi_stats$Source <- rownames(phi_stats)
  rownames(phi_stats) <- NULL
  names(phi_stats)[1] <- "Phi"
  
  result <- full_join(components, phi_stats, by = "Source") %>%
    mutate(
      Model = model_label,
      N_individuals_used = adegenet::nInd(gi_use),
      N_groups_used = dplyr::n_distinct(strata_df$Site),
      Grouping_source = grouping_source,
      Euclidean_correction_note = "If non-euclidean, poppr::poppr.amova internally applies correction (e.g., quasieuclid)."
    ) %>%
    select(Model, Source, everything())
  
  if (is.null(randtest_fit)) {
    rand <- data.frame(
      model = model_label,
      test = "AMOVA_randtest",
      component = NA_character_,
      statistic = NA_real_,
      p_value = NA_real_,
      permutations = 999,
      grouping_source = grouping_source,
      n_individuals_used = adegenet::nInd(gi_use),
      n_groups_used = dplyr::n_distinct(strata_df$Site),
      dropped_for_missing_group = as.integer(n_drop),
      note = "randtest_failed",
      stringsAsFactors = FALSE
    )
  } else {
    obs_vec <- as.numeric(randtest_fit$obs)
    p_vec <- as.numeric(randtest_fit$pvalue)
    n_comp <- max(length(obs_vec), length(p_vec), 1)
    if (length(obs_vec) == 0) obs_vec <- rep(NA_real_, n_comp)
    if (length(p_vec) == 0) p_vec <- rep(NA_real_, n_comp)
    if (length(obs_vec) < n_comp) obs_vec <- c(obs_vec, rep(NA_real_, n_comp - length(obs_vec)))
    if (length(p_vec) < n_comp) p_vec <- c(p_vec, rep(NA_real_, n_comp - length(p_vec)))
    
    comp_names <- names(randtest_fit$obs)
    if (is.null(comp_names) || length(comp_names) == 0) comp_names <- paste0("component_", seq_len(n_comp))
    if (length(comp_names) < n_comp) comp_names <- c(comp_names, paste0("component_", (length(comp_names) + 1):n_comp))
    
    rand <- data.frame(
      model = model_label,
      test = "AMOVA_randtest",
      component = comp_names[seq_len(n_comp)],
      statistic = obs_vec[seq_len(n_comp)],
      p_value = p_vec[seq_len(n_comp)],
      permutations = 999,
      grouping_source = grouping_source,
      n_individuals_used = adegenet::nInd(gi_use),
      n_groups_used = dplyr::n_distinct(strata_df$Site),
      dropped_for_missing_group = as.integer(n_drop),
      note = "ok",
      stringsAsFactors = FALSE
    )
  }
  
  list(result = result, rand = rand)
}

# ------------------------------------------------------------
# 1) Resolve and validate grouping variables
# ------------------------------------------------------------
df_ids_cols <- resolve_df_ids_columns(df_ids, context = "[04_amova]", require = TRUE)
id_col <- df_ids_cols$id_col
site_col <- df_ids_cols$site_col
region_col <- resolve_region_col(df_ids)

id_to_site <- setNames(as.character(df_ids[[site_col]]), normalize_id(df_ids[[id_col]]))
id_to_region <- if (!is.na(region_col)) setNames(as.character(df_ids[[region_col]]), normalize_id(df_ids[[id_col]])) else NULL

inds <- adegenet::indNames(gi_mll)
site_from_dfids <- id_to_site[normalize_id(inds)]
site_from_pop <- as.character(adegenet::pop(gi_mll))

if (all(!is.na(site_from_dfids))) {
  group_site <- site_from_dfids
  grouping_source <- "df_ids$Site"
} else {
  group_site <- site_from_pop
  grouping_source <- "pop(gi_mll)"
}

group_region <- if (!is.null(id_to_region)) id_to_region[normalize_id(inds)] else rep(NA_character_, length(inds))

valid <- !is.na(group_site) & nzchar(group_site)
n_drop <- sum(!valid)

if (n_drop > 0) message("[04_amova] Dropping ", n_drop, " individuals with missing Site labels.")

gi_use <- gi_mll[valid, , drop = FALSE]
site_use <- as.factor(group_site[valid])
region_use <- as.factor(group_region[valid])

if (adegenet::nInd(gi_use) < 2) stop("[04_amova] Fewer than 2 individuals remain after grouping alignment.")

tab_groups <- table(site_use)
keep_groups <- names(tab_groups)[tab_groups >= 2]
if (length(keep_groups) < 2) stop("[04_amova] Need at least two sites with >=2 individuals each for AMOVA.")

keep_idx <- site_use %in% keep_groups
gi_use <- gi_use[keep_idx, , drop = FALSE]
site_use <- droplevels(site_use[keep_idx])
region_use <- droplevels(region_use[keep_idx])

adegenet::pop(gi_use) <- site_use

message("[04_amova] Individuals used: ", adegenet::nInd(gi_use))
message("[04_amova] Site groups included: ", nlevels(site_use), " -> ", paste(levels(site_use), collapse = ", "))

# ------------------------------------------------------------
# 2) Run Site-only AMOVA (main)
# ------------------------------------------------------------
strata_site <- data.frame(
  pop = site_use,
  Site = site_use,
  row.names = adegenet::indNames(gi_use),
  stringsAsFactors = TRUE
)

site_fit <- run_amova_model(
  gi_use = gi_use,
  strata_df = strata_site,
  model_formula = ~pop,
  model_label = "Site_only",
  grouping_source = grouping_source,
  n_drop = n_drop
)

amova_results <- site_fit$result
amova_randtest_summary <- site_fit$rand

# ------------------------------------------------------------
# 3) Optional hierarchical AMOVA (Region/Site)
# ------------------------------------------------------------
has_region <- !all(is.na(region_use) | !nzchar(as.character(region_use)))
if (has_region && nlevels(droplevels(region_use)) >= 2) {
  message("[04_amova] Region column detected (", region_col, "); running hierarchical AMOVA: ~Region/Site")
  
  keep_region <- !is.na(region_use) & nzchar(as.character(region_use))
  gi_h <- gi_use[keep_region, , drop = FALSE]
  site_h <- droplevels(site_use[keep_region])
  region_h <- droplevels(region_use[keep_region])
  adegenet::pop(gi_h) <- site_h
  
  strata_h <- data.frame(
    Region = region_h,
    Site = site_h,
    pop = site_h,
    row.names = adegenet::indNames(gi_h),
    stringsAsFactors = TRUE
  )
  
  h_fit <- run_amova_model(
    gi_use = gi_h,
    strata_df = strata_h,
    model_formula = ~Region/Site,
    model_label = "Region_Site_hierarchical",
    grouping_source = paste0(grouping_source, " + df_ids$", region_col),
    n_drop = n_drop
  )
  
  amova_results <- bind_rows(amova_results, h_fit$result)
  amova_randtest_summary <- bind_rows(amova_randtest_summary, h_fit$rand)
} else {
  message("[04_amova] Region variable not available with >=2 levels; keeping Site-only AMOVA.")
}

amova_file <- file.path(TABLES_DIR, "amova_results.csv")
write.csv(amova_results, amova_file, row.names = FALSE)
message("[04_amova] Saved: ", amova_file)

supp_file <- file.path(TABLES_SUPP_DIR, "amova_randtest_summary.csv")
write.csv(amova_randtest_summary, supp_file, row.names = FALSE)
message("[04_amova] Saved: ", supp_file)