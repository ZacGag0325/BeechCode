# scripts/04_amova.R
############################################################
# AMOVA (clone-corrected: gi_mll)
# Why gi_mll here?
# - AMOVA is a population-level analysis and assumes independent genotypes.
# - Clone-corrected genotypes reduce bias from repeated ramets.
#
# This script explicitly aligns individuals and grouping labels to avoid
# strata/data mismatches such as:
#   strata: pop
#   Data:   Site
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

# ------------------------------------------------------------
# 1) Resolve and validate grouping variable (Site)
# ------------------------------------------------------------
resolve_col <- function(df, choices) {
  nms <- names(df)
  idx <- match(TRUE, tolower(nms) %in% tolower(choices), nomatch = 0)
  if (idx == 0) return(NA_character_)
  nms[idx]
}

id_col <- resolve_col(df_ids, c("ind", "individual", "sample", "sampleid", "id"))
site_col <- resolve_col(df_ids, c("site", "population", "pop"))

if (is.na(id_col) || is.na(site_col)) {
  stop("[04_amova] df_ids must contain individual ID and Site columns.")
}

id_to_site <- setNames(as.character(df_ids[[site_col]]), as.character(df_ids[[id_col]]))
inds <- adegenet::indNames(gi_mll)
site_from_dfids <- id_to_site[inds]
site_from_pop <- as.character(adegenet::pop(gi_mll))

# Prefer mapped Site from df_ids if complete; fallback to pop(gi_mll)
if (all(!is.na(site_from_dfids))) {
  group_site <- site_from_dfids
  grouping_source <- "df_ids$Site"
} else {
  group_site <- site_from_pop
  grouping_source <- "pop(gi_mll)"
}

valid <- !is.na(group_site) & nzchar(group_site)
n_drop <- sum(!valid)

if (n_drop > 0) {
  message("[04_amova] Dropping ", n_drop, " individuals with missing group labels.")
}

gi_use <- gi_mll[valid, , drop = FALSE]
site_use <- as.factor(group_site[valid])

if (adegenet::nInd(gi_use) < 2) {
  stop("[04_amova] Fewer than 2 individuals remain after grouping alignment.")
}

# Keep only groups with >= 2 individuals for defensible AMOVA variance estimation
tab_groups <- table(site_use)
keep_groups <- names(tab_groups)[tab_groups >= 2]
if (length(keep_groups) < 2) {
  stop("[04_amova] Need at least two groups with >=2 individuals each for AMOVA.")
}

keep_idx <- site_use %in% keep_groups
gi_use <- gi_use[keep_idx, , drop = FALSE]
site_use <- droplevels(site_use[keep_idx])

message("[04_amova] Individuals used: ", adegenet::nInd(gi_use))
message("[04_amova] Grouping variable used: Site (source: ", grouping_source, ")")
message("[04_amova] Groups included: ", nlevels(site_use), " -> ", paste(levels(site_use), collapse = ", "))
message("[04_amova] AMOVA formula: ~pop")

# ------------------------------------------------------------
# 2) Enforce consistent pop + strata naming for poppr.amova
# ------------------------------------------------------------
adegenet::pop(gi_use) <- site_use
strata_df <- data.frame(
  pop = site_use,
  Site = site_use,
  row.names = adegenet::indNames(gi_use),
  stringsAsFactors = TRUE
)
adegenet::strata(gi_use) <- strata_df

if (!"pop" %in% names(adegenet::strata(gi_use))) {
  stop("[04_amova] strata(gi_use) does not contain column 'pop'.")
}
if (length(adegenet::pop(gi_use)) != adegenet::nInd(gi_use)) {
  stop("[04_amova] pop vector length does not match number of individuals.")
}
if (!identical(as.character(adegenet::pop(gi_use)), as.character(adegenet::strata(gi_use)$pop))) {
  stop("[04_amova] pop(gi_use) and strata(gi_use)$pop are not aligned.")
}

# ------------------------------------------------------------
# 3) Run AMOVA and permutation test
# ------------------------------------------------------------
amova_fit <- poppr::poppr.amova(gi_use, ~pop)

# Run randtest robustly and do not force incompatible vector lengths into mutate
randtest_fit <- tryCatch(
  ade4::randtest(amova_fit, nrepet = 999),
  error = function(e) {
    warning("[04_amova] randtest failed: ", conditionMessage(e))
    NULL
  }
)

# Components and Phi table
components <- as.data.frame(amova_fit$componentsofcovariance)
components$Source <- rownames(components)
rownames(components) <- NULL
names(components)[1] <- "Sigma"

phi_stats <- as.data.frame(amova_fit$statphi)
phi_stats$Source <- rownames(phi_stats)
rownames(phi_stats) <- NULL
names(phi_stats)[1] <- "Phi"

amova_results <- full_join(components, phi_stats, by = "Source") %>%
  mutate(
    N_individuals_used = adegenet::nInd(gi_use),
    N_groups_used = nlevels(site_use),
    Grouping_variable = "Site",
    Grouping_source = grouping_source,
    Euclidean_correction_note = "If non-euclidean, poppr::poppr.amova internally applies correction (e.g., quasieuclid)."
  ) %>%
  select(Source, everything())

amova_file <- file.path(TABLES_DIR, "amova_results.csv")
write.csv(amova_results, amova_file, row.names = FALSE)
message("[04_amova] Saved: ", amova_file)

# ------------------------------------------------------------
# 4) Supplementary randtest summary (structurally safe)
# ------------------------------------------------------------
if (is.null(randtest_fit)) {
  amova_randtest_summary <- data.frame(
    test = "AMOVA_randtest",
    component = NA_character_,
    statistic = NA_real_,
    p_value = NA_real_,
    permutations = 999,
    grouping_variable = "Site",
    grouping_source = grouping_source,
    n_individuals_used = adegenet::nInd(gi_use),
    n_groups_used = nlevels(site_use),
    groups_included = paste(levels(site_use), collapse = ";"),
    dropped_for_missing_group = as.integer(n_drop),
    note = "randtest_failed",
    stringsAsFactors = FALSE
  )
} else {
  obs_vec <- as.numeric(randtest_fit$obs)
  p_vec <- as.numeric(randtest_fit$pvalue)
  
  # Harmonize lengths safely
  n_comp <- max(length(obs_vec), length(p_vec), 1)
  if (length(obs_vec) == 0) obs_vec <- rep(NA_real_, n_comp)
  if (length(p_vec) == 0) p_vec <- rep(NA_real_, n_comp)
  if (length(obs_vec) < n_comp) obs_vec <- c(obs_vec, rep(NA_real_, n_comp - length(obs_vec)))
  if (length(p_vec) < n_comp) p_vec <- c(p_vec, rep(NA_real_, n_comp - length(p_vec)))
  
  comp_names <- names(randtest_fit$obs)
  if (is.null(comp_names) || length(comp_names) == 0) comp_names <- paste0("component_", seq_len(n_comp))
  if (length(comp_names) < n_comp) comp_names <- c(comp_names, paste0("component_", (length(comp_names)+1):n_comp))
  
  amova_randtest_summary <- data.frame(
    test = "AMOVA_randtest",
    component = comp_names[seq_len(n_comp)],
    statistic = obs_vec[seq_len(n_comp)],
    p_value = p_vec[seq_len(n_comp)],
    permutations = 999,
    grouping_variable = "Site",
    grouping_source = grouping_source,
    n_individuals_used = adegenet::nInd(gi_use),
    n_groups_used = nlevels(site_use),
    groups_included = paste(levels(site_use), collapse = ";"),
    dropped_for_missing_group = as.integer(n_drop),
    note = "ok",
    stringsAsFactors = FALSE
  )
}

supp_file <- file.path(TABLES_SUPP_DIR, "amova_randtest_summary.csv")
write.csv(amova_randtest_summary, supp_file, row.names = FALSE)
message("[04_amova] Saved: ", supp_file)