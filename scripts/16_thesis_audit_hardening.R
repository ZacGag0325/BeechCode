############################################################
# scripts/16_thesis_audit_hardening.R
# Thesis-defense hardening: QC, differentiation, IBD/IBE, robustness
############################################################

options(repos = c(CRAN = "https://cloud.r-project.org"))
pkgs <- c("adegenet", "poppr", "hierfstat", "pegas", "dplyr", "tidyr", "ggplot2", "vegan", "reshape2")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
if (!requireNamespace("mmod", quietly = TRUE)) install.packages("mmod")

suppressPackageStartupMessages({
  library(adegenet)
  library(poppr)
  library(hierfstat)
  library(pegas)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(vegan)
  library(reshape2)
  library(mmod)
})

source(file.path("scripts", "_load_objects.R"))

`%||%` <- function(a, b) if (!is.null(a)) a else b

# ---- Reproducible output structure ----
DIRS <- c(
  file.path(RUN_OUT, "outputs", "tables"),
  file.path(RUN_OUT, "outputs", "qc"),
  file.path(RUN_OUT, "outputs", "robustness"),
  file.path(RUN_OUT, "outputs", "ibd_ibe"),
  file.path(RUN_OUT, "figures", "qc"),
  file.path(RUN_OUT, "figures", "differentiation"),
  file.path(RUN_OUT, "figures", "ibd_ibe"),
  file.path(RUN_OUT, "outputs", "logs")
)
for (d in DIRS) dir.create(d, recursive = TRUE, showWarnings = FALSE)

LOG_FILE <- file.path(RUN_OUT, "outputs", "logs", paste0("16_thesis_audit_hardening_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
log_con <- file(LOG_FILE, open = "wt")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message", append = TRUE)
on.exit({
  sink(type = "message")
  sink(type = "output")
  close(log_con)
}, add = TRUE)

cat("Running thesis-defense hardening script\n")
cat("Log:", LOG_FILE, "\n")

# ---- Helpers ----
read_matrix_csv <- function(path) {
  m <- as.matrix(read.csv(path, row.names = 1, check.names = FALSE))
  storage.mode(m) <- "numeric"
  m <- (m + t(m)) / 2
  diag(m) <- 0
  m
}

get_site_meta <- function(df_ids, metadata_csv = file.path("inputs", "site_metadata.csv")) {
  x <- as_tibble(df_ids)
  cand_site <- intersect(c("Site", "site", "Population", "population"), names(x))[1]
  if (is.na(cand_site)) stop("No site column in df_ids.")
  
  # Try site-level metadata file first (supports elevation + lat/lon)
  md <- NULL
  if (file.exists(metadata_csv)) {
    md <- suppressWarnings(read.csv(metadata_csv, stringsAsFactors = FALSE, comment.char = "#", check.names = FALSE))
    names(md) <- trimws(names(md))
  }
  
  if (!is.null(md) && nrow(md) > 0) {
    site_col <- names(md)[tolower(names(md)) %in% c("site", "population", "pop")][1]
    lat_col <- names(md)[tolower(names(md)) %in% c("latitude", "lat")][1]
    lon_col <- names(md)[tolower(names(md)) %in% c("longitude", "lon", "long")][1]
    elev_col <- names(md)[tolower(names(md)) %in% c("elevation_m", "elevation", "elev_m", "altitude")][1]
    if (!anyNA(c(site_col, lat_col, lon_col, elev_col))) {
      out <- md %>%
        transmute(
          Site = as.character(.data[[site_col]]),
          latitude = suppressWarnings(as.numeric(.data[[lat_col]])),
          longitude = suppressWarnings(as.numeric(.data[[lon_col]])),
          elevation_m = suppressWarnings(as.numeric(.data[[elev_col]]))
        ) %>%
        group_by(Site) %>%
        summarise(across(c(latitude, longitude, elevation_m), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
      return(out)
    }
  }
  
  # Fallback to columns already present in df_ids
  lat_col <- intersect(c("Latitude", "latitude", "Lat"), names(x))[1]
  lon_col <- intersect(c("Longitude", "longitude", "Lon", "Long"), names(x))[1]
  elev_col <- intersect(c("Elevation", "elevation", "elevation_m", "Elev_m"), names(x))[1]
  
  if (anyNA(c(lat_col, lon_col, elev_col))) {
    stop("Could not find complete site-level latitude/longitude/elevation metadata.")
  }
  
  x %>%
    transmute(
      Site = as.character(.data[[cand_site]]),
      latitude = suppressWarnings(as.numeric(.data[[lat_col]])),
      longitude = suppressWarnings(as.numeric(.data[[lon_col]])),
      elevation_m = suppressWarnings(as.numeric(.data[[elev_col]]))
    ) %>%
    group_by(Site) %>%
    summarise(across(c(latitude, longitude, elevation_m), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
}

compare_id_sets <- function(reference_ids, external_ids, reference_name = "reference", external_name = "external") {
  reference_ids <- unique(as.character(reference_ids))
  external_ids <- unique(as.character(external_ids))
  missing_in_external <- setdiff(reference_ids, external_ids)
  extra_in_external <- setdiff(external_ids, reference_ids)
  out <- tibble(
    reference_name = reference_name,
    external_name = external_name,
    n_reference = length(reference_ids),
    n_external = length(external_ids),
    n_missing_in_external = length(missing_in_external),
    n_extra_in_external = length(extra_in_external),
    perfect_match = length(missing_in_external) == 0 && length(extra_in_external) == 0,
    missing_ids = paste(head(missing_in_external, 25), collapse = ";"),
    extra_ids = paste(head(extra_in_external, 25), collapse = ";")
  )
  
  if (!out$perfect_match) {
    warning(
      sprintf(
        "ID mismatch (%s vs %s): missing=%d, extra=%d",
        reference_name, external_name, out$n_missing_in_external, out$n_extra_in_external
      )
    )
  }
  out
}

export_structure_input_from_gi_mll <- function(gi_mll, file_out = file.path(RUN_OUT, "outputs", "tables", "structure_input_from_gi_mll.tsv")) {
  gdf <- adegenet::genind2df(gi_mll, sep = "/")
  loci <- adegenet::locNames(gi_mll)
  popv <- as.character(pop(gi_mll))
  
  parse_two <- function(x) {
    x <- as.character(x)
    x[is.na(x) | x == "" | x == "NA/NA" | x == "0/0"] <- NA_character_
    sp <- strsplit(x, "/", fixed = TRUE)
    a1 <- vapply(sp, function(z) if (length(z) >= 1) z[1] else NA_character_, character(1))
    a2 <- vapply(sp, function(z) if (length(z) >= 2) z[2] else NA_character_, character(1))
    a1[is.na(x)] <- "-9"
    a2[is.na(x)] <- "-9"
    list(a1 = a1, a2 = a2)
  }
  
  out <- data.frame(ID = indNames(gi_mll), POP = popv, stringsAsFactors = FALSE)
  for (L in loci) {
    pz <- parse_two(gdf[[L]])
    out[[paste0(L, "_1")]] <- pz$a1
    out[[paste0(L, "_2")]] <- pz$a2
  }
  write.table(out, file = file_out, sep = "\t", quote = FALSE, row.names = FALSE)
  file_out
}

# ---- Consistency checks ----
consistency_tbl <- tibble(
  object = c("gi", "gi_mll"),
  nInd = c(nInd(gi), nInd(gi_mll)),
  nLoc = c(nLoc(gi), nLoc(gi_mll)),
  nPop = c(nPop(gi), nPop(gi_mll))
)
write.csv(consistency_tbl, file.path(RUN_OUT, "outputs", "qc", "dataset_consistency_counts.csv"), row.names = FALSE)

# Optional STRUCTURE ID check from common input files
structure_candidates <- c(
  file.path("inputs", "structure_input.txt"),
  file.path("inputs", "structure_input.tsv"),
  file.path("inputs", "structure_ids.txt"),
  file.path(RUN_OUT, "structure", "structure_input.txt")
)
structure_file <- structure_candidates[file.exists(structure_candidates)][1]
if (!is.na(structure_file)) {
  st <- read.table(structure_file, header = FALSE, stringsAsFactors = FALSE)
  ext_ids <- as.character(st[[1]])
  id_cmp <- compare_id_sets(indNames(gi_mll), ext_ids, "gi_mll", basename(structure_file))
} else {
  id_cmp <- tibble(
    reference_name = "gi_mll", external_name = "none_found",
    n_reference = length(indNames(gi_mll)), n_external = 0,
    n_missing_in_external = length(indNames(gi_mll)), n_extra_in_external = 0,
    perfect_match = FALSE, missing_ids = paste(head(indNames(gi_mll), 25), collapse = ";"), extra_ids = ""
  )
}
write.csv(id_cmp, file.path(RUN_OUT, "outputs", "qc", "structure_id_consistency_check.csv"), row.names = FALSE)

struct_path <- export_structure_input_from_gi_mll(gi_mll)
cat("STRUCTURE helper export written to:", struct_path, "\n")

# ---- HWE + multiple testing + locus flags ----
gi_site <- seppop(gi_mll)
hwe_rows <- list()
for (s in names(gi_site)) {
  g <- gi_site[[s]]
  if (nInd(g) < 4) next
  loci_df <- adegenet::genind2df(g, sep = "/")
  loci_df <- loci_df[, locNames(g), drop = FALSE]
  loci_obj <- tryCatch(pegas::as.loci(loci_df), error = function(e) NULL)
  if (is.null(loci_obj)) next
  hw <- tryCatch(pegas::hw.test(loci_obj, B = 1000), error = function(e) NULL)
  if (is.null(hw)) next
  pvals <- if (is.list(hw) && !is.null(hw$p.value)) hw$p.value else as.numeric(hw[, ncol(hw)])
  locs <- names(pvals) %||% rownames(hw)
  hwe_rows[[s]] <- tibble(Site = s, Locus = as.character(locs), p_raw = as.numeric(pvals))
}
hwe_tbl <- bind_rows(hwe_rows)
if (nrow(hwe_tbl) > 0) {
  hwe_tbl <- hwe_tbl %>%
    mutate(
      p_bh_global = p.adjust(p_raw, method = "BH"),
      p_bonf_global = p.adjust(p_raw, method = "bonferroni")
    ) %>%
    group_by(Site) %>%
    mutate(
      p_bh_site = p.adjust(p_raw, method = "BH"),
      bonf_threshold_site = 0.05 / n(),
      sig_bh_site = p_bh_site < 0.05,
      sig_bh_global = p_bh_global < 0.05,
      sig_bonf_site = p_raw < bonf_threshold_site
    ) %>%
    ungroup()
}
write.csv(hwe_tbl, file.path(RUN_OUT, "outputs", "qc", "hwe_site_locus_raw_and_adjusted.csv"), row.names = FALSE)

locus_sig_summary <- hwe_tbl %>%
  group_by(Locus) %>%
  summarise(
    n_sites_tested = n(),
    n_sig_bh_site = sum(sig_bh_site, na.rm = TRUE),
    n_sig_bh_global = sum(sig_bh_global, na.rm = TRUE),
    n_sig_bonf_site = sum(sig_bonf_site, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_sig_bh_site), desc(n_sig_bonf_site))
write.csv(locus_sig_summary, file.path(RUN_OUT, "outputs", "qc", "hwe_locus_significance_summary.csv"), row.names = FALSE)

flagged_loci <- locus_sig_summary %>% filter(n_sig_bh_site >= 2) %>% pull(Locus)
write.csv(data.frame(flagged_loci = flagged_loci), file.path(RUN_OUT, "outputs", "qc", "problem_loci_list.csv"), row.names = FALSE)

# Proxy null-allele estimator from heterozygote deficit
loc_stats <- hierfstat::basic.stats(genind2hierfstat(gi_mll))
ho <- colMeans(loc_stats$Ho, na.rm = TRUE)
he <- colMeans(loc_stats$Hs, na.rm = TRUE)
null_proxy <- pmax(0, pmin(1, (he - ho) / (1 + he)))
null_tbl <- tibble(Locus = names(he), Ho = as.numeric(ho), He = as.numeric(he), null_freq_proxy = as.numeric(null_proxy))
write.csv(null_tbl, file.path(RUN_OUT, "outputs", "qc", "null_allele_proxy_by_locus.csv"), row.names = FALSE)

# ---- Pairwise differentiation (clone-corrected only) ----
calc_pairwise <- function(gobj) {
  hf <- genind2hierfstat(gobj)
  pw_fst <- hierfstat::pairwise.WCfst(hf[,-1, drop = FALSE], hf[,1])
  pw_jost <- mmod::pairwise_D(gobj, linearized = FALSE)
  list(fst = as.matrix(pw_fst), jost = as.matrix(pw_jost))
}

pw_full <- calc_pairwise(gi_mll)
write.csv(pw_full$fst, file.path(RUN_OUT, "outputs", "tables", "pairwise_fst_wc_full.csv"))
write.csv(pw_full$jost, file.path(RUN_OUT, "outputs", "tables", "pairwise_jostd_full.csv"))

plot_heat <- function(mat, out_file, title_txt) {
  d <- reshape2::melt(mat, varnames = c("Site1", "Site2"), value.name = "value")
  p <- ggplot(d, aes(Site1, Site2, fill = value)) +
    geom_tile() +
    scale_fill_viridis_c(option = "magma") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = title_txt, x = NULL, y = NULL, fill = "Distance")
  ggsave(out_file, p, width = 8, height = 7, dpi = 300)
}
plot_heat(pw_full$fst, file.path(RUN_OUT, "figures", "differentiation", "heatmap_pairwise_fst_wc_full.png"), "Pairwise FST (Weir-Cockerham, full loci)")
plot_heat(pw_full$jost, file.path(RUN_OUT, "figures", "differentiation", "heatmap_pairwise_jostd_full.png"), "Pairwise Jost's D (full loci)")

# CI/permutation (if feasible)
fst_ci <- tryCatch(hierfstat::boot.ppfst(genind2hierfstat(gi_mll), nboot = 200), error = function(e) NULL)
if (!is.null(fst_ci)) {
  sink(file.path(RUN_OUT, "outputs", "tables", "pairwise_fst_bootstrap_summary.txt"))
  print(fst_ci)
  sink()
}

# ---- Sensitivity: exclude flagged loci ----
keep_loci <- setdiff(locNames(gi_mll), flagged_loci)
if (length(keep_loci) >= 3) {
  gi_filtered <- gi_mll[, keep_loci]
} else {
  gi_filtered <- gi_mll
}

# AMOVA
amova_full <- poppr::poppr.amova(gi_mll, ~Site)
amova_filt <- poppr::poppr.amova(gi_filtered, ~Site)
capture.output(amova_full, file = file.path(RUN_OUT, "outputs", "robustness", "amova_full.txt"))
capture.output(amova_filt, file = file.path(RUN_OUT, "outputs", "robustness", "amova_filtered_problem_loci_removed.txt"))

# PCA
x_full <- tab(gi_mll, freq = TRUE, NA.method = "mean")
x_filt <- tab(gi_filtered, freq = TRUE, NA.method = "mean")
pca_full <- prcomp(x_full, center = TRUE, scale. = FALSE)
pca_filt <- prcomp(x_filt, center = TRUE, scale. = FALSE)

pca_var <- tibble(
  dataset = c("full", "filtered"),
  pc1_var = c(summary(pca_full)$importance[2, 1], summary(pca_filt)$importance[2, 1]),
  pc2_var = c(summary(pca_full)$importance[2, 2], summary(pca_filt)$importance[2, 2])
)
write.csv(pca_var, file.path(RUN_OUT, "outputs", "robustness", "pca_variance_full_vs_filtered.csv"), row.names = FALSE)

pw_filtered <- calc_pairwise(gi_filtered)
write.csv(pw_filtered$fst, file.path(RUN_OUT, "outputs", "robustness", "pairwise_fst_wc_filtered.csv"))
write.csv(pw_filtered$jost, file.path(RUN_OUT, "outputs", "robustness", "pairwise_jostd_filtered.csv"))

robust_tbl <- tibble(
  metric = c("n_loci", "mean_pairwise_fst", "mean_pairwise_jostd", "amova_sigma_between", "amova_sigma_within", "pca_pc1_var", "pca_pc2_var"),
  full = c(
    nLoc(gi_mll),
    mean(pw_full$fst[upper.tri(pw_full$fst)], na.rm = TRUE),
    mean(pw_full$jost[upper.tri(pw_full$jost)], na.rm = TRUE),
    amova_full$componentsofcovariance[1, 1],
    amova_full$componentsofcovariance[nrow(amova_full$componentsofcovariance), 1],
    pca_var$pc1_var[pca_var$dataset == "full"],
    pca_var$pc2_var[pca_var$dataset == "full"]
  ),
  filtered = c(
    nLoc(gi_filtered),
    mean(pw_filtered$fst[upper.tri(pw_filtered$fst)], na.rm = TRUE),
    mean(pw_filtered$jost[upper.tri(pw_filtered$jost)], na.rm = TRUE),
    amova_filt$componentsofcovariance[1, 1],
    amova_filt$componentsofcovariance[nrow(amova_filt$componentsofcovariance), 1],
    pca_var$pc1_var[pca_var$dataset == "filtered"],
    pca_var$pc2_var[pca_var$dataset == "filtered"]
  )
) %>% mutate(delta_filtered_minus_full = filtered - full)
write.csv(robust_tbl, file.path(RUN_OUT, "outputs", "robustness", "ROBUSTNESS_SUMMARY.csv"), row.names = FALSE)

# ---- IBD + IBE + partial Mantel using site matrices ----
site_meta <- get_site_meta(df_ids)
if (anyNA(site_meta$latitude) || anyNA(site_meta$longitude) || anyNA(site_meta$elevation_m)) {
  stop("Site metadata has missing latitude/longitude/elevation values required for IBD/IBE.")
}

# geographic distance (km)
coords <- as.matrix(site_meta[, c("longitude", "latitude")])
rownames(coords) <- site_meta$Site
geo_mat <- as.matrix(dist(coords)) * 111 # approx km, robust fallback when geosphere unavailable
if (requireNamespace("geosphere", quietly = TRUE)) {
  geo_mat <- geosphere::distm(coords, fun = geosphere::distHaversine) / 1000
}
rownames(geo_mat) <- site_meta$Site
colnames(geo_mat) <- site_meta$Site

# elevation distance
elev_mat <- abs(outer(site_meta$elevation_m, site_meta$elevation_m, `-`))
rownames(elev_mat) <- site_meta$Site
colnames(elev_mat) <- site_meta$Site

# genetic matrices: Nei + FST + JostD (aligned to site names)
genetic_mats <- list()
nei_path <- file.path(RUN_OUT, "matrix_genetic_distance_nei.csv")
if (file.exists(nei_path)) genetic_mats$Nei <- read_matrix_csv(nei_path)
genetic_mats$FST <- pw_full$fst
genetic_mats$JostD <- pw_full$jost

mantel_results <- list()
pair_dfs <- list()
for (nm in names(genetic_mats)) {
  gm <- genetic_mats[[nm]]
  common <- Reduce(intersect, list(rownames(gm), rownames(geo_mat), rownames(elev_mat)))
  common <- sort(common)
  if (length(common) < 3) next
  g <- gm[common, common, drop = FALSE]
  gd <- geo_mat[common, common, drop = FALSE]
  ed <- elev_mat[common, common, drop = FALSE]
  
  m_ibd <- vegan::mantel(as.dist(g), as.dist(gd), permutations = 9999, method = "pearson")
  m_ibe <- vegan::mantel(as.dist(g), as.dist(ed), permutations = 9999, method = "pearson")
  m_p_ibd <- vegan::mantel.partial(as.dist(g), as.dist(gd), as.dist(ed), permutations = 9999, method = "pearson")
  m_p_ibe <- vegan::mantel.partial(as.dist(g), as.dist(ed), as.dist(gd), permutations = 9999, method = "pearson")
  
  mantel_results[[nm]] <- tibble(
    genetic_matrix = nm,
    test = c("mantel_genetic_geo", "mantel_genetic_elev", "partial_genetic_geo_given_elev", "partial_genetic_elev_given_geo"),
    r = c(m_ibd$statistic, m_ibe$statistic, m_p_ibd$statistic, m_p_ibe$statistic),
    p = c(m_ibd$signif, m_ibe$signif, m_p_ibd$signif, m_p_ibe$signif),
    permutations = 9999,
    n_sites = length(common),
    n_pairs = choose(length(common), 2)
  )
  
  idx <- which(upper.tri(g), arr.ind = TRUE)
  pair_dfs[[nm]] <- tibble(
    genetic_matrix = nm,
    Site1 = rownames(g)[idx[, 1]], Site2 = colnames(g)[idx[, 2]],
    genetic_distance = g[idx], geographic_km = gd[idx], elevation_distance_m = ed[idx]
  )
}

mantel_tbl <- bind_rows(mantel_results)
pairs_tbl <- bind_rows(pair_dfs)
write.csv(mantel_tbl, file.path(RUN_OUT, "outputs", "ibd_ibe", "mantel_ibd_ibe_partial_results.csv"), row.names = FALSE)
write.csv(pairs_tbl, file.path(RUN_OUT, "outputs", "ibd_ibe", "pairwise_site_distances_for_scatter.csv"), row.names = FALSE)
write.csv(geo_mat, file.path(RUN_OUT, "outputs", "ibd_ibe", "matrix_geographic_km.csv"))
write.csv(elev_mat, file.path(RUN_OUT, "outputs", "ibd_ibe", "matrix_elevation_absdiff_m.csv"))

plot_scatter <- function(df, xcol, out, ttl) {
  p <- ggplot(df, aes(.data[[xcol]], genetic_distance)) +
    geom_point(size = 2, alpha = 0.85, color = "#1f78b4") +
    geom_smooth(method = "lm", se = TRUE, color = "#d7301f", fill = "#fcbba1") +
    theme_minimal(base_size = 12) +
    labs(title = ttl, x = xcol, y = "Genetic distance")
  ggsave(out, p, width = 8, height = 6, dpi = 300)
}

if (nrow(pairs_tbl) > 0) {
  d_nei <- pairs_tbl %>% filter(genetic_matrix == "Nei")
  if (nrow(d_nei) > 0) {
    plot_scatter(d_nei, "geographic_km", file.path(RUN_OUT, "figures", "ibd_ibe", "scatter_nei_vs_geographic.png"), "Genetic vs geographic distance")
    plot_scatter(d_nei, "elevation_distance_m", file.path(RUN_OUT, "figures", "ibd_ibe", "scatter_nei_vs_elevation.png"), "Genetic vs elevation distance")
  }
}

cat("Hardening analyses complete.\n")