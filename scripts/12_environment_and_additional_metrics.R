############################################################
# scripts/12_environment_and_additional_metrics.R
# IBE (elevation), FST, private alleles, heterozygosity, HWE corrections
############################################################

options(repos = c(CRAN = "https://cloud.r-project.org"))
pkgs <- c("adegenet", "poppr", "hierfstat", "dplyr", "tidyr", "tibble", "ggplot2", "vegan")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)

suppressPackageStartupMessages({
  library(adegenet)
  library(poppr)
  library(hierfstat)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(vegan)
})

set.seed(20250218)

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

PROJECT_ROOT <- find_project_root()
setwd(PROJECT_ROOT)
source(file.path("scripts", "_load_objects.R"))

INPUTS_DIR <- file.path(PROJECT_ROOT, "inputs")
dir.create(INPUTS_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(RUN_OUT, showWarnings = FALSE, recursive = TRUE)

resolve_file <- function(candidates, label, required = TRUE) {
  existing <- candidates[file.exists(candidates)]
  if (length(existing) == 0) {
    if (required) {
      stop("Missing required ", label, ". Checked:\n- ", paste(candidates, collapse = "\n- "))
    } else {
      message("Optional input missing: ", label)
      return(NA_character_)
    }
  }
  if (length(existing) > 1) {
    info <- file.info(existing)$mtime
    pick <- existing[order(info, decreasing = TRUE)][1]
    message("Multiple matches found for ", label, ":\n- ", paste(existing, collapse = "\n- "), "\nSelected newest: ", pick)
    return(pick)
  }
  message("Resolved ", label, " from: ", existing[1])
  existing[1]
}

# Note: some downstream analyses are optional. Missing optional files are handled
# with graceful skip messages instead of stopping the entire pipeline.
if (!file.exists(file.path(RUN_OUT, "hwe_by_site_by_locus.csv"))) {
  stop("Missing required input: outputs/v1/hwe_by_site_by_locus.csv")
}

# Defensive diploid handling with scalar logical for control flow
ploidy_vec <- tryCatch(as.numeric(adegenet::ploidy(gi_mll)), error = function(e) rep(NA_real_, adegenet::nInd(gi_mll)))
diploid_vec <- (ploidy_vec == 2)
diploid_all <- isTRUE(all(diploid_vec, na.rm = TRUE))
if (anyNA(diploid_vec)) {
  message("Found NA ploidy values in gi_mll; diploid_all computed with na.rm = TRUE.")
}
message("Diploid check (gi_mll): all individuals diploid = ", diploid_all)

# Keep clone-corrected data as default working dataset
# and only run diploid-only metrics when diploid_all is TRUE.
gi_work <- gi_mll
run_diploid_only <- diploid_all
if (!run_diploid_only) {
  warning("Detected non-diploid/mixed ploidy in gi_mll. Diploid-only steps will be skipped.")
}

read_matrix_csv <- function(path) {
  mat_df <- read.csv(path, check.names = FALSE, row.names = 1)
  mat <- as.matrix(mat_df)
  storage.mode(mat) <- "numeric"
  if (nrow(mat) != ncol(mat)) stop("Matrix is not square: ", path)
  if (is.null(rownames(mat)) || is.null(colnames(mat))) stop("Matrix lacks row/column names: ", path)
  if (!identical(sort(rownames(mat)), sort(colnames(mat)))) stop("Row/column names differ in matrix: ", path)
  mat <- mat[sort(rownames(mat)), sort(colnames(mat)), drop = FALSE]
  mat <- (mat + t(mat)) / 2
  diag(mat) <- 0
  mat
}

nei_path <- resolve_file(c(file.path(RUN_OUT, "matrix_genetic_distance_nei.csv")), "Nei distance matrix", required = FALSE)
geo_path <- resolve_file(c(file.path(RUN_OUT, "matrix_geographic_distance_km.csv")), "geographic distance matrix", required = FALSE)

has_dist_mats <- !is.na(nei_path) && !is.na(geo_path)
if (has_dist_mats) {
  nei_mat <- read_matrix_csv(nei_path)
  geo_mat <- read_matrix_csv(geo_path)
  
  if (!identical(rownames(nei_mat), rownames(geo_mat))) {
    stop("Site names do not match between Nei and geographic distance matrices in outputs/v1.")
  }
  site_names <- rownames(nei_mat)
} else {
  message("Distance matrices not found; skipping Mantel/IBE distance-based steps.")
  nei_mat <- NULL
  geo_mat <- NULL
  site_names <- sort(unique(as.character(pop(gi_work))))
}

meta_candidates <- c(
  file.path(INPUTS_DIR, "site_metadata.csv"),
  file.path(PROJECT_ROOT, "data", "processed", "site_metadata.csv"),
  file.path(PROJECT_ROOT, "data", "raw", "site_metadata.csv")
)
meta_existing <- meta_candidates[file.exists(meta_candidates)]
meta_path <- if (length(meta_existing) > 0) resolve_file(meta_candidates, "site metadata") else NA_character_

if (is.na(meta_path) || !nzchar(meta_path)) {
  meta_path <- file.path(INPUTS_DIR, "site_metadata.csv")
  template <- data.frame(
    Site = site_names,
    latitude = NA_real_,
    longitude = NA_real_,
    elevation_m = NA_real_,
    Region_NS = ifelse(site_names %in% c("ML1", "ML2", "ML3", "CPF", "PLI", "LDF"), "North", "South"),
    stringsAsFactors = FALSE
  )
  writeLines(
    c(
      "# Required columns for downstream analyses:",
      "# Site, latitude, longitude, elevation_m, Region_NS",
      "# Region_NS values should be North or South.",
      "# Fill elevation_m to enable elevation-based analyses (IBE and elevation regressions)."
    ),
    con = meta_path
  )
  suppressWarnings(write.table(template, file = meta_path, sep = ",", row.names = FALSE, col.names = TRUE, append = TRUE))
  message("Created metadata template at: ", meta_path)
}

site_meta <- read.csv(meta_path, comment.char = "#", stringsAsFactors = FALSE, check.names = FALSE)
required_meta_cols <- c("Site", "latitude", "longitude", "elevation_m", "Region_NS")
if (!all(required_meta_cols %in% names(site_meta))) {
  stop(
    "Metadata file is missing required columns: ",
    paste(setdiff(required_meta_cols, names(site_meta)), collapse = ", "),
    "\nFile: ", meta_path
  )
}

site_meta <- site_meta %>%
  mutate(
    Site = trimws(as.character(Site)),
    elevation_m = suppressWarnings(as.numeric(elevation_m)),
    latitude = suppressWarnings(as.numeric(latitude)),
    longitude = suppressWarnings(as.numeric(longitude)),
    Region_NS = trimws(as.character(Region_NS))
  )

if (anyDuplicated(site_meta$Site) > 0) {
  dup_sites <- unique(site_meta$Site[duplicated(site_meta$Site)])
  stop("Duplicate Site entries found in metadata: ", paste(dup_sites, collapse = ", "))
}
if (!setequal(site_meta$Site, site_names)) {
  stop(
    "Site names in metadata do not exactly match matrix sites.\nMissing in metadata: ",
    paste(setdiff(site_names, site_meta$Site), collapse = ", "),
    "\nExtra in metadata: ",
    paste(setdiff(site_meta$Site, site_names), collapse = ", ")
  )
}
site_meta <- site_meta[match(site_names, site_meta$Site), , drop = FALSE]

find_site_order <- function(default_sites) {
  order_candidates <- list.files(RUN_OUT, pattern = "(?i)(order|south.*north|north.*south).*\\.csv$", full.names = TRUE)
  if (length(order_candidates) == 0) return(sort(default_sites))
  for (f in order_candidates) {
    ord <- tryCatch(read.csv(f, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
    if (is.null(ord) || nrow(ord) == 0) next
    site_col <- names(ord)[tolower(names(ord)) %in% c("site", "site_code", "population", "pop")]
    if (length(site_col) == 0) {
      if (ncol(ord) >= 1) {
        vec <- trimws(as.character(ord[[1]]))
      } else {
        next
      }
    } else {
      vec <- trimws(as.character(ord[[site_col[1]]]))
    }
    vec <- vec[nzchar(vec)]
    if (setequal(vec, default_sites) && length(vec) == length(default_sites)) return(vec)
  }
  sort(default_sites)
}

site_order <- find_site_order(site_names)

# A) Isolation-by-Environment (Elevation)
run_ibe <- function() {
  if (anyNA(site_meta$elevation_m)) {
    missing_sites <- site_meta$Site[is.na(site_meta$elevation_m)]
    stop(
      "Elevation-based analyses require non-missing elevation_m for all sites. Missing elevation_m for: ",
      paste(missing_sites, collapse = ", "),
      ". Fill inputs/site_metadata.csv and rerun."
    )
  }
  
  elev_diff <- abs(outer(site_meta$elevation_m, site_meta$elevation_m, "-"))
  rownames(elev_diff) <- site_meta$Site
  colnames(elev_diff) <- site_meta$Site
  elev_diff <- (elev_diff + t(elev_diff)) / 2
  diag(elev_diff) <- 0
  write.csv(elev_diff, file.path(RUN_OUT, "matrix_elevation_diff_m.csv"), row.names = TRUE)
  
  mantel_elev <- vegan::mantel(as.dist(nei_mat), as.dist(elev_diff), method = "pearson", permutations = 9999)
  writeLines(
    c(
      "Mantel test: Nei genetic distance vs elevation difference",
      paste0("Statistic r: ", as.numeric(mantel_elev$statistic)),
      paste0("p-value: ", as.numeric(mantel_elev$signif)),
      paste0("Permutations: ", as.integer(mantel_elev$permutations))
    ),
    con = file.path(RUN_OUT, "mantel_elevation_results.txt")
  )
  
  pair_idx <- which(upper.tri(nei_mat, diag = FALSE), arr.ind = TRUE)
  ibe_df <- data.frame(
    Site1 = rownames(nei_mat)[pair_idx[, 1]],
    Site2 = colnames(nei_mat)[pair_idx[, 2]],
    Nei_distance = nei_mat[pair_idx],
    Elevation_difference_m = elev_diff[pair_idx],
    stringsAsFactors = FALSE
  )
  
  p_ibe <- ggplot(ibe_df, aes(x = Elevation_difference_m, y = Nei_distance)) +
    geom_point(size = 2.5, alpha = 0.9, color = "#2c7fb8") +
    geom_smooth(method = "lm", se = TRUE, color = "#d95f0e", fill = "#fdd0a2", linewidth = 0.9) +
    theme_minimal(base_size = 12) +
    labs(
      title = "Isolation by Environment: Nei distance vs elevation difference",
      x = "Absolute elevation difference (m)",
      y = "Nei genetic distance"
    )
  ggsave(file.path(RUN_OUT, "ibe_scatter_nei_vs_elev.jpeg"), p_ibe, width = 8.5, height = 6.5, dpi = 350)
  
  partial_out <- vegan::mantel.partial(
    as.dist(nei_mat),
    as.dist(elev_diff),
    as.dist(geo_mat),
    method = "pearson",
    permutations = 9999
  )
  writeLines(
    c(
      "Partial Mantel: Nei genetic distance vs elevation difference | geographic distance",
      paste0("Statistic r: ", as.numeric(partial_out$statistic)),
      paste0("p-value: ", as.numeric(partial_out$signif)),
      paste0("Permutations: ", as.integer(partial_out$permutations))
    ),
    con = file.path(RUN_OUT, "mantel_partial_elev_given_geo.txt")
  )
}

if (has_dist_mats) {
  ibe_error <- tryCatch({
    run_ibe()
    NULL
  }, error = function(e) e)
  if (inherits(ibe_error, "error")) {
    writeLines(
      c("IBE skipped due to missing/invalid elevation metadata:", conditionMessage(ibe_error)),
      con = file.path(RUN_OUT, "mantel_elevation_results.txt")
    )
    writeLines(
      c("IBE partial Mantel skipped due to missing/invalid elevation metadata:", conditionMessage(ibe_error)),
      con = file.path(RUN_OUT, "mantel_partial_elev_given_geo.txt")
    )
  }
} else {
  writeLines(
    "IBE skipped: matrix_genetic_distance_nei.csv and/or matrix_geographic_distance_km.csv not found.",
    con = file.path(RUN_OUT, "mantel_elevation_results.txt")
  )
  writeLines(
    "Partial Mantel skipped: distance matrix inputs not found.",
    con = file.path(RUN_OUT, "mantel_partial_elev_given_geo.txt")
  )
}

# B) Allelic richness vs elevation
ar_path <- resolve_file(c(file.path(RUN_OUT, "allelic_richness_site_summary.csv")), "allelic richness site summary", required = FALSE)
if (is.na(ar_path)) {
  writeLines(
    "Skipped ar_vs_elevation: allelic_richness_site_summary.csv not found.",
    con = file.path(RUN_OUT, "ar_vs_elevation_stats.txt")
  )
} else {
  ar_site <- read.csv(ar_path, stringsAsFactors = FALSE, check.names = FALSE)
  ar_name_candidates <- names(ar_site)[tolower(names(ar_site)) %in% c("mean_allelic_richness", "allelic_richness", "ar_mean", "mean_ar")]
  if (length(ar_name_candidates) == 0) {
    numeric_cols <- names(ar_site)[sapply(ar_site, is.numeric)]
    numeric_cols <- setdiff(numeric_cols, c("n", "n_ind", "n_individuals"))
    ar_col <- numeric_cols[1]
  } else {
    ar_col <- ar_name_candidates[1]
  }
  if (is.na(ar_col) || !nzchar(ar_col) || !(ar_col %in% names(ar_site))) {
    stop("Could not identify allelic richness column in outputs/v1/allelic_richness_site_summary.csv")
  }
  if (!("Site" %in% names(ar_site))) stop("allelic_richness_site_summary.csv must include a Site column.")
  
  ar_df <- ar_site %>%
    mutate(Site = as.character(Site)) %>%
    left_join(site_meta %>% select(Site, elevation_m), by = "Site") %>%
    mutate(allelic_richness = suppressWarnings(as.numeric(.data[[ar_col]])))
  
  if (all(is.na(ar_df$elevation_m))) {
    writeLines(
      "Skipped ar_vs_elevation: elevation_m is missing for all sites in metadata.",
      con = file.path(RUN_OUT, "ar_vs_elevation_stats.txt")
    )
  } else {
    use_ar <- ar_df %>% filter(!is.na(elevation_m), !is.na(allelic_richness))
    if (nrow(use_ar) < 3) {
      writeLines(
        "Skipped ar_vs_elevation: fewer than 3 sites with both elevation_m and allelic richness.",
        con = file.path(RUN_OUT, "ar_vs_elevation_stats.txt")
      )
    } else {
      p_ar <- ggplot(use_ar, aes(x = elevation_m, y = allelic_richness)) +
        geom_point(size = 2.5, alpha = 0.9, color = "#1b9e77") +
        geom_smooth(method = "lm", se = TRUE, color = "#d95f02", fill = "#fdcc8a", linewidth = 0.9) +
        theme_minimal(base_size = 12) +
        labs(x = "Elevation (m)", y = "Allelic richness", title = "Allelic richness vs elevation")
      ggsave(file.path(RUN_OUT, "ar_vs_elevation.jpeg"), p_ar, width = 8, height = 6, dpi = 350)
      
      pear <- cor.test(use_ar$allelic_richness, use_ar$elevation_m, method = "pearson")
      spear <- cor.test(use_ar$allelic_richness, use_ar$elevation_m, method = "spearman", exact = FALSE)
      writeLines(
        c(
          "Allelic richness vs elevation correlations",
          paste0("n_sites_used: ", nrow(use_ar)),
          paste0("Pearson r: ", unname(pear$estimate), "; p-value: ", pear$p.value),
          paste0("Spearman rho: ", unname(spear$estimate), "; p-value: ", spear$p.value)
        ),
        con = file.path(RUN_OUT, "ar_vs_elevation_stats.txt")
      )
    }
  }
}

# C) Pairwise FST (Weir & Cockerham)
hf <- NULL
if (run_diploid_only) {
  hf <- hierfstat::genind2hierfstat(gi_work)
  fst_raw <- hierfstat::pairwise.WCfst(hf[, -1, drop = FALSE], hf[, 1])
  fst_mat <- as.matrix(fst_raw)
  fst_mat <- (fst_mat + t(fst_mat)) / 2
  diag(fst_mat) <- 0
  
  if (!setequal(site_names, intersect(rownames(fst_mat), colnames(fst_mat)))) {
    stop("Pairwise FST matrix site names do not match Nei/geographic matrix site names.")
  }
  fst_mat <- fst_mat[site_names, site_names, drop = FALSE]
  write.csv(fst_mat, file.path(RUN_OUT, "pairwise_fst_matrix.csv"), row.names = TRUE)
  
  fst_long <- as.data.frame(fst_mat) %>%
    tibble::rownames_to_column("Site1") %>%
    pivot_longer(-Site1, names_to = "Site2", values_to = "fst") %>%
    filter(Site1 < Site2)
  write.csv(fst_long, file.path(RUN_OUT, "pairwise_fst_long.csv"), row.names = FALSE)
  
  fst_plot_df <- as.data.frame(fst_mat) %>%
    tibble::rownames_to_column("Site1") %>%
    pivot_longer(-Site1, names_to = "Site2", values_to = "fst") %>%
    mutate(
      Site1 = factor(Site1, levels = site_order),
      Site2 = factor(Site2, levels = site_order)
    )
  
  p_fst <- ggplot(fst_plot_df, aes(x = Site1, y = Site2, fill = fst)) +
    geom_tile(color = "white", linewidth = 0.1) +
    scale_fill_gradient(low = "#f7fbff", high = "#08306b", na.value = "grey85") +
    coord_equal() +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank()) +
    labs(title = "Pairwise FST (Weir & Cockerham)", x = "Site", y = "Site", fill = "FST")
  ggsave(file.path(RUN_OUT, "pairwise_fst_heatmap.jpeg"), p_fst, width = 8.5, height = 7.5, dpi = 350)
} else {
  message("Skipping diploid-only FST calculations (non-diploid ploidy detected).")
  writeLines("Skipped pairwise FST: gi_mll is not fully diploid.", con = file.path(RUN_OUT, "pairwise_fst_notes.txt"))
}

# D) Private alleles per site
private_counts <- NULL
if ("private_alleles" %in% getNamespaceExports("poppr")) {
  pa_try <- tryCatch(poppr::private_alleles(gi_work, count.alleles = TRUE), error = function(e) NULL)
  if (is.data.frame(pa_try) && all(c("population", "count") %in% names(pa_try))) {
    private_counts <- pa_try %>%
      transmute(Site = as.character(population), private_alleles_total = as.numeric(count))
  }
}

if (is.null(private_counts)) {
  allele_tab <- adegenet::tab(gi_work, freq = FALSE, NA.method = "zero")
  loci_map <- as.character(adegenet::locFac(gi_work))
  site_fac <- as.character(pop(gi_work))
  allele_presence <- sapply(split(seq_len(ncol(allele_tab)), loci_map), function(cols) {
    loc_dat <- allele_tab[, cols, drop = FALSE]
    present <- sapply(unique(site_fac), function(s) {
      colSums(loc_dat[site_fac == s, , drop = FALSE], na.rm = TRUE) > 0
    })
    if (is.vector(present)) {
      present <- matrix(present, ncol = 1)
      colnames(present) <- unique(site_fac)
    }
    rownames(present) <- colnames(loc_dat)
    present
  }, simplify = FALSE)
  
  sites_u <- sort(unique(site_fac))
  per_locus_counts <- matrix(0, nrow = length(sites_u), ncol = length(allele_presence), dimnames = list(sites_u, names(allele_presence)))
  for (loc in names(allele_presence)) {
    pres <- allele_presence[[loc]]
    pres <- pres[, sites_u, drop = FALSE]
    present_n <- rowSums(pres)
    private_alleles <- which(present_n == 1)
    if (length(private_alleles) > 0) {
      owner_idx <- apply(pres[private_alleles, , drop = FALSE], 1, function(v) which(v)[1])
      owner_sites <- sites_u[owner_idx]
      per_locus_counts[owner_sites, loc] <- per_locus_counts[owner_sites, loc] + 1
    }
  }
  private_counts <- data.frame(
    Site = rownames(per_locus_counts),
    private_alleles_total = rowSums(per_locus_counts),
    private_alleles_per_locus_mean = rowMeans(per_locus_counts),
    stringsAsFactors = FALSE
  )
}

site_n <- as.data.frame(table(as.character(pop(gi_work))), stringsAsFactors = FALSE)
names(site_n) <- c("Site", "n_individuals")
private_out <- site_n %>%
  left_join(private_counts, by = "Site") %>%
  mutate(
    private_alleles_total = replace_na(as.numeric(private_alleles_total), 0),
    private_alleles_per_locus_mean = if ("private_alleles_per_locus_mean" %in% names(.)) replace_na(as.numeric(private_alleles_per_locus_mean), 0) else NA_real_
  ) %>%
  arrange(match(Site, site_order))

write.csv(private_out, file.path(RUN_OUT, "private_alleles_by_site.csv"), row.names = FALSE)

p_priv <- ggplot(private_out, aes(x = factor(Site, levels = site_order), y = private_alleles_total)) +
  geom_col(fill = "#3182bd") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Site", y = "Private alleles (total)", title = "Private alleles by site")
ggsave(file.path(RUN_OUT, "private_alleles_barplot.jpeg"), p_priv, width = 8.5, height = 6.5, dpi = 350)

top_private <- private_out %>% arrange(desc(private_alleles_total), Site)
writeLines(
  c(
    "Private alleles summary by site",
    paste0("Top site(s): ", paste(top_private$Site[top_private$private_alleles_total == max(top_private$private_alleles_total, na.rm = TRUE)], collapse = ", ")),
    paste0("Maximum private alleles observed: ", max(top_private$private_alleles_total, na.rm = TRUE)),
    "Interpretation note: higher private-allele counts can indicate stronger site-specific differentiation and/or sampling effects."
  ),
  con = file.path(RUN_OUT, "private_alleles_summary.txt")
)

# E) Heterozygosity by site
# Estimator: hierfstat::basic.stats; Ho = observed heterozygosity, He = Hs within-pop expected heterozygosity.
if (run_diploid_only && !is.null(hf)) {
  bs <- hierfstat::basic.stats(hf)
  ho_mat <- bs$Ho
  he_mat <- bs$Hs
  fis_mat <- bs$Fis
  
  het_df <- data.frame(
    Site = colnames(ho_mat),
    Ho_mean = colMeans(ho_mat, na.rm = TRUE),
    He_mean = colMeans(he_mat, na.rm = TRUE),
    FIS_mean = colMeans(fis_mat, na.rm = TRUE),
    stringsAsFactors = FALSE
  ) %>%
    left_join(site_n, by = "Site") %>%
    arrange(match(Site, site_order))
  
  write.csv(het_df, file.path(RUN_OUT, "heterozygosity_by_site.csv"), row.names = FALSE)
  
  het_long <- het_df %>%
    select(Site, Ho_mean, He_mean) %>%
    pivot_longer(cols = c(Ho_mean, He_mean), names_to = "Metric", values_to = "Value") %>%
    mutate(Site = factor(Site, levels = site_order))
  
  p_het <- ggplot(het_long, aes(x = Site, y = Value, fill = Metric)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    scale_fill_manual(values = c(Ho_mean = "#1b9e77", He_mean = "#7570b3")) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Observed and expected heterozygosity by site", x = "Site", y = "Mean heterozygosity")
  ggsave(file.path(RUN_OUT, "Ho_He_by_site.jpeg"), p_het, width = 9, height = 6.5, dpi = 350)
  
  if (all(is.na(site_meta$elevation_m))) {
    writeLines("Skipped He_vs_elevation: elevation_m is missing for all sites in metadata.", con = file.path(RUN_OUT, "He_vs_elevation.txt"))
  } else {
    he_elev <- het_df %>% left_join(site_meta %>% select(Site, elevation_m), by = "Site") %>% filter(!is.na(elevation_m), !is.na(He_mean))
    if (nrow(he_elev) >= 3) {
      p_he <- ggplot(he_elev, aes(x = elevation_m, y = He_mean)) +
        geom_point(size = 2.5, alpha = 0.9, color = "#386cb0") +
        geom_smooth(method = "lm", se = TRUE, color = "#ef3b2c", fill = "#fcbba1", linewidth = 0.9) +
        theme_minimal(base_size = 12) +
        labs(x = "Elevation (m)", y = "Expected heterozygosity (He)", title = "He vs elevation")
      ggsave(file.path(RUN_OUT, "He_vs_elevation.jpeg"), p_he, width = 8, height = 6, dpi = 350)
    } else {
      writeLines("Skipped He_vs_elevation: fewer than 3 sites with both He and elevation.", con = file.path(RUN_OUT, "He_vs_elevation.txt"))
    }
  }
} else {
  message("Skipping diploid-only heterozygosity calculations (non-diploid ploidy detected).")
  writeLines("Skipped heterozygosity_by_site and He_vs_elevation: gi_mll is not fully diploid.", con = file.path(RUN_OUT, "heterozygosity_notes.txt"))
}

# F) HWE p-value adjustments
hwe_path <- resolve_file(c(file.path(RUN_OUT, "hwe_by_site_by_locus.csv")), "HWE by-site-by-locus table")
hwe_df <- read.csv(hwe_path, stringsAsFactors = FALSE, check.names = FALSE)
adjusted_path <- file.path(RUN_OUT, "hwe_by_site_by_locus_adjusted.csv")
notes_path <- file.path(RUN_OUT, "hwe_adjustment_notes.txt")

adj_cols <- c("p_bonf_site", "p_fdr_site", "p_bonf_global", "p_fdr_global")
if (all(adj_cols %in% names(hwe_df))) {
  write.csv(hwe_df, adjusted_path, row.names = FALSE)
  writeLines("Adjusted p-value columns already existed in hwe_by_site_by_locus.csv; copied without recalculation.", con = notes_path)
} else {
  p_candidates <- names(hwe_df)[tolower(names(hwe_df)) %in% c("p", "p_value", "pvalue", "pval", "p_raw", "pvalue_raw")]
  if (length(p_candidates) == 0) {
    p_candidates <- names(hwe_df)[grepl("^p", tolower(names(hwe_df)))]
  }
  if (length(p_candidates) == 0) stop("Could not identify raw p-value column in hwe_by_site_by_locus.csv")
  p_col <- p_candidates[1]
  
  if (!("Site" %in% names(hwe_df))) stop("hwe_by_site_by_locus.csv must include a Site column for per-site adjustments.")
  p_raw <- suppressWarnings(as.numeric(hwe_df[[p_col]]))
  
  hwe_adj <- hwe_df %>%
    mutate(.p_raw = p_raw) %>%
    group_by(Site) %>%
    mutate(
      p_bonf_site = p.adjust(.p_raw, method = "bonferroni"),
      p_fdr_site = p.adjust(.p_raw, method = "BH")
    ) %>%
    ungroup() %>%
    mutate(
      p_bonf_global = p.adjust(.p_raw, method = "bonferroni"),
      p_fdr_global = p.adjust(.p_raw, method = "BH")
    ) %>%
    select(-.p_raw)
  
  write.csv(hwe_adj, adjusted_path, row.names = FALSE)
  writeLines(
    c(
      "HWE multiple-testing correction notes",
      paste0("Raw p-value column used: ", p_col),
      "Per-site corrections were computed across loci within each Site (Bonferroni and BH/FDR).",
      "Global corrections were computed across all Site x Locus tests (Bonferroni and BH/FDR).",
      "Original hwe_by_site_by_locus.csv was not overwritten."
    ),
    con = notes_path
  )
}

cat("DONE environment + additional metrics. Outputs in: ", RUN_OUT, "\n", sep = "")
