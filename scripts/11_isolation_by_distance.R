############################################################
# scripts/11_isolation_by_distance.R
# Isolation by Distance (population-level)
############################################################

options(repos = c(CRAN = "https://cloud.r-project.org"))
pkgs <- c("adegenet", "dplyr", "tidyr", "tibble", "ggplot2", "vegan", "geosphere", "hierfstat")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)

suppressPackageStartupMessages({
  library(adegenet)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(vegan)
  library(geosphere)
  library(hierfstat)
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

setwd(find_project_root())
source(file.path("scripts", "_load_objects.R"))

OUTDIR <- file.path(RUN_OUT, "ibd_population")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

normalize_site <- function(x) toupper(trimws(as.character(x)))

site_levels <- levels(pop(gi))
if (length(site_levels) < 2) stop("Need at least 2 populations/sites in pop(gi) to run IBD.")

site_lookup <- tibble(
  Site = as.character(site_levels),
  Site_norm = normalize_site(site_levels)
)

extract_site_coords <- function(df, name = "unknown") {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) return(NULL)
  nms <- names(df)
  site_col <- nms[tolower(nms) %in% c("site", "population", "pop", "site_code", "code_site")]
  lat_col <- nms[tolower(nms) %in% c("latitude", "lat")]
  lon_col <- nms[tolower(nms) %in% c("longitude", "lon", "long", "lng")]
  if (length(site_col) == 0 || length(lat_col) == 0 || length(lon_col) == 0) return(NULL)
  
  out <- tibble(
    Site_raw = as.character(df[[site_col[1]]]),
    Latitude = suppressWarnings(as.numeric(df[[lat_col[1]]])),
    Longitude = suppressWarnings(as.numeric(df[[lon_col[1]]])),
    source = name
  ) %>%
    mutate(
      Site_norm = normalize_site(Site_raw)
    )
  
  reg_col <- nms[tolower(nms) %in% c("region", "north_south", "group", "zone")]
  if (length(reg_col) > 0) {
    out$region <- as.character(df[[reg_col[1]]])
  }
  
  out
}

meta_coords <- extract_site_coords(meta, "meta")
df_ids_coords <- extract_site_coords(df_ids, "df_ids")

coord_candidates <- bind_rows(meta_coords, df_ids_coords)
if (nrow(coord_candidates) == 0) {
  stop("Could not find a metadata table with Site + Latitude + Longitude columns in loaded objects (meta/df_ids).")
}

coord_candidates <- coord_candidates %>%
  filter(!is.na(Site_norm) & nzchar(Site_norm))

if (nrow(coord_candidates) == 0) {
  stop("No usable Site rows in metadata after normalization.")
}

coord_matched <- coord_candidates %>%
  inner_join(site_lookup, by = "Site_norm") %>%
  mutate(Site = Site.y) %>%
  select(Site, Site_norm, Latitude, Longitude, region, source)

if (nrow(coord_matched) == 0) {
  stop("No metadata Site codes matched pop(gi) Site names after trim/case normalization.")
}

coord_summary <- coord_matched %>%
  group_by(Site, Site_norm) %>%
  summarise(
    n_rows = n(),
    n_non_missing_coords = sum(!is.na(Latitude) & !is.na(Longitude)),
    n_unique_coord_pairs = n_distinct(paste(Latitude, Longitude)),
    Latitude = dplyr::first(na.omit(Latitude)),
    Longitude = dplyr::first(na.omit(Longitude)),
    region = dplyr::first(na.omit(region)),
    source = paste(sort(unique(source)), collapse = ","),
    .groups = "drop"
  )

if (any(coord_summary$n_unique_coord_pairs > 1, na.rm = TRUE)) {
  warning(
    "Some sites have multiple coordinate pairs in metadata; using first non-missing pair: ",
    paste(coord_summary$Site[coord_summary$n_unique_coord_pairs > 1], collapse = ", ")
  )
}

missing_coord_sites <- coord_summary$Site[is.na(coord_summary$Latitude) | is.na(coord_summary$Longitude)]
if (length(missing_coord_sites) > 0) {
  stop(
    "Missing Latitude/Longitude for matched sites: ",
    paste(missing_coord_sites, collapse = ", "),
    ". Fix metadata before running IBD."
  )
}

coord_summary <- site_lookup %>%
  left_join(coord_summary, by = c("Site", "Site_norm"))

if (any(is.na(coord_summary$Latitude) | is.na(coord_summary$Longitude))) {
  stop(
    "At least one pop(gi) site has no coordinates after metadata matching: ",
    paste(coord_summary$Site[is.na(coord_summary$Latitude) | is.na(coord_summary$Longitude)], collapse = ", ")
  )
}

coords <- coord_summary %>%
  select(Site, Latitude, Longitude, region) %>%
  distinct(Site, .keep_all = TRUE)

# A) Genetic distances (population-level)
nei_dist <- tryCatch(
  adegenet::dist.genpop(adegenet::genind2genpop(gi), method = 1),
  error = function(e) e
)
if (inherits(nei_dist, "error")) {
  stop("Could not compute Nei genetic distances: ", nei_dist$message)
}
nei_mat <- as.matrix(nei_dist)
nei_mat <- (nei_mat + t(nei_mat)) / 2
diag(nei_mat) <- 0

if (!identical(sort(rownames(nei_mat)), sort(site_levels))) {
  stop("Nei matrix site names do not align with pop(gi) Site names.")
}
nei_mat <- nei_mat[site_levels, site_levels, drop = FALSE]

fst_mat <- NULL
fst_status <- "not_computed"
fst_try <- tryCatch({
  hf <- hierfstat::genind2hierfstat(gi)
  pairwise <- hierfstat::pairwise.WCfst(hf[,-1, drop = FALSE], hf[,1])
  pairwise
}, error = function(e) e)

if (!inherits(fst_try, "error") && is.matrix(fst_try)) {
  fst_mat <- (fst_try + t(fst_try)) / 2
  diag(fst_mat) <- 0
  common <- intersect(site_levels, intersect(rownames(fst_mat), colnames(fst_mat)))
  if (length(common) >= 2) {
    aligned <- matrix(NA_real_, nrow = length(site_levels), ncol = length(site_levels), dimnames = list(site_levels, site_levels))
    diag(aligned) <- 0
    aligned[common, common] <- fst_mat[common, common, drop = FALSE]
    fst_mat <- aligned
    fst_status <- "computed_hierfstat_pairwise.WCfst"
  }
} else if (inherits(fst_try, "error")) {
  fst_status <- paste0("failed: ", fst_try$message)
}

# B) Geographic distances (great-circle, km)
geo_m <- geosphere::distm(
  x = as.matrix(coords[, c("Longitude", "Latitude")]),
  fun = geosphere::distHaversine
)
geo_mat <- geo_m / 1000
rownames(geo_mat) <- coords$Site
colnames(geo_mat) <- coords$Site
geo_mat <- (geo_mat + t(geo_mat)) / 2
diag(geo_mat) <- 0

if (!identical(rownames(geo_mat), site_levels)) {
  geo_mat <- geo_mat[site_levels, site_levels, drop = FALSE]
}

if (any(is.na(geo_mat[upper.tri(geo_mat)]))) {
  stop("Geographic distance matrix has NA off-diagonal values; likely caused by invalid coordinates.")
}
if (any(is.na(nei_mat[upper.tri(nei_mat)]))) {
  warning("Nei distance matrix contains NA off-diagonal values (e.g., monomorphic loci). Mantel will use complete cases only.")
}

# C) Mantel test
nei_dist_obj <- stats::as.dist(nei_mat)
geo_dist_obj <- stats::as.dist(geo_mat)
mantel_out <- vegan::mantel(nei_dist_obj, geo_dist_obj, method = "pearson", permutations = 9999)

mantel_tbl <- tibble(
  analysis = "mantel_nei_vs_geo",
  statistic_r = as.numeric(mantel_out$statistic),
  p_value = as.numeric(mantel_out$signif),
  permutations = as.integer(mantel_out$permutations)
)

if ("region" %in% names(coords) && sum(!is.na(coords$region)) >= 2) {
  reg_clean <- normalize_site(coords$region)
  if (length(unique(reg_clean[!is.na(reg_clean) & nzchar(reg_clean)])) >= 2) {
    reg_dist <- as.matrix(dist(as.numeric(factor(reg_clean))))
    rownames(reg_dist) <- coords$Site
    colnames(reg_dist) <- coords$Site
    reg_dist <- reg_dist[site_levels, site_levels, drop = FALSE]
    partial_out <- tryCatch(
      vegan::mantel.partial(nei_dist_obj, geo_dist_obj, stats::as.dist(reg_dist), method = "pearson", permutations = 9999),
      error = function(e) e
    )
    if (!inherits(partial_out, "error")) {
      mantel_tbl <- bind_rows(
        mantel_tbl,
        tibble(
          analysis = "partial_mantel_nei_vs_geo_given_region",
          statistic_r = as.numeric(partial_out$statistic),
          p_value = as.numeric(partial_out$signif),
          permutations = as.integer(partial_out$permutations)
        )
      )
    } else {
      mantel_tbl <- bind_rows(
        mantel_tbl,
        tibble(
          analysis = "partial_mantel_nei_vs_geo_given_region",
          statistic_r = NA_real_,
          p_value = NA_real_,
          permutations = 9999L
        )
      )
      warning("Partial Mantel skipped: ", partial_out$message)
    }
  }
}

# D + E) Pairwise tidy table and figures
pair_idx <- which(upper.tri(nei_mat, diag = FALSE), arr.ind = TRUE)
pairs_tbl <- tibble(
  site1 = rownames(nei_mat)[pair_idx[, 1]],
  site2 = colnames(nei_mat)[pair_idx[, 2]],
  geo_km = geo_mat[pair_idx],
  gen_dist = nei_mat[pair_idx]
)

if (!is.null(fst_mat) && is.matrix(fst_mat)) {
  fst_mat <- fst_mat[site_levels, site_levels, drop = FALSE]
  pairs_tbl$fst <- fst_mat[pair_idx]
}

if (nrow(pairs_tbl %>% filter(!is.na(geo_km), !is.na(gen_dist))) < 3) {
  stop("Not enough complete pairwise comparisons for plotting/regression after filtering NAs.")
}

scatter_df <- pairs_tbl %>% filter(!is.na(geo_km), !is.na(gen_dist))

p_scatter <- ggplot(scatter_df, aes(x = geo_km, y = gen_dist)) +
  geom_point(size = 2, alpha = 0.9, color = "#1f78b4") +
  geom_smooth(method = "lm", se = TRUE, color = "#d7301f", fill = "#fcae91", linewidth = 0.9) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Isolation by Distance (Population-Level)",
    x = "Geographic distance (km)",
    y = "Nei genetic distance"
  )

scatter_file <- file.path(OUTDIR, "ibd_scatter_nei_vs_geo.png")
ggsave(scatter_file, p_scatter, width = 8, height = 6, dpi = 300)

heat_df_nei <- as.data.frame(nei_mat) %>%
  tibble::rownames_to_column("Site1") %>%
  tidyr::pivot_longer(-Site1, names_to = "Site2", values_to = "value")
heat_df_geo <- as.data.frame(geo_mat) %>%
  tibble::rownames_to_column("Site1") %>%
  tidyr::pivot_longer(-Site1, names_to = "Site2", values_to = "value")

p_nei <- ggplot(heat_df_nei, aes(Site1, Site2, fill = value)) +
  geom_tile(color = "white", linewidth = 0.1) +
  scale_fill_gradient(low = "#f7fbff", high = "#08306b", na.value = "grey85") +
  coord_equal() +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank()) +
  labs(title = "Nei genetic distance matrix", x = "Site", y = "Site", fill = "Nei D")

p_geo <- ggplot(heat_df_geo, aes(Site1, Site2, fill = value)) +
  geom_tile(color = "white", linewidth = 0.1) +
  scale_fill_gradient(low = "#fff5f0", high = "#67000d", na.value = "grey85") +
  coord_equal() +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank()) +
  labs(title = "Geographic distance matrix (km)", x = "Site", y = "Site", fill = "km")

ggsave(file.path(OUTDIR, "heatmap_nei_distance.png"), p_nei, width = 8, height = 7, dpi = 300)
ggsave(file.path(OUTDIR, "heatmap_geographic_distance_km.png"), p_geo, width = 8, height = 7, dpi = 300)

write.csv(nei_mat, file.path(OUTDIR, "matrix_nei_distance.csv"), row.names = TRUE)
write.csv(geo_mat, file.path(OUTDIR, "matrix_geographic_distance_km.csv"), row.names = TRUE)
if (!is.null(fst_mat) && is.matrix(fst_mat)) {
  write.csv(fst_mat, file.path(OUTDIR, "matrix_pairwise_fst_wc.csv"), row.names = TRUE)
}

write.csv(coords, file.path(OUTDIR, "site_coordinates_used.csv"), row.names = FALSE)
write.csv(pairs_tbl, file.path(OUTDIR, "ibd_pairwise_table.csv"), row.names = FALSE)
write.csv(mantel_tbl, file.path(OUTDIR, "mantel_results.csv"), row.names = FALSE)

status_tbl <- tibble(
  metric = c("n_sites", "n_pairs", "fst_status"),
  value = c(as.character(length(site_levels)), as.character(nrow(pairs_tbl)), fst_status)
)
write.csv(status_tbl, file.path(OUTDIR, "ibd_status.csv"), row.names = FALSE)

cat("DONE population-level IBD. Outputs in: ", OUTDIR, "\n", sep = "")
