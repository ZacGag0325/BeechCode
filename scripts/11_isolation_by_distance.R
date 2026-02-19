
############################################################
# scripts/11_isolation_by_distance.R
# Isolation by Distance (population-level)
############################################################

options(repos = c(CRAN = "https://cloud.r-project.org"))
pkgs <- c("adegenet", "dplyr", "tidyr", "tibble", "ggplot2", "vegan", "geosphere")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)

suppressPackageStartupMessages({
  library(adegenet)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(vegan)
  library(geosphere)
})

set.seed(123)

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

normalize_site <- function(x) toupper(trimws(as.character(x)))

extract_site_coords <- function(df, source_name = "unknown") {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) return(NULL)
  
  nms <- names(df)
  nms_low <- tolower(nms)
  
  site_col <- nms[match(TRUE, nms_low %in% c("site", "population", "pop", "site_code", "code_site"), nomatch = 0)]
  lat_col <- nms[match(TRUE, nms_low %in% c("latitude", "lat"), nomatch = 0)]
  lon_col <- nms[match(TRUE, nms_low %in% c("longitude", "lon", "long", "lng"), nomatch = 0)]
  
  if (length(site_col) == 0 || length(lat_col) == 0 || length(lon_col) == 0 || !nzchar(site_col) || !nzchar(lat_col) || !nzchar(lon_col)) {
    return(NULL)
  }
  
  # optional grouping column for partial Mantel
  region_col <- nms[match(TRUE, nms_low %in% c("region_ns", "region", "north_south", "group", "zone"), nomatch = 0)]
  
  out <- tibble(
    Site_raw = as.character(df[[site_col]]),
    Latitude = suppressWarnings(as.numeric(df[[lat_col]])),
    Longitude = suppressWarnings(as.numeric(df[[lon_col]])),
    source = source_name
  )
  
  out$Region_NS <- if (length(region_col) > 0 && nzchar(region_col)) as.character(df[[region_col]]) else NA_character_
  out$Site_norm <- normalize_site(out$Site_raw)
  out
}

setwd(find_project_root())
source(file.path("scripts", "_load_objects.R"))

OUTDIR <- file.path(RUN_OUT, "ibd_only")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

if (!inherits(gi, "genind")) stop("Object 'gi' must be a genind object.")
site_levels <- normalize_site(levels(pop(gi)))
levels(pop(gi)) <- site_levels

if (length(site_levels) < 2) stop("Need at least 2 populations/sites in pop(gi) to run IBD.")

site_lookup <- tibble(Site = site_levels, Site_norm = site_levels)

meta_coords <- extract_site_coords(meta, "meta")
df_ids_coords <- extract_site_coords(df_ids, "df_ids")
coord_candidates <- bind_rows(meta_coords, df_ids_coords)

if (nrow(coord_candidates) == 0) {
  stop("Could not find Site + Latitude + Longitude columns in loaded objects (meta/df_ids).")
}

coord_summary <- coord_candidates %>%
  filter(!is.na(Site_norm) & nzchar(Site_norm)) %>%
  inner_join(site_lookup, by = "Site_norm") %>%
  filter(!is.na(Latitude), !is.na(Longitude)) %>%
  group_by(Site) %>%
  summarise(
    Latitude = mean(Latitude, na.rm = TRUE),
    Longitude = mean(Longitude, na.rm = TRUE),
    Region_NS = dplyr::first(Region_NS[!is.na(Region_NS) & nzchar(trimws(Region_NS))]),
    n_rows = n(),
    source = paste(sort(unique(source)), collapse = ","),
    .groups = "drop"
  )

coord_summary <- site_lookup %>%
  select(Site) %>%
  left_join(coord_summary, by = "Site")

missing_sites <- coord_summary$Site[is.na(coord_summary$Latitude) | is.na(coord_summary$Longitude)]
if (length(missing_sites) > 0) {
  stop("Missing Latitude/Longitude for sites: ", paste(missing_sites, collapse = ", "))
}

# optional derived Region_NS if absent
if (all(is.na(coord_summary$Region_NS) | !nzchar(trimws(coord_summary$Region_NS)))) {
  med_lat <- stats::median(coord_summary$Latitude, na.rm = TRUE)
  coord_summary <- coord_summary %>%
    mutate(Region_NS = ifelse(Latitude >= med_lat, "NORTH", "SOUTH"))
  message("Region_NS not found in metadata; derived from median latitude.")
}

coords <- coord_summary %>%
  select(Site, Latitude, Longitude, Region_NS, n_rows, source)

# Genetic distances (population-level Nei)
nei_obj <- tryCatch(
  adegenet::dist.genpop(adegenet::genind2genpop(gi), method = 1),
  error = function(e) e
)
if (inherits(nei_obj, "error")) stop("Could not compute genetic distances (Nei): ", nei_obj$message)

gen_mat <- as.matrix(nei_obj)
rownames(gen_mat) <- normalize_site(rownames(gen_mat))
colnames(gen_mat) <- normalize_site(colnames(gen_mat))

gen_aligned <- matrix(NA_real_, nrow = length(site_levels), ncol = length(site_levels), dimnames = list(site_levels, site_levels))
diag(gen_aligned) <- 0
common <- intersect(site_levels, intersect(rownames(gen_mat), colnames(gen_mat)))
if (length(common) >= 2) {
  gen_aligned[common, common] <- gen_mat[common, common, drop = FALSE]
}
gen_mat <- (gen_aligned + t(gen_aligned)) / 2
diag(gen_mat) <- 0

# Geographic great-circle distances (km)
coords <- coords %>% arrange(match(Site, site_levels))
geo_m <- geosphere::distm(as.matrix(coords[, c("Longitude", "Latitude")]), fun = geosphere::distHaversine)
geo_mat <- geo_m / 1000
rownames(geo_mat) <- coords$Site
colnames(geo_mat) <- coords$Site
geo_mat <- (geo_mat + t(geo_mat)) / 2
diag(geo_mat) <- 0
geo_mat <- geo_mat[site_levels, site_levels, drop = FALSE]

if (any(is.na(geo_mat[upper.tri(geo_mat)]))) {
  stop("Geographic distance matrix contains NA off-diagonal values.")
}

pair_idx <- which(upper.tri(gen_mat, diag = FALSE), arr.ind = TRUE)
pairs_tbl <- tibble(
  site1 = rownames(gen_mat)[pair_idx[, 1]],
  site2 = colnames(gen_mat)[pair_idx[, 2]],
  geo_km = geo_mat[pair_idx],
  gen_dist = gen_mat[pair_idx]
)

scatter_df <- pairs_tbl %>% filter(is.finite(geo_km), is.finite(gen_dist))
if (nrow(scatter_df) < 3) stop("Not enough complete pairwise comparisons for IBD (need at least 3).")

mantel_out <- vegan::mantel(
  xdis = stats::as.dist(gen_mat),
  ydis = stats::as.dist(geo_mat),
  method = "pearson",
  permutations = 9999
)

mantel_tbl <- tibble(
  analysis = "mantel_nei_vs_geo",
  statistic_r = as.numeric(mantel_out$statistic),
  p_value = as.numeric(mantel_out$signif),
  permutations = as.integer(mantel_out$permutations),
  n_pairs_used = nrow(scatter_df)
)

# Optional partial Mantel (only if >=2 non-empty groups)
region_clean <- toupper(trimws(coords$Region_NS))
region_ok <- !is.na(region_clean) & nzchar(region_clean)
if (sum(region_ok) >= 2 && length(unique(region_clean[region_ok])) >= 2) {
  reg_num <- as.numeric(factor(region_clean, levels = unique(region_clean)))
  reg_dist <- as.matrix(dist(reg_num, diag = TRUE, upper = TRUE))
  rownames(reg_dist) <- coords$Site
  colnames(reg_dist) <- coords$Site
  reg_dist <- reg_dist[site_levels, site_levels, drop = FALSE]
  
  partial_out <- tryCatch(
    vegan::mantel.partial(
      xdis = stats::as.dist(gen_mat),
      ydis = stats::as.dist(geo_mat),
      zdis = stats::as.dist(reg_dist),
      method = "pearson",
      permutations = 9999
    ),
    error = function(e) e
  )
  
  if (!inherits(partial_out, "error")) {
    mantel_tbl <- bind_rows(
      mantel_tbl,
      tibble(
        analysis = "partial_mantel_nei_vs_geo_given_region",
        statistic_r = as.numeric(partial_out$statistic),
        p_value = as.numeric(partial_out$signif),
        permutations = as.integer(partial_out$permutations),
        n_pairs_used = nrow(scatter_df)
      )
    )
  } else {
    message("Partial Mantel skipped: ", partial_out$message)
  }
} else {
  message("Partial Mantel skipped: Region_NS has <2 groups.")
}

p_scatter <- ggplot(scatter_df, aes(x = geo_km, y = gen_dist)) +
  geom_point(size = 2, alpha = 0.9, color = "#1f78b4") +
  geom_smooth(method = "lm", se = TRUE, color = "#d7301f", fill = "#fcae91", linewidth = 0.9) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Isolation by Distance (Population-Level)",
    x = "Geographic distance (km)",
    y = "Nei genetic distance"
  )

ggsave(file.path(OUTDIR, "ibd_scatter_nei_vs_geo.png"), p_scatter, width = 8, height = 6, dpi = 300)

write.csv(coords, file.path(OUTDIR, "site_coordinates_used.csv"), row.names = FALSE)
write.csv(geo_mat, file.path(OUTDIR, "matrix_geographic_distance_km.csv"), row.names = TRUE)
write.csv(gen_mat, file.path(OUTDIR, "matrix_genetic_distance_nei.csv"), row.names = TRUE)
write.csv(pairs_tbl, file.path(OUTDIR, "ibd_pairwise_table.csv"), row.names = FALSE)
write.csv(mantel_tbl, file.path(OUTDIR, "mantel_results.csv"), row.names = FALSE)

writeLines(
  c(
    "Isolation-by-Distance (population-level)",
    paste0("Sites: ", length(site_levels)),
    paste0("Pairs used: ", nrow(scatter_df)),
    paste0("Mantel r: ", signif(mantel_tbl$statistic_r[1], 4)),
    paste0("Mantel p: ", signif(mantel_tbl$p_value[1], 4))
  ),
  con = file.path(OUTDIR, "mantel_results.txt")
)

cat("DONE population-level IBD. Outputs in: ", OUTDIR, "\n", sep = "")
