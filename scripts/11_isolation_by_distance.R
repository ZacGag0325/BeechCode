# scripts/11_isolation_by_distance.R
############################################################
# Isolation by Distance (Mantel test)
# Genetic distance input: pairwise Jost's D
# Rationale: For highly polymorphic microsatellites, Jost's D is
# treated as primary differentiation metric and therefore used for IBD.
# Outputs:
# - outputs/tables/mantel_test_results.csv
# - outputs/figures/isolation_by_distance_plot.jpeg
############################################################

suppressPackageStartupMessages({
  library(vegan)
  library(geosphere)
  library(ggplot2)
  library(dplyr)
})

source("scripts/_load_objects.R")

message("[11_isolation_by_distance] Running Mantel test: Jost's D vs geographic distance...")

# ----------------------------
# 1) Load genetic distance matrix (Jost's D)
# ----------------------------
jost_file <- file.path(MATRICES_DIR, "pairwise_jostD.csv")
if (!file.exists(jost_file)) {
  stop("Missing pairwise_jostD.csv. Run scripts/06_distance_matrices.R first.")
}
pairwise_jostD <- as.matrix(read.csv(jost_file, row.names = 1, check.names = FALSE))

# ----------------------------
# 2) Load site coordinates
# ----------------------------
site_meta_file <- file.path(PROJECT_ROOT, "inputs", "site_metadata.csv")
if (!file.exists(site_meta_file)) stop("Missing required file: inputs/site_metadata.csv")

site_meta <- read.csv(site_meta_file, stringsAsFactors = FALSE)
resolve_col <- function(df, patterns) {
  nms <- names(df)
  nms_low <- tolower(nms)
  idx <- match(TRUE, nms_low %in% patterns, nomatch = 0)
  if (idx > 0) return(nms[idx])
  NA_character_
}

site_col <- resolve_col(site_meta, c("site", "population", "pop"))
lat_col <- names(site_meta)[match(TRUE, grepl("lat", tolower(names(site_meta))), nomatch = 0)]
lon_col <- names(site_meta)[match(TRUE, grepl("lon|long", tolower(names(site_meta))), nomatch = 0)]

if (is.na(site_col) || is.na(lat_col) || is.na(lon_col)) {
  stop("site_metadata.csv must contain Site, Latitude, and Longitude columns.")
}

site_meta <- site_meta %>%
  transmute(
    Site = as.character(.data[[site_col]]),
    Latitude = as.numeric(.data[[lat_col]]),
    Longitude = as.numeric(.data[[lon_col]])
  ) %>%
  filter(!is.na(Site), !is.na(Latitude), !is.na(Longitude))

# ----------------------------
# 3) Align sites and compute geographic distance matrix
# ----------------------------
shared_sites <- intersect(rownames(pairwise_jostD), site_meta$Site)
if (length(shared_sites) < 3) {
  stop("Need at least 3 shared sites between pairwise_jostD.csv and site_metadata.csv")
}

pairwise_jostD <- pairwise_jostD[shared_sites, shared_sites, drop = FALSE]
site_meta <- site_meta[match(shared_sites, site_meta$Site), , drop = FALSE]

coords <- as.matrix(site_meta[, c("Longitude", "Latitude")])
geographic_km <- geosphere::distm(coords, fun = geosphere::distHaversine) / 1000
rownames(geographic_km) <- site_meta$Site
colnames(geographic_km) <- site_meta$Site

# ----------------------------
# 4) Mantel test
# ----------------------------
mantel_fit <- vegan::mantel(
  as.dist(pairwise_jostD),
  as.dist(geographic_km),
  method = "pearson",
  permutations = 9999
)

mantel_results <- data.frame(
  test = "Mantel_JostD_vs_Geographic",
  statistic_r = as.numeric(mantel_fit$statistic),
  p_value = as.numeric(mantel_fit$signif),
  permutations = as.integer(mantel_fit$permutations),
  stringsAsFactors = FALSE
)

mantel_file <- file.path(TABLES_DIR, "mantel_test_results.csv")
write.csv(mantel_results, mantel_file, row.names = FALSE)
message("[11_isolation_by_distance] Saved: ", mantel_file)

# ----------------------------
# 5) IBD plot with Mantel annotation
# ----------------------------
upper_idx <- upper.tri(pairwise_jostD)
plot_df <- data.frame(
  Geographic_km = geographic_km[upper_idx],
  JostD = pairwise_jostD[upper_idx]
)

ann_text <- sprintf("Mantel r = %.3f\np = %.4f", as.numeric(mantel_fit$statistic), as.numeric(mantel_fit$signif))

ibd_plot <- ggplot(plot_df, aes(Geographic_km, JostD)) +
  geom_point(size = 2.2, alpha = 0.8, color = "#2C3E50") +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
  annotate("text", x = Inf, y = Inf, label = ann_text, hjust = 1.05, vjust = 1.3, size = 4.1) +
  theme_bw(base_size = 12) +
  labs(
    title = "Isolation by distance",
    x = "Geographic distance (km)",
    y = "Pairwise Jost's D"
  )

ibd_plot_file <- file.path(FIGURES_DIR, "isolation_by_distance_plot.jpeg")
ggsave(ibd_plot_file, plot = ibd_plot, width = 8, height = 6, dpi = 320)
message("[11_isolation_by_distance] Saved: ", ibd_plot_file)