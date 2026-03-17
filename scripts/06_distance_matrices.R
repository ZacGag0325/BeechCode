# scripts/06_distance_matrices.R
############################################################
# Pairwise population differentiation matrices (gi_mll)
#
# IMPORTANT (microsatellite context):
# - Jost's D is treated as the primary differentiation metric.
# - FST is retained as a secondary comparison metric.
#
# Outputs:
# - outputs/matrices/pairwise_jostD.csv
# - outputs/matrices/pairwise_fst.csv
# - outputs/tables/pairwise_jostD_long.csv
# - outputs/tables/pairwise_fst_long.csv
# - outputs/figures/pairwise_jostD_heatmap.jpeg
# - outputs/figures/pairwise_fst_heatmap.jpeg
# - outputs/tables/site_differentiation_summary.csv
############################################################

suppressPackageStartupMessages({
  library(adegenet)
  library(mmod)
  library(hierfstat)
  library(ggplot2)
  library(dplyr)
})

source("scripts/_load_objects.R")

message("[06_distance_matrices] Using object: gi_mll (clone-corrected)")
message("[06_distance_matrices] Jost's D function: mmod::pairwise_D")
message("[06_distance_matrices] FST function: hierfstat::pairwise.WCfst (fallback: hierfstat::pairwise.neifst)")

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
to_numeric_matrix <- function(x, name = "object") {
  out <- NULL
  
  if (inherits(x, "dist")) {
    out <- as.matrix(x)
  } else if (is.matrix(x)) {
    out <- x
  } else if (is.data.frame(x)) {
    out <- as.matrix(x)
  } else {
    out <- tryCatch(as.matrix(x), error = function(e) NULL)
  }
  
  if (is.null(out) || !is.matrix(out)) {
    stop("[06_distance_matrices] Could not coerce ", name, " to matrix.")
  }
  
  suppressWarnings(storage.mode(out) <- "numeric")
  if (!is.numeric(out)) {
    stop("[06_distance_matrices] ", name, " matrix is not numeric after coercion.")
  }
  
  out
}

ensure_square_named <- function(m, fallback_names, name = "matrix") {
  if (!is.matrix(m) || nrow(m) != ncol(m)) {
    stop("[06_distance_matrices] ", name, " is not a square matrix.")
  }
  
  rn <- rownames(m)
  cn <- colnames(m)
  
  if (is.null(rn) || is.null(cn) || any(rn == "") || any(cn == "") || !setequal(rn, cn)) {
    if (length(fallback_names) == nrow(m)) {
      rownames(m) <- fallback_names
      colnames(m) <- fallback_names
    } else {
      stop("[06_distance_matrices] ", name, " lacks valid dimnames and fallback names are incompatible.")
    }
  }
  
  m <- m[rownames(m), rownames(m), drop = FALSE]
  m
}

make_long_pairwise <- function(m, value_name) {
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
  list(all = long_all, unique = long_unique %>% select(Site1, Site2, all_of(value_name)))
}

# ------------------------------------------------------------
# 1) Pairwise Jost's D (primary)
# ------------------------------------------------------------
jost_raw <- mmod::pairwise_D(gi_mll, linearized = FALSE)
pairwise_jostD <- to_numeric_matrix(jost_raw, name = "Jost's D")

site_levels <- sort(unique(as.character(adegenet::pop(gi_mll))))
pairwise_jostD <- ensure_square_named(pairwise_jostD, fallback_names = site_levels, name = "Jost's D")
diag(pairwise_jostD) <- 0
message("[06_distance_matrices] Set Jost's D diagonal to 0 for downstream distance compatibility.")

message("[06_distance_matrices] Jost's D matrix dimensions: ", nrow(pairwise_jostD), " x ", ncol(pairwise_jostD))

jost_file <- file.path(MATRICES_DIR, "pairwise_jostD.csv")
write.csv(pairwise_jostD, jost_file)
message("[06_distance_matrices] Saved: ", jost_file)

jost_long <- make_long_pairwise(pairwise_jostD, "JostD")
jost_long_file <- file.path(TABLES_DIR, "pairwise_jostD_long.csv")
write.csv(jost_long$unique, jost_long_file, row.names = FALSE)
message("[06_distance_matrices] Saved: ", jost_long_file)

# ------------------------------------------------------------
# 2) Pairwise FST (secondary)
# ------------------------------------------------------------
# Use hierfstat-formatted data.frame directly.
# NOTE: pairwise.WCfst signature expects a single hierfstat data object;
# passing a population vector as second argument is interpreted as diploid
# and can trigger 'condition has length > 1'.
hf <- hierfstat::genind2hierfstat(gi_mll)
if (!is.data.frame(hf) || ncol(hf) < 2) {
  stop("[06_distance_matrices] genind2hierfstat returned invalid object for FST.")
}

# Ensure first column is numeric population coding
if (!is.numeric(hf[[1]])) {
  hf[[1]] <- as.integer(factor(hf[[1]]))
}

fst_raw <- tryCatch(
  hierfstat::pairwise.WCfst(hf),
  error = function(e) {
    warning("[06_distance_matrices] pairwise.WCfst failed: ", conditionMessage(e),
            " | Falling back to hierfstat::pairwise.neifst")
    hierfstat::pairwise.neifst(hf)
  }
)

pairwise_fst <- to_numeric_matrix(fst_raw, name = "pairwise FST")

# Name rows/cols with site labels
site_lookup <- tapply(as.character(adegenet::pop(gi_mll)), hf[[1]], function(x) names(sort(table(x), decreasing = TRUE))[1])
site_lookup <- as.character(site_lookup)
if (length(site_lookup) == nrow(pairwise_fst)) {
  rownames(pairwise_fst) <- site_lookup
  colnames(pairwise_fst) <- site_lookup
}

pairwise_fst <- ensure_square_named(pairwise_fst, fallback_names = sort(unique(as.character(adegenet::pop(gi_mll)))), name = "pairwise FST")
diag(pairwise_fst) <- 0
message("[06_distance_matrices] Set FST diagonal to 0 for downstream distance compatibility.")

message("[06_distance_matrices] FST matrix dimensions: ", nrow(pairwise_fst), " x ", ncol(pairwise_fst))

fst_file <- file.path(MATRICES_DIR, "pairwise_fst.csv")
write.csv(pairwise_fst, fst_file)
message("[06_distance_matrices] Saved: ", fst_file)

fst_long <- make_long_pairwise(pairwise_fst, "FST")
fst_long_file <- file.path(TABLES_DIR, "pairwise_fst_long.csv")
write.csv(fst_long$unique, fst_long_file, row.names = FALSE)
message("[06_distance_matrices] Saved: ", fst_long_file)

# ------------------------------------------------------------
# 3) Heatmaps (publication-ready)
# ------------------------------------------------------------
heat_theme <- theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold")
  )

jost_heatmap <- ggplot(jost_long$all, aes(Site1, Site2, fill = JostD)) +
  geom_tile(color = "grey90", linewidth = 0.2) +
  scale_fill_viridis_c(option = "C", na.value = "white") +
  coord_equal() +
  heat_theme +
  labs(title = "Pairwise Jost's D (primary metric)", x = NULL, y = NULL, fill = "Jost's D")

jost_fig <- file.path(FIGURES_DIR, "pairwise_jostD_heatmap.jpeg")
ggsave(jost_fig, plot = jost_heatmap, width = 7.2, height = 6.2, dpi = 320)
message("[06_distance_matrices] Saved: ", jost_fig)

fst_heatmap <- ggplot(fst_long$all, aes(Site1, Site2, fill = FST)) +
  geom_tile(color = "grey90", linewidth = 0.2) +
  scale_fill_viridis_c(option = "B", na.value = "white") +
  coord_equal() +
  heat_theme +
  labs(title = "Pairwise FST (secondary metric)", x = NULL, y = NULL, fill = "FST")

fst_fig <- file.path(FIGURES_DIR, "pairwise_fst_heatmap.jpeg")
ggsave(fst_fig, plot = fst_heatmap, width = 7.2, height = 6.2, dpi = 320)
message("[06_distance_matrices] Saved: ", fst_fig)

# ------------------------------------------------------------
# 4) Site-level differentiation summary
# ------------------------------------------------------------
common_sites <- intersect(rownames(pairwise_jostD), rownames(pairwise_fst))
if (length(common_sites) < 2) {
  stop("[06_distance_matrices] Fewer than 2 shared sites between Jost's D and FST matrices.")
}

jost_use <- pairwise_jostD[common_sites, common_sites, drop = FALSE]
fst_use <- pairwise_fst[common_sites, common_sites, drop = FALSE]

site_differentiation_summary <- data.frame(
  Site = common_sites,
  mean_pairwise_JostD = vapply(common_sites, function(s) mean(jost_use[s, ], na.rm = TRUE), numeric(1)),
  mean_pairwise_FST = vapply(common_sites, function(s) mean(fst_use[s, ], na.rm = TRUE), numeric(1)),
  stringsAsFactors = FALSE
) %>%
  arrange(Site)

summary_file <- file.path(TABLES_DIR, "site_differentiation_summary.csv")
write.csv(site_differentiation_summary, summary_file, row.names = FALSE)
message("[06_distance_matrices] Saved: ", summary_file)