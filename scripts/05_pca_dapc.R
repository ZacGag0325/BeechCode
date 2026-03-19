# scripts/05_pca_dapc.R
############################################################
# PCA + DAPC (clone-corrected: gi_mll)
#
# Why gi_mll here?
# - Population-level structure analyses assume approximately independent
#   genotypes. Repeated ramets can dominate ordinations and inflate apparent
#   clustering, so clone-corrected MLL representatives are used.
#
# Why use Site as the DAPC grouping variable?
# - In this workflow, Site is the biologically interpretable a priori group
#   used throughout population-level summaries (HWE, differentiation, AMOVA,
#   IBD). Using Site here keeps the DAPC directly comparable to the rest of
#   the thesis pipeline.
# - STRUCTURE was run externally on the full dataset and remains separate.
#
# Outputs:
# - outputs/tables/pca_scores.csv
# - outputs/tables/pca_eigenvalues.csv
# - outputs/figures/pca_plot.jpeg (PC1 vs PC2)
# - outputs/figures/supplementary/PC1_vs_PC3_plot.jpeg
# - outputs/figures/supplementary/PC2_vs_PC3_plot.jpeg
# - outputs/tables/dapc_coordinates.csv
# - outputs/tables/dapc_group_centroids.csv
# - outputs/tables/dapc_model_metadata.csv
# - outputs/figures/dapc_plot.jpeg
############################################################

suppressPackageStartupMessages({
  library(adegenet)
  library(ggplot2)
  library(dplyr)
})

source("scripts/_load_objects.R")

message("[05_pca_dapc] Running PCA and DAPC on gi_mll...")

if (!inherits(gi_mll, "genind")) {
  stop("[05_pca_dapc] gi_mll must be a genind object. Run scripts/00_master_pipeline.R first.")
}

validate_columns(df_ids_mll, c("ind_id", "Site"), df_name = "[05_pca_dapc] df_ids_mll")
if (!all(adegenet::indNames(gi_mll) == df_ids_mll$ind_id)) {
  stop("[05_pca_dapc] gi_mll and df_ids_mll are not aligned.")
}
if (!identical(as.character(adegenet::pop(gi_mll)), as.character(df_ids_mll$Site))) {
  stop("[05_pca_dapc] pop(gi_mll) must match df_ids_mll$Site for consistent grouping.")
}

site_factor <- droplevels(factor(as.character(df_ids_mll$Site)))
site_counts <- table(site_factor)

if (nlevels(site_factor) < 2) {
  stop("[05_pca_dapc] DAPC requires at least two Site groups.")
}
if (any(site_counts < 2)) {
  warning(
    "[05_pca_dapc] Some Site groups contain fewer than 2 clone-corrected individuals: ",
    paste(names(site_counts)[site_counts < 2], collapse = ", "),
    ". These sites will still appear in PCA; DAPC interpretation for those groups should be cautious."
  )
}

# ------------------------------------------------------------
# 1) PCA
# ------------------------------------------------------------
# adegenet::tab gives allele frequencies/coded counts from the genind object.
# Missing data are replaced by mean allele frequencies so the PCA can proceed
# without discarding entire individuals or loci.
X <- adegenet::tab(gi_mll, freq = TRUE, NA.method = "mean")

if (!is.matrix(X) || nrow(X) < 2 || ncol(X) < 2) {
  stop("[05_pca_dapc] Allele-frequency table from gi_mll is too small for PCA.")
}

keep_cols <- apply(X, 2, function(v) stats::var(v, na.rm = TRUE) > 0)
if (!any(keep_cols)) {
  stop("[05_pca_dapc] All PCA columns are invariant after NA imputation.")
}
X <- X[, keep_cols, drop = FALSE]

pca_fit <- stats::prcomp(X, center = TRUE, scale. = FALSE)
var_exp <- summary(pca_fit)$importance[2, ] * 100
scores_mat <- pca_fit$x
n_pc <- ncol(scores_mat)

pca_scores <- data.frame(
  Individual = rownames(scores_mat),
  Site = as.character(site_factor),
  PC1 = if (n_pc >= 1) scores_mat[, 1] else NA_real_,
  PC2 = if (n_pc >= 2) scores_mat[, 2] else NA_real_,
  PC3 = if (n_pc >= 3) scores_mat[, 3] else NA_real_,
  stringsAsFactors = FALSE
) %>%
  filter(!is.na(Site), nzchar(Site))

pca_scores_file <- file.path(TABLES_DIR, "pca_scores.csv")
write.csv(pca_scores, pca_scores_file, row.names = FALSE)
message("[05_pca_dapc] Saved: ", pca_scores_file)

pca_eigenvalues <- data.frame(
  PC = paste0("PC", seq_along(pca_fit$sdev)),
  Eigenvalue = pca_fit$sdev^2,
  Percent_Variance = (pca_fit$sdev^2 / sum(pca_fit$sdev^2)) * 100,
  stringsAsFactors = FALSE
)

pca_eig_file <- file.path(TABLES_DIR, "pca_eigenvalues.csv")
write.csv(pca_eigenvalues, pca_eig_file, row.names = FALSE)
message("[05_pca_dapc] Saved: ", pca_eig_file)

plot_theme <- theme_bw(base_size = 12) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

plot_12_df <- pca_scores %>% filter(is.finite(PC1), is.finite(PC2))
if (nrow(plot_12_df) == 0) {
  stop("[05_pca_dapc] No valid rows for PC1 vs PC2 plot.")
}

pca_plot_12 <- ggplot(plot_12_df, aes(PC1, PC2, color = Site)) +
  geom_point(size = 2.3, alpha = 0.9) +
  plot_theme +
  labs(
    title = "PCA (PC1 vs PC2)",
    x = sprintf("PC1 (%.1f%%)", ifelse(length(var_exp) >= 1, var_exp[1], NA_real_)),
    y = sprintf("PC2 (%.1f%%)", ifelse(length(var_exp) >= 2, var_exp[2], NA_real_))
  )

pca_plot_file <- file.path(FIGURES_DIR, "pca_plot.jpeg")
ggsave(pca_plot_file, plot = pca_plot_12, width = 8, height = 6, dpi = 320)
message("[05_pca_dapc] Saved: ", pca_plot_file)

if (n_pc >= 3) {
  plot_13_df <- pca_scores %>% filter(is.finite(PC1), is.finite(PC3))
  plot_23_df <- pca_scores %>% filter(is.finite(PC2), is.finite(PC3))
  
  if (nrow(plot_13_df) > 0) {
    pca_plot_13 <- ggplot(plot_13_df, aes(PC1, PC3, color = Site)) +
      geom_point(size = 2.3, alpha = 0.9) +
      plot_theme +
      labs(
        title = "PCA (PC1 vs PC3)",
        x = sprintf("PC1 (%.1f%%)", var_exp[1]),
        y = sprintf("PC3 (%.1f%%)", var_exp[3])
      )
    
    pca13_file <- file.path(FIGURES_SUPP_DIR, "PC1_vs_PC3_plot.jpeg")
    ggsave(pca13_file, plot = pca_plot_13, width = 8, height = 6, dpi = 320)
    message("[05_pca_dapc] Saved: ", pca13_file)
  } else {
    message("[05_pca_dapc] Skipped PC1 vs PC3 plot (no finite rows).")
  }
  
  if (nrow(plot_23_df) > 0) {
    pca_plot_23 <- ggplot(plot_23_df, aes(PC2, PC3, color = Site)) +
      geom_point(size = 2.3, alpha = 0.9) +
      plot_theme +
      labs(
        title = "PCA (PC2 vs PC3)",
        x = sprintf("PC2 (%.1f%%)", var_exp[2]),
        y = sprintf("PC3 (%.1f%%)", var_exp[3])
      )
    
    pca23_file <- file.path(FIGURES_SUPP_DIR, "PC2_vs_PC3_plot.jpeg")
    ggsave(pca23_file, plot = pca_plot_23, width = 8, height = 6, dpi = 320)
    message("[05_pca_dapc] Saved: ", pca23_file)
  } else {
    message("[05_pca_dapc] Skipped PC2 vs PC3 plot (no finite rows).")
  }
} else {
  message("[05_pca_dapc] PC3 unavailable; supplementary PC3 plots skipped.")
}

# ------------------------------------------------------------
# 2) DAPC using Site as the a priori group
# ------------------------------------------------------------
# We keep the retained number of PCs moderate and explicit. For a thesis
# workflow, this is preferable to a hidden default because it makes the model
# choice reportable and reproducible. The number of discriminant axes is bound
# by (number of groups - 1).
set.seed(123)
max_n_pca <- min(50L, adegenet::nInd(gi_mll) - 1L, ncol(X))
if (max_n_pca < 1) {
  stop("[05_pca_dapc] Not enough informative dimensions to fit DAPC.")
}

suggested_n_pca <- min(
  max_n_pca,
  max(10L, floor(adegenet::nInd(gi_mll) / 3), nlevels(site_factor) + 1L)
)
suggested_n_pca <- max(1L, suggested_n_pca)
n_da <- min(nlevels(site_factor) - 1L, 10L)

if (n_da < 1) {
  stop("[05_pca_dapc] DAPC requires at least two Site groups with non-empty membership.")
}

message("[05_pca_dapc] DAPC grouping variable: Site")
message("[05_pca_dapc] Retained principal components (n.pca): ", suggested_n_pca)
message("[05_pca_dapc] Retained discriminant axes (n.da): ", n_da)

dapc_fit <- adegenet::dapc(
  x = gi_mll,
  pop = site_factor,
  n.pca = suggested_n_pca,
  n.da = n_da
)

dapc_coords <- as.data.frame(dapc_fit$ind.coord)
if (ncol(dapc_coords) < 1) {
  stop("[05_pca_dapc] DAPC returned no individual coordinates.")
}
if (ncol(dapc_coords) < 2) {
  dapc_coords$LD2 <- 0
}

colnames(dapc_coords)[1:2] <- c("LD1", "LD2")

dapc_coordinates <- data.frame(
  Individual = adegenet::indNames(gi_mll),
  Site = as.character(site_factor),
  LD1 = as.numeric(dapc_coords$LD1),
  LD2 = as.numeric(dapc_coords$LD2),
  stringsAsFactors = FALSE
) %>%
  filter(!is.na(Site), nzchar(Site), is.finite(LD1), is.finite(LD2))

if (nrow(dapc_coordinates) == 0) {
  stop("[05_pca_dapc] No valid rows available for DAPC outputs after filtering.")
}

dapc_coord_file <- file.path(TABLES_DIR, "dapc_coordinates.csv")
write.csv(dapc_coordinates, dapc_coord_file, row.names = FALSE)
message("[05_pca_dapc] Saved: ", dapc_coord_file)

dapc_centroids <- dapc_coordinates %>%
  group_by(Site) %>%
  summarise(
    n_clone_corrected = n(),
    LD1_centroid = mean(LD1, na.rm = TRUE),
    LD2_centroid = mean(LD2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Site)

centroids_file <- file.path(TABLES_DIR, "dapc_group_centroids.csv")
write.csv(dapc_centroids, centroids_file, row.names = FALSE)
message("[05_pca_dapc] Saved: ", centroids_file)

dapc_metadata <- data.frame(
  grouping_variable = "Site",
  genotype_object = "gi_mll",
  n_sites = nlevels(site_factor),
  n_individuals = adegenet::nInd(gi_mll),
  n_pca_retained = suggested_n_pca,
  n_da_retained = n_da,
  stringsAsFactors = FALSE
)

metadata_file <- file.path(TABLES_DIR, "dapc_model_metadata.csv")
write.csv(dapc_metadata, metadata_file, row.names = FALSE)
message("[05_pca_dapc] Saved: ", metadata_file)

ellipse_sites <- dapc_coordinates %>%
  count(Site, name = "n_site") %>%
  filter(n_site >= 3) %>%
  pull(Site)

if (length(ellipse_sites) == 0) {
  warning("[05_pca_dapc] No Site group has >=3 clone-corrected individuals; confidence ellipses cannot be drawn.")
}

# The DAPC plot uses site-level 95% ellipses where possible, which makes group
# overlap and outliers much easier to interpret in a publication figure.
dapc_plot <- ggplot(dapc_coordinates, aes(LD1, LD2, color = Site)) +
  {
    if (length(ellipse_sites) > 0) {
      stat_ellipse(
        data = subset(dapc_coordinates, Site %in% ellipse_sites),
        aes(fill = Site),
        geom = "polygon",
        alpha = 0.12,
        level = 0.95,
        linewidth = 0.35,
        show.legend = FALSE,
        type = "t"
      )
    }
  } +
  {
    if (length(ellipse_sites) > 0) {
      stat_ellipse(
        data = subset(dapc_coordinates, Site %in% ellipse_sites),
        linewidth = 0.55,
        level = 0.95,
        show.legend = FALSE,
        type = "t"
      )
    }
  } +
  geom_point(size = 2.6, alpha = 0.92, shape = 16) +
  geom_point(
    data = dapc_centroids,
    aes(x = LD1_centroid, y = LD2_centroid),
    inherit.aes = FALSE,
    shape = 4,
    stroke = 1,
    size = 3,
    color = "black"
  ) +
  plot_theme +
  labs(
    title = "DAPC by Site (clone-corrected microsatellite data)",
    subtitle = "Points are clone-corrected MLL representatives; ellipses show 95% site-level dispersion when n ≥ 3",
    x = "Discriminant axis 1",
    y = "Discriminant axis 2"
  )

dapc_plot_file <- file.path(FIGURES_DIR, "dapc_plot.jpeg")
ggsave(dapc_plot_file, plot = dapc_plot, width = 8.6, height = 6.4, dpi = 320)
message("[05_pca_dapc] Saved: ", dapc_plot_file)