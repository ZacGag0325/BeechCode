# scripts/05_pca_dapc.R
############################################################
# PCA + DAPC (clone-corrected: gi_mll)
# Why gi_mll here?
# - Population-level structure methods assume independent genotypes.
# - Clone-corrected data reduce bias from repeated ramets.
# Outputs:
# - outputs/tables/pca_scores.csv
# - outputs/tables/pca_eigenvalues.csv
# - outputs/figures/pca_plot.jpeg (PC1 vs PC2)
# - outputs/figures/supplementary/PC1_vs_PC3_plot.jpeg
# - outputs/figures/supplementary/PC2_vs_PC3_plot.jpeg
# - outputs/tables/dapc_coordinates.csv
# - outputs/figures/dapc_plot.jpeg
############################################################

suppressPackageStartupMessages({
  library(adegenet)
  library(ggplot2)
  library(dplyr)
})

source("scripts/_load_objects.R")

message("[05_pca_dapc] Running PCA and DAPC on gi_mll...")

# ----------------------------
# 1) PCA
# ----------------------------
X <- adegenet::tab(gi_mll, freq = TRUE, NA.method = "mean")
keep_cols <- apply(X, 2, function(v) stats::var(v, na.rm = TRUE) > 0)
X <- X[, keep_cols, drop = FALSE]

pca_fit <- stats::prcomp(X, center = TRUE, scale. = FALSE)
var_exp <- summary(pca_fit)$importance[2, ] * 100

pca_scores <- data.frame(
  Individual = rownames(pca_fit$x),
  Site = as.character(adegenet::pop(gi_mll)),
  PC1 = pca_fit$x[, 1],
  PC2 = pca_fit$x[, 2],
  PC3 = if (ncol(pca_fit$x) >= 3) pca_fit$x[, 3] else NA_real_,
  stringsAsFactors = FALSE
)

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
    plot.title = element_text(face = "bold")
  )

pca_plot_12 <- ggplot(pca_scores, aes(PC1, PC2, color = Site)) +
  geom_point(size = 2.2, alpha = 0.9) +
  plot_theme +
  labs(
    title = "PCA (PC1 vs PC2)",
    x = sprintf("PC1 (%.1f%%)", var_exp[1]),
    y = sprintf("PC2 (%.1f%%)", var_exp[2])
  )

pca_plot_file <- file.path(FIGURES_DIR, "pca_plot.jpeg")
ggsave(pca_plot_file, plot = pca_plot_12, width = 8, height = 6, dpi = 320)
message("[05_pca_dapc] Saved: ", pca_plot_file)

if (!all(is.na(pca_scores$PC3))) {
  pca_plot_13 <- ggplot(pca_scores, aes(PC1, PC3, color = Site)) +
    geom_point(size = 2.2, alpha = 0.9) +
    plot_theme +
    labs(
      title = "PCA (PC1 vs PC3)",
      x = sprintf("PC1 (%.1f%%)", var_exp[1]),
      y = sprintf("PC3 (%.1f%%)", var_exp[3])
    )
  
  pca_plot_23 <- ggplot(pca_scores, aes(PC2, PC3, color = Site)) +
    geom_point(size = 2.2, alpha = 0.9) +
    plot_theme +
    labs(
      title = "PCA (PC2 vs PC3)",
      x = sprintf("PC2 (%.1f%%)", var_exp[2]),
      y = sprintf("PC3 (%.1f%%)", var_exp[3])
    )
  
  pca13_file <- file.path(FIGURES_SUPP_DIR, "PC1_vs_PC3_plot.jpeg")
  pca23_file <- file.path(FIGURES_SUPP_DIR, "PC2_vs_PC3_plot.jpeg")
  ggsave(pca13_file, plot = pca_plot_13, width = 8, height = 6, dpi = 320)
  ggsave(pca23_file, plot = pca_plot_23, width = 8, height = 6, dpi = 320)
  message("[05_pca_dapc] Saved: ", pca13_file)
  message("[05_pca_dapc] Saved: ", pca23_file)
}

# ----------------------------
# 2) DAPC
# ----------------------------
set.seed(123)
max_n_pca <- min(50, adegenet::nInd(gi_mll) - 1)
cluster_fit <- adegenet::find.clusters(
  gi_mll,
  max.n.clust = 20,
  n.pca = max_n_pca,
  choose.n.clust = FALSE
)
cluster_assign <- factor(cluster_fit$grp)

n_da <- min(length(unique(cluster_assign)) - 1, 10)
dapc_fit <- adegenet::dapc(gi_mll, pop = cluster_assign, n.pca = max_n_pca, n.da = n_da)

dapc_coords <- as.data.frame(dapc_fit$ind.coord)
if (ncol(dapc_coords) < 2) dapc_coords$LD2 <- 0

dapc_coordinates <- data.frame(
  Individual = adegenet::indNames(gi_mll),
  Site = as.character(adegenet::pop(gi_mll)),
  Cluster = as.character(cluster_assign),
  LD1 = dapc_coords[, 1],
  LD2 = dapc_coords[, 2],
  stringsAsFactors = FALSE
)

dapc_coord_file <- file.path(TABLES_DIR, "dapc_coordinates.csv")
write.csv(dapc_coordinates, dapc_coord_file, row.names = FALSE)
message("[05_pca_dapc] Saved: ", dapc_coord_file)

dapc_plot <- ggplot(dapc_coordinates, aes(LD1, LD2, color = Site, shape = Cluster)) +
  geom_point(size = 2.4, alpha = 0.9) +
  plot_theme +
  labs(
    title = "DAPC (clone-corrected microsatellite data)",
    x = "Discriminant axis 1",
    y = "Discriminant axis 2"
  )

dapc_plot_file <- file.path(FIGURES_DIR, "dapc_plot.jpeg")
ggsave(dapc_plot_file, plot = dapc_plot, width = 8, height = 6, dpi = 320)
message("[05_pca_dapc] Saved: ", dapc_plot_file)