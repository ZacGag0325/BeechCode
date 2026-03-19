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
#   used throughout the population-level summaries (HWE, differentiation,
#   AMOVA, and site-level IBD). Using Site here keeps DAPC directly comparable
#   to the rest of the pipeline.
# - STRUCTURE remains a separate, externally run analysis on the full dataset.
#
# Biological interpretation notes:
# - PCA is unsupervised. It summarizes multivariate genetic variation without
#   using site labels to build axes.
# - Site-level confidence ellipses in PCA are added ONLY as a visualization
#   aid to help readers see overlap among sampling localities; they do not
#   imply that discrete clusters truly exist.
# - DAPC is supervised by the Site grouping variable. It is therefore useful
#   for describing separation among known sampling localities, but it should
#   be interpreted alongside unsupervised PCA and the distance-based analyses.
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
# - outputs/tables/dapc_assignment_summary.csv
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

plot_theme <- theme_bw(base_size = 12) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

eligible_ellipse_groups <- function(df, x_col, y_col, group_col = "Site") {
  df %>%
    group_by(.data[[group_col]]) %>%
    summarise(
      n_group = n(),
      var_x = stats::var(.data[[x_col]], na.rm = TRUE),
      var_y = stats::var(.data[[y_col]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(
      n_group >= 3,
      is.finite(var_x),
      is.finite(var_y),
      var_x > 0,
      var_y > 0
    ) %>%
    pull(.data[[group_col]])
}

# ------------------------------------------------------------
# 1) PCA (unsupervised)
# ------------------------------------------------------------
# adegenet::tab converts the genind object into allele-frequency / coded-count
# data for multivariate analysis. Missing values are replaced by mean allele
# frequencies so that individuals and loci are not discarded unnecessarily.
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

plot_12_df <- pca_scores %>% filter(is.finite(PC1), is.finite(PC2))
if (nrow(plot_12_df) == 0) {
  stop("[05_pca_dapc] No valid rows for PC1 vs PC2 plot.")
}

pca_ellipse_sites <- eligible_ellipse_groups(plot_12_df, x_col = "PC1", y_col = "PC2", group_col = "Site")
if (length(pca_ellipse_sites) == 0) {
  warning("[05_pca_dapc] No Site group has enough spread for PCA confidence ellipses; PCA plot will show points only.")
}

# PCA remains unsupervised. Ellipses are added beneath the points with low
# alpha so they help readability without obscuring observations.
pca_plot_12 <- ggplot(plot_12_df, aes(PC1, PC2, color = Site)) +
  {
    if (length(pca_ellipse_sites) > 0) {
      stat_ellipse(
        data = subset(plot_12_df, Site %in% pca_ellipse_sites),
        aes(fill = Site),
        geom = "polygon",
        alpha = 0.12,
        level = 0.95,
        linewidth = 0.25,
        show.legend = FALSE,
        type = "t"
      )
    }
  } +
  {
    if (length(pca_ellipse_sites) > 0) {
      stat_ellipse(
        data = subset(plot_12_df, Site %in% pca_ellipse_sites),
        linewidth = 0.5,
        alpha = 0.8,
        level = 0.95,
        show.legend = FALSE,
        type = "t"
      )
    }
  } +
  geom_point(size = 2.4, alpha = 0.92) +
  plot_theme +
  labs(
    title = "PCA (PC1 vs PC2)",
    subtitle = "PCA is unsupervised; ellipses are for visualization only and do not imply discrete clusters",
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
        subtitle = "PCA is unsupervised; points only are shown for supplementary axes",
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
        subtitle = "PCA is unsupervised; points only are shown for supplementary axes",
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
# We retain the existing Site-group DAPC, but choose n.pca more defensibly.
# First, a broad candidate model is fit. Then optim.a.score is attempted to
# identify a conservative number of retained PCs. If optimisation is unstable,
# the workflow falls back to a reproducible variance-based rule.
set.seed(123)
max_n_pca <- min(50L, adegenet::nInd(gi_mll) - 1L, ncol(X))
if (max_n_pca < 1) {
  stop("[05_pca_dapc] Not enough informative dimensions to fit DAPC.")
}

n_da_max <- min(nlevels(site_factor) - 1L, nrow(X) - 1L, max_n_pca)
if (n_da_max < 1) {
  stop("[05_pca_dapc] DAPC requires at least one valid discriminant axis.")
}

n_da_final <- min(n_da_max, 10L)
if (n_da_final < 1) {
  stop("[05_pca_dapc] DAPC requires at least two Site groups with non-empty membership.")
}

cumvar_threshold_idx <- which(cumsum(pca_eigenvalues$Percent_Variance) >= 80)
fallback_n_pca <- min(
  max_n_pca,
  max(
    n_da_final + 1L,
    nlevels(site_factor) + 1L,
    ifelse(length(cumvar_threshold_idx) > 0, min(cumvar_threshold_idx), max_n_pca)
  )
)
fallback_n_pca <- max(1L, as.integer(fallback_n_pca))

message("[05_pca_dapc] DAPC grouping variable: Site")
message("[05_pca_dapc] Candidate maximum retained principal components for optimisation: ", max_n_pca)
message("[05_pca_dapc] Final discriminant axes retained (n.da): ", n_da_final)

initial_dapc <- adegenet::dapc(
  x = gi_mll,
  pop = site_factor,
  n.pca = max_n_pca,
  n.da = n_da_final,
  pca.loadings = TRUE,
  var.contrib = FALSE,
  var.loadings = FALSE
)

optim_used <- FALSE
optim_best_n_pca <- NA_integer_
optim_note <- "optim.a.score_not_attempted"
optim_error_message <- NA_character_

optim_result <- tryCatch({
  adegenet::optim.a.score(
    initial_dapc,
    n.pca = seq_len(max_n_pca),
    smart = TRUE,
    plot = FALSE,
    n = 25,
    n.sim = 30,
    n.da = n_da_final
  )
}, error = function(e) {
  optim_error_message <<- conditionMessage(e)
  NULL
})

if (!is.null(optim_result) && !is.null(optim_result$best)) {
  optim_candidate <- suppressWarnings(as.integer(round(optim_result$best[1])))
  if (length(optim_candidate) == 1 && is.finite(optim_candidate) && !is.na(optim_candidate)) {
    optim_best_n_pca <- max(1L, min(max_n_pca, optim_candidate))
    optim_used <- TRUE
    optim_note <- "optim.a.score_success"
  } else {
    optim_note <- "optim.a.score_returned_non_finite_best"
  }
} else if (!is.null(optim_error_message) && nzchar(optim_error_message)) {
  optim_note <- paste0("optim.a.score_failed: ", optim_error_message)
}

n_pca_final <- if (optim_used && !is.na(optim_best_n_pca)) {
  max(n_da_final + 1L, optim_best_n_pca)
} else {
  fallback_n_pca
}
n_pca_final <- min(max_n_pca, as.integer(n_pca_final))
n_pca_final <- max(1L, n_pca_final)

if (n_pca_final <= n_da_final) {
  n_pca_final <- min(max_n_pca, n_da_final + 1L)
}
if (n_pca_final <= n_da_final) {
  n_da_final <- max(1L, n_pca_final - 1L)
}
if (n_da_final > nlevels(site_factor) - 1L) {
  n_da_final <- nlevels(site_factor) - 1L
}

if (optim_used) {
  n_pca_justification <- paste0(
    "n.pca selected from optim.a.score (best=", optim_best_n_pca,
    ") and constrained to remain > n.da for stable DAPC interpretation."
  )
} else {
  n_pca_justification <- paste0(
    "n.pca selected by reproducible fallback rule: retain at least the number required to exceed n.da and preserve >=80% cumulative PCA variance when possible (fallback=", fallback_n_pca, ")."
  )
}

message("[05_pca_dapc] Final retained principal components (n.pca): ", n_pca_final)
message("[05_pca_dapc] Final retained discriminant axes (n.da): ", n_da_final)
message("[05_pca_dapc] n.pca selection note: ", optim_note)
message("[05_pca_dapc] n.pca justification: ", n_pca_justification)

dapc_fit <- adegenet::dapc(
  x = gi_mll,
  pop = site_factor,
  n.pca = n_pca_final,
  n.da = n_da_final,
  pca.loadings = TRUE,
  var.contrib = FALSE,
  var.loadings = FALSE
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

posterior_df <- if (!is.null(dapc_fit$posterior)) {
  as.data.frame(dapc_fit$posterior, stringsAsFactors = FALSE)
} else {
  data.frame(matrix(nrow = adegenet::nInd(gi_mll), ncol = 0))
}

assignment_summary <- data.frame(
  Individual = adegenet::indNames(gi_mll),
  observed_site = as.character(site_factor),
  assigned_site = if (!is.null(dapc_fit$assign)) as.character(dapc_fit$assign) else NA_character_,
  max_posterior = if (ncol(posterior_df) > 0) apply(as.matrix(posterior_df), 1, max, na.rm = TRUE) else NA_real_,
  stringsAsFactors = FALSE
)

if (ncol(posterior_df) > 0) {
  colnames(posterior_df) <- paste0("posterior_", make.names(colnames(posterior_df), unique = TRUE))
  assignment_summary <- cbind(assignment_summary, posterior_df)
}

assignment_file <- file.path(TABLES_DIR, "dapc_assignment_summary.csv")
write.csv(assignment_summary, assignment_file, row.names = FALSE)
message("[05_pca_dapc] Saved: ", assignment_file)

constraint_notes <- c(
  paste0("max_n_pca=min(50, nInd-1, informative_allele_columns) => ", max_n_pca),
  paste0("n.da capped at min(n_groups-1, n_individuals-1, max_n_pca, 10) => ", n_da_final),
  if (optim_used) {
    paste0("optim.a.score_best=", optim_best_n_pca, "; final_n_pca adjusted to remain > n.da if needed")
  } else {
    paste0("fallback_n_pca=", fallback_n_pca, " based on cumulative_PCA_variance>=80% and group-count safeguards")
  }
)

dapc_metadata <- data.frame(
  grouping_variable = "Site",
  genotype_object = "gi_mll",
  n_sites = nlevels(site_factor),
  n_individuals = adegenet::nInd(gi_mll),
  n_pca_retained = n_pca_final,
  n_da_retained = n_da_final,
  optim_a_score_used = optim_used,
  optim_a_score_best_n_pca = ifelse(is.na(optim_best_n_pca), NA_integer_, optim_best_n_pca),
  optimisation_note = optim_note,
  n_pca_justification = n_pca_justification,
  constraint_notes = paste(constraint_notes, collapse = " | "),
  stringsAsFactors = FALSE
)

metadata_file <- file.path(TABLES_DIR, "dapc_model_metadata.csv")
write.csv(dapc_metadata, metadata_file, row.names = FALSE)
message("[05_pca_dapc] Saved: ", metadata_file)

dapc_ellipse_sites <- eligible_ellipse_groups(dapc_coordinates, x_col = "LD1", y_col = "LD2", group_col = "Site")
if (length(dapc_ellipse_sites) == 0) {
  warning("[05_pca_dapc] No Site group has enough spread for DAPC confidence ellipses; DAPC plot will show points and centroids only.")
}

# The DAPC plot shows supervised separation among Site groups. Filled 95%
# confidence ellipses summarize within-site dispersion and are drawn below the
# points. Centroids are added to make the position of each site's average DAPC
# score explicit for publication-ready interpretation.
dapc_plot <- ggplot(dapc_coordinates, aes(LD1, LD2, color = Site)) +
  {
    if (length(dapc_ellipse_sites) > 0) {
      stat_ellipse(
        data = subset(dapc_coordinates, Site %in% dapc_ellipse_sites),
        aes(fill = Site),
        geom = "polygon",
        alpha = 0.14,
        level = 0.95,
        linewidth = 0.35,
        show.legend = FALSE,
        type = "t"
      )
    }
  } +
  {
    if (length(dapc_ellipse_sites) > 0) {
      stat_ellipse(
        data = subset(dapc_coordinates, Site %in% dapc_ellipse_sites),
        linewidth = 0.7,
        level = 0.95,
        show.legend = FALSE,
        type = "t"
      )
    }
  } +
  geom_point(size = 2.6, alpha = 0.9, shape = 16) +
  geom_point(
    data = dapc_centroids,
    aes(x = LD1_centroid, y = LD2_centroid),
    inherit.aes = FALSE,
    shape = 4,
    stroke = 1.2,
    size = 3.4,
    color = "black"
  ) +
  geom_text(
    data = dapc_centroids,
    aes(x = LD1_centroid, y = LD2_centroid, label = Site),
    inherit.aes = FALSE,
    color = "black",
    size = 3,
    fontface = "bold",
    vjust = -0.7,
    check_overlap = TRUE
  ) +
  plot_theme +
  labs(
    title = "DAPC by Site (clone-corrected microsatellite data)",
    subtitle = paste0(
      "Points are clone-corrected MLL representatives; 95% confidence ellipses summarize within-site dispersion when estimable. ",
      "Centroids are labelled to support comparison among sites."
    ),
    x = "Discriminant axis 1",
    y = "Discriminant axis 2"
  )

dapc_plot_file <- file.path(FIGURES_DIR, "dapc_plot.jpeg")
ggsave(dapc_plot_file, plot = dapc_plot, width = 8.4, height = 6.4, dpi = 320)
message("[05_pca_dapc] Saved: ", dapc_plot_file)