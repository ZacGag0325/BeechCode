############################################################
# scripts/05_pca_dapc.R
# PCA + DAPC (non-interactive) to compare with STRUCTURE
# Fix: handle constant allele columns before PCA/DAPC
############################################################

options(repos = c(CRAN = "https://cloud.r-project.org"))

pkgs <- c("adegenet", "poppr", "ggplot2")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)

suppressPackageStartupMessages(library(adegenet))
suppressPackageStartupMessages(library(poppr))
suppressPackageStartupMessages(library(ggplot2))

source("scripts/_load_objects.R")

OUTDIR <- file.path(RUN_OUT, "pca_dapc_only")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# Prefer clone-corrected object when available
if (exists("gi_mll", inherits = FALSE)) {
  gi_use <- gi_mll
  gi_label <- "gi_mll"
} else {
  gi_use <- gi
  gi_label <- "gi"
}

# Ensure pop is Site (consistent)
pop(gi_use) <- as.factor(df_ids$Site[match(indNames(gi_use), df_ids$ind)])

save_plot <- function(p, filename, w = 7, h = 5) {
  ggsave(file.path(OUTDIR, filename), p, width = w, height = h, dpi = 200)
}

# Helper: return 1 or 2 LD columns safely
get_ld_df <- function(dapc_obj) {
  M <- dapc_obj$ind.coord
  if (is.null(M)) stop("dapc object has no ind.coord")
  M <- as.matrix(M)
  
  # ensure column names
  if (is.null(colnames(M))) {
    colnames(M) <- paste0("LD", seq_len(ncol(M)))
  }
  
  if (ncol(M) == 0) stop("dapc ind.coord has 0 columns")
  if (ncol(M) == 1) {
    out <- data.frame(LD1 = M[, 1], LD2 = NA_real_, stringsAsFactors = FALSE)
  } else {
    out <- data.frame(LD1 = M[, 1], LD2 = M[, 2], stringsAsFactors = FALSE)
  }
  out$ind <- rownames(M)
  out
}


# Helper: subset rows to groups with enough non-missing 2D points for ellipses
ellipse_subset <- function(df, group_col, x_col, y_col, min_n = 3) {
  keep <- !is.na(df[[group_col]]) & !is.na(df[[x_col]]) & !is.na(df[[y_col]])
  df_ok <- df[keep, , drop = FALSE]
  if (nrow(df_ok) == 0) return(df_ok)
  
  grp_counts <- table(df_ok[[group_col]])
  valid_groups <- names(grp_counts)[grp_counts >= min_n]
  df_ok[df_ok[[group_col]] %in% valid_groups, , drop = FALSE]
}

############################################################
# 1) PCA (allele table -> prcomp)
############################################################
cat("\n[PCA] Building PCA from allele table using ", gi_label, "...\n", sep = "")

X <- adegenet::tab(gi_use, NA.method = "mean")

# Remove constant columns (zero variance) because prcomp cannot rescale them.
# This directly fixes: "cannot rescale a constant/zero column to unit variance".
keep_cols <- apply(X, 2, function(col) stats::var(col, na.rm = TRUE) > 0)
message("Removing ", sum(!keep_cols), " constant columns out of ", length(keep_cols))
X <- X[, keep_cols, drop = FALSE]

write.csv(
  data.frame(removed_cols = sum(!keep_cols)),
  file.path(RUN_OUT, "pca_removed_constant_columns.csv"),
  row.names = FALSE
)

# Impute any residual NA with column means for stable PCA/DAPC computations.
if (anyNA(X)) {
  message("Imputing remaining NA values with column means")
  for (j in seq_len(ncol(X))) {
    if (anyNA(X[, j])) {
      X[is.na(X[, j]), j] <- mean(X[, j], na.rm = TRUE)
    }
  }
}

stopifnot(ncol(X) > 1)

pca_res <- prcomp(X, center = TRUE, scale. = FALSE)
message("PCA completed: ", nrow(X), " individuals × ", ncol(X), " variables")

pc <- as.data.frame(pca_res$x[, 1:2, drop = FALSE])
pc$ind  <- rownames(pc)
pc$Site <- pop(gi_use)

keep_meta_cols <- intersect(names(df_ids), c("ind", "Latitude", "Longitude"))
if (length(keep_meta_cols) > 0) {
  pc <- merge(pc, df_ids[, keep_meta_cols, drop = FALSE], by = "ind", all.x = TRUE)
}

pc_ell <- ellipse_subset(pc, group_col = "Site", x_col = "PC1", y_col = "PC2", min_n = 3)

p_pca12 <- ggplot(pc, aes(x = PC1, y = PC2, color = Site)) +
  stat_ellipse(
    data = pc_ell,
    aes(x = PC1, y = PC2, group = Site, fill = Site),
    type = "norm", level = 0.95,
    geom = "polygon", alpha = 0.25,
    color = NA, inherit.aes = FALSE
  ) +
  stat_ellipse(
    data = pc_ell,
    aes(x = PC1, y = PC2, group = Site, color = Site),
    type = "norm", level = 0.95,
    linewidth = 0.9, inherit.aes = FALSE
  ) +
  geom_point(size = 2, alpha = 0.85) +
  theme_minimal(base_size = 12) +
  labs(title = "PCA (PC1 vs PC2)", subtitle = "Individuals colored by Site")

save_plot(p_pca12, "PCA_PC1_PC2_bySite.png", w = 7.5, h = 5.5)
ggsave(file.path(OUTDIR, "PCA_PC1_PC2_bySite.pdf"), p_pca12, width = 7.5, height = 5.5)

if ("Latitude" %in% names(pc) && any(!is.na(pc$Latitude))) {
  p_pca_lat <- ggplot(pc, aes(x = Latitude, y = PC1, color = Site)) +
    geom_point(size = 2, alpha = 0.85) +
    theme_minimal(base_size = 12) +
    labs(title = "PCA: PC1 vs Latitude", subtitle = "Check cline / IBD pattern")
  
  save_plot(p_pca_lat, "PCA_PC1_vs_Latitude.png", w = 7.5, h = 5.5)
  ggsave(file.path(OUTDIR, "PCA_PC1_vs_Latitude.pdf"), p_pca_lat, width = 7.5, height = 5.5)
}

write.csv(pc, file.path(OUTDIR, "PCA_scores_PC1_PC2.csv"), row.names = FALSE)

############################################################
# 2) DAPC supervised by Site
############################################################
cat("\n[DAPC] Running DAPC supervised by Site...\n")

grp_site <- as.factor(pop(gi_use))
if (any(is.na(grp_site))) {
  stop("DAPC grouping factor contains NA values in pop(gi_use).")
}

max_pca <- min(60, nrow(X) - 1, ncol(X) - 1)
n_pca <- min(30, max_pca)
n_da <- min(nlevels(grp_site) - 1, 10)

dapc_site <- adegenet::dapc(X, grp_site, n.pca = n_pca, n.da = n_da, scale = FALSE)

dapc_df <- get_ld_df(dapc_site)
dapc_df$Site <- grp_site[match(dapc_df$ind, rownames(X))]

# If LD2 exists, do 2D; otherwise do 1D
if (all(is.na(dapc_df$LD2))) {
  p_dapc <- ggplot(dapc_df, aes(x = Site, y = LD1, color = Site)) +
    geom_jitter(width = 0.15, height = 0, size = 2, alpha = 0.85) +
    theme_minimal(base_size = 12) +
    labs(
      title = "DAPC (supervised by Site)",
      subtitle = paste0("Only LD1 available | n.pca=", n_pca, " | n.da=", n_da),
      x = NULL, y = "LD1"
    )
} else {
  dapc_ell <- ellipse_subset(dapc_df, group_col = "Site", x_col = "LD1", y_col = "LD2", min_n = 3)
  
  p_dapc <- ggplot(dapc_df, aes(x = LD1, y = LD2, color = Site)) +
    stat_ellipse(
      data = dapc_ell,
      aes(x = LD1, y = LD2, group = Site, fill = Site),
      type = "norm", level = 0.95,
      geom = "polygon", alpha = 0.25,
      color = NA, inherit.aes = FALSE
    ) +
    stat_ellipse(
      data = dapc_ell,
      aes(x = LD1, y = LD2, group = Site, color = Site),
      type = "norm", level = 0.95,
      linewidth = 0.9, inherit.aes = FALSE
    ) +
    geom_point(size = 2, alpha = 0.85) +
    theme_minimal(base_size = 12) +
    labs(
      title = "DAPC (supervised by Site)",
      subtitle = paste0("n.pca=", n_pca, " | n.da=", n_da)
    )
}

save_plot(p_dapc, "DAPC_supervised_bySite.png", w = 7.5, h = 5.5)
ggsave(file.path(OUTDIR, "DAPC_supervised_bySite.pdf"), p_dapc, width = 7.5, height = 5.5)

write.csv(dapc_df, file.path(OUTDIR, "DAPC_supervised_ind_coords.csv"), row.names = FALSE)
if (!is.null(dapc_site$posterior)) {
  post <- as.data.frame(dapc_site$posterior)
  post$ind <- rownames(post)
  write.csv(post, file.path(OUTDIR, "DAPC_supervised_posterior.csv"), row.names = FALSE)
}

############################################################
# 3) Unsupervised DAPC with K=2 (no prompts)
############################################################
K_BEST <- 2
N_PCA_FC <- min(30, max_pca)

cat("\n[DAPC] Unsupervised find.clusters with K_BEST=", K_BEST, "...\n", sep = "")

fc <- tryCatch(
  adegenet::find.clusters(
    X,
    max.n.clust = 20,
    n.pca = N_PCA_FC,
    choose.n.clust = FALSE,
    n.clust = K_BEST,
    scale = FALSE
  ),
  error = function(e) {
    tryCatch(
      adegenet::find.clusters(
        X,
        max.n.clust = 20,
        n.pca = N_PCA_FC,
        choose.n.clust = FALSE,
        k = K_BEST,
        scale = FALSE
      ),
      error = function(e2) {
        message("[DAPC] find.clusters failed (signature mismatch). Skipping unsupervised DAPC.\n  ",
                conditionMessage(e), "\n  ", conditionMessage(e2))
        NULL
      }
    )
  }
)

if (!is.null(fc)) {
  # For K=2, n.da max is 1
  n_da_k <- min(K_BEST - 1, 10)
  
  dapc_k <- adegenet::dapc(X, fc$grp, n.pca = N_PCA_FC, n.da = n_da_k, scale = FALSE)
  
  dk <- get_ld_df(dapc_k)
  dk$Cluster <- as.factor(fc$grp)[match(dk$ind, rownames(X))]
  dk$Site <- grp_site[match(dk$ind, rownames(X))]
  
  # Plot: if only LD1 exists (K=2), make a 1D jitter plot
  if (all(is.na(dk$LD2))) {
    p_dapc_k <- ggplot(dk, aes(x = Cluster, y = LD1, color = Site)) +
      geom_jitter(width = 0.15, height = 0, size = 2, alpha = 0.85) +
      theme_minimal(base_size = 12) +
      labs(
        title = paste0("DAPC (unsupervised clusters, K=", K_BEST, ")"),
        subtitle = paste0("Only LD1 available | find.clusters n.pca=", N_PCA_FC),
        x = "Cluster", y = "LD1"
      )
  } else {
    dk_ell <- ellipse_subset(dk, group_col = "Cluster", x_col = "LD1", y_col = "LD2", min_n = 3)
    
    p_dapc_k <- ggplot(dk, aes(x = LD1, y = LD2, color = Cluster, shape = Site)) +
      stat_ellipse(
        data = dk_ell,
        aes(x = LD1, y = LD2, group = Cluster, fill = Cluster),
        type = "norm", level = 0.95,
        geom = "polygon", alpha = 0.25,
        color = NA, inherit.aes = FALSE
      ) +
      stat_ellipse(
        data = dk_ell,
        aes(x = LD1, y = LD2, group = Cluster, color = Cluster),
        type = "norm", level = 0.95,
        linewidth = 0.9, inherit.aes = FALSE
      ) +
      geom_point(size = 2, alpha = 0.85) +
      theme_minimal(base_size = 12) +
      labs(
        title = paste0("DAPC (unsupervised clusters, K=", K_BEST, ")"),
        subtitle = paste0("find.clusters n.pca=", N_PCA_FC)
      )
  }
  
  save_plot(p_dapc_k, paste0("DAPC_unsupervised_K", K_BEST, ".png"), w = 8.5, h = 6)
  ggsave(
    file.path(OUTDIR, paste0("DAPC_unsupervised_K", K_BEST, ".pdf")),
    p_dapc_k, width = 8.5, height = 6
  )
  
  write.csv(dk, file.path(OUTDIR, paste0("DAPC_unsupervised_K", K_BEST, "_ind_coords.csv")), row.names = FALSE)
  
  tab_cs <- with(dk, table(Site, Cluster))
  write.csv(as.data.frame.matrix(tab_cs),
            file.path(OUTDIR, paste0("DAPC_unsupervised_K", K_BEST, "_cluster_counts_bySite.csv")))
}

cat("\nDONE PCA + DAPC. Outputs in: ", OUTDIR, "\n", sep = "")
