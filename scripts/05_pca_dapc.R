############################################################
# scripts/05_pca_dapc.R
# PCA + DAPC (non-interactive) to compare with STRUCTURE
# Fix: handle K=2 (only 1 LD axis) without crashing
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

# Ensure pop is Site (consistent)
pop(gi) <- as.factor(df_ids$Site[match(indNames(gi), df_ids$ind)])

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

############################################################
# 1) PCA (allele freq table -> prcomp)
############################################################
cat("\n[PCA] Building PCA from allele frequency table...\n")

X <- adegenet::tab(gi, freq = TRUE, NA.method = "mean")  # individuals x alleles
X <- scale(X)

pca <- prcomp(X, center = TRUE, scale. = FALSE)

pc <- as.data.frame(pca$x[, 1:2, drop = FALSE])
pc$ind  <- rownames(pc)
pc$Site <- pop(gi)

keep_cols <- intersect(names(df_ids), c("ind", "Latitude", "Longitude"))
if (length(keep_cols) > 0) {
  pc <- merge(pc, df_ids[, keep_cols, drop = FALSE], by = "ind", all.x = TRUE)
}

p_pca12 <- ggplot(pc, aes(x = PC1, y = PC2, color = Site)) +
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

max_pca <- min(60, nInd(gi) - 1)
n_pca  <- min(30, max_pca)
n_da   <- min(nPop(gi) - 1, 10)

dapc_site <- adegenet::dapc(gi, pop(gi), n.pca = n_pca, n.da = n_da)

dapc_df <- get_ld_df(dapc_site)
dapc_df$Site <- pop(gi)[match(dapc_df$ind, indNames(gi))]

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
  p_dapc <- ggplot(dapc_df, aes(x = LD1, y = LD2, color = Site)) +
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
    gi,
    max.n.clust = 20,
    n.pca = N_PCA_FC,
    choose.n.clust = FALSE,
    n.clust = K_BEST
  ),
  error = function(e) {
    tryCatch(
      adegenet::find.clusters(
        gi,
        max.n.clust = 20,
        n.pca = N_PCA_FC,
        choose.n.clust = FALSE,
        k = K_BEST
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
  
  dapc_k <- adegenet::dapc(gi, fc$grp, n.pca = N_PCA_FC, n.da = n_da_k)
  
  dk <- get_ld_df(dapc_k)
  dk$Cluster <- as.factor(fc$grp)[match(dk$ind, indNames(gi))]
  dk$Site    <- pop(gi)[match(dk$ind, indNames(gi))]
  
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
    p_dapc_k <- ggplot(dk, aes(x = LD1, y = LD2, color = Cluster, shape = Site)) +
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
