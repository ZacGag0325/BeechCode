# scripts/06a_bruvo_distance.R
############################################################
# Individual-level Bruvo distance (full dataset: gi)
# Why Bruvo here?
# - Bruvo's distance is appropriate for microsatellite loci and stepwise
#   mutation differences, so we use it for individual-level exploratory work.
# - This is supplementary and NOT the main site-level differentiation metric.
#
# Outputs:
# - outputs/matrices/bruvo_individual_distance.csv
# - outputs/tables/supplementary/bruvo_pairwise_long.csv
# - outputs/tables/supplementary/bruvo_near_clone_pairs.csv
# - outputs/figures/supplementary/bruvo_distance_distribution.jpeg
# - outputs/figures/supplementary/bruvo_distance_heatmap.jpeg
############################################################

suppressPackageStartupMessages({
  library(poppr)
  library(adegenet)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

source("scripts/_load_objects.R")

message("[06a_bruvo_distance] Computing individual-level Bruvo distances on gi...")

bd <- poppr::bruvo.dist(gi)
bd_mat <- as.matrix(bd)
storage.mode(bd_mat) <- "numeric"
diag(bd_mat) <- 0

if (!isTRUE(all.equal(bd_mat, t(bd_mat), tolerance = 1e-10, check.attributes = FALSE))) {
  stop("[06a_bruvo_distance] Bruvo distance matrix is not symmetric.")
}

bruvo_file <- file.path(MATRICES_DIR, "bruvo_individual_distance.csv")
write.csv(bd_mat, bruvo_file, row.names = TRUE)
message("[06a_bruvo_distance] Saved: ", bruvo_file)

upper <- upper.tri(bd_mat)
bruvo_long <- data.frame(
  Individual1 = rownames(bd_mat)[row(bd_mat)[upper]],
  Individual2 = colnames(bd_mat)[col(bd_mat)[upper]],
  BruvoDistance = as.numeric(bd_mat[upper]),
  stringsAsFactors = FALSE
)

pair_file <- file.path(TABLES_SUPP_DIR, "bruvo_pairwise_long.csv")
write.csv(bruvo_long, pair_file, row.names = FALSE)
message("[06a_bruvo_distance] Saved: ", pair_file)

# Near-clone diagnostics: flag very small Bruvo distances.
q01 <- as.numeric(stats::quantile(bruvo_long$BruvoDistance, probs = 0.01, na.rm = TRUE))
near_clone_cutoff <- max(0, min(0.05, q01))

near_clone_pairs <- bruvo_long %>%
  filter(is.finite(BruvoDistance), BruvoDistance <= near_clone_cutoff) %>%
  arrange(BruvoDistance)

near_clone_file <- file.path(TABLES_SUPP_DIR, "bruvo_near_clone_pairs.csv")
write.csv(near_clone_pairs, near_clone_file, row.names = FALSE)
message("[06a_bruvo_distance] Saved: ", near_clone_file)
message("[06a_bruvo_distance] Near-clone threshold used: ", signif(near_clone_cutoff, 4),
        " (min(0.05, 1st percentile))")

hist_plot <- ggplot(bruvo_long, aes(BruvoDistance)) +
  geom_histogram(bins = 40, fill = "#3E7CB1", color = "white") +
  geom_vline(xintercept = near_clone_cutoff, linetype = "dashed", color = "#C0392B", linewidth = 0.8) +
  theme_bw(base_size = 12) +
  labs(
    title = "Distribution of pairwise Bruvo distances (individual-level, supplementary)",
    x = "Bruvo distance",
    y = "Pair count"
  )

hist_file <- file.path(FIGURES_SUPP_DIR, "bruvo_distance_distribution.jpeg")
ggsave(hist_file, plot = hist_plot, width = 8, height = 6, dpi = 320)
message("[06a_bruvo_distance] Saved: ", hist_file)

heat_df <- as.data.frame(as.table(bd_mat), stringsAsFactors = FALSE)
names(heat_df) <- c("Individual1", "Individual2", "BruvoDistance")

heat_plot <- ggplot(heat_df, aes(Individual1, Individual2, fill = BruvoDistance)) +
  geom_tile() +
  scale_fill_viridis_c(option = "C") +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  ) +
  labs(
    title = "Bruvo distance heatmap (individual-level, supplementary)",
    x = "Individuals",
    y = "Individuals"
  )

heat_file <- file.path(FIGURES_SUPP_DIR, "bruvo_distance_heatmap.jpeg")
ggsave(heat_file, plot = heat_plot, width = 8, height = 7, dpi = 320)
message("[06a_bruvo_distance] Saved: ", heat_file)