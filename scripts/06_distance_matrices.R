# scripts/06_distance_matrices.R
############################################################
# Pairwise population differentiation matrices (gi_mll)
# IMPORTANT:
# - Jost's D is the primary differentiation metric for highly
#   polymorphic microsatellite loci.
# - FST is included as a secondary comparison metric.
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
  library(mmod)
  library(hierfstat)
  library(ggplot2)
  library(dplyr)
})

source("scripts/_load_objects.R")

message("[06_distance_matrices] Calculating Jost's D (primary) and FST (secondary)...")

# ----------------------------
# 1) Pairwise Jost's D
# ----------------------------
pairwise_jostD <- as.matrix(mmod::pairwise_D(gi_mll, linearized = FALSE))
diag(pairwise_jostD) <- NA_real_

jost_file <- file.path(MATRICES_DIR, "pairwise_jostD.csv")
write.csv(pairwise_jostD, jost_file)
message("[06_distance_matrices] Saved: ", jost_file)

jost_long_all <- as.data.frame(as.table(pairwise_jostD), stringsAsFactors = FALSE)
names(jost_long_all) <- c("Site1", "Site2", "JostD")
jost_long <- jost_long_all %>%
  filter(Site1 != Site2) %>%
  mutate(site_pair = paste(pmin(Site1, Site2), pmax(Site1, Site2), sep = "__")) %>%
  group_by(site_pair) %>%
  summarise(Site1 = first(pmin(Site1, Site2)), Site2 = first(pmax(Site1, Site2)), JostD = mean(JostD, na.rm = TRUE), .groups = "drop") %>%
  select(Site1, Site2, JostD)

jost_long_file <- file.path(TABLES_DIR, "pairwise_jostD_long.csv")
write.csv(jost_long, jost_long_file, row.names = FALSE)
message("[06_distance_matrices] Saved: ", jost_long_file)

# ----------------------------
# 2) Pairwise FST
# ----------------------------
hf <- hierfstat::genind2hierfstat(gi_mll)
pairwise_fst <- as.matrix(hierfstat::pairwise.WCfst(hf[, -1, drop = FALSE], hf[, 1]))
diag(pairwise_fst) <- NA_real_

fst_file <- file.path(MATRICES_DIR, "pairwise_fst.csv")
write.csv(pairwise_fst, fst_file)
message("[06_distance_matrices] Saved: ", fst_file)

fst_long_all <- as.data.frame(as.table(pairwise_fst), stringsAsFactors = FALSE)
names(fst_long_all) <- c("Site1", "Site2", "FST")
fst_long <- fst_long_all %>%
  filter(Site1 != Site2) %>%
  mutate(site_pair = paste(pmin(Site1, Site2), pmax(Site1, Site2), sep = "__")) %>%
  group_by(site_pair) %>%
  summarise(Site1 = first(pmin(Site1, Site2)), Site2 = first(pmax(Site1, Site2)), FST = mean(FST, na.rm = TRUE), .groups = "drop") %>%
  select(Site1, Site2, FST)

fst_long_file <- file.path(TABLES_DIR, "pairwise_fst_long.csv")
write.csv(fst_long, fst_long_file, row.names = FALSE)
message("[06_distance_matrices] Saved: ", fst_long_file)

# ----------------------------
# 3) Heatmaps (publication-ready)
# ----------------------------
heat_theme <- theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold")
  )

jost_heatmap <- ggplot(jost_long_all, aes(Site1, Site2, fill = JostD)) +
  geom_tile(color = "grey90", linewidth = 0.2) +
  scale_fill_viridis_c(option = "C", na.value = "white") +
  coord_equal() +
  heat_theme +
  labs(title = "Pairwise Jost's D (primary metric)", x = NULL, y = NULL, fill = "Jost's D")

jost_fig <- file.path(FIGURES_DIR, "pairwise_jostD_heatmap.jpeg")
ggsave(jost_fig, plot = jost_heatmap, width = 7.2, height = 6.2, dpi = 320)
message("[06_distance_matrices] Saved: ", jost_fig)

fst_heatmap <- ggplot(fst_long_all, aes(Site1, Site2, fill = FST)) +
  geom_tile(color = "grey90", linewidth = 0.2) +
  scale_fill_viridis_c(option = "B", na.value = "white") +
  coord_equal() +
  heat_theme +
  labs(title = "Pairwise FST (secondary metric)", x = NULL, y = NULL, fill = "FST")

fst_fig <- file.path(FIGURES_DIR, "pairwise_fst_heatmap.jpeg")
ggsave(fst_fig, plot = fst_heatmap, width = 7.2, height = 6.2, dpi = 320)
message("[06_distance_matrices] Saved: ", fst_fig)

# ----------------------------
# 4) Site-level differentiation summary
# ----------------------------
site_names <- rownames(pairwise_jostD)

site_differentiation_summary <- data.frame(
  Site = site_names,
  mean_pairwise_JostD = vapply(site_names, function(s) mean(pairwise_jostD[s, ], na.rm = TRUE), numeric(1)),
  mean_pairwise_FST = vapply(site_names, function(s) mean(pairwise_fst[s, ], na.rm = TRUE), numeric(1)),
  stringsAsFactors = FALSE
) %>% arrange(Site)

summary_file <- file.path(TABLES_DIR, "site_differentiation_summary.csv")
write.csv(site_differentiation_summary, summary_file, row.names = FALSE)
message("[06_distance_matrices] Saved: ", summary_file)