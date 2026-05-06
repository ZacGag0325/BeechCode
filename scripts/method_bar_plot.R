# ============================================================================
# Publication-ready methodological review bar plots (EN + FR)
# Project: beechcode
# ============================================================================

# 1) Set working directory to project root ------------------------------------
# IMPORTANT: adjust only if needed
project_dir <- "~/Desktop/BeechCode"
if (dir.exists(path.expand(project_dir))) {
  setwd(path.expand(project_dir))
} else {
  message("Using current working directory: ", getwd())
}

# 2) Load libraries ------------------------------------------------------------
library(tidyverse)
library(readxl)
library(stringr)

# 3) Define file path ----------------------------------------------------------
input_file <- "data/raw/scopus_review_database.xlsx"
input_sheet <- "Data extraction"
output_dir <- "figures"

# Check file exists
if (!file.exists(input_file)) {
  stop("File not found at: ", input_file)
}

# Create figures folder if needed
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# 4) Read Excel ----------------------------------------------------------------
df <- read_excel(input_file, sheet = input_sheet)

# 5) Clean data ----------------------------------------------------------------
# Standardize spaces and convert key columns to character
df <- df %>%
  mutate(
    Method_Category = str_squish(as.character(Method_Category)),
    Assumed_or_Tested = str_squish(as.character(Assumed_or_Tested)),
    First_Author_Year = str_squish(as.character(First_Author_Year)),
    Stade_Development = str_squish(as.character(Stade_Development)),
    Sourced_from = str_squish(as.character(Sourced_from))
  )

# 6) Recode method categories --------------------------------------------------
df_method <- df %>%
  mutate(
    Method_Category_clean_lower = str_to_lower(str_squish(Method_Category)),
    Method_Category_clean = case_when(
      str_detect(Method_Category_clean_lower, "excavation") &
        str_detect(Method_Category_clean_lower, "morphologie") &
        str_detect(Method_Category_clean_lower, "lien racinaire") ~ "Excavation - Morphologie du collet et lien racinaire",
      str_detect(Method_Category_clean_lower, "excavation") &
        str_detect(Method_Category_clean_lower, "morphologie") ~ "Excavation - Morphologie du collet",
      str_detect(Method_Category_clean_lower, "excavation") &
        str_detect(Method_Category_clean_lower, "lien racinaire") ~ "Excavation - Lien racinaire entre individus",
      str_detect(Method_Category_clean_lower, "excavation") &
        str_detect(Method_Category_clean_lower, "non explicite") ~ "Excavation - Non explicite",
      str_detect(Method_Category_clean_lower, "lien racinaire") &
        str_detect(Method_Category_clean_lower, "racine de surface") ~ "Lien racinaire - Racine de surface",
      str_detect(Method_Category_clean_lower, "proximité|proximite") ~ "Proximité des individus",
      str_detect(Method_Category_clean_lower, "identification génétique|identification genetique|génétique|genetique") ~ "Identification génétique",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Method_Category_clean))

# 7) Define EN/FR dictionaries -------------------------------------------------
method_levels_en <- c(
  "Excavation - Morphologie du collet et lien racinaire",
  "Excavation - Morphologie du collet",
  "Excavation - Lien racinaire entre individus",
  "Excavation - Non explicite",
  "Lien racinaire - Racine de surface",
  "Proximité des individus",
  "Identification génétique"
)

method_labels_fr <- c(
  "Excavation - Morphologie du collet et lien racinaire" = "Excavation - Morphologie du collet et lien racinaire",
  "Excavation - Morphologie du collet" = "Excavation - Morphologie du collet",
  "Excavation - Lien racinaire entre individus" = "Excavation - Lien racinaire entre individus",
  "Excavation - Non explicite" = "Excavation - Non explicite",
  "Lien racinaire - Racine de surface" = "Lien racinaire - Racine de surface",
  "Proximité des individus" = "Proximité des individus",
  "Identification génétique" = "Identification génétique"
)

# Apply ordering for method categories
df_method <- df_method %>%
  mutate(Method_Category_clean = factor(Method_Category_clean, levels = method_levels_en))

# 8) Count data for method plot ------------------------------------------------
summary_df <- df_method %>%
  count(Method_Category_clean)

print(summary_df)

# Shared presentation settings for PNG exports.
# This style is matched to nearest_neighbor_sampling_check.R so all BeechCode
# presentation bar plots use consistent text sizing, ticks, and grid styling.
bar_fill <- "#2E8B57"
count_label_size <- 7

theme_pub <- theme_bw(base_size = 22) +
  theme(
    plot.title = element_text(size = 26, face = "bold"),
    axis.title = element_text(size = 26, face = "bold"),
    axis.text = element_text(size = 22),
    axis.text.x = element_text(size = 22),
    axis.text.y = element_text(size = 22),
    panel.grid.minor = element_blank(),
    # Enlarged tick marks for presentation readability
    axis.ticks = element_line(linewidth = 1.0),
    axis.ticks.length = grid::unit(0.28, "cm")
  )

# 9) Method stacked bar plot (EN + FR) ----------------------------------------
# English (main plot without legend)
p_method_en <- ggplot(summary_df, aes(x = Method_Category_clean, y = n)) +
  geom_col(fill = bar_fill, color = NA, width = 0.75) +
  geom_text(aes(label = n), hjust = -0.15, size = count_label_size, color = "black") +
  coord_flip() +
  expand_limits(y = max(summary_df$n, na.rm = TRUE) * 1.15 + 0.3) +
  labs(
    x = "Method category",
    y = "Number of studies",
    title = NULL
  ) +
  theme_pub

print(p_method_en)

ggsave(
  file.path(output_dir, "method_barplot_stacked.png"),
  p_method_en,
  width = 12,
  height = 8,
  dpi = 300
)

# Separate method legends are omitted because all method bars use the same presentation color.

# French (translate visible plotted labels)
summary_df_fr <- summary_df %>%
  mutate(
    Method_Category_clean_fr = recode(as.character(Method_Category_clean), !!!method_labels_fr)
  )

method_levels_fr <- unname(method_labels_fr[method_levels_en])

summary_df_fr <- summary_df_fr %>%
  mutate(
    Method_Category_clean_fr = factor(Method_Category_clean_fr, levels = method_levels_fr)
  )

p_method_fr <- ggplot(summary_df_fr, aes(x = Method_Category_clean_fr, y = n)) +
  geom_col(fill = bar_fill, color = NA, width = 0.75) +
  geom_text(aes(label = n), hjust = -0.15, size = count_label_size, color = "black") +
  coord_flip() +
  expand_limits(y = max(summary_df_fr$n, na.rm = TRUE) * 1.15 + 0.3) +
  labs(
    x = "Catégorie de méthode",
    y = "Nombre d’études",
    title = NULL
  ) +
  theme_pub

print(p_method_fr)

ggsave(
  file.path(output_dir, "method_barplot_stacked_fr.png"),
  p_method_fr,
  width = 12,
  height = 8,
  dpi = 300
)

# 10) Assumed vs Tested (EN + FR) ---------------------------------------------
assumed_tested_df <- df %>%
  count(Assumed_or_Tested) %>%
  filter(!is.na(Assumed_or_Tested), Assumed_or_Tested != "")

# English
p2_en <- ggplot(assumed_tested_df, aes(x = reorder(Assumed_or_Tested, n), y = n)) +
  geom_col(fill = bar_fill, color = NA, width = 0.75) +
  geom_text(aes(label = n), hjust = -0.15, size = count_label_size, color = "black") +
  coord_flip() +
  expand_limits(y = max(assumed_tested_df$n, na.rm = TRUE) * 1.15 + 0.3) +
  labs(title = NULL, x = "", y = "Count") +
  theme_pub

print(p2_en)

ggsave(
  file.path(output_dir, "assumed_tested_barplot.png"),
  p2_en,
  width = 10,
  height = 6,
  dpi = 300
)

# French
assumed_tested_df_fr <- assumed_tested_df %>%
  mutate(
    Assumed_or_Tested_fr = case_when(
      str_detect(str_to_lower(Assumed_or_Tested), "assum") ~ "Présumé",
      str_detect(str_to_lower(Assumed_or_Tested), "test") ~ "Vérifié",
      TRUE ~ Assumed_or_Tested
    )
  )

p2_fr <- ggplot(assumed_tested_df_fr, aes(x = reorder(Assumed_or_Tested_fr, n), y = n)) +
  geom_col(fill = bar_fill, color = NA, width = 0.75) +
  geom_text(aes(label = n), hjust = -0.15, size = count_label_size, color = "black") +
  coord_flip() +
  expand_limits(y = max(assumed_tested_df_fr$n, na.rm = TRUE) * 1.15 + 0.3) +
  labs(title = NULL, x = "", y = "Nombre") +
  theme_pub

print(p2_fr)

ggsave(
  file.path(output_dir, "assumed_tested_barplot_fr.png"),
  p2_fr,
  width = 10,
  height = 6,
  dpi = 300
)

# 11) Stade development --------------------------------------------------------
# Clean and standardize developmental stage values using pattern matching.
df <- df %>%
  mutate(
    Stade_Development_clean = str_to_lower(str_squish(Stade_Development)),
    Stade_Development_clean = str_replace_all(Stade_Development_clean, "[;|/]", ","),
    
    Stade_Development_std = case_when(
      is.na(Stade_Development_clean) | Stade_Development_clean == "" ~ "Unspecified",
      str_detect(Stade_Development_clean, "mixed stages|stades? mixtes|juvenile tree|arbre juvenile|\\bother\\b|\\bautre\\b") ~ NA_character_,
      str_detect(Stade_Development_clean, "^sapling\\s*(and|&)\\s*seedling|^saplings\\s*(and|&)\\s*seedlings|^gaule\\s*(et|&)\\s*semis|^gaule\\s*(et|&)\\s*semis") ~ "Sapling and Seedling",
      str_detect(Stade_Development_clean, "all trees|all stages|all developmental stages|all individuals|all size classes") ~ "All trees",
      str_detect(Stade_Development_clean, "^(na|n/a|nd|none)$|non renseign|non préc|non prec|not specified|not stated|unspecified|unknown") ~ "Unspecified",
      str_detect(Stade_Development_clean, "sapling|saplings|gaule|gaule") &
        str_detect(Stade_Development_clean, "mature|adult|arbre mature") ~ "Sapling and Mature",
      str_detect(Stade_Development_clean, "seedling|seedlings|semis") &
        str_detect(Stade_Development_clean, "sapling|saplings|gaule|gaule") ~ "Seedling and Sapling",
      str_detect(Stade_Development_clean, "seedling|seedlings|semis") ~ "Seedling",
      str_detect(Stade_Development_clean, "sapling|saplings|gaule|gaule") ~ "Sapling",
      str_detect(Stade_Development_clean, "mature|adult|arbre mature") |
        (str_detect(Stade_Development_clean, "tree|trees|arbre|arbres") &
           !str_detect(Stade_Development_clean, "young|juvenile")) ~ "Mature tree",
      TRUE ~ NA_character_
    )
  )

stage_order_en <- c(
  "Seedling",
  "Sapling",
  "Mature tree",
  "Seedling and Sapling",
  "Sapling and Seedling",
  "Sapling and Mature",
  "All trees",
  "Unspecified"
)

stage_labels_fr <- c(
  "Seedling" = "Semis",
  "Sapling" = "Gaule",
  "Mature tree" = "Arbre mature",
  "Seedling and Sapling" = "Semis et gaule",
  "Sapling and Seedling" = "Gaule et semis",
  "Sapling and Mature" = "Gaule et arbre mature",
  "All trees" = "Tous les arbres",
  "Unspecified" = "Non précisé"
)

df <- df %>%
  filter(!is.na(Stade_Development_std), Stade_Development_std %in% stage_order_en) %>%
  mutate(Stade_Development_std = factor(Stade_Development_std, levels = stage_order_en))

stage_summary <- df %>%
  count(Stade_Development_std, name = "n")

print(stage_summary)

# Development stage plot - English
p3_en <- ggplot(stage_summary, aes(x = Stade_Development_std, y = n)) +
  geom_col(fill = bar_fill, color = NA, width = 0.75) +
  geom_text(aes(label = n), hjust = -0.15, size = count_label_size, color = "black") +
  coord_flip() +
  expand_limits(y = max(stage_summary$n, na.rm = TRUE) * 1.15 + 0.3) +
  labs(
    title = NULL,
    x = "",
    y = "Number of studies"
  ) +
  theme_pub

print(p3_en)

ggsave(
  file.path(output_dir, "stade_development_barplot.png"),
  p3_en,
  width = 10.5,
  height = 7,
  dpi = 300
)

# Development stage plot - French (translated category labels in plotted data)
stage_summary_fr <- stage_summary %>%
  mutate(Stade_Development_std_fr = recode(as.character(Stade_Development_std), !!!stage_labels_fr))

stage_levels_fr <- unname(stage_labels_fr[stage_order_en])
stage_summary_fr <- stage_summary_fr %>%
  mutate(Stade_Development_std_fr = factor(Stade_Development_std_fr, levels = stage_levels_fr))

p3_fr <- ggplot(stage_summary_fr, aes(x = Stade_Development_std_fr, y = n)) +
  geom_col(fill = bar_fill, color = NA, width = 0.75) +
  geom_text(aes(label = n), hjust = -0.15, size = count_label_size, color = "black") +
  coord_flip() +
  expand_limits(y = max(stage_summary_fr$n, na.rm = TRUE) * 1.15 + 0.3) +
  labs(
    title = NULL,
    x = "",
    y = "Nombre d’études"
  ) +
  theme_pub

print(p3_fr)

ggsave(
  file.path(output_dir, "stade_development_barplot_fr.png"),
  p3_fr,
  width = 10.5,
  height = 7,
  dpi = 300
)

# 12) Source of studies (EN + FR) ---------------------------------------------
source_df <- df %>%
  count(Sourced_from) %>%
  filter(!is.na(Sourced_from), Sourced_from != "")

# English
p4_en <- ggplot(source_df, aes(x = reorder(Sourced_from, n), y = n)) +
  geom_col(fill = bar_fill, color = NA, width = 0.75) +
  geom_text(aes(label = n), hjust = -0.15, size = count_label_size, color = "black") +
  coord_flip() +
  expand_limits(y = max(source_df$n, na.rm = TRUE) * 1.15 + 0.3) +
  labs(title = NULL, x = "", y = "Count") +
  theme_pub

print(p4_en)

ggsave(
  file.path(output_dir, "sourced_from_barplot.png"),
  p4_en,
  width = 10.5,
  height = 6.5,
  dpi = 300
)

# French
p4_fr <- ggplot(source_df, aes(x = reorder(Sourced_from, n), y = n)) +
  geom_col(fill = bar_fill, color = NA, width = 0.75) +
  geom_text(aes(label = n), hjust = -0.15, size = count_label_size, color = "black") +
  coord_flip() +
  expand_limits(y = max(source_df$n, na.rm = TRUE) * 1.15 + 0.3) +
  labs(title = NULL, x = "", y = "Nombre") +
  theme_pub

print(p4_fr)

ggsave(
  file.path(output_dir, "sourced_from_barplot_fr.png"),
  p4_fr,
  width = 10.5,
  height = 6.5,
  dpi = 300
)

# 13) Final confirmation -------------------------------------------------------
cat("English and French plots were created successfully in the 'figures' folder ✅\n")