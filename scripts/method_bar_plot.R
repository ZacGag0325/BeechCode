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
# Keep original logic: Excavation pooled as one main category,
# while preserving stacked subcategories.
df <- df %>%
  mutate(
    Method_Main = case_when(
      str_detect(Method_Category, "Excavation") ~ "Excavation",
      str_detect(Method_Category, "Lien racinaire") &
        !str_detect(Method_Category, "Excavation") ~ "Root connection",
      str_detect(Method_Category, "Proximité") ~ "Spatial proximity",
      str_detect(Method_Category, "Identification génétique") ~ "Genetic identification",
      TRUE ~ "Other"
    ),
    
    Method_Sub = case_when(
      str_detect(Method_Category, "Morphologie du collet") &
        str_detect(Method_Category, "Lien racinaire") ~ "Morphology + root link",
      str_detect(Method_Category, "Morphologie du collet") ~ "Morphology",
      str_detect(Method_Category, "Lien racinaire") ~ "Root link",
      str_detect(Method_Category, "Non explicite") ~ "Not specified",
      TRUE ~ "Other"
    )
  )

# Force subcategory for non-excavation
df <- df %>%
  mutate(Method_Sub = ifelse(Method_Main != "Excavation", Method_Main, Method_Sub))

# 7) Define EN/FR dictionaries -------------------------------------------------
method_main_levels_en <- c(
  "Excavation",
  "Genetic identification",
  "Spatial proximity",
  "Root connection",
  "Other"
)

method_main_labels_fr <- c(
  "Excavation" = "Excavation",
  "Genetic identification" = "Identification génétique",
  "Spatial proximity" = "Proximité spatiale",
  "Root connection" = "Connexion racinaire",
  "Other" = "Autre"
)

method_sub_labels_fr <- c(
  "Morphology + root link" = "Morphologie + lien racinaire",
  "Morphology" = "Morphologie",
  "Root link" = "Lien racinaire",
  "Not specified" = "Non précisé",
  "Excavation" = "Excavation",
  "Genetic identification" = "Identification génétique",
  "Spatial proximity" = "Proximité spatiale",
  "Root connection" = "Connexion racinaire",
  "Other" = "Autre"
)

# Apply ordering for method main categories
df <- df %>%
  mutate(Method_Main = factor(Method_Main, levels = method_main_levels_en))

# 8) Count data for method plot ------------------------------------------------
summary_df <- df %>%
  count(Method_Main, Method_Sub)

print(summary_df)

# 9) Method stacked bar plot (EN + FR) ----------------------------------------
# English
p_method_en <- ggplot(summary_df, aes(x = Method_Main, y = n, fill = Method_Sub)) +
  geom_col(width = 0.75) +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5), size = 3.2, color = "black") +
  coord_flip() +
  labs(
    x = "Method category",
    y = "Number of studies",
    fill = "Method detail",
    title = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

print(p_method_en)

ggsave(
  file.path(output_dir, "method_barplot_stacked.png"),
  p_method_en,
  width = 8,
  height = 6,
  dpi = 300
)

# French (translate visible plotted labels)
summary_df_fr <- summary_df %>%
  mutate(
    Method_Main_fr = recode(as.character(Method_Main), !!!method_main_labels_fr),
    Method_Sub_fr = recode(as.character(Method_Sub), !!!method_sub_labels_fr)
  )

method_main_levels_fr <- unname(method_main_labels_fr[method_main_levels_en])

summary_df_fr <- summary_df_fr %>%
  mutate(
    Method_Main_fr = factor(Method_Main_fr, levels = method_main_levels_fr),
    Method_Sub_fr = factor(Method_Sub_fr)
  )

p_method_fr <- ggplot(summary_df_fr, aes(x = Method_Main_fr, y = n, fill = Method_Sub_fr)) +
  geom_col(width = 0.75) +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5), size = 3.2, color = "black") +
  coord_flip() +
  labs(
    x = "Catégorie de méthode",
    y = "Nombre d’études",
    fill = "Détail de méthode",
    title = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

print(p_method_fr)

ggsave(
  file.path(output_dir, "method_barplot_stacked_fr.png"),
  p_method_fr,
  width = 8,
  height = 6,
  dpi = 300
)

# 10) Assumed vs Tested (EN + FR) ---------------------------------------------
assumed_tested_df <- df %>%
  count(Assumed_or_Tested) %>%
  filter(!is.na(Assumed_or_Tested), Assumed_or_Tested != "")

# English
p2_en <- ggplot(assumed_tested_df, aes(x = reorder(Assumed_or_Tested, n), y = n)) +
  geom_col(fill = "#2C7FB8", width = 0.75) +
  geom_text(aes(label = n), hjust = -0.15, size = 3.2) +
  coord_flip() +
  labs(title = NULL, x = "", y = "Count") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

print(p2_en)

ggsave(
  file.path(output_dir, "assumed_tested_barplot.png"),
  p2_en,
  width = 8,
  height = 5,
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
  geom_col(fill = "#2C7FB8", width = 0.75) +
  geom_text(aes(label = n), hjust = -0.15, size = 3.2) +
  coord_flip() +
  labs(title = NULL, x = "", y = "Nombre") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

print(p2_fr)

ggsave(
  file.path(output_dir, "assumed_tested_barplot_fr.png"),
  p2_fr,
  width = 8,
  height = 5,
  dpi = 300
)

# 11) Stade development --------------------------------------------------------
# Clean and standardize developmental stage values using pattern matching.
df <- df %>%
  mutate(
    Stade_Development_clean = str_to_lower(str_squish(Stade_Development)),
    Stade_Development_clean = str_replace_all(Stade_Development_clean, "[;|/]", ","),
    
    Stade_Development_std = case_when(
      # Missing / unspecified
      is.na(Stade_Development_clean) |
        Stade_Development_clean == "" |
        str_detect(Stade_Development_clean, "^na$|^n/a$|^nd$|non renseign|non préc|not specified|not stated|unspecified|unknown|none") ~ "Not specified",
      
      # Explicitly mixed / multiple stages
      str_detect(Stade_Development_clean, "mixed|multiple|plusieurs|all stages|all developmental|diverse") |
        str_detect(Stade_Development_clean, "seedling") & str_detect(Stade_Development_clean, "sapling|juvenile|mature|adult|tree") |
        str_detect(Stade_Development_clean, "sapling|gaulis") & str_detect(Stade_Development_clean, "juvenile|mature|adult|tree") |
        str_detect(Stade_Development_clean, "juvenile|young tree|young stem") & str_detect(Stade_Development_clean, "mature|adult|tree|arbre mature") |
        str_detect(Stade_Development_clean, "et|and|,") &
        str_detect(Stade_Development_clean, "seedling|semis|sapling|gaulis|juvenile|young tree|young stem|mature|adult|arbre mature|tree") ~ "Mixed stages",
      
      # Seedlings
      str_detect(Stade_Development_clean, "seedling|seedlings|semis") ~ "Seedling",
      
      # Saplings
      str_detect(Stade_Development_clean, "sapling|saplings|gaulis") ~ "Sapling",
      
      # Juvenile trees
      str_detect(Stade_Development_clean, "juvenile|young tree|young trees|young stem|young stems") ~ "Juvenile tree",
      
      # Mature trees
      str_detect(Stade_Development_clean, "mature|adult|arbre mature") |
        (str_detect(Stade_Development_clean, "tree|trees|arbre|arbres") &
           !str_detect(Stade_Development_clean, "young|juvenile|seedling|sapling|gaulis|semis")) ~ "Mature tree",
      
      TRUE ~ "Other"
    )
  )

stage_order_en <- c(
  "Seedling",
  "Sapling",
  "Juvenile tree",
  "Mature tree",
  "Mixed stages",
  "Not specified",
  "Other"
)

stage_labels_fr <- c(
  "Seedling" = "Semis",
  "Sapling" = "Gaulis",
  "Juvenile tree" = "Jeune arbre",
  "Mature tree" = "Arbre mature",
  "Mixed stages" = "Stades mixtes",
  "Not specified" = "Non précisé",
  "Other" = "Autre"
)

df <- df %>%
  mutate(Stade_Development_std = factor(Stade_Development_std, levels = stage_order_en))

stage_summary <- df %>%
  count(Stade_Development_std, name = "n") %>%
  complete(Stade_Development_std = factor(stage_order_en, levels = stage_order_en), fill = list(n = 0))

print(stage_summary)

# Development stage plot - English
p3_en <- ggplot(stage_summary, aes(x = Stade_Development_std, y = n)) +
  geom_col(fill = "#2C7FB8", width = 0.75) +
  geom_text(aes(label = n), hjust = -0.15, size = 3.3) +
  coord_flip() +
  expand_limits(y = max(stage_summary$n, na.rm = TRUE) * 1.12 + 0.2) +
  labs(
    title = NULL,
    x = "",
    y = "Number of studies"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

print(p3_en)

ggsave(
  file.path(output_dir, "stade_development_barplot.png"),
  p3_en,
  width = 8.5,
  height = 5.5,
  dpi = 300
)

# Development stage plot - French (translated category labels in plotted data)
stage_summary_fr <- stage_summary %>%
  mutate(Stade_Development_std_fr = recode(as.character(Stade_Development_std), !!!stage_labels_fr))

stage_levels_fr <- unname(stage_labels_fr[stage_order_en])
stage_summary_fr <- stage_summary_fr %>%
  mutate(Stade_Development_std_fr = factor(Stade_Development_std_fr, levels = stage_levels_fr))

p3_fr <- ggplot(stage_summary_fr, aes(x = Stade_Development_std_fr, y = n)) +
  geom_col(fill = "#2C7FB8", width = 0.75) +
  geom_text(aes(label = n), hjust = -0.15, size = 3.3) +
  coord_flip() +
  expand_limits(y = max(stage_summary_fr$n, na.rm = TRUE) * 1.12 + 0.2) +
  labs(
    title = NULL,
    x = "",
    y = "Nombre d’études"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

print(p3_fr)

ggsave(
  file.path(output_dir, "stade_development_barplot_fr.png"),
  p3_fr,
  width = 8.5,
  height = 5.5,
  dpi = 300
)

# 12) Source of studies (EN + FR) ---------------------------------------------
source_df <- df %>%
  count(Sourced_from) %>%
  filter(!is.na(Sourced_from), Sourced_from != "")

# English
p4_en <- ggplot(source_df, aes(x = reorder(Sourced_from, n), y = n)) +
  geom_col(fill = "#2C7FB8", width = 0.75) +
  geom_text(aes(label = n), hjust = -0.15, size = 3.2) +
  coord_flip() +
  labs(title = NULL, x = "", y = "Count") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

print(p4_en)

ggsave(
  file.path(output_dir, "sourced_from_barplot.png"),
  p4_en,
  width = 8,
  height = 5,
  dpi = 300
)

# French
p4_fr <- ggplot(source_df, aes(x = reorder(Sourced_from, n), y = n)) +
  geom_col(fill = "#2C7FB8", width = 0.75) +
  geom_text(aes(label = n), hjust = -0.15, size = 3.2) +
  coord_flip() +
  labs(title = NULL, x = "", y = "Nombre") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

print(p4_fr)

ggsave(
  file.path(output_dir, "sourced_from_barplot_fr.png"),
  p4_fr,
  width = 8,
  height = 5,
  dpi = 300
)

# 13) Final confirmation -------------------------------------------------------
cat("English and French plots were created successfully in the 'figures' folder ✅\n")