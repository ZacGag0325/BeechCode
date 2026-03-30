# ============================================================================
# Publication-ready methodological review bar plots
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

# 6) Recode categories ---------------------------------------------------------
# Keep your original logic: Excavation as one pooled main category,
# while still showing stacked subcategories in Method_Sub.
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
  mutate(
    Method_Sub = ifelse(Method_Main != "Excavation", Method_Main, Method_Sub)
  )

# 7) Order categories ----------------------------------------------------------
df$Method_Main <- factor(
  df$Method_Main,
  levels = c("Excavation", "Genetic identification", "Spatial proximity", "Root connection", "Other")
)

# 8) Count data ----------------------------------------------------------------
summary_df <- df %>%
  count(Method_Main, Method_Sub)

print(summary_df)

# 9) Stacked bar plot ----------------------------------------------------------
# Publication-ready stacked plot for methods
p <- ggplot(summary_df, aes(x = Method_Main, y = n, fill = Method_Sub)) +
  geom_col(width = 0.75) +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5), size = 3.2, color = "black") +
  coord_flip() +
  labs(
    x = "Method category",
    y = "Number of studies",
    fill = "Method detail",
    title = "Methods used to determine regeneration origin"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

print(p)

ggsave(file.path(output_dir, "method_barplot_stacked.png"), p, width = 8, height = 6, dpi = 300)

# 10) Assumed vs Tested --------------------------------------------------------
p2 <- df %>%
  count(Assumed_or_Tested) %>%
  filter(!is.na(Assumed_or_Tested), Assumed_or_Tested != "") %>%
  ggplot(aes(x = reorder(Assumed_or_Tested, n), y = n)) +
  geom_col(fill = "#2C7FB8", width = 0.75) +
  geom_text(aes(label = n), hjust = -0.15, size = 3.2) +
  coord_flip() +
  labs(title = "Assumed vs Tested", x = "", y = "Count") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

print(p2)

ggsave(file.path(output_dir, "assumed_tested_barplot.png"), p2, width = 8, height = 5, dpi = 300)

# 11) Stade development --------------------------------------------------------
# Clean and standardize developmental stage values using pattern matching.
# Goal categories:
# - Seedling
# - Sapling
# - Juvenile tree
# - Mature tree
# - Mixed stages
# - Not specified
# - Other

df <- df %>%
  mutate(
    Stade_Development_clean = str_to_lower(str_squish(Stade_Development)),
    
    # Normalize common separators and punctuation for robust matching
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

stage_order <- c(
  "Seedling",
  "Sapling",
  "Juvenile tree",
  "Mature tree",
  "Mixed stages",
  "Not specified",
  "Other"
)

df <- df %>%
  mutate(Stade_Development_std = factor(Stade_Development_std, levels = stage_order))

stage_summary <- df %>%
  count(Stade_Development_std, name = "n") %>%
  complete(Stade_Development_std = factor(stage_order, levels = stage_order), fill = list(n = 0))

print(stage_summary)

# Horizontal publication-ready bar plot for developmental stage
p3 <- ggplot(stage_summary, aes(x = Stade_Development_std, y = n)) +
  geom_col(fill = "#2C7FB8", width = 0.75) +
  geom_text(aes(label = n), hjust = -0.15, size = 3.3) +
  coord_flip() +
  expand_limits(y = max(stage_summary$n, na.rm = TRUE) * 1.12 + 0.2) +
  labs(
    title = "Developmental stage of trees/regeneration used in studies",
    x = "",
    y = "Number of studies"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

print(p3)

ggsave(file.path(output_dir, "stade_development_barplot.png"), p3, width = 8.5, height = 5.5, dpi = 300)

# 12) Source -------------------------------------------------------------------
p4 <- df %>%
  count(Sourced_from) %>%
  filter(!is.na(Sourced_from), Sourced_from != "") %>%
  ggplot(aes(x = reorder(Sourced_from, n), y = n)) +
  geom_col(fill = "#2C7FB8", width = 0.75) +
  geom_text(aes(label = n), hjust = -0.15, size = 3.2) +
  coord_flip() +
  labs(title = "Source of studies", x = "", y = "Count") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

print(p4)

ggsave(file.path(output_dir, "sourced_from_barplot.png"), p4, width = 8, height = 5, dpi = 300)

# 13) Final confirmation -------------------------------------------------------
cat("All plots were created successfully in the 'figures' folder ✅\n")