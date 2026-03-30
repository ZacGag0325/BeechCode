# ============================================================================
# Publication-ready methodological review bar plots
# Project: beechcode
# ============================================================================

# 1) Set working directory to project root ------------------------------------
# IMPORTANT: adjust only if needed
setwd("~/Desktop/BeechCode")

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
df <- df %>%
  mutate(
    Method_Category = str_squish(as.character(Method_Category)),
    Assumed_or_Tested = str_squish(as.character(Assumed_or_Tested)),
    First_Author_Year = str_squish(as.character(First_Author_Year)),
    Stade_Development = str_squish(as.character(Stade_Development)),
    Sourced_from = str_squish(as.character(Sourced_from))
  )

# 6) Recode categories ---------------------------------------------------------
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
df$Method_Main <- factor(df$Method_Main,
                         levels = c("Excavation", "Genetic identification", "Spatial proximity", "Root connection", "Other")
)

# 8) Count data ----------------------------------------------------------------
summary_df <- df %>%
  count(Method_Main, Method_Sub)

print(summary_df)

# 9) Stacked bar plot ----------------------------------------------------------
p <- ggplot(summary_df, aes(x = Method_Main, y = n, fill = Method_Sub)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5), size = 3) +
  coord_flip() +
  labs(
    x = "Method category",
    y = "Number of studies",
    fill = "Method detail",
    title = "Methods used to determine regeneration origin"
  ) +
  theme_minimal()

print(p)

ggsave("figures/method_barplot_stacked.png", p, width = 8, height = 6, dpi = 300)

# 10) Assumed vs Tested --------------------------------------------------------
p2 <- df %>%
  count(Assumed_or_Tested) %>%
  ggplot(aes(x = reorder(Assumed_or_Tested, n), y = n)) +
  geom_bar(stat = "identity", fill = "#2C7FB8") +
  geom_text(aes(label = n), hjust = -0.2) +
  coord_flip() +
  theme_minimal() +
  labs(title = "Assumed vs Tested", x = "", y = "Count")

ggsave("figures/assumed_tested_barplot.png", p2, width = 8, height = 5)

# 11) Stade development --------------------------------------------------------
p3 <- df %>%
  count(Stade_Development) %>%
  ggplot(aes(x = reorder(Stade_Development, n), y = n)) +
  geom_bar(stat = "identity", fill = "#2C7FB8") +
  geom_text(aes(label = n), hjust = -0.2) +
  coord_flip() +
  theme_minimal() +
  labs(title = "Stage of development", x = "", y = "Count")

ggsave("figures/stade_development_barplot.png", p3, width = 8, height = 5)

# 12) Source -------------------------------------------------------------------
p4 <- df %>%
  count(Sourced_from) %>%
  ggplot(aes(x = reorder(Sourced_from, n), y = n)) +
  geom_bar(stat = "identity", fill = "#2C7FB8") +
  geom_text(aes(label = n), hjust = -0.2) +
  coord_flip() +
  theme_minimal() +
  labs(title = "Source of studies", x = "", y = "Count")

ggsave("figures/sourced_from_barplot.png", p4, width = 8, height = 5)

cat("All plots successfully created in /figures folder ✅\n")