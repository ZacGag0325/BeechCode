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
site_lookup_sheet <- "site_lookup"
output_dir <- "figures"

# Check file exists
if (!file.exists(input_file)) {
  stop("File not found at: ", input_file)
}

# Create figures folder if needed
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# 4) Site lookup helpers -------------------------------------------------------
read_site_lookup <- function(path_in, sheet = "site_lookup") {
  sheets <- readxl::excel_sheets(path_in)
  if (!(sheet %in% sheets)) {
    stop(
      paste0(
        "Required site lookup sheet '", sheet, "' was not found in workbook: ", path_in, "\n",
        "Available sheets: ", paste(sheets, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
  lookup_raw <- read_excel(path_in, sheet = sheet)
  required_cols <- c("Site", "Site_label", "Region", "Site_order")
  missing_cols <- setdiff(required_cols, names(lookup_raw))
  if (length(missing_cols) > 0) {
    stop(
      paste0(
        "Sheet '", sheet, "' is missing required columns: ",
        paste(missing_cols, collapse = ", "), "\n",
        "Required columns are: ", paste(required_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
  lookup <- lookup_raw %>%
    transmute(
      Site = str_squish(as.character(.data$Site)),
      Site_label = str_squish(as.character(.data$Site_label)),
      Region = str_squish(as.character(.data$Region)),
      Site_order = suppressWarnings(as.numeric(.data$Site_order))
    ) %>%
    filter(!is.na(Site), Site != "")
  
  duplicated_sites <- lookup %>%
    count(Site, name = "n") %>%
    filter(n > 1)
  if (nrow(duplicated_sites) > 0) {
    stop(
      paste0(
        "Duplicate Site values found in sheet '", sheet, "': ",
        paste(duplicated_sites$Site, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
  invalid_lookup <- lookup %>%
    filter(is.na(Site_label) | Site_label == "" | is.na(Region) | Region == "" | is.na(Site_order))
  if (nrow(invalid_lookup) > 0) {
    stop(
      paste0(
        "Site lookup rows must have non-missing Site_label, Region, and numeric Site_order. ",
        "Problem Site values: ", paste(invalid_lookup$Site, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
  lookup %>%
    arrange(Site_order, Site_label)
}

attach_site_lookup <- function(data, lookup, data_name = "data frame") {
  if (!("Site" %in% names(data))) {
    return(data)
  }
  
  data <- data %>%
    mutate(Site = str_squish(as.character(.data$Site)))
  
  missing_sites <- setdiff(
    sort(unique(data$Site[!is.na(data$Site) & data$Site != ""])),
    lookup$Site
  )
  if (length(missing_sites) > 0) {
    stop(
      paste0(
        "The following Site value(s) in ", data_name, " are missing from site_lookup: ",
        paste(missing_sites, collapse = ", "), "\n",
        "Add them to the '", site_lookup_sheet, "' sheet with Site_label, Region, and Site_order."
      ),
      call. = FALSE
    )
  }
  
  site_label_levels <- lookup %>%
    arrange(Site_order, Site_label) %>%
    pull(Site_label)
  
  data %>%
    select(-any_of(c("Site_label", "Region", "Site_order"))) %>%
    left_join(lookup, by = "Site") %>%
    mutate(Site_label = factor(Site_label, levels = site_label_levels)) %>%
    arrange(Site_order, Site_label)
}

# 5) Read Excel ----------------------------------------------------------------
df <- read_excel(input_file, sheet = input_sheet)
site_lookup <- read_site_lookup(input_file, site_lookup_sheet)

# 6) Clean data ----------------------------------------------------------------
# Standardize spaces and convert key columns to character
df <- df %>%
  mutate(
    Method_Category = str_squish(as.character(Method_Category)),
    Assumed_or_Tested = str_squish(as.character(Assumed_or_Tested)),
    First_Author_Year = str_squish(as.character(First_Author_Year)),
    Stade_Development = str_squish(as.character(Stade_Development)),
    Sourced_from = str_squish(as.character(Sourced_from))
  ) %>%
  attach_site_lookup(site_lookup, "Data extraction sheet")

# 7) Recode method categories --------------------------------------------------
# Raw method labels are normalized into stable internal keys first. Plot-specific
# English and French labels are then applied from dictionaries so factor levels
# cannot mix languages.
method_dictionary <- tribble(
  ~method_key, ~method_label_en, ~method_label_fr,
  "excavation_collar_root", "Excavation + collar morphology + root connection", "Excavation + morphologie du collet + lien racinaire",
  "excavation_collar", "Excavation + collar morphology", "Excavation + morphologie du collet",
  "excavation_root", "Excavation + root connection between individuals", "Excavation + lien racinaire entre individus",
  "excavation_unspecified", "Excavation, method not specified", "Excavation, méthode non explicite",
  "root_surface", "Root connection / surface root", "Lien racinaire / racine de surface",
  "spatial_proximity", "Spatial proximity between individuals", "Proximité entre individus",
  "genetic_identification", "Genetic identification", "Identification génétique"
)

method_levels <- method_dictionary$method_key
method_labels_en <- setNames(method_dictionary$method_label_en, method_dictionary$method_key)
method_labels_fr <- setNames(method_dictionary$method_label_fr, method_dictionary$method_key)

df_method <- df %>%
  mutate(
    Method_Category_clean_lower = str_to_lower(str_squish(Method_Category)),
    Method_Category_key = case_when(
      str_detect(Method_Category_clean_lower, "excavation") &
        str_detect(Method_Category_clean_lower, "morphologie|morphology|collet|collar") &
        str_detect(Method_Category_clean_lower, "lien racinaire|root connection|root link|racinaire") ~ "excavation_collar_root",
      str_detect(Method_Category_clean_lower, "excavation") &
        str_detect(Method_Category_clean_lower, "morphologie|morphology|collet|collar") ~ "excavation_collar",
      str_detect(Method_Category_clean_lower, "excavation") &
        str_detect(Method_Category_clean_lower, "lien racinaire|root connection|root link|racinaire") ~ "excavation_root",
      str_detect(Method_Category_clean_lower, "excavation") &
        str_detect(Method_Category_clean_lower, "non explicite|not specified|unspecified|method not specified") ~ "excavation_unspecified",
      str_detect(Method_Category_clean_lower, "lien racinaire|root connection|root link|racinaire") &
        str_detect(Method_Category_clean_lower, "racine de surface|surface root") ~ "root_surface",
      str_detect(Method_Category_clean_lower, "proximité|proximite|proximity|spatial") ~ "spatial_proximity",
      str_detect(Method_Category_clean_lower, "identification génétique|identification genetique|génétique|genetique|genetic") ~ "genetic_identification",
      TRUE ~ NA_character_
    ),
    Method_Category_key = factor(Method_Category_key, levels = method_levels)
  ) %>%
  filter(!is.na(Method_Category_key)) %>%
  attach_site_lookup(site_lookup, "method category data")

# 9) Count data for method plot ------------------------------------------------
summary_df <- df_method %>%
  count(Method_Category_key, name = "n")

print(summary_df)

# Shared presentation settings for PNG exports.
# This style is matched to nearest_neighbor_sampling_check.R so all BeechCode
# presentation bar plots use consistent text sizing, ticks, and grid styling.
bar_fill <- "#2E8B57"
count_label_size <- 7

theme_pub <- theme_classic(base_size = 22) +
  theme(
    plot.title = element_text(size = 26, face = "bold"),
    axis.title = element_text(size = 26, face = "bold"),
    axis.text = element_text(size = 22),
    axis.text.x = element_text(size = 22),
    axis.text.y = element_text(size = 22),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    plot.background = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(linewidth = 0.4, color = "black"),
    # Enlarged tick marks for presentation readability
    axis.ticks = element_line(linewidth = 1.0),
    axis.ticks.length = grid::unit(0.28, "cm")
  )

# 10) Method stacked bar plot (EN + FR) ----------------------------------------
summary_df_en <- summary_df %>%
  mutate(
    Method_Category_label = recode(as.character(Method_Category_key), !!!method_labels_en),
    Method_Category_label = factor(Method_Category_label, levels = unname(method_labels_en[method_levels]))
  )

# English (main plot without legend)
p_method_en <- ggplot(summary_df_en, aes(x = Method_Category_label, y = n)) +
  geom_col(fill = bar_fill, color = NA, width = 0.75) +
  geom_text(aes(label = n), hjust = -0.15, size = count_label_size, color = "black") +
  coord_flip() +
  expand_limits(y = max(summary_df_en$n, na.rm = TRUE) * 1.15 + 0.3) +
  labs(
    x = "Method used to identify clones",
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
  dpi = 300,
  bg = "white"
)

# Separate method legends are omitted because all method bars use the same presentation color.

summary_df_fr <- summary_df %>%
  mutate(
    Method_Category_label = recode(as.character(Method_Category_key), !!!method_labels_fr),
    Method_Category_label = factor(Method_Category_label, levels = unname(method_labels_fr[method_levels]))
  )

p_method_fr <- ggplot(summary_df_fr, aes(x = Method_Category_label, y = n)) +
  geom_col(fill = bar_fill, color = NA, width = 0.75) +
  geom_text(aes(label = n), hjust = -0.15, size = count_label_size, color = "black") +
  coord_flip() +
  expand_limits(y = max(summary_df_fr$n, na.rm = TRUE) * 1.15 + 0.3) +
  labs(
    x = "Méthode utilisée pour identifier les clones",
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
  dpi = 300,
  bg = "white"
)

# 11) Assumed vs Tested (EN + FR) ---------------------------------------------
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
  dpi = 300,
  bg = "white"
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
  dpi = 300,
  bg = "white"
)

# 12) Stade development --------------------------------------------------------
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
  dpi = 300,
  bg = "white"
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
  dpi = 300,
  bg = "white"
)

# 13) Source of studies (EN + FR) ---------------------------------------------
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
  dpi = 300,
  bg = "white"
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
  dpi = 300,
  bg = "white"
)

# 14) Final confirmation -------------------------------------------------------
cat("English and French plots were created successfully in the 'figures' folder ✅\n")