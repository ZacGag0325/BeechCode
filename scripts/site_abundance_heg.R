# ============================================================
# BEECH DENSITY BY SITE
# Full script from start to end
# ============================================================

# -----------------------------
# 1) Load packages
# -----------------------------
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(writexl)
})

# -----------------------------
# 2) User settings
# -----------------------------
# Where your raw data file is stored
# (based on your note: file is in data/raw)
data_dir <- "data/raw"

# Excel file name inside data/raw
input_filename <- "donnees_modifiees_west_summer2024 copie.xlsx"

# Full input path
file_path <- file.path(data_dir, input_filename)

# Beech species code in your dataset
# Change this if needed (example: "HEG", "FASY", etc.)
beech_code <- "HEG"

# Regeneration subplot area in m²
# IMPORTANT:
# Replace this with the real area used in your field protocol if different.
# Example:
# - 4 m² subplot  -> 4
# - 1 m² subplot  -> 1
# - 2 m radius circular plot -> pi * 2^2 = 12.56637
regen_plot_area_m2 <- 4

# Output file
output_excel <- "beech_density_by_site_results.xlsx"

# -----------------------------
# 3) Input checks
# -----------------------------
if (!file.exists(file_path)) {
  available_xlsx <- list.files(data_dir, pattern = "\\.xlsx$", full.names = FALSE)
  
  msg <- paste0(
    "Input file not found: ", normalizePath(file_path, winslash = "/", mustWork = FALSE), "\n",
    "\nFiles currently found in '", data_dir, "':\n",
    if (length(available_xlsx) == 0) "(none)" else paste0("- ", available_xlsx, collapse = "\n"),
    "\n\nUpdate 'input_filename' or move your file into ", data_dir, "."
  )
  
  stop(msg)
}

if (!is.numeric(regen_plot_area_m2) || regen_plot_area_m2 <= 0) {
  stop("'regen_plot_area_m2' must be a positive numeric value.")
}

# -----------------------------
# 4) Read sheets
# -----------------------------
arbre <- read_excel(file_path, sheet = "arbre")
gaule <- read_excel(file_path, sheet = "gaule")
regeneration <- read_excel(file_path, sheet = "regeneration")

# -----------------------------
# 5) Basic cleaning
# -----------------------------
# Make sure text columns are character, trimmed and standardized
clean_species <- function(df) {
  df %>%
    mutate(
      Site = str_trim(as.character(Site)),
      Espece = str_to_upper(str_trim(as.character(Espece)))
    )
}

arbre <- clean_species(arbre)
gaule <- clean_species(gaule)
regeneration <- clean_species(regeneration)

beech_code <- str_to_upper(str_trim(beech_code))

# -----------------------------
# 6) MATURE TREES (sheet: arbre)
# -----------------------------
# Assumption:
# Nb_tiges_ha is already the per-hectare contribution for each tree record.
# So for each site, we sum Nb_tiges_ha for beech.

beech_arbre <- arbre %>%
  filter(Espece == beech_code)

mature_density_site <- beech_arbre %>%
  group_by(Site) %>%
  summarise(
    n_records_beech = n(),
    mature_density_ha = sum(Nb_tiges_ha, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mature_density_ha))

# -----------------------------
# 7) GAULE (sheet: gaule)
# -----------------------------
# Assumption:
# Each row = one individual counted inside a circular plot of radius R_utilise (m).
# Density per point = count * 10000 / sampled_area
# sampled_area = pi * R_utilise^2
#
# We first count beech individuals per point, then convert to stems/ha,
# then average point densities within each site.

beech_gaule <- gaule %>%
  filter(Espece == beech_code)

gaule_point_counts <- beech_gaule %>%
  group_by(Site, Point_prisme, R_utilise) %>%
  summarise(
    n_beech_gaule = n(),
    .groups = "drop"
  ) %>%
  mutate(
    sampled_area_m2 = pi * (R_utilise^2),
    gaule_density_ha_point = (n_beech_gaule / sampled_area_m2) * 10000
  )

# Include points with zero beech in the site average
all_gaule_points <- gaule %>%
  distinct(Site, Point_prisme, R_utilise)

gaule_point_counts_complete <- all_gaule_points %>%
  left_join(gaule_point_counts, by = c("Site", "Point_prisme", "R_utilise")) %>%
  mutate(
    n_beech_gaule = replace_na(n_beech_gaule, 0),
    sampled_area_m2 = pi * (R_utilise^2),
    gaule_density_ha_point = replace_na(gaule_density_ha_point, 0)
  )

gaule_density_site <- gaule_point_counts_complete %>%
  group_by(Site) %>%
  summarise(
    n_points = n(),
    total_beech_gaule = sum(n_beech_gaule, na.rm = TRUE),
    mean_gaule_density_ha = mean(gaule_density_ha_point, na.rm = TRUE),
    sd_gaule_density_ha = sd(gaule_density_ha_point, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_gaule_density_ha))

# -----------------------------
# 8) REGENERATION (sheet: regeneration)
# -----------------------------
# Assumption:
# "Compte" = number of stems in each regeneration subplot.
# Density per subplot = count / subplot_area * 10000
# Then average subplot densities by site.

beech_regen <- regeneration %>%
  filter(Espece == beech_code)

regen_subplot <- beech_regen %>%
  group_by(Site, Id_spa) %>%
  summarise(
    regen_count_beech = sum(Compte, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    regen_density_ha_subplot = (regen_count_beech / regen_plot_area_m2) * 10000
  )

# Include subplots with zero beech in the site average
all_regen_subplots <- regeneration %>%
  distinct(Site, Id_spa)

regen_subplot_complete <- all_regen_subplots %>%
  left_join(regen_subplot, by = c("Site", "Id_spa")) %>%
  mutate(
    regen_count_beech = replace_na(regen_count_beech, 0),
    regen_density_ha_subplot = replace_na(regen_density_ha_subplot, 0)
  )

regen_density_site <- regen_subplot_complete %>%
  group_by(Site) %>%
  summarise(
    n_subplots = n(),
    total_beech_regen = sum(regen_count_beech, na.rm = TRUE),
    mean_regen_density_ha = mean(regen_density_ha_subplot, na.rm = TRUE),
    sd_regen_density_ha = sd(regen_density_ha_subplot, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_regen_density_ha))

# -----------------------------
# 9) Top site in each stratum
# -----------------------------
top_mature <- mature_density_site %>% slice(1)
top_gaule <- gaule_density_site %>% slice(1)
top_regen <- regen_density_site %>% slice(1)

cat("\n=====================================\n")
cat("TOP SITE FOR BEECH DENSITY BY STRATUM\n")
cat("=====================================\n\n")

cat("MATURE TREES:\n")
print(top_mature)

cat("\nGAULE:\n")
print(top_gaule)

cat("\nREGENERATION:\n")
print(top_regen)

# -----------------------------
# 10) Combined ranking table
# -----------------------------
combined_ranking <- mature_density_site %>%
  select(Site, mature_density_ha) %>%
  full_join(gaule_density_site %>% select(Site, mean_gaule_density_ha), by = "Site") %>%
  full_join(regen_density_site %>% select(Site, mean_regen_density_ha), by = "Site") %>%
  arrange(Site)

cat("\n========================\n")
cat("COMBINED SITE TABLE\n")
cat("========================\n\n")
print(combined_ranking)

# -----------------------------
# 11) Optional: site ranks
# -----------------------------
mature_ranked <- mature_density_site %>%
  mutate(rank_mature = row_number(desc(mature_density_ha)))

gaule_ranked <- gaule_density_site %>%
  mutate(rank_gaule = row_number(desc(mean_gaule_density_ha)))

regen_ranked <- regen_density_site %>%
  mutate(rank_regen = row_number(desc(mean_regen_density_ha)))

site_ranks <- combined_ranking %>%
  left_join(mature_ranked %>% select(Site, rank_mature), by = "Site") %>%
  left_join(gaule_ranked %>% select(Site, rank_gaule), by = "Site") %>%
  left_join(regen_ranked %>% select(Site, rank_regen), by = "Site")

cat("\n========================\n")
cat("SITE RANKS\n")
cat("========================\n\n")
print(site_ranks)

# -----------------------------
# 12) Save outputs to Excel
# -----------------------------
output_list <- list(
  mature_density_site = mature_density_site,
  gaule_density_site = gaule_density_site,
  regen_density_site = regen_density_site,
  combined_ranking = combined_ranking,
  site_ranks = site_ranks,
  gaule_point_counts_complete = gaule_point_counts_complete,
  regen_subplot_complete = regen_subplot_complete
)

write_xlsx(output_list, output_excel)

cat("\nResults saved to:\n")
cat(output_excel, "\n")

# -----------------------------
# 13) Optional quick interpretation
# -----------------------------
cat("\n=====================================\n")
cat("QUICK INTERPRETATION\n")
cat("=====================================\n\n")

cat(
  "Highest mature beech density site: ",
  top_mature$Site,
  " (", round(top_mature$mature_density_ha, 2), " stems/ha)\n",
  sep = ""
)

cat(
  "Highest gaule beech density site: ",
  top_gaule$Site,
  " (", round(top_gaule$mean_gaule_density_ha, 2), " stems/ha)\n",
  sep = ""
)

cat(
  "Highest regeneration beech density site: ",
  top_regen$Site,
  " (", round(top_regen$mean_regen_density_ha, 2), " stems/ha)\n",
  sep = ""
)

cat("\nDone.\n")