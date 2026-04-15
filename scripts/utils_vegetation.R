# Utility functions for vegetation / stand structure descriptive analyses

suppressPackageStartupMessages({
  library(dplyr)
  library(janitor)
  library(readxl)
  library(stringr)
  library(tidyr)
  library(purrr)
})

find_workbook_path <- function(preferred_path = "donnees_modifiees_west_summer2024 copie.xlsx") {
  candidate_paths <- c(
    preferred_path,
    file.path("data", preferred_path),
    file.path("data", "raw", preferred_path),
    file.path("data", "input", preferred_path)
  )
  
  candidate_paths <- unique(candidate_paths[file.exists(candidate_paths)])
  
  if (length(candidate_paths) > 0) {
    return(candidate_paths[[1]])
  }
  
  discovered <- list.files(
    path = ".",
    pattern = "donnees_modifiees_west_summer2024.*\\.xlsx$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )
  
  if (length(discovered) > 0) {
    return(discovered[[1]])
  }
  
  return(NA_character_)
}

safe_read_sheet <- function(path, sheet) {
  out <- tryCatch(
    readxl::read_excel(path, sheet = sheet),
    error = function(e) {
      message(sprintf("Skipping sheet '%s' (read error): %s", sheet, e$message))
      return(NULL)
    }
  )
  
  if (is.null(out)) {
    return(NULL)
  }
  
  out %>%
    janitor::clean_names() %>%
    mutate(across(where(is.character), ~ str_squish(.x)))
}

normalize_species <- function(x) {
  x_std <- x %>%
    as.character() %>%
    str_squish() %>%
    str_to_lower()
  
  dplyr::case_when(
    is.na(x_std) ~ NA_character_,
    x_std %in% c("fagus grandifolia", "fagus_grandifolia", "american beech", "hetre", "hêtre", "hetre americain", "hêtre américain") ~ "Fagus grandifolia",
    TRUE ~ str_to_sentence(x_std)
  )
}

is_beech <- function(species_vec) {
  sp <- species_vec %>% as.character() %>% str_to_lower()
  str_detect(sp, "fagus") | str_detect(sp, "beech") | str_detect(sp, "h[eê]tre")
}

pick_first_col <- function(df, patterns) {
  nms <- names(df)
  idx <- purrr::detect_index(
    nms,
    ~ any(str_detect(.x, regex(paste(patterns, collapse = "|"), ignore_case = TRUE)))
  )
  
  if (idx == 0) {
    return(NA_character_)
  }
  
  nms[[idx]]
}

coerce_site_species <- function(df) {
  site_col <- pick_first_col(df, c("^site$", "site_", "station", "plot", "parcelle", "stand"))
  species_col <- pick_first_col(df, c("espece", "species", "taxon", "sp$", "nom"))
  
  if (!is.na(site_col)) {
    df <- df %>% mutate(site_std = .data[[site_col]] %>% as.character() %>% str_squish() %>% str_to_title())
  }
  
  if (!is.na(species_col)) {
    df <- df %>% mutate(species_std = normalize_species(.data[[species_col]]))
  }
  
  list(
    data = df,
    site_col = site_col,
    species_col = species_col
  )
}

sheet_stratum <- function(sheet_name) {
  nm <- sheet_name %>% str_to_lower() %>% str_squish()
  
  dplyr::case_when(
    str_detect(nm, "arbre|tree") ~ "arbre",
    str_detect(nm, "gaule|sapl") ~ "gaule",
    str_detect(nm, "regen") ~ "regeneration",
    str_detect(nm, "veget") ~ "vegetation",
    TRUE ~ nm
  )
}

infer_count_col <- function(df) {
  ccol <- pick_first_col(df, c("^n$", "count", "abond", "nombre", "effectif", "tige", "stem", "individu"))
  
  if (is.na(ccol)) {
    return(NULL)
  }
  
  parsed <- suppressWarnings(as.numeric(df[[ccol]]))
  
  if (all(is.na(parsed))) {
    return(NULL)
  }
  
  ccol
}

summarise_site_beech <- function(df, stratum_label) {
  if (!("site_std" %in% names(df))) {
    return(NULL)
  }
  
  count_col <- infer_count_col(df)
  
  if (is.null(count_col)) {
    base_df <- df %>% mutate(.count = 1)
  } else {
    base_df <- df %>% mutate(.count = coalesce(suppressWarnings(as.numeric(.data[[count_col]])), 0))
  }
  
  if (!("species_std" %in% names(base_df))) {
    return(
      base_df %>%
        group_by(site_std) %>%
        summarise(
          stratum = stratum_label,
          n_records = n(),
          total_abundance = sum(.count, na.rm = TRUE),
          .groups = "drop"
        )
    )
  }
  
  base_df %>%
    mutate(is_beech = is_beech(species_std)) %>%
    group_by(site_std) %>%
    summarise(
      stratum = stratum_label,
      n_records = n(),
      total_abundance = sum(.count, na.rm = TRUE),
      beech_abundance = sum(.count[is_beech], na.rm = TRUE),
      beech_prop = ifelse(total_abundance > 0, beech_abundance / total_abundance, NA_real_),
      .groups = "drop"
    )
}

calc_diversity <- function(df, stratum_label) {
  if (!("site_std" %in% names(df)) || !("species_std" %in% names(df))) {
    return(NULL)
  }
  
  count_col <- infer_count_col(df)
  
  work_df <- if (is.null(count_col)) {
    df %>% mutate(.count = 1)
  } else {
    df %>% mutate(.count = coalesce(suppressWarnings(as.numeric(.data[[count_col]])), 0))
  }
  
  agg <- work_df %>%
    filter(!is.na(species_std), species_std != "") %>%
    group_by(site_std, species_std) %>%
    summarise(n = sum(.count, na.rm = TRUE), .groups = "drop")
  
  if (nrow(agg) == 0) {
    return(NULL)
  }
  
  agg %>%
    group_by(site_std) %>%
    mutate(p = n / sum(n)) %>%
    summarise(
      stratum = stratum_label,
      species_richness = n_distinct(species_std),
      shannon = -sum(ifelse(p > 0, p * log(p), 0), na.rm = TRUE),
      evenness = ifelse(species_richness > 1, shannon / log(species_richness), NA_real_),
      .groups = "drop"
    )
}

calc_basal_area <- function(df) {
  if (!("site_std" %in% names(df))) {
    return(NULL)
  }
  
  diameter_col <- pick_first_col(df, c("dbh", "dhp", "diam", "d_hp", "diametre", "diameter"))
  
  if (is.na(diameter_col)) {
    return(NULL)
  }
  
  out <- df %>%
    mutate(
      diameter_raw = suppressWarnings(as.numeric(.data[[diameter_col]])),
      diameter_cm = ifelse(diameter_raw > 3, diameter_raw, diameter_raw * 10),
      basal_area_m2 = pi * (diameter_cm / 200)^2,
      is_beech = if ("species_std" %in% names(df)) is_beech(species_std) else FALSE
    ) %>%
    filter(!is.na(diameter_cm), diameter_cm > 0)
  
  if (nrow(out) == 0) {
    return(NULL)
  }
  
  out %>%
    group_by(site_std) %>%
    summarise(
      stems_for_basal_area = n(),
      total_basal_area_m2 = sum(basal_area_m2, na.rm = TRUE),
      beech_basal_area_m2 = sum(basal_area_m2[is_beech], na.rm = TRUE),
      beech_basal_area_prop = ifelse(total_basal_area_m2 > 0, beech_basal_area_m2 / total_basal_area_m2, NA_real_),
      diameter_source_col = diameter_col,
      diameter_assumption = "Values >3 treated as cm; <=3 treated as dm and converted to cm",
      .groups = "drop"
    )
}

write_csv_safe <- function(df, path) {
  if (is.null(df) || nrow(df) == 0) {
    return(invisible(FALSE))
  }
  readr::write_csv(df, path)
  invisible(TRUE)
}