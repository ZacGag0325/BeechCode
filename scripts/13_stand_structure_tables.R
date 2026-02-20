# scripts/13_stand_structure_tables.R
# ------------------------------------------------------------
# Site-level stand structure summaries from tree/stem field data
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
})

message_info <- function(...) message("[INFO] ", ...)
message_warn <- function(...) message("[WARN] ", ...)

normalize_name <- function(x) {
  x %>%
    stringr::str_to_lower() %>%
    stringr::str_replace_all("[^a-z0-9]+", "_") %>%
    stringr::str_replace_all("_+", "_") %>%
    stringr::str_replace("^_", "") %>%
    stringr::str_replace("_$", "")
}

find_best_column <- function(df, candidates) {
  nms <- names(df)
  nms_norm <- normalize_name(nms)
  
  for (cand in candidates) {
    idx <- which(nms_norm == normalize_name(cand))
    if (length(idx) > 0) return(nms[idx[1]])
  }
  
  for (cand in candidates) {
    pat <- normalize_name(cand)
    idx <- which(stringr::str_detect(nms_norm, paste0("(^|_)", pat, "(_|$)|", pat)))
    if (length(idx) > 0) return(nms[idx[1]])
  }
  
  NULL
}

read_tabular <- function(path, sheet = NULL) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "csv") {
    readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
  } else if (ext %in% c("tsv", "tab")) {
    readr::read_tsv(path, show_col_types = FALSE, progress = FALSE)
  } else if (ext == "txt") {
    readr::read_delim(path, delim = "\t", show_col_types = FALSE, progress = FALSE)
  } else if (ext %in% c("xlsx", "xls")) {
    if (!requireNamespace("readxl", quietly = TRUE)) {
      stop("Excel file found but package 'readxl' is not installed.")
    }
    if (is.null(sheet)) {
      readxl::read_excel(path) %>% as_tibble()
    } else {
      readxl::read_excel(path, sheet = sheet) %>% as_tibble()
    }
  } else {
    stop("Unsupported extension: ", ext)
  }
}

excel_candidate_score <- function(df) {
  if (nrow(df) == 0 || ncol(df) == 0) return(-Inf)
  nms <- normalize_name(names(df))
  score <- 0
  score <- score + as.integer(any(nms %in% c("site", "population", "pop", "site_id", "stand", "plot"))) * 5
  score <- score + as.integer(any(stringr::str_detect(nms, "dbh|diameter|dhp|d130|cap"))) * 8
  score <- score + as.integer(any(stringr::str_detect(nms, "height|ht")))
  score <- score + as.integer(any(stringr::str_detect(nms, "species|espece|taxon")))
  score
}

read_best_excel_sheet <- function(path) {
  sheets <- readxl::excel_sheets(path)
  if (length(sheets) == 0) stop("Excel file has no sheets: ", path)
  
  best_score <- -Inf
  best_sheet <- sheets[1]
  best_df <- NULL
  
  for (sh in sheets) {
    candidate <- tryCatch(read_tabular(path, sheet = sh), error = function(e) NULL)
    if (is.null(candidate)) next
    sc <- excel_candidate_score(candidate)
    if (sc > best_score) {
      best_score <- sc
      best_sheet <- sh
      best_df <- candidate
    }
  }
  
  if (is.null(best_df)) stop("Could not read any Excel sheet from: ", path)
  message_info("Excel detected; selected sheet: ", best_sheet, " (score=", best_score, ")")
  best_df
}

choose_dataset <- function(paths) {
  preferred <- c("field", "stand", "structure", "tree", "stem", "inventory", "plot", "donnees", "genetique")
  score <- purrr::map_int(paths, ~ {
    nm <- normalize_name(basename(.x))
    sum(stringr::str_detect(nm, preferred))
  })
  paths[order(score, file.info(paths)$size, decreasing = TRUE)[1]]
}

compute_density <- function(df) {
  if (!("Area" %in% names(df))) return(tibble(Site = levels(df$Site), Stem_density = NA_real_))
  
  stems <- df %>% count(Site, name = "N_individuals")
  area <- df %>%
    group_by(Site) %>%
    summarize(area_value = suppressWarnings(first(Area[!is.na(Area) & Area > 0])), .groups = "drop")
  
  stems %>%
    left_join(area, by = "Site") %>%
    mutate(Stem_density = ifelse(!is.na(area_value) & area_value > 0, N_individuals / area_value, NA_real_)) %>%
    select(Site, Stem_density)
}

# 1) Detect dataset -----------------------------------------------------------
search_dirs <- c("inputs", "data")
candidates <- c()
for (d in search_dirs) {
  if (dir.exists(d)) {
    candidates <- c(candidates, list.files(
      d,
      pattern = "\\.(csv|tsv|tab|txt|xlsx|xls)$",
      full.names = TRUE,
      recursive = TRUE,
      ignore.case = TRUE
    ))
  }
}

if (length(candidates) == 0) {
  stop("No dataset detected in inputs/ or data/.")
}

selected_file <- choose_dataset(candidates)
message_info("Detected ", length(candidates), " candidate file(s).")
message_info("Selected dataset: ", selected_file)

selected_ext <- tolower(tools::file_ext(selected_file))
if (selected_ext %in% c("xlsx", "xls")) {
  raw_df <- read_best_excel_sheet(selected_file)
} else {
  raw_df <- read_tabular(selected_file)
}
if (nrow(raw_df) == 0) stop("Selected dataset is empty: ", selected_file)

# 2) Detect columns + clean ---------------------------------------------------
site_col <- find_best_column(raw_df, c("site", "population", "pop", "site_id", "population_id", "stand", "plot"))
dbh_col <- find_best_column(raw_df, c("dbh", "dbh_cm", "diameter", "diameter_cm", "d130", "cap", "dhp", "dhp_tige"))
height_col <- find_best_column(raw_df, c("height", "ht", "height_m", "tree_height"))
status_col <- find_best_column(raw_df, c("status", "alive_dead", "vital_status", "state"))
species_col <- find_best_column(raw_df, c("species", "taxon", "sp", "species_name", "espece"))
regen_col <- find_best_column(raw_df, c("regeneration_class", "regen", "size_class", "stage", "life_stage", "class"))
area_col <- find_best_column(raw_df, c("plot_area", "area", "area_ha", "sample_area", "site_area"))

if (is.null(site_col)) stop("No Site-like column detected.")
if (is.null(dbh_col)) stop("No DBH-like column detected.")

message_info("Detected columns:")
message_info("  Site: ", site_col)
message_info("  DBH: ", dbh_col)
message_info("  Height: ", ifelse(is.null(height_col), "<missing>", height_col))
message_info("  Status: ", ifelse(is.null(status_col), "<missing>", status_col))
message_info("  Species: ", ifelse(is.null(species_col), "<missing>", species_col))
message_info("  Regeneration class: ", ifelse(is.null(regen_col), "<missing>", regen_col))
message_info("  Area: ", ifelse(is.null(area_col), "<missing>", area_col))

df <- raw_df %>%
  mutate(
    Site = as.factor(.data[[site_col]]),
    DBH = suppressWarnings(readr::parse_number(as.character(.data[[dbh_col]]))),
    DBH_invalid = is.na(DBH) | DBH < 0
  )

if (!is.null(height_col)) df <- df %>% mutate(Height = suppressWarnings(readr::parse_number(as.character(.data[[height_col]]))))
if (!is.null(status_col)) df <- df %>% mutate(Status = as.character(.data[[status_col]]))
if (!is.null(species_col)) df <- df %>% mutate(Species = as.character(.data[[species_col]]))
if (!is.null(regen_col)) df <- df %>% mutate(Regeneration_Class = as.character(.data[[regen_col]]))
if (!is.null(area_col)) df <- df %>% mutate(Area = suppressWarnings(readr::parse_number(as.character(.data[[area_col]]))))

invalid_n <- sum(df$DBH_invalid, na.rm = TRUE)
if (invalid_n > 0) {
  message_warn(invalid_n, " row(s) have invalid/missing DBH and will be excluded from DBH summaries.")
}

df_valid_dbh <- df %>% filter(!DBH_invalid)

# 3) Build summary tables -----------------------------------------------------
out_dir <- file.path("outputs", "v1")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

basic_summary <- df %>%
  group_by(Site) %>%
  summarize(
    N_individuals = n(),
    Mean_DBH = mean(DBH[!DBH_invalid], na.rm = TRUE),
    SD_DBH = sd(DBH[!DBH_invalid], na.rm = TRUE),
    Median_DBH = median(DBH[!DBH_invalid], na.rm = TRUE),
    Min_DBH = suppressWarnings(min(DBH[!DBH_invalid], na.rm = TRUE)),
    Max_DBH = suppressWarnings(max(DBH[!DBH_invalid], na.rm = TRUE)),
    pct_missing_DBH = 100 * mean(DBH_invalid, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Min_DBH = ifelse(is.infinite(Min_DBH), NA_real_, Min_DBH),
    Max_DBH = ifelse(is.infinite(Max_DBH), NA_real_, Max_DBH)
  ) %>%
  left_join(compute_density(df), by = "Site")

readr::write_csv(basic_summary, file.path(out_dir, "stand_structure_summary_by_site.csv"))

dbh_labels <- c("0-5", "5-10", "10-20", "20-40", "40+")
dbh_dist <- df_valid_dbh %>%
  mutate(DBH_class = cut(DBH, breaks = c(0, 5, 10, 20, 40, Inf), labels = dbh_labels, right = FALSE, include.lowest = TRUE)) %>%
  group_by(Site, DBH_class) %>%
  summarize(n_individuals = n(), .groups = "drop") %>%
  complete(Site, DBH_class = factor(dbh_labels, levels = dbh_labels), fill = list(n_individuals = 0)) %>%
  group_by(Site) %>%
  mutate(prop_percent = ifelse(sum(n_individuals) > 0, 100 * n_individuals / sum(n_individuals), NA_real_)) %>%
  ungroup()

readr::write_csv(dbh_dist, file.path(out_dir, "dbh_class_distribution_by_site.csv"))

if (!is.null(regen_col)) {
  regen_tbl <- df %>%
    mutate(
      Regen_group = normalize_name(Regeneration_Class),
      Regen_group = case_when(
        stringr::str_detect(Regen_group, "seed") ~ "seedling",
        stringr::str_detect(Regen_group, "sap") ~ "sapling",
        stringr::str_detect(Regen_group, "adult|mature|tree") ~ "adult",
        TRUE ~ Regen_group
      )
    ) %>%
    group_by(Site, Regen_group) %>%
    summarize(n_individuals = n(), .groups = "drop") %>%
    group_by(Site) %>%
    mutate(prop_percent = 100 * n_individuals / sum(n_individuals)) %>%
    ungroup()
} else {
  message_warn("No regeneration class detected. Using DBH < 10 cm as regeneration.")
  regen_tbl <- df_valid_dbh %>%
    mutate(Regen_group = ifelse(DBH < 10, "regeneration", "adult")) %>%
    group_by(Site, Regen_group) %>%
    summarize(n_individuals = n(), .groups = "drop") %>%
    group_by(Site) %>%
    mutate(prop_percent = 100 * n_individuals / sum(n_individuals)) %>%
    ungroup()
}

readr::write_csv(regen_tbl, file.path(out_dir, "regeneration_summary_by_site.csv"))

if (!is.null(species_col)) {
  species_tbl <- df %>%
    mutate(
      Species_clean = stringr::str_squish(as.character(Species)),
      Fagus_group = ifelse(
        stringr::str_detect(stringr::str_to_lower(Species_clean), "fagus\\s+grandifolia|f\\.?\\s*grandifolia"),
        "Fagus grandifolia",
        "Other species"
      )
    ) %>%
    group_by(Site, Species_clean, Fagus_group) %>%
    summarize(n_individuals = n(), .groups = "drop") %>%
    group_by(Site) %>%
    mutate(relative_abundance_percent = 100 * n_individuals / sum(n_individuals)) %>%
    ungroup()
  
  readr::write_csv(species_tbl, file.path(out_dir, "species_composition_by_site.csv"))
} else {
  message_warn("Species column missing; skipping species composition output.")
}

regen_wide <- regen_tbl %>%
  select(Site, Regen_group, prop_percent) %>%
  mutate(metric = paste0("regen_prop_", normalize_name(Regen_group), "_pct")) %>%
  select(-Regen_group) %>%
  pivot_wider(names_from = metric, values_from = prop_percent, values_fill = 0)

pub_ready <- basic_summary %>%
  mutate(mean_dbh_sd = ifelse(is.na(Mean_DBH), NA_character_, paste0(round(Mean_DBH, 2), " ± ", round(SD_DBH, 2)))) %>%
  left_join(regen_wide, by = "Site") %>%
  select(Site, N_individuals, mean_dbh_sd, Median_DBH, Min_DBH, Max_DBH, pct_missing_DBH, Stem_density, starts_with("regen_prop_"))

readr::write_csv(pub_ready, file.path(out_dir, "publication_ready_site_table.csv"))

# 4) Short suggestions linking to genetic results ----------------------------
message_info("Suggested links with genetics:")
message_info(" - Join site stand metrics with clonality (MLG/MLL richness) by Site.")
message_info(" - Model He / allelic richness as a function of regeneration structure.")
message_info(" - Compare Fagus dominance with pairwise Fst or within-site kinship patterns.")

message_info("Done. Outputs written in: ", out_dir)
