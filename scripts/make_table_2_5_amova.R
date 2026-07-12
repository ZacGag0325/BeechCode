# scripts/make_table_2_5_amova.R
############################################################
# Thesis-ready Table 2.5: AMOVA results for American beech
# microsatellite dataset.
#
# This script formats existing AMOVA output files only. It does
# not rerun AMOVA or permutation tests.
#
# Inputs:
# - outputs/tables/amova_results.csv
# - outputs/tables/supplementary/amova_randtest_summary.csv
#   fallback: outputs/tables/amova_randtest_summary.csv
#
# Outputs:
# - outputs/tables/table_2_5_amova.csv
# - outputs/tables/table_2_5_amova.docx
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(flextable)
  library(officer)
  library(readr)
  library(stringr)
  library(tibble)
})

# ----------------------------
# Project-root and path helpers
# ----------------------------
find_project_root <- function() {
  candidates <- c(getwd(), normalizePath(file.path(getwd(), ".."), mustWork = FALSE))
  cmd_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
  if (length(cmd_file) > 0 && nzchar(cmd_file[1])) {
    candidates <- c(candidates, dirname(normalizePath(cmd_file[1], mustWork = FALSE)))
  }
  
  for (start in unique(candidates)) {
    cur <- normalizePath(start, mustWork = FALSE)
    repeat {
      if (file.exists(file.path(cur, "scripts", "make_table_2_5_amova.R"))) return(cur)
      parent <- dirname(cur)
      if (identical(parent, cur)) break
      cur <- parent
    }
  }
  
  stop(
    "[make_table_2_5_amova] Cannot find project root containing scripts/make_table_2_5_amova.R.",
    call. = FALSE
  )
}

PROJECT_ROOT <- find_project_root()
setwd(PROJECT_ROOT)

TABLES_DIR <- file.path(PROJECT_ROOT, "outputs", "tables")
AMOVA_RESULTS_FILE <- file.path(TABLES_DIR, "amova_results.csv")
RANDTEST_CANDIDATES <- c(
  file.path(TABLES_DIR, "supplementary", "amova_randtest_summary.csv"),
  file.path(TABLES_DIR, "amova_randtest_summary.csv")
)
OUTPUT_CSV <- file.path(TABLES_DIR, "table_2_5_amova.csv")
OUTPUT_DOCX <- file.path(TABLES_DIR, "table_2_5_amova.docx")

message_table <- function(...) {
  message("[make_table_2_5_amova] ", ...)
}

require_file <- function(path, description) {
  if (!file.exists(path)) {
    stop(
      "[make_table_2_5_amova] Missing ", description, ": ", path,
      call. = FALSE
    )
  }
  path
}

first_existing_file <- function(paths, description) {
  path <- paths[file.exists(paths)][1]
  if (is.na(path)) {
    stop(
      "[make_table_2_5_amova] Missing ", description, ". Tried:\n- ",
      paste(paths, collapse = "\n- "),
      call. = FALSE
    )
  }
  path
}

require_columns <- function(df, required, df_name) {
  missing <- setdiff(required, names(df))
  if (length(missing) > 0) {
    stop(
      "[make_table_2_5_amova] ", df_name, " is missing required column(s): ",
      paste(missing, collapse = ", "),
      ". Available columns are: ", paste(names(df), collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

find_column <- function(df, candidates, df_name, purpose) {
  hit <- candidates[candidates %in% names(df)][1]
  if (is.na(hit)) {
    stop(
      "[make_table_2_5_amova] Could not find ", purpose, " column in ", df_name,
      ". Tried: ", paste(candidates, collapse = ", "),
      ". Available columns are: ", paste(names(df), collapse = ", "),
      call. = FALSE
    )
  }
  hit
}

format_p_value <- function(x) {
  ifelse(is.na(x), NA_character_, sprintf("%.3f", round(x, 3)))
}

# ----------------------------
# Read and validate inputs
# ----------------------------
amova_results_path <- require_file(AMOVA_RESULTS_FILE, "AMOVA results CSV")
randtest_path <- first_existing_file(RANDTEST_CANDIDATES, "AMOVA randtest summary CSV")

message_table("Reading AMOVA results: ", amova_results_path)
amova_results <- readr::read_csv(amova_results_path, show_col_types = FALSE)

message_table("Reading AMOVA randtest summary: ", randtest_path)
randtest_summary <- readr::read_csv(randtest_path, show_col_types = FALSE)

require_columns(amova_results, c("Model", "Source"), "amova_results.csv")
require_columns(randtest_summary, c("model", "component", "p_value"), "amova_randtest_summary.csv")

variance_col <- find_column(
  amova_results,
  candidates = c("Sigma", "Variance component", "Variance_component", "variance_component", "Variance"),
  df_name = "amova_results.csv",
  purpose = "variance component"
)
percent_col <- find_column(
  amova_results,
  candidates = c("%", "Variation (%)", "Variation_percent", "variation_percent", "Percent", "Percentage", "percent"),
  df_name = "amova_results.csv",
  purpose = "variation percentage"
)

# ----------------------------
# Clean AMOVA rows and attach permutation p-values
# ----------------------------
p_value_lookup <- tribble(
  ~Model, ~`Source of variation`, ~component,
  "Site_only", "Among sites", "component_3",
  "Site_only", "Among individuals within sites", "component_2",
  "Site_only", "Within individuals", "component_1",
  "NorthSouth_Site_hierarchical", "Between regions", "component_4",
  "NorthSouth_Site_hierarchical", "Among sites within regions", "component_3",
  "NorthSouth_Site_hierarchical", "Among individuals within sites", "component_2",
  "NorthSouth_Site_hierarchical", "Within individuals", "component_1"
)

model_labels <- c(
  Site_only = "Site-level AMOVA",
  NorthSouth_Site_hierarchical = "North–south hierarchical AMOVA"
)

cleaned_amova <- amova_results %>%
  filter(
    !str_detect(Source, fixed("Total variations", ignore_case = TRUE)),
    !str_detect(Source, fixed("Phi", ignore_case = TRUE)),
    Model %in% names(model_labels)
  ) %>%
  mutate(
    `Source of variation` = case_when(
      Model == "Site_only" & Source == "Variations  Between pop" ~ "Among sites",
      Model == "Site_only" & Source == "Variations  Between samples Within pop" ~ "Among individuals within sites",
      Model == "Site_only" & Source == "Variations  Within samples" ~ "Within individuals",
      Model == "NorthSouth_Site_hierarchical" & Source == "Variations  Between Region" ~ "Between regions",
      Model == "NorthSouth_Site_hierarchical" & Source == "Variations  Between Site Within Region" ~ "Among sites within regions",
      Model == "NorthSouth_Site_hierarchical" & Source == "Variations  Between samples Within Site" ~ "Among individuals within sites",
      Model == "NorthSouth_Site_hierarchical" & Source == "Variations  Within samples" ~ "Within individuals",
      TRUE ~ NA_character_
    )
  )

unmapped_sources <- cleaned_amova %>%
  filter(is.na(`Source of variation`)) %>%
  distinct(Model, Source)
if (nrow(unmapped_sources) > 0) {
  stop(
    "[make_table_2_5_amova] Found AMOVA variance row(s) with unmapped Source labels:\n",
    paste(sprintf("- Model: %s | Source: %s", unmapped_sources$Model, unmapped_sources$Source), collapse = "\n"),
    call. = FALSE
  )
}

p_values <- randtest_summary %>%
  transmute(
    Model = model,
    component = component,
    p_value = as.numeric(p_value)
  )

final_table <- cleaned_amova %>%
  left_join(p_value_lookup, by = c("Model", "Source of variation")) %>%
  left_join(p_values, by = c("Model", "component")) %>%
  mutate(
    Model = recode(Model, !!!model_labels),
    `Variance component` = round(as.numeric(.data[[variance_col]]), 3),
    `Variation (%)` = round(as.numeric(.data[[percent_col]]), 2),
    `p-value` = format_p_value(p_value)
  ) %>%
  select(Model, `Source of variation`, `Variance component`, `Variation (%)`, `p-value`)

if (any(is.na(final_table$`p-value`))) {
  missing_p <- final_table %>% filter(is.na(`p-value`))
  stop(
    "[make_table_2_5_amova] Missing p-value(s) after joining randtest summary. Affected row(s):\n",
    paste(sprintf("- %s | %s", missing_p$Model, missing_p$`Source of variation`), collapse = "\n"),
    call. = FALSE
  )
}

expected_n <- nrow(p_value_lookup)
if (nrow(final_table) != expected_n) {
  stop(
    "[make_table_2_5_amova] Expected ", expected_n, " AMOVA table row(s), but found ",
    nrow(final_table), ". Check input AMOVA Source labels and models.",
    call. = FALSE
  )
}

# ----------------------------
# Write CSV and Word table
# ----------------------------
dir.create(TABLES_DIR, recursive = TRUE, showWarnings = FALSE)

readr::write_csv(final_table, OUTPUT_CSV, na = "")

note_text <- paste(
  "Note. AMOVA was performed on the clone-corrected dataset (N = 268).",
  "The site-level model tested differentiation among the 12 study sites.",
  "The hierarchical model tested differentiation between southern and northern regions, among sites within regions, among individuals within sites, and within individuals.",
  "Significance was assessed using 999 permutations.",
  "Negative variance components can occur when differentiation among groups is effectively zero and were interpreted as no regional genetic structure."
)

title_text <- "Table 2.5. Analysis of molecular variance (AMOVA) of clone-corrected American beech (Fagus grandifolia Ehrh.) genotypes."

ft <- flextable(final_table) %>%
  theme_booktabs() %>%
  fontsize(size = 10, part = "all") %>%
  font(fontname = "Times New Roman", part = "all") %>%
  align(align = "left", part = "all") %>%
  align(j = c("Variance component", "Variation (%)", "p-value"), align = "right", part = "body") %>%
  bold(part = "header") %>%
  autofit()

doc <- read_docx() %>%
  body_add_par(title_text, style = "Normal") %>%
  body_add_flextable(ft) %>%
  body_add_par(note_text, style = "Normal")

print(doc, target = OUTPUT_DOCX)

# ----------------------------
# Console output
# ----------------------------
print(final_table)
message_table("Saved CSV: ", OUTPUT_CSV)
message_table("Saved Word table: ", OUTPUT_DOCX)