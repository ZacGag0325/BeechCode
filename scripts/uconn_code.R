#!/usr/bin/env Rscript

# Prepare a clean, shareable metadata workbook for potential UConn/American beech
# landscape genomics collaboration.
#
# This script intentionally excludes thesis-sensitive genetic/clonality/root-
# connection results and exports only basic site/sample metadata for the sites
# present in the 24-site workbook but absent from the MSc/thesis workbook.

required_packages <- c("readxl", "dplyr", "tidyr", "stringr", "stringi", "purrr", "writexl", "janitor", "tibble", "readr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  message("Installing missing R packages: ", paste(missing_packages, collapse = ", "))
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}
invisible(lapply(required_packages, library, character.only = TRUE))

find_script_path <- function() {
  file_arg <- commandArgs(trailingOnly = FALSE)
  file_arg <- file_arg[stringr::str_detect(file_arg, "^--file=")]
  if (length(file_arg) > 0) return(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE))
  
  # Works when sourced with source("scripts/prepare_uconn_leaf_metadata.R").
  sys_files <- vapply(sys.frames(), function(frame) {
    if (!is.null(frame$ofile)) frame$ofile else NA_character_
  }, character(1))
  sys_files <- stats::na.omit(sys_files)
  if (length(sys_files) > 0) return(normalizePath(sys_files[[length(sys_files)]], mustWork = FALSE))
  
  NA_character_
}

required_raw_files <- c(
  "donnees_F.grandifolia_2024.xlsx",
  "donnees_modifiees_west_summer2024 copie.xlsx"
)

has_required_raw_files <- function(path) {
  all(file.exists(file.path(path, "data", "raw", required_raw_files)))
}

has_data_raw_folder <- function(path) {
  dir.exists(file.path(path, "data", "raw"))
}

candidate_parent_first <- function(path) {
  path <- normalizePath(path, mustWork = FALSE)
  # If R/RStudio is currently inside a folder named BeechCode.Rproj, the actual
  # BeechCode project/data folder may be its parent. Prefer that parent first.
  if (stringr::str_detect(basename(path), stringr::regex("\\.Rproj$", ignore_case = TRUE))) {
    return(c(dirname(path), path))
  }
  c(path, dirname(path))
}

find_project_root <- function() {
  script_path <- find_script_path()
  seed_paths <- unique(stats::na.omit(c(
    candidate_parent_first(getwd()),
    if (!is.na(script_path)) candidate_parent_first(dirname(script_path)) else NA_character_,
    if (!is.na(script_path)) dirname(dirname(script_path)) else NA_character_
  )))
  
  searched <- character()
  roots_with_git <- character()
  
  # First pass: prefer the directory that actually contains the required raw files.
  for (start in seed_paths) {
    current <- normalizePath(start, mustWork = FALSE)
    repeat {
      searched <- unique(c(searched, current))
      if (has_required_raw_files(current)) return(current)
      if (dir.exists(file.path(current, ".git"))) roots_with_git <- unique(c(roots_with_git, current))
      parent <- dirname(current)
      if (identical(parent, current)) break
      current <- parent
    }
  }
  
  # Second pass: if files have not been downloaded yet, prefer a parent that has data/raw.
  data_roots <- searched[vapply(searched, has_data_raw_folder, logical(1))]
  if (length(data_roots) > 0) return(data_roots[[1]])
  
  # Third pass: fall back to a git root, but avoid choosing a nested *.Rproj folder
  # when the actual project root is likely its parent.
  if (length(roots_with_git) > 0) {
    non_rproj_roots <- roots_with_git[!stringr::str_detect(basename(roots_with_git), stringr::regex("\\.Rproj$", ignore_case = TRUE))]
    if (length(non_rproj_roots) > 0) return(non_rproj_roots[[1]])
    return(dirname(roots_with_git[[1]]))
  }
  
  normalizePath(getwd(), mustWork = FALSE)
}

project_root <- find_project_root()
message("Using project root: ", project_root)

path_from_root <- function(...) file.path(project_root, ...)

raw_24_path <- path_from_root("data", "raw", required_raw_files[[1]])
thesis_path <- path_from_root("data", "raw", required_raw_files[[2]])
map_path <- path_from_root("data", "raw", "map_csv_east_sites_qgis_points.csv")
out_dir <- path_from_root("outputs", "collaboration_uconn")
out_xlsx <- file.path(out_dir, "uconn_additional_12_sites_metadata.xlsx")
out_site_csv <- file.path(out_dir, "site_summary_additional_12_sites.csv")
out_sample_csv <- file.path(out_dir, "sample_metadata_additional_12_sites.csv")

stop_if_missing <- function(path, label) {
  if (!file.exists(path)) {
    stop(
      label, " is missing: ", path, "\n",
      "Current working directory: ", getwd(), "\n",
      "Detected project root: ", project_root, "\n",
      "Place the required workbook under data/raw/ in the project root, run the script from the BeechCode root, ",
      "or source the script from its repository path so the project root can be detected.",
      call. = FALSE
    )
  }
}
stop_if_missing(raw_24_path, "24-site field workbook")
stop_if_missing(thesis_path, "MSc/thesis workbook")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

normalize_name <- function(x) {
  x |>
    stringi::stri_trans_general("Latin-ASCII") |>
    stringr::str_to_lower() |>
    stringr::str_replace_all("[^a-z0-9]+", "_") |>
    stringr::str_replace_all("^_|_$", "")
}

normalize_site_value <- function(x) {
  x <- as.character(x)
  x <- stringr::str_squish(x)
  x <- stringi::stri_trans_general(x, "Latin-ASCII")
  x <- stringr::str_to_upper(x)
  x[x == "" | is.na(x)] <- NA_character_
  x
}

safe_read_sheet <- function(path, sheet) {
  readxl::read_excel(path, sheet = sheet, .name_repair = "unique") |>
    janitor::remove_empty(c("rows", "cols"))
}

sheet_exists <- function(path, sheet) sheet %in% readxl::excel_sheets(path)

find_col <- function(dat, patterns, required = FALSE, label = "column") {
  if (is.null(dat) || ncol(dat) == 0) {
    if (required) stop("Cannot detect ", label, " because the data frame has no columns.", call. = FALSE)
    return(NA_character_)
  }
  original <- names(dat)
  normalized <- normalize_name(original)
  for (pat in patterns) {
    hit <- which(stringr::str_detect(normalized, pat))
    if (length(hit) > 0) return(original[hit[1]])
  }
  if (required) {
    stop("Could not detect ", label, ". Available columns: ", paste(original, collapse = ", "), call. = FALSE)
  }
  NA_character_
}

find_site_col <- function(dat, required = FALSE) {
  find_col(
    dat,
    c("^site$", "^sites$", "^code_site$", "^site_code$", "^codesite$", "^nom_site$", "^site_name$", "^station$", "^localite$", "^location$"),
    required = required,
    label = "site column"
  )
}

id_col_patterns <- c("id", "ident", "tree", "arbre", "sample", "echant", "lab", "code", "tag", "plant", "individu")

sensitive_patterns <- c(
  "genotype", "microsat", "allele", "locus", "^ssr", "bruvo", "clone", "clonal", "clon", "mlg", "mll",
  "structure", "pca", "amova", "qvalue", "q_value", "admixture", "kinship", "genetic", "genetique",
  "excav", "root", "racine", "connect", "connexion", "connected", "unconnected", "apparentement"
)

safe_col_patterns <- c(
  "^site$", "code_site", "site_code", "nom_site", "site_name", "station", "localite", "location",
  "sample", "echant", "lab", "tree", "arbre", "id", "tag", "code", "latitude", "longitude", "lat", "lon", "lng",
  "elevation", "altitude", "region", "cluster", "secteur", "zone", "note", "comment", "remarque", "date", "year", "annee",
  "stade", "stage", "dbh", "dhp", "diam", "diametre", "sante", "health", "maladie", "disease", "symptom", "canopy", "hauteur", "height"
)

is_sensitive_col <- function(cols) {
  norm <- normalize_name(cols)
  purrr::map_lgl(norm, ~ any(stringr::str_detect(.x, sensitive_patterns)))
}

select_safe_columns <- function(dat, force_cols = character()) {
  if (is.null(dat) || ncol(dat) == 0) return(dat)
  cols <- names(dat)
  norm <- normalize_name(cols)
  keep_safe <- purrr::map_lgl(norm, ~ any(stringr::str_detect(.x, safe_col_patterns)))
  keep <- (cols %in% force_cols) | keep_safe
  keep <- keep & !is_sensitive_col(cols)
  # If a forced site/sample column has a broad name like genetique-safe ID, retain it.
  keep[cols %in% force_cols] <- TRUE
  dat[, keep, drop = FALSE]
}

extract_sites_from_sheet <- function(path, sheet) {
  dat <- safe_read_sheet(path, sheet)
  site_col <- find_site_col(dat)
  if (is.na(site_col)) return(character())
  normalize_site_value(dat[[site_col]]) |> unique() |> stats::na.omit() |> as.character()
}

extract_thesis_sites <- function() {
  sheets <- readxl::excel_sheets(thesis_path)
  candidates <- intersect(c("genetique", "site_lookup"), sheets)
  sites <- purrr::map(candidates, ~ extract_sites_from_sheet(thesis_path, .x)) |> unlist(use.names = FALSE)
  sites <- unique(normalize_site_value(sites)) |> stats::na.omit() |> as.character()
  if (length(sites) == 0) {
    stop("Could not detect thesis sites. Expected a site column in the thesis workbook genetique and/or site_lookup sheet.", call. = FALSE)
  }
  sites
}

field_sheets <- c("arbre", "gaule", "regeneration", "vegetation", "genetique", "sol")
available_field_sheets <- intersect(field_sheets, readxl::excel_sheets(raw_24_path))
if (length(available_field_sheets) == 0) {
  stop("None of the expected field sheets were found in the 24-site workbook: ", paste(field_sheets, collapse = ", "), call. = FALSE)
}

thesis_sites <- extract_thesis_sites()

# For a leaf/genomics collaboration, prefer the 24-site genetique sheet as the
# authoritative site list when it is available. Some non-genetic field sheets can
# contain extra site-like codes or inconsistent capitalization that should not
# expand the leaf-sample export set.
if ("genetique" %in% available_field_sheets) {
  all_sites <- extract_sites_from_sheet(raw_24_path, "genetique") |>
    normalize_site_value() |>
    unique() |>
    stats::na.omit() |>
    as.character()
  site_source_note <- "24-site genetique sheet"
} else {
  all_sites <- purrr::map(available_field_sheets, ~ extract_sites_from_sheet(raw_24_path, .x)) |>
    unlist(use.names = FALSE) |>
    normalize_site_value() |>
    unique() |>
    stats::na.omit() |>
    as.character()
  site_source_note <- paste("all available field sheets:", paste(available_field_sheets, collapse = ", "))
}

if (length(all_sites) == 0) {
  stop("Could not detect any sites in the 24-site workbook. Expected a recognizable site column in at least one of: ", paste(available_field_sheets, collapse = ", "), call. = FALSE)
}

extra_sites <- setdiff(all_sites, thesis_sites) |> sort()
if (length(extra_sites) == 0) {
  stop("Automatic detection found no additional/non-thesis sites. Check that site columns use comparable codes in both workbooks.", call. = FALSE)
}
if (length(extra_sites) != 12) {
  stop(
    "Detected ", length(extra_sites), " additional/non-thesis sites, not exactly 12.\n",
    "Site source used: ", site_source_note, "\n",
    "Thesis sites detected (", length(thesis_sites), "): ", paste(sort(thesis_sites), collapse = ", "), "\n",
    "Total sites detected in 24-site workbook (", length(all_sites), "): ", paste(sort(all_sites), collapse = ", "), "\n",
    "Additional/non-thesis sites detected (", length(extra_sites), "): ", paste(extra_sites, collapse = ", "), "\n",
    "This usually means a site code is misspelled/inconsistently coded or a non-leaf site code was included. Fix the site codes or review the detected site column before exporting.",
    call. = FALSE
  )
}

message("Site source used: ", site_source_note)
message("Thesis sites detected (", length(thesis_sites), "): ", paste(sort(thesis_sites), collapse = ", "))
message("Total sites in 24-site workbook (", length(all_sites), "): ", paste(sort(all_sites), collapse = ", "))
message("Additional/non-thesis sites (", length(extra_sites), "): ", paste(extra_sites, collapse = ", "))

filter_extra_sites <- function(dat, site_col) {
  dat |>
    dplyr::mutate(.site_filter_value = normalize_site_value(.data[[site_col]])) |>
    dplyr::filter(.site_filter_value %in% extra_sites) |>
    dplyr::select(-.site_filter_value)
}

read_field_extra <- function(sheet) {
  dat <- safe_read_sheet(raw_24_path, sheet)
  site_col <- find_site_col(dat, required = TRUE)
  list(sheet = sheet, data = filter_extra_sites(dat, site_col), site_col = site_col)
}
field_data <- purrr::map(available_field_sheets, read_field_extra)
names(field_data) <- available_field_sheets

count_records <- function(sheet_name, label) {
  if (!sheet_name %in% names(field_data)) return(tibble::tibble(site = extra_sites, !!label := 0L))
  item <- field_data[[sheet_name]]
  item$data |>
    dplyr::transmute(site = normalize_site_value(.data[[item$site_col]])) |>
    dplyr::count(site, name = label) |>
    dplyr::right_join(tibble::tibble(site = extra_sites), by = "site") |>
    dplyr::mutate(dplyr::across(dplyr::all_of(label), ~ tidyr::replace_na(.x, 0L)))
}

site_summary <- tibble::tibble(site = extra_sites)
for (pair in list(c("genetique", "leaf_genetic_sample_count"), c("arbre", "tree_record_count"), c("gaule", "sapling_record_count"), c("regeneration", "regeneration_record_count"), c("vegetation", "vegetation_record_count"), c("sol", "soil_record_count"))) {
  site_summary <- dplyr::left_join(site_summary, count_records(pair[1], pair[2]), by = "site")
}

# Add safe site-level fields discovered in field sheets, taking the first non-missing value per site.
site_level_patterns <- c("latitude", "longitude", "^lat$", "^lon$", "^lng$", "elevation", "altitude", "region", "cluster", "secteur", "zone", "note", "comment", "remarque")
for (item in field_data) {
  safe <- select_safe_columns(item$data, force_cols = item$site_col)
  safe_cols <- setdiff(names(safe), item$site_col)
  site_cols <- safe_cols[purrr::map_lgl(normalize_name(safe_cols), ~ any(stringr::str_detect(.x, site_level_patterns)))]
  if (length(site_cols) == 0) next
  tmp <- safe |>
    dplyr::mutate(site = normalize_site_value(.data[[item$site_col]])) |>
    dplyr::select(site, dplyr::all_of(site_cols)) |>
    dplyr::group_by(site) |>
    dplyr::summarise(dplyr::across(dplyr::everything(), ~ dplyr::first(stats::na.omit(as.character(.x)), default = NA_character_)), .groups = "drop")
  names(tmp)[-1] <- paste0(item$sheet, "_", names(tmp)[-1])
  site_summary <- dplyr::left_join(site_summary, tmp, by = "site")
}

site_summary <- site_summary |>
  dplyr::mutate(ecological_reserve_inferred = stringr::str_detect(stringr::str_to_lower(site), "reserve|réserve|ecological|ecologique|écologique|r[eé]serve"))

if ("genetique" %in% names(field_data)) {
  genetic_item <- field_data[["genetique"]]
  sample_metadata <- select_safe_columns(genetic_item$data, force_cols = genetic_item$site_col)
} else {
  sample_metadata <- tibble::tibble(site = character())
}

field_metadata_available <- purrr::map_dfr(field_data, function(item) {
  dat <- item$data
  id_cols <- names(dat)[purrr::map_lgl(normalize_name(names(dat)), ~ any(stringr::str_detect(.x, id_col_patterns)))]
  tibble::tibble(
    sheet_name = item$sheet,
    rows_retained = nrow(dat),
    detected_site_column = item$site_col,
    detected_id_columns = paste(id_cols, collapse = ", "),
    notes = "Filtered to additional/non-thesis sites; sensitive result columns excluded from exports."
  )
})

readme <- tibble::tibble(
  field = c("date_generated", "purpose", "site_scope", "language_note", "data_use_note", "review_note"),
  value = c(
    as.character(Sys.Date()),
    "Clean preliminary metadata workbook for potential UConn/American beech landscape genomics collaboration.",
    "Includes only additional sites detected in the 24-site workbook that are not included in the MSc/thesis site list.",
    "Metadata were originally recorded in French and may require review/translation before external sharing.",
    "No genetic results, microsatellite genotypes, MLG/MLL, Bruvo distances, clone calls, STRUCTURE/PCA/AMOVA outputs, excavation data, or root-connection data are included.",
    "Review all sheets before sharing externally."
  )
)

if (sheet_exists(raw_24_path, "dictionnaire")) {
  dictionary_original <- safe_read_sheet(raw_24_path, "dictionnaire")
} else {
  dictionary_original <- tibble::tibble(note = "No dictionnaire sheet found in 24-site workbook.")
}

dict_names <- names(dictionary_original)
field_col <- find_col(dictionary_original, c("field", "champ", "variable", "colonne", "nom"))
desc_col <- find_col(dictionary_original, c("description", "definition", "desc", "details"))
if (!is.na(field_col)) {
  dictionary_translation_template <- tibble::tibble(
    original_field_name = as.character(dictionary_original[[field_col]]),
    original_description_fr = if (!is.na(desc_col)) as.character(dictionary_original[[desc_col]]) else NA_character_
  )
} else {
  dictionary_translation_template <- tibble::tibble(
    original_field_name = dict_names,
    original_description_fr = NA_character_
  )
}
suggestions <- c(site = "Site", latitude = "Latitude", longitude = "Longitude", dhp = "DBH", stade = "Life_stage", date = "Date")
dictionary_translation_template <- dictionary_translation_template |>
  dplyr::mutate(
    suggested_field_name_en = dplyr::case_when(
      normalize_name(original_field_name) %in% names(suggestions) ~ unname(suggestions[normalize_name(original_field_name)]),
      TRUE ~ NA_character_
    ),
    suggested_description_en = NA_character_,
    notes = NA_character_
  )

workbook <- list(
  README = readme,
  site_summary = site_summary,
  sample_metadata = sample_metadata,
  field_metadata_available = field_metadata_available,
  dictionary_original = dictionary_original,
  dictionary_translation_template = dictionary_translation_template
)

if (file.exists(map_path)) {
  map_dat <- read.csv(map_path, check.names = FALSE, stringsAsFactors = FALSE)
  map_site_col <- find_site_col(map_dat)
  if (!is.na(map_site_col)) {
    map_points <- filter_extra_sites(map_dat, map_site_col) |> select_safe_columns(force_cols = map_site_col)
    workbook$map_points <- map_points
  } else {
    message("Optional map CSV exists but no site column was detected; map_points sheet skipped.")
  }
}

writexl::write_xlsx(workbook, out_xlsx)
readr::write_csv(site_summary, out_site_csv)
readr::write_csv(sample_metadata, out_sample_csv)

message("Output Excel workbook: ", out_xlsx)
message("Output site summary CSV: ", out_site_csv)
message("Output sample metadata CSV: ", out_sample_csv)
message("Number of sites exported: ", length(extra_sites))
message("Number of samples exported: ", nrow(sample_metadata))
message("Reminder: review the output carefully before sharing externally.")