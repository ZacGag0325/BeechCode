# Select one mature American beech sample per site for Michelle's landscape project.
# Rules: DBH > 10 cm, at least 8 km between selected trees, one tree per site.

required_packages <- c("readxl", "dplyr", "writexl")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Install the missing package(s) before running this script: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages(library(dplyr))

# source() does not change the working directory. Find the folder containing
# this script so the workbook can be stored beside the script if desired.
get_sourced_script_path <- function() {
  frame_paths <- vapply(
    sys.frames(),
    function(frame) {
      if (is.null(frame$ofile)) NA_character_ else as.character(frame$ofile)[1]
    },
    character(1)
  )
  frame_paths <- frame_paths[!is.na(frame_paths) & nzchar(frame_paths)]
  if (length(frame_paths) == 0) {
    NA_character_
  } else {
    normalizePath(tail(frame_paths, 1), winslash = "/", mustWork = FALSE)
  }
}

script_path <- get_sourced_script_path()
script_dir <- if (is.na(script_path)) getwd() else dirname(script_path)
project_root <- if (basename(script_dir) == "scripts") dirname(script_dir) else script_dir

possible_input_paths <- unique(c(
  Sys.getenv("BEECH_METADATA_FILE", unset = ""),
  file.path(script_dir, "donnees_F.grandifolia_2024.xlsx"),
  file.path(project_root, "donnees_F.grandifolia_2024.xlsx"),
  file.path(project_root, "genetique", "donnees_F.grandifolia_2024.xlsx"),
  file.path(project_root, "data", "donnees_F.grandifolia_2024.xlsx"),
  file.path(getwd(), "donnees_F.grandifolia_2024.xlsx"),
  file.path(getwd(), "genetique", "donnees_F.grandifolia_2024.xlsx"),
  file.path(getwd(), "data", "donnees_F.grandifolia_2024.xlsx")
))
possible_input_paths <- possible_input_paths[nzchar(possible_input_paths)]
matching_input_paths <- possible_input_paths[file.exists(possible_input_paths)]

if (length(matching_input_paths) == 0) {
  stop(
    "Could not find donnees_F.grandifolia_2024.xlsx.\n",
    "The script searched:\n- ",
    paste(possible_input_paths, collapse = "\n- "),
    "\n\nPut the workbook beside this script or set BEECH_METADATA_FILE."
  )
}

input_file <- normalizePath(
  matching_input_paths[1],
  winslash = "/",
  mustWork = TRUE
)

output_file <- file.path(
  project_root,
  "outputs",
  "michelle_beech_sample_selection.xlsx"
)
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

dbh_threshold_cm <- 10
minimum_spacing_km <- 8

as_numeric_clean <- function(x) {
  suppressWarnings(as.numeric(gsub(",", ".", as.character(x), fixed = TRUE)))
}

median_or_na <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) NA_real_ else median(x)
}

add_flag <- function(current, condition, flag) {
  ifelse(
    condition,
    ifelse(current == "", flag, paste(current, flag, sep = "; ")),
    current
  )
}

haversine_km <- function(lat1, lon1, lat2, lon2) {
  radius_km <- 6371.0088
  to_radians <- pi / 180
  dlat <- (lat2 - lat1) * to_radians
  dlon <- (lon2 - lon1) * to_radians
  a <- sin(dlat / 2)^2 +
    cos(lat1 * to_radians) * cos(lat2 * to_radians) * sin(dlon / 2)^2
  2 * radius_km * asin(sqrt(a))
}

genetique <- readxl::read_excel(input_file, sheet = "genetique") |>
  mutate(
    id_tige = as_numeric_clean(id_tige),
    dhp_tige = as_numeric_clean(dhp_tige),
    Latitude = as_numeric_clean(Latitude),
    Longitude = as_numeric_clean(Longitude),
    Elevation = as_numeric_clean(Elevation),
    valid_coordinate = !is.na(Latitude) &
      !is.na(Longitude) &
      dplyr::between(Latitude, 40, 85) &
      dplyr::between(Longitude, -145, -45),
    name_prefix_matches = startsWith(as.character(Name), paste0(site, "-")),
    name_suffix = ifelse(
      grepl("-[0-9]+$", as.character(Name)),
      as_numeric_clean(sub(".*-([0-9]+)$", "\\1", as.character(Name))),
      NA_real_
    ),
    qa_flag = ""
  )

genetique$qa_flag <- add_flag(
  genetique$qa_flag,
  is.na(genetique$dhp_tige),
  "Missing/non-numeric DBH"
)
genetique$qa_flag <- add_flag(
  genetique$qa_flag,
  !genetique$valid_coordinate,
  "Coordinate outside expected Canadian bounds"
)
genetique$qa_flag <- add_flag(
  genetique$qa_flag,
  !genetique$name_prefix_matches,
  "Site and tree-name prefix do not match"
)
genetique$qa_flag <- add_flag(
  genetique$qa_flag,
  !is.na(genetique$name_suffix) &
    !is.na(genetique$id_tige) &
    genetique$name_suffix != genetique$id_tige,
  "Tree-name number and id_tige do not match"
)

all_candidates <- genetique |>
  filter(!is.na(dhp_tige), dhp_tige > dbh_threshold_cm) |>
  arrange(site, desc(dhp_tige), id_tige) |>
  group_by(site) |>
  mutate(rank_within_site = row_number()) |>
  ungroup()

candidate_counts <- all_candidates |>
  group_by(site) |>
  summarise(
    trees_dbh_gt_10 = n(),
    valid_spatial_candidates = sum(valid_coordinate),
    maximum_dbh_cm = max(dhp_tige),
    .groups = "drop"
  )

best_choices <- all_candidates |>
  filter(valid_coordinate) |>
  group_by(site) |>
  arrange(desc(dhp_tige), id_tige, .by_group = TRUE) |>
  slice_head(n = 1) |>
  ungroup() |>
  left_join(
    candidate_counts |> select(site, trees_dbh_gt_10),
    by = "site"
  )

n_choices <- nrow(best_choices)
conflict_matrix <- matrix(FALSE, n_choices, n_choices)

if (n_choices > 1) {
  for (i in seq_len(n_choices - 1)) {
    for (j in (i + 1):n_choices) {
      distance <- haversine_km(
        best_choices$Latitude[i],
        best_choices$Longitude[i],
        best_choices$Latitude[j],
        best_choices$Longitude[j]
      )
      conflict_matrix[i, j] <- distance < minimum_spacing_km
      conflict_matrix[j, i] <- conflict_matrix[i, j]
    }
  }
}

# Exact branch-and-bound search: maximize the number of compatible sites.
# If several solutions use the same number of sites, prefer the larger total DBH.
best_indices <- integer(0)
best_total_dbh <- -Inf

search_compatible_sets <- function(position, chosen) {
  remaining <- n_choices - position + 1
  if (length(chosen) + max(remaining, 0) < length(best_indices)) return(invisible(NULL))
  
  if (position > n_choices) {
    total_dbh <- sum(best_choices$dhp_tige[chosen])
    if (
      length(chosen) > length(best_indices) ||
      (length(chosen) == length(best_indices) && total_dbh > best_total_dbh)
    ) {
      best_indices <<- chosen
      best_total_dbh <<- total_dbh
    }
    return(invisible(NULL))
  }
  
  compatible <- length(chosen) == 0 || all(!conflict_matrix[position, chosen])
  if (compatible) search_compatible_sets(position + 1, c(chosen, position))
  search_compatible_sets(position + 1, chosen)
  invisible(NULL)
}

if (n_choices > 0) search_compatible_sets(1, integer(0))
selected <- best_choices[best_indices, , drop = FALSE] |>
  arrange(site)

selected$key <- paste(selected$site, selected$Name, selected$id_tige, sep = "|")
selected_keys <- selected$key
selected_site_codes <- selected$site
all_candidates$key <- paste(
  all_candidates$site,
  all_candidates$Name,
  all_candidates$id_tige,
  sep = "|"
)

if (nrow(selected) > 1) {
  selected_distance_matrix <- outer(
    seq_len(nrow(selected)),
    seq_len(nrow(selected)),
    Vectorize(function(i, j) {
      if (i == j) return(0)
      haversine_km(
        selected$Latitude[i],
        selected$Longitude[i],
        selected$Latitude[j],
        selected$Longitude[j]
      )
    })
  )
  dimnames(selected_distance_matrix) <- list(selected$site, selected$site)
  
  selected$nearest_selected_site <- vapply(
    seq_len(nrow(selected)),
    function(i) {
      distances <- selected_distance_matrix[i, ]
      distances[i] <- Inf
      names(which.min(distances))
    },
    character(1)
  )
  selected$nearest_distance_km <- vapply(
    seq_len(nrow(selected)),
    function(i) {
      distances <- selected_distance_matrix[i, ]
      distances[i] <- Inf
      min(distances)
    },
    numeric(1)
  )
} else {
  selected_distance_matrix <- matrix(
    0,
    nrow = nrow(selected),
    ncol = nrow(selected),
    dimnames = list(selected$site, selected$site)
  )
  selected$nearest_selected_site <- NA_character_
  selected$nearest_distance_km <- NA_real_
}

selected_samples <- selected |>
  transmute(
    site,
    selected_tree = Name,
    id_tige,
    dbh_cm = dhp_tige,
    latitude = Latitude,
    longitude = Longitude,
    elevation_m = Elevation,
    candidate_count_at_site = trees_dbh_gt_10,
    nearest_selected_site,
    nearest_distance_km,
    spacing_pass = is.na(nearest_distance_km) |
      nearest_distance_km >= minimum_spacing_km,
    qa_flag
  )

all_candidates_output <- all_candidates |>
  mutate(
    selected = key %in% .env$selected_keys,
    selection_note = case_when(
      selected ~ "Largest valid-DBH tree at selected site",
      !valid_coordinate ~ "Coordinate must be corrected before spatial selection",
      TRUE ~ "Alternate candidate at same site"
    )
  ) |>
  transmute(
    site,
    rank_within_site,
    id_tige,
    tree_name = Name,
    dbh_cm = dhp_tige,
    latitude = Latitude,
    longitude = Longitude,
    elevation_m = Elevation,
    valid_coordinate,
    selected,
    selection_note,
    qa_flag
  )

site_summary <- genetique |>
  group_by(site) |>
  summarise(
    sampled_trees = n(),
    site_latitude_median = median_or_na(Latitude[valid_coordinate]),
    site_longitude_median = median_or_na(Longitude[valid_coordinate]),
    .groups = "drop"
  ) |>
  left_join(candidate_counts, by = "site") |>
  mutate(
    across(
      c(trees_dbh_gt_10, valid_spatial_candidates),
      ~ ifelse(is.na(.x), 0L, .x)
    ),
    selected = site %in% .env$selected_site_codes,
    reason = case_when(
      selected ~ "Selected",
      trees_dbh_gt_10 == 0 ~ "No sampled tree with DBH >10 cm",
      TRUE ~ "No valid-coordinate candidate or excluded by spacing"
    )
  ) |>
  select(
    site,
    sampled_trees,
    trees_dbh_gt_10,
    valid_spatial_candidates,
    maximum_dbh_cm,
    site_latitude_median,
    site_longitude_median,
    selected,
    reason
  ) |>
  arrange(site)

spacing_matrix <- if (nrow(selected) > 0) {
  data.frame(
    site = rownames(selected_distance_matrix),
    round(selected_distance_matrix, 3),
    check.names = FALSE,
    row.names = NULL
  )
} else {
  data.frame(site = character(0))
}

qa_flags <- genetique |>
  filter(qa_flag != "") |>
  transmute(
    site,
    id_tige,
    tree_name = Name,
    dbh_cm = dhp_tige,
    latitude = Latitude,
    longitude = Longitude,
    qa_flag,
    dbh_candidate = !is.na(dhp_tige) & dhp_tige > dbh_threshold_cm
  )

readme <- data.frame(
  item = c(
    "Source",
    "Worksheet",
    "DBH threshold",
    "Spacing threshold",
    "Selection rule",
    "Sites in source",
    "Trees with DBH >10 cm",
    "Selected sites",
    "Closest selected pair (km)",
    "Important limitation"
  ),
  value = c(
    input_file,
    "genetique",
    "Strictly greater than 10 cm",
    "At least 8 km",
    "One tree per site; maximize sites, then prefer larger DBH",
    dplyr::n_distinct(genetique$site),
    nrow(all_candidates),
    nrow(selected_samples),
    if (nrow(selected_samples) > 1) {
      round(min(selected_samples$nearest_distance_km), 3)
    } else {
      NA
    },
    "Spacing was checked only among these proposed trees; Michelle's existing coordinates were not supplied."
  ),
  stringsAsFactors = FALSE
)

writexl::write_xlsx(
  list(
    README = readme,
    Selected_samples = selected_samples,
    All_candidates = all_candidates_output,
    Site_summary = site_summary,
    Spacing_matrix = spacing_matrix,
    QA_flags = qa_flags
  ),
  path = output_file
)

message("Saved: ", normalizePath(output_file, winslash = "/", mustWork = FALSE))
message("Selected sites: ", nrow(selected_samples))
message("All DBH >10 cm candidates: ", nrow(all_candidates_output))