#!/usr/bin/env Rscript

# Select spatially separated Fagus grandifolia trees for Michelle's samples.
# Run with Rscript from the BeechCode project root, or source it from anywhere.
# The project can also be specified with BEECHCODE_ROOT=/path/to/beechcode.

required_packages <- c("readxl", "openxlsx")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Install required package(s): ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

workbook_relative_path <- file.path("data", "raw", "donnees_F.grandifolia_2024.xlsx")
input_sheet <- "genetique"
minimum_dbh_cm <- 10
minimum_distance_km <- 8

ancestor_paths <- function(path) {
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  out <- character()
  repeat {
    out <- c(out, path)
    parent <- dirname(path)
    if (identical(parent, path)) break
    path <- parent
  }
  out
}

source_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
source_dir <- if (!is.null(source_file) && length(source_file) == 1L && nzchar(source_file)) {
  dirname(normalizePath(source_file, winslash = "/", mustWork = FALSE))
} else character()
desktop_dir <- file.path(path.expand("~"), "Desktop")
desktop_projects <- if (dir.exists(desktop_dir)) {
  children <- list.dirs(desktop_dir, recursive = FALSE, full.names = TRUE)
  children[tolower(basename(children)) == "beechcode"]
} else character()

configured_root <- Sys.getenv("BEECHCODE_ROOT", unset = "")
if (nzchar(configured_root)) {
  root_candidates <- configured_root
} else {
  # Search in a documented order. In particular, check ~/Desktop/beechcode
  # before the working directory so source("~/Desktop/...") works even when R
  # was launched from another project.
  root_candidates <- c(
    desktop_projects,
    if (length(source_dir)) ancestor_paths(source_dir) else character(),
    ancestor_paths(getwd()),
    file.path(getwd(), c("beechcode", "BeechCode"))
  )
}
root_candidates <- root_candidates[nzchar(root_candidates)]
root_candidates <- unique(normalizePath(root_candidates, winslash = "/", mustWork = FALSE))
workbook_candidates <- file.path(root_candidates, workbook_relative_path)
matches <- which(file.exists(workbook_candidates))

if (!length(matches)) {
  stop(
    "Input workbook not found. Expected 'data/raw/donnees_F.grandifolia_2024.xlsx' ",
    "under the BeechCode project. Checked:\n- ",
    paste(workbook_candidates, collapse = "\n- "),
    "\nSet BEECHCODE_ROOT to the project directory if it is elsewhere.",
    call. = FALSE
  )
}

# The first match follows the precedence above; the resolved path is printed at
# the end so the selected input is always visible rather than silent.
match_index <- matches[1L]
project_root <- normalizePath(root_candidates[match_index], winslash = "/", mustWork = TRUE)
input_path <- normalizePath(workbook_candidates[match_index], winslash = "/", mustWork = TRUE)
output_dir <- file.path(project_root, "genetique", "outputs", "michelle_beech_samples")
output_path <- file.path(output_dir, "michelle_beech_sample_selection.xlsx")

sheets <- readxl::excel_sheets(input_path)
if (!(input_sheet %in% sheets)) {
  stop(
    "Required worksheet '", input_sheet, "' not found. Available worksheets: ",
    paste(sheets, collapse = ", "), call. = FALSE
  )
}

raw <- as.data.frame(
  readxl::read_excel(input_path, sheet = input_sheet, .name_repair = "minimal"),
  stringsAsFactors = FALSE, check.names = FALSE
)

# Normalization is used only to match column headings. Data values are preserved.
heading_key <- function(x) {
  x <- iconv(x, from = "", to = "ASCII//TRANSLIT")
  tolower(gsub("[^a-z0-9]+", "", x))
}

resolve_column <- function(aliases, label) {
  keys <- heading_key(names(raw))
  hits <- which(keys %in% heading_key(aliases))
  if (length(hits) != 1L) {
    stop(
      "Expected exactly one ", label, " column; matched ", length(hits),
      ". Available columns: ", paste(names(raw), collapse = ", "),
      call. = FALSE
    )
  }
  names(raw)[hits]
}

site_col <- resolve_column(c("site", "site_id", "code_site", "population", "pop"), "site")
id_col <- resolve_column(
  c("id", "id_tige", "identifiant", "individual", "individual_id", "individu",
    "numero", "numero_arbre", "numero_individu", "sample", "sample_id", "tree"),
  "tree identifier"
)
dbh_col <- resolve_column(c("Dhp_tige", "dbh", "dbh_cm", "dhp", "diametre", "diametre_cm"), "DBH")
lat_col <- resolve_column(c("latitude", "lat", "y_wgs84"), "latitude")
lon_col <- resolve_column(c("longitude", "lon", "long", "x_wgs84"), "longitude")

numeric_exact <- function(x) suppressWarnings(as.numeric(as.character(x)))
site_value <- as.character(raw[[site_col]])
id_value <- as.character(raw[[id_col]])
dbh_value <- numeric_exact(raw[[dbh_col]])
lat_value <- numeric_exact(raw[[lat_col]])
lon_value <- numeric_exact(raw[[lon_col]])

flag_rows <- function(issue, test, severity = "warning") {
  rows <- which(test %in% TRUE)
  if (!length(rows)) return(NULL)
  data.frame(
    source_row = rows + 1L,
    tree_identifier = id_value[rows],
    site = site_value[rows],
    severity = severity,
    issue = issue,
    stringsAsFactors = FALSE
  )
}

qa_flags <- do.call(rbind, Filter(Negate(is.null), list(
  flag_rows("missing site", is.na(site_value) | site_value == "", "error"),
  flag_rows("missing tree identifier", is.na(id_value) | id_value == "", "error"),
  flag_rows("identifier has leading/trailing whitespace (not corrected)",
            !is.na(id_value) & id_value != trimws(id_value)),
  flag_rows("identifier contains a control character (not corrected)",
            !is.na(id_value) & grepl("[[:cntrl:]]", id_value)),
  flag_rows("duplicate tree identifier (not corrected)",
            !is.na(id_value) & (duplicated(id_value) | duplicated(id_value, fromLast = TRUE))),
  flag_rows("DBH is missing or non-numeric", is.na(dbh_value)),
  flag_rows("latitude is missing/non-numeric or outside [-90, 90] (not corrected)",
            is.na(lat_value) | lat_value < -90 | lat_value > 90, "error"),
  flag_rows("longitude is missing/non-numeric or outside [-180, 180] (not corrected)",
            is.na(lon_value) | lon_value < -180 | lon_value > 180, "error")
)))
if (is.null(qa_flags)) {
  qa_flags <- data.frame(source_row = integer(), tree_identifier = character(),
                         site = character(), severity = character(), issue = character())
}

candidate_index <- which(!is.na(dbh_value) & dbh_value > minimum_dbh_cm)
candidates <- raw[candidate_index, , drop = FALSE]
candidates$source_row <- candidate_index + 1L
candidates$.site <- site_value[candidate_index]
candidates$.tree_identifier <- id_value[candidate_index]
candidates$.dbh_cm <- dbh_value[candidate_index]
candidates$.latitude <- lat_value[candidate_index]
candidates$.longitude <- lon_value[candidate_index]

if (!nrow(candidates)) stop("No rows have DBH strictly greater than 10 cm.", call. = FALSE)

# Retain every DBH-qualified row in all_candidates. Rows whose exact, unaltered
# site/coordinate values cannot support spacing are retained there and flagged,
# but cannot be considered for selection.
selectable <- candidates[
  !is.na(candidates$.site) & candidates$.site != "" &
    !is.na(candidates$.latitude) & candidates$.latitude >= -90 & candidates$.latitude <= 90 &
    !is.na(candidates$.longitude) & candidates$.longitude >= -180 & candidates$.longitude <= 180,
  , drop = FALSE
]
if (!nrow(selectable)) stop("No DBH-qualified rows have a site and valid coordinates.", call. = FALSE)

# Pick the largest valid-DBH tree at each site. Exact source order is the explicit
# deterministic tie-break, without changing identifiers or coordinates.
candidate_order <- order(selectable$.site, -selectable$.dbh_cm, selectable$source_row)
ranked <- selectable[candidate_order, , drop = FALSE]
representatives <- ranked[!duplicated(ranked$.site), , drop = FALSE]
row.names(representatives) <- NULL

haversine_km <- function(lat1, lon1, lat2, lon2) {
  rad <- pi / 180
  dlat <- (lat2 - lat1) * rad
  dlon <- (lon2 - lon1) * rad
  a <- sin(dlat / 2)^2 + cos(lat1 * rad) * cos(lat2 * rad) * sin(dlon / 2)^2
  6371.0088 * 2 * atan2(sqrt(a), sqrt(pmax(0, 1 - a)))
}

n_sites <- nrow(representatives)
distances <- matrix(0, n_sites, n_sites,
                    dimnames = list(representatives$.site, representatives$.site))
if (n_sites > 1L) {
  for (i in seq_len(n_sites - 1L)) for (j in (i + 1L):n_sites) {
    distances[i, j] <- distances[j, i] <- haversine_km(
      representatives$.latitude[i], representatives$.longitude[i],
      representatives$.latitude[j], representatives$.longitude[j]
    )
  }
}
conflicts <- distances < minimum_distance_km & row(distances) != col(distances)

# Exact branch-and-bound maximum independent set. Thus a conflict never causes
# a greedy, order-dependent loss of a compatible site.
best <- integer()
search_independent_set <- function(available, chosen = integer()) {
  if (length(chosen) + length(available) <= length(best)) return(invisible(NULL))
  if (!length(available)) {
    if (length(chosen) > length(best)) best <<- chosen
    return(invisible(NULL))
  }
  vertex <- available[1L]
  search_independent_set(available[!conflicts[vertex, available] & available != vertex],
                         c(chosen, vertex))
  search_independent_set(available[-1L], chosen)
  invisible(NULL)
}
search_independent_set(seq_len(n_sites))
selected <- representatives[best, , drop = FALSE]

selected_sites <- selected$.site
site_summary <- data.frame(
  site = sort(unique(candidates$.site)),
  stringsAsFactors = FALSE
)
site_summary$candidate_trees <- as.integer(table(candidates$.site)[site_summary$site])
site_summary$largest_valid_dbh_cm <- vapply(site_summary$site, function(s) {
  max(candidates$.dbh_cm[candidates$.site == s])
}, numeric(1))
site_summary$preferred_tree_identifier <- representatives$.tree_identifier[
  match(site_summary$site, representatives$.site)
]
site_summary$selected <- site_summary$site %in% selected_sites

spacing_matrix <- data.frame(site = rownames(distances), round(distances, 3),
                             check.names = FALSE, row.names = NULL)
drop_helpers <- function(x) x[, !startsWith(names(x), "."), drop = FALSE]

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
openxlsx::write.xlsx(
  list(
    selected_samples = drop_helpers(selected),
    all_candidates = drop_helpers(candidates),
    site_summary = site_summary,
    spacing_km = spacing_matrix,
    qa_flags = qa_flags
  ),
  file = output_path, overwrite = TRUE
)

selected_minimum <- if (nrow(selected) > 1L) {
  min(distances[best, best][upper.tri(distances[best, best])])
} else {
  NA_real_
}
stopifnot(!any(duplicated(selected$.site)))
stopifnot(is.na(selected_minimum) || selected_minimum >= minimum_distance_km)

message("Input: ", input_path, " [worksheet: ", input_sheet, "]")
message("Candidate trees (DBH > ", minimum_dbh_cm, " cm): ", nrow(candidates))
message("Eligible sites: ", n_sites)
message("Selected trees/sites: ", nrow(selected))
message("Minimum selected-tree distance: ", round(selected_minimum, 3), " km")
message("QA flags (values were not corrected): ", nrow(qa_flags))
message("Output: ", output_path)