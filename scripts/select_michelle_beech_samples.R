#!/usr/bin/env Rscript

# Select samples for Michelle's landscape-genomics project without changing the
# source workbook. The source is discovered by worksheet name, not filename.

required_packages <- c("readxl", "openxlsx")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Missing required package(s): ", paste(missing_packages, collapse = ", "),
    ". Install them explicitly, then rerun; this script never installs packages.",
    call. = FALSE
  )
}

DBH_THRESHOLD_CM <- 10
MINIMUM_SPACING_KM <- 8
OUTPUT_NAME <- "michelle_beech_sample_selection.xlsx"

script_file <- function() {
  frame_files <- vapply(sys.frames(), function(x) {
    if (is.null(x$ofile)) NA_character_ else as.character(x$ofile)[1]
  }, character(1))
  frame_files <- frame_files[!is.na(frame_files) & nzchar(frame_files)]
  if (length(frame_files)) {
    return(normalizePath(
      tail(frame_files, 1),
      winslash = "/",
      mustWork = FALSE
    ))
  }
  args <- commandArgs(trailingOnly = FALSE)
  hit <- sub("^--file=", "", args[grepl("^--file=", args)])
  if (length(hit)) {
    normalizePath(hit[1], winslash = "/", mustWork = FALSE)
  } else {
    NA_character_
  }
}

find_project_root <- function() {
  starts <- unique(c(dirname(script_file()), getwd()))
  starts <- starts[!is.na(starts) & nzchar(starts)]
  
  for (start in starts) {
    probe <- normalizePath(start, winslash = "/", mustWork = FALSE)
    
    repeat {
      if (dir.exists(file.path(probe, ".git"))) {
        return(probe)
      }
      
      parent <- dirname(probe)
      if (identical(parent, probe)) {
        break
      }
      
      probe <- parent
    }
  }
  
  stop(
    "Could not locate the R-project root (the directory containing .git).",
    call. = FALSE
  )
}

project_root <- find_project_root()

normalize_name <- function(x) {
  x <- iconv(as.character(x), from = "", to = "ASCII//TRANSLIT")
  x <- tolower(x)
  gsub("(^_+|_+$)", "", gsub("[^a-z0-9]+", "_", x))
}

find_source_workbook <- function(root) {
  override <- Sys.getenv("BEECH_METADATA_FILE", unset = "")
  
  if (nzchar(override)) {
    override <- normalizePath(
      override,
      winslash = "/",
      mustWork = TRUE
    )
    
    sheets <- readxl::excel_sheets(override)
    hit <- which(normalize_name(sheets) == "genetique")
    
    if (!length(hit)) {
      stop(
        "BEECH_METADATA_FILE has no 'genetique' worksheet: ",
        override,
        call. = FALSE
      )
    }
    
    return(list(
      path = override,
      sheet = sheets[hit[1]]
    ))
  }
  
  files <- list.files(
    root,
    pattern = "\\.(xlsx|xls|xlsm)$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )
  
  files <- files[!grepl("(^|/)\\.git(/|$)", files)]
  
  matches <- lapply(files, function(path) {
    sheets <- tryCatch(
      readxl::excel_sheets(path),
      error = function(e) character()
    )
    
    hit <- which(normalize_name(sheets) == "genetique")
    
    if (length(hit)) {
      list(
        path = normalizePath(
          path,
          winslash = "/",
          mustWork = TRUE
        ),
        sheet = sheets[hit[1]]
      )
    } else {
      NULL
    }
  })
  
  matches <- Filter(Negate(is.null), matches)
  
  if (!length(matches)) {
    stop(
      "No Excel workbook containing a 'genetique' worksheet was found under ",
      root,
      ". Supply it in the project or set BEECH_METADATA_FILE.",
      call. = FALSE
    )
  }
  
  if (length(matches) > 1) {
    stop(
      "More than one workbook contains 'genetique'; set BEECH_METADATA_FILE:\n- ",
      paste(
        vapply(matches, `[[`, character(1), "path"),
        collapse = "\n- "
      ),
      call. = FALSE
    )
  }
  
  matches[[1]]
}

source_info <- find_source_workbook(project_root)

find_output_dir <- function(root) {
  preferred <- c("results", "outputs", "output")
  hit <- preferred[dir.exists(file.path(root, preferred))]
  
  if (!length(hit)) {
    stop(
      "No existing output/results folder was found at the project root.",
      call. = FALSE
    )
  }
  
  file.path(root, hit[1])
}

output_file <- file.path(
  find_output_dir(project_root),
  OUTPUT_NAME
)

pick_column <- function(df, aliases, label, required = TRUE) {
  nms <- names(df)
  normalized <- normalize_name(nms)
  hits <- which(normalized %in% normalize_name(aliases))
  
  if (!length(hits)) {
    if (!required) {
      return(NA_character_)
    }
    
    stop(
      "Could not identify the ",
      label,
      " column. Available columns: ",
      paste(nms, collapse = ", "),
      call. = FALSE
    )
  }
  
  if (length(hits) > 1) {
    stop(
      "Ambiguous ",
      label,
      " columns: ",
      paste(nms[hits], collapse = ", "),
      call. = FALSE
    )
  }
  
  nms[hits]
}

raw <- as.data.frame(
  readxl::read_excel(
    source_info$path,
    sheet = source_info$sheet,
    .name_repair = "minimal"
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

columns <- c(
  site = pick_column(
    raw,
    c("site", "site_id", "code_site", "population", "pop"),
    "site"
  ),
  tree_name = pick_column(
    raw,
    c(
      "name",
      "tree_name",
      "nom",
      "nom_tige",
      "sample_name",
      "individual"
    ),
    "tree/sample name",
    FALSE
  ),
  tree_id = pick_column(
    raw,
    c(
      "id_tige",
      "tree_id",
      "sample_id",
      "id_arbre",
      "numero_tige",
      "id"
    ),
    "tree/sample ID"
  ),
  dbh = pick_column(
    raw,
    c(
      "dhp_tige",
      "dbh",
      "dbh_cm",
      "dhp",
      "diametre",
      "diameter"
    ),
    "DBH"
  ),
  latitude = pick_column(
    raw,
    c("latitude", "lat", "y_wgs84"),
    "latitude"
  ),
  longitude = pick_column(
    raw,
    c("longitude", "lon", "long", "x_wgs84"),
    "longitude"
  ),
  elevation = pick_column(
    raw,
    c(
      "elevation",
      "elevation_m",
      "altitude",
      "altitude_m",
      "elev"
    ),
    "elevation",
    FALSE
  )
)

value <- function(name, default = NA_character_) {
  col <- columns[[name]]
  
  if (is.na(col)) {
    rep(default, nrow(raw))
  } else {
    raw[[col]]
  }
}

number <- function(x) {
  suppressWarnings(
    as.numeric(
      gsub(
        ",",
        ".",
        trimws(as.character(x)),
        fixed = TRUE
      )
    )
  )
}

text <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | x == ""] <- NA_character_
  x
}

append_flag <- function(x, condition, label) {
  condition[is.na(condition)] <- FALSE
  
  x[condition] <- ifelse(
    is.na(x[condition]) | x[condition] == "",
    label,
    paste(x[condition], label, sep = "; ")
  )
  
  x
}

dat <- data.frame(
  source_row = seq_len(nrow(raw)) + 1L,
  site = text(value("site")),
  tree_name = text(value("tree_name")),
  tree_id = text(value("tree_id")),
  dbh_cm = number(value("dbh")),
  latitude = number(value("latitude")),
  longitude = number(value("longitude")),
  elevation_m = number(value("elevation")),
  stringsAsFactors = FALSE
)

dat$missing_coordinate <- is.na(dat$latitude) |
  is.na(dat$longitude)

dat$suspicious_coordinate <- !dat$missing_coordinate & (
  !is.finite(dat$latitude) |
    !is.finite(dat$longitude) |
    abs(dat$latitude) > 90 |
    abs(dat$longitude) > 180 |
    dat$latitude == 0 |
    dat$longitude == 0
)

# Study sites are in eastern Canada. Values outside broad Canadian bounds are
# retained but flagged rather than silently removed or corrected.
outside_expected_bounds <- !dat$missing_coordinate &
  !(
    dat$latitude >= 40 &
      dat$latitude <= 85 &
      dat$longitude >= -145 &
      dat$longitude <= -45
  )

dat$suspicious_coordinate <- dat$suspicious_coordinate |
  outside_expected_bounds

dat$valid_coordinate <- !dat$missing_coordinate &
  !dat$suspicious_coordinate

site_key <- normalize_name(dat$site)
name_key <- normalize_name(dat$tree_name)

dat$site_name_disagreement <- !is.na(dat$tree_name) &
  !is.na(dat$site) &
  !startsWith(name_key, paste0(site_key, "_")) &
  name_key != site_key

name_number <- ifelse(
  grepl("[0-9]+$", name_key),
  sub(".*?([0-9]+)$", "\\1", name_key),
  NA_character_
)

id_number <- ifelse(
  grepl("[0-9]+$", normalize_name(dat$tree_id)),
  sub(
    ".*?([0-9]+)$",
    "\\1",
    normalize_name(dat$tree_id)
  ),
  NA_character_
)

dat$name_id_disagreement <- !is.na(name_number) &
  !is.na(id_number) &
  suppressWarnings(
    as.numeric(name_number) != as.numeric(id_number)
  )

coord_key <- ifelse(
  dat$valid_coordinate,
  paste(
    round(dat$latitude, 6),
    round(dat$longitude, 6)
  ),
  NA_character_
)

coord_sites <- ave(
  dat$site,
  coord_key,
  FUN = function(x) {
    length(unique(x[!is.na(x)]))
  }
)

duplicate_across_sites <- !is.na(coord_key) &
  coord_sites > 1

dat$qa_flag <- NA_character_

dat$qa_flag <- append_flag(
  dat$qa_flag,
  is.na(dat$site),
  "Missing site"
)

dat$qa_flag <- append_flag(
  dat$qa_flag,
  is.na(dat$tree_id),
  "Missing tree/sample ID"
)

dat$qa_flag <- append_flag(
  dat$qa_flag,
  is.na(dat$dbh_cm),
  "Missing/non-numeric DBH"
)

dat$qa_flag <- append_flag(
  dat$qa_flag,
  dat$missing_coordinate,
  "Missing coordinate"
)

dat$qa_flag <- append_flag(
  dat$qa_flag,
  dat$suspicious_coordinate,
  "Suspicious coordinate"
)

dat$qa_flag <- append_flag(
  dat$qa_flag,
  duplicate_across_sites,
  "Identical coordinate used by different sites"
)

dat$qa_flag <- append_flag(
  dat$qa_flag,
  dat$site_name_disagreement,
  "Site and tree name disagree"
)

dat$qa_flag <- append_flag(
  dat$qa_flag,
  dat$name_id_disagreement,
  "Tree name and tree ID disagree"
)

candidates <- dat[
  !is.na(dat$site) &
    !is.na(dat$dbh_cm) &
    dat$dbh_cm > DBH_THRESHOLD_CM,
]

candidates <- candidates[
  order(
    candidates$site,
    -candidates$dbh_cm,
    candidates$tree_id,
    candidates$source_row
  ),
]

candidates$rank_within_site <- ave(
  seq_len(nrow(candidates)),
  candidates$site,
  FUN = seq_along
)

haversine_km <- function(lat1, lon1, lat2, lon2) {
  rad <- pi / 180
  dlat <- (lat2 - lat1) * rad
  dlon <- (lon2 - lon1) * rad
  
  a <- sin(dlat / 2)^2 +
    cos(lat1 * rad) *
    cos(lat2 * rad) *
    sin(dlon / 2)^2
  
  2 * 6371.0088 * atan2(
    sqrt(a),
    sqrt(pmax(0, 1 - a))
  )
}

# Choose among all coordinate-valid candidates. This exact branch-and-bound
# search maximizes represented sites first, then total DBH (therefore preferring
# each site's largest tree whenever that does not reduce the number of sites).
spatial <- candidates[candidates$valid_coordinate, ]

sites <- sort(unique(candidates$site))

choices <- lapply(
  sites,
  function(s) {
    which(spatial$site == s)
  }
)

names(choices) <- sites

best <- integer()
best_dbh <- -Inf

search <- function(position, chosen) {
  if (
    length(chosen) +
    length(sites) -
    position +
    1L <
    length(best)
  ) {
    return(invisible(NULL))
  }
  
  if (position > length(sites)) {
    score <- sum(spatial$dbh_cm[chosen])
    
    if (
      length(chosen) > length(best) ||
      (
        length(chosen) == length(best) &&
        score > best_dbh
      )
    ) {
      best <<- chosen
      best_dbh <<- score
    }
    
    return(invisible(NULL))
  }
  
  for (candidate in choices[[position]]) {
    compatible <- !length(chosen) ||
      all(
        haversine_km(
          spatial$latitude[candidate],
          spatial$longitude[candidate],
          spatial$latitude[chosen],
          spatial$longitude[chosen]
        ) >= MINIMUM_SPACING_KM
      )
    
    if (compatible) {
      search(
        position + 1L,
        c(chosen, candidate)
      )
    }
  }
  
  search(position + 1L, chosen)
  invisible(NULL)
}

if (length(sites)) {
  search(1L, integer())
}

selected <- spatial[best, , drop = FALSE]

selected <- selected[
  order(selected$site),
]

sample_key <- function(x) {
  paste(
    x$source_row,
    x$site,
    x$tree_id,
    sep = "|"
  )
}

selected_keys <- sample_key(selected)

candidates$selected <- sample_key(candidates) %in%
  selected_keys

candidates$selection_note <- ifelse(
  candidates$selected,
  "Selected",
  ifelse(
    !candidates$valid_coordinate,
    "Eligible alternate; coordinate requires verification",
    ifelse(
      candidates$site %in% selected$site,
      "Eligible alternate at selected site",
      "Eligible candidate at site excluded by 8-km optimization"
    )
  )
)

distance_matrix <- matrix(
  numeric(),
  nrow = nrow(selected),
  ncol = nrow(selected),
  dimnames = list(
    selected$site,
    selected$site
  )
)

if (nrow(selected)) {
  for (i in seq_len(nrow(selected))) {
    for (j in seq_len(nrow(selected))) {
      distance_matrix[i, j] <- if (i == j) {
        0
      } else {
        haversine_km(
          selected$latitude[i],
          selected$longitude[i],
          selected$latitude[j],
          selected$longitude[j]
        )
      }
    }
  }
}

if (nrow(selected) > 1) {
  nearest <- apply(
    distance_matrix + diag(Inf, nrow(selected)),
    1,
    which.min
  )
  
  selected$nearest_selected_site <-
    colnames(distance_matrix)[nearest]
  
  selected$nearest_distance_km <-
    distance_matrix[
      cbind(
        seq_len(nrow(selected)),
        nearest
      )
    ]
} else {
  selected$nearest_selected_site <- NA_character_
  selected$nearest_distance_km <- NA_real_
}

selected_samples <- selected[
  ,
  c(
    "site",
    "tree_name",
    "tree_id",
    "dbh_cm",
    "latitude",
    "longitude",
    "elevation_m",
    "nearest_selected_site",
    "nearest_distance_km",
    "qa_flag"
  )
]

names(selected_samples)[2] <- "selected_tree"

selected_samples$spacing_pass <-
  is.na(selected_samples$nearest_distance_km) |
  selected_samples$nearest_distance_km >= MINIMUM_SPACING_KM

all_candidates <- candidates[
  ,
  c(
    "source_row",
    "site",
    "rank_within_site",
    "tree_id",
    "tree_name",
    "dbh_cm",
    "latitude",
    "longitude",
    "elevation_m",
    "valid_coordinate",
    "selected",
    "selection_note",
    "qa_flag"
  )
]

site_summary <- do.call(
  rbind,
  lapply(sites, function(s) {
    x <- candidates[candidates$site == s, ]
    
    data.frame(
      site = s,
      eligible_trees = nrow(x),
      valid_spatial_candidates = sum(x$valid_coordinate),
      maximum_dbh_cm = max(x$dbh_cm),
      selected = s %in% selected$site,
      reason = if (s %in% selected$site) {
        "Selected"
      } else if (!any(x$valid_coordinate)) {
        "No valid-coordinate eligible tree"
      } else {
        "Excluded by exact 8-km spacing optimization"
      }
    )
  })
)

if (is.null(site_summary)) {
  site_summary <- data.frame(
    site = character(),
    eligible_trees = integer(),
    valid_spatial_candidates = integer(),
    maximum_dbh_cm = numeric(),
    selected = logical(),
    reason = character()
  )
}

spacing_matrix <- if (nrow(selected)) {
  data.frame(
    site = rownames(distance_matrix),
    round(distance_matrix, 3),
    check.names = FALSE,
    row.names = NULL
  )
} else {
  data.frame(site = character())
}

qa_flags <- dat[
  !is.na(dat$qa_flag),
  c(
    "source_row",
    "site",
    "tree_name",
    "tree_id",
    "dbh_cm",
    "latitude",
    "longitude",
    "elevation_m",
    "qa_flag"
  )
]

qa_flags$eligible_dbh <- !is.na(qa_flags$dbh_cm) &
  qa_flags$dbh_cm > DBH_THRESHOLD_CM

closest <- if (nrow(selected) > 1) {
  z <- distance_matrix
  diag(z) <- Inf
  ij <- arrayInd(which.min(z), dim(z))
  
  paste(
    rownames(z)[ij[1]],
    "-",
    colnames(z)[ij[2]],
    sprintf("(%.3f km)", z[ij])
  )
} else {
  "Not applicable"
}

readme <- data.frame(
  item = c(
    "Source workbook",
    "Source worksheet",
    "Detected columns",
    "DBH rule",
    "Spacing rule",
    "Optimization",
    "Eligible trees",
    "Eligible sites",
    "Selected trees",
    "Closest selected pair",
    "QA scope",
    "Original data"
  ),
  value = c(
    source_info$path,
    source_info$sheet,
    paste(
      names(columns),
      columns,
      sep = "=",
      collapse = "; "
    ),
    "DBH strictly greater than 10 cm",
    "At least 8 km great-circle distance between selected coordinates",
    paste(
      "Exact search: maximize sites, then maximize total DBH;",
      "maximum one tree per site"
    ),
    nrow(candidates),
    length(sites),
    nrow(selected),
    closest,
    paste(
      "Flags missing/suspicious coordinates, cross-site coordinate",
      "duplicates, and site/name/ID disagreements"
    ),
    "Read-only; the source workbook is never written or modified"
  ),
  stringsAsFactors = FALSE
)

openxlsx::write.xlsx(
  list(
    Selected_samples = selected_samples,
    All_candidates = all_candidates,
    Site_summary = site_summary,
    Spacing_matrix = spacing_matrix,
    QA_flags = qa_flags,
    README = readme
  ),
  file = output_file,
  overwrite = TRUE,
  asTable = TRUE
)

message("Input file: ", source_info$path)

message(
  "Output file: ",
  normalizePath(
    output_file,
    winslash = "/",
    mustWork = TRUE
  )
)

message(
  "Eligible trees: ",
  nrow(candidates),
  "; eligible sites: ",
  length(sites),
  "; selected trees: ",
  nrow(selected)
)

message("Closest selected pair: ", closest)