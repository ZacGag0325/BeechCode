#!/usr/bin/env Rscript

# Download and summarize 1981-2010 Environment and Climate Change Canada
# Climate Normals for two American beech study-area climate clusters in Quebec.
#
# The original target climate IDs are used only to locate the study-area cluster
# centroids. The script then selects the nearest Quebec station with 1981-2010
# Climate Normals and downloads normals with weathercan::normals_dl(). It never
# falls back to incomplete daily records.

required_packages <- c(
  "weathercan",
  "dplyr",
  "readr",
  "stringr",
  "lubridate",
  "tibble",
  "tidyr"
)

# weathercan is distributed through the rOpenSci R-universe rather than CRAN
# for some R versions, so include both repositories when installing packages.
package_repos <- c(
  ropensci = "https://ropensci.r-universe.dev",
  CRAN = "https://cloud.r-project.org"
)

install_if_missing <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    message("Installing missing package: ", package)
    install.packages(package, repos = package_repos)
  }
  
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(
      "Package '", package, "' could not be installed. ",
      "For weathercan, try running: install.packages('weathercan', repos = c('https://ropensci.r-universe.dev', 'https://cloud.r-project.org'))",
      call. = FALSE
    )
  }
}

invisible(lapply(required_packages, install_if_missing))
invisible(lapply(required_packages, library, character.only = TRUE))

normals_period <- "1981-2010"

clusters <- tibble::tribble(
  ~cluster, ~centroid_source_station, ~centroid_source_climate_id,
  "northern", "STE ANNE DU LAC", "7036855",
  "southern", "NOTRE DAME DE LA PAIX", "7035666"
)

raw_dir <- "data/raw/climate"
output_dir <- "outputs/climate"
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

message("Updating weathercan station inventory with weathercan::stations_dl()...")
tryCatch(
  weathercan::stations_dl(verbose = FALSE, quiet = FALSE),
  error = function(e) {
    stop(
      "Could not update the weathercan station inventory with stations_dl(): ",
      conditionMessage(e),
      call. = FALSE
    )
  }
)

station_inventory <- weathercan::stations() |>
  dplyr::mutate(
    climate_id = as.character(.data$climate_id),
    prov = as.character(.data$prov),
    station_name = as.character(.data$station_name),
    normals_1981_2010 = as.logical(.data$normals_1981_2010)
  )

haversine_km <- function(lat1, lon1, lat2, lon2) {
  earth_radius_km <- 6371.0088
  to_radians <- pi / 180
  
  lat1 <- lat1 * to_radians
  lon1 <- lon1 * to_radians
  lat2 <- lat2 * to_radians
  lon2 <- lon2 * to_radians
  
  dlat <- lat2 - lat1
  dlon <- lon2 - lon1
  
  a <- sin(dlat / 2)^2 + cos(lat1) * cos(lat2) * sin(dlon / 2)^2
  2 * earth_radius_km * atan2(sqrt(a), sqrt(1 - a))
}

nearest_alternatives <- function(cluster_row, n = 10) {
  station_inventory |>
    dplyr::filter(
      .data$prov == "QC",
      .data$normals_1981_2010,
      !is.na(.data$lat),
      !is.na(.data$lon)
    ) |>
    dplyr::distinct(.data$climate_id, .keep_all = TRUE) |>
    dplyr::mutate(
      distance_from_centroid_km = haversine_km(
        cluster_row$centroid_lat,
        cluster_row$centroid_lon,
        .data$lat,
        .data$lon
      )
    ) |>
    dplyr::arrange(.data$distance_from_centroid_km) |>
    dplyr::slice_head(n = n)
}

get_cluster_centroids <- function() {
  dplyr::bind_rows(lapply(seq_len(nrow(clusters)), function(i) {
    cluster_row <- clusters[i, ]
    
    centroid_station <- station_inventory |>
      dplyr::filter(.data$climate_id == cluster_row$centroid_source_climate_id) |>
      dplyr::filter(!is.na(.data$lat), !is.na(.data$lon)) |>
      dplyr::arrange(.data$start) |>
      dplyr::slice(1)
    
    if (nrow(centroid_station) == 0) {
      stop(
        "Could not locate centroid source climate ID ",
        cluster_row$centroid_source_climate_id,
        " (", cluster_row$centroid_source_station, ") in the station inventory.",
        call. = FALSE
      )
    }
    
    dplyr::bind_cols(
      cluster_row,
      tibble::tibble(
        centroid_lat = centroid_station$lat,
        centroid_lon = centroid_station$lon,
        centroid_elev_m = centroid_station$elev
      )
    )
  }))
}

select_nearest_normals_station <- function(cluster_row) {
  alternatives <- nearest_alternatives(cluster_row, n = 10)
  
  if (nrow(alternatives) == 0) {
    stop(
      "No Quebec stations with 1981-2010 Climate Normals were found near the ",
      cluster_row$cluster,
      " cluster.",
      call. = FALSE
    )
  }
  
  selected <- alternatives |>
    dplyr::slice(1)
  
  if (!isTRUE(selected$normals_1981_2010)) {
    alternative_text <- alternatives |>
      dplyr::transmute(
        alternative = stringr::str_glue(
          "{station_name} ({climate_id}), {round(distance_from_centroid_km, 1)} km"
        )
      ) |>
      dplyr::pull(.data$alternative) |>
      paste(collapse = "; ")
    
    stop(
      "The nearest selected station for the ", cluster_row$cluster,
      " cluster does not have 1981-2010 Climate Normals. Nearest alternatives with normals: ",
      alternative_text,
      call. = FALSE
    )
  }
  
  dplyr::bind_cols(
    cluster_row,
    selected |>
      dplyr::select(
        selected_station_name = .data$station_name,
        selected_climate_id = .data$climate_id,
        selected_station_id = .data$station_id,
        selected_lat = .data$lat,
        selected_lon = .data$lon,
        selected_elev_m = .data$elev,
        distance_from_centroid_km = .data$distance_from_centroid_km,
        normals_1981_2010 = .data$normals_1981_2010
      )
  )
}

cluster_centroids <- get_cluster_centroids()
selected_stations <- dplyr::bind_rows(lapply(seq_len(nrow(cluster_centroids)), function(i) {
  select_nearest_normals_station(cluster_centroids[i, ])
})) |>
  dplyr::mutate(normals_period = normals_period)

if (any(is.na(selected_stations$selected_climate_id))) {
  stop("At least one cluster did not receive a valid selected Climate Normals station.", call. = FALSE)
}

if (dplyr::n_distinct(selected_stations$selected_climate_id) < nrow(selected_stations)) {
  duplicate_station_text <- selected_stations |>
    dplyr::count(.data$selected_climate_id, .data$selected_station_name) |>
    dplyr::filter(.data$n > 1) |>
    dplyr::transmute(station = stringr::str_glue("{selected_station_name} ({selected_climate_id})")) |>
    dplyr::pull(.data$station) |>
    paste(collapse = "; ")
  
  message(
    "Note: more than one cluster intentionally selected the same nearest normals station: ",
    duplicate_station_text
  )
}

message("Selected 1981-2010 Climate Normals stations:")
for (i in seq_len(nrow(selected_stations))) {
  message(
    stringr::str_glue(
      "- {selected_stations$cluster[i]} cluster: {selected_stations$selected_station_name[i]} ",
      "(Climate ID {selected_stations$selected_climate_id[i]}), ",
      "lat {round(selected_stations$selected_lat[i], 4)}, lon {round(selected_stations$selected_lon[i], 4)}, ",
      "elev {round(selected_stations$selected_elev_m[i], 1)} m, ",
      "{round(selected_stations$distance_from_centroid_km[i], 1)} km from centroid, ",
      "normals period {selected_stations$normals_period[i]}"
    )
  )
}

readr::write_csv(
  selected_stations,
  file.path(raw_dir, "selected_climate_normals_stations.csv"),
  na = ""
)

message("Downloading 1981-2010 Climate Normals with weathercan::normals_dl()...")
normals_raw <- weathercan::normals_dl(
  climate_ids = unique(selected_stations$selected_climate_id),
  normals_years = normals_period,
  verbose = FALSE,
  quiet = FALSE
)

normals_flat <- normals_raw |>
  dplyr::select(-dplyr::any_of("frost")) |>
  tidyr::unnest(.data$normals)

readr::write_csv(
  normals_flat,
  file.path(raw_dir, "climate_normals_1981_2010_raw.csv"),
  na = ""
)

find_column <- function(data, candidates, required = TRUE) {
  available <- names(data)
  exact <- candidates[candidates %in% available]
  if (length(exact) > 0) {
    return(exact[1])
  }
  
  for (pattern in candidates) {
    regex_matches <- available[stringr::str_detect(available, stringr::regex(pattern, ignore_case = TRUE))]
    if (length(regex_matches) > 0) {
      return(regex_matches[1])
    }
  }
  
  if (required) {
    stop(
      "Could not find any of these columns in the downloaded normals data: ",
      paste(candidates, collapse = ", "),
      call. = FALSE
    )
  }
  
  NA_character_
}

month_column <- find_column(normals_flat, c("month", "period", "date"), required = FALSE)
if (is.na(month_column)) {
  stop("Could not identify the month/period column in the downloaded normals data.", call. = FALSE)
}

annual_normals <- normals_flat |>
  dplyr::mutate(normals_row_label = as.character(.data[[month_column]])) |>
  dplyr::filter(stringr::str_detect(stringr::str_to_lower(.data$normals_row_label), "annual|year"))

if (nrow(annual_normals) == 0) {
  stop(
    "The downloaded normals did not include an annual/year row. ",
    "Inspect data/raw/climate/climate_normals_1981_2010_raw.csv for available rows.",
    call. = FALSE
  )
}

temperature_col <- find_column(
  annual_normals,
  c("mean_temp", "mean_temperature", "daily_mean_temp", "^mean.*temp", "temp.*mean")
)
precipitation_col <- find_column(
  annual_normals,
  c("total_precip", "total_precipitation", "precipitation", "^total.*precip", "precip.*total")
)
snowfall_col <- find_column(
  annual_normals,
  c("total_snow", "snowfall", "snow_precip", "snow_precipitation", "^total.*snow", "snow.*precip"),
  required = FALSE
)
degree_days_col <- find_column(
  annual_normals,
  c("degree_days_above_0", "degree_days_above_0_c", "growing_degree_days", "degree.*above.*0", "above.*0.*degree"),
  required = FALSE
)

value_or_na <- function(data, column) {
  if (is.na(column)) {
    return(rep(NA_real_, nrow(data)))
  }
  as.numeric(data[[column]])
}

summary_by_cluster <- selected_stations |>
  dplyr::left_join(
    annual_normals,
    by = c("selected_climate_id" = "climate_id")
  ) |>
  dplyr::mutate(
    mean_annual_temperature_c = as.numeric(.data[[temperature_col]]),
    mean_annual_precipitation_mm = as.numeric(.data[[precipitation_col]]),
    mean_annual_snowfall_cm = value_or_na(dplyr::cur_data_all(), snowfall_col),
    mean_degree_days_above_0_c = value_or_na(dplyr::cur_data_all(), degree_days_col),
    snowfall_variable_used = ifelse(is.na(snowfall_col), NA_character_, snowfall_col),
    degree_days_variable_used = ifelse(is.na(degree_days_col), NA_character_, degree_days_col)
  ) |>
  dplyr::select(
    .data$cluster,
    .data$centroid_source_station,
    .data$centroid_source_climate_id,
    .data$centroid_lat,
    .data$centroid_lon,
    .data$selected_station_name,
    .data$selected_climate_id,
    .data$selected_station_id,
    .data$selected_lat,
    .data$selected_lon,
    .data$selected_elev_m,
    .data$distance_from_centroid_km,
    .data$normals_period,
    .data$mean_annual_temperature_c,
    .data$mean_annual_precipitation_mm,
    .data$mean_annual_snowfall_cm,
    .data$mean_degree_days_above_0_c,
    .data$snowfall_variable_used,
    .data$degree_days_variable_used
  ) |>
  dplyr::arrange(.data$cluster)

if (any(is.na(summary_by_cluster$mean_annual_temperature_c)) ||
    any(is.na(summary_by_cluster$mean_annual_precipitation_mm))) {
  stop(
    "At least one selected station is missing annual temperature or precipitation normals; ",
    "inspect the raw normals CSV before using these summaries.",
    call. = FALSE
  )
}

if (nrow(summary_by_cluster) != nrow(clusters)) {
  stop("The summary table does not have exactly one row per cluster.", call. = FALSE)
}

if (dplyr::n_distinct(summary_by_cluster$selected_climate_id) == nrow(summary_by_cluster) &&
    any(duplicated(summary_by_cluster[, c("mean_annual_temperature_c", "mean_annual_precipitation_mm")])) ) {
  stop(
    "Distinct selected stations produced duplicate summary values. ",
    "Check that the downloaded normals data were not accidentally reused across clusters.",
    call. = FALSE
  )
}

format_optional_value <- function(value, suffix, missing_text) {
  if (is.na(value)) {
    return(missing_text)
  }
  paste0(round(value, 0), suffix)
}

sentences <- summary_by_cluster |>
  dplyr::mutate(
    snowfall_text = vapply(
      .data$mean_annual_snowfall_cm,
      format_optional_value,
      character(1),
      suffix = " cm",
      missing_text = "not available in the downloaded normals"
    ),
    degree_days_text = vapply(
      .data$mean_degree_days_above_0_c,
      format_optional_value,
      character(1),
      suffix = "",
      missing_text = "not available in the downloaded normals"
    ),
    sentence = stringr::str_glue(
      "For the {cluster} cluster, the nearest Quebec station with {normals_period} ",
      "Climate Normals was {selected_station_name} (Climate ID {selected_climate_id}; ",
      "{round(distance_from_centroid_km, 1)} km from the cluster centroid; ",
      "lat {round(selected_lat, 4)}, lon {round(selected_lon, 4)}, ",
      "elevation {round(selected_elev_m, 1)} m). Mean annual temperature was ",
      "{round(mean_annual_temperature_c, 1)}°C, mean annual precipitation was ",
      "{round(mean_annual_precipitation_mm, 0)} mm, mean annual snowfall was ",
      "{snowfall_text}, and degree-days above 0°C were {degree_days_text}."
    )
  ) |>
  dplyr::pull(.data$sentence)

thesis_paragraph <- paste(sentences, collapse = " ")

readr::write_csv(
  summary_by_cluster,
  file.path(output_dir, "climate_summary_by_cluster.csv"),
  na = ""
)
readr::write_lines(thesis_paragraph, file.path(output_dir, "climate_sentences.txt"))

cat(thesis_paragraph, "\n", sep = "")