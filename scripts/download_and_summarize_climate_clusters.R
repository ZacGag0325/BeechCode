#!/usr/bin/env Rscript

# Download and summarize Environment and Climate Change Canada daily climate data
# for the two climate-station clusters used in the American beech study area.

required_packages <- c(
  "weathercan",
  "dplyr",
  "readr",
  "stringr",
  "lubridate",
  "tibble"
)

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}

invisible(lapply(required_packages, library, character.only = TRUE))

clusters <- tibble::tribble(
  ~cluster, ~station_label, ~climate_id, ~raw_file,
  "northern", "Sainte-Anne-du-Lac", "7036855", "data/raw/climate/st_anne_du_lac_daily.csv",
  "southern", "Notre-Dame-de-la-Paix", "7035666", "data/raw/climate/notre_dame_de_la_paix_daily.csv"
)

preferred_start <- lubridate::ymd("1981-01-01")
preferred_end <- lubridate::ymd("2010-12-31")

raw_dir <- "data/raw/climate"
output_dir <- "outputs/climate"
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

message("Searching the weathercan station inventory for target climate IDs...")
station_inventory <- weathercan::stations()

get_daily_station <- function(climate_id) {
  matches <- station_inventory |>
    dplyr::mutate(climate_id = as.character(.data$climate_id)) |>
    dplyr::filter(.data$climate_id == as.character(climate_id), .data$interval == "day") |>
    dplyr::arrange(.data$start)
  
  if (nrow(matches) == 0) {
    stop("No daily station_id found for climate ID ", climate_id, call. = FALSE)
  }
  
  matches |>
    dplyr::mutate(
      covers_preferred = .data$start <= preferred_start & .data$end >= preferred_end,
      record_length_days = as.numeric(.data$end - .data$start)
    ) |>
    dplyr::arrange(dplyr::desc(.data$covers_preferred), dplyr::desc(.data$record_length_days)) |>
    dplyr::slice(1)
}

daily_stations <- dplyr::bind_rows(lapply(seq_len(nrow(clusters)), function(i) {
  cluster_row <- clusters[i, ]
  station_info <- get_daily_station(cluster_row$climate_id)
  
  dplyr::bind_cols(
    cluster_row,
    station_info |>
      dplyr::rename_with(~ paste0("station_info_", .x))
  )
}))

download_cluster_daily <- function(cluster_row) {
  station_start <- lubridate::as_date(cluster_row$station_info_start)
  station_end <- lubridate::as_date(cluster_row$station_info_end)
  
  use_start <- max(preferred_start, station_start, na.rm = TRUE)
  use_end <- min(preferred_end, station_end, na.rm = TRUE)
  
  if (use_start > use_end || station_start > preferred_start || station_end < preferred_end) {
    use_start <- station_start
    use_end <- station_end
    message(
      stringr::str_glue(
        "Climate ID {cluster_row$climate_id} does not fully cover 1981-2010; ",
        "using full available daily record ({lubridate::year(use_start)}-{lubridate::year(use_end)})."
      )
    )
  } else {
    message(
      stringr::str_glue(
        "Climate ID {cluster_row$climate_id} covers 1981-2010; using the 1981-2010 normal period."
      )
    )
  }
  
  daily_data <- weathercan::weather_dl(
    station_ids = cluster_row$station_info_station_id,
    start = use_start,
    end = use_end,
    interval = "day"
  ) |>
    dplyr::mutate(
      cluster = cluster_row$cluster,
      cluster_station = cluster_row$station_label,
      target_climate_id = cluster_row$climate_id,
      station_id = cluster_row$station_info_station_id,
      year = lubridate::year(.data$date)
    )
  
  readr::write_csv(daily_data, cluster_row$raw_file, na = "")
  
  actual_years <- range(daily_data$year, na.rm = TRUE)
  message(
    stringr::str_glue(
      "Saved {nrow(daily_data)} daily records for the {cluster_row$cluster} cluster ",
      "to {cluster_row$raw_file}; actual years used: {actual_years[1]}-{actual_years[2]}."
    )
  )
  
  daily_data
}

all_daily <- dplyr::bind_rows(
  lapply(seq_len(nrow(daily_stations)), function(i) download_cluster_daily(daily_stations[i, ]))
)

annual_by_cluster <- all_daily |>
  dplyr::group_by(.data$cluster, .data$cluster_station, .data$target_climate_id, .data$station_id, .data$year) |>
  dplyr::summarise(
    mean_annual_temperature_c = mean(.data$mean_temp, na.rm = TRUE),
    annual_precipitation_mm = sum(.data$total_precip, na.rm = TRUE),
    annual_snowfall_cm = sum(.data$total_snow, na.rm = TRUE),
    degree_days_above_0_c = sum(pmax(.data$mean_temp, 0), na.rm = TRUE),
    daily_records = dplyr::n(),
    .groups = "drop"
  ) |>
  dplyr::arrange(.data$cluster, .data$year)

summary_by_cluster <- annual_by_cluster |>
  dplyr::group_by(.data$cluster, .data$cluster_station, .data$target_climate_id, .data$station_id) |>
  dplyr::summarise(
    mean_annual_temperature_c = mean(.data$mean_annual_temperature_c, na.rm = TRUE),
    mean_annual_precipitation_mm = mean(.data$annual_precipitation_mm, na.rm = TRUE),
    mean_annual_snowfall_cm = mean(.data$annual_snowfall_cm, na.rm = TRUE),
    mean_degree_days_above_0_c = mean(.data$degree_days_above_0_c, na.rm = TRUE),
    first_year_used = min(.data$year, na.rm = TRUE),
    last_year_used = max(.data$year, na.rm = TRUE),
    number_of_years_used = dplyr::n_distinct(.data$year),
    .groups = "drop"
  ) |>
  dplyr::arrange(.data$cluster)

sentences <- summary_by_cluster |>
  dplyr::mutate(
    sentence = stringr::str_glue(
      "For the {cluster} cluster, represented by the {cluster_station} station, ",
      "mean annual temperature was {round(mean_annual_temperature_c, 1)}°C, ",
      "mean annual precipitation was {round(mean_annual_precipitation_mm, 0)} mm, ",
      "mean annual snowfall was {round(mean_annual_snowfall_cm, 0)} cm, ",
      "and the average number of degree-days above 0°C was {round(mean_degree_days_above_0_c, 0)} ",
      "between {first_year_used} and {last_year_used}."
    )
  ) |>
  dplyr::pull(.data$sentence)

readr::write_csv(annual_by_cluster, file.path(output_dir, "climate_annual_by_cluster.csv"), na = "")
readr::write_csv(summary_by_cluster, file.path(output_dir, "climate_summary_by_cluster.csv"), na = "")
readr::write_lines(sentences, file.path(output_dir, "climate_sentences.txt"))

cat(paste(sentences, collapse = "\n"), "\n", sep = "")