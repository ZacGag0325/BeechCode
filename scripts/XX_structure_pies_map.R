#!/usr/bin/env Rscript

############################################################
# scripts/XX_structure_pies_map.R
# Publication-style STRUCTURE geographic pie maps (ALL K)
############################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readxl)
  library(stringr)
  library(stringi)
  library(scales)
  library(sf)
  library(rnaturalearth)
})

if (!requireNamespace("scatterpie", quietly = TRUE)) {
  stop("MISSING FROM YOUR SIDE: Package 'scatterpie' is required. Install it with install.packages('scatterpie').")
}

has_ggspatial <- requireNamespace("ggspatial", quietly = TRUE)
has_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)

# -----------------------------
# User options
# -----------------------------
pad_lon <- 2.0
pad_lat <- 1.5
pie_r_min <- 0.04
pie_r_max <- 0.10
save_per_k_pdf <- TRUE

# overlap reduction (stable deterministic offsets)
jitter_lon <- 0.03
jitter_lat <- 0.02

# optional labels
add_site_labels <- FALSE

# optional tile background
use_tiles <- FALSE

# -----------------------------
# Strict paths
# -----------------------------
project_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
q_dir <- file.path(project_root, "outputs", "v1", "structure_runs", "Q_extracted")
meta_path <- file.path(project_root, "data", "raw", "donnees_modifiees_west_summer2024 copie.xlsx")
out_dir <- file.path(project_root, "outputs", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(q_dir)) {
  stop("MISSING FROM YOUR SIDE: Q directory not found at outputs/v1/structure_runs/Q_extracted/.")
}
if (!file.exists(meta_path)) {
  stop("MISSING FROM YOUR SIDE: Missing metadata file data/raw/donnees_modifiees_west_summer2024 copie.xlsx.")
}

# -----------------------------
# Helpers
# -----------------------------
to_utf8 <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  if (!is.character(x)) return(x)
  y <- iconv(x, from = "", to = "UTF-8", sub = "")
  y[is.na(y)] <- ""
  y
}

clean_names <- function(x) {
  x <- to_utf8(x)
  x <- stringi::stri_trans_general(x, "Latin-ASCII")
  x <- tolower(x)
  x <- gsub("\\s+", "_", x)
  x <- gsub("[^a-z0-9_]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

apply_utf8_df <- function(df) {
  names(df) <- to_utf8(names(df))
  for (j in seq_along(df)) {
    if (is.character(df[[j]]) || is.factor(df[[j]])) {
      df[[j]] <- to_utf8(df[[j]])
    }
  }
  df
}

remove_empty_auto_cols <- function(df) {
  if (ncol(df) == 0) return(df)
  nm <- names(df)
  nm[is.na(nm)] <- ""
  keep_nm <- !grepl("^\\.\\.\\.[0-9]+$", nm)
  df <- df[, keep_nm, drop = FALSE]
  if (ncol(df) == 0) return(df)
  
  non_empty <- vapply(df, function(col) {
    z <- trimws(as.character(col))
    any(!is.na(z) & z != "")
  }, logical(1))
  df[, non_empty, drop = FALSE]
}

pick_first <- function(nms, candidates) {
  hit <- intersect(candidates, nms)
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

extract_k <- function(path, default_k) {
  b <- basename(path)
  m <- stringr::str_match(b, "(?i)k\\s*[-_=]?\\s*(\\d+)")
  if (!is.na(m[1, 2])) return(as.integer(m[1, 2]))
  as.integer(default_k)
}

extract_rep <- function(path) {
  b <- basename(path)
  pats <- c(
    "(?i)rep(?:licate)?[-_ ]?(\\d+)",
    "(?i)run[-_ ]?(\\d+)",
    "(?i)k\\s*\\d+[-_](\\d+)",
    "(?i)k\\s*\\d+\\D+(\\d+)"
  )
  for (pt in pats) {
    m <- stringr::str_match(b, pt)
    if (!is.na(m[1, 2])) return(as.integer(m[1, 2]))
  }
  1L
}

is_num_like <- function(x) {
  suppressWarnings(y <- as.numeric(as.character(x)))
  mean(!is.na(y)) >= 0.95
}

load_states_layer <- function() {
  if (requireNamespace("rnaturalearthhires", quietly = TRUE)) {
    st <- tryCatch(
      rnaturalearth::ne_states(
        country = c("canada", "united states of america"),
        returnclass = "sf"
      ),
      error = function(e) NULL
    )
    if (!is.null(st) && nrow(st) > 0) return(list(data = st, method = "rnaturalearth::ne_states"))
  }
  
  st <- tryCatch(
    rnaturalearth::ne_download(
      scale = 50,
      type = "admin_1_states_provinces_lines",
      category = "cultural",
      returnclass = "sf"
    ),
    error = function(e) NULL
  )
  
  if (is.null(st) || nrow(st) == 0) {
    st <- tryCatch(
      rnaturalearth::ne_download(
        scale = 50,
        type = "admin_1_states_provinces",
        category = "cultural",
        returnclass = "sf"
      ),
      error = function(e) NULL
    )
  }
  
  if (is.null(st) || nrow(st) == 0) return(list(data = NULL, method = "none"))
  
  nm <- names(st)
  country_col <- intersect(c("admin", "adm0_name", "geonunit", "sr_adm0_a3", "adm0_a3"), nm)
  if (length(country_col) > 0) {
    cc <- country_col[1]
    keep <- grepl("canada|united states", tolower(as.character(st[[cc]]))) |
      as.character(st[[cc]]) %in% c("CAN", "USA")
    st <- st[keep, , drop = FALSE]
  }
  
  if (nrow(st) == 0) return(list(data = NULL, method = "none"))
  list(data = st, method = "rnaturalearth::ne_download")
}

# -----------------------------
# 1) Read metadata
# -----------------------------
sheets <- readxl::excel_sheets(meta_path)
if (length(sheets) == 0) {
  stop("MISSING FROM YOUR SIDE: No readable sheet in metadata Excel.")
}

target_sheet <- NA_character_
meta_raw <- NULL
for (sh in sheets) {
  tmp <- tryCatch(readxl::read_excel(meta_path, sheet = sh), error = function(e) NULL)
  if (is.null(tmp) || nrow(tmp) == 0 || ncol(tmp) == 0) next
  tmp <- apply_utf8_df(tmp)
  names(tmp) <- clean_names(names(tmp))
  tmp <- remove_empty_auto_cols(tmp)
  if (ncol(tmp) == 0) next
  
  if ("nom_labo_echantillons" %in% names(tmp) || "nom_labo_echantillon" %in% names(tmp)) {
    target_sheet <- sh
    meta_raw <- tmp
    break
  }
}

if (is.na(target_sheet) || is.null(meta_raw)) {
  stop("MISSING FROM YOUR SIDE: Could not find a sheet containing Nom_Labo_Échantillons in metadata Excel.")
}

id_col <- pick_first(names(meta_raw), c("nom_labo_echantillons", "nom_labo_echantillon"))
site_col <- pick_first(names(meta_raw), c("numero_population", "num_population"))
lon_col <- pick_first(names(meta_raw), c("lon", "long", "longitude", "x"))
lat_col <- pick_first(names(meta_raw), c("lat", "latitude", "y"))

if (is.na(id_col) || is.na(site_col)) {
  stop(sprintf(
    "MISSING FROM YOUR SIDE: Could not detect required ID/SITE columns. Columns found: %s",
    paste(names(meta_raw), collapse = ", ")
  ))
}
if (is.na(lon_col) || is.na(lat_col)) {
  stop("MISSING FROM YOUR SIDE: Longitude and latitude columns not found in genetic sheet.")
}

meta_ind <- meta_raw %>%
  transmute(
    id_ind = to_utf8(as.character(.data[[id_col]])),
    site = to_utf8(as.character(.data[[site_col]])),
    lon = suppressWarnings(as.numeric(.data[[lon_col]])),
    lat = suppressWarnings(as.numeric(.data[[lat_col]])),
    row_idx = row_number()
  )

if (nrow(meta_ind) == 0) {
  stop("MISSING FROM YOUR SIDE: Genetic sheet has 0 data rows.")
}
if (any(is.na(meta_ind$site) | meta_ind$site == "")) {
  stop("MISSING FROM YOUR SIDE: Some individuals have missing Numéro_Population in genetic sheet.")
}
if (any(is.na(meta_ind$lon) | is.na(meta_ind$lat))) {
  stop("MISSING FROM YOUR SIDE: Some individuals have missing longitude/latitude values in genetic sheet.")
}
if (any(abs(meta_ind$lon) > 180, na.rm = TRUE) || any(abs(meta_ind$lat) > 90, na.rm = TRUE)) {
  stop(
    paste0(
      "MISSING FROM YOUR SIDE: Coordinates are not valid lon/lat degrees (WGS84). ",
      "Values suggest projected coordinates (e.g., UTM)."
    )
  )
}

# -----------------------------
# 2) Read Q files
# -----------------------------
read_q_file <- function(path) {
  dat <- tryCatch(
    utils::read.table(path, header = FALSE, fill = TRUE, stringsAsFactors = FALSE,
                      check.names = FALSE, comment.char = "", sep = ""),
    error = function(e) NULL
  )
  if (is.null(dat) || nrow(dat) == 0 || ncol(dat) == 0) return(NULL)
  
  dat <- apply_utf8_df(dat)
  names(dat) <- clean_names(names(dat))
  dat <- remove_empty_auto_cols(dat)
  if (ncol(dat) < 2) return(NULL)
  
  num_like <- vapply(dat, is_num_like, logical(1))
  if (sum(num_like) < 2) return(NULL)
  
  q <- dat[, num_like, drop = FALSE]
  q <- as.data.frame(lapply(q, function(z) suppressWarnings(as.numeric(as.character(z)))))
  
  q <- q[stats::complete.cases(q), , drop = FALSE]
  if (nrow(q) == 0) return(NULL)
  
  if (any(as.matrix(q) < -1e-8 | as.matrix(q) > 1 + 1e-8, na.rm = TRUE)) {
    stop(sprintf("MISSING FROM YOUR SIDE: Q values outside [0,1] in file %s.", basename(path)))
  }
  
  rs <- rowSums(q)
  if (mean(abs(rs - 1) <= 0.02) < 0.95) {
    stop(sprintf("MISSING FROM YOUR SIDE: Q rowSums are not approximately 1 in file %s.", basename(path)))
  }
  
  k <- extract_k(path, ncol(q))
  if (k > ncol(q)) {
    stop(sprintf("MISSING FROM YOUR SIDE: Detected K=%s but only %s numeric columns in %s.",
                 k, ncol(q), basename(path)))
  }
  
  q <- q[, seq_len(k), drop = FALSE]
  names(q) <- paste0("Q", seq_len(k))
  
  list(file = path, K = k, replicate = extract_rep(path), q = q)
}

q_files <- list.files(q_dir, full.names = TRUE)
q_files <- q_files[file.info(q_files)$isdir %in% FALSE]
if (length(q_files) == 0) {
  stop("MISSING FROM YOUR SIDE: No files found in outputs/v1/structure_runs/Q_extracted/.")
}

q_list_all <- lapply(q_files, read_q_file)
q_list_all <- q_list_all[!vapply(q_list_all, is.null, logical(1))]
if (length(q_list_all) == 0) {
  stop("MISSING FROM YOUR SIDE: No valid Q files could be parsed from outputs/v1/structure_runs/Q_extracted/.")
}

k_detected_all <- sort(unique(vapply(q_list_all, function(x) x$K, integer(1))))
k_plot <- sort(intersect(k_detected_all, 2:12))
if (length(k_plot) == 0) {
  stop(sprintf(
    "MISSING FROM YOUR SIDE: No K in 2:12 detected in Q files. Available K: %s",
    paste(k_detected_all, collapse = ", ")
  ))
}

q_list <- q_list_all[vapply(q_list_all, function(x) x$K %in% k_plot, logical(1))]

# -----------------------------
# 3) Join by order + aggregate
# -----------------------------
long_all <- list()
for (i in seq_along(q_list)) {
  qi <- q_list[[i]]
  q_df <- qi$q
  
  if (nrow(q_df) != nrow(meta_ind)) {
    stop(sprintf(
      paste0(
        "MISSING FROM YOUR SIDE: Q rows (%s) do not match number of individuals in genetic sheet (%s). ",
        "Ensure STRUCTURE used the same individual order."
      ),
      nrow(q_df), nrow(meta_ind)
    ))
  }
  
  merged <- dplyr::bind_cols(meta_ind, q_df)
  long_i <- merged %>%
    pivot_longer(cols = starts_with("Q"), names_to = "cluster", values_to = "prop") %>%
    mutate(K = qi$K, replicate = qi$replicate)
  
  long_all[[i]] <- long_i
}
long_all <- bind_rows(long_all)

site_rep <- long_all %>%
  group_by(site, K, replicate, cluster) %>%
  summarise(prop = mean(prop, na.rm = TRUE), .groups = "drop")

site_k <- site_rep %>%
  group_by(site, K, cluster) %>%
  summarise(prop = mean(prop, na.rm = TRUE), .groups = "drop")

site_meta <- meta_ind %>%
  group_by(site) %>%
  summarise(
    n_ind = n(),
    lon = mean(lon, na.rm = TRUE),
    lat = mean(lat, na.rm = TRUE),
    .groups = "drop"
  )

plot_df <- site_k %>% left_join(site_meta, by = "site")
plot_wide <- plot_df %>%
  select(site, K, cluster, prop, n_ind, lon, lat) %>%
  pivot_wider(names_from = cluster, values_from = prop, values_fill = 0) %>%
  filter(K %in% k_plot)

# bounded radii
plot_wide <- plot_wide %>%
  mutate(
    base_r = sqrt(n_ind),
    r = scales::rescale(base_r, to = c(pie_r_min, pie_r_max))
  )

# deterministic offsets to reduce overlap
site_offset <- site_meta %>%
  arrange(site) %>%
  mutate(
    idx = row_number(),
    ang = 2 * pi * (idx - 1) / pmax(n(), 1),
    dx = jitter_lon * cos(ang),
    dy = jitter_lat * sin(ang),
    lon_plot = lon + dx,
    lat_plot = lat + dy
  ) %>%
  select(site, lon_plot, lat_plot)

plot_wide <- plot_wide %>% left_join(site_offset, by = "site")

# -----------------------------
# 4) Basemap + extent
# -----------------------------
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
states_info <- load_states_layer()
states <- states_info$data
basemap_method <- states_info$method

if (is.null(world) || nrow(world) == 0) {
  stop("MISSING FROM YOUR SIDE: Could not load world basemap from rnaturalearth.")
}
if (is.null(states) || nrow(states) == 0) {
  stop(
    paste0(
      "MISSING FROM YOUR SIDE: Could not load provinces/states boundaries. ",
      "Install 'rnaturalearthhires' OR allow Natural Earth downloads, then rerun."
    )
  )
}

lon_range <- range(site_meta$lon, na.rm = TRUE)
lat_range <- range(site_meta$lat, na.rm = TRUE)
xlim <- lon_range + c(-pad_lon, pad_lon)
ylim <- lat_range + c(-pad_lat, pad_lat)

# -----------------------------
# 5) Diagnostics
# -----------------------------
rep_by_k <- dplyr::bind_rows(lapply(q_list_all, function(x) {
  data.frame(K = x$K, replicate = x$replicate)
})) %>%
  count(K, name = "n_runs")

cat("\n==== STRUCTURE map diagnostics ====\n")
cat("Individuals:", nrow(meta_ind), "\n")
cat("Sites:", dplyr::n_distinct(meta_ind$site), "\n")
cat("K detected:", paste(k_detected_all, collapse = ", "), "\n")
cat("K plotted (2..12 present):", paste(k_plot, collapse = ", "), "\n")
cat("Replicates detected:", ifelse(any(rep_by_k$n_runs > 1), "YES", "NO"), "\n")
cat("Coordinates found in genetic sheet: YES\n")
cat("Sheet used:", target_sheet, "\n")
cat("Lon range:", paste(signif(lon_range, 6), collapse = " to "), "\n")
cat("Lat range:", paste(signif(lat_range, 6), collapse = " to "), "\n")
cat("xlim used:", paste(signif(xlim, 6), collapse = " to "), "\n")
cat("ylim used:", paste(signif(ylim, 6), collapse = " to "), "\n")
cat("Jitter settings (lon, lat):", paste0(jitter_lon, ", ", jitter_lat), "\n")
cat("Pie radius range:", sprintf("min=%.3f max=%.3f", min(plot_wide$r), max(plot_wide$r)), "\n")
cat("Basemap method:", basemap_method, "\n")

# -----------------------------
# 6) Plot objects
# -----------------------------
base_pal <- c(
  "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e",
  "#e6ab02", "#a6761d", "#666666", "#1f78b4", "#b2df8a",
  "#fb9a99", "#cab2d6"
)
max_k <- max(k_plot)
fill_vals <- setNames(base_pal[seq_len(max_k)], paste0("Q", seq_len(max_k)))

ref_n <- as.numeric(stats::quantile(site_meta$n_ind, probs = c(0, 0.5, 1), na.rm = TRUE, type = 1))
if (any(!is.finite(ref_n))) {
  ref_n <- c(min(site_meta$n_ind, na.rm = TRUE), stats::median(site_meta$n_ind, na.rm = TRUE), max(site_meta$n_ind, na.rm = TRUE))
}
if (diff(range(ref_n, na.rm = TRUE)) == 0) {
  ref_n <- c(ref_n[1], ref_n[1] + 1, ref_n[1] + 2)
}

x0 <- xlim[2] - 0.22 * diff(xlim)
y0 <- ylim[1] + 0.16 * diff(ylim)
leg_df <- data.frame(
  n = ref_n,
  lon_plot = rep(x0, 3),
  lat_plot = y0 + c(0, 0.25, 0.5)
)

p_all <- ggplot()

if (isTRUE(use_tiles) && has_ggspatial) {
  cache_dir <- file.path(project_root, "outputs", "cache_tiles")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  tile_added <- FALSE
  p_try <- tryCatch({
    p_all + ggspatial::annotation_map_tile(type = "osm", cachedir = cache_dir)
  }, error = function(e) NULL)
  if (!is.null(p_try)) {
    p_all <- p_try
    tile_added <- TRUE
    basemap_method <- "ggspatial::annotation_map_tile(osm)"
  }
  if (!tile_added) {
    p_all <- ggplot()
  }
}

# Colored basemap fallback / default
if (!grepl("annotation_map_tile", basemap_method)) {
  p_all <- p_all +
    geom_sf(data = world, fill = "#dbe9f6", color = NA) +
    geom_sf(data = states, fill = "#dfead6", color = "grey65", linewidth = 0.2)
}

p_all <- p_all +
  scatterpie::geom_scatterpie(
    data = plot_wide,
    aes(x = lon_plot, y = lat_plot, r = r),
    cols = paste0("Q", seq_len(max_k)),
    color = "black",
    linewidth = 0.2,
    alpha = 0.95
  ) +
  scale_fill_manual(values = fill_vals, name = "Ancestry cluster") +
  geom_point(
    data = leg_df,
    aes(x = lon_plot, y = lat_plot, size = n),
    shape = 21,
    fill = "white",
    color = "black",
    stroke = 0.25,
    inherit.aes = FALSE
  ) +
  scale_size_continuous(
    name = "Samples (n)",
    breaks = ref_n,
    range = c(pie_r_min * 90, pie_r_max * 90),
    guide = guide_legend(order = 2)
  ) +
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
  facet_wrap(~K, ncol = 4) +
  labs(
    title = "STRUCTURE mean Q by site (K = 2..12)",
    subtitle = "Pie size proportional to sqrt(n individuals per site)",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "#dbe9f6", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    legend.position = "right",
    plot.title = element_text(face = "bold")
  )

if (isTRUE(add_site_labels) && has_ggrepel) {
  labels_df <- plot_wide %>%
    group_by(site) %>%
    summarise(lon_plot = first(lon_plot), lat_plot = first(lat_plot), .groups = "drop")
  p_all <- p_all +
    ggrepel::geom_text_repel(
      data = labels_df,
      aes(x = lon_plot, y = lat_plot, label = site),
      size = 2.5,
      min.segment.length = 0,
      box.padding = 0.15,
      point.padding = 0.05,
      inherit.aes = FALSE
    )
}

if (has_ggspatial) {
  p_all <- p_all +
    ggspatial::annotation_scale(location = "bl") +
    ggspatial::annotation_north_arrow(
      location = "tr",
      style = ggspatial::north_arrow_fancy_orienteering
    )
}

# -----------------------------
# 7) Export ALL-K
# -----------------------------
out_pdf_all <- file.path(out_dir, "structure_map_pies_allK.pdf")
out_png_all <- file.path(out_dir, "structure_map_pies_allK.png")

ggplot2::ggsave(out_pdf_all, p_all, width = 14, height = 9)
ggplot2::ggsave(out_png_all, p_all, width = 14, height = 9, dpi = 400)

if (!file.exists(out_pdf_all)) {
  stop("MISSING FROM YOUR SIDE: Failed to produce combined ALL-K PDF output.")
}

# -----------------------------
# 8) Optional per-K PDFs
# -----------------------------
if (isTRUE(save_per_k_pdf)) {
  for (k in k_plot) {
    d_k <- plot_wide %>% filter(K == k)
    cols_k <- paste0("Q", seq_len(k))
    
    p_k <- ggplot()
    if (!grepl("annotation_map_tile", basemap_method)) {
      p_k <- p_k +
        geom_sf(data = world, fill = "#dbe9f6", color = NA) +
        geom_sf(data = states, fill = "#dfead6", color = "grey65", linewidth = 0.2)
    }
    
    p_k <- p_k +
      scatterpie::geom_scatterpie(
        data = d_k,
        aes(x = lon_plot, y = lat_plot, r = r),
        cols = cols_k,
        color = "black",
        linewidth = 0.2,
        alpha = 0.95
      ) +
      scale_fill_manual(
        values = setNames(base_pal[seq_len(k)], cols_k),
        name = "Ancestry cluster"
      ) +
      geom_point(
        data = leg_df,
        aes(x = lon_plot, y = lat_plot, size = n),
        shape = 21,
        fill = "white",
        color = "black",
        stroke = 0.25,
        inherit.aes = FALSE
      ) +
      scale_size_continuous(
        name = "Samples (n)",
        breaks = ref_n,
        range = c(pie_r_min * 90, pie_r_max * 90),
        guide = guide_legend(order = 2)
      ) +
      coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
      labs(
        title = paste0("STRUCTURE mean Q by site (K = ", k, ")"),
        subtitle = "Pie size proportional to sqrt(n individuals per site)",
        x = "Longitude",
        y = "Latitude"
      ) +
      theme_bw(base_size = 12) +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "#dbe9f6", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        legend.position = "right",
        plot.title = element_text(face = "bold")
      )
    
    if (isTRUE(add_site_labels) && has_ggrepel) {
      labels_df <- d_k %>% distinct(site, lon_plot, lat_plot)
      p_k <- p_k +
        ggrepel::geom_text_repel(
          data = labels_df,
          aes(x = lon_plot, y = lat_plot, label = site),
          size = 2.5,
          min.segment.length = 0,
          box.padding = 0.15,
          point.padding = 0.05,
          inherit.aes = FALSE
        )
    }
    
    if (has_ggspatial) {
      p_k <- p_k +
        ggspatial::annotation_scale(location = "bl") +
        ggspatial::annotation_north_arrow(
          location = "tr",
          style = ggspatial::north_arrow_fancy_orienteering
        )
    }
    
    out_pdf_k <- file.path(out_dir, paste0("structure_map_pies_K", k, ".pdf"))
    ggplot2::ggsave(out_pdf_k, p_k, width = 12, height = 8)
  }
}

cat("[structure_pies_map] Output written:", out_pdf_all, "\n")
cat("[structure_pies_map] Output written:", out_png_all, "\n")
if (isTRUE(save_per_k_pdf)) {
  cat("[structure_pies_map] Per-K PDFs written for K:", paste(k_plot, collapse = ", "), "\n")
}
cat("[structure_pies_map] Done.\n")