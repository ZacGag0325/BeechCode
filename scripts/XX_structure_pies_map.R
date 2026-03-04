#!/usr/bin/env Rscript

############################################################
# scripts/XX_structure_pies_map.R
# Publication-style STRUCTURE geographic pie maps
# - all-K overview facets
# - per-K overview + zoom panels
############################################################

ensure_pkg <- function(pkg, required = TRUE) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  if (!requireNamespace(pkg, quietly = TRUE) && required) {
    stop(sprintf("MISSING FROM YOUR SIDE: Package '%s' is required but could not be loaded.", pkg))
  }
  requireNamespace(pkg, quietly = TRUE)
}

invisible(lapply(
  c("ggplot2", "dplyr", "tidyr", "readxl", "stringr", "stringi", "scales", "sf", "rnaturalearth", "scatterpie", "patchwork"),
  ensure_pkg
))

has_ggspatial <- ensure_pkg("ggspatial", required = FALSE)
has_ggrepel <- ensure_pkg("ggrepel", required = FALSE)

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
  library(scatterpie)
  library(patchwork)
})

# -----------------------------
# User options
# -----------------------------
pad_lon_over <- 2.0
pad_lat_over <- 1.5

zoom_xlim <- c(NA, NA)
zoom_ylim <- c(NA, NA)
zoom_n_sites <- 8
pad_lon_zoom <- 0.25
pad_lat_zoom <- 0.20

jitter_lon <- 0.02
jitter_lat <- 0.015

pie_r_min <- 0.04
pie_r_max <- 0.10
save_per_k_pdf <- TRUE

add_site_labels <- FALSE
add_map_decorations <- TRUE

# -----------------------------
# Strict paths
# -----------------------------
get_script_path <- function() {
  ofile <- NULL
  nfr <- sys.nframe()
  if (nfr >= 1) for (i in rev(seq_len(nfr))) {
    f <- sys.frame(i)
    if (!is.null(f$ofile)) {
      ofile <- f$ofile
      break
    }
  }
  if (is.null(ofile)) {
    args <- commandArgs(trailingOnly = FALSE)
    idx <- grep("^--file=", args)
    if (length(idx) > 0) ofile <- sub("^--file=", "", args[idx[1]])
  }
  if (is.null(ofile) || !file.exists(ofile)) return(NULL)
  normalizePath(ofile, winslash = "/", mustWork = TRUE)
}

script_path <- get_script_path()
project_root_candidates <- c(
  if (!is.null(script_path)) normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = FALSE) else NULL,
  normalizePath(".", winslash = "/", mustWork = TRUE)
)
project_root_candidates <- unique(project_root_candidates[file.exists(project_root_candidates)])

project_root <- NULL
for (cand in project_root_candidates) {
  if (dir.exists(file.path(cand, "outputs", "v1", "structure_runs", "Q_extracted")) &&
      file.exists(file.path(cand, "data", "raw", "donnees_modifiees_west_summer2024 copie.xlsx"))) {
    project_root <- cand
    break
  }
}
if (is.null(project_root)) {
  stop(
    paste0(
      "MISSING FROM YOUR SIDE: Could not resolve project root containing both required inputs.\n",
      "Expected:\n",
      "- outputs/v1/structure_runs/Q_extracted/\n",
      "- data/raw/donnees_modifiees_west_summer2024 copie.xlsx\n"
    )
  )
}

q_dir <- file.path(project_root, "outputs", "v1", "structure_runs", "Q_extracted")
meta_path <- file.path(project_root, "data", "raw", "donnees_modifiees_west_summer2024 copie.xlsx")
out_dir <- file.path(project_root, "outputs", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

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
    if (is.character(df[[j]]) || is.factor(df[[j]])) df[[j]] <- to_utf8(df[[j]])
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
  pats <- c("(?i)rep(?:licate)?[-_ ]?(\\d+)", "(?i)run[-_ ]?(\\d+)", "(?i)k\\s*\\d+[-_](\\d+)", "(?i)k\\s*\\d+\\D+(\\d+)")
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
      rnaturalearth::ne_states(country = c("canada", "united states of america"), returnclass = "sf"),
      error = function(e) NULL
    )
    if (!is.null(st) && nrow(st) > 0) return(list(data = st, method = "rnaturalearth::ne_states"))
  }
  
  st <- tryCatch(
    rnaturalearth::ne_download(scale = 50, type = "admin_1_states_provinces_lines", category = "cultural", returnclass = "sf"),
    error = function(e) NULL
  )
  if (is.null(st) || nrow(st) == 0) {
    st <- tryCatch(
      rnaturalearth::ne_download(scale = 50, type = "admin_1_states_provinces", category = "cultural", returnclass = "sf"),
      error = function(e) NULL
    )
  }
  if (is.null(st) || nrow(st) == 0) return(list(data = NULL, method = "none"))
  
  nm <- names(st)
  country_col <- intersect(c("admin", "adm0_name", "geonunit", "sr_adm0_a3", "adm0_a3"), nm)
  if (length(country_col) > 0) {
    cc <- country_col[1]
    keep <- grepl("canada|united states", tolower(as.character(st[[cc]]))) | as.character(st[[cc]]) %in% c("CAN", "USA")
    st <- st[keep, , drop = FALSE]
  }
  if (nrow(st) == 0) return(list(data = NULL, method = "none"))
  list(data = st, method = "rnaturalearth::ne_download")
}

compute_zoom_bbox <- function(site_df, zoom_xlim, zoom_ylim, n_near = 8, pad_lon = 0.25, pad_lat = 0.20) {
  manual_ok <- all(is.finite(zoom_xlim)) && all(is.finite(zoom_ylim))
  if (manual_ok) {
    return(list(
      xlim = zoom_xlim,
      ylim = zoom_ylim,
      center_site = NA_character_,
      n_sites = nrow(site_df),
      method = "manual"
    ))
  }
  
  coords <- as.matrix(site_df[, c("lon", "lat")])
  dmat <- as.matrix(stats::dist(coords, method = "euclidean"))
  dsum <- rowSums(dmat)
  center_idx <- which.min(dsum)
  center_row <- site_df[center_idx, , drop = FALSE]
  zoom_center_site <- as.character(center_row$site[[1]])
  
  ord <- order(dmat[center_idx, ])
  keep_idx <- ord[seq_len(min(n_near, nrow(site_df)))]
  near_df <- site_df[keep_idx, , drop = FALSE]
  
  xlim <- range(near_df$lon, na.rm = TRUE) + c(-pad_lon, pad_lon)
  ylim <- range(near_df$lat, na.rm = TRUE) + c(-pad_lat, pad_lat)
  
  list(
    xlim = xlim,
    ylim = ylim,
    center_site = zoom_center_site,
    n_sites = nrow(near_df),
    method = "auto_dense"
  )
}

# -----------------------------
# 1) Read metadata
# -----------------------------
sheets <- readxl::excel_sheets(meta_path)
if (length(sheets) == 0) stop("MISSING FROM YOUR SIDE: No readable sheet in metadata Excel.")

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

if (is.na(id_col) || is.na(site_col)) stop("MISSING FROM YOUR SIDE: Could not detect required ID/SITE columns.")
if (is.na(lon_col) || is.na(lat_col)) stop("MISSING FROM YOUR SIDE: Longitude and latitude columns not found.")

meta_ind <- meta_raw %>%
  transmute(
    id_ind = to_utf8(as.character(.data[[id_col]])),
    site = to_utf8(as.character(.data[[site_col]])),
    lon = suppressWarnings(as.numeric(.data[[lon_col]])),
    lat = suppressWarnings(as.numeric(.data[[lat_col]])),
    row_idx = row_number()
  )

if (nrow(meta_ind) == 0) stop("MISSING FROM YOUR SIDE: Genetic sheet has 0 data rows.")
if (any(is.na(meta_ind$site) | meta_ind$site == "")) stop("MISSING FROM YOUR SIDE: Some individuals have missing site.")
if (any(is.na(meta_ind$lon) | is.na(meta_ind$lat))) stop("MISSING FROM YOUR SIDE: Some individuals have missing coordinates.")

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
  
  k <- extract_k(path, ncol(q))
  if (k > ncol(q)) stop(sprintf("MISSING FROM YOUR SIDE: Detected K=%s but only %s numeric columns in %s.", k, ncol(q), basename(path)))
  q <- q[, seq_len(k), drop = FALSE]
  names(q) <- paste0("Q", seq_len(k))
  
  list(file = path, K = k, replicate = extract_rep(path), q = q)
}

q_files <- list.files(q_dir, full.names = TRUE)
q_files <- q_files[file.info(q_files)$isdir %in% FALSE]
if (length(q_files) == 0) stop("MISSING FROM YOUR SIDE: No files found in outputs/v1/structure_runs/Q_extracted/.")

q_list_all <- lapply(q_files, read_q_file)
q_list_all <- q_list_all[!vapply(q_list_all, is.null, logical(1))]
if (length(q_list_all) == 0) stop("MISSING FROM YOUR SIDE: No valid Q files parsed.")

k_detected_all <- sort(unique(vapply(q_list_all, function(x) x$K, integer(1))))
k_plot <- sort(intersect(k_detected_all, 2:12))
if (length(k_plot) == 0) stop(sprintf("MISSING FROM YOUR SIDE: No K in 2:12 detected. Available K: %s", paste(k_detected_all, collapse = ", ")))

q_list <- q_list_all[vapply(q_list_all, function(x) x$K %in% k_plot, logical(1))]

# -----------------------------
# 3) Join by order + aggregate
# -----------------------------
long_all <- list()
for (i in seq_along(q_list)) {
  qi <- q_list[[i]]
  q_df <- qi$q
  if (nrow(q_df) != nrow(meta_ind)) {
    stop(sprintf("MISSING FROM YOUR SIDE: Q rows (%s) do not match number of individuals in genetic sheet (%s).", nrow(q_df), nrow(meta_ind)))
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

plot_wide <- site_k %>%
  left_join(site_meta, by = "site") %>%
  select(site, K, cluster, prop, n_ind, lon, lat) %>%
  pivot_wider(names_from = cluster, values_from = prop, values_fill = 0) %>%
  filter(K %in% k_plot) %>%
  mutate(
    base_r = sqrt(n_ind),
    r = scales::rescale(base_r, to = c(pie_r_min, pie_r_max))
  )

# Stable zoom-only jitter (seeded)
set.seed(1)
site_offset <- site_meta %>%
  arrange(site) %>%
  mutate(
    idx = row_number(),
    ang = 2 * pi * (idx - 1) / pmax(n(), 1),
    dx = jitter_lon * cos(ang),
    dy = jitter_lat * sin(ang),
    lon_zoom = lon + dx,
    lat_zoom = lat + dy
  ) %>%
  select(site, lon_zoom, lat_zoom)

plot_wide <- plot_wide %>%
  left_join(site_offset, by = "site") %>%
  mutate(lon_plot = lon, lat_plot = lat)

# -----------------------------
# 4) Basemap + extents
# -----------------------------
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
states_info <- load_states_layer()
states <- states_info$data
basemap_method <- states_info$method
if (is.null(world) || nrow(world) == 0) stop("MISSING FROM YOUR SIDE: Could not load world basemap from rnaturalearth.")
if (is.null(states) || nrow(states) == 0) stop("MISSING FROM YOUR SIDE: Could not load provinces/states boundaries.")

lon_range <- range(site_meta$lon, na.rm = TRUE)
lat_range <- range(site_meta$lat, na.rm = TRUE)
over_xlim <- lon_range + c(-pad_lon_over, pad_lon_over)
over_ylim <- lat_range + c(-pad_lat_over, pad_lat_over)

zoom_info <- compute_zoom_bbox(
  site_df = site_meta,
  zoom_xlim = zoom_xlim,
  zoom_ylim = zoom_ylim,
  n_near = zoom_n_sites,
  pad_lon = pad_lon_zoom,
  pad_lat = pad_lat_zoom
)

zoom_xlim_use <- zoom_info$xlim
zoom_ylim_use <- zoom_info$ylim
zoom_sites <- site_meta %>%
  filter(lon >= zoom_xlim_use[1], lon <= zoom_xlim_use[2], lat >= zoom_ylim_use[1], lat <= zoom_ylim_use[2])

# -----------------------------
# 5) Diagnostics
# -----------------------------
cat("\n==== STRUCTURE map diagnostics ====\n")
cat("K detected:", paste(k_detected_all, collapse = ", "), "\n")
cat("K plotted (2..12 present):", paste(k_plot, collapse = ", "), "\n")
cat("Overview bbox xlim:", paste(signif(over_xlim, 6), collapse = " to "), "\n")
cat("Overview bbox ylim:", paste(signif(over_ylim, 6), collapse = " to "), "\n")
cat("Zoom bbox xlim:", paste(signif(zoom_xlim_use, 6), collapse = " to "), "\n")
cat("Zoom bbox ylim:", paste(signif(zoom_ylim_use, 6), collapse = " to "), "\n")
cat("Zoom method:", zoom_info$method, "\n")
if (!is.na(zoom_info$center_site)) cat("Zoom center site:", zoom_info$center_site, "\n")
cat("Number of sites in zoom bbox:", nrow(zoom_sites), "\n")
cat("Pie radius range:", sprintf("min=%.3f max=%.3f", min(plot_wide$r), max(plot_wide$r)), "\n")
cat("Basemap method:", basemap_method, "\n")

# -----------------------------
# 6) Plot builders
# -----------------------------
base_pal <- c(
  "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02",
  "#a6761d", "#666666", "#1f78b4", "#b2df8a", "#fb9a99", "#cab2d6"
)
max_k <- max(k_plot)
fill_vals <- setNames(base_pal[seq_len(max_k)], paste0("Q", seq_len(max_k)))

ref_n <- as.numeric(stats::quantile(site_meta$n_ind, probs = c(0, 0.5, 1), na.rm = TRUE, type = 1))
if (diff(range(ref_n, na.rm = TRUE)) == 0) ref_n <- c(ref_n[1], ref_n[1] + 1, ref_n[1] + 2)

leg_df <- data.frame(
  n = ref_n,
  lon_plot = over_xlim[2] - 0.16 * diff(over_xlim),
  lat_plot = over_ylim[1] + c(0.12, 0.23, 0.36) * diff(over_ylim)
)

build_map_base <- function() {
  ggplot() +
    geom_sf(data = world, fill = "#dfead6", color = "grey55", linewidth = 0.2) +
    geom_sf(data = states, color = "grey65", linewidth = 0.2, fill = NA) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "#dbe9f6", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(face = "bold")
    )
}

add_pies <- function(p, d, cols_k, use_zoom_coords = FALSE) {
  if (use_zoom_coords) {
    p + scatterpie::geom_scatterpie(
      data = d,
      aes(x = lon_zoom, y = lat_zoom, r = r),
      cols = cols_k,
      color = "black",
      linewidth = 0.2,
      alpha = 0.95
    )
  } else {
    p + scatterpie::geom_scatterpie(
      data = d,
      aes(x = lon, y = lat, r = r),
      cols = cols_k,
      color = "black",
      linewidth = 0.2,
      alpha = 0.95
    )
  }
}

safe_save_pdf <- function(filename, plot_obj, width, height) {
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  tryCatch({
    ggsave(filename, plot_obj, width = width, height = height, device = cairo_pdf)
  }, error = function(e) {
    ggsave(filename, plot_obj, width = width, height = height)
  })
}

safe_save_png <- function(filename, plot_obj, width, height, dpi = 400) {
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  tryCatch({
    ggsave(filename, plot_obj, width = width, height = height, dpi = dpi, type = "cairo")
  }, error = function(e) {
    ggsave(filename, plot_obj, width = width, height = height, dpi = dpi)
  })
}

# -----------------------------
# 7) Export all-K overview facets
# -----------------------------
p_all_base <- build_map_base()
p_all <- add_pies(p_all_base, plot_wide, paste0("Q", seq_len(max_k)), use_zoom_coords = FALSE) +
  scale_fill_manual(values = fill_vals, name = "Ancestry cluster") +
  geom_point(
    data = leg_df,
    aes(x = lon_plot, y = lat_plot, size = n),
    shape = 21, fill = "white", color = "black", stroke = 0.25, inherit.aes = FALSE
  ) +
  scale_size_continuous(
    name = "Samples (n)",
    breaks = ref_n,
    range = c(pie_r_min * 90, pie_r_max * 90),
    guide = guide_legend(order = 2)
  ) +
  coord_sf(xlim = over_xlim, ylim = over_ylim, expand = FALSE) +
  facet_wrap(~K, ncol = 4) +
  labs(
    title = "STRUCTURE mean Q by site (K = 2..12)",
    subtitle = "Pie size proportional to sqrt(n individuals per site)",
    x = "Longitude", y = "Latitude"
  ) +
  theme(legend.position = "right")

if (has_ggspatial && isTRUE(add_map_decorations)) {
  p_all <- p_all +
    ggspatial::annotation_scale(location = "bl") +
    ggspatial::annotation_north_arrow(location = "tr", style = ggspatial::north_arrow_fancy_orienteering)
}

out_pdf_all_over <- file.path(out_dir, "structure_map_pies_allK_overview.pdf")
out_png_all_over <- file.path(out_dir, "structure_map_pies_allK_overview.png")
# legacy names kept
out_pdf_all_legacy <- file.path(out_dir, "structure_map_pies_allK.pdf")
out_png_all_legacy <- file.path(out_dir, "structure_map_pies_allK.png")

safe_save_pdf(out_pdf_all_over, p_all, width = 14, height = 9)
safe_save_png(out_png_all_over, p_all, width = 14, height = 9, dpi = 400)
safe_save_pdf(out_pdf_all_legacy, p_all, width = 14, height = 9)
safe_save_png(out_png_all_legacy, p_all, width = 14, height = 9, dpi = 400)

# -----------------------------
# 8) Per-K outputs (legacy + overview/zoom)
# -----------------------------
for (k in k_plot) {
  pies_df <- plot_wide %>% filter(K == k)
  site_df <- site_meta
  cols_k <- paste0("Q", seq_len(k))
  
  message("[DEBUG] Objects available: ", paste(ls(), collapse = ", "))
  message("[DEBUG] head(site_df):")
  print(head(site_df))
  message("[DEBUG] head(pies_df):")
  print(head(pies_df))
  
  plot_pair <- tryCatch({
    # overview panel
    p_overview_base <- build_map_base()
    p_overview <- add_pies(p_overview_base, pies_df, cols_k, use_zoom_coords = FALSE) +
      scale_fill_manual(values = setNames(base_pal[seq_len(k)], cols_k), name = "Ancestry cluster") +
      geom_point(
        data = leg_df,
        aes(x = lon_plot, y = lat_plot, size = n),
        shape = 21, fill = "white", color = "black", stroke = 0.25, inherit.aes = FALSE
      ) +
      scale_size_continuous(
        name = "Samples (n)",
        breaks = ref_n,
        range = c(pie_r_min * 90, pie_r_max * 90),
        guide = guide_legend(order = 2)
      ) +
      coord_sf(xlim = over_xlim, ylim = over_ylim, expand = FALSE) +
      labs(
        title = paste0("STRUCTURE mean Q by site (K = ", k, ")"),
        subtitle = "Pie size proportional to sqrt(n individuals per site)",
        x = "Longitude", y = "Latitude"
      ) +
      theme(legend.position = "right")
    
    if (has_ggspatial && isTRUE(add_map_decorations)) {
      p_overview <- p_overview +
        ggspatial::annotation_scale(location = "bl") +
        ggspatial::annotation_north_arrow(location = "tr", style = ggspatial::north_arrow_fancy_orienteering)
    }
    
    # zoom panel (jittered pie positions only)
    p_zoom_base <- build_map_base()
    p_zoom <- add_pies(p_zoom_base, pies_df, cols_k, use_zoom_coords = TRUE) +
      scale_fill_manual(values = setNames(base_pal[seq_len(k)], cols_k), name = "Ancestry cluster") +
      coord_sf(xlim = zoom_xlim_use, ylim = zoom_ylim_use, expand = FALSE) +
      labs(title = paste0("Zoom (K = ", k, ")"), x = "Longitude", y = "Latitude") +
      theme(legend.position = "none")
    
    if (isTRUE(add_site_labels) && has_ggrepel) {
      labels_df <- pies_df %>% distinct(site, lon_zoom, lat_zoom)
      p_zoom <- p_zoom +
        ggrepel::geom_text_repel(
          data = labels_df,
          aes(x = lon_zoom, y = lat_zoom, label = site),
          size = 2.5, min.segment.length = 0, box.padding = 0.15, point.padding = 0.05,
          inherit.aes = FALSE
        )
    }
    
    p_k_pair <- p_overview + p_zoom + patchwork::plot_layout(ncol = 2, widths = c(1.25, 1))
    
    list(p_overview = p_overview, p_k_pair = p_k_pair)
  }, error = function(e) {
    message("[DEBUG] failed at per-K plotting block")
    stop(e)
  })
  
  p_overview <- plot_pair$p_overview
  p_k_pair <- plot_pair$p_k_pair
  
  # New requested outputs
  out_pdf_k_pair <- file.path(out_dir, paste0("structure_map_pies_K", k, "_overview_zoom.pdf"))
  out_png_k_pair <- file.path(out_dir, paste0("structure_map_pies_K", k, "_overview_zoom.png"))
  safe_save_pdf(out_pdf_k_pair, p_k_pair, width = 14, height = 7.5)
  safe_save_png(out_png_k_pair, p_k_pair, width = 14, height = 7.5, dpi = 400)
  
  # Legacy per-K PDF kept
  if (isTRUE(save_per_k_pdf)) {
    out_pdf_k_legacy <- file.path(out_dir, paste0("structure_map_pies_K", k, ".pdf"))
    safe_save_pdf(out_pdf_k_legacy, p_overview, width = 12, height = 8)
  }
}

cat("[structure_pies_map] Output written:", out_pdf_all_over, "\n")
cat("[structure_pies_map] Output written:", out_png_all_over, "\n")
cat("[structure_pies_map] Per-K overview+zoom written for K:", paste(k_plot, collapse = ", "), "\n")
if (isTRUE(save_per_k_pdf)) cat("[structure_pies_map] Legacy per-K PDFs preserved for K:", paste(k_plot, collapse = ", "), "\n")
cat("[structure_pies_map] Done.\n")