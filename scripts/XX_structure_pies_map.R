############################################################
# scripts/XX_structure_pies_map.R
# STRUCTURE pies by site (map + non-spatial fallback)
############################################################

options(repos = c(CRAN = "https://cloud.r-project.org"))

suppressPackageStartupMessages(library(here))

FALLBACK_ROOT <- file.path(path.expand("~"), "Desktop", "BeechCode")
candidate_root <- here::here()
if (dir.exists(file.path(candidate_root, "data"))) {
  PROJECT_ROOT <- candidate_root
} else if (dir.exists(file.path(FALLBACK_ROOT, "data"))) {
  PROJECT_ROOT <- FALLBACK_ROOT
  setwd(PROJECT_ROOT)
} else {
  stop("Can't find project root. Open BeechCode.Rproj or set FALLBACK_ROOT.")
}

OUT_DIR <- file.path(PROJECT_ROOT, "outputs", "figures")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
cat("\n[structure_pies_map] Project root:", PROJECT_ROOT, "\n")

req <- c("dplyr", "tidyr", "ggplot2", "readr", "stringi", "scales")
for (p in req) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(stringi)
})

use_scatterpie <- requireNamespace("scatterpie", quietly = TRUE)
if (!use_scatterpie) {
  cat("[structure_pies_map] 'scatterpie' missing, trying install...\n")
  try(install.packages("scatterpie"), silent = TRUE)
  use_scatterpie <- requireNamespace("scatterpie", quietly = TRUE)
}
use_sf <- requireNamespace("sf", quietly = TRUE)
use_rnaturalearth <- requireNamespace("rnaturalearth", quietly = TRUE)
use_ggspatial <- requireNamespace("ggspatial", quietly = TRUE)

q_dirs <- c("outputs/v1/structure_runs/Q_extracted")

clean_names <- function(x) {
  x <- iconv(as.character(x), from = "", to = "UTF-8", sub = "")
  x[is.na(x)] <- ""
  x <- trimws(x)
  x <- stringi::stri_trans_general(x, "Latin-ASCII")
  x <- gsub("\\s+", "_", x)
  x <- gsub("[^a-zA-Z0-9_]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  tolower(x)
}

safe_basename <- function(x) {
  b <- basename(x)
  b <- iconv(as.character(b), from = "", to = "UTF-8", sub = "")
  b[is.na(b)] <- ""
  b
}

read_tabular_clean <- function(path) {
  ext <- tolower(tools::file_ext(path))
  dat <- tryCatch({
    if (ext == "csv") {
      suppressWarnings(readr::read_csv(path, show_col_types = FALSE, guess_max = 5000, locale = readr::locale(encoding = "UTF-8"), comment = "#", name_repair = "minimal"))
    } else if (ext %in% c("tsv", "txt", "dat", "pop")) {
      utils::read.table(path, header = FALSE, fill = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = "", comment.char = "")
    } else {
      suppressWarnings(readr::read_csv(path, show_col_types = FALSE, guess_max = 5000, locale = readr::locale(encoding = "UTF-8"), comment = "#", name_repair = "minimal"))
    }
  }, error = function(e) NULL)
  if (is.null(dat)) return(NULL)
  
  dat <- as.data.frame(dat, stringsAsFactors = FALSE, check.names = FALSE)
  nms <- names(dat)
  if (is.null(nms)) nms <- paste0("v", seq_len(ncol(dat)))
  names(dat) <- clean_names(nms)
  
  keep <- !grepl("^\\.\\.\\.[0-9]+$", names(dat)) & !grepl("^\\.\\.[0-9]+$", names(dat))
  dat <- dat[, keep, drop = FALSE]
  non_empty <- vapply(dat, function(col) {
    v <- trimws(as.character(col)); any(!is.na(v) & v != "")
  }, logical(1))
  dat <- dat[, non_empty, drop = FALSE]
  names(dat) <- make.unique(names(dat), sep = "_")
  dat
}

choose_col <- function(nms, candidates) {
  i <- which(nms %in% candidates)
  if (length(i) > 0) nms[i[1]] else NA_character_
}

extract_k <- function(path) {
  nm <- safe_basename(path)
  pats <- c("(?i)k\\s*=\\s*(\\d+)", "(?i)-k(\\d+)\\b", "(?i)_k(\\d+)\\b", "(?i)\\bk(\\d+)\\b")
  for (pt in pats) {
    m <- regexec(pt, nm, perl = TRUE)
    g <- regmatches(nm, m)[[1]]
    if (length(g) >= 2) return(as.integer(g[2]))
  }
  NA_integer_
}

extract_rep <- function(path) {
  nm <- safe_basename(path)
  pats <- c(
    "(?i)-k\\d+-(\\d+)_f\\.q$",
    "(?i)-k\\d+[-_](\\d+)\\b",
    "(?i)rep(?:licate)?[-_ ]?(\\d+)\\b",
    "(?i)run[-_ ]?(\\d+)\\b"
  )
  for (pt in pats) {
    m <- regexec(pt, nm, perl = TRUE)
    g <- regmatches(nm, m)[[1]]
    if (length(g) >= 2) return(as.integer(g[2]))
  }
  1L
}

is_numeric_like <- function(x) {
  suppressWarnings(y <- as.numeric(as.character(x)))
  mean(!is.na(y)) >= 0.95
}

site_candidates <- c("site", "population", "pop", "numero_population", "num_population", "station", "sampling_site", "localite", "location", "site_id", "pop_id", "id_site", "v1")
id_candidates <- c("id", "ind_id", "ind", "nom_labo_echantillons", "sample", "sample_id", "individual", "individu", "v2")
lon_candidates <- c("lon", "long", "longitude", "x")
lat_candidates <- c("lat", "latitude", "y")

meta_site_path <- file.path(PROJECT_ROOT, "inputs", "site_metadata.csv")
if (!file.exists(meta_site_path)) stop("Missing inputs/site_metadata.csv")
meta_site <- read_tabular_clean(meta_site_path)
if (is.null(meta_site) || nrow(meta_site) == 0) stop("Could not read inputs/site_metadata.csv")
cat("[structure_pies_map] meta_site columns:\n")
print(names(meta_site))

site_col_site <- choose_col(names(meta_site), site_candidates)
if (is.na(site_col_site)) stop("inputs/site_metadata.csv must contain a site/population column.")
meta_site <- meta_site %>% mutate(site = as.character(.data[[site_col_site]]))

coords_paths <- c(
  file.path(PROJECT_ROOT, "inputs", "site_coordinates.csv"),
  file.path(PROJECT_ROOT, "inputs", "site_coords.csv"),
  file.path(PROJECT_ROOT, "inputs", "sites_latlon.csv"),
  file.path(PROJECT_ROOT, "inputs", "coordinates.csv")
)
coords_tbl <- NULL
for (p in coords_paths) {
  if (!file.exists(p)) next
  x <- read_tabular_clean(p)
  if (is.null(x)) next
  sc <- choose_col(names(x), site_candidates)
  lc <- choose_col(names(x), lon_candidates)
  tc <- choose_col(names(x), lat_candidates)
  if (is.na(sc) || is.na(lc) || is.na(tc)) next
  coords_tbl <- x %>% transmute(site = as.character(.data[[sc]]), lon = suppressWarnings(as.numeric(.data[[lc]])), lat = suppressWarnings(as.numeric(.data[[tc]])))
  cat("[structure_pies_map] Using coords file:", p, "\n")
  break
}
if (!is.null(coords_tbl)) {
  meta_site <- meta_site %>% left_join(coords_tbl, by = "site")
} else {
  lc <- choose_col(names(meta_site), lon_candidates)
  tc <- choose_col(names(meta_site), lat_candidates)
  if (!is.na(lc) && !is.na(tc)) {
    meta_site <- meta_site %>% mutate(lon = suppressWarnings(as.numeric(.data[[lc]])), lat = suppressWarnings(as.numeric(.data[[tc]])))
  } else {
    meta_site <- meta_site %>% mutate(lon = NA_real_, lat = NA_real_)
  }
}
cat("[structure_pies_map] site coords: with=", sum(!is.na(meta_site$lon) & !is.na(meta_site$lat)), " missing=", sum(is.na(meta_site$lon) | is.na(meta_site$lat)), "\n", sep = "")

q_abs <- unique(file.path(PROJECT_ROOT, q_dirs))
q_abs <- q_abs[dir.exists(q_abs)]
if (length(q_abs) == 0) stop("q_dirs not found: ", paste(q_dirs, collapse = ", "))

q_files <- unlist(lapply(q_abs, function(d) list.files(d, pattern = "\\.(Q|q|txt|csv)$", recursive = TRUE, full.names = TRUE)), use.names = FALSE)
q_files <- unique(q_files)
q_files <- q_files[!grepl("dapc_supervised_posterior", safe_basename(q_files), ignore.case = TRUE)]
q_files <- q_files[!grepl("harvest|summary|log|evanno|delta|trace|plot|meta", safe_basename(q_files), ignore.case = TRUE)]
if (length(q_files) == 0) stop("No Q files found in q_dirs.")

parse_q <- function(path) {
  ext <- tolower(tools::file_ext(path))
  dat <- tryCatch({
    if (ext == "csv") suppressWarnings(readr::read_csv(path, show_col_types = FALSE, guess_max = 5000, locale = readr::locale(encoding = "UTF-8"), comment = "#", name_repair = "minimal"))
    else utils::read.table(path, header = FALSE, fill = TRUE, stringsAsFactors = FALSE, check.names = FALSE, comment.char = "")
  }, error = function(e) NULL)
  if (is.null(dat)) return(NULL)
  
  dat <- as.data.frame(dat, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(dat) < 3 || ncol(dat) < 2) return(NULL)
  
  nms <- names(dat); if (is.null(nms)) nms <- paste0("v", seq_len(ncol(dat)))
  names(dat) <- clean_names(nms)
  names(dat)[names(dat) == ""] <- paste0("col", seq_len(sum(names(dat) == "")))
  
  keep <- !grepl("^\\.\\.\\.[0-9]+$", names(dat))
  dat <- dat[, keep, drop = FALSE]
  
  num_like <- vapply(dat, is_numeric_like, logical(1))
  id_col <- NA_character_
  if (!num_like[1] && sum(num_like[-1]) >= 2) id_col <- names(dat)[1]
  
  q_cols <- if (!is.na(id_col)) names(dat)[-1] else names(dat)[num_like]
  if (length(q_cols) < 2) return(NULL)
  
  q <- as.data.frame(lapply(dat[, q_cols, drop = FALSE], function(z) suppressWarnings(as.numeric(as.character(z)))))
  keep_rows <- rowSums(is.na(q)) == 0
  q <- q[keep_rows, , drop = FALSE]
  if (nrow(q) < 3) return(NULL)
  
  if (!all(q >= -1e-8 & q <= 1 + 1e-8, na.rm = TRUE)) return(NULL)
  if (mean(abs(rowSums(q) - 1) <= 0.02) < 0.95) return(NULL)
  
  k <- extract_k(path); if (is.na(k)) k <- ncol(q)
  if (k > 20 || k == nrow(q)) return(NULL)
  
  ids <- NULL
  if (!is.na(id_col)) {
    ids <- as.character(dat[[id_col]])
    ids <- ids[keep_rows]
  }
  
  colnames(q) <- paste0("Q", seq_len(ncol(q)))
  list(path = path, k = k, rep = extract_rep(path), q = q, ids = ids, n = nrow(q), p = ncol(q))
}

q_list <- list()
for (f in q_files) {
  x <- parse_q(f)
  if (is.null(x)) next
  q_list[[f]] <- x
  cat("[structure_pies_map] Q retained:", f, "| K=", x$k, "| rep=", x$rep, "| n=", x$n, "| p=", x$p, "\n", sep = "")
}
if (length(q_list) == 0) stop("No valid Q files retained.")

find_ids_q <- function(target_n) {
  base_dir <- file.path(PROJECT_ROOT, "outputs", "v1", "structure_runs")
  if (!dir.exists(base_dir)) return(NULL)
  
  all_files <- list.files(base_dir, recursive = TRUE, full.names = TRUE)
  all_files <- all_files[file.exists(all_files)]
  all_files <- all_files[!dir.exists(all_files)]
  all_files <- all_files[grepl("ind|inds|id|ids|order|sample|individual|popfile|mainparams", safe_basename(all_files), ignore.case = TRUE)]
  if (length(all_files) == 0) return(NULL)
  
  best <- NULL; best_score <- -Inf
  for (f in all_files) {
    d <- read_tabular_clean(f)
    
    # line-based fallback for weird/no-extension files
    if (is.null(d)) {
      lines <- tryCatch(readLines(f, warn = FALSE), error = function(e) character(0))
      if (length(lines) > 0) {
        spl <- strsplit(trimws(lines), "\\s+")
        lens <- lengths(spl)
        keep <- lens >= 2
        if (sum(keep) > 0) {
          m <- do.call(rbind, lapply(spl[keep], function(v) c(v[1], v[2])))
          d <- data.frame(v1 = m[, 1], v2 = m[, 2], stringsAsFactors = FALSE)
        }
      }
    }
    
    if (is.null(d)) next
    nr <- nrow(d)
    if (nr < 200 || nr > 400) next
    
    idc <- choose_col(names(d), id_candidates)
    if (is.na(idc) && ncol(d) >= 2) idc <- names(d)[2]
    if (is.na(idc)) next
    
    score <- 0
    score <- score + ifelse(abs(nr - target_n) <= 2, 8, ifelse(abs(nr - target_n) <= 10, 4, 0))
    score <- score + ifelse(grepl("popfile|order|ind_ids|individual", safe_basename(f), ignore.case = TRUE), 3, 0)
    
    if (score > best_score) {
      best_score <- score
      best <- list(path = f, data = d, id_col = idc)
    }
  }
  
  if (is.null(best)) return(NULL)
  out <- best$data %>% transmute(ind_id = as.character(.data[[best$id_col]]))
  cat("[structure_pies_map] Using ids_q file:", best$path, "(n=", nrow(out), ", id_col=", best$id_col, ")\n", sep = "")
  out
}

find_meta_ind <- function(target_n, ids_q = NULL) {
  scan_dirs <- c(file.path(PROJECT_ROOT, "outputs", "v1", "structure_runs"), file.path(PROJECT_ROOT, "outputs", "v1"), file.path(PROJECT_ROOT, "inputs"))
  files <- unlist(lapply(scan_dirs, function(d) {
    if (!dir.exists(d)) return(character(0))
    list.files(d, recursive = TRUE, full.names = TRUE)
  }), use.names = FALSE)
  files <- unique(files)
  files <- files[file.exists(files)]
  files <- files[!dir.exists(files)]
  files <- files[!grepl("site_metadata\\.csv$", files, ignore.case = TRUE)]
  
  best <- NULL; best_score <- -Inf
  for (f in files) {
    d <- read_tabular_clean(f)
    if (is.null(d)) next
    nr <- nrow(d)
    if (nr < 250 || nr > 320) next
    
    idc <- choose_col(names(d), id_candidates)
    sc <- choose_col(names(d), site_candidates)
    if (is.na(sc) && "v1" %in% names(d)) sc <- "v1"
    if (is.na(idc) && "v2" %in% names(d)) idc <- "v2"
    if (is.na(idc) || is.na(sc)) next
    
    score <- 0
    score <- score + ifelse(abs(nr - target_n) <= 1, 12, ifelse(abs(nr - target_n) <= 3, 7, ifelse(abs(nr - target_n) <= 10, 3, 0)))
    score <- score + ifelse(grepl("structure_runs", f, ignore.case = TRUE), 4, 0)
    score <- score + ifelse(grepl("popfile|sample|individual|meta", safe_basename(f), ignore.case = TRUE), 2, 0)
    
    if (!is.null(ids_q) && nrow(ids_q) == nr) {
      tmp_ids <- as.character(d[[idc]])
      overlap <- mean(ids_q$ind_id %in% tmp_ids)
      score <- score + 8 * overlap
    }
    
    if (score > best_score) {
      best_score <- score
      best <- list(path = f, n = nr, id_col = idc, site_col = sc, data = d)
    }
  }
  
  if (is.null(best)) return(NULL)
  
  # hard guard: if no ids_q and selected meta_ind does not match target_n exactly, reject
  if (is.null(ids_q) && best$n != target_n) return(NULL)
  
  cat("[structure_pies_map] Using individual metadata:", best$path, "(n=", best$n, ", id_col=", best$id_col, ", site_col=", best$site_col, ")\n", sep = "")
  best$data %>% transmute(ind_id = as.character(.data[[best$id_col]]), site = as.character(.data[[best$site_col]]))
}

target_n <- max(vapply(q_list, function(x) x$n, numeric(1)))
ids_q <- find_ids_q(target_n)
meta_ind <- find_meta_ind(target_n, ids_q)
if (is.null(meta_ind)) {
  stop("No compatible individual metadata found for STRUCTURE cohort (n=", target_n, ").\nProvide an IDs/order file in outputs/v1/structure_runs/ (e.g., popfile/ind_ids).")
}

site_rep_list <- list()
for (nm in names(q_list)) {
  x <- q_list[[nm]]
  qdf <- x$q %>% mutate(row_index = row_number())
  cat("[structure_pies_map] merge for", basename(nm), "| Q n=", nrow(qdf), "| meta_ind n=", nrow(meta_ind), "\n")
  
  if (!is.null(x$ids)) {
    qdf$ind_id <- as.character(x$ids)
  } else if (!is.null(ids_q)) {
    if (nrow(ids_q) != nrow(qdf)) {
      stop("Impossible to map Q rows to individuals: provide an IDs file in outputs/v1/structure_runs/ (ind_id order) OR use the same meta_ind cohort as STRUCTURE (n=", nrow(qdf), ").")
    }
    qdf$ind_id <- ids_q$ind_id
  }
  
  merged <- NULL
  if ("ind_id" %in% names(qdf)) {
    merged <- qdf %>% left_join(meta_ind, by = "ind_id")
    matched <- sum(!is.na(merged$site))
    unmatched <- nrow(merged) - matched
    bad <- unique(merged$ind_id[is.na(merged$site)]); bad <- bad[!is.na(bad)]
    cat("[structure_pies_map] matched=", matched, " unmatched=", unmatched, " bad_ids=", ifelse(length(bad)==0,"none",paste(head(bad,10), collapse=", ")), "\n", sep = "")
    merged <- merged %>% filter(!is.na(site))
  } else {
    if (nrow(qdf) != nrow(meta_ind)) {
      stop("Impossible to map Q rows to individuals: provide an IDs file in outputs/v1/structure_runs/ (ind_id order) OR use the same meta_ind cohort as STRUCTURE (n=", nrow(qdf), ").")
    }
    merged <- qdf %>% left_join(meta_ind %>% mutate(row_index = row_number()), by = "row_index")
  }
  
  if (nrow(merged) == 0) next
  
  site_rep <- merged %>%
    group_by(site) %>%
    summarise(n_ind = n(), across(starts_with("Q"), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
    mutate(K = x$k, rep = x$rep, q_file = basename(nm))
  
  site_rep_list[[nm]] <- site_rep
}

if (length(site_rep_list) == 0) stop("No valid K datasets after merge checks.")

site_rep_df <- bind_rows(site_rep_list)
cluster_cols_all <- grep("^Q[0-9]+$", names(site_rep_df), value = TRUE)

site_k_df <- site_rep_df %>%
  group_by(site, K) %>%
  summarise(
    n_ind = mean(n_ind, na.rm = TRUE),
    n_reps = n_distinct(rep),
    across(all_of(cluster_cols_all), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  left_join(meta_site, by = "site")

if (nrow(site_k_df) == 0) stop("No site-level rows after aggregation.")

max_k <- max(site_k_df$K, na.rm = TRUE)
cluster_cols <- paste0("Q", seq_len(max_k))
for (cc in cluster_cols) if (!(cc %in% names(site_k_df))) site_k_df[[cc]] <- 0
pal <- setNames(scales::hue_pal()(max_k), cluster_cols)

has_coords <- any(!is.na(site_k_df$lon) & !is.na(site_k_df$lat))
subtitle_reps <- paste0("Pies = mean Q by site; replicate-averaged per K (n_reps range: ", min(site_k_df$n_reps, na.rm = TRUE), "-", max(site_k_df$n_reps, na.rm = TRUE), ")")

if (has_coords) {
  miss <- site_k_df %>% filter(is.na(lon) | is.na(lat)) %>% distinct(site)
  if (nrow(miss) > 0) warning("Sites excluded from spatial plot (missing lon/lat): ", paste(miss$site, collapse = ", "))
  plot_df <- site_k_df %>% filter(!is.na(lon), !is.na(lat)) %>% mutate(K_label = paste0("K=", K), radius_raw = sqrt(n_ind))
  if (nrow(plot_df) == 0) has_coords <- FALSE
}

if (has_coords) {
  span <- max(diff(range(plot_df$lon, na.rm = TRUE)), diff(range(plot_df$lat, na.rm = TRUE)))
  if (!is.finite(span) || span <= 0) span <- 1
  plot_df <- plot_df %>% mutate(radius = radius_raw * (0.03 * span / max(radius_raw, na.rm = TRUE)))
  
  p <- ggplot()
  if (use_sf && use_rnaturalearth) {
    base <- tryCatch({
      w <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
      bb <- sf::st_bbox(c(xmin = min(plot_df$lon) - 1, xmax = max(plot_df$lon) + 1, ymin = min(plot_df$lat) - 1, ymax = max(plot_df$lat) + 1), crs = sf::st_crs(4326))
      geom_sf(data = sf::st_crop(w, bb), fill = "grey97", color = "grey70", linewidth = 0.2)
    }, error = function(e) NULL)
    if (!is.null(base)) p <- p + base
  }
  
  if (use_scatterpie) {
    p <- p + scatterpie::geom_scatterpie(data = plot_df, aes(x = lon, y = lat, r = radius), cols = cluster_cols, color = "grey25", linewidth = 0.2, alpha = 0.95)
  } else {
    long <- plot_df %>%
      pivot_longer(cols = all_of(cluster_cols), names_to = "cluster", values_to = "q") %>%
      group_by(site, K, K_label, lon, lat, n_ind, radius_raw) %>%
      slice_max(q, n = 1, with_ties = FALSE) %>%
      ungroup()
    p <- p + geom_point(data = long, aes(x = lon, y = lat, size = radius_raw, fill = cluster), shape = 21, color = "grey20")
  }
  
  bks <- pretty(plot_df$n_ind, n = 3); bks <- bks[bks > 0]
  if (length(bks) == 0) bks <- sort(unique(plot_df$n_ind))[1]
  
  p <- p +
    geom_point(data = plot_df, aes(x = lon, y = lat, size = n_ind), alpha = 0, color = NA) +
    scale_fill_manual(values = pal, name = "Clusters") +
    scale_size_continuous(name = "Samples (n)", breaks = bks) +
    facet_wrap(~K_label, ncol = ifelse(length(unique(plot_df$K)) <= 4, 2, 3)) +
    coord_equal() +
    labs(title = "STRUCTURE mean Q by site", subtitle = subtitle_reps, x = "Longitude", y = "Latitude") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold"), strip.text = element_text(face = "bold"), panel.grid.minor = element_blank(), legend.position = "right")
  
  if (use_ggspatial) {
    p <- p +
      ggspatial::annotation_scale(location = "bl", width_hint = 0.25) +
      ggspatial::annotation_north_arrow(location = "tl", which_north = "true", style = ggspatial::north_arrow_fancy_orienteering)
  }
  
  pdf_path <- file.path(OUT_DIR, "structure_pies_map.pdf")
  png_path <- file.path(OUT_DIR, "structure_pies_map.png")
  ggsave(pdf_path, p, width = 12, height = 8, device = grDevices::cairo_pdf)
  ggsave(png_path, p, width = 12, height = 8, dpi = 400)
  cat("[structure_pies_map] Figure exported:\n - ", pdf_path, "\n - ", png_path, "\n", sep = "")
  
} else {
  cat("[structure_pies_map] No lon/lat found: generated non-spatial pies. To enable map, add inputs/site_coordinates.csv with columns site, lon, lat.\n")
  
  plot_df <- site_k_df %>% mutate(site = factor(site, levels = unique(site[order(region_ns, site)])), x_site = as.numeric(site), y_site = 0, K_label = paste0("K=", K), radius_raw = sqrt(n_ind))
  plot_df <- plot_df %>% mutate(radius = radius_raw * (0.35 / max(radius_raw, na.rm = TRUE)))
  
  p <- ggplot()
  if (use_scatterpie) {
    p <- p + scatterpie::geom_scatterpie(data = plot_df, aes(x = x_site, y = y_site, r = radius), cols = cluster_cols, color = "grey25", linewidth = 0.2, alpha = 0.95)
  } else {
    long <- plot_df %>%
      pivot_longer(cols = all_of(cluster_cols), names_to = "cluster", values_to = "q") %>%
      group_by(site, K, K_label, x_site, y_site, n_ind, radius_raw) %>%
      slice_max(q, n = 1, with_ties = FALSE) %>%
      ungroup()
    p <- p + geom_point(data = long, aes(x = x_site, y = y_site, size = radius_raw, fill = cluster), shape = 21, color = "grey20")
  }
  
  bks <- pretty(plot_df$n_ind, n = 3); bks <- bks[bks > 0]
  if (length(bks) == 0) bks <- sort(unique(plot_df$n_ind))[1]
  
  p <- p +
    geom_point(data = plot_df, aes(x = x_site, y = y_site, size = n_ind), alpha = 0, color = NA) +
    geom_text(data = plot_df %>% distinct(site, x_site), aes(x = x_site, y = -0.75, label = site), angle = 45, hjust = 1, size = 3) +
    scale_fill_manual(values = pal, name = "Clusters") +
    scale_size_continuous(name = "Samples (n)", breaks = bks) +
    facet_wrap(~K_label, ncol = ifelse(length(unique(plot_df$K)) <= 4, 2, 3)) +
    scale_x_continuous(breaks = NULL) +
    scale_y_continuous(breaks = NULL) +
    labs(title = "STRUCTURE mean Q by site (non-spatial fallback)", subtitle = subtitle_reps, x = NULL, y = NULL) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold"), strip.text = element_text(face = "bold"), panel.grid = element_blank(), legend.position = "right")
  
  pdf_path <- file.path(OUT_DIR, "structure_pies_no_coords.pdf")
  png_path <- file.path(OUT_DIR, "structure_pies_no_coords.png")
  ggsave(pdf_path, p, width = 12, height = 8, device = grDevices::cairo_pdf)
  ggsave(png_path, p, width = 12, height = 8, dpi = 400)
  cat("[structure_pies_map] Figure exported:\n - ", pdf_path, "\n - ", png_path, "\n", sep = "")
}
