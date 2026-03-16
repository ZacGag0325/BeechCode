# scripts/03_structure_all_plots_S2N_bySite.R
############################################################
# STRUCTURE helper outputs (final individual barplots only)
# - Reads externally generated STRUCTURE Q files
# - Builds ONE final individual plot per K (no per-run figures)
# - Individuals ordered SOUTH -> NORTH by site latitude
# Outputs:
# - outputs/figures/structure_individual_barplot_K{K}.jpeg
# - outputs/tables/supplementary/structure_run_inventory.csv
# - outputs/tables/supplementary/structure_selected_runs.csv
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

source("scripts/_load_objects.R")

message("[03_structure] Preparing final STRUCTURE individual barplots...")

resolve_col <- function(df, choices) {
  nms <- names(df)
  idx <- match(TRUE, tolower(nms) %in% tolower(choices), nomatch = 0)
  if (idx == 0) return(NA_character_)
  nms[idx]
}

# ----------------------------
# 1) Individual IDs and Site map
# ----------------------------
id_col_dfids <- resolve_col(df_ids, c("ind", "individual", "sample", "sampleid", "id"))
site_col_dfids <- resolve_col(df_ids, c("site", "population", "pop"))
if (is.na(id_col_dfids) || is.na(site_col_dfids)) {
  stop("[03_structure] df_ids must contain individual ID and Site columns.")
}

ids_default <- as.character(df_ids[[id_col_dfids]])
site_map <- setNames(as.character(df_ids[[site_col_dfids]]), ids_default)

load_individual_order <- function(default_ids) {
  ids_path <- file.path(PROJECT_ROOT, "data", "structure", "ids_order_from_raw.csv")
  
  if (!file.exists(ids_path)) {
    return(default_ids)
  }
  
  ids_df <- tryCatch(
    read.csv(ids_path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL
  )
  
  if (is.null(ids_df) || nrow(ids_df) == 0) {
    warning("[03_structure] Could not read ids_order_from_raw.csv; using df_ids order.")
    return(default_ids)
  }
  
  id_col <- resolve_col(ids_df, c("ind", "individual", "sample", "sampleid", "id"))
  if (is.na(id_col)) id_col <- names(ids_df)[1]
  
  ordered_ids <- trimws(as.character(ids_df[[id_col]]))
  ordered_ids <- ordered_ids[nzchar(ordered_ids)]
  ordered_ids <- unique(ordered_ids)
  
  matched <- ordered_ids[ordered_ids %in% default_ids]
  if (length(matched) == 0) {
    warning("[03_structure] ids_order_from_raw.csv did not match df_ids IDs; using df_ids order.")
    return(default_ids)
  }
  
  missing <- setdiff(default_ids, matched)
  c(matched, missing)
}

id_reference <- load_individual_order(ids_default)
expected_n <- length(id_reference)

# ----------------------------
# 2) Site latitude map (for SOUTH -> NORTH ordering)
# ----------------------------
load_site_latitude <- function(meta_df) {
  site_col <- resolve_col(meta_df, c("site", "population", "pop"))
  lat_col <- resolve_col(meta_df, c("latitude", "lat"))
  
  if (is.na(site_col) || is.na(lat_col)) {
    fallback <- file.path(PROJECT_ROOT, "inputs", "site_metadata.csv")
    if (!file.exists(fallback)) {
      stop("[03_structure] Could not find Site/Latitude in meta or inputs/site_metadata.csv")
    }
    meta_df <- read.csv(fallback, stringsAsFactors = FALSE, check.names = FALSE)
    site_col <- resolve_col(meta_df, c("site", "population", "pop"))
    lat_col <- resolve_col(meta_df, c("latitude", "lat"))
    if (is.na(site_col) || is.na(lat_col)) {
      stop("[03_structure] site_metadata.csv must contain Site and Latitude columns.")
    }
  }
  
  out <- data.frame(
    Site = trimws(as.character(meta_df[[site_col]])),
    Latitude = suppressWarnings(as.numeric(meta_df[[lat_col]])),
    stringsAsFactors = FALSE
  ) %>%
    filter(nzchar(Site), !is.na(Latitude)) %>%
    group_by(Site) %>%
    summarise(Latitude = mean(Latitude, na.rm = TRUE), .groups = "drop")
  
  if (nrow(out) == 0) {
    stop("[03_structure] No valid site latitude information available.")
  }
  
  setNames(out$Latitude, out$Site)
}

site_lat_map <- load_site_latitude(meta)

# ----------------------------
# 3) Q file discovery
# ----------------------------
find_q_files <- function() {
  q_root <- file.path(PROJECT_ROOT, "data", "structure", "Q_Files")
  if (!dir.exists(q_root)) {
    q_root_alt <- file.path(PROJECT_ROOT, "data", "structure", "Q_files")
    if (dir.exists(q_root_alt)) q_root <- q_root_alt
  }
  
  message("[03_structure] Searching folder: ", q_root)
  
  if (!dir.exists(q_root)) {
    return(list(root = q_root, files = character(0), ignored = character(0)))
  }
  
  files <- list.files(
    q_root,
    recursive = TRUE,
    full.names = TRUE,
    pattern = "(?i)\\.(q|txt|csv)$"
  )
  
  files <- files[file.info(files)$isdir %in% FALSE]
  b <- basename(files)
  
  ignore_idx <- grepl("(?i)summary|extracted|evanno|likelihood|readme", b)
  ignored <- unique(files[ignore_idx])
  kept <- unique(files[!ignore_idx])
  
  list(root = q_root, files = kept, ignored = ignored)
}

extract_k_run <- function(path) {
  b <- basename(path)
  k_hit <- regmatches(b, regexpr("(?i)K[0-9]+", b, perl = TRUE))
  K <- if (length(k_hit) == 0 || !nzchar(k_hit)) NA_integer_ else as.integer(gsub("[^0-9]", "", k_hit))
  
  run <- suppressWarnings(as.integer(sub("(?i).*K[0-9]+[-_]?([0-9]+).*$", "\\1", b, perl = TRUE)))
  if (is.na(run)) {
    run <- suppressWarnings(as.integer(sub("(?i).*[_-]run[_-]?([0-9]+).*$", "\\1", b, perl = TRUE)))
  }
  
  list(K = K, run = run)
}

# ----------------------------
# 4) Robust Q parser
# ----------------------------
numeric_columns <- function(df) {
  num_df <- as.data.frame(lapply(df, function(col) {
    x <- trimws(as.character(col))
    x[x == ""] <- NA_character_
    suppressWarnings(as.numeric(x))
  }), check.names = FALSE, stringsAsFactors = FALSE)
  
  keep <- vapply(seq_len(ncol(num_df)), function(j) {
    original <- trimws(as.character(df[[j]]))
    original[original == ""] <- NA_character_
    converted <- num_df[[j]]
    all(is.na(original) | !is.na(converted))
  }, logical(1))
  
  num_df[, keep, drop = FALSE]
}

read_q_matrix <- function(path) {
  readers <- list(
    list(label = "whitespace_no_header", fun = function() read.table(path, header = FALSE, sep = "", stringsAsFactors = FALSE, check.names = FALSE, fill = TRUE, comment.char = "", quote = "")),
    list(label = "whitespace_header",    fun = function() read.table(path, header = TRUE,  sep = "", stringsAsFactors = FALSE, check.names = FALSE, fill = TRUE, comment.char = "", quote = "")),
    list(label = "tab_no_header",        fun = function() read.table(path, header = FALSE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, fill = TRUE, comment.char = "", quote = "")),
    list(label = "tab_header",           fun = function() read.table(path, header = TRUE,  sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, fill = TRUE, comment.char = "", quote = "")),
    list(label = "csv_header",           fun = function() read.csv(path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)),
    list(label = "csv_no_header",        fun = function() read.csv(path, header = FALSE, stringsAsFactors = FALSE, check.names = FALSE))
  )
  
  for (r in readers) {
    dat <- tryCatch(r$fun(), error = function(e) NULL)
    if (is.null(dat) || nrow(dat) == 0 || ncol(dat) == 0) next
    
    q <- numeric_columns(dat)
    if (ncol(q) < 1) next
    
    if (ncol(q) >= 2) {
      first_col <- q[[1]]
      if (!anyNA(first_col) && all(abs(first_col - seq_len(length(first_col))) < 1e-8)) {
        q <- q[, -1, drop = FALSE]
      }
    }
    
    q <- q[rowSums(!is.na(q)) > 0, , drop = FALSE]
    if (nrow(q) == 0 || ncol(q) < 1) next
    
    mat <- as.matrix(q)
    storage.mode(mat) <- "numeric"
    colnames(mat) <- paste0("Q", seq_len(ncol(mat)))
    
    return(list(ok = TRUE, matrix = mat, reader = r$label))
  }
  
  list(ok = FALSE, matrix = NULL, reader = NA_character_)
}

# ----------------------------
# 5) Parse all runs and choose ONE final run per K
#    (no across-run averaging because cluster-label alignment is not implemented)
# ----------------------------
scan <- find_q_files()
if (length(scan$ignored) > 0) {
  message("[03_structure] Ignored helper files: ", length(scan$ignored))
}

if (length(scan$files) == 0) {
  message("[03_structure] No external STRUCTURE Q files found; checked: ", scan$root)
  message("[03_structure] Done.")
} else {
  parsed_runs <- list()
  run_inventory <- list()
  
  for (f in scan$files) {
    meta_name <- extract_k_run(f)
    parsed <- read_q_matrix(f)
    
    if (!parsed$ok) {
      run_inventory[[length(run_inventory) + 1]] <- data.frame(
        file = basename(f),
        full_path = f,
        K = meta_name$K,
        run = meta_name$run,
        reader = parsed$reader,
        n_rows = NA_integer_,
        n_cols = NA_integer_,
        mean_abs_row_sum_deviation = NA_real_,
        status = "unreadable",
        stringsAsFactors = FALSE
      )
      next
    }
    
    q_mat <- parsed$matrix
    K_detected <- ncol(q_mat)
    K_final <- if (!is.na(meta_name$K)) meta_name$K else K_detected
    
    if (!is.na(meta_name$K) && meta_name$K != K_detected) {
      K_final <- K_detected
    }
    
    if (K_final == 1) {
      run_inventory[[length(run_inventory) + 1]] <- data.frame(
        file = basename(f),
        full_path = f,
        K = K_final,
        run = meta_name$run,
        reader = parsed$reader,
        n_rows = nrow(q_mat),
        n_cols = ncol(q_mat),
        mean_abs_row_sum_deviation = mean(abs(rowSums(q_mat, na.rm = TRUE) - 1), na.rm = TRUE),
        status = "skipped_k1",
        stringsAsFactors = FALSE
      )
      next
    }
    
    if (nrow(q_mat) != expected_n) {
      run_inventory[[length(run_inventory) + 1]] <- data.frame(
        file = basename(f),
        full_path = f,
        K = K_final,
        run = meta_name$run,
        reader = parsed$reader,
        n_rows = nrow(q_mat),
        n_cols = ncol(q_mat),
        mean_abs_row_sum_deviation = mean(abs(rowSums(q_mat, na.rm = TRUE) - 1), na.rm = TRUE),
        status = "row_mismatch_skipped",
        stringsAsFactors = FALSE
      )
      next
    }
    
    mad_row_sum <- mean(abs(rowSums(q_mat, na.rm = TRUE) - 1), na.rm = TRUE)
    
    run_inventory[[length(run_inventory) + 1]] <- data.frame(
      file = basename(f),
      full_path = f,
      K = K_final,
      run = meta_name$run,
      reader = parsed$reader,
      n_rows = nrow(q_mat),
      n_cols = ncol(q_mat),
      mean_abs_row_sum_deviation = mad_row_sum,
      status = "parsed",
      stringsAsFactors = FALSE
    )
    
    parsed_runs[[length(parsed_runs) + 1]] <- list(
      file = f,
      file_base = basename(f),
      K = K_final,
      run = meta_name$run,
      reader = parsed$reader,
      q = q_mat,
      mad = mad_row_sum
    )
  }
  
  if (length(run_inventory) > 0) {
    inv_df <- bind_rows(run_inventory) %>% arrange(K, run, file)
    inv_file <- file.path(TABLES_SUPP_DIR, "structure_run_inventory.csv")
    write.csv(inv_df, inv_file, row.names = FALSE)
    message("[03_structure] Saved run inventory: ", inv_file)
  }
  
  if (length(parsed_runs) == 0) {
    message("[03_structure] No usable K>=2 Q runs after validation.")
    message("[03_structure] Done.")
  } else {
    runs_by_k <- split(parsed_runs, sapply(parsed_runs, function(x) x$K))
    selected_rows <- list()
    
    for (k_name in names(runs_by_k)) {
      k_runs <- runs_by_k[[k_name]]
      k_num <- as.integer(k_name)
      
      score_df <- bind_rows(lapply(k_runs, function(x) {
        data.frame(
          file = x$file,
          file_base = x$file_base,
          K = x$K,
          run = ifelse(is.na(x$run), Inf, x$run),
          run_report = x$run,
          reader = x$reader,
          mad = x$mad,
          stringsAsFactors = FALSE
        )
      })) %>%
        arrange(mad, run, file_base)
      
      best <- k_runs[[match(score_df$file[1], sapply(k_runs, function(x) x$file))]]
      
      message("[03_structure] Found ", length(k_runs), " usable run(s) for K=", k_num)
      message("[03_structure] Using best run for K=", k_num,
              " (run=", ifelse(is.na(best$run), "NA", best$run),
              ", reader=", best$reader,
              ", mean|rowSum-1|=", signif(best$mad, 4), ")")
      
      q_df <- data.frame(
        Individual = id_reference,
        best$q,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
      q_df$Site <- site_map[q_df$Individual]
      q_df$SiteLat <- as.numeric(site_lat_map[q_df$Site])
      
      if (any(is.na(q_df$Site))) {
        warning("[03_structure] Some individuals have missing Site labels for K=", k_num)
      }
      
      if (any(is.na(q_df$SiteLat))) {
        warning("[03_structure] Some Site latitudes are missing for K=", k_num, "; these individuals will be placed at end.")
      }
      
      q_df <- q_df %>%
        mutate(
          Site = as.character(Site),
          SiteLat = ifelse(is.na(SiteLat), Inf, SiteLat)
        ) %>%
        arrange(SiteLat, Site, Individual)
      
      q_df$PlotIndex <- seq_len(nrow(q_df))
      q_cols <- grep("^Q", names(q_df), value = TRUE)
      
      plot_df <- q_df %>%
        select(PlotIndex, Site, all_of(q_cols)) %>%
        pivot_longer(cols = all_of(q_cols), names_to = "Cluster", values_to = "Q")
      
      site_blocks <- q_df %>%
        count(Site, name = "n") %>%
        mutate(
          xmax = cumsum(n),
          xmin = xmax - n + 1,
          xmid = (xmin + xmax) / 2
        )
      
      separators <- site_blocks$xmax[-nrow(site_blocks)] + 0.5
      
      p <- ggplot(plot_df, aes(x = PlotIndex, y = Q, fill = Cluster)) +
        geom_col(width = 1) +
        geom_vline(xintercept = separators, linewidth = 0.25, color = "grey25") +
        scale_x_continuous(
          breaks = site_blocks$xmid,
          labels = site_blocks$Site,
          expand = c(0, 0)
        ) +
        scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
        theme_bw(base_size = 11) +
        theme(
          panel.grid = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
          axis.ticks.x = element_blank(),
          legend.position = "right",
          plot.title = element_text(face = "bold")
        ) +
        labs(
          title = sprintf("STRUCTURE ancestry barplot (K=%d)", k_num),
          subtitle = "Individuals ordered south to north by site latitude",
          x = "Site blocks (south -> north)",
          y = "Ancestry proportion"
        )
      
      out_fig <- file.path(FIGURES_DIR, sprintf("structure_individual_barplot_K%d.jpeg", k_num))
      ggsave(out_fig, p, width = 13, height = 5.5, dpi = 320)
      message("[03_structure] Saved: ", out_fig)
      
      selected_rows[[length(selected_rows) + 1]] <- data.frame(
        K = k_num,
        selected_file = best$file_base,
        selected_path = best$file,
        selected_run = best$run,
        reader = best$reader,
        mean_abs_row_sum_deviation = best$mad,
        n_runs_available = length(k_runs),
        method = "best_single_run_no_cluster_relabeling",
        stringsAsFactors = FALSE
      )
    }
    
    if (length(selected_rows) > 0) {
      selected_df <- bind_rows(selected_rows) %>% arrange(K)
      selected_file <- file.path(TABLES_SUPP_DIR, "structure_selected_runs.csv")
      write.csv(selected_df, selected_file, row.names = FALSE)
      message("[03_structure] Saved selected-run summary: ", selected_file)
    }
    
    message("[03_structure] Done.")
  }
}