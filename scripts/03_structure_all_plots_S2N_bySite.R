# scripts/03_structure_all_plots_S2N_bySite.R
############################################################
# STRUCTURE helper outputs only (NO STRUCTURE reruns)
# - Uses externally generated Q tables
# - Produces mean Q by site for each available K
# - Produces site-level dominant cluster summary
# - Produces individual and site-level STRUCTURE barplots
# Outputs (examples):
# - outputs/tables/supplementary/structure_meanQ_by_site_K2.csv
# - outputs/tables/supplementary/structure_meanQ_by_site_K2_run1.csv
# - outputs/tables/supplementary/structure_parsed_files_summary.csv
# - outputs/figures/supplementary/structure_individual_barplot_K2_run1.jpeg
# - outputs/figures/supplementary/structure_site_mean_barplot_K2_run1.jpeg
# - outputs/figures/supplementary/structure_site_mean_barplot_K2.jpeg
# - outputs/tables/supplementary/structure_site_cluster_summary.csv
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

source("scripts/_load_objects.R")

message("[03_structure] Preparing STRUCTURE helper outputs from external Q files...")

# STRUCTURE was run externally and was not clone-corrected.
# We do not rerun STRUCTURE here; this script only summarizes existing Q outputs.

resolve_col <- function(df, choices) {
  nms <- names(df)
  idx <- match(TRUE, tolower(nms) %in% tolower(choices), nomatch = 0)
  if (idx == 0) return(NA_character_)
  nms[idx]
}

id_col_dfids <- resolve_col(df_ids, c("ind", "individual", "sample", "sampleid", "id"))
site_col_dfids <- resolve_col(df_ids, c("site", "population", "pop"))
if (is.na(id_col_dfids) || is.na(site_col_dfids)) {
  stop("[03_structure] df_ids must contain individual ID and Site columns.")
}

site_map <- setNames(as.character(df_ids[[site_col_dfids]]), as.character(df_ids[[id_col_dfids]]))

load_individual_order <- function(default_ids) {
  ids_path <- file.path(PROJECT_ROOT, "data", "structure", "ids_order_from_raw.csv")
  if (!file.exists(ids_path)) {
    message("[03_structure] ids_order_from_raw.csv not found; using df_ids order.")
    return(default_ids)
  }
  
  ids_df <- tryCatch(
    read.csv(ids_path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL
  )
  if (is.null(ids_df) || nrow(ids_df) == 0) {
    warning("[03_structure] Could not read IDs order file; using df_ids order.")
    return(default_ids)
  }
  
  id_col <- resolve_col(ids_df, c("ind", "individual", "sample", "sampleid", "id"))
  if (is.na(id_col)) {
    id_col <- names(ids_df)[1]
  }
  
  ordered_ids <- trimws(as.character(ids_df[[id_col]]))
  ordered_ids <- ordered_ids[nzchar(ordered_ids)]
  ordered_ids <- unique(ordered_ids)
  
  if (length(ordered_ids) == 0) {
    warning("[03_structure] IDs order file has no valid IDs; using df_ids order.")
    return(default_ids)
  }
  
  matched <- ordered_ids[ordered_ids %in% default_ids]
  missing <- setdiff(default_ids, matched)
  
  if (length(matched) == 0) {
    warning("[03_structure] IDs order file did not match df_ids IDs; using df_ids order.")
    return(default_ids)
  }
  
  if (length(missing) > 0) {
    message("[03_structure] IDs order file missing ", length(missing), " IDs from df_ids; appending missing IDs at end.")
    matched <- c(matched, missing)
  }
  
  message("[03_structure] Loaded individual order from: ", ids_path)
  matched
}

id_reference <- load_individual_order(as.character(df_ids[[id_col_dfids]]))
expected_n <- length(id_reference)
message("[03_structure] Expected individuals for Q alignment: ", expected_n)

extract_k_run_from_filename <- function(path) {
  b <- basename(path)
  
  k_hit <- regmatches(b, regexpr("(?i)K[0-9]+", b, perl = TRUE))
  K <- if (length(k_hit) == 0 || !nzchar(k_hit)) NA_integer_ else as.integer(gsub("[^0-9]", "", k_hit))
  
  # Example supported: STRUCTURE_ZG_HEG-K2-1_f.Q -> run 1
  run <- suppressWarnings(as.integer(sub("(?i).*K[0-9]+[-_]?([0-9]+).*$", "\\1", b, perl = TRUE)))
  if (is.na(run)) {
    run <- suppressWarnings(as.integer(sub("(?i).*[_-]run[_-]?([0-9]+).*$", "\\1", b, perl = TRUE)))
  }
  
  list(K = K, run = run)
}

find_external_q_files <- function() {
  q_root <- file.path(PROJECT_ROOT, "data", "structure", "Q_Files")
  if (!dir.exists(q_root)) {
    q_root_alt <- file.path(PROJECT_ROOT, "data", "structure", "Q_files")
    if (dir.exists(q_root_alt)) q_root <- q_root_alt
  }
  
  message("[03_structure] Searching for STRUCTURE Q files in: ", q_root)
  
  if (!dir.exists(q_root)) {
    return(list(root = q_root, files = character(0), ignored = character(0)))
  }
  
  all_candidates <- list.files(
    q_root,
    recursive = TRUE,
    full.names = TRUE,
    pattern = "(?i)\\.(q|txt|csv)$"
  )
  all_candidates <- all_candidates[file.info(all_candidates)$isdir %in% FALSE]
  
  bnames <- basename(all_candidates)
  
  # Intentionally ignore helper/summary/non-matrix files
  ignore_idx <- grepl("(?i)summary|extracted|evanno|likelihood|readme", bnames)
  ignored <- unique(all_candidates[ignore_idx])
  q_files <- unique(all_candidates[!ignore_idx])
  
  list(root = q_root, files = q_files, ignored = ignored)
}

select_numeric_columns <- function(df) {
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

attempt_readers <- function(path) {
  list(
    list(label = "whitespace_no_header", fun = function() read.table(path, header = FALSE, sep = "", stringsAsFactors = FALSE, check.names = FALSE, fill = TRUE, comment.char = "", quote = "")),
    list(label = "whitespace_header",    fun = function() read.table(path, header = TRUE,  sep = "", stringsAsFactors = FALSE, check.names = FALSE, fill = TRUE, comment.char = "", quote = "")),
    list(label = "tab_no_header",        fun = function() read.table(path, header = FALSE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, fill = TRUE, comment.char = "", quote = "")),
    list(label = "tab_header",           fun = function() read.table(path, header = TRUE,  sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, fill = TRUE, comment.char = "", quote = "")),
    list(label = "csv_header",           fun = function() read.csv(path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)),
    list(label = "csv_no_header",        fun = function() read.csv(path, header = FALSE, stringsAsFactors = FALSE, check.names = FALSE))
  )
}

read_q_file <- function(path) {
  parsed <- NULL
  reader_label <- NA_character_
  
  for (reader in attempt_readers(path)) {
    dat <- tryCatch(reader$fun(), error = function(e) NULL)
    if (is.null(dat) || nrow(dat) == 0 || ncol(dat) == 0) next
    
    q_num <- select_numeric_columns(dat)
    if (ncol(q_num) < 1) next
    
    # Remove index-like leading column (1..N) if present and there are other columns
    if (ncol(q_num) >= 2) {
      first_col <- q_num[[1]]
      if (!anyNA(first_col) && all(abs(first_col - seq_len(length(first_col))) < 1e-8)) {
        q_num <- q_num[, -1, drop = FALSE]
      }
    }
    
    # Drop completely empty rows
    keep_rows <- rowSums(!is.na(q_num)) > 0
    q_num <- q_num[keep_rows, , drop = FALSE]
    
    if (ncol(q_num) < 1 || nrow(q_num) == 0) next
    
    parsed <- q_num
    reader_label <- reader$label
    break
  }
  
  if (is.null(parsed)) {
    return(list(ok = FALSE, reason = "unreadable/non-tabular", data = NULL, reader = NA_character_))
  }
  
  K_detected <- ncol(parsed)
  names(parsed) <- paste0("Q", seq_len(K_detected))
  
  # Row alignment by known individual order
  n_rows <- nrow(parsed)
  if (n_rows == expected_n) {
    individuals <- id_reference
  } else if (n_rows < expected_n) {
    warning("[03_structure] ", basename(path), ": row count (", n_rows, ") < expected individuals (", expected_n, "). Using first ", n_rows, " IDs from order.")
    individuals <- id_reference[seq_len(n_rows)]
  } else {
    warning("[03_structure] ", basename(path), ": row count (", n_rows, ") > expected individuals (", expected_n, "). Truncating to expected size.")
    parsed <- parsed[seq_len(expected_n), , drop = FALSE]
    individuals <- id_reference
  }
  
  q_df <- data.frame(Individual = individuals, parsed, stringsAsFactors = FALSE, check.names = FALSE)
  q_df$Site <- site_map[q_df$Individual]
  
  q_cols <- grep("^Q", names(q_df), value = TRUE)
  row_sums <- rowSums(q_df[, q_cols, drop = FALSE], na.rm = TRUE)
  
  list(
    ok = TRUE,
    reason = "parsed",
    data = q_df,
    reader = reader_label,
    K_detected = K_detected,
    n_rows = nrow(q_df),
    n_cols = length(q_cols),
    row_sum_min = suppressWarnings(min(row_sums, na.rm = TRUE)),
    row_sum_max = suppressWarnings(max(row_sums, na.rm = TRUE)),
    row_sums_ok = all(abs(row_sums - 1) <= 0.05, na.rm = TRUE)
  )
}

write_site_mean_q <- function(q_df, k, run = NA_integer_) {
  out <- q_df %>%
    filter(!is.na(Site)) %>%
    group_by(Site) %>%
    summarise(across(starts_with("Q"), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
    arrange(Site)
  
  if (!is.na(run)) {
    out_file <- file.path(TABLES_SUPP_DIR, sprintf("structure_meanQ_by_site_K%d_run%d.csv", k, run))
  } else {
    out_file <- file.path(TABLES_SUPP_DIR, sprintf("structure_meanQ_by_site_K%d.csv", k))
  }
  
  write.csv(out, out_file, row.names = FALSE)
  message("[03_structure] Saved: ", out_file)
  out
}

plot_individual_q <- function(q_df, k, run = NA_integer_) {
  q_cols <- grep("^Q", names(q_df), value = TRUE)
  if (length(q_cols) < 2) return(invisible(NULL))
  
  plot_df <- q_df %>%
    mutate(Individual = factor(Individual, levels = unique(Individual))) %>%
    select(Individual, all_of(q_cols)) %>%
    pivot_longer(cols = all_of(q_cols), names_to = "Cluster", values_to = "Q")
  
  p <- ggplot(plot_df, aes(x = Individual, y = Q, fill = Cluster)) +
    geom_col(width = 1) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "right"
    ) +
    labs(
      title = if (is.na(run)) sprintf("STRUCTURE ancestry proportions (K=%d)", k) else sprintf("STRUCTURE ancestry proportions (K=%d, run=%d)", k, run),
      x = "Individuals (ordered)",
      y = "Ancestry proportion"
    )
  
  out_file <- if (is.na(run)) {
    file.path(FIGURES_SUPP_DIR, sprintf("structure_individual_barplot_K%d.jpeg", k))
  } else {
    file.path(FIGURES_SUPP_DIR, sprintf("structure_individual_barplot_K%d_run%d.jpeg", k, run))
  }
  
  ggsave(out_file, p, width = 12, height = 4.5, dpi = 320)
  message("[03_structure] Saved: ", out_file)
}

plot_site_mean_q <- function(site_df, k, run = NA_integer_) {
  q_cols <- grep("^Q", names(site_df), value = TRUE)
  if (length(q_cols) < 2) return(invisible(NULL))
  
  plot_df <- site_df %>%
    mutate(Site = factor(Site, levels = sort(unique(Site)))) %>%
    pivot_longer(cols = all_of(q_cols), names_to = "Cluster", values_to = "Q")
  
  p <- ggplot(plot_df, aes(x = Site, y = Q, fill = Cluster)) +
    geom_col() +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      title = if (is.na(run)) sprintf("Site mean ancestry proportions (K=%d)", k) else sprintf("Site mean ancestry proportions (K=%d, run=%d)", k, run),
      x = "Site",
      y = "Mean ancestry proportion"
    )
  
  out_file <- if (is.na(run)) {
    file.path(FIGURES_SUPP_DIR, sprintf("structure_site_mean_barplot_K%d.jpeg", k))
  } else {
    file.path(FIGURES_SUPP_DIR, sprintf("structure_site_mean_barplot_K%d_run%d.jpeg", k, run))
  }
  
  ggsave(out_file, p, width = 8, height = 5.5, dpi = 320)
  message("[03_structure] Saved: ", out_file)
}

file_scan <- find_external_q_files()
q_files_found <- file_scan$files

if (length(file_scan$ignored) > 0) {
  message("[03_structure] Ignoring ", length(file_scan$ignored), " helper file(s):")
  for (f in file_scan$ignored) {
    message("[03_structure]  - ", basename(f))
  }
}

if (length(q_files_found) == 0) {
  message("[03_structure] No external STRUCTURE Q files found; checked folder: ", file_scan$root)
  message("[03_structure] Done.")
} else {
  message("[03_structure] Found ", length(q_files_found), " candidate STRUCTURE Q file(s).")
  
  parsed_by_k <- list()
  parsed_log <- list()
  
  for (f in q_files_found) {
    meta <- extract_k_run_from_filename(f)
    parsed <- read_q_file(f)
    
    if (!parsed$ok) {
      message("[03_structure] Skipping unreadable/non-tabular Q file: ", basename(f))
      parsed_log[[length(parsed_log) + 1]] <- data.frame(
        file = basename(f),
        full_path = f,
        K_from_filename = meta$K,
        run_from_filename = meta$run,
        K_detected = NA_integer_,
        rows = NA_integer_,
        cols = NA_integer_,
        reader = parsed$reader,
        status = "unreadable",
        note = parsed$reason,
        row_sum_min = NA_real_,
        row_sum_max = NA_real_,
        stringsAsFactors = FALSE
      )
      next
    }
    
    q_df <- parsed$data
    k_detected <- parsed$K_detected
    run_id <- meta$run
    k_file <- meta$K
    
    k_final <- if (!is.na(k_file)) k_file else k_detected
    if (!is.na(k_file) && k_file != k_detected) {
      warning("[03_structure] K mismatch in ", basename(f), ": filename K=", k_file, " but detected K=", k_detected, ". Using detected K.")
      k_final <- k_detected
    }
    
    if (!parsed$row_sums_ok) {
      warning("[03_structure] ", basename(f), ": row sums are not approximately 1 (tolerance 0.05).")
    }
    
    message(
      "[03_structure] Parsed: ", basename(f),
      " | rows=", parsed$n_rows,
      " cols=", parsed$n_cols,
      " | K=", k_final,
      " | run=", ifelse(is.na(run_id), "NA", run_id),
      " | reader=", parsed$reader
    )
    
    if (k_final == 1) {
      message("[03_structure] Skipping K=1 because single-cluster STRUCTURE output is not informative for admixture plots: ", basename(f))
      parsed_log[[length(parsed_log) + 1]] <- data.frame(
        file = basename(f),
        full_path = f,
        K_from_filename = k_file,
        run_from_filename = run_id,
        K_detected = k_detected,
        rows = parsed$n_rows,
        cols = parsed$n_cols,
        reader = parsed$reader,
        status = "skipped_k1",
        note = "valid K=1, intentionally skipped for admixture plotting",
        row_sum_min = parsed$row_sum_min,
        row_sum_max = parsed$row_sum_max,
        stringsAsFactors = FALSE
      )
      next
    }
    
    # Save per-run tables + plots
    run_site_summary <- write_site_mean_q(q_df, k = k_final, run = run_id)
    plot_individual_q(q_df, k = k_final, run = run_id)
    plot_site_mean_q(run_site_summary, k = k_final, run = run_id)
    
    k_name <- paste0("K", k_final)
    parsed_by_k[[k_name]] <- append(parsed_by_k[[k_name]], list(run_site_summary))
    
    parsed_log[[length(parsed_log) + 1]] <- data.frame(
      file = basename(f),
      full_path = f,
      K_from_filename = k_file,
      run_from_filename = run_id,
      K_detected = k_detected,
      rows = parsed$n_rows,
      cols = parsed$n_cols,
      reader = parsed$reader,
      status = "parsed",
      note = "ok",
      row_sum_min = parsed$row_sum_min,
      row_sum_max = parsed$row_sum_max,
      stringsAsFactors = FALSE
    )
  }
  
  # Write parsing audit table
  if (length(parsed_log) > 0) {
    parsed_summary <- bind_rows(parsed_log)
    summary_file <- file.path(TABLES_SUPP_DIR, "structure_parsed_files_summary.csv")
    write.csv(parsed_summary, summary_file, row.names = FALSE)
    message("[03_structure] Saved: ", summary_file)
  }
  
  # Aggregate across runs per K
  all_summaries <- list()
  for (k_name in names(parsed_by_k)) {
    run_tables <- parsed_by_k[[k_name]]
    if (length(run_tables) == 0) next
    
    if (length(run_tables) == 1) {
      agg <- run_tables[[1]]
    } else {
      agg <- bind_rows(run_tables, .id = "run") %>%
        group_by(Site) %>%
        summarise(across(starts_with("Q"), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
        arrange(Site)
    }
    
    all_summaries[[k_name]] <- agg
    
    k_num <- as.integer(gsub("[^0-9]", "", k_name))
    out_file <- file.path(TABLES_SUPP_DIR, sprintf("structure_meanQ_by_site_K%d.csv", k_num))
    write.csv(agg, out_file, row.names = FALSE)
    message("[03_structure] Saved aggregated-by-K summary: ", out_file)
    
    plot_site_mean_q(agg, k = k_num, run = NA_integer_)
  }
  
  # Dominant cluster summary across available K summaries
  if (length(all_summaries) > 0) {
    dom_list <- lapply(names(all_summaries), function(k_name) {
      dat <- all_summaries[[k_name]]
      q_cols <- grep("^Q", names(dat), value = TRUE)
      if (length(q_cols) == 0) return(NULL)
      
      dom_idx <- apply(dat[, q_cols, drop = FALSE], 1, which.max)
      dom_val <- apply(dat[, q_cols, drop = FALSE], 1, max)
      
      data.frame(
        K = gsub("[^0-9]", "", k_name),
        Site = dat$Site,
        Dominant_Cluster = paste0("Q", dom_idx),
        Dominant_Q = dom_val,
        stringsAsFactors = FALSE
      )
    })
    
    structure_site_cluster_summary <- bind_rows(dom_list)
    out_dom <- file.path(TABLES_SUPP_DIR, "structure_site_cluster_summary.csv")
    write.csv(structure_site_cluster_summary, out_dom, row.names = FALSE)
    message("[03_structure] Saved: ", out_dom)
  }
  
  message("[03_structure] Done.")
}