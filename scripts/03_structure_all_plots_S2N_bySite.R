# scripts/03_structure_all_plots_S2N_bySite.R
############################################################
# STRUCTURE helper outputs (final individual barplots only)
# - Reads externally generated STRUCTURE Q files
# - Builds ONE final individual plot per K (no per-run figures)
# - Builds ONE combined all-K figure
# - Individuals ordered SOUTH -> NORTH by site latitude
# Outputs:
# - outputs/figures/structure_individual_barplot_K{K}.jpeg
# - outputs/figures/structure_individual_barplot_allK.jpeg
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

clamp01 <- function(x, tol = 1e-6) {
  x[x < 0 & x >= -tol] <- 0
  x[x > 1 & x <= 1 + tol] <- 1
  x
}

fail_validation <- function(k, file_base, run_id, msg) {
  stop(
    paste0(
      "[03_structure] Validation failed for K=", k,
      " | file=", file_base,
      " | run=", ifelse(is.na(run_id), "NA", run_id),
      "\n", msg
    ),
    call. = FALSE
  )
}

# ----------------------------
# 1) Individual IDs and Site map
# ----------------------------
id_col_dfids <- resolve_col(df_ids, c("ind", "individual", "sample", "sampleid", "id"))
site_col_dfids <- resolve_col(df_ids, c("site", "population", "pop"))
if (is.na(id_col_dfids) || is.na(site_col_dfids)) {
  stop("[03_structure] df_ids must contain individual ID and Site columns.")
}

ids_dfids <- trimws(as.character(df_ids[[id_col_dfids]]))
ids_dfids <- ids_dfids[nzchar(ids_dfids)]
site_map_dfids <- setNames(as.character(df_ids[[site_col_dfids]]), ids_dfids)

load_ids_order_from_raw <- function() {
  ids_path <- file.path(PROJECT_ROOT, "data", "structure", "ids_order_from_raw.csv")
  if (!file.exists(ids_path)) {
    return(list(ids = character(0), source = "ids_order_from_raw.csv (missing)"))
  }
  
  ids_df <- tryCatch(
    read.csv(ids_path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL
  )
  if (is.null(ids_df) || nrow(ids_df) == 0) {
    return(list(ids = character(0), source = "ids_order_from_raw.csv (unreadable)"))
  }
  
  id_col <- resolve_col(ids_df, c("ind", "individual", "sample", "sampleid", "id", "ind_id"))
  if (is.na(id_col)) id_col <- names(ids_df)[1]
  
  ids <- trimws(as.character(ids_df[[id_col]]))
  ids <- ids[nzchar(ids)]
  ids <- unique(ids)
  
  list(ids = ids, source = "data/structure/ids_order_from_raw.csv")
}

extract_structure_order_from_results <- function() {
  res_dir <- file.path(PROJECT_ROOT, "data", "structure", "STRUCTURE_ZG_HEG-results")
  if (!dir.exists(res_dir)) {
    return(list(ids = character(0), pop = integer(0), source = "STRUCTURE results (missing)"))
  }
  
  files <- list.files(res_dir, full.names = TRUE)
  files <- files[file.info(files)$isdir %in% FALSE]
  if (length(files) == 0) {
    return(list(ids = character(0), pop = integer(0), source = "STRUCTURE results (empty)"))
  }
  
  # use first K>=2 file if available
  b <- basename(files)
  k_vals <- suppressWarnings(as.integer(gsub(".*K([0-9]+).*", "\\1", b, perl = TRUE)))
  ord <- order(ifelse(is.na(k_vals) | k_vals < 2, Inf, k_vals), b)
  
  for (i in ord) {
    f <- files[i]
    txt <- tryCatch(readLines(f, warn = FALSE), error = function(e) character(0))
    if (length(txt) == 0) next
    
    anc_start <- grep("^Inferred ancestry of individuals:", txt)
    if (length(anc_start) == 0) next
    block <- txt[(anc_start[1] + 1):length(txt)]
    
    keep <- grepl("^\\s*[0-9]+\\s+\\S+\\s+\\([^)]*\\)\\s+[0-9]+\\s*:\\s+", block)
    lines <- block[keep]
    if (length(lines) == 0) next
    
    labels <- sub("^\\s*[0-9]+\\s+(\\S+)\\s+\\([^)]*\\)\\s+[0-9]+\\s*:.*$", "\\1", lines, perl = TRUE)
    pops <- suppressWarnings(as.integer(sub("^\\s*[0-9]+\\s+\\S+\\s+\\([^)]*\\)\\s+([0-9]+)\\s*:.*$", "\\1", lines, perl = TRUE)))
    labels <- trimws(labels)
    labels <- labels[nzchar(labels)]
    
    if (length(labels) > 0) {
      return(list(ids = labels, pop = pops[seq_along(labels)], source = paste0("", f)))
    }
  }
  
  list(ids = character(0), pop = integer(0), source = "STRUCTURE results (no ancestry block parsed)")
}

ids_order_raw <- load_ids_order_from_raw()
ids_results <- extract_structure_order_from_results()

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
  
  if (nrow(out) == 0) stop("[03_structure] No valid site latitude information available.")
  setNames(out$Latitude, out$Site)
}

site_lat_map <- load_site_latitude(meta)

build_base_order <- function(ids_vec, site_map_final) {
  out <- data.frame(
    Individual = ids_vec,
    Site = site_map_final[ids_vec],
    stringsAsFactors = FALSE
  )
  out$SiteLat <- as.numeric(site_lat_map[out$Site])
  out <- out %>%
    mutate(
      Site = as.character(Site),
      SiteLat = ifelse(is.na(SiteLat), Inf, SiteLat)
    ) %>%
    arrange(SiteLat, Site, Individual)
  out$PlotIndex <- seq_len(nrow(out))
  out
}

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
  
  files <- list.files(q_root, recursive = TRUE, full.names = TRUE, pattern = "(?i)\\.(q|txt|csv)$")
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
  if (is.na(run)) run <- suppressWarnings(as.integer(sub("(?i).*[_-]run[_-]?([0-9]+).*$", "\\1", b, perl = TRUE)))
  
  list(K = K, run = run)
}

# ----------------------------
# 4) Robust Q parser
# ----------------------------
numeric_columns <- function(df) {
  as.data.frame(lapply(df, function(col) {
    x <- trimws(as.character(col))
    x[x == ""] <- NA_character_
    suppressWarnings(as.numeric(x))
  }), check.names = FALSE, stringsAsFactors = FALSE)
}

window_stats <- function(mat, tol = 1e-6) {
  na_count <- sum(is.na(mat))
  lo <- sum(mat < -tol, na.rm = TRUE)
  hi <- sum(mat > 1 + tol, na.rm = TRUE)
  valid_rows <- complete.cases(mat)
  if (any(valid_rows)) {
    rs <- rowSums(mat[valid_rows, , drop = FALSE])
    row_min <- min(rs)
    row_max <- max(rs)
    row_bad <- sum(abs(rs - 1) > 1e-3)
    row_dev <- mean(abs(rs - 1))
  } else {
    row_min <- NA_real_
    row_max <- NA_real_
    row_bad <- Inf
    row_dev <- Inf
  }
  list(na = na_count, lo = lo, hi = hi, row_min = row_min, row_max = row_max, row_bad = row_bad, row_dev = row_dev)
}

choose_q_columns <- function(num_df, k_hint = NA_integer_, tol = 1e-6) {
  m <- ncol(num_df)
  if (m < 2) return(NULL)
  
  get_candidates <- function(k_values) {
    out <- list()
    for (k in k_values) {
      if (k < 2 || k > m) next
      for (s in seq_len(m - k + 1)) {
        idx <- s:(s + k - 1)
        mat <- as.matrix(num_df[, idx, drop = FALSE])
        storage.mode(mat) <- "numeric"
        st <- window_stats(mat, tol = tol)
        score <- st$na * 1e7 + (st$lo + st$hi) * 1e6 + st$row_bad * 1e3 + st$row_dev * 100
        out[[length(out) + 1]] <- list(idx = idx, k = k, q = mat, stats = st, score = score)
      }
    }
    out
  }
  
  cands <- list()
  if (!is.na(k_hint)) cands <- get_candidates(k_hint)
  if (length(cands) == 0) cands <- get_candidates(2:m)
  if (length(cands) == 0) return(NULL)
  
  best_i <- which.min(vapply(cands, function(z) z$score, numeric(1)))
  cands[[best_i]]
}

inspect_q_file_lines <- function(path, expected_k = NA_integer_) {
  txt <- tryCatch(readLines(path, warn = FALSE), error = function(e) character(0))
  if (length(txt) == 0) {
    return(list(total = 0L, blank = 0L, numeric_like = 0L, wrong_token_count = 0L, duplicate_lines = 0L))
  }
  
  raw <- trimws(txt)
  nonblank <- raw[nzchar(raw)]
  blank_n <- sum(!nzchar(raw))
  
  token_n <- vapply(strsplit(nonblank, "\\s+"), length, integer(1))
  numeric_like <- vapply(nonblank, function(line) {
    toks <- strsplit(line, "\\s+")[[1]]
    all(!is.na(suppressWarnings(as.numeric(toks))))
  }, logical(1))
  
  if (is.na(expected_k)) {
    numeric_rows <- sum(numeric_like)
    wrong_tok <- NA_integer_
  } else {
    numeric_rows <- sum(numeric_like & token_n == expected_k)
    wrong_tok <- sum(numeric_like & token_n != expected_k)
  }
  
  list(
    total = length(txt),
    blank = blank_n,
    numeric_like = numeric_rows,
    wrong_token_count = wrong_tok,
    duplicate_lines = sum(duplicated(nonblank))
  )
}

read_q_matrix <- function(path, k_hint = NA_integer_, tol = 1e-6) {
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
    
    num_df <- numeric_columns(dat)
    if (nrow(num_df) == 0 || ncol(num_df) < 2) next
    
    keep_non_empty <- vapply(num_df, function(x) sum(!is.na(x)) > 0, logical(1))
    num_df <- num_df[, keep_non_empty, drop = FALSE]
    if (ncol(num_df) < 2) next
    
    first_col <- num_df[[1]]
    if (!anyNA(first_col) && length(first_col) > 1 && all(abs(first_col - seq_len(length(first_col))) < 1e-8)) {
      num_df <- num_df[, -1, drop = FALSE]
      if (ncol(num_df) < 2) next
    }
    
    sel <- choose_q_columns(num_df, k_hint = k_hint, tol = tol)
    if (is.null(sel)) next
    
    q <- sel$q
    q <- q[rowSums(!is.na(q)) > 0, , drop = FALSE]
    colnames(q) <- paste0("Q", seq_len(ncol(q)))
    
    line_diag <- inspect_q_file_lines(path, expected_k = ncol(q))
    
    return(list(ok = TRUE, matrix = q, reader = r$label, q_indices = sel$idx, stats = sel$stats, line_diag = line_diag))
  }
  
  list(ok = FALSE, matrix = NULL, reader = NA_character_, q_indices = integer(0), stats = NULL, line_diag = inspect_q_file_lines(path, expected_k = k_hint))
}

summarize_candidate <- function(file_base, k_infer, parsed, reason) {
  if (!parsed$ok) {
    reason_print <- if (grepl("^REJECTED:", reason)) reason else paste0("REJECTED: ", reason)
    message("[03_structure] Candidate: ", file_base,
            " | inferred K=", ifelse(is.na(k_infer), "NA", k_infer),
            " | ", reason_print)
    return(invisible(NULL))
  }
  
  q <- parsed$matrix
  st <- parsed$stats
  ld <- parsed$line_diag
  message("[03_structure] Candidate: ", file_base,
          " | inferred K=", ifelse(is.na(k_infer), "NA", k_infer),
          " | rows=", nrow(q),
          " | cols=", ncol(q),
          " | q_col_indices=", paste(parsed$q_indices, collapse = ","),
          " | numeric=TRUE",
          " | min=", signif(min(q, na.rm = TRUE), 6),
          " | max=", signif(max(q, na.rm = TRUE), 6),
          " | rowSumRange=", signif(st$row_min, 6), "-", signif(st$row_max, 6),
          " | lines(total/blank/numeric_like)=", ld$total, "/", ld$blank, "/", ld$numeric_like,
          " | duplicate_nonblank_rows=", ld$duplicate_lines,
          " | ", reason)
}

validate_selected_run <- function(k_num, run_info, q_df, q_cols, id_reference,
                                  tol = 1e-6, row_tol = 5e-3,
                                  hard_row_tol = 5e-2,
                                  renormalize_rows = TRUE) {
  q_mat <- as.matrix(q_df[, q_cols, drop = FALSE])
  storage.mode(q_mat) <- "numeric"
  
  na_count <- sum(is.na(q_mat))
  lt_count <- sum(q_mat < -tol, na.rm = TRUE)
  gt_count <- sum(q_mat > 1 + tol, na.rm = TRUE)
  
  q_mat <- clamp01(q_mat, tol = tol)
  row_sums <- rowSums(q_mat)
  abs_dev <- abs(row_sums - 1)
  max_abs_dev <- max(abs_dev)
  bad_row_sums <- sum(abs_dev > row_tol)
  hard_bad_row_sums <- sum(abs_dev > hard_row_tol)
  
  duplicated_ids <- sum(duplicated(q_df$Individual))
  missing_ids <- sum(!(id_reference %in% q_df$Individual))
  
  message(sprintf("[03_structure] K=%d validation details:", k_num))
  message("  selected file/run: ", run_info$file_base, " / ", ifelse(is.na(run_info$run), "NA", run_info$run))
  message("  number of individuals: ", nrow(q_df))
  message("  number of Q columns: ", length(q_cols))
  message("  Q column names: ", paste(q_cols, collapse = ", "))
  message("  all Q columns numeric: TRUE")
  message("  total NA in Q columns: ", na_count)
  message("  Q min/max: ", signif(min(q_mat, na.rm = TRUE), 6), " / ", signif(max(q_mat, na.rm = TRUE), 6))
  message("  row-sum range: ", signif(min(row_sums), 6), " - ", signif(max(row_sums), 6))
  message("  max |rowSum-1|: ", signif(max_abs_dev, 6))
  message("  rows with row sums not close to 1: ", bad_row_sums)
  message("  duplicated individual IDs after merging: ", duplicated_ids)
  message("  missing IDs after matching to reference order: ", missing_ids)
  
  if (na_count > 0) fail_validation(k_num, run_info$file_base, run_info$run, paste0("NA values in Q columns: ", na_count))
  if (lt_count > 0 || gt_count > 0) {
    fail_validation(k_num, run_info$file_base, run_info$run,
                    paste0("Out-of-range Q values detected. <0 count=", lt_count, ", >1 count=", gt_count))
  }
  if (duplicated_ids > 0 || missing_ids > 0) {
    fail_validation(k_num, run_info$file_base, run_info$run,
                    paste0("ID matching issue: duplicated_ids=", duplicated_ids, ", missing_ids=", missing_ids))
  }
  if (hard_bad_row_sums > 0) {
    fail_validation(k_num, run_info$file_base, run_info$run,
                    paste0("Row sums show true invalidity for ", hard_bad_row_sums,
                           " rows (hard tol=", hard_row_tol,
                           "; max |rowSum-1|=", signif(max_abs_dev, 6), ")."))
  }
  
  if (renormalize_rows) {
    renorm_rows <- sum(row_sums > 0 & abs_dev > 1e-12)
    if (renorm_rows > 0) {
      message("  renormalizing ", renorm_rows,
              " rows to sum exactly 1 (prevents stacked-bar clipping at y=1).")
      q_mat <- q_mat / row_sums
      q_mat <- clamp01(q_mat, tol = tol)
      row_sums_after <- rowSums(q_mat)
      message("  post-renormalization row-sum range: ",
              signif(min(row_sums_after), 6), " - ", signif(max(row_sums_after), 6))
    }
  }
  
  q_df[, q_cols] <- as.data.frame(q_mat, stringsAsFactors = FALSE)
  q_df
}

# ----------------------------
# 5) Parse all runs and choose ONE final run per K
# ----------------------------
scan <- find_q_files()
if (length(scan$ignored) > 0) message("[03_structure] Ignored helper files: ", length(scan$ignored))

if (length(scan$files) == 0) {
  message("[03_structure] No external STRUCTURE Q files found; checked: ", scan$root)
  message("[03_structure] Done.")
} else {
  parsed_candidates <- list()
  run_inventory <- list()
  
  for (f in scan$files) {
    meta_name <- extract_k_run(f)
    parsed <- read_q_matrix(f, k_hint = meta_name$K, tol = 1e-6)
    
    if (!parsed$ok) {
      summarize_candidate(basename(f), meta_name$K, parsed, "could not parse a numeric Q matrix")
      run_inventory[[length(run_inventory) + 1]] <- data.frame(
        file = basename(f), full_path = f, K = meta_name$K, run = meta_name$run,
        reader = parsed$reader, n_rows = NA_integer_, n_cols = NA_integer_,
        q_col_indices = NA_character_, mean_abs_row_sum_deviation = NA_real_,
        status = "unreadable", reason = "could_not_parse_numeric_q_matrix",
        stringsAsFactors = FALSE
      )
      next
    }
    
    q_mat <- parsed$matrix
    K_detected <- ncol(q_mat)
    K_final <- if (!is.na(meta_name$K)) meta_name$K else K_detected
    if (!is.na(meta_name$K) && meta_name$K != K_detected) K_final <- K_detected
    
    row_dev <- mean(abs(rowSums(clamp01(q_mat, 1e-6), na.rm = TRUE) - 1), na.rm = TRUE)
    
    reject_reason <- NULL
    status <- "parsed_candidate"
    if (K_final == 1) {
      reject_reason <- "K1 skipped"
      status <- "skipped_k1"
    }
    
    ld <- parsed$line_diag
    cause_msg <- "ACCEPTED as parsed candidate"
    if (ld$blank > 0 && ld$numeric_like == nrow(q_mat)) {
      cause_msg <- paste0(cause_msg, " (blank lines ignored; no extra non-data rows kept)")
    }
    if (!is.na(ld$wrong_token_count) && ld$wrong_token_count > 0) {
      cause_msg <- paste0(cause_msg, " (non-Q token-count rows excluded)")
    }
    
    summarize_candidate(
      basename(f),
      meta_name$K,
      parsed,
      ifelse(is.null(reject_reason), cause_msg, paste0("REJECTED: ", reject_reason))
    )
    
    run_inventory[[length(run_inventory) + 1]] <- data.frame(
      file = basename(f), full_path = f, K = K_final, run = meta_name$run,
      reader = parsed$reader, n_rows = nrow(q_mat), n_cols = ncol(q_mat),
      q_col_indices = paste(parsed$q_indices, collapse = ";"),
      mean_abs_row_sum_deviation = row_dev,
      status = status,
      reason = ifelse(is.null(reject_reason), "candidate", reject_reason),
      stringsAsFactors = FALSE
    )
    
    if (is.null(reject_reason)) {
      parsed_candidates[[length(parsed_candidates) + 1]] <- list(
        file = f,
        file_base = basename(f),
        K = K_final,
        run = meta_name$run,
        reader = parsed$reader,
        q = q_mat,
        mad = row_dev
      )
    }
  }
  
  if (length(run_inventory) > 0) {
    inv_df <- bind_rows(run_inventory) %>% arrange(K, run, file)
    inv_file <- file.path(TABLES_SUPP_DIR, "structure_run_inventory.csv")
    write.csv(inv_df, inv_file, row.names = FALSE)
    message("[03_structure] Saved run inventory: ", inv_file)
  }
  
  if (length(parsed_candidates) == 0) {
    message("[03_structure] No usable K>=2 Q runs after parsing.")
    message("[03_structure] Done.")
  } else {
    parsed_runs <- Filter(function(x) x$K >= 2, parsed_candidates)
    row_counts <- sort(unique(vapply(parsed_runs, function(x) nrow(x$q), integer(1))))
    if (length(row_counts) != 1) {
      stop("[03_structure] Parsed Q files have inconsistent row counts across runs: ", paste(row_counts, collapse = ", "))
    }
    q_n <- row_counts[1]
    
    # choose reference source with exact match to parsed Q row count
    ref_candidates <- list(
      list(name = ids_order_raw$source, ids = ids_order_raw$ids),
      list(name = ids_results$source, ids = unique(ids_results$ids)),
      list(name = "df_ids object order", ids = unique(ids_dfids))
    )
    
    ref_lengths <- vapply(ref_candidates, function(x) length(x$ids), integer(1))
    message("[03_structure] Reference lengths (for matching parsed n=", q_n, "): ",
            paste0(vapply(ref_candidates, function(x) x$name, character(1)), "=", ref_lengths, collapse = " | "))
    
    idx_match <- which(ref_lengths == q_n)
    if (length(idx_match) == 0) {
      stop("[03_structure] No reference ID source matches parsed Q row count n=", q_n,
           ". Check ids_order_from_raw.csv and STRUCTURE input export order.")
    }
    
    # preference order already encoded above: raw IDs -> results labels -> df_ids
    ref_pick <- idx_match[1]
    id_reference <- ref_candidates[[ref_pick]]$ids
    ref_source_used <- ref_candidates[[ref_pick]]$name
    message("[03_structure] Using reference ID source: ", ref_source_used, " (n=", length(id_reference), ")")
    
    # diagnose mismatch vs df_ids
    missing_in_df <- setdiff(id_reference, ids_dfids)
    extra_in_df <- setdiff(ids_dfids, id_reference)
    if (length(missing_in_df) > 0 || length(extra_in_df) > 0) {
      message("[03_structure] ID mismatch diagnostic:")
      message("  IDs in reference but not df_ids: ", length(missing_in_df),
              ifelse(length(missing_in_df) > 0, paste0(" (e.g., ", paste(head(missing_in_df, 5), collapse = ", "), ")"), ""))
      message("  IDs in df_ids but not reference: ", length(extra_in_df),
              ifelse(length(extra_in_df) > 0, paste0(" (e.g., ", paste(head(extra_in_df, 5), collapse = ", "), ")"), ""))
    }
    
    # build best-available site map for all reference IDs
    site_map_final <- site_map_dfids
    missing_site_before <- sum(is.na(site_map_final[id_reference]))
    
    if (length(ids_results$ids) > 0 && length(ids_results$pop) == length(ids_results$ids)) {
      lab_pop <- data.frame(Individual = ids_results$ids, PopIdx = ids_results$pop, stringsAsFactors = FALSE)
      known <- lab_pop %>%
        mutate(Site = site_map_dfids[Individual]) %>%
        filter(!is.na(Site), !is.na(PopIdx))
      
      if (nrow(known) > 0) {
        pop_to_site <- known %>%
          count(PopIdx, Site, name = "n") %>%
          group_by(PopIdx) %>%
          arrange(desc(n), Site, .by_group = TRUE) %>%
          slice(1) %>%
          ungroup() %>%
          select(PopIdx, Site)
        
        inferred <- lab_pop %>%
          left_join(pop_to_site, by = "PopIdx")
        inferred_map <- setNames(inferred$Site, inferred$Individual)
        
        need_fill <- id_reference[is.na(site_map_final[id_reference])]
        if (length(need_fill) > 0) {
          fill_vals <- inferred_map[need_fill]
          site_map_final[need_fill] <- fill_vals
        }
      }
    }
    
    missing_site_after <- sum(is.na(site_map_final[id_reference]))
    recovered_site_n <- missing_site_before - missing_site_after
    if (missing_site_before > 0) {
      message("[03_structure] Site metadata recovery diagnostic:")
      message("  missing Site labels before recovery: ", missing_site_before)
      message("  recovered via STRUCTURE results metadata: ", recovered_site_n)
      message("  still unmatched after recovery: ", missing_site_after,
              " (these are kept as unmatched and plotted in 'Unknown').")
    }
    
    base_order_df <- build_base_order(id_reference, site_map_final)
    if (nrow(base_order_df) != q_n) {
      stop("[03_structure] Internal alignment error: reference order n=", nrow(base_order_df), " but parsed Q n=", q_n)
    }
    
    missing_site_n <- sum(is.na(base_order_df$Site))
    if (missing_site_n > 0) {
      warning("[03_structure] ", missing_site_n,
              " individuals have missing Site labels after site-map recovery; they will be plotted at end.")
    }
    
    runs_by_k <- split(parsed_runs, sapply(parsed_runs, function(x) x$K))
    selected_rows <- list()
    allk_plot_data <- list()
    validations_passed <- 0L
    
    site_blocks <- base_order_df %>%
      mutate(Site = ifelse(is.na(Site), "Unknown", Site)) %>%
      count(Site, name = "n") %>%
      mutate(xmax = cumsum(n), xmin = xmax - n + 1, xmid = (xmin + xmax) / 2)
    separators <- site_blocks$xmax[-nrow(site_blocks)] + 0.5
    
    k_levels <- sort(unique(as.integer(names(runs_by_k))))
    
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
      })) %>% arrange(mad, run, file_base)
      
      best <- k_runs[[match(score_df$file[1], sapply(k_runs, function(x) x$file))]]
      
      message("[03_structure] Found ", length(k_runs), " usable run(s) for K=", k_num)
      message("[03_structure] Using best run for K=", k_num,
              " (run=", ifelse(is.na(best$run), "NA", best$run),
              ", reader=", best$reader,
              ", mean|rowSum-1|=", signif(best$mad, 4), ")")
      
      q_df <- cbind(base_order_df, as.data.frame(best$q, stringsAsFactors = FALSE, check.names = FALSE))
      q_cols <- grep("^Q", names(q_df), value = TRUE)
      
      q_df <- validate_selected_run(
        k_num = k_num,
        run_info = best,
        q_df = q_df,
        q_cols = q_cols,
        id_reference = base_order_df$Individual,
        tol = 1e-6,
        row_tol = 5e-3,
        hard_row_tol = 5e-2,
        renormalize_rows = TRUE
      )
      
      plot_df <- q_df %>%
        mutate(Site = ifelse(is.na(Site), "Unknown", Site)) %>%
        select(PlotIndex, Site, all_of(q_cols)) %>%
        pivot_longer(cols = all_of(q_cols), names_to = "Cluster", values_to = "Q")
      
      if (!is.numeric(plot_df$Q)) fail_validation(k_num, best$file_base, best$run, "Long-format Q column is not numeric after pivot.")
      if (anyNA(plot_df$Q)) fail_validation(k_num, best$file_base, best$run, paste0("Long-format Q contains NA: ", sum(is.na(plot_df$Q))))
      if (any(plot_df$Q < -1e-6 | plot_df$Q > 1 + 1e-6)) {
        fail_validation(k_num, best$file_base, best$run, "Long-format Q contains values outside [0,1].")
      }
      
      plot_df$Q <- clamp01(plot_df$Q, tol = 1e-6)
      plot_df$K <- k_num
      plot_df$KLabel <- factor(paste0("K=", k_num), levels = paste0("K=", k_levels))
      allk_plot_data[[length(allk_plot_data) + 1]] <- plot_df
      
      p <- ggplot(plot_df, aes(x = PlotIndex, y = Q, fill = Cluster)) +
        geom_col(width = 1) +
        geom_vline(xintercept = separators, linewidth = 0.25, color = "grey25") +
        scale_x_continuous(breaks = site_blocks$xmid, labels = site_blocks$Site, expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        coord_cartesian(ylim = c(0, 1), expand = FALSE, clip = "on") +
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
      
      validations_passed <- validations_passed + 1L
    }
    
    if (length(allk_plot_data) > 0) {
      combined_plot_df <- bind_rows(allk_plot_data)
      
      p_all <- ggplot(combined_plot_df, aes(x = PlotIndex, y = Q, fill = Cluster)) +
        geom_col(width = 1) +
        geom_vline(xintercept = separators, linewidth = 0.2, color = "grey25") +
        scale_x_continuous(breaks = site_blocks$xmid, labels = site_blocks$Site, expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        coord_cartesian(ylim = c(0, 1), expand = FALSE, clip = "on") +
        facet_grid(rows = vars(KLabel)) +
        theme_bw(base_size = 11) +
        theme(
          panel.grid = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
          axis.ticks.x = element_blank(),
          legend.position = "right",
          strip.text.y = element_text(face = "bold")
        ) +
        labs(
          title = "STRUCTURE ancestry barplots across K",
          subtitle = "Individuals ordered south to north by site latitude",
          x = "Site blocks (south -> north)",
          y = "Ancestry proportion"
        )
      
      out_all <- file.path(FIGURES_DIR, "structure_individual_barplot_allK.jpeg")
      ggsave(out_all, p_all, width = 13, height = max(5.5, 2.2 * length(unique(combined_plot_df$K))), dpi = 320)
      message("[03_structure] Saved: ", out_all)
    }
    
    if (length(selected_rows) > 0) {
      selected_df <- bind_rows(selected_rows) %>% arrange(K)
      selected_file <- file.path(TABLES_SUPP_DIR, "structure_selected_runs.csv")
      write.csv(selected_df, selected_file, row.names = FALSE)
      message("[03_structure] Saved selected-run summary: ", selected_file)
    }
    
    message("[03_structure] Validation summary: all K plots passed validation (", validations_passed, " K values).")
    message("[03_structure] Validation summary: no rows dropped in plotting.")
    message("[03_structure] Validation summary: combined figure built successfully.")
    message("[03_structure] Done.")
  }
}