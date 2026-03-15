# scripts/03_structure_all_plots_S2N_bySite.R
############################################################
# STRUCTURE helper outputs only (NO STRUCTURE reruns)
# - Uses externally generated Q tables
# - Produces mean Q by site for each available K
# - Produces site-level dominant cluster summary
# Outputs (examples):
# - outputs/tables/supplementary/structure_meanQ_by_site_K2.csv
# - outputs/tables/supplementary/structure_meanQ_by_site_K3.csv
# - outputs/tables/supplementary/structure_meanQ_by_site_K4.csv
# - outputs/tables/supplementary/structure_site_cluster_summary.csv
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
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

read_q_file <- function(path) {
  q <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  id_col <- resolve_col(q, c("ind", "individual", "sample", "sampleid", "id"))
  if (is.na(id_col)) return(NULL)
  
  q_cols <- names(q)[sapply(q, is.numeric)]
  q_cols <- setdiff(q_cols, id_col)
  if (length(q_cols) < 2) return(NULL)
  
  # Normalize names to Q1..QK for readability
  q_sub <- q[, c(id_col, q_cols), drop = FALSE]
  names(q_sub) <- c("Individual", paste0("Q", seq_along(q_cols)))
  q_sub$Site <- site_map[as.character(q_sub$Individual)]
  q_sub
}

extract_k_from_filename <- function(path, k_default) {
  b <- basename(path)
  hit <- regmatches(b, regexpr("K[0-9]+", b, ignore.case = TRUE))
  if (length(hit) == 0 || !nzchar(hit)) return(k_default)
  as.integer(gsub("[^0-9]", "", hit))
}

write_site_mean_q <- function(q_df, k) {
  out <- q_df %>%
    filter(!is.na(Site)) %>%
    group_by(Site) %>%
    summarise(across(starts_with("Q"), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
    arrange(Site)
  
  out_file <- file.path(TABLES_SUPP_DIR, sprintf("structure_meanQ_by_site_K%d.csv", k))
  write.csv(out, out_file, row.names = FALSE)
  message("[03_structure] Saved: ", out_file)
  out
}

# Candidate single-table files and multi-K directory
q_candidates <- c(
  file.path(PROJECT_ROOT, "inputs", "structure_Q.csv"),
  file.path(PROJECT_ROOT, "inputs", "Q_matrix.csv"),
  file.path(PROJECT_ROOT, "inputs", "structure", "structure_Q.csv")
)

structure_dir <- file.path(PROJECT_ROOT, "inputs", "structure")
q_multi_files <- if (dir.exists(structure_dir)) {
  list.files(structure_dir, pattern = "(?i)(^Q_?K[0-9]+.*\\.csv$|structure.*K[0-9]+.*\\.csv$)", full.names = TRUE)
} else {
  character(0)
}

all_summaries <- list()

if (length(q_multi_files) > 0) {
  message("[03_structure] Detected multi-K STRUCTURE files: ", length(q_multi_files))
  for (f in q_multi_files) {
    q_df <- read_q_file(f)
    if (is.null(q_df)) {
      message("[03_structure] Skipping unreadable Q file: ", f)
      next
    }
    k <- extract_k_from_filename(f, k_default = ncol(q_df) - 2)
    all_summaries[[paste0("K", k)]] <- write_site_mean_q(q_df, k)
  }
} else {
  # Fallback to one Q file if available
  q_file <- q_candidates[file.exists(q_candidates)][1]
  if (is.na(q_file) || length(q_file) == 0) {
    message("[03_structure] No external STRUCTURE Q files found; skipping STRUCTURE helper outputs.")
  } else {
    q_df <- read_q_file(q_file)
    if (is.null(q_df)) {
      stop("[03_structure] Could not parse external Q file: ", q_file)
    }
    k <- ncol(q_df) - 2
    all_summaries[[paste0("K", k)]] <- write_site_mean_q(q_df, k)
  }
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