# scripts/00_run_all.R
############################################################
# Master runner: BeechCode genetics pipeline
# - Executes scripts in fixed order
# - Stops on first error
# - Anchors paths at project root
# - Runs each script in a clean R session (--vanilla)
# - Audits warnings emitted during the run without suppressing them
############################################################

find_project_root <- function() {
  candidates <- c(getwd(), normalizePath(file.path(getwd(), ".."), mustWork = FALSE))
  cmd_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
  if (length(cmd_file) > 0 && nzchar(cmd_file[1])) {
    candidates <- c(candidates, dirname(normalizePath(cmd_file[1], mustWork = FALSE)))
  }
  for (start in unique(candidates)) {
    cur <- normalizePath(start, mustWork = FALSE)
    repeat {
      if (file.exists(file.path(cur, "scripts", "00_run_all.R"))) return(cur)
      parent <- dirname(cur)
      if (identical(parent, cur)) break
      cur <- parent
    }
  }
  stop("[00_run_all] Cannot find project root containing scripts/00_run_all.R")
}

PROJECT_ROOT <- find_project_root()
setwd(PROJECT_ROOT)

cat("[00_run_all] Project root:", PROJECT_ROOT, "\n")

validate_file_exists <- function(path, context) {
  if (!file.exists(path)) stop("[00_run_all] Missing expected file after ", context, ": ", path)
}

classify_warning_text <- function(x) {
  x <- tolower(trimws(paste(x, collapse = " ")))
  if (!nzchar(x)) return("other")
  if (grepl("duplicate|duplicated|name repair|repaired names|make.unique", x)) return("duplicate_names")
  if (grepl("coerc|n[ao] introduced by coercion|converted|parsing", x)) return("coercion_or_parsing")
  if (grepl("read|import|column|sheet|csv|excel|metadata|workbook", x)) return("import_or_format")
  if (grepl("ellipse|too few|group has|centroid", x)) return("small_group_or_plotting")
  "other"
}

extract_warning_blocks <- function(lines) {
  if (length(lines) == 0) {
    return(data.frame(
      warning_index = integer(0),
      line_start = integer(0),
      line_end = integer(0),
      warning_text = character(0),
      warning_category = character(0),
      stringsAsFactors = FALSE
    ))
  }
  
  blocks <- list()
  i <- 1L
  k <- 1L
  while (i <= length(lines)) {
    this_line <- lines[i]
    is_warning_start <- grepl("^(Warning|warning)", this_line) || grepl("^Warning message", this_line)
    if (!is_warning_start) {
      i <- i + 1L
      next
    }
    
    start_i <- i
    block_lines <- this_line
    i <- i + 1L
    while (i <= length(lines)) {
      next_line <- lines[i]
      if (!nzchar(trimws(next_line))) {
        break
      }
      if (grepl("^--- ", next_line) || grepl("^\\[[0-9]{2}_", next_line) || grepl("^(Error|Execution halted)", next_line)) {
        break
      }
      if (grepl("^(Warning|warning)", next_line) && !grepl("^\\s", next_line)) {
        break
      }
      block_lines <- c(block_lines, next_line)
      i <- i + 1L
    }
    
    blocks[[k]] <- data.frame(
      warning_index = k,
      line_start = start_i,
      line_end = start_i + length(block_lines) - 1L,
      warning_text = paste(trimws(block_lines), collapse = " | "),
      warning_category = classify_warning_text(block_lines),
      stringsAsFactors = FALSE
    )
    k <- k + 1L
  }
  
  if (length(blocks) == 0) {
    return(data.frame(
      warning_index = integer(0),
      line_start = integer(0),
      line_end = integer(0),
      warning_text = character(0),
      warning_category = character(0),
      stringsAsFactors = FALSE
    ))
  }
  
  do.call(rbind, blocks)
}

run_script <- function(script_name, log_dir) {
  path <- file.path(PROJECT_ROOT, "scripts", script_name)
  if (!file.exists(path)) stop("[00_run_all] Missing script: ", path)
  
  log_file <- file.path(log_dir, paste0(sub("\\.R$", "", script_name), ".log"))
  if (file.exists(log_file)) file.remove(log_file)
  
  cat("\n--- Running in clean session:", script_name, "\n")
  res <- system2(
    command = "Rscript",
    args = c("--vanilla", path),
    stdout = log_file,
    stderr = log_file
  )
  
  log_lines <- if (file.exists(log_file)) readLines(log_file, warn = FALSE) else character(0)
  if (length(log_lines) > 0) {
    cat(paste0(log_lines, collapse = "\n"), "\n")
  }
  
  warning_blocks <- extract_warning_blocks(log_lines)
  if (nrow(warning_blocks) > 0) {
    warning_blocks$script <- script_name
    warning_blocks$log_file <- log_file
  }
  
  if (!identical(res, 0L)) {
    stop("[00_run_all] Script failed: ", script_name, " (exit code ", res, ")")
  }
  cat("--- Completed:", script_name, "\n")
  
  list(log_file = log_file, warnings = warning_blocks)
}

for (d in c(
  file.path(PROJECT_ROOT, "outputs"),
  file.path(PROJECT_ROOT, "outputs", "tables"),
  file.path(PROJECT_ROOT, "outputs", "figures"),
  file.path(PROJECT_ROOT, "outputs", "matrices"),
  file.path(PROJECT_ROOT, "outputs", "tables", "supplementary"),
  file.path(PROJECT_ROOT, "outputs", "figures", "supplementary")
)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

warning_log_dir <- file.path(PROJECT_ROOT, "outputs", "tables", "supplementary", "pipeline_logs")
dir.create(warning_log_dir, recursive = TRUE, showWarnings = FALSE)
warning_audit_file <- file.path(PROJECT_ROOT, "outputs", "tables", "supplementary", "pipeline_warning_audit.csv")
warning_summary_file <- file.path(PROJECT_ROOT, "outputs", "tables", "supplementary", "pipeline_warning_summary.csv")

all_warning_blocks <- list()
run_index <- 1L

master_run <- run_script("00_master_pipeline.R", warning_log_dir)
all_warning_blocks[[run_index]] <- master_run$warnings
run_index <- run_index + 1L

obj_dir <- file.path(PROJECT_ROOT, "outputs", "v1", "objects")
required_obj <- c("gi.rds", "gi_mll.rds", "df_ids.rds", "meta.rds")
missing_obj <- required_obj[!file.exists(file.path(obj_dir, required_obj))]
if (length(missing_obj) > 0) {
  stop("[00_run_all] Object build incomplete. Missing: ", paste(missing_obj, collapse = ", "))
}
cat("[00_run_all] Verified canonical objects in", obj_dir, "\n")
validate_file_exists(file.path(obj_dir, "df_ids.rds"), "00_master_pipeline.R")

for (script_name in c(
  "01_clonality.R",
  "07_allelic_richness.R",
  "02_hwe.R",
  "05_pca_dapc.R",
  "04_amova.R",
  "06_distance_matrices.R",
  "11_isolation_by_distance.R",
  "03_structure_all_plots_S2N_bySite.R"
)) {
  run_out <- run_script(script_name, warning_log_dir)
  all_warning_blocks[[run_index]] <- run_out$warnings
  run_index <- run_index + 1L
  
  if (identical(script_name, "01_clonality.R")) {
    validate_file_exists(file.path(PROJECT_ROOT, "outputs", "tables", "clonality_summary.csv"), "01_clonality.R")
  }
}

warning_audit <- do.call(
  rbind,
  c(
    Filter(Negate(is.null), all_warning_blocks),
    list(data.frame(
      warning_index = integer(0),
      line_start = integer(0),
      line_end = integer(0),
      warning_text = character(0),
      warning_category = character(0),
      script = character(0),
      log_file = character(0),
      stringsAsFactors = FALSE
    ))
  )
)

write.csv(warning_audit, warning_audit_file, row.names = FALSE)

warning_summary <- if (nrow(warning_audit) > 0) {
  aggregate(
    list(n_warnings = warning_audit$warning_text),
    by = list(script = warning_audit$script, warning_category = warning_audit$warning_category),
    FUN = length
  )
} else {
  data.frame(
    script = character(0),
    warning_category = character(0),
    n_warnings = integer(0),
    stringsAsFactors = FALSE
  )
}
write.csv(warning_summary, warning_summary_file, row.names = FALSE)

cat("\n[00_run_all] Pipeline completed successfully.\n")
cat(
  "[00_run_all] Warning audit summary:",
  nrow(warning_audit),
  "warning block(s) logged in",
  warning_audit_file,
  "with per-script counts in",
  warning_summary_file,
  "and full logs in",
  warning_log_dir,
  "\n"
)