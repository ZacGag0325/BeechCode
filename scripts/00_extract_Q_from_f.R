############################################################
# 00_extract_Q_from_f.R  (BeechCode)
# Batch extract Q matrices from STRUCTURE *_f files (mac-safe)
# Input:  configured STRUCTURE run folder, e.g. data/structure/HEG_ZG_run2/**/_f
# Output: data/structure/Q_Files_2/*.Q  + summary CSV
############################################################

### 0) PROJECT ROOT + PATHS (ROBUST) ----
suppressPackageStartupMessages(library(here))

FALLBACK_ROOT <- file.path(path.expand("~"), "Desktop", "BeechCode")

normalize_root_candidate <- function(path) {
  path <- normalizePath(path, mustWork = FALSE)
  # If R/here accidentally gives the .Rproj path as the root, step back to the
  # containing project directory before appending data/ or outputs/ paths.
  if (grepl("\\.Rproj$", basename(path), ignore.case = TRUE)) {
    path <- dirname(path)
  }
  path
}

find_project_root_for_extractor <- function() {
  source_file <- tryCatch(normalizePath(sys.frames()[[1]]$ofile, mustWork = FALSE), error = function(e) NA_character_)
  candidates <- c(
    here::here(),
    getwd(),
    if (!is.na(source_file)) dirname(dirname(source_file)) else character(0),
    FALLBACK_ROOT
  )
  candidates <- unique(normalize_root_candidate(candidates))
  
  for (candidate in candidates) {
    if (dir.exists(file.path(candidate, "data", "raw")) ||
        file.exists(file.path(candidate, "scripts", "00_extract_Q_from_f.R")) ||
        file.exists(file.path(candidate, "scripts", "00_extracted_Q_from_f.R"))) {
      return(candidate)
    }
  }
  
  stop(
    "Can't find project root.\n",
    "Expected a folder containing data/raw or the BeechCode scripts directory.\n\n",
    "Tried:\n  - ", paste(candidates, collapse = "\n  - "), "\n\n",
    "Fix: open your BeechCode.Rproj, source this script from inside the project, ",
    "or change FALLBACK_ROOT to your actual folder.",
    call. = FALSE
  )
}

PROJECT_ROOT <- find_project_root_for_extractor()
setwd(PROJECT_ROOT)

OUTPUTS_DIR <- file.path(PROJECT_ROOT, "outputs")
RUN_TAG <- "v1"
RUN_OUT <- file.path(OUTPUTS_DIR, RUN_TAG)

STRUCTURE_RUN_NAME <- "HEG_ZG_run2"
Q_FOLDER_NAME <- "Q_Files_2"

STRUCTURE_RUNS_DIR <- file.path(RUN_OUT, "structure_runs")

find_structure_results_dir <- function(structure_run_name = STRUCTURE_RUN_NAME) {
  candidates <- c(
    file.path(PROJECT_ROOT, "data", "structure", structure_run_name),
    file.path(PROJECT_ROOT, "outputs", "v1", "structure_runs", structure_run_name)
  )
  hit <- candidates[dir.exists(candidates)]
  if (length(hit) == 0) {
    stop(
      "Could not find configured STRUCTURE results folder for run '", structure_run_name, "'.\n",
      "Looked for:\n  - ", paste(candidates, collapse = "\n  - "),
      call. = FALSE
    )
  }
  hit[1]
}

message("Project root: ", PROJECT_ROOT)
message("Configured STRUCTURE_RUN_NAME: ", STRUCTURE_RUN_NAME)
message("Configured Q_FOLDER_NAME: ", Q_FOLDER_NAME)

RESULTS_DIR <- find_structure_results_dir()
message("Using raw STRUCTURE folder: ", RESULTS_DIR)

# Where extracted Q files will go
Q_OUT_DIR <- file.path(PROJECT_ROOT, "data", "structure", Q_FOLDER_NAME)
dir.create(Q_OUT_DIR, showWarnings = FALSE, recursive = TRUE)
message("Q files will be written to: ", Q_OUT_DIR)

### 1) Helper functions ----

# Parse a numeric line -> numeric vector (or NULL)
parse_num_line <- function(ln) {
  ln <- trimws(ln)
  if (ln == "") return(NULL)
  # Allow: digits, spaces, decimal points, minus, scientific notation e/E
  if (!grepl("^[-0-9eE\\.\\s]+$", ln)) return(NULL)
  
  parts <- strsplit(ln, "\\s+")[[1]]
  suppressWarnings(v <- as.numeric(parts))
  if (length(v) < 2) return(NULL)
  if (any(!is.finite(v))) return(NULL)
  v
}

# Extract the best candidate Q block from a STRUCTURE *_f file
extract_Q_from_f <- function(f) {
  lines <- readLines(f, warn = FALSE)
  
  blocks <- list()
  cur <- list()
  
  flush_block <- function() {
    if (length(cur) > 0) {
      blocks[[length(blocks) + 1]] <<- cur
      cur <<- list()
    }
  }
  
  for (ln in lines) {
    v <- parse_num_line(ln)
    if (is.null(v)) {
      flush_block()
    } else {
      cur[[length(cur) + 1]] <- v
    }
  }
  flush_block()
  
  if (length(blocks) == 0) return(NULL)
  
  score_block <- function(b) {
    maxc <- max(vapply(b, length, integer(1)))
    M <- matrix(NA_real_, nrow = length(b), ncol = maxc)
    for (i in seq_along(b)) M[i, seq_along(b[[i]])] <- b[[i]]
    
    # drop empty columns
    keep <- colSums(!is.na(M)) > 0
    M <- M[, keep, drop = FALSE]
    
    if (nrow(M) < 10 || ncol(M) < 2) return(list(score = -Inf, M = NULL))
    if (min(M, na.rm = TRUE) < -1e-6) return(list(score = -Inf, M = NULL))
    if (max(M, na.rm = TRUE) >  1 + 1e-6) return(list(score = -Inf, M = NULL))
    
    rs <- rowSums(M, na.rm = TRUE)
    prop_good <- mean(abs(rs - 1) < 0.05)
    
    # Require at least 70% rows summing ~1, then score by quality + size
    if (!is.finite(prop_good) || prop_good < 0.70) return(list(score = -Inf, M = NULL))
    
    score <- (prop_good * 1000) + nrow(M) + (ncol(M) / 10)
    list(score = score, M = M)
  }
  
  scored <- lapply(blocks, score_block)
  scores <- vapply(scored, `[[`, numeric(1), "score")
  
  if (all(!is.finite(scores))) return(NULL)
  
  best <- scored[[which.max(scores)]]$M
  
  # final drop all-NA columns (safety)
  keep <- colSums(!is.na(best)) > 0
  best <- best[, keep, drop = FALSE]
  best
}

# Try to guess K from filename (optional; not required)
guess_K_from_name <- function(f) {
  m <- regmatches(f, regexpr("K\\d+", f, ignore.case = TRUE))
  if (length(m) == 0 || is.na(m)) return(NA_integer_)
  as.integer(gsub("[^0-9]", "", m))
}

### 2) Find *_f files ----
f_files <- list.files(RESULTS_DIR, recursive = TRUE, full.names = TRUE, pattern = "_f$")
message("Found ", length(f_files), " *_f files in raw STRUCTURE folder")

if (length(f_files) == 0) {
  stop("No files ending in '_f' found under:\n  ", RESULTS_DIR)
}

### 3) Extract Q from each file ----
results <- data.frame(
  f_file = character(0),
  q_file = character(0),
  K_detected = integer(0),
  n_inds = integer(0),
  K_from_name = integer(0),
  ok = logical(0),
  stringsAsFactors = FALSE
)

ok_count <- 0
fail_count <- 0

for (f in f_files) {
  Q <- tryCatch(extract_Q_from_f(f), error = function(e) NULL)
  
  if (is.null(Q)) {
    fail_count <- fail_count + 1
    results <- rbind(results, data.frame(
      f_file = f,
      q_file = NA_character_,
      K_detected = NA_integer_,
      n_inds = NA_integer_,
      K_from_name = guess_K_from_name(f),
      ok = FALSE
    ))
    next
  }
  
  Kdet <- ncol(Q)
  nind <- nrow(Q)
  
  # output name: keep base filename, add .Q
  base <- basename(f)
  q_out <- file.path(Q_OUT_DIR, paste0(base, ".Q"))
  
  write.table(Q, file = q_out, quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  ok_count <- ok_count + 1
  results <- rbind(results, data.frame(
    f_file = f,
    q_file = q_out,
    K_detected = Kdet,
    n_inds = nind,
    K_from_name = guess_K_from_name(f),
    ok = TRUE
  ))
}

message("Done extracting. OK: ", ok_count, " | Failed: ", fail_count)

# Save summary
summary_csv <- file.path(Q_OUT_DIR, "extracted_Q_summary.csv")
write.csv(results, summary_csv, row.names = FALSE)
message("Saved summary CSV: ", summary_csv)

# Print counts per detected K
if (ok_count > 0) {
  tab <- sort(table(results$K_detected[results$ok]))
  message("Counts per K (detected from Q matrices):")
  print(tab)
  row_tab <- sort(table(results$n_inds[results$ok]))
  message("Detected Q row counts across successful files:")
  print(row_tab)
}

message("All extracted Q files are in: ", Q_OUT_DIR)