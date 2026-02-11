############################################################
# 00_extract_Q_from_f.R  (BeechCode)
# Batch extract Q matrices from STRUCTURE *_f files (mac-safe)
# Input:  outputs/v1/structure_runs/STRUCTURE_ZG_HEG-results/**/_f
# Output: outputs/v1/structure_runs/Q_extracted/*.Q  + summary CSV
############################################################

### 0) PROJECT ROOT + PATHS (ROBUST) ----
suppressPackageStartupMessages(library(here))

FALLBACK_ROOT <- file.path(path.expand("~"), "Desktop", "BeechCode")
candidate_root <- here::here()

if (dir.exists(file.path(candidate_root, "data", "raw"))) {
  PROJECT_ROOT <- candidate_root
} else if (dir.exists(file.path(FALLBACK_ROOT, "data", "raw"))) {
  PROJECT_ROOT <- FALLBACK_ROOT
  setwd(PROJECT_ROOT)
} else {
  stop(
    "Can't find project root.\n",
    "Expected to find a folder with data/raw.\n\n",
    "Tried:\n  - ", candidate_root, "\n  - ", FALLBACK_ROOT, "\n\n",
    "Fix: open your BeechCode.Rproj OR change FALLBACK_ROOT to your actual folder."
  )
}

OUTPUTS_DIR <- file.path(PROJECT_ROOT, "outputs")
RUN_TAG <- "v1"
RUN_OUT <- file.path(OUTPUTS_DIR, RUN_TAG)

STRUCTURE_RUNS_DIR <- file.path(RUN_OUT, "structure_runs")
RESULTS_DIR <- file.path(STRUCTURE_RUNS_DIR, "STRUCTURE_ZG_HEG-results")

if (!dir.exists(RESULTS_DIR)) {
  stop("RESULTS_DIR not found:\n  ", RESULTS_DIR,
       "\n\nCheck spelling/case in Finder and update RESULTS_DIR if needed.")
}

message("Project root: ", PROJECT_ROOT)
message("Using RESULTS_DIR: ", RESULTS_DIR)

# Where extracted Q files will go
Q_OUT_DIR <- file.path(STRUCTURE_RUNS_DIR, "Q_extracted")
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
message("Found ", length(f_files), " *_f files")

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
}

message("All extracted Q files are in: ", Q_OUT_DIR)
