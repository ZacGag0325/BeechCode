############################################################
# 00_extract_Q_from_f.R  (BeechCode)
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
dir.create(STRUCTURE_RUNS_DIR, showWarnings = FALSE, recursive = TRUE)

message("Project root: ", PROJECT_ROOT)
message("RUN_OUT: ", RUN_OUT)
message("STRUCTURE_RUNS_DIR: ", STRUCTURE_RUNS_DIR)

### 1) Pick correct results folder name ----
cand1 <- file.path(STRUCTURE_RUNS_DIR, "STRUCTURE_ZG_HEG-RESULTS")
cand2 <- file.path(STRUCTURE_RUNS_DIR, "STRUCTURE_ZG_HEG-results")

if (dir.exists(cand1)) {
  RESULTS_DIR <- cand1
} else if (dir.exists(cand2)) {
  RESULTS_DIR <- cand2
} else {
  stop(
    "Could not find your STRUCTURE results folder.\n",
    "Looked for:\n  - ", cand1, "\n  - ", cand2
  )
}

message("Using RESULTS_DIR: ", RESULTS_DIR)

### 2) Output folder ----
Q_OUT_DIR <- file.path(STRUCTURE_RUNS_DIR, "Q_extracted")
dir.create(Q_OUT_DIR, showWarnings = FALSE, recursive = TRUE)
message("Q files will be written to: ", Q_OUT_DIR)

### 3) Helpers ----

# Extract ALL numeric tokens from a string (regex), return numeric vector (or NULL)
extract_numbers_regex <- function(s) {
  s <- trimws(s)
  if (s == "") return(NULL)
  
  # match numbers like 0.123, .123, 1, 1.0, 1e-3, -2.5E+02
  m <- gregexpr("[-+]?(?:\\d+\\.?\\d*|\\.\\d+)(?:[eE][-+]?\\d+)?", s, perl = TRUE)
  hits <- regmatches(s, m)[[1]]
  if (length(hits) == 0) return(NULL)
  
  suppressWarnings(v <- as.numeric(hits))
  v <- v[is.finite(v)]
  if (length(v) == 0) return(NULL)
  
  v
}

# Guess K from filename like "...-K11-4_f"
guess_K_from_name <- function(f) {
  m <- regmatches(f, regexpr("K\\d+", f, ignore.case = TRUE))
  if (length(m) == 0 || is.na(m)) return(NA_integer_)
  as.integer(gsub("[^0-9]", "", m))
}

# Extract Q matrix from a STRUCTURE _f file
extract_Q_from_f <- function(f) {
  lines <- readLines(f, warn = FALSE)
  
  # Find ancestry section start
  start <- grep("Inferred ancestry of individuals", lines, ignore.case = TRUE)
  if (length(start) == 0) {
    start <- grep("Inferred ancestry|Membership coefficients|Estimated membership", lines, ignore.case = TRUE)
  }
  if (length(start) == 0) return(list(Q = NULL, reason = "no_ancestry_header"))
  
  i0 <- start[1]
  q_rows <- list()
  started <- FALSE
  
  for (j in (i0 + 1):length(lines)) {
    ln <- lines[j]
    
    # Once collecting, stop at blank line or new section
    if (started) {
      if (trimws(ln) == "") break
      if (grepl("^\\s*Run\\s+parameters|^\\s*Estimated|^\\s*Allele|^\\s*Log\\s|^\\s*Average", ln, ignore.case = TRUE)) break
    }
    
    # Must contain ":" somewhere (STRUCTURE lines usually have it)
    if (!grepl(":", ln, fixed = TRUE)) next
    
    # Use text after the LAST ":" on the line
    pos <- gregexpr(":", ln, fixed = TRUE)[[1]]
    if (pos[1] == -1) next
    last_colon <- tail(pos, 1)
    rhs <- substr(ln, last_colon + 1, nchar(ln))
    
    nums <- extract_numbers_regex(rhs)
    if (is.null(nums)) next
    
    # Keep only plausible Q values in [0,1] (tolerate tiny float error)
    q <- nums[nums >= -1e-6 & nums <= 1 + 1e-6]
    if (length(q) < 1) next
    
    # Sum sanity (tolerate rounding / formatting)
    s <- sum(q)
    if (abs(s - 1) > 0.30) next
    
    q_rows[[length(q_rows) + 1]] <- q
    started <- TRUE
  }
  
  if (length(q_rows) == 0) return(list(Q = NULL, reason = "no_q_rows_found"))
  
  # Determine modal K length in this file (most common row length)
  lens <- vapply(q_rows, length, integer(1))
  K_mode <- as.integer(names(sort(table(lens), decreasing = TRUE)[1]))
  
  # Keep only rows with modal length
  keep_idx <- which(lens == K_mode)
  q_rows <- q_rows[keep_idx]
  
  if (length(q_rows) < 10) return(list(Q = NULL, reason = "too_few_rows_after_modalK"))
  
  # Build matrix
  M <- do.call(rbind, lapply(q_rows, function(v) as.numeric(v[1:K_mode])))
  
  # Final sanity
  if (ncol(M) < 1 || nrow(M) < 10) return(list(Q = NULL, reason = "bad_dimensions"))
  
  list(Q = M, reason = "ok")
}

### 4) Find all *_f files ----
f_files <- list.files(RESULTS_DIR, recursive = TRUE, full.names = TRUE, pattern = "_f$")
message("Found ", length(f_files), " *_f files in RESULTS_DIR")
if (length(f_files) == 0) stop("No *_f files found under: ", RESULTS_DIR)

### 5) Extract from all files ----
results <- data.frame(
  f_file = character(0),
  q_file = character(0),
  K_detected = integer(0),
  n_inds = integer(0),
  K_from_name = integer(0),
  ok = logical(0),
  fail_reason = character(0),
  stringsAsFactors = FALSE
)

ok_count <- 0
fail_count <- 0

for (f in f_files) {
  out <- tryCatch(extract_Q_from_f(f), error = function(e) list(Q = NULL, reason = "error"))
  Q <- out$Q
  reason <- out$reason
  
  if (is.null(Q)) {
    fail_count <- fail_count + 1
    results <- rbind(results, data.frame(
      f_file = f,
      q_file = NA_character_,
      K_detected = NA_integer_,
      n_inds = NA_integer_,
      K_from_name = guess_K_from_name(f),
      ok = FALSE,
      fail_reason = reason
    ))
    next
  }
  
  base <- basename(f)
  q_out <- file.path(Q_OUT_DIR, paste0(base, ".Q"))
  write.table(Q, file = q_out, quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  ok_count <- ok_count + 1
  results <- rbind(results, data.frame(
    f_file = f,
    q_file = q_out,
    K_detected = ncol(Q),
    n_inds = nrow(Q),
    K_from_name = guess_K_from_name(f),
    ok = TRUE,
    fail_reason = "ok"
  ))
}

message("Extraction finished. OK: ", ok_count, " | Failed: ", fail_count)

summary_csv <- file.path(Q_OUT_DIR, "extracted_Q_summary.csv")
write.csv(results, summary_csv, row.names = FALSE)
message("Saved summary CSV: ", summary_csv)

if (ok_count > 0) {
  tab <- sort(table(results$K_detected[results$ok]))
  message("Counts per K (detected from extracted Q):")
  print(tab)
}

message("All extracted Q files are in: ", Q_OUT_DIR)
