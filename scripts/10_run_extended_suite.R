# scripts/10_run_extended_suite.R
############################################################
# scripts/10_run_extended_suite.R
# Runner for optional extended analyses
############################################################

suppressPackageStartupMessages(library(here))

FALLBACK_ROOT <- file.path(path.expand("~"), "Desktop", "BeechCode")
candidate_root <- here::here()

if (dir.exists(file.path(candidate_root, "scripts"))) {
  PROJECT_ROOT <- candidate_root
} else if (dir.exists(file.path(FALLBACK_ROOT, "scripts"))) {
  PROJECT_ROOT <- FALLBACK_ROOT
  setwd(PROJECT_ROOT)
} else {
  stop("Can't find project root. Open BeechCode.Rproj or set FALLBACK_ROOT.")
}

cat("Project root:", PROJECT_ROOT, "\n")

run_one <- function(path) {
  f <- file.path(PROJECT_ROOT, path)
  out <- list(path = path, status = "ok", msg = "")
  
  if (!file.exists(f)) {
    out$status <- "skipped_missing"
    out$msg <- "file not found"
    cat("SKIP (missing):", path, "\n")
    return(out)
  }
  
  cat("\n========================================================================\n")
  cat("RUNNING:", basename(path), "\n")
  cat("========================================================================\n")
  
  res <- tryCatch({
    source(f, local = FALSE, encoding = "UTF-8")
    NULL
  }, error = function(e) e)
  
  if (inherits(res, "error")) {
    out$status <- "failed"
    out$msg <- conditionMessage(res)
    cat("FAILED:", path, "\n", "  -> ", out$msg, "\n", sep = "")
  }
  
  out
}

extended_scripts <- c(
  "scripts/08_ibs_ibd.R"
)

results <- lapply(extended_scripts, run_one)
status <- vapply(results, `[[`, character(1), "status")

cat("\nExtended suite complete.\n")
cat(" - ran ok:      ", sum(status == "ok"), "\n", sep = "")
cat(" - skipped:     ", sum(status == "skipped_missing"), "\n", sep = "")
cat(" - failed:      ", sum(status == "failed"), "\n", sep = "")

if (any(status == "failed")) {
  cat("\nFailures:\n")
  for (i in which(status == "failed")) {
    cat(" - ", results[[i]]$path, ": ", results[[i]]$msg, "\n", sep = "")
  }
}
