#!/usr/bin/env Rscript

# list_package_usage.R
# Scan all R scripts in scripts/ and summarize R package usage.

options(stringsAsFactors = FALSE)

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
scripts_dir <- file.path(project_root, "scripts")
output_dir <- file.path(project_root, "outputs", "tables")

if (!dir.exists(scripts_dir)) {
  stop("Could not find scripts directory: ", scripts_dir)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

r_files <- list.files(scripts_dir, pattern = "\\.[Rr]$", full.names = TRUE)
r_files <- r_files[basename(r_files) != "list_package_usage.R"]
r_files <- sort(r_files)

if (length(r_files) == 0) {
  stop("No .R files found in ", scripts_dir)
}

extract_matches <- function(pattern, text, group = 1L) {
  m <- gregexpr(pattern, text, perl = TRUE)
  reg <- regmatches(text, m)
  out <- character(0)
  for (hits in reg) {
    if (length(hits) == 0) next
    groups <- regexec(pattern, hits, perl = TRUE)
    cap <- regmatches(hits, groups)
    for (x in cap) {
      if (length(x) >= group && nzchar(x[group])) out <- c(out, x[group])
    }
  }
  unique(out)
}

`%||%` <- function(x, y) if (length(x) == 0 || all(is.na(x)) || !nzchar(x[1])) y else x

infer_purpose <- function(lines, file_name) {
  top <- lines[seq_len(min(length(lines), 80))]
  comments <- trimws(gsub("^\\s*#+'?\\s?", "", top))
  comments <- comments[nzchar(comments)]
  comments <- comments[!grepl("^(----+|====+)$", comments)]
  
  if (length(comments) > 0) {
    candidate <- comments[1]
    return(candidate)
  }
  
  stem <- sub("\\.[Rr]$", "", file_name)
  stem <- gsub("[_\\-]+", " ", stem)
  stem <- gsub("\\s+", " ", stem)
  stem <- trimws(stem)
  paste("Analysis script:", stem)
}

# Heuristic mapping from common function names to package roles.
role_dictionary <- list(
  dplyr = "data wrangling and summaries",
  tidyr = "data reshaping",
  ggplot2 = "data visualization",
  readr = "read/write delimited data",
  readxl = "read Excel files",
  writexl = "write Excel files",
  openxlsx = "write Excel files",
  stringr = "string cleaning and parsing",
  purrr = "functional iteration/mapping",
  tibble = "tibble/data-frame helpers",
  rlang = "tidy evaluation helpers",
  adegenet = "genetic marker data structures and analyses",
  poppr = "clonal/genotypic population genetics",
  hierfstat = "F-statistics and allelic richness",
  pegas = "HWE and population genetics tests",
  mmod = "genetic differentiation metrics (e.g., Jost's D)",
  vegan = "Mantel tests and ecological distance statistics",
  geosphere = "geographic distance calculations",
  ade4 = "multivariate stats and randtest",
  ape = "phylogenetic distance/tree tools",
  ggtree = "phylogenetic tree visualization/annotation",
  cowplot = "plot composition/legend extraction",
  tidyverse = "meta-package for tidy data workflow",
  here = "project-root path handling",
  stats = "base statistical functions",
  tools = "file/path utilities",
  utils = "base utility functions"
)

build_reason <- function(pkg, ns_funcs, loaded_only = FALSE) {
  role <- role_dictionary[[pkg]] %||% "general package functionality used in this script"
  
  if (length(ns_funcs) > 0) {
    funcs <- sort(unique(ns_funcs))
    funcs_show <- paste(utils::head(funcs, 6), collapse = ", ")
    if (length(funcs) > 6) funcs_show <- paste0(funcs_show, ", ...")
    return(sprintf("%s (namespace calls: %s)", role, funcs_show))
  }
  
  if (loaded_only) {
    return(sprintf("%s (loaded via library/require/requireNamespace; unqualified calls may be used)", role))
  }
  
  sprintf("%s", role)
}

script_rows <- list()
package_rows <- list()

for (f in r_files) {
  lines <- readLines(f, warn = FALSE, encoding = "UTF-8")
  txt <- paste(lines, collapse = "\n")
  
  lib_pkgs <- c(
    extract_matches("(?:library|require)\\s*\\(\\s*([A-Za-z][A-Za-z0-9._]*)\\s*\\)", txt, group = 1L),
    extract_matches("(?:library|require)\\s*\\(\\s*['\"]([A-Za-z][A-Za-z0-9._]*)['\"]\\s*\\)", txt, group = 1L)
  )
  reqns_pkgs <- c(
    extract_matches("requireNamespace\\s*\\(\\s*\"([A-Za-z][A-Za-z0-9._]*)\"", txt, group = 1L),
    extract_matches("requireNamespace\\s*\\(\\s*'([A-Za-z][A-Za-z0-9._]*)'", txt, group = 1L)
  )
  
  ns_calls <- unique(unlist(regmatches(
    txt,
    gregexpr("[A-Za-z][A-Za-z0-9._]*:::{0,1}[A-Za-z][A-Za-z0-9._]*", txt, perl = TRUE)
  )))
  
  ns_pkg <- character(0)
  ns_fun <- character(0)
  if (length(ns_calls) > 0) {
    ns_pkg <- sub("::.*$", "", ns_calls)
    ns_fun <- sub("^.*:::{0,1}", "", ns_calls)
  }
  
  all_pkgs <- sort(unique(c(lib_pkgs, reqns_pkgs, ns_pkg)))
  
  pkg_reason_vec <- character(0)
  
  for (pkg in all_pkgs) {
    idx <- which(ns_pkg == pkg)
    funcs <- if (length(idx) > 0) ns_fun[idx] else character(0)
    reason <- build_reason(pkg, ns_funcs = funcs, loaded_only = (length(funcs) == 0))
    pkg_reason_vec <- c(pkg_reason_vec, sprintf("%s: %s", pkg, reason))
    
    package_rows[[length(package_rows) + 1L]] <- data.frame(
      package_name = pkg,
      script_name = basename(f),
      reason_in_script = reason,
      stringsAsFactors = FALSE
    )
  }
  
  script_rows[[length(script_rows) + 1L]] <- data.frame(
    script_name = basename(f),
    analysis_purpose = infer_purpose(lines, basename(f)),
    packages_used = if (length(all_pkgs) > 0) paste(all_pkgs, collapse = "; ") else "",
    package_reasons = if (length(pkg_reason_vec) > 0) paste(pkg_reason_vec, collapse = " | ") else "",
    stringsAsFactors = FALSE
  )
}

by_script <- do.call(rbind, script_rows)
by_script <- by_script[order(by_script$script_name), , drop = FALSE]

if (length(package_rows) > 0) {
  pkg_long <- do.call(rbind, package_rows)
  
  package_split <- split(pkg_long, pkg_long$package_name)
  summary_rows <- lapply(package_split, function(df) {
    scripts_used <- sort(unique(df$script_name))
    reasons <- sort(unique(df$reason_in_script))
    
    data.frame(
      package_name = df$package_name[1],
      scripts_used = paste(scripts_used, collapse = "; "),
      overall_role = paste(utils::head(reasons, 3), collapse = " | "),
      stringsAsFactors = FALSE
    )
  })
  package_summary <- do.call(rbind, summary_rows)
  package_summary <- package_summary[order(package_summary$package_name), , drop = FALSE]
} else {
  package_summary <- data.frame(
    package_name = character(0),
    scripts_used = character(0),
    overall_role = character(0),
    stringsAsFactors = FALSE
  )
}

path_by_script <- file.path(output_dir, "package_usage_by_script.csv")
path_summary <- file.path(output_dir, "package_usage_summary.csv")
path_xlsx <- file.path(output_dir, "package_usage_by_script.xlsx")

utils::write.csv(by_script, path_by_script, row.names = FALSE, fileEncoding = "UTF-8")
utils::write.csv(package_summary, path_summary, row.names = FALSE, fileEncoding = "UTF-8")

if (requireNamespace("openxlsx", quietly = TRUE)) {
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "by_script")
  openxlsx::writeData(wb, "by_script", by_script)
  openxlsx::addWorksheet(wb, "package_summary")
  openxlsx::writeData(wb, "package_summary", package_summary)
  openxlsx::saveWorkbook(wb, path_xlsx, overwrite = TRUE)
} else if (requireNamespace("writexl", quietly = TRUE)) {
  writexl::write_xlsx(list(by_script = by_script, package_summary = package_summary), path = path_xlsx)
} else {
  message("Skipping XLSX export (neither openxlsx nor writexl is available).")
}

cat("\n=== Package usage by script ===\n")
print(by_script, row.names = FALSE, right = FALSE)

cat("\n=== Package usage summary ===\n")
print(package_summary, row.names = FALSE, right = FALSE)

cat("\nSaved files:\n")
cat(" - ", path_by_script, "\n", sep = "")
cat(" - ", path_summary, "\n", sep = "")
if (file.exists(path_xlsx)) {
  cat(" - ", path_xlsx, "\n", sep = "")
}