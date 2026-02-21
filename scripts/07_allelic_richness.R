# filename: scripts/07_allelic_richness.R
############################################################
# scripts/07_allelic_richness.R
# Per-site allelic richness summary (robust + non-crashing)
############################################################

options(repos = c(CRAN = "https://cloud.r-project.org"))
pkgs <- c("adegenet", "hierfstat", "dplyr")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
suppressPackageStartupMessages({library(adegenet); library(hierfstat); library(dplyr)})

find_project_root <- function() {
  candidates <- c(getwd(), normalizePath(file.path(getwd(), ".."), mustWork = FALSE))
  for (start in unique(candidates)) {
    cur <- normalizePath(start, mustWork = FALSE)
    repeat {
      if (file.exists(file.path(cur, "scripts", "_load_objects.R"))) return(cur)
      parent <- dirname(cur); if (identical(parent, cur)) break; cur <- parent
    }
  }
  stop("Cannot find project root containing scripts/_load_objects.R")
}

compute_ar <- function(gen) {
  if (!inherits(gen, "genind")) stop("Input must be a genind object.")
  if (nPop(gen) < 2) {
    warning("Need at least two sites for allelic richness summaries.")
    return(NULL)
  }
  
  pop_sizes <- table(pop(gen))
  min_n <- min(pop_sizes)
  if (is.na(min_n) || min_n < 2) {
    warning("At least one site has fewer than 2 individuals; allelic richness cannot be robustly computed.")
    return(NULL)
  }
  
  hf <- hierfstat::genind2hierfstat(gen)
  ar <- hierfstat::allelic.richness(hf, min.n = min_n)
  
  ar_mat <- ar$Ar
  if (is.null(ar_mat) || ncol(ar_mat) < 1) {
    warning("allelic.richness returned no per-site values.")
    return(NULL)
  }
  
  data.frame(
    Site = colnames(ar_mat),
    n = as.integer(pop_sizes[colnames(ar_mat)]),
    allelic_richness = as.numeric(colMeans(ar_mat, na.rm = TRUE)),
    allelic_richness_se = as.numeric(apply(ar_mat, 2, sd, na.rm = TRUE) / sqrt(nrow(ar_mat))),
    stringsAsFactors = FALSE
  )
}

PROJECT_ROOT <- find_project_root()
setwd(PROJECT_ROOT)
source(file.path("scripts", "_load_objects.R"))

RUN_TAG <- if (exists("RUN_TAG", inherits = TRUE)) get("RUN_TAG", inherits = TRUE) else "v1"
RUN_OUT <- if (exists("RUN_OUT", inherits = TRUE)) get("RUN_OUT", inherits = TRUE) else file.path(PROJECT_ROOT, "outputs", RUN_TAG)
DIVERSITY_DIR <- file.path(PROJECT_ROOT, "outputs", RUN_TAG, "diversity")
dir.create(file.path(PROJECT_ROOT, "outputs"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(PROJECT_ROOT, "outputs", RUN_TAG), showWarnings = FALSE, recursive = TRUE)
dir.create(DIVERSITY_DIR, showWarnings = FALSE, recursive = TRUE)

ar_site <- tryCatch(compute_ar(gi_mll), error = function(e) {
  warning("Allelic richness computation failed: ", conditionMessage(e))
  NULL
})

csv_path <- file.path(DIVERSITY_DIR, "allelic_richness_by_site.csv")
rds_path <- file.path(DIVERSITY_DIR, "allelic_richness_by_site.rds")
legacy_csv <- file.path(RUN_OUT, "allelic_richness_site_summary.csv")

if (is.null(ar_site) || nrow(ar_site) == 0) {
  warning("No allelic richness output produced (insufficient data).")
} else {
  write.csv(ar_site, csv_path, row.names = FALSE)
  saveRDS(ar_site, rds_path)
  write.csv(ar_site, legacy_csv, row.names = FALSE)
  cat("Created allelic richness summary CSV: ", csv_path, "\n", sep = "")
  cat("Created allelic richness summary RDS: ", rds_path, "\n", sep = "")
  cat("Created legacy compatibility CSV: ", legacy_csv, "\n", sep = "")
}

cat("Allelic-richness step complete. Output folder: ", DIVERSITY_DIR, "\n", sep = "")
