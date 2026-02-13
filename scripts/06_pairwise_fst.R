############################################################
# scripts/06_pairwise_fst.R
# Pairwise Weir & Cockerham FST by Site
############################################################

options(repos = c(CRAN = "https://cloud.r-project.org"))
pkgs <- c("adegenet", "hierfstat", "dplyr", "tidyr", "tibble")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)

suppressPackageStartupMessages({
  library(adegenet)
  library(hierfstat)
  library(dplyr)
  library(tidyr)
  library(tibble)
})

find_project_root <- function() {
  candidates <- c(getwd(), normalizePath(file.path(getwd(), ".."), mustWork = FALSE))
  cmd_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
  if (length(cmd_file) > 0 && nzchar(cmd_file[1])) {
    candidates <- c(candidates, dirname(normalizePath(cmd_file[1], mustWork = FALSE)))
  }
  for (start in unique(candidates)) {
    cur <- normalizePath(start, mustWork = FALSE)
    repeat {
      if (file.exists(file.path(cur, "scripts", "_load_objects.R"))) return(cur)
      parent <- dirname(cur)
      if (identical(parent, cur)) break
      cur <- parent
    }
  }
  stop("Cannot find project root containing scripts/_load_objects.R. Open BeechCode project first.")
}

setwd(find_project_root())
source(file.path("scripts", "_load_objects.R"))

OUTDIR <- file.path(RUN_OUT, "fst_only")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

if (nPop(gi) < 2) stop("Need at least 2 populations (Site levels) to compute pairwise FST.")

# IMPORTANT FIX: this function is from hierfstat, not adegenet
hf <- hierfstat::genind2hierfstat(gi)
fst_mat <- hierfstat::pairwise.WCfst(hf)

if (!is.matrix(fst_mat) || nrow(fst_mat) == 0 || ncol(fst_mat) == 0) {
  stop("pairwise.WCfst returned an empty result.")
}

write.csv(fst_mat, file.path(OUTDIR, "pairwise_fst_matrix.csv"), row.names = TRUE)

fst_long <- as.data.frame(fst_mat) %>%
  tibble::rownames_to_column("Site_1") %>%
  tidyr::pivot_longer(-Site_1, names_to = "Site_2", values_to = "FST") %>%
  dplyr::filter(Site_1 != Site_2)

write.csv(fst_long, file.path(OUTDIR, "pairwise_fst_long.csv"), row.names = FALSE)
cat("DONE pairwise FST. Outputs in: ", OUTDIR, "\n", sep = "")
