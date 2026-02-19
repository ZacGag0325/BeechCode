# filename: scripts/07_allelic_richness.R
############################################################
# scripts/07_allelic_richness.R
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

compute_ar <- function(gen, dataset_label) {
  if (!inherits(gen, "genind")) stop("Input must be a genind object.")
  if (nPop(gen) < 2) stop("Need at least two sites for allelic richness summaries.")
  
  hf <- hierfstat::genind2hierfstat(gen)
  min_n <- min(table(pop(gen)))
  ar <- hierfstat::allelic.richness(hf, min.n = min_n)
  
  ar_locus_site <- as.data.frame(ar$Ar)
  ar_locus_site$Locus <- rownames(ar_locus_site)
  ar_locus_site$dataset <- dataset_label
  
  ar_site_summary <- data.frame(
    Site = colnames(ar$Ar),
    mean_allelic_richness = colMeans(ar$Ar, na.rm = TRUE),
    dataset = dataset_label,
    stringsAsFactors = FALSE
  )
  
  list(locus = ar_locus_site, site = ar_site_summary)
}

setwd(find_project_root())
source(file.path("scripts", "_load_objects.R"))

OUTDIR <- file.path(RUN_OUT, "allelic_richness_only")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

ar_raw <- compute_ar(gi, "raw")
ar_cc <- compute_ar(gi_mll, "clone_corrected")

write.csv(bind_rows(ar_raw$locus, ar_cc$locus),
          file.path(OUTDIR, "allelic_richness_locus_by_site_raw_vs_clone_corrected.csv"),
          row.names = FALSE)
write.csv(bind_rows(ar_raw$site, ar_cc$site),
          file.path(OUTDIR, "allelic_richness_site_summary_raw_vs_clone_corrected.csv"),
          row.names = FALSE)

# Backward-compatible filename now clone-corrected by default
write.csv(ar_cc$locus, file.path(OUTDIR, "allelic_richness_locus_by_site.csv"), row.names = FALSE)
write.csv(ar_cc$site, file.path(OUTDIR, "allelic_richness_site_summary.csv"), row.names = FALSE)

cat("DONE allelic richness (raw + clone-corrected). Outputs in: ", OUTDIR, "\n", sep = "")
