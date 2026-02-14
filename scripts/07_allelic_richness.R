############################################################
# scripts/07_allelic_richness.R
############################################################
options(repos = c(CRAN = "https://cloud.r-project.org"))
pkgs <- c("adegenet","hierfstat","dplyr")
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
setwd(find_project_root())
source(file.path("scripts","_load_objects.R"))
OUTDIR <- file.path(RUN_OUT, "allelic_richness_only")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
hf <- hierfstat::genind2hierfstat(gi)
min_n <- min(table(pop(gi)))
ar <- hierfstat::allelic.richness(hf, min.n = min_n)
ar_locus_site <- as.data.frame(ar$Ar); ar_locus_site$Locus <- rownames(ar_locus_site)
ar_site_summary <- data.frame(Site = colnames(ar$Ar), mean_allelic_richness = colMeans(ar$Ar, na.rm = TRUE))
write.csv(ar_locus_site, file.path(OUTDIR, "allelic_richness_locus_by_site.csv"), row.names = FALSE)
write.csv(ar_site_summary, file.path(OUTDIR, "allelic_richness_site_summary.csv"), row.names = FALSE)
cat("DONE allelic richness. Outputs in: ", OUTDIR, "\n", sep = "")
