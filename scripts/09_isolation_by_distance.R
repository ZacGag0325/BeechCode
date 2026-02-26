# scripts/09_isolation_by_distance.R
############################################################
# scripts/09_isolation_by_distance.R
############################################################

options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
suppressPackageStartupMessages(library(dplyr))

source(file.path("scripts", "_load_objects.R"))

OUTDIR <- file.path(RUN_OUT, "ibd_only")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

infile <- file.path(RUN_OUT, "ibs_ibd_only", "pairwise_individual_distances.csv")
if (!file.exists(infile)) source(file.path("scripts", "08_ibs_ibs.R"))

if (!file.exists(infile)) stop("Expected pairwise_individual_distances.csv was not created.")

d <- read.csv(infile, stringsAsFactors = FALSE) %>%
  dplyr::filter(!is.na(genetic_distance), !is.na(geographic_km))

if (nrow(d) < 3) stop("Not enough complete pairs for IBD regression.")

fit <- stats::lm(genetic_distance ~ geographic_km, data = d)
ct <- suppressWarnings(stats::cor.test(d$genetic_distance, d$geographic_km, method = "spearman"))

summary_tbl <- data.frame(
  n_pairs = nrow(d),
  slope = unname(coef(fit)["geographic_km"]),
  intercept = unname(coef(fit)["(Intercept)"]),
  r_squared = summary(fit)$r.squared,
  spearman_rho = unname(ct$estimate),
  spearman_p = ct$p.value
)

write.csv(summary_tbl, file.path(OUTDIR, "ibd_summary.csv"), row.names = FALSE)
write.csv(d, file.path(OUTDIR, "ibd_pairs_used.csv"), row.names = FALSE)
cat("DONE isolation-by-distance summary. Outputs in: ", OUTDIR, "\n", sep = "")
