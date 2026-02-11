############################################################
# scripts/04_amova.R
# AMOVA among Sites (+ permutation test)
############################################################

options(repos = c(CRAN = "https://cloud.r-project.org"))

pkgs <- c("adegenet", "poppr", "ade4")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)

suppressPackageStartupMessages(library(adegenet))
suppressPackageStartupMessages(library(poppr))
suppressPackageStartupMessages(library(ade4))

source("scripts/_load_objects.R")

OUTDIR <- file.path(RUN_OUT, "amova_only")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# Ensure pop is Site (consistent)
pop(gi) <- as.factor(df_ids$Site[match(indNames(gi), df_ids$ind)])

# AMOVA (ade4-style object returned by poppr.amova)
am <- poppr::poppr.amova(gi, ~Site, within = TRUE)
print(am)

# Permutation test (THIS is the correct randtest)
set.seed(123)
am_test <- ade4::randtest(am, nrepet = 999)
print(am_test)

# Save outputs
capture.output(am,      file = file.path(OUTDIR, "amova_site.txt"))
capture.output(am_test, file = file.path(OUTDIR, "amova_site_randtest.txt"))

# Optional: save a quick plot of permutation distribution
pdf(file.path(OUTDIR, "amova_site_randtest.pdf"), width = 7, height = 5)
plot(am_test)
dev.off()

cat("DONE AMOVA. Outputs in: ", OUTDIR, "\n", sep = "")
