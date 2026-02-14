# scripts/08_ibs_ibd.R
############################################################
# scripts/08_ibs_ibd.R
# Pairwise individual genetic/geographic distances
############################################################

options(repos = c(CRAN = "https://cloud.r-project.org"))
pkgs <- c("adegenet", "dplyr")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
suppressPackageStartupMessages({library(adegenet); library(dplyr)})

source(file.path("scripts", "_load_objects.R"))

OUTDIR <- file.path(RUN_OUT, "ibs_ibd_only")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

haversine_km <- function(lat1, lon1, lat2, lon2) {
  rad <- pi / 180
  dlat <- (lat2 - lat1) * rad
  dlon <- (lon2 - lon1) * rad
  a <- sin(dlat / 2)^2 + cos(lat1 * rad) * cos(lat2 * rad) * sin(dlon / 2)^2
  2 * 6371 * asin(pmin(1, sqrt(a)))
}

X <- adegenet::tab(gi, freq = TRUE, NA.method = "mean")
gen_dist <- as.matrix(dist(X, method = "euclidean"))
inds <- rownames(gen_dist)

ids <- df_ids %>%
  dplyr::select(ind, Site, Latitude, Longitude) %>%
  dplyr::distinct(ind, .keep_all = TRUE)

pair_rows <- vector("list", length = max(1, (length(inds) * (length(inds) - 1)) / 2))
k <- 1L
for (i in seq_len(length(inds) - 1)) {
  for (j in (i + 1):length(inds)) {
    id1 <- inds[i]
    id2 <- inds[j]
    r1 <- ids[ids$ind == id1, , drop = FALSE]
    r2 <- ids[ids$ind == id2, , drop = FALSE]
    if (nrow(r1) == 0 || nrow(r2) == 0) next
    geo <- if (all(!is.na(c(r1$Latitude, r1$Longitude, r2$Latitude, r2$Longitude)))) {
      haversine_km(r1$Latitude, r1$Longitude, r2$Latitude, r2$Longitude)
    } else {
      NA_real_
    }
    
    pair_rows[[k]] <- data.frame(
      ind_1 = id1,
      ind_2 = id2,
      site_1 = as.character(r1$Site),
      site_2 = as.character(r2$Site),
      genetic_distance = gen_dist[i, j],
      geographic_km = geo,
      stringsAsFactors = FALSE
    )
    k <- k + 1L
  }
}

pairs_df <- dplyr::bind_rows(pair_rows)
if (nrow(pairs_df) == 0) {
  stop("No pairwise rows were generated in scripts/08_ibs_ibs.R")
}

write.csv(pairs_df, file.path(OUTDIR, "pairwise_individual_distances.csv"), row.names = FALSE)
cat("DONE IBS/IBD distance extraction. Outputs in: ", OUTDIR, "\n", sep = "")
