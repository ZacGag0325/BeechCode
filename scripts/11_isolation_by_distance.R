############################################################
# scripts/11_isolation_by_distance.R
# Mantel / IBE analyses from genetic + geographic distances
############################################################

options(repos = c(CRAN = "https://cloud.r-project.org"))
pkgs <- c("vegan", "ggplot2", "dplyr", "readr")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)

suppressPackageStartupMessages({
  library(vegan)
  library(ggplot2)
  library(dplyr)
  library(readr)
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
  stop("Cannot find project root containing scripts/_load_objects.R")
}

read_matrix <- function(path) {
  mat_df <- read.csv(path, row.names = 1, check.names = FALSE)
  mat <- as.matrix(mat_df)
  storage.mode(mat) <- "numeric"
  mat <- (mat + t(mat)) / 2
  diag(mat) <- 0
  mat
}

load_or_build_distance_matrices <- function(project_root, run_tag = "v1") {
  dist_dir <- file.path(project_root, "outputs", run_tag, "distances")
  dir.create(dist_dir, showWarnings = FALSE, recursive = TRUE)
  
  nei_csv <- file.path(dist_dir, "genetic_nei_dist.csv")
  geo_csv <- file.path(dist_dir, "geographic_dist.csv")
  
  if (!file.exists(nei_csv) || !file.exists(geo_csv)) {
    message("Distance matrices missing; computing on the fly via scripts/06_distance_matrices.R")
    source(file.path(project_root, "scripts", "06_distance_matrices.R"), local = FALSE)
  }
  
  nei_ok <- file.exists(nei_csv)
  geo_ok <- file.exists(geo_csv)
  
  if (!nei_ok || !geo_ok) {
    warning("Distance matrices still unavailable after attempted build.")
    return(list(nei = NULL, geo = NULL, nei_path = nei_csv, geo_path = geo_csv))
  }
  
  nei <- read_matrix(nei_csv)
  geo <- read_matrix(geo_csv)
  list(nei = nei, geo = geo, nei_path = nei_csv, geo_path = geo_csv)
}

align_matrices <- function(nei, geo) {
  common <- intersect(rownames(nei), rownames(geo))
  common <- sort(common)
  if (length(common) < 3) {
    warning("Need at least 3 shared sites to run Mantel. Shared sites: ", length(common))
    return(NULL)
  }
  nei2 <- nei[common, common, drop = FALSE]
  geo2 <- geo[common, common, drop = FALSE]
  list(nei = nei2, geo = geo2, sites = common)
}

PROJECT_ROOT <- find_project_root()
setwd(PROJECT_ROOT)
source(file.path("scripts", "_load_objects.R"))

RUN_TAG <- if (exists("RUN_TAG", inherits = TRUE)) get("RUN_TAG", inherits = TRUE) else "v1"
RUN_OUT <- if (exists("RUN_OUT", inherits = TRUE)) get("RUN_OUT", inherits = TRUE) else file.path(PROJECT_ROOT, "outputs", RUN_TAG)
MANTEL_DIR <- file.path(PROJECT_ROOT, "outputs", RUN_TAG, "mantel")
DIST_DIR <- file.path(PROJECT_ROOT, "outputs", RUN_TAG, "distances")

for (d in c(file.path(PROJECT_ROOT, "outputs"), file.path(PROJECT_ROOT, "outputs", RUN_TAG), DIST_DIR, MANTEL_DIR)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

mat_obj <- load_or_build_distance_matrices(PROJECT_ROOT, RUN_TAG)
if (is.null(mat_obj$nei) || is.null(mat_obj$geo)) {
  write.csv(
    data.frame(analysis = "mantel_nei_vs_geographic", status = "skipped", note = "Distance matrices unavailable", stringsAsFactors = FALSE),
    file.path(MANTEL_DIR, "mantel_results.csv"),
    row.names = FALSE
  )
  message("Mantel/IBE step skipped cleanly; distance matrices unavailable.")
} else {
  
  aligned <- align_matrices(mat_obj$nei, mat_obj$geo)
  if (is.null(aligned)) {
    write.csv(
      data.frame(analysis = "mantel_nei_vs_geographic", status = "skipped", note = "Too few shared sites", stringsAsFactors = FALSE),
      file.path(MANTEL_DIR, "mantel_results.csv"),
      row.names = FALSE
    )
    message("Mantel/IBE step skipped cleanly; too few shared sites.")
  } else {
    
    nei_mat <- aligned$nei
    geo_mat <- aligned$geo
    permutations_n <- 9999
    
    mantel_main <- vegan::mantel(as.dist(nei_mat), as.dist(geo_mat), method = "pearson", permutations = permutations_n)
    
    env_dist_path <- file.path(DIST_DIR, "environmental_dist.csv")
    partial_res <- NULL
    if (file.exists(env_dist_path)) {
      env_mat <- read_matrix(env_dist_path)
      common3 <- Reduce(intersect, list(rownames(nei_mat), rownames(geo_mat), rownames(env_mat)))
      if (length(common3) >= 3) {
        partial_res <- tryCatch(
          vegan::mantel.partial(
            as.dist(nei_mat[common3, common3, drop = FALSE]),
            as.dist(env_mat[common3, common3, drop = FALSE]),
            as.dist(geo_mat[common3, common3, drop = FALSE]),
            permutations = permutations_n,
            method = "pearson"
          ),
          error = function(e) {
            warning("Partial Mantel failed: ", conditionMessage(e))
            NULL
          }
        )
      } else {
        warning("Environmental distance matrix found, but fewer than 3 shared sites for partial Mantel.")
      }
    } else {
      message("Optional environmental distance matrix not found; skipping partial Mantel.")
    }
    
    pair_idx <- which(upper.tri(nei_mat, diag = FALSE), arr.ind = TRUE)
    scatter_df <- data.frame(
      Site1 = rownames(nei_mat)[pair_idx[, 1]],
      Site2 = colnames(nei_mat)[pair_idx[, 2]],
      geographic_km = geo_mat[pair_idx],
      genetic_nei = nei_mat[pair_idx],
      stringsAsFactors = FALSE
    )
    
    plot_obj <- ggplot(scatter_df, aes(x = geographic_km, y = genetic_nei)) +
      geom_point(size = 2, alpha = 0.85, color = "#1f78b4") +
      geom_smooth(method = "lm", se = TRUE, color = "#d7301f", fill = "#fcbba1", linewidth = 0.9) +
      theme_minimal(base_size = 12) +
      labs(
        title = "Isolation by distance (Mantel)",
        subtitle = paste0("r = ", signif(as.numeric(mantel_main$statistic), 4), ", p = ", signif(as.numeric(mantel_main$signif), 4)),
        x = "Geographic distance (km)",
        y = "Nei genetic distance"
      )
    
    ggsave(file.path(MANTEL_DIR, "genetic_vs_geographic_scatter.png"), plot_obj, width = 8, height = 6, dpi = 300)
    
    results <- data.frame(
      analysis = "mantel_nei_vs_geographic",
      statistic_r = as.numeric(mantel_main$statistic),
      p_value = as.numeric(mantel_main$signif),
      permutations = as.integer(mantel_main$permutations),
      n_sites = nrow(nei_mat),
      n_pairs = nrow(scatter_df),
      status = "ok",
      note = NA_character_,
      stringsAsFactors = FALSE
    )
    
    if (!is.null(partial_res)) {
      results <- bind_rows(
        results,
        data.frame(
          analysis = "partial_mantel_genetic_vs_environment_given_geographic",
          statistic_r = as.numeric(partial_res$statistic),
          p_value = as.numeric(partial_res$signif),
          permutations = as.integer(partial_res$permutations),
          n_sites = length(attr(as.dist(nei_mat), "Size")),
          n_pairs = nrow(scatter_df),
          status = "ok",
          note = "environmental_dist.csv used",
          stringsAsFactors = FALSE
        )
      )
    }
    
    write.csv(results, file.path(MANTEL_DIR, "mantel_results.csv"), row.names = FALSE)
    write.csv(scatter_df, file.path(MANTEL_DIR, "mantel_pairwise_points.csv"), row.names = FALSE)
    
    message("Created Mantel results table: ", file.path(MANTEL_DIR, "mantel_results.csv"))
    message("Created Mantel scatter data: ", file.path(MANTEL_DIR, "mantel_pairwise_points.csv"))
    message("Created Mantel scatter plot: ", file.path(MANTEL_DIR, "genetic_vs_geographic_scatter.png"))
    message("Mantel/IBE step complete. Output folder: ", MANTEL_DIR)
    
  }
}
