############################################################
# scripts/06_pairwise_fst.R
# Pairwise Jost's D by Site (microsatellite data)
############################################################

options(repos = c(CRAN = "https://cloud.r-project.org"))
pkgs <- c("adegenet", "mmod", "dplyr", "tidyr", "tibble", "ggplot2")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)

suppressPackageStartupMessages({
  library(adegenet)
  library(mmod)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
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

if (nPop(gi) < 2) stop("Need at least 2 populations (Site levels) to compute pairwise Jost's D.")

# Keep only populations with non-empty labels in pop(gi), preserving existing pipeline object usage.
site_factor <- as.factor(pop(gi))
site_levels <- levels(site_factor)

# Enforce a minimum sample size per site for pairwise comparisons.
# Pairs involving sites below this threshold are kept as NA (with warning) instead of crashing.
min_n_per_site <- 2L
site_counts <- table(site_factor)
low_n_sites <- names(site_counts)[site_counts < min_n_per_site]
if (length(low_n_sites) > 0) {
  warning(
    "Sites with n < ", min_n_per_site,
    " will return NA in pairwise comparisons: ",
    paste(low_n_sites, collapse = ", ")
  )
}

# Initialize symmetric matrix with NA, force diagonal to 0.
jost_mat <- matrix(
  NA_real_,
  nrow = length(site_levels),
  ncol = length(site_levels),
  dimnames = list(site_levels, site_levels)
)
diag(jost_mat) <- 0

# Preferred method for microsatellite pairwise Jost's D:
# 1) Try mmod::pairwise_D(gi), which directly computes a full pairwise matrix from genind.
# 2) If that is incompatible in this runtime/data setup, fallback to per-pair mmod::D(two-pop genind).
used_method <- ""
pw_try <- tryCatch(mmod::pairwise_D(gi), error = function(e) e)

if (is.matrix(pw_try) && nrow(pw_try) > 0 && ncol(pw_try) > 0) {
  pw <- pw_try
  
  # Align returned matrix to site levels used by pop(gi)
  common_sites <- intersect(site_levels, intersect(rownames(pw), colnames(pw)))
  if (length(common_sites) >= 2) {
    jost_mat[common_sites, common_sites] <- pw[common_sites, common_sites, drop = FALSE]
    used_method <- "mmod::pairwise_D(genind)"
  }
}

if (used_method == "") {
  message("mmod::pairwise_D was not usable for this object; falling back to mmod::D on each site pair.")
  used_method <- "mmod::D(two-pop genind per pair)"
  
  for (i in seq_along(site_levels)) {
    for (j in i:length(site_levels)) {
      s1 <- site_levels[i]
      s2 <- site_levels[j]
      
      if (i == j) {
        jost_mat[s1, s2] <- 0
        next
      }
      
      n1 <- as.integer(site_counts[s1])
      n2 <- as.integer(site_counts[s2])
      if (is.na(n1) || is.na(n2) || n1 < min_n_per_site || n2 < min_n_per_site) {
        jost_mat[s1, s2] <- NA_real_
        jost_mat[s2, s1] <- NA_real_
        next
      }
      
      sub_gi <- gi[site_factor %in% c(s1, s2), , drop = FALSE]
      
      # Ensure sub-object population factor has exactly the two compared sites.
      pop(sub_gi) <- droplevels(pop(sub_gi))
      
      d_val <- tryCatch({
        as.numeric(mmod::D(sub_gi))
      }, error = function(e) {
        warning("Could not compute Jost's D for pair ", s1, " vs ", s2, ": ", e$message)
        NA_real_
      })
      
      jost_mat[s1, s2] <- d_val
      jost_mat[s2, s1] <- d_val
    }
  }
}

# Final matrix cleanup for robustness: symmetric + zero diagonal.
jost_mat <- (jost_mat + t(jost_mat)) / 2
diag(jost_mat) <- 0

if (!is.matrix(jost_mat) || nrow(jost_mat) == 0 || ncol(jost_mat) == 0) {
  stop("Jost's D matrix is empty.")
}

matrix_file <- file.path(OUTDIR, "pairwise_jostd_matrix.csv")
write.csv(jost_mat, matrix_file, row.names = TRUE)

jost_long <- as.data.frame(jost_mat) %>%
  tibble::rownames_to_column("Site1") %>%
  tidyr::pivot_longer(-Site1, names_to = "Site2", values_to = "JostD") %>%
  dplyr::filter(Site1 < Site2)

long_file <- file.path(OUTDIR, "pairwise_jostd_long.csv")
write.csv(jost_long, long_file, row.names = FALSE)

# Optional (preferred): heatmap figure.
heatmap_file <- file.path(OUTDIR, "pairwise_jostd_heatmap.png")
heat_df <- as.data.frame(jost_mat) %>%
  tibble::rownames_to_column("Site1") %>%
  tidyr::pivot_longer(-Site1, names_to = "Site2", values_to = "JostD")

p <- ggplot(heat_df, aes(x = Site1, y = Site2, fill = JostD)) +
  geom_tile(color = "white", linewidth = 0.1) +
  scale_fill_gradient(low = "#f7fbff", high = "#08306b", na.value = "grey85") +
  coord_equal() +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(title = "Pairwise Jost's D by Site", x = "Site", y = "Site", fill = "Jost's D")

ggsave(filename = heatmap_file, plot = p, width = 8.5, height = 7.5, dpi = 300)

# Console summary requested.
vals <- as.numeric(jost_mat[upper.tri(jost_mat, diag = FALSE)])
vals <- vals[!is.na(vals)]

site_n <- nrow(jost_mat)
if (length(vals) > 0) {
  summary_msg <- paste0(
    "DONE pairwise Jost's D (", used_method, ").\n",
    "Number of sites: ", site_n, "\n",
    "Jost's D (excluding NA) min/median/max: ",
    format(min(vals), digits = 4), " / ",
    format(stats::median(vals), digits = 4), " / ",
    format(max(vals), digits = 4), "\n",
    "Saved files:\n",
    "- ", matrix_file, "\n",
    "- ", long_file, "\n",
    "- ", heatmap_file
  )
} else {
  summary_msg <- paste0(
    "DONE pairwise Jost's D (", used_method, ").\n",
    "Number of sites: ", site_n, "\n",
    "No non-NA pairwise values were available after filtering.\n",
    "Saved files:\n",
    "- ", matrix_file, "\n",
    "- ", long_file, "\n",
    "- ", heatmap_file
  )
}

cat(summary_msg, "\n", sep = "")
