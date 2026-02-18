############################################################
# scripts/06_pairwise_fst.R
# Pairwise Jost's D by Site (with permutation p-values)
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

set.seed(20250218)

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

site_factor <- as.factor(pop(gi))
site_levels <- levels(site_factor)

min_n_per_site <- 2L
n_perm <- 999L

site_counts <- table(site_factor)
low_n_sites <- names(site_counts)[site_counts < min_n_per_site]
if (length(low_n_sites) > 0) {
  warning(
    "Sites with n < ", min_n_per_site,
    " will return NA in pairwise comparisons: ",
    paste(low_n_sites, collapse = ", ")
  )
}

jost_mat <- matrix(
  NA_real_,
  nrow = length(site_levels),
  ncol = length(site_levels),
  dimnames = list(site_levels, site_levels)
)
diag(jost_mat) <- 0

pair_rows <- list()
row_id <- 1L

compute_pair_D <- function(genind_obj, s1, s2) {
  sub_gi <- genind_obj[pop(genind_obj) %in% c(s1, s2), , drop = FALSE]
  pop(sub_gi) <- droplevels(pop(sub_gi))
  as.numeric(mmod::D(sub_gi))
}

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
      pair_rows[[row_id]] <- data.frame(
        Site1 = s1,
        Site2 = s2,
        n_site1 = n1,
        n_site2 = n2,
        JostD = NA_real_,
        p_value = NA_real_,
        reason = paste0("too_few_individuals (min_n_per_site=", min_n_per_site, ")"),
        stringsAsFactors = FALSE
      )
      row_id <- row_id + 1L
      next
    }
    
    d_obs <- tryCatch(compute_pair_D(gi, s1, s2), error = function(e) e)
    if (inherits(d_obs, "error") || !is.finite(d_obs)) {
      jost_mat[s1, s2] <- NA_real_
      jost_mat[s2, s1] <- NA_real_
      pair_rows[[row_id]] <- data.frame(
        Site1 = s1,
        Site2 = s2,
        n_site1 = n1,
        n_site2 = n2,
        JostD = NA_real_,
        p_value = NA_real_,
        reason = paste0("D_error: ", if (inherits(d_obs, "error")) d_obs$message else "non-finite D"),
        stringsAsFactors = FALSE
      )
      row_id <- row_id + 1L
      next
    }
    
    jost_mat[s1, s2] <- d_obs
    jost_mat[s2, s1] <- d_obs
    
    pair_gi <- gi[site_factor %in% c(s1, s2), , drop = FALSE]
    perm_vals <- rep(NA_real_, n_perm)
    for (b in seq_len(n_perm)) {
      perm_gi <- pair_gi
      pop(perm_gi) <- sample(pop(pair_gi), size = nInd(pair_gi), replace = FALSE)
      perm_vals[b] <- tryCatch({
        as.numeric(mmod::D(perm_gi))
      }, error = function(e) NA_real_)
    }
    
    valid_perm <- perm_vals[is.finite(perm_vals)]
    p_raw <- if (length(valid_perm) > 0) {
      (sum(valid_perm >= d_obs) + 1) / (length(valid_perm) + 1)
    } else {
      NA_real_
    }
    
    pair_rows[[row_id]] <- data.frame(
      Site1 = s1,
      Site2 = s2,
      n_site1 = n1,
      n_site2 = n2,
      JostD = d_obs,
      p_value = p_raw,
      reason = ifelse(is.na(p_raw), "permutation_failed", NA_character_),
      stringsAsFactors = FALSE
    )
    row_id <- row_id + 1L
  }
}

jost_mat <- (jost_mat + t(jost_mat)) / 2
diag(jost_mat) <- 0

if (!is.matrix(jost_mat) || nrow(jost_mat) == 0 || ncol(jost_mat) == 0) {
  stop("Jost's D matrix is empty.")
}

jost_long <- bind_rows(pair_rows) %>%
  mutate(
    p_adj_bh = p.adjust(p_value, method = "BH"),
    p_adj_bonf = p.adjust(p_value, method = "bonferroni")
  )

matrix_file <- file.path(OUTDIR, "pairwise_jostd_matrix.csv")
write.csv(jost_mat, matrix_file, row.names = TRUE)

long_file <- file.path(OUTDIR, "pairwise_jostd_long.csv")
write.csv(jost_long, long_file, row.names = FALSE)

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

vals <- as.numeric(jost_mat[upper.tri(jost_mat, diag = FALSE)])
vals <- vals[!is.na(vals)]

site_n <- nrow(jost_mat)
if (length(vals) > 0) {
  summary_msg <- paste0(
    "DONE pairwise Jost's D with permutation p-values (n_perm=", n_perm, ").\n",
    "Number of sites: ", site_n, "\n",
    "Jost's D (excluding NA) min/median/max: ",
    format(min(vals), digits = 4), " / ",
    format(stats::median(vals), digits = 4), " / ",
    format(max(vals), digits = 4), "\n",
    "Pairs with non-NA raw p-values: ", sum(!is.na(jost_long$p_value)), " / ", nrow(jost_long), "\n",
    "Saved files:\n",
    "- ", matrix_file, "\n",
    "- ", long_file, "\n",
    "- ", heatmap_file
  )
} else {
  summary_msg <- paste0(
    "DONE pairwise Jost's D with permutation p-values (n_perm=", n_perm, ").\n",
    "Number of sites: ", site_n, "\n",
    "No non-NA pairwise values were available after filtering.\n",
    "Saved files:\n",
    "- ", matrix_file, "\n",
    "- ", long_file, "\n",
    "- ", heatmap_file
  )
}

cat(summary_msg, "\n", sep = "")
