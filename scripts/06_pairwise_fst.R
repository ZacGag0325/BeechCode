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

set.seed(123)

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

normalize_site <- function(x) toupper(trimws(as.character(x)))

setwd(find_project_root())
source(file.path("scripts", "_load_objects.R"))

OUTDIR <- file.path(RUN_OUT, "fst_only")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

if (!inherits(gi, "genind")) stop("Object 'gi' must be a genind object.")
if (nPop(gi) < 2) stop("Need at least 2 populations (sites) to compute pairwise Jost's D.")

site_factor <- as.factor(pop(gi))
levels(site_factor) <- normalize_site(levels(site_factor))
pop(gi) <- site_factor
site_levels <- levels(site_factor)

min_n_per_site <- 2L
n_perm <- 999L
site_counts <- table(site_factor)

if (any(site_counts < min_n_per_site)) {
  warning(
    "Sites with n < ", min_n_per_site,
    " may return NA pairwise values: ",
    paste(names(site_counts)[site_counts < min_n_per_site], collapse = ", ")
  )
}

get_pair_reason <- function(genind_pair, s1, s2, min_n = min_n_per_site) {
  n1 <- sum(pop(genind_pair) == s1)
  n2 <- sum(pop(genind_pair) == s2)
  if (n1 < min_n || n2 < min_n) {
    return(paste0("too_few_individuals (", s1, "=", n1, ", ", s2, "=", n2, ")"))
  }
  
  pair_tab <- adegenet::tab(genind_pair, NA.method = "mean")
  if (!is.matrix(pair_tab) || nrow(pair_tab) == 0 || ncol(pair_tab) == 0) {
    return("empty_genotype_table")
  }
  
  loc_fac <- adegenet::locFac(genind_pair)
  loci <- levels(loc_fac)
  if (length(loci) == 0) return("no_loci")
  
  n_poly <- sum(vapply(loci, function(loc) {
    idx <- which(loc_fac == loc)
    if (length(idx) == 0) return(FALSE)
    alleles_present <- colSums(pair_tab[, idx, drop = FALSE], na.rm = TRUE) > 0
    sum(alleles_present) >= 2
  }, logical(1)))
  
  if (n_poly == 0) return("all_loci_monomorphic")
  NA_character_
}

# A) Core pairwise Jost's D matrix using mmod::pairwise_D
jost_mat_raw <- tryCatch(
  mmod::pairwise_D(gi),
  error = function(e) e
)
if (inherits(jost_mat_raw, "error")) {
  stop("Failed to compute pairwise Jost's D with mmod::pairwise_D: ", jost_mat_raw$message)
}

jost_mat <- as.matrix(jost_mat_raw)
if (is.null(rownames(jost_mat)) || is.null(colnames(jost_mat))) {
  stop("pairwise_D did not return named matrix dimensions.")
}

rownames(jost_mat) <- normalize_site(rownames(jost_mat))
colnames(jost_mat) <- normalize_site(colnames(jost_mat))

# align to pop(gi) site levels and guard against missing rows/cols
aligned <- matrix(
  NA_real_,
  nrow = length(site_levels),
  ncol = length(site_levels),
  dimnames = list(site_levels, site_levels)
)
diag(aligned) <- 0
common_sites <- intersect(site_levels, intersect(rownames(jost_mat), colnames(jost_mat)))
if (length(common_sites) >= 2) {
  aligned[common_sites, common_sites] <- jost_mat[common_sites, common_sites, drop = FALSE]
}
jost_mat <- aligned
jost_mat <- (jost_mat + t(jost_mat)) / 2
diag(jost_mat) <- 0

# B) Long table with permutation p-values + reasons
pair_idx <- which(upper.tri(jost_mat, diag = FALSE), arr.ind = TRUE)
if (nrow(pair_idx) == 0) stop("No site pairs found to summarize.")

pair_rows <- vector("list", nrow(pair_idx))

for (k in seq_len(nrow(pair_idx))) {
  s1 <- rownames(jost_mat)[pair_idx[k, 1]]
  s2 <- colnames(jost_mat)[pair_idx[k, 2]]
  d_obs <- jost_mat[s1, s2]
  
  pair_gi <- gi[pop(gi) %in% c(s1, s2), , drop = FALSE]
  pop(pair_gi) <- droplevels(as.factor(pop(pair_gi)))
  
  reason <- get_pair_reason(pair_gi, s1, s2)
  
  p_raw <- NA_real_
  if (is.na(reason) && is.finite(d_obs)) {
    perm_vals <- rep(NA_real_, n_perm)
    for (b in seq_len(n_perm)) {
      perm_gi <- pair_gi
      pop(perm_gi) <- sample(pop(pair_gi), size = nInd(pair_gi), replace = FALSE)
      perm_vals[b] <- tryCatch(as.numeric(mmod::D(perm_gi)), error = function(e) NA_real_)
    }
    valid_perm <- perm_vals[is.finite(perm_vals)]
    if (length(valid_perm) > 0) {
      p_raw <- (sum(valid_perm >= d_obs) + 1) / (length(valid_perm) + 1)
    } else {
      reason <- "permutation_failed"
    }
  }
  
  if (!is.finite(d_obs) && is.na(reason)) {
    reason <- "pairwise_D_returned_NA"
  }
  
  pair_rows[[k]] <- tibble(
    site1 = s1,
    site2 = s2,
    D = as.numeric(d_obs),
    p_value = as.numeric(p_raw),
    reason = reason
  )
}

jost_long <- bind_rows(pair_rows) %>%
  mutate(
    p_adj_bh = p.adjust(p_value, method = "BH"),
    p_adj_bonf = p.adjust(p_value, method = "bonferroni")
  ) %>%
  select(site1, site2, D, p_value, p_adj_bh, p_adj_bonf, reason)

# C) Outputs
matrix_file <- file.path(OUTDIR, "pairwise_jostd_matrix.csv")
write.csv(jost_mat, matrix_file, row.names = TRUE)

long_file <- file.path(OUTDIR, "pairwise_jostd_long.csv")
write.csv(jost_long, long_file, row.names = FALSE)

heatmap_png <- file.path(OUTDIR, "pairwise_jostd_heatmap.png")
heatmap_pdf <- file.path(OUTDIR, "pairwise_jostd_heatmap.pdf")

heat_df <- as.data.frame(jost_mat) %>%
  rownames_to_column("site1") %>%
  pivot_longer(-site1, names_to = "site2", values_to = "D")

p_heat <- ggplot(heat_df, aes(x = site1, y = site2, fill = D)) +
  geom_tile(color = "white", linewidth = 0.1) +
  scale_fill_gradient(low = "#f7fbff", high = "#08306b", na.value = "grey85") +
  coord_equal() +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank()) +
  labs(title = "Pairwise Jost's D by Site", x = "Site", y = "Site", fill = "Jost's D")

ggsave(heatmap_png, p_heat, width = 8.5, height = 7.5, dpi = 300)
ggsave(heatmap_pdf, p_heat, width = 8.5, height = 7.5)

vals <- as.numeric(jost_mat[upper.tri(jost_mat, diag = FALSE)])
vals_ok <- vals[is.finite(vals)]

cat(
  "DONE pairwise Jost's D (n_perm=", n_perm, ").\n",
  "Sites: ", length(site_levels), "\n",
  "Pairs with finite D: ", length(vals_ok), " / ", length(vals), "\n",
  "Pairs with finite p-values: ", sum(is.finite(jost_long$p_value)), " / ", nrow(jost_long), "\n",
  "Outputs:\n",
  "- ", matrix_file, "\n",
  "- ", long_file, "\n",
  "- ", heatmap_png, "\n",
  "- ", heatmap_pdf, "\n",
  sep = ""
)
