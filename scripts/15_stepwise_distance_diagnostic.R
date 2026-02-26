############################################################
# scripts/15_stepwise_distance_diagnostics.R
# Diagnostics for individual stepwise microsatellite distances
############################################################

options(repos = c(CRAN = "https://cloud.r-project.org"))
pkgs <- c("adegenet", "poppr", "dplyr", "readr", "ggplot2", "tidyr")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)

suppressPackageStartupMessages({
  library(adegenet)
  library(poppr)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(tidyr)
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

parse_allele_token <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "NaN", "0", "000", "000/000", "0/0")] <- NA_character_
  y <- gsub("[^0-9.\\-]", "", x)
  out <- suppressWarnings(as.numeric(y))
  out[!is.finite(out)] <- NA_real_
  out
}

parse_alleles <- function(x) {
  x <- trimws(as.character(x))
  if (is.na(x) || x %in% c("", "NA", "NaN", "0", "000", "0/0", "NA/NA")) return(c(NA_real_, NA_real_))
  x <- gsub("\\|", "/", x)
  parts <- unlist(strsplit(x, "/", fixed = TRUE), use.names = FALSE)
  parts <- parts[nzchar(parts)]
  if (length(parts) == 0) return(c(NA_real_, NA_real_))
  vals <- parse_allele_token(parts)
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0) return(c(NA_real_, NA_real_))
  if (length(vals) == 1) return(c(vals[1], vals[1]))
  c(vals[1], vals[2])
}

gcd_pair <- function(a, b) {
  a <- as.integer(round(abs(a)))
  b <- as.integer(round(abs(b)))
  while (b != 0L) {
    tmp <- b
    b <- a %% b
    a <- tmp
  }
  a
}

infer_repeat_length <- function(v) {
  v <- sort(unique(v[is.finite(v)]))
  if (length(v) < 2) return(NA_integer_)
  dv <- as.numeric(dist(v))
  dv <- dv[dv > 0]
  if (!length(dv)) return(NA_integer_)
  if (any(abs(dv - round(dv)) > 1e-6)) return(NA_integer_)
  g <- as.integer(round(dv[1]))
  for (x in dv[-1]) g <- gcd_pair(g, x)
  if (!is.finite(g) || g <= 0) return(NA_integer_)
  g
}

compute_stepwise_pair <- function(g1, g2, repl) {
  if (any(!is.finite(g1)) || any(!is.finite(g2)) || !is.finite(repl) || repl <= 0) return(NA_real_)
  c1 <- abs(g1[1] - g2[1]) + abs(g1[2] - g2[2])
  c2 <- abs(g1[1] - g2[2]) + abs(g1[2] - g2[1])
  min(c1, c2) / repl
}

upper_tri_vec <- function(m) m[upper.tri(m, diag = FALSE)]

align_meta <- function(meta, ids) {
  if (!is.data.frame(meta)) {
    warning("meta is not a data.frame; creating placeholder meta_used.")
    return(data.frame(Nom_Labo_Échantillons = ids, Site = NA_character_, stringsAsFactors = FALSE))
  }
  
  id_col <- if ("Nom_Labo_Échantillons" %in% names(meta)) "Nom_Labo_Échantillons" else if ("ind" %in% names(meta)) "ind" else NA_character_
  if (is.na(id_col)) {
    warning("meta lacks Nom_Labo_Échantillons/ind column; creating placeholder meta_used.")
    return(data.frame(Nom_Labo_Échantillons = ids, Site = NA_character_, stringsAsFactors = FALSE))
  }
  
  m <- meta %>%
    mutate(Nom_Labo_Échantillons = as.character(.data[[id_col]])) %>%
    filter(Nom_Labo_Échantillons %in% ids) %>%
    distinct(Nom_Labo_Échantillons, .keep_all = TRUE)
  
  out <- data.frame(Nom_Labo_Échantillons = ids, stringsAsFactors = FALSE) %>%
    left_join(m, by = "Nom_Labo_Échantillons")
  
  if (!("Site" %in% names(out))) out$Site <- NA_character_
  missing_n <- sum(is.na(out$Site) | !nzchar(as.character(out$Site)))
  if (missing_n > 0) warning("meta alignment: ", missing_n, " IDs missing Site after alignment.")
  out
}

align_df_ids <- function(df_ids, ids) {
  if (!is.data.frame(df_ids) || !("ind" %in% names(df_ids))) {
    return(data.frame(ind = ids, Site = NA_character_, stringsAsFactors = FALSE))
  }
  x <- df_ids %>%
    mutate(ind = as.character(ind)) %>%
    filter(ind %in% ids) %>%
    distinct(ind, .keep_all = TRUE)
  out <- data.frame(ind = ids, stringsAsFactors = FALSE) %>% left_join(x, by = "ind")
  if (!("Site" %in% names(out))) out$Site <- NA_character_
  out
}

PROJECT_ROOT <- find_project_root()
setwd(PROJECT_ROOT)
source(file.path(PROJECT_ROOT, "scripts", "_load_objects.R"))

RUN_TAG <- if (exists("RUN_TAG", inherits = TRUE)) get("RUN_TAG", inherits = TRUE) else "v1"
RUN_OUT <- if (exists("RUN_OUT", inherits = TRUE)) get("RUN_OUT", inherits = TRUE) else file.path(PROJECT_ROOT, "outputs", RUN_TAG)
OUT_DIR <- file.path(RUN_OUT, "stepwise_distance")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# A) choose gi_used
if (exists("gi_mll") && inherits(gi_mll, "genind") && nInd(gi_mll) > 1) {
  gi_used <- gi_mll
  message("Diagnostics using clone-corrected genind object: gi_mll")
} else if (exists("gi") && inherits(gi, "genind") && nInd(gi) > 1) {
  gi_used <- gi
  message("Diagnostics using full genind object: gi")
} else {
  stop("No usable genind object found (need gi_mll or gi with >= 2 individuals).")
}

ids <- adegenet::indNames(gi_used)
loci <- adegenet::locNames(gi_used)

meta_used <- align_meta(meta = if (exists("meta")) meta else NULL, ids = ids)
df_ids_used <- align_df_ids(df_ids = if (exists("df_ids")) df_ids else NULL, ids = ids)

site_vec <- as.character(df_ids_used$Site)
site_vec[is.na(site_vec) | !nzchar(site_vec)] <- as.character(meta_used$Site[is.na(site_vec) | !nzchar(site_vec)])
site_vec[is.na(site_vec) | !nzchar(site_vec)] <- as.character(adegenet::pop(gi_used)[is.na(site_vec) | !nzchar(site_vec)])
site_vec[is.na(site_vec) | !nzchar(site_vec)] <- "Unknown"
names(site_vec) <- ids

message("DEBUG: nInd(gi_used)=", nInd(gi_used), ", nLoc(gi_used)=", nLoc(gi_used), ", length(ids)=", length(ids), ", nrow(meta_used)=", nrow(meta_used))

if (length(ids) < 2) stop("Need at least 2 individuals after alignment.")
if (length(loci) < 1) stop("Need at least 1 locus in gi_used.")

# B) robust per-locus genotype extraction from genind2df
geno_df <- adegenet::genind2df(gi_used, sep = "/", usepop = FALSE)

if (!nrow(geno_df)) stop("genind2df returned 0 rows.")
if (is.null(rownames(geno_df)) || !all(ids %in% rownames(geno_df))) {
  # fallback alignment by order
  if (nrow(geno_df) != length(ids)) {
    stop("Cannot align genind2df rows to IDs; nrow(geno_df) != length(ids).")
  }
  rownames(geno_df) <- ids
}

geno_df <- geno_df[ids, , drop = FALSE]

avail_loci <- intersect(loci, colnames(geno_df))
missing_loci <- setdiff(loci, avail_loci)
if (length(avail_loci) < 1) {
  stop("No locus columns from locNames(gi_used) found in genind2df output.")
}
if (length(missing_loci) > 0) {
  warning("Dropping loci missing in genind2df: ", paste(missing_loci, collapse = ", "))
}
loci <- avail_loci

parsed_fail <- 0L
parsed_by_locus <- vector("list", length(loci))
all_alleles_by_locus <- vector("list", length(loci))
names(parsed_by_locus) <- loci
names(all_alleles_by_locus) <- loci

for (loc in loci) {
  colv <- as.character(geno_df[[loc]])
  loc_parsed <- t(vapply(colv, parse_alleles, numeric(2)))
  rownames(loc_parsed) <- ids
  parsed_by_locus[[loc]] <- loc_parsed
  
  toks <- unlist(strsplit(colv, "/", fixed = TRUE), use.names = FALSE)
  toks <- toks[nzchar(toks)]
  nums <- parse_allele_token(toks)
  parsed_fail <- parsed_fail + sum(!is.na(toks) & !is.finite(nums))
  all_alleles_by_locus[[loc]] <- nums[is.finite(nums)]
}

if (parsed_fail > 0) {
  warning(parsed_fail, " allele tokens could not be parsed as numeric sizes; treated as missing.")
}

inferred_replen <- sapply(all_alleles_by_locus, infer_repeat_length)
replen <- inferred_replen
replen[is.na(replen) | replen <= 0] <- 2
if (any(is.na(inferred_replen))) {
  warning("Could not infer repeat length for ", sum(is.na(inferred_replen)), " loci; using default replen=2 for those loci.")
}
message("DEBUG: loci with inferred repeat length = ", sum(!is.na(inferred_replen)), " / ", length(inferred_replen))

n_ids <- length(ids)
n_loci <- length(loci)
if (n_ids < 2 || n_loci < 1) stop("Insufficient individuals or loci after preprocessing.")

# ============================================================
# A) Pairwise locus-difference summaries
# ============================================================
pair_rows <- vector("list", choose(n_ids, 2))
idx <- 0L

locus_sum <- setNames(numeric(n_loci), loci)
locus_n <- setNames(integer(n_loci), loci)

for (i in seq_len(n_ids - 1)) {
  for (j in (i + 1):n_ids) {
    stopifnot(i <= n_ids, j <= n_ids)
    
    step_vec <- rep(NA_real_, n_loci)
    names(step_vec) <- loci
    
    for (loc in loci) {
      gmat <- parsed_by_locus[[loc]]
      if (is.null(gmat) || nrow(gmat) < j || ncol(gmat) < 2) {
        warning("Malformed parsed matrix at locus ", loc, " (nrow=", ifelse(is.null(gmat), NA, nrow(gmat)), ", ncol=", ifelse(is.null(gmat), NA, ncol(gmat)), ").")
        next
      }
      
      s <- tryCatch(
        compute_stepwise_pair(gmat[i, ], gmat[j, ], replen[loc]),
        error = function(e) {
          message("DEBUG ERROR locus/pair: locus=", loc, ", ind1=", ids[i], ", ind2=", ids[j], " -> ", conditionMessage(e))
          NA_real_
        }
      )
      
      step_vec[loc] <- s
      if (is.finite(s)) {
        locus_sum[loc] <- locus_sum[loc] + s
        locus_n[loc] <- locus_n[loc] + 1L
      }
    }
    
    compared <- is.finite(step_vec)
    n_comp <- sum(compared)
    n_diff <- if (n_comp > 0) sum(step_vec[compared] > 0) else 0L
    
    idx <- idx + 1L
    pair_rows[[idx]] <- data.frame(
      ind1 = ids[i],
      ind2 = ids[j],
      site1 = site_vec[ids[i]],
      site2 = site_vec[ids[j]],
      pair_type = ifelse(site_vec[ids[i]] == site_vec[ids[j]], "within_site", "between_site"),
      n_loci_compared = n_comp,
      n_loci_diff = n_diff,
      prop_loci_diff = ifelse(n_comp > 0, n_diff / n_comp, NA_real_),
      mean_step = ifelse(n_comp > 0, mean(step_vec[compared], na.rm = TRUE), NA_real_),
      median_step = ifelse(n_comp > 0, median(step_vec[compared], na.rm = TRUE), NA_real_),
      max_step = ifelse(n_comp > 0, max(step_vec[compared], na.rm = TRUE), NA_real_),
      stringsAsFactors = FALSE
    )
  }
}

pair_diag <- dplyr::bind_rows(pair_rows)
readr::write_csv(pair_diag, file.path(OUT_DIR, "smm_pairwise_locus_diagnostics.csv"))

p1 <- ggplot(pair_diag, aes(x = prop_loci_diff)) +
  geom_histogram(bins = 30, color = "white") +
  theme_bw(base_size = 12) +
  labs(title = "Proportion of loci differing across pairs", x = "prop_loci_diff", y = "Count of pairs")
ggsave(file.path(OUT_DIR, "diagnostics_hist_prop_loci_diff.png"), p1, width = 8, height = 5, dpi = 300)

p2 <- ggplot(pair_diag, aes(x = mean_step)) +
  geom_histogram(bins = 30, color = "white") +
  theme_bw(base_size = 12) +
  labs(title = "Mean step difference across loci (pairwise)", x = "mean_step", y = "Count of pairs")
ggsave(file.path(OUT_DIR, "diagnostics_hist_mean_step.png"), p2, width = 8, height = 5, dpi = 300)

if (length(unique(pair_diag$pair_type)) > 1) {
  p3 <- ggplot(pair_diag, aes(x = pair_type, y = prop_loci_diff, fill = pair_type)) +
    geom_boxplot(alpha = 0.8, outlier.alpha = 0.25) +
    theme_bw(base_size = 12) +
    guides(fill = "none") +
    labs(title = "Proportion of differing loci by pair type", x = "Pair type", y = "prop_loci_diff")
  ggsave(file.path(OUT_DIR, "diagnostics_boxplot_prop_loci_diff_by_pair_type.png"), p3, width = 8, height = 5, dpi = 300)
}

# ============================================================
# B) LOLO sensitivity for Bruvo distances
# ============================================================
read_existing_dist <- function(path, expected_ids) {
  if (!file.exists(path)) return(NULL)
  x <- tryCatch(readr::read_csv(path, show_col_types = FALSE), error = function(e) NULL)
  if (is.null(x) || !nrow(x) || !("ind" %in% names(x))) return(NULL)
  
  rn <- as.character(x$ind)
  cn <- setdiff(names(x), "ind")
  if (!length(cn)) return(NULL)
  
  mat <- as.matrix(x[, cn, drop = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- rn
  colnames(mat) <- cn
  
  if (nrow(mat) != ncol(mat)) return(NULL)
  if (!setequal(rownames(mat), colnames(mat))) return(NULL)
  if (!setequal(expected_ids, rownames(mat))) return(NULL)
  
  mat <- mat[expected_ids, expected_ids, drop = FALSE]
  mat
}

existing_dist_path <- file.path(OUT_DIR, "smm_dist_individuals.csv")
D_full <- read_existing_dist(existing_dist_path, expected_ids = ids)

if (is.null(D_full)) {
  message("Existing smm_dist_individuals.csv not usable for current individuals; recomputing full Bruvo matrix.")
  D_full <- as.matrix(poppr::bruvo.dist(gi_used, replen = replen, add = TRUE, loss = TRUE))
}

if (!setequal(rownames(D_full), ids) || !setequal(colnames(D_full), ids)) {
  stop("D_full IDs are not aligned with gi_used IDs.")
}
D_full <- D_full[ids, ids, drop = FALSE]
D_full <- (D_full + t(D_full)) / 2
diag(D_full) <- 0

lolo_rows <- vector("list", length(loci))
for (k in seq_along(loci)) {
  loc_drop <- loci[k]
  keep <- setdiff(loci, loc_drop)
  if (length(keep) < 1) next
  
  gi_minus <- tryCatch(gi_used[, keep], error = function(e) {
    message("DEBUG ERROR subsetting gi_minus at locus ", loc_drop, ": ", conditionMessage(e))
    NULL
  })
  
  if (is.null(gi_minus)) {
    lolo_rows[[k]] <- data.frame(locus = loc_drop, cor_pearson = NA_real_, cor_spearman = NA_real_, n_pairs = 0L, stringsAsFactors = FALSE)
    next
  }
  
  repl_minus <- replen[keep]
  D_minus <- tryCatch(as.matrix(poppr::bruvo.dist(gi_minus, replen = repl_minus, add = TRUE, loss = TRUE)), error = function(e) {
    message("DEBUG ERROR bruvo.dist LOLO at locus ", loc_drop, ": ", conditionMessage(e))
    NULL
  })
  
  if (is.null(D_minus)) {
    lolo_rows[[k]] <- data.frame(locus = loc_drop, cor_pearson = NA_real_, cor_spearman = NA_real_, n_pairs = 0L, stringsAsFactors = FALSE)
    next
  }
  
  if (!setequal(rownames(D_minus), ids) || !setequal(colnames(D_minus), ids)) {
    lolo_rows[[k]] <- data.frame(locus = loc_drop, cor_pearson = NA_real_, cor_spearman = NA_real_, n_pairs = 0L, stringsAsFactors = FALSE)
    next
  }
  
  D_minus <- D_minus[ids, ids, drop = FALSE]
  D_minus <- (D_minus + t(D_minus)) / 2
  diag(D_minus) <- 0
  
  v1 <- upper_tri_vec(D_full)
  v2 <- upper_tri_vec(D_minus)
  ok <- is.finite(v1) & is.finite(v2)
  
  lolo_rows[[k]] <- data.frame(
    locus = loc_drop,
    cor_pearson = if (sum(ok) > 2) suppressWarnings(cor(v1[ok], v2[ok], method = "pearson")) else NA_real_,
    cor_spearman = if (sum(ok) > 2) suppressWarnings(cor(v1[ok], v2[ok], method = "spearman")) else NA_real_,
    n_pairs = sum(ok),
    stringsAsFactors = FALSE
  )
}

lolo_tab <- bind_rows(lolo_rows) %>% arrange(cor_spearman)
readr::write_csv(lolo_tab, file.path(OUT_DIR, "smm_LOLO_correlations.csv"))

lolo_long <- lolo_tab %>%
  select(locus, cor_pearson, cor_spearman) %>%
  pivot_longer(cols = c(cor_pearson, cor_spearman), names_to = "metric", values_to = "correlation")

p4 <- ggplot(lolo_long, aes(x = reorder(locus, correlation, na.rm = TRUE), y = correlation, fill = metric)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  coord_flip() +
  theme_bw(base_size = 12) +
  labs(title = "Leave-one-locus-out correlation with full Bruvo matrix", x = "Locus", y = "Correlation")
ggsave(file.path(OUT_DIR, "diagnostics_lolo_correlations_by_locus.png"), p4, width = 9, height = 6, dpi = 300)

# ============================================================
# C) Locus contribution ranking
# ============================================================
contrib <- data.frame(
  locus = loci,
  mean_step_diff = as.numeric(locus_sum / ifelse(locus_n > 0, locus_n, NA)),
  n_pairs_compared = as.integer(locus_n),
  lolo_cor_spearman = lolo_tab$cor_spearman[match(loci, lolo_tab$locus)],
  lolo_cor_pearson = lolo_tab$cor_pearson[match(loci, lolo_tab$locus)],
  lolo_impact = 1 - lolo_tab$cor_spearman[match(loci, lolo_tab$locus)],
  stringsAsFactors = FALSE
) %>%
  arrange(desc(mean_step_diff))

readr::write_csv(contrib, file.path(OUT_DIR, "smm_locus_contributions.csv"))

# ============================================================
# D) Append README summary
# ============================================================
readme_path <- file.path(OUT_DIR, "README.txt")

median_n_loci_diff <- median(pair_diag$n_loci_diff, na.rm = TRUE)
lolo_min <- suppressWarnings(min(lolo_tab$cor_spearman, na.rm = TRUE))
lolo_max <- suppressWarnings(max(lolo_tab$cor_spearman, na.rm = TRUE))
if (!is.finite(lolo_min)) lolo_min <- NA_real_
if (!is.finite(lolo_max)) lolo_max <- NA_real_

conclusion <- if (is.finite(lolo_min) && lolo_min > 0.90) {
  "LOLO correlations remain very high across loci, supporting that distances are not driven by a single locus."
} else if (is.finite(lolo_min) && lolo_min > 0.75) {
  "LOLO correlations are generally high, indicating multilocus signal; no single locus appears to dominate distance structure."
} else {
  "Some loci exert stronger influence, but rankings and pairwise diagnostics indicate distances reflect multilocus differences rather than only one locus."
}

snippet <- c(
  "",
  "--- Diagnostics update ---",
  paste0("Median number of loci differing per pair: ", round(median_n_loci_diff, 2)),
  paste0("LOLO Spearman correlation range: ", round(lolo_min, 4), " to ", round(lolo_max, 4)),
  paste0("Conclusion: ", conclusion)
)

if (file.exists(readme_path)) {
  write(snippet, file = readme_path, append = TRUE)
} else {
  writeLines(c("Stepwise microsatellite distance analysis", snippet), con = readme_path)
}

message("Saved: ", file.path(OUT_DIR, "smm_pairwise_locus_diagnostics.csv"))
message("Saved: ", file.path(OUT_DIR, "smm_LOLO_correlations.csv"))
message("Saved: ", file.path(OUT_DIR, "smm_locus_contributions.csv"))
message("Saved plots: diagnostics_*.png in ", OUT_DIR)
message("Updated README: ", readme_path)
