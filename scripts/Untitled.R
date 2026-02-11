############################################################
# FULL PIPELINE (QC + MLG/MLL + HWE + clonality by SITE)
# + HWE FDR correction + HWE summary table
# + NJ tree (Bruvo distance) + “boxed”/highlighted MLL clusters (like the paper figure)
# + NJ tree BY SITE (Nei + optional chord distance)
# + Optional spatial genetic structure (Mantel) if coords + packages exist
#
# ADD-ON (robust when clonality is rare):
# + Clone-pair geographic distances (if Latitude/Longitude)
# + Permutation test: are clone pairs closer than random pairs?
# + (Exploratory) clonality ~ latitude correlation at site level
############################################################

### 0) PACKAGES ----
suppressPackageStartupMessages({
  library(adegenet)
  library(poppr)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(hierfstat)
  library(pegas)     # HWE exact tests
  library(ape)       # NJ tree
})

# Optional packages
suppressPackageStartupMessages({
  ok_vegan  <- requireNamespace("vegan", quietly = TRUE)
  ok_geo    <- requireNamespace("geosphere", quietly = TRUE)
  ok_ggtree <- requireNamespace("ggtree", quietly = TRUE)  # for “boxes”/highlights
  ok_ggplot <- requireNamespace("ggplot2", quietly = TRUE)
})

### 0B) SETTINGS ----
save_outputs <- TRUE
missing_threshold <- 0.35
chosen_thresh <- 0.09               # MLL Bruvo threshold

# HWE options
HWE_B <- 10000
FDR_scope <- "site"

# NJ tree plot options
tree_tip_cex <- 0.4
min_group_for_box <- 2
tree_out_basic_pdf <- "NJ_tree_Bruvo_basic.pdf"
tree_out_boxed_pdf <- "NJ_tree_Bruvo_MLL_boxes.pdf"

# NEW: NJ by site outputs
tree_site_nei_pdf   <- "NJ_tree_by_site_Nei.pdf"
tree_site_nei_nwk   <- "NJ_tree_by_site_Nei.nwk"
tree_site_nei_csv   <- "NJ_by_site_Nei_distance_matrix.csv"

tree_site_chord_pdf <- "NJ_tree_by_site_Chord.pdf"
tree_site_chord_nwk <- "NJ_tree_by_site_Chord.nwk"
tree_site_chord_csv <- "NJ_by_site_Chord_distance_matrix.csv"

# NEW: rare-clonality robust analyses
n_perm_clone_distance <- 999
min_clone_pairs_for_test <- 3     # if fewer, we only report descriptives
min_sites_for_latitude_test <- 3  # Spearman across sites

# Files
geno_file   <- "poppr.xlsx"
field_file  <- "donnees_modifiees_west_summer2024 copie.xlsx"
field_sheet <- "genetique"

### 1) DOSSIER DE TRAVAIL ----
setwd("/Users/zacharygagnon/Desktop/R data workshop")
getwd()
list.files()

### 2) READ GENOTYPE DATA ----
geno_raw <- read_xlsx(geno_file)

geno_raw <- geno_raw %>%
  select(where(~ !all(is.na(.)))) %>%
  select(-matches("^Unnamed"))

names(geno_raw) <- trimws(names(geno_raw))

req_geno <- c("Nom_Labo_Échantillons", "Numéro_Population")
if (!all(req_geno %in% names(geno_raw))) {
  stop(paste0("poppr.xlsx must contain: ", paste(req_geno, collapse = ", ")))
}

geno_raw$Nom_Labo_Échantillons <- as.character(geno_raw$Nom_Labo_Échantillons)
geno <- geno_raw

cat("\nPop counts (before QC):\n")
print(table(geno$Numéro_Population))

### 3) READ FIELD DATA (LAB ID + SITE) ----
field <- read_xlsx(field_file, sheet = field_sheet)
names(field) <- trimws(names(field))

req_field <- c("Nom_Labo_Échantillons", "Site")
if (!all(req_field %in% names(field))) {
  stop(paste0(
    "Field sheet '", field_sheet, "' must contain: ",
    paste(req_field, collapse = ", "),
    "\nColumns found: ", paste(names(field), collapse = ", ")
  ))
}

field <- field %>%
  mutate(
    Nom_Labo_Échantillons = as.character(Nom_Labo_Échantillons),
    Site = as.character(Site)
  )

# Optional columns
has_lat <- "Latitude"  %in% names(field)
has_lon <- "Longitude" %in% names(field)
has_dbh <- "Dhp_tige"  %in% names(field)

if (has_lat) field$Latitude  <- as.numeric(field$Latitude)
if (has_lon) field$Longitude <- as.numeric(field$Longitude)
if (has_dbh) field$Dhp_tige  <- as.numeric(field$Dhp_tige)

### 4) JOIN SITE INFO INTO GENO ----
field_keep <- c("Nom_Labo_Échantillons", "Site")
if (has_lat) field_keep <- c(field_keep, "Latitude")
if (has_lon) field_keep <- c(field_keep, "Longitude")
if (has_dbh) field_keep <- c(field_keep, "Dhp_tige")

field2 <- field %>%
  select(all_of(field_keep)) %>%
  distinct(Nom_Labo_Échantillons, .keep_all = TRUE)

geno <- geno %>%
  left_join(field2, by = "Nom_Labo_Échantillons")

n_missing_site <- sum(is.na(geno$Site))
if (n_missing_site > 0) {
  cat("\nWARNING: Some samples in poppr.xlsx did not match the field file (Site is NA).\n")
  cat("Count missing Site:", n_missing_site, "\n")
  cat("Examples:\n")
  print(head(geno$Nom_Labo_Échantillons[is.na(geno$Site)], 10))
  stop("Fix the IDs so Nom_Labo_Échantillons matches in BOTH files (exact same text).")
}

cat("\nSites found in joined data:\n")
print(table(geno$Site))

### 5) DEFINE ALLELE COLUMNS (_1/_2) ----
allele_cols <- which(grepl("_(1|2)$", names(geno)))
if (length(allele_cols) == 0) stop("No allele columns detected. Must end with _1 and _2.")
if (length(allele_cols) %% 2 != 0) stop("Allele columns not even. Check _1/_2 headers.")

locus_base <- sub("_(1|2)$", "", names(geno)[allele_cols])
bad <- table(locus_base)[table(locus_base) != 2]
if (length(bad) > 0) {
  print(bad)
  stop("At least one locus does not have exactly 2 columns (_1 and _2). Fix Excel headers.")
}

locus_starts <- sort(allele_cols[grepl("_1$", names(geno)[allele_cols])])

### 6) BUILD 1-COLUMN-PER-LOCUS TABLE ----
geno_locus <- data.frame(row.names = geno$Nom_Labo_Échantillons)

for (i in seq_along(locus_starts)) {
  col1_idx <- locus_starts[i]
  col2_idx <- col1_idx + 1
  
  locus_name <- sub("_1$", "", names(geno)[col1_idx])
  
  a1 <- as.character(geno[[col1_idx]])
  a2 <- as.character(geno[[col2_idx]])
  
  geno_locus[[locus_name]] <- paste(a1, a2, sep = "/")
}

geno_locus[] <- lapply(geno_locus, function(col) {
  col[grepl("NA", col, fixed = TRUE)] <- NA
  col
})

names(geno_locus) <- gsub("\\.+", "_", names(geno_locus))

cat("\nLoci detected:\n")
print(names(geno_locus))

### 7) QC FILTER (>35% missing loci removed) ----
missing_prop <- rowMeans(is.na(geno_locus))
cat("\nMissingness summary:\n")
print(summary(missing_prop))

keep_inds <- missing_prop <= missing_threshold

cat("\nIndividuals before QC:", nrow(geno_locus), "\n")
cat("Individuals kept (<= ", missing_threshold * 100, "% missing): ", sum(keep_inds), "\n", sep = "")

geno_locus <- geno_locus[keep_inds, , drop = FALSE]
geno       <- geno[keep_inds, , drop = FALSE]
stopifnot(nrow(geno_locus) == nrow(geno))

cat("\nSite counts (after QC):\n")
print(table(geno$Site))

### 8) BUILD genind ----
gen <- df2genind(
  geno_locus,
  ploidy    = 2,
  ind.names = rownames(geno_locus),
  pop       = as.factor(geno$Numéro_Population),
  type      = "codom",
  sep       = "/",
  ncode     = 3
)

cat("\nGenind summary:\n")
print(gen)
print(summary(gen))

# Bruvo requires numeric alleles
aa <- alleles(gen)
bad_alleles <- lapply(aa, function(v) v[is.na(as.numeric(v))])
bad_alleles <- bad_alleles[sapply(bad_alleles, length) > 0]
if (length(bad_alleles) > 0) {
  print(bad_alleles)
  stop("Some alleles are not numeric. Fix allele values (Bruvo requires numeric).")
} else {
  message("OK: all alleles are numeric (Bruvo compatible).")
}

### 9) MLG ----
gc_mlg <- as.genclone(gen)
geno$MLG <- as.integer(mlg.vector(gc_mlg))

cat("\nUnique MLGs:", length(unique(mlg.vector(gc_mlg))), "out of", nInd(gen), "individuals\n")
cat("\nClonality stats (MLG):\n")
print(poppr(gc_mlg))

### 10) MLL via BRUVO DISTANCE ----
replen <- rep(2, nLoc(gen)); names(replen) <- locNames(gen)

gc_mll <- gc_mlg
mlg.filter(
  gc_mll,
  distance  = bruvo.dist,
  replen    = replen,
  algorithm = "farthest_neighbor"
) <- chosen_thresh

geno$MLL <- as.integer(mll(gc_mll))

cat("\nMLL threshold used:", chosen_thresh, "\n")
cat("nmll before filter:", nmll(gc_mlg), "\n")
cat("nmll after  filter:", nmll(gc_mll), "\n")

cat("\nClonality stats (MLL):\n")
print(poppr(gc_mll))

### 11) CLONE LISTS (MLG & MLL) ----
mlg_full_list <- split(
  data.frame(Individual = indNames(gc_mlg), Population = pop(gc_mlg)),
  mlg.vector(gc_mlg)
)
mlg_clone_list <- mlg_full_list[sapply(mlg_full_list, nrow) > 1]

mll_full_list <- split(
  data.frame(Individual = indNames(gc_mll), Population = pop(gc_mlg)),
  mll(gc_mll)
)
mll_clone_list <- mll_full_list[sapply(mll_full_list, nrow) > 1]

cat("\n# shared clones (MLG):", length(mlg_clone_list), "\n")
cat("# shared clones (MLL):", length(mll_clone_list), "\n")

### 12) HWE PER SITE × LOCUS (EXACT TEST) + FDR + SUMMARY ----
subset_genind_by_inds <- function(gen_obj, inds) gen_obj[inds, , drop = TRUE]

hw_to_df <- function(hw_obj) {
  if (is.numeric(hw_obj) && !is.null(names(hw_obj)) && length(hw_obj) > 0) {
    return(data.frame(Locus = names(hw_obj), P_value = as.numeric(hw_obj), row.names = NULL))
  }
  if (is.matrix(hw_obj) || is.data.frame(hw_obj)) {
    hw_tbl <- as.data.frame(hw_obj)
    if (!is.null(rownames(hw_tbl)) && nrow(hw_tbl) > 0) {
      hw_tbl$Locus <- rownames(hw_tbl)
      rownames(hw_tbl) <- NULL
    }
    pcol_candidates <- c("p.value", "pvalue", "p", "pr.exact", "pr.exact.")
    pcol <- names(hw_tbl)[tolower(names(hw_tbl)) %in% pcol_candidates]
    if (length(pcol) == 0) {
      num_cols <- names(hw_tbl)[sapply(hw_tbl, is.numeric)]
      if (length(num_cols) == 1) {
        return(hw_tbl %>% transmute(Locus = Locus, P_value = .data[[num_cols]]))
      }
      vals <- as.numeric(unlist(hw_tbl))
      return(data.frame(Locus = seq_along(vals), P_value = vals))
    } else {
      pcol <- pcol[1]
      return(hw_tbl %>% transmute(
        Locus = if ("Locus" %in% names(hw_tbl)) Locus else seq_len(nrow(hw_tbl)),
        P_value = .data[[pcol]]
      ))
    }
  }
  vals <- tryCatch(as.numeric(hw_obj), error = function(e) numeric(0))
  if (length(vals) == 0) return(data.frame(Locus = character(0), P_value = numeric(0)))
  data.frame(Locus = seq_along(vals), P_value = vals)
}

sites <- sort(unique(geno$Site))

hwe_one_site <- function(site_code) {
  inds <- geno$Nom_Labo_Échantillons[geno$Site == site_code]
  gsite <- subset_genind_by_inds(gen, inds)
  
  if (nInd(gsite) < 3) {
    return(data.frame(Site = site_code, Locus = NA_character_, P_value = NA_real_))
  }
  
  hw <- hw.test(gsite, B = HWE_B)
  hw_df <- hw_to_df(hw)
  
  if (nrow(hw_df) == 0) {
    return(data.frame(Site = site_code, Locus = NA_character_, P_value = NA_real_))
  }
  
  hw_df$Site <- site_code
  hw_df %>% select(Site, Locus, P_value)
}

hwe_results <- bind_rows(lapply(sites, hwe_one_site)) %>%
  mutate(Significant_0.05 = P_value < 0.05) %>%
  arrange(Site, Locus)

# FDR correction (BH)
if (FDR_scope == "site") {
  hwe_results <- hwe_results %>%
    group_by(Site) %>%
    mutate(
      P_FDR = p.adjust(P_value, method = "BH"),
      Significant_FDR_0.05 = P_FDR < 0.05
    ) %>%
    ungroup()
} else if (FDR_scope == "global") {
  hwe_results <- hwe_results %>%
    mutate(
      P_FDR = p.adjust(P_value, method = "BH"),
      Significant_FDR_0.05 = P_FDR < 0.05
    )
}

cat("\nHWE results (head):\n")
print(head(hwe_results, 30))

hwe_summary <- hwe_results %>%
  group_by(Site) %>%
  summarise(
    n_tests = sum(!is.na(P_value)),
    n_sig_005 = sum(P_value < 0.05, na.rm = TRUE),
    prop_sig_005 = n_sig_005 / n_tests,
    n_sig_FDR_005 = sum(Significant_FDR_0.05, na.rm = TRUE),
    prop_sig_FDR_005 = n_sig_FDR_005 / n_tests,
    .groups = "drop"
  ) %>%
  arrange(desc(prop_sig_FDR_005), desc(prop_sig_005))

cat("\nHWE summary by site:\n")
print(hwe_summary)

### 13) CLONALITY SUMMARY BY SITE ----
clonality_summary <- geno %>%
  group_by(Site) %>%
  summarise(
    N = n(),
    nMLG = n_distinct(MLG),
    nMLL = n_distinct(MLL),
    clonal_fraction_MLG = 1 - (nMLG / N),
    clonal_fraction_MLL = 1 - (nMLL / N),
    n_shared_MLG = sum(table(MLG) >= 2),
    n_shared_MLL = sum(table(MLL) >= 2),
    Lat_mean = if (has_lat) mean(Latitude, na.rm = TRUE) else NA_real_,
    Lon_mean = if (has_lon) mean(Longitude, na.rm = TRUE) else NA_real_,
    mean_DBH = if (has_dbh) mean(Dhp_tige, na.rm = TRUE) else NA_real_,
    .groups = "drop"
  ) %>%
  arrange(Site)

cat("\nClonality summary by Site:\n")
print(clonality_summary)

### 13B) RARE-CLONALITY ROBUST ANALYSES (pairs + distances + permutation) ----
clone_pairs <- NULL
clone_distance_test <- NULL
clonality_latitude_test <- NULL

# Helper: all unique pairs within a vector
make_pairs <- function(x) {
  x <- unique(x)
  if (length(x) < 2) return(NULL)
  combn(x, 2, simplify = FALSE)
}

# Use MLL as the "clone" definition (more conservative than exact MLG if you want)
if (has_lat && has_lon && ok_geo) {
  library(geosphere)
  
  # Identify clone groups (MLL) with >=2 individuals
  mll_sizes <- table(geno$MLL)
  shared_mll <- as.integer(names(mll_sizes[mll_sizes >= 2]))
  
  if (length(shared_mll) == 0) {
    cat("\n[INFO] No shared MLL detected -> skipping clone-pair distance analyses.\n")
  } else {
    # Build pair table
    pair_rows <- list()
    
    for (g in shared_mll) {
      members <- geno %>% filter(MLL == g)
      
      # Option: keep within-site only (prevents weird long-distance “same clone” if any ID issues)
      members_by_site <- split(members, members$Site)
      
      for (s in names(members_by_site)) {
        msub <- members_by_site[[s]]
        pr <- make_pairs(msub$Nom_Labo_Échantillons)
        if (is.null(pr)) next
        
        for (k in seq_along(pr)) {
          id1 <- pr[[k]][1]
          id2 <- pr[[k]][2]
          
          a <- msub %>% filter(Nom_Labo_Échantillons == id1) %>% slice(1)
          b <- msub %>% filter(Nom_Labo_Échantillons == id2) %>% slice(1)
          
          d_m <- geosphere::distHaversine(
            c(a$Longitude, a$Latitude),
            c(b$Longitude, b$Latitude)
          )
          
          pair_rows[[length(pair_rows) + 1]] <- data.frame(
            clone_def = "MLL",
            clone_id = g,
            Site = s,
            ind1 = id1,
            ind2 = id2,
            dist_m = as.numeric(d_m),
            stringsAsFactors = FALSE
          )
        }
      }
    }
    
    clone_pairs <- bind_rows(pair_rows)
    
    if (nrow(clone_pairs) == 0) {
      cat("\n[INFO] Shared MLL exist, but no within-site clone pairs -> skipping distance test.\n")
    } else {
      cat("\nClone pair distances (first rows):\n")
      print(head(clone_pairs, 20))
      
      # Permutation test: are clone pairs closer than random pairs within same site?
      # Statistic: median distance of clone pairs (robust in small N)
      obs <- median(clone_pairs$dist_m, na.rm = TRUE)
      
      # Build a function to draw random pairs matching number of clone pairs per site
      set.seed(1)
      perm_stats <- numeric(n_perm_clone_distance)
      
      pairs_per_site <- clone_pairs %>% count(Site, name = "n_pairs")
      
      for (p in seq_len(n_perm_clone_distance)) {
        rand_d <- c()
        
        for (i in seq_len(nrow(pairs_per_site))) {
          s <- pairs_per_site$Site[i]
          np <- pairs_per_site$n_pairs[i]
          
          sub <- geno %>% filter(Site == s) %>% filter(!is.na(Latitude), !is.na(Longitude))
          if (nrow(sub) < 2) next
          
          # Sample 2*np individuals with replacement if needed, then pair them up
          # (If you want strictly without replacement, we can tighten it, but this is stable for small n.)
          ids <- sample(sub$Nom_Labo_Échantillons, size = 2 * np, replace = TRUE)
          
          for (k in seq(1, length(ids), by = 2)) {
            a <- sub %>% filter(Nom_Labo_Échantillons == ids[k]) %>% slice(1)
            b <- sub %>% filter(Nom_Labo_Échantillons == ids[k + 1]) %>% slice(1)
            
            d_m <- geosphere::distHaversine(
              c(a$Longitude, a$Latitude),
              c(b$Longitude, b$Latitude)
            )
            rand_d <- c(rand_d, as.numeric(d_m))
          }
        }
        
        perm_stats[p] <- median(rand_d, na.rm = TRUE)
      }
      
      # one-sided p-value: clone pairs are SMALLER than random
      pval <- (sum(perm_stats <= obs) + 1) / (length(perm_stats) + 1)
      
      clone_distance_test <- data.frame(
        statistic = "median_distance_m",
        observed = obs,
        perm_mean = mean(perm_stats, na.rm = TRUE),
        perm_median = median(perm_stats, na.rm = TRUE),
        p_one_sided = pval,
        n_clone_pairs = nrow(clone_pairs),
        n_perm = n_perm_clone_distance
      )
      
      cat("\nPermutation test (clone distances vs random within-site):\n")
      print(clone_distance_test)
    }
  }
} else {
  cat("\n[INFO] Skipping clone-pair distance analyses (need Latitude+Longitude + geosphere).\n")
}

# Exploratory: clonality ~ latitude at site level (Spearman)
if (has_lat) {
  tmp <- clonality_summary %>%
    filter(!is.na(Lat_mean), !is.na(clonal_fraction_MLL))
  
  if (nrow(tmp) >= min_sites_for_latitude_test) {
    ct <- suppressWarnings(cor.test(tmp$Lat_mean, tmp$clonal_fraction_MLL, method = "spearman", exact = FALSE))
    clonality_latitude_test <- data.frame(
      method = "spearman",
      n_sites = nrow(tmp),
      rho = as.numeric(ct$estimate),
      p_value = as.numeric(ct$p.value)
    )
    cat("\nExploratory clonality~latitude (Spearman across sites):\n")
    print(clonality_latitude_test)
  } else {
    cat("\n[INFO] Not enough sites with latitude for clonality~latitude test.\n")
  }
}

### 14) OPTIONAL: SPATIAL GENETIC STRUCTURE (MANTEL) ----
mantel_results <- NULL

if (has_lat && has_lon && ok_vegan && ok_geo) {
  library(vegan)
  library(geosphere)
  
  Dgen_all <- bruvo.dist(gen, replen = replen)
  Dgen_all <- as.matrix(Dgen_all)
  
  mantel_one_site <- function(site_code) {
    sub <- geno %>% filter(Site == site_code)
    if (nrow(sub) < 4) {
      return(data.frame(Site = site_code, N = nrow(sub), mantel_r = NA_real_, mantel_p = NA_real_))
    }
    
    inds <- sub$Nom_Labo_Échantillons
    idx  <- match(inds, rownames(Dgen_all))
    gen_mat <- Dgen_all[idx, idx]
    
    geo_mat <- geosphere::distm(
      sub[, c("Longitude", "Latitude")],
      fun = geosphere::distHaversine
    )
    
    m <- vegan::mantel(as.dist(gen_mat), as.dist(geo_mat),
                       method = "spearman", permutations = 999)
    
    data.frame(
      Site = site_code,
      N = nrow(sub),
      mantel_r = as.numeric(m$statistic),
      mantel_p = as.numeric(m$signif)
    )
  }
  
  mantel_results <- bind_rows(lapply(sites, mantel_one_site)) %>%
    arrange(Site)
  
  cat("\nMantel results (genetic vs geographic distance):\n")
  print(mantel_results)
  
} else {
  cat("\n[INFO] Skipping Mantel spatial test (need Latitude+Longitude + packages vegan/geosphere).\n")
}

### 14B) NJ TREE (Bruvo) + “boxes” for MLL clusters ----
cat("\nBuilding NJ tree (Bruvo distance)...\n")
D_bruvo <- bruvo.dist(gen, replen = replen)
tree_nj <- ape::nj(as.dist(D_bruvo))

# Basic PDF (always)
pdf(tree_out_basic_pdf, width = 10, height = 12)
plot(tree_nj, cex = tree_tip_cex, no.margin = TRUE)
title(main = "Neighbor-Joining tree (Bruvo distance)")
dev.off()

cat("Saved basic NJ tree to:", file.path(getwd(), tree_out_basic_pdf), "\n")

# Boxed/highlighted version (requires ggtree)
if (ok_ggtree && ok_ggplot) {
  library(ggtree)
  library(ggplot2)
  
  mll_vec <- geno$MLL
  names(mll_vec) <- geno$Nom_Labo_Échantillons
  
  tip_site <- geno$Site
  names(tip_site) <- geno$Nom_Labo_Échantillons
  
  tip_df <- data.frame(label = tree_nj$tip.label,
                       Site  = tip_site[tree_nj$tip.label],
                       MLL   = mll_vec[tree_nj$tip.label],
                       stringsAsFactors = FALSE)
  
  p <- ggtree(tree_nj) %<+% tip_df +
    geom_tiplab(aes(label = label), size = 2) +
    ggtitle("Neighbor-Joining tree (Bruvo distance) with MLL clusters highlighted")
  
  tab_mll <- table(mll_vec)
  big_groups <- names(tab_mll[tab_mll >= min_group_for_box])
  
  for (g in big_groups) {
    tips <- names(mll_vec[mll_vec == as.integer(g)])
    tips <- intersect(tips, tree_nj$tip.label)
    if (length(tips) >= min_group_for_box) {
      node <- tryCatch(ggtree::MRCA(tree_nj, tips), error = function(e) NA_integer_)
      if (!is.na(node)) {
        p <- p + geom_hilight(node = node, alpha = 0.25)
      }
    }
  }
  
  ggsave(tree_out_boxed_pdf, plot = p, width = 10, height = 12)
  cat("Saved boxed NJ tree to:", file.path(getwd(), tree_out_boxed_pdf), "\n")
  
} else {
  cat("\n[INFO] ggtree not installed, so only basic NJ PDF was created.\n")
  cat("To install: BiocManager::install('ggtree')\n")
}

### 14C) NJ TREE BY SITE (population-level NJ) ----
cat("\nBuilding NJ tree BY SITE (population-level)...\n")

# Set SITE as pop for genind, then collapse -> genpop
gen_site <- gen
pop(gen_site) <- as.factor(geno$Site)

# Convert individuals -> site allele frequencies
gp_site <- genind2genpop(gen_site)

# A) Nei distance among sites
D_site_nei <- poppr::nei.dist(gp_site)
tree_site_nei <- ape::nj(D_site_nei)

# Save + plot
write.csv(as.matrix(D_site_nei), tree_site_nei_csv)
write.tree(tree_site_nei, file = tree_site_nei_nwk)

pdf(tree_site_nei_pdf, width = 10, height = 6)
plot(tree_site_nei, cex = 0.9, no.margin = TRUE)
title(main = "Neighbor-Joining tree BY SITE (Nei distance)")
dev.off()

cat("Saved NJ by site (Nei) to:\n",
    "  ", file.path(getwd(), tree_site_nei_pdf), "\n",
    "  ", file.path(getwd(), tree_site_nei_nwk), "\n",
    "  ", file.path(getwd(), tree_site_nei_csv), "\n", sep = "")

# B) Optional chord distance (often good for microsats) if supported by your adegenet version
D_site_chord <- tryCatch(
  adegenet::dist.genpop(gp_site, method = "CSE"),
  error = function(e) NULL
)

if (!is.null(D_site_chord)) {
  tree_site_chord <- ape::nj(D_site_chord)
  
  write.csv(as.matrix(D_site_chord), tree_site_chord_csv)
  write.tree(tree_site_chord, file = tree_site_chord_nwk)
  
  pdf(tree_site_chord_pdf, width = 10, height = 6)
  plot(tree_site_chord, cex = 0.9, no.margin = TRUE)
  title(main = "Neighbor-Joining tree BY SITE (Chord/CSE distance)")
  dev.off()
  
  cat("Saved NJ by site (Chord/CSE) to:\n",
      "  ", file.path(getwd(), tree_site_chord_pdf), "\n",
      "  ", file.path(getwd(), tree_site_chord_nwk), "\n",
      "  ", file.path(getwd(), tree_site_chord_csv), "\n", sep = "")
} else {
  cat("[INFO] Chord/CSE distance not available in your adegenet version; kept Nei NJ by site only.\n")
}

### 15) SAVE OUTPUTS ----
if (save_outputs) {
  write.csv(geno, "genotypes_avec_MLG_MLL_QC_SITE.csv", row.names = FALSE)
  write.csv(clonality_summary, "clonality_summary_by_site.csv", row.names = FALSE)
  write.csv(hwe_results, "HWE_per_site_exact.csv", row.names = FALSE)
  write.csv(hwe_summary, "HWE_summary_by_site.csv", row.names = FALSE)
  
  if (!is.null(mantel_results)) {
    write.csv(mantel_results, "mantel_genetic_vs_geo_by_site.csv", row.names = FALSE)
  }
  
  # NEW outputs
  if (!is.null(clone_pairs)) {
    write.csv(clone_pairs, "clone_pairs_MLL_within_site_geodist.csv", row.names = FALSE)
  }
  if (!is.null(clone_distance_test)) {
    write.csv(clone_distance_test, "clone_distance_permutation_test.csv", row.names = FALSE)
  }
  if (!is.null(clonality_latitude_test)) {
    write.csv(clonality_latitude_test, "clonality_latitude_spearman.csv", row.names = FALSE)
  }
  
  cat("\nSaved outputs to:\n")
  cat("  ", getwd(), "\n", sep = "")
}

cat("\nDONE ✅\n")
