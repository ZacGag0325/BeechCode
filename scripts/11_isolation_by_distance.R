############################################################
# scripts/11_isolation_by_distance.R
# Isolation by Distance (population-level Mantel: Nei vs geographic)
############################################################

options(repos = c(CRAN = "https://cloud.r-project.org"))
pkgs <- c("readxl", "dplyr", "tibble", "tidyr", "stringr", "adegenet", "poppr", "vegan", "geosphere", "ggplot2")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(stringr)
  library(adegenet)
  library(poppr)
  library(vegan)
  library(geosphere)
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

normalize_site <- function(x) {
  x <- str_to_upper(str_squish(as.character(x)))
  x[nchar(x) == 0] <- NA_character_
  x
}

read_required_sheet <- function(path, sheet, required_cols) {
  out <- readxl::read_xlsx(path, sheet = sheet)
  miss <- setdiff(required_cols, names(out))
  if (length(miss) > 0) {
    stop("Sheet '", sheet, "' in ", basename(path), " is missing columns: ", paste(miss, collapse = ", "))
  }
  out
}

find_existing_nei_matrix <- function(path) {
  sh <- readxl::excel_sheets(path)
  nei_sheet <- sh[str_detect(str_to_lower(sh), "nei")]
  if (length(nei_sheet) > 0) {
    return(list(found = TRUE, reason = paste0("Sheet name suggests Nei matrix: ", paste(nei_sheet, collapse = ", "))))
  }
  list(found = FALSE, reason = "No sheet name indicates a precomputed Nei matrix in the input workbooks.")
}

PROJECT_ROOT <- find_project_root()
setwd(PROJECT_ROOT)
source(file.path("scripts", "_load_objects.R"))

RUN_OUT <- if (exists("RUN_OUT", inherits = TRUE)) get("RUN_OUT", inherits = TRUE) else file.path(PROJECT_ROOT, "outputs", "v1")
OUTDIR <- file.path(RUN_OUT, "ibd_only")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# Input paths requested by user
field_xlsx <- "/mnt/data/donnees_modifiees_west_summer2024 copie.xlsx"
geno_xlsx <- "/mnt/data/poppr.xlsx"

if (!file.exists(field_xlsx)) stop("Field workbook not found: ", field_xlsx)
if (!file.exists(geno_xlsx)) stop("Genotype workbook not found: ", geno_xlsx)

# -----------------------------------------------------------------------------
# Variable dictionary (French -> English)
# -----------------------------------------------------------------------------
# point_prisme: prism point orientation (N, SE, or SO)
# valeur_prisme: prism inclusion class (2 = tree inside plot; 1 = tree on boundary)
# espece: 3-letter botanical species code
# DHP_cm: DBH (diameter at breast height) in cm
# etat: tree status (V = alive, M = dead)
# MCH: beech bark disease presence/absence for Fagus grandifolia (P/A; NA for others)
# p_rec_MCH: percent stem covered by MCH (class)
# dep_MCH: decline class (1/2/3)
# mo: organic layer thickness
# surface: charcoal presence/absence (P/A)
# r_utilise: radius used around prism point (e.g., 3.57 m)
# compte: regeneration/sapling count per species in subplot
# id_spa: regeneration subplot id (starts North, then clockwise)
# id_tige: sampled stem id (1-24 per site)
# dhp_tige: DBH for each sampled stem
# point_gps: GPS coordinates for each sampled stem
# strate: canopy stratum (inf/moy/sup)

# -----------------------------------------------------------------------------
# 1) Read and standardize required sheets/columns
# -----------------------------------------------------------------------------
field_genetique <- read_required_sheet(
  field_xlsx,
  "genetique",
  c("Site", "Id_tige", "Dhp_tige", "Name_terrain", "Nom_Labo_Échantillons", "Geometrytype", "Latitude", "Longitude", "Elevation", "Geoidoffset")
) %>%
  transmute(
    Site = as.character(Site),
    Nom_Labo_Échantillons = as.character(Nom_Labo_Échantillons),
    Latitude = suppressWarnings(as.numeric(Latitude)),
    Longitude = suppressWarnings(as.numeric(Longitude))
  )

# Read other sheets (structure check + cleaning readiness)
field_arbre <- read_required_sheet(field_xlsx, "arbre", c("Site", "Point_prisme", "Valeur_prisme", "Espece", "Dhp_cm", "Etat", "Mch", "P_rec_mch", "Dep_mch", "St_tige", "Nb_tiges_ha"))
field_gaule <- read_required_sheet(field_xlsx, "gaule", c("Site", "Point_prisme", "R_utilise", "Espece", "Dhp_cm", "Etat", "Mch", "P_rec_mch", "Dep_mch"))
field_regen <- read_required_sheet(field_xlsx, "regeneration", c("Site", "Id_spa", "Espece", "Compte"))
field_veg <- read_required_sheet(field_xlsx, "vegetation", c("Site", "Espece", "Nombre"))
field_sol <- read_required_sheet(field_xlsx, "sol", c("Site", "Id_spa", "Mo", "Surface"))
field_dict <- read_required_sheet(field_xlsx, "dictionnaire", names(readxl::read_xlsx(field_xlsx, sheet = "dictionnaire", n_max = 0)))

geno_ind <- read_required_sheet(geno_xlsx, "Genotypes-Indiv.", c("Nom_Labo_Échantillons", "Numéro_Population"))
legend_sites <- read_required_sheet(geno_xlsx, "Légende-Sites", c("Numéro_Population", "Nom_Terrain"))

# -----------------------------------------------------------------------------
# 2) Build one coordinate per Site (centroids from genetique sheet)
# -----------------------------------------------------------------------------
site_centroids <- field_genetique %>%
  mutate(Site = normalize_site(Site)) %>%
  filter(!is.na(Site)) %>%
  group_by(Site) %>%
  summarise(
    Latitude = mean(Latitude, na.rm = TRUE),
    Longitude = mean(Longitude, na.rm = TRUE),
    n_stems_with_gps = sum(!is.na(Latitude) & !is.na(Longitude)),
    .groups = "drop"
  ) %>%
  mutate(
    Latitude = ifelse(is.nan(Latitude), NA_real_, Latitude),
    Longitude = ifelse(is.nan(Longitude), NA_real_, Longitude)
  )

# -----------------------------------------------------------------------------
# 3) Build genind from poppr-style Excel and compute Nei distance by population
# -----------------------------------------------------------------------------
legend_clean <- legend_sites %>%
  transmute(
    Numéro_Population = as.character(Numéro_Population),
    Site = normalize_site(Nom_Terrain)
  ) %>%
  filter(!is.na(Numéro_Population), !is.na(Site)) %>%
  distinct(Numéro_Population, .keep_all = TRUE)

geno_clean <- geno_ind %>%
  mutate(
    Nom_Labo_Échantillons = as.character(Nom_Labo_Échantillons),
    Numéro_Population = as.character(`Numéro_Population`)
  ) %>%
  left_join(legend_clean, by = "Numéro_Population")

locus_cols <- names(geno_clean)[str_detect(names(geno_clean), "_[12]$")]
if (length(locus_cols) == 0) {
  stop("No locus allele columns found (expected names like sfc_0036_1, sfc_0036_2).")
}

loci_base <- unique(str_remove(locus_cols, "_[12]$"))
allele_df <- data.frame(row.names = geno_clean$Nom_Labo_Échantillons, stringsAsFactors = FALSE)
for (loc in loci_base) {
  a1 <- paste0(loc, "_1")
  a2 <- paste0(loc, "_2")
  if (!all(c(a1, a2) %in% names(geno_clean))) next
  allele_df[[loc]] <- paste0(as.character(geno_clean[[a1]]), "/", as.character(geno_clean[[a2]]))
}

if (ncol(allele_df) == 0) stop("No paired loci (_1/_2) available to build genind object.")

if (anyDuplicated(rownames(allele_df)) > 0) {
  keep <- !duplicated(rownames(allele_df))
  allele_df <- allele_df[keep, , drop = FALSE]
  geno_clean <- geno_clean[keep, , drop = FALSE]
}

gi_from_excel <- adegenet::df2genind(
  X = allele_df,
  sep = "/",
  ploidy = 2,
  ncode = 3,
  ind.names = rownames(allele_df),
  pop = as.factor(geno_clean$Site),
  NA.char = c("NA", "NA/NA", "<NA>/<NA>", "-/-", "0/0")
)

# Drop individuals with missing population label after legend join
missing_pop_inds <- which(is.na(pop(gi_from_excel)) | pop(gi_from_excel) == "NA")
if (length(missing_pop_inds) > 0) gi_from_excel <- gi_from_excel[-missing_pop_inds]

if (nInd(gi_from_excel) < 2) stop("Not enough individuals with mapped population labels.")

gp <- adegenet::genind2genpop(gi_from_excel)
nei_obj <- tryCatch(poppr::nei.dist(gp), error = function(e) e)
if (inherits(nei_obj, "error")) {
  nei_obj <- tryCatch(adegenet::dist.genpop(gp, method = 1), error = function(e) e)
}
if (inherits(nei_obj, "error")) stop("Could not compute Nei distance matrix: ", nei_obj$message)

gen_mat <- as.matrix(nei_obj)
rownames(gen_mat) <- normalize_site(rownames(gen_mat))
colnames(gen_mat) <- normalize_site(colnames(gen_mat))
gen_mat <- (gen_mat + t(gen_mat)) / 2
diag(gen_mat) <- 0

# Optional Bruvo at individual level (not used for population Mantel)
bruvo_available <- FALSE
if (inherits(gi_from_excel, "genind") && nInd(gi_from_excel) > 1) {
  bruvo_available <- TRUE
}

# -----------------------------------------------------------------------------
# 4) Geographic matrix from centroid coordinates (great-circle, km)
# -----------------------------------------------------------------------------
coord_ok <- site_centroids %>% filter(!is.na(Latitude), !is.na(Longitude))

# -----------------------------------------------------------------------------
# 5) Align both matrices on the exact same site set/order + report dropped sites
# -----------------------------------------------------------------------------
sites_in_gen <- rownames(gen_mat)
sites_in_geo <- coord_ok$Site

all_sites <- sort(unique(c(sites_in_gen, sites_in_geo)))
drop_tbl <- tibble(Site = all_sites) %>%
  mutate(
    has_genetic = Site %in% sites_in_gen,
    has_coordinates = Site %in% sites_in_geo,
    drop_reason = case_when(
      has_genetic & has_coordinates ~ NA_character_,
      !has_genetic & has_coordinates ~ "Missing genotypes/population mapping",
      has_genetic & !has_coordinates ~ "Missing or invalid coordinates",
      TRUE ~ "Missing both genetics and coordinates"
    )
  )

sites_keep <- drop_tbl %>% filter(is.na(drop_reason)) %>% pull(Site)

if (length(sites_keep) < 3) {
  stop("Need at least 3 sites after alignment. Kept: ", length(sites_keep))
}

coord_keep <- coord_ok %>% filter(Site %in% sites_keep) %>% arrange(match(Site, sites_keep))
geo_mat <- geosphere::distm(as.matrix(coord_keep[, c("Longitude", "Latitude")]), fun = geosphere::distHaversine) / 1000
rownames(geo_mat) <- coord_keep$Site
colnames(geo_mat) <- coord_keep$Site
geo_mat <- (geo_mat + t(geo_mat)) / 2
diag(geo_mat) <- 0

gen_mat <- gen_mat[sites_keep, sites_keep, drop = FALSE]
geo_mat <- geo_mat[sites_keep, sites_keep, drop = FALSE]

if (!identical(rownames(gen_mat), rownames(geo_mat))) stop("Matrix row names are not aligned after filtering.")
if (!identical(colnames(gen_mat), colnames(geo_mat))) stop("Matrix column names are not aligned after filtering.")
if (any(is.na(gen_mat[upper.tri(gen_mat)]))) stop("Genetic matrix contains NA off-diagonal values.")
if (any(is.na(geo_mat[upper.tri(geo_mat)]))) stop("Geographic matrix contains NA off-diagonal values.")

# -----------------------------------------------------------------------------
# 6) Mantel test + outputs
# -----------------------------------------------------------------------------
mantel_out <- vegan::mantel(
  xdis = as.dist(gen_mat),
  ydis = as.dist(geo_mat),
  method = "pearson",
  permutations = 9999
)

pair_idx <- which(upper.tri(gen_mat), arr.ind = TRUE)
scatter_df <- tibble(
  Site1 = rownames(gen_mat)[pair_idx[, 1]],
  Site2 = colnames(gen_mat)[pair_idx[, 2]],
  Geographic_km = geo_mat[pair_idx],
  Nei_distance = gen_mat[pair_idx]
)

p <- ggplot(scatter_df, aes(x = Geographic_km, y = Nei_distance)) +
  geom_point(size = 2, alpha = 0.9, color = "#1f78b4") +
  geom_smooth(method = "lm", se = TRUE, color = "#d7301f", fill = "#fcae91", linewidth = 0.9) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Mantel test at population level",
    subtitle = paste0("r = ", signif(as.numeric(mantel_out$statistic), 4), ", p = ", signif(as.numeric(mantel_out$signif), 4)),
    x = "Geographic distance (km)",
    y = "Nei genetic distance"
  )

ggsave(file.path(OUTDIR, "ibd_scatter_nei_vs_geo.png"), p, width = 8, height = 6, dpi = 300)

mantel_tbl <- tibble(
  analysis = "mantel_nei_vs_geographic_km",
  statistic_r = as.numeric(mantel_out$statistic),
  p_value = as.numeric(mantel_out$signif),
  permutations = as.integer(mantel_out$permutations),
  n_sites = length(sites_keep),
  n_pairs = nrow(scatter_df)
)

nei_check_field <- find_existing_nei_matrix(field_xlsx)
nei_check_geno <- find_existing_nei_matrix(geno_xlsx)
nei_status <- tibble(
  input_file = c(basename(field_xlsx), basename(geno_xlsx)),
  nei_matrix_found = c(nei_check_field$found, nei_check_geno$found),
  note = c(nei_check_field$reason, nei_check_geno$reason)
)

write.csv(site_centroids, file.path(OUTDIR, "site_centroids_from_genetique.csv"), row.names = FALSE)
write.csv(geo_mat, file.path(OUTDIR, "matrix_geographic_distance_km.csv"), row.names = TRUE)
write.csv(gen_mat, file.path(OUTDIR, "matrix_genetic_distance_nei.csv"), row.names = TRUE)
write.csv(scatter_df, file.path(OUTDIR, "ibd_pairwise_table.csv"), row.names = FALSE)
write.csv(mantel_tbl, file.path(OUTDIR, "mantel_results.csv"), row.names = FALSE)
write.csv(drop_tbl %>% filter(!is.na(drop_reason)), file.path(OUTDIR, "dropped_sites_summary.csv"), row.names = FALSE)
write.csv(nei_status, file.path(OUTDIR, "nei_matrix_input_check.csv"), row.names = FALSE)

summary_lines <- c(
  "Population-level Mantel workflow completed.",
  paste0("Sites kept (both matrices): ", length(sites_keep)),
  paste0("Sites dropped: ", sum(!is.na(drop_tbl$drop_reason))),
  if (sum(!is.na(drop_tbl$drop_reason)) > 0) {
    paste0("Dropped details: ", paste0(drop_tbl$Site[!is.na(drop_tbl$drop_reason)], " [", drop_tbl$drop_reason[!is.na(drop_tbl$drop_reason)], "]", collapse = "; "))
  } else {
    "Dropped details: none"
  },
  paste0("Mantel r = ", signif(as.numeric(mantel_out$statistic), 5), "; p = ", signif(as.numeric(mantel_out$signif), 5)),
  paste0("Nei matrix pre-existing in inputs? ", ifelse(any(nei_status$nei_matrix_found), "Possibly yes (see nei_matrix_input_check.csv)", "No, computed from genotype loci.")),
  paste0("Bruvo (individual-level) option available in script: ", bruvo_available, " (not used for this Mantel).")
)
writeLines(summary_lines, con = file.path(OUTDIR, "mantel_summary.txt"))

cat(paste(summary_lines, collapse = "\n"), "\n")
cat("Outputs written to: ", OUTDIR, "\n", sep = "")
