############################################################
# POPPR / ADEGENET PIPELINE (MLG + MLL via Bruvo) – VERSION FIX
# + Liste des clones MLG/MLL avec individus + populations
############################################################

### 0. PACKAGES ----
library(adegenet)
library(poppr)
library(readxl)
library(dplyr)
library(tidyr)
library(hierfstat)

### 0B. OPTION : SAUVEGARDES SUR DISQUE ? ----
save_outputs <- FALSE   # mets TRUE si tu veux save les table.csv

### 1. DOSSIER DE TRAVAIL ----
setwd("/Users/zacharygagnon/Desktop/R data workshop")
getwd()
list.files()

### 2. LECTURE DU FICHIER D’ALLÈLES ----
geno_raw <- read_xlsx("poppr_avec_dup_E1&2.xlsx")

# 2A) Nettoyage: enlever les colonnes complètement vides + colonnes "Unnamed"
geno_raw <- geno_raw %>%
  select(where(~ !all(is.na(.)))) %>%      # drop colonnes 100% NA
  select(-matches("^Unnamed"))             # drop colonnes Unnamed: xx (Excel)

# (optionnel) enlever espaces dans les noms de colonnes
names(geno_raw) <- trimws(names(geno_raw))

# 2B) ID en caractère
geno_raw$Nom_Labo_Échantillons <- as.character(geno_raw$Nom_Labo_Échantillons)

# On garde toutes les populations (incl. ML1 si elle existe)
geno <- geno_raw

# Vérifier les codes de population
print(table(geno$Numéro_Population))

### 3. DÉFINITION DES COLONNES D’ALLÈLES (ROBUSTE) ----
# On prend uniquement les colonnes dont le nom finit par _1 ou _2
allele_cols <- which(grepl("_(1|2)$", names(geno)))

if (length(allele_cols) == 0) {
  stop("Aucune colonne d'allèles détectée. Assure-toi que tes colonnes finissent par _1 et _2.")
}

# Sécurité : vérifier qu’on a bien un nombre pair de colonnes d’allèles
if (length(allele_cols) %% 2 != 0) {
  stop("Le nombre de colonnes d'allèles détectées (_1/_2) n'est pas pair. Vérifie les noms de colonnes.")
}

# Vérifier que chaque locus a exactement 2 colonnes (une _1 et une _2)
locus_base <- sub("_(1|2)$", "", names(geno)[allele_cols])
bad <- table(locus_base)[table(locus_base) != 2]
if (length(bad) > 0) {
  print(bad)
  stop("Au moins un locus n'a pas exactement 2 colonnes (une _1 et une _2). Corrige les entêtes dans Excel.")
}

# Indices des premières colonnes de chaque locus (= celles qui finissent par _1)
locus_starts <- allele_cols[grepl("_1$", names(geno)[allele_cols])]
locus_starts <- sort(locus_starts)

### 4. CONSTRUIRE UN TABLEAU "1 COLONNE PAR LOCUS" ----
# Chaque colonne = génotype "allèle1/allèle2" (texte)
geno_locus <- data.frame(row.names = geno$Nom_Labo_Échantillons)

for (i in seq_along(locus_starts)) {
  col1_idx <- locus_starts[i]
  col2_idx <- col1_idx + 1
  
  col1_name <- names(geno)[col1_idx]
  col2_name <- names(geno)[col2_idx]
  
  # Nom du locus = base sans _1/_2
  locus_name <- sub("_1$", "", col1_name)
  
  # Construire le génotype allélique "A1/A2"
  a1 <- as.character(geno[[col1_name]])
  a2 <- as.character(geno[[col2_name]])
  
  geno_locus[[locus_name]] <- paste(a1, a2, sep = "/")
}

# Remplacer les génotypes contenant "NA" par NA complet
geno_locus[] <- lapply(geno_locus, function(col) {
  col[grepl("NA", col, fixed = TRUE)] <- NA
  col
})

print(head(geno_locus))
print(str(geno_locus))

### 4B. NETTOYER LES NOMS DES COLONNES ----
clean_names <- names(geno_locus)
clean_names <- gsub("\\.+", "_", clean_names)
names(geno_locus) <- clean_names
print(names(geno_locus))

### 5. CRÉER L’OBJET genind ----
gen <- df2genind(
  geno_locus,
  ploidy    = 2,
  ind.names = rownames(geno_locus),
  pop       = as.factor(geno$Numéro_Population),
  type      = "codom",
  sep       = "/",
  ncode     = 3
)

print(gen)
print(summary(gen))

### (CHECK) : les allèles doivent être numériques pour Bruvo ----
aa <- alleles(gen)
bad_alleles <- lapply(aa, function(v) v[is.na(as.numeric(v))])
bad_alleles <- bad_alleles[sapply(bad_alleles, length) > 0]

if (length(bad_alleles) > 0) {
  print(bad_alleles)
  stop("Certains allèles ne sont pas numériques. Corrige les valeurs non-numériques dans le fichier.")
} else {
  message("OK: tous les allèles sont numériques (Bruvo compatible).")
}

### 6. PASSER EN genclone ET CALCULER LES MLG (naïfs) ----
gc_mlg <- as.genclone(gen)
print(gc_mlg)

mlg_vec <- mlg.vector(gc_mlg)
print(head(mlg_vec))
print(length(unique(mlg_vec)))

mlg_tab <- mlg.table(gc_mlg)
print(mlg_tab)

# Liste des individus par clone (MLG)
mlg_list <- mlg.id(gc_mlg)

### 7. AJOUTER LES MLG AU TABLEAU GÉNOTYPIQUE ----
geno$MLG <- as.integer(mlg_vec)
print(head(geno[, c("Nom_Labo_Échantillons", "Numéro_Population", "MLG")]))

if (save_outputs) {
  write.csv(geno, "genotypes_avec_MLG.csv", row.names = FALSE)
}

### 8. STATISTIQUES CLONALES PAR POPULATION (MLG NAÏFS) ----
clone_stats_mlg <- poppr(gc_mlg)
print(clone_stats_mlg)

if (save_outputs) {
  write.csv(clone_stats_mlg, "clone_stats_par_population_MLG.csv", row.names = TRUE)
}

### 9. RICHESSE ALLÉLIQUE RAREFIÉE (N = 14) ----
hf <- genind2hierfstat(gen)
ar <- hierfstat::allelic.richness(hf, min.n = 14)
print(ar)

if (save_outputs) {
  write.csv(ar$Ar, "allelic_richness_rarefied_min14.csv", row.names = TRUE)
}

###############################
#  10. MLL AVEC BRUVO.DIST   #
###############################

### 10A. Longueur de répétition par locus ----
replen <- rep(2, nLoc(gen))
names(replen) <- locNames(gen)
print(replen)

### 10B. Explorer les seuils possibles avec filter_stats ----
fs <- filter_stats(
  gc_mlg,
  distance  = bruvo.dist,
  replen    = replen,
  plot      = TRUE
)

farthest_thresh <- cutoff_predictor(fs$farthest$THRESHOLDS)
print(farthest_thresh)

### 10C. ESSAI D'UNE GRILLE DE SEUILS BRUVO ----
thresh_seq <- c(0, 0.005, 0.01, 0.02, 0.03, 0.05, 0.1)

mll_summary <- data.frame(
  threshold       = thresh_seq,
  nmll            = NA_integer_,
  mean_clone_size = NA_real_
)

for (i in seq_along(thresh_seq)) {
  th <- thresh_seq[i]
  
  gc_tmp <- gc_mlg
  mlg.filter(
    gc_tmp,
    distance  = bruvo.dist,
    replen    = replen,
    algorithm = "farthest_neighbor"
  ) <- th
  
  mll_vec_tmp <- mll(gc_tmp)
  mll_tab_tmp <- table(mll_vec_tmp)
  
  mll_summary$nmll[i]            <- length(mll_tab_tmp)
  mll_summary$mean_clone_size[i] <- mean(mll_tab_tmp)
}

print(mll_summary)

if (save_outputs) {
  write.csv(mll_summary, "mll_summary_thresholds.csv", row.names = FALSE)
}

### 10D. CHOISIR UN SEUIL POUR DÉFINIR LES MLL ----
chosen_thresh <- 0.09

gc_mll <- gc_mlg
mlg.filter(
  gc_mll,
  distance  = bruvo.dist,
  replen    = replen,
  algorithm = "farthest_neighbor"
) <- chosen_thresh

print(gc_mll)

print(nmll(gc_mlg))
print(nmll(gc_mll))

mll_tab <- mlg.table(gc_mll)
print(mll_tab)

mll_vec <- mll(gc_mll)
print(head(mll_vec))

geno$MLL <- as.integer(mll_vec)
print(head(geno[, c("Nom_Labo_Échantillons", "Numéro_Population", "MLG", "MLL")]))

if (save_outputs) {
  write.csv(geno, "genotypes_avec_MLG_MLL.csv", row.names = FALSE)
}

### 10E. Statistiques clonalité par population avec MLL ----
clone_stats_mll <- poppr(gc_mll)
print(clone_stats_mll)

if (save_outputs) {
  write.csv(clone_stats_mll, "clone_stats_par_population_MLL.csv", row.names = TRUE)
}

############################################################
# 11. LISTES CLONALES (INDIVIDUS + POPULATION) – MLG & MLL
############################################################

### 11A) MLG: liste style "mlg.id" mais avec Population aussi
mlg_full_list <- split(
  data.frame(Individual = indNames(gc_mlg),
             Population = pop(gc_mlg)),
  mlg.vector(gc_mlg)
)

# Garder seulement les clones (MLG partagés par ≥ 2 individus)
mlg_clone_list <- mlg_full_list[sapply(mlg_full_list, nrow) > 1]

mlg_clone_list   # <-- ta liste de clones MLG avec pop

### 11B) MLL: liste style "mlg.id" mais groupée par MLL + Population
mll_full_list <- split(
  data.frame(Individual = indNames(gc_mll),
             Population = pop(gc_mll)),
  mll(gc_mll)
)

# Garder seulement les clones (MLL partagés par ≥ 2 individus)
mll_clone_list <- mll_full_list[sapply(mll_full_list, nrow) > 1]

mll_clone_list   # <-- ta liste de clones MLL avec pop

### (OPTION) Sortie compacte, une ligne par clone:
# lapply(mll_clone_list, function(x) paste(x$Individual, "(", x$Population, ")", collapse = ", "))
