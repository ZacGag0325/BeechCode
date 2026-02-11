###############################
#   Clonality / MLG - Fagus  #
#   Zach - nov 2025          #
###############################

### 0. PACKAGES ----
library(adegenet)
library(poppr)
library(readxl)
library(dplyr)
library(tidyr)
library(hierfstat)

### 0B. OPTION : SAUVEGARDES SUR DISQUE ? ----
save_outputs <- FALSE   # mets TRUE si tu veux écrire les .csv

### 1. DOSSIER DE TRAVAIL ----
setwd("/Users/zacharygagnon/Desktop/R data workshop")
getwd()
list.files()  # tu devrais voir tbd.xlsx

### 2. LECTURE DU FICHIER D’ALLÈLES ----
geno_raw <- read_xlsx("tbd.xlsx")

# ID en caractère
geno_raw$Nom_Labo_Échantillons <- as.character(geno_raw$Nom_Labo_Échantillons)

head(geno_raw)
str(geno_raw)

### 3. RETIRER LA POP ML1 ----
# Vérifie les codes de population
table(geno_raw$Numéro_Population)

# Ici tu as mis ML1 = 11
ML1_code <- 11

geno <- geno_raw %>%
  filter(Numéro_Population != ML1_code)

table(geno$Numéro_Population)

### 4. COLONNES D’ALLÈLES ----
# col 1 = Nom_Labo_Échantillons
# col 2 = Numéro_Population
# col 3+ = allèles

allele_cols <- 3:ncol(geno)

# Sécurité : vérifier qu’on a bien un nombre pair de colonnes d’allèles
if (length(allele_cols) %% 2 != 0) {
  stop("Le nombre de colonnes d'allèles n'est pas pair. Vérifie ton fichier.")
}

# Indices des premières colonnes de chaque locus (1ère allèle)
locus_starts <- allele_cols[seq(1, length(allele_cols), by = 2)]

### 5. CONSTRUIRE UN TABLEAU "1 COLONNE PAR LOCUS" ----
# Chaque colonne = génotype "allèle1/allèle2" (texte)

geno_locus <- data.frame(row.names = geno$Nom_Labo_Échantillons)

for (i in seq_along(locus_starts)) {
  col1_idx <- locus_starts[i]
  col2_idx <- col1_idx + 1
  
  col1_name <- names(geno)[col1_idx]
  col2_name <- names(geno)[col2_idx]
  
  # Nom du locus = nom de la première colonne, sans le ".1" éventuel
  locus_name <- gsub("\\.1$", "", col1_name)
  
  # Construire le génotype allélique "A1/A2"
  a1 <- as.character(geno[[col1_name]])
  a2 <- as.character(geno[[col2_name]])
  
  geno_locus[[locus_name]] <- paste(a1, a2, sep = "/")
}

# Vérifier la structure
head(geno_locus)
str(geno_locus)

### 5C. NETTOYER LES NOMS DES COLONNES ----
clean_names <- names(geno_locus)
clean_names <- gsub("\\.+", "_", clean_names)  # remplace les "." par "_"
names(geno_locus) <- clean_names
names(geno_locus)

### 6. CRÉER L’OBJET genind ----
gen <- df2genind(
  geno_locus,
  ploidy    = 2,
  ind.names = rownames(geno_locus),
  pop       = as.factor(geno$Numéro_Population),
  type      = "codom",
  sep       = "/",   # séparateur entre les deux allèles
  ncode     = 3      # longueur des codes alléliques (115, 152, etc.)
)

gen
summary(gen)

### 7. PASSER EN genclone ET CALCULER LES MLG (naïfs) ----
gc_mlg <- as.genclone(gen)
gc_mlg

# IMPORTANT : vector d’assignation des MLG par individu
mlg_vec <- mlg.vector(gc_mlg)
head(mlg_vec)
length(unique(mlg_vec))    # nombre de MLG uniques

# Table des tailles de clones (MLG naïfs)
mlg_tab <- mlg.table(gc_mlg)
mlg_tab

# Liste des individus par clone
mlg_list <- mlg.id(gc_mlg)

### 8. AJOUTER LES MLG AU TABLEAU GÉNOTYPIQUE ----
geno$MLG <- as.integer(mlg_vec)

head(geno[, c("Nom_Labo_Échantillons", "Numéro_Population", "MLG")])

if (save_outputs) {
  write.csv(
    geno,
    "genotypes_avec_MLG_sans_ML1.csv",
    row.names = FALSE
  )
}

### 9. STATISTIQUES CLONALES PAR POPULATION (MLG NAÏFS) ----
clone_stats_mlg <- poppr(gc_mlg)
clone_stats_mlg

if (save_outputs) {
  write.csv(
    clone_stats_mlg,
    "clone_stats_par_population_sans_ML1_MLG.csv",
    row.names = TRUE
  )
}

### 10. RICHESSE ALLÉLIQUE RAREFIÉE (N = 14) ----
hf <- genind2hierfstat(gen)
ar <- hierfstat::allelic.richness(hf, min.n = 14)
ar

if (save_outputs) {
  write.csv(
    ar$Ar,
    "allelic_richness_rarefied_min14.csv",
    row.names = TRUE
  )
}

###############################
#   11. MLL AVEC BRUVO.DIST  #
###############################

# ICI on veut regrouper des MLG très proches en lignées clonales (MLL)
# via Bruvo (stepwise mutation), donc parfait pour tes microsats.

# 11A. Vecteur des longueurs de répétition par locus ----
# A ADAPTER selon le motif de chaque locus.
# Si (temporairement) tu assumes que tous sont di-nucléotidiques :
replen <- rep(2, nLoc(gen))
names(replen) <- locNames(gen)
replen

# (Plus tard, tu pourras mettre par ex. c(2,2,3,2,4,...) selon tes repeats.)

# 11B. Explorer les seuils possibles avec filter_stats ----
# Ça calcule des distances et regarde à quels seuils des MLG se fusionnent.
fs <- filter_stats(
  gc_mlg,
  distance  = bruvo.dist,
  replen    = replen,
  plot      = TRUE  # donne un graphique interactif dans R
)

# Proposer un seuil automatique pour l’algo "farthest neighbor"
farthest_thresh <- cutoff_predictor(fs$farthest$THRESHOLDS)
farthest_thresh  # regarde la valeur, ajuste si besoin

# 11C. Créer un objet cloné pour les MLL ----
gc_mll <- gc_mlg

# On applique le filtrage des MLG par distance Bruvo pour définir les MLL
mlg.filter(
  gc_mll,
  distance  = bruvo.dist,
  replen    = replen,
  algorithm = "farthest_neighbor"  # conservatif; tu peux tester "nearest" aussi
) <- farthest_thresh

gc_mll  # devrait dire "contracted multilocus genotypes"

# 11D. Regarder combien de MLL tu as maintenant ----
nmll(gc_mlg)   # nombre de MLG naïfs
nmll(gc_mll)   # nombre de MLL (après collapse)

# Table des MLL
mll_tab <- mlg.table(gc_mll)
mll_tab

# Assignation des MLL par individu
mll_vec <- mll(gc_mll)     # vecteur de lignées clonales
head(mll_vec)

# Ajouter au tableau génotypique
geno$MLL <- as.integer(mll_vec)

head(geno[, c("Nom_Labo_Échantillons", "Numéro_Population", "MLG", "MLL")])

if (save_outputs) {
  write.csv(
    geno,
    "genotypes_avec_MLG_MLL_sans_ML1.csv",
    row.names = FALSE
  )
}

# 11E. Statistiques clonalité par population avec MLL ----
clone_stats_mll <- poppr(gc_mll)
clone_stats_mll

if (save_outputs) {
  write.csv(
    clone_stats_mll,
    "clone_stats_par_population_sans_ML1_MLL.csv",
    row.names = TRUE
  )
}

###############################
# Fin du script MLG / MLL
###############################
