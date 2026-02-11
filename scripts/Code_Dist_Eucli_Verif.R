##############################################
# Détection totalement indépendante des clones avec la distance euclidienne
##############################################
#fais une fois
library(readxl)
library(dplyr)
library(tibble)

##############################################
# 1. Lire ton fichier 
##############################################

### 1. DOSSIER DE TRAVAIL ----
setwd("/Users/zacharygagnon/Desktop/R data workshop")

# Lire le fichier réel
geno_raw <- read_xlsx("tbd.xlsx", sheet = "Genotypes-Indiv.")#à changer le nom

# Vérifier les colonnes si c'est OK
cat("Colonnes trouvées :", ncol(geno_raw), "\n")
print(names(geno_raw))

##############################################
# 2. Extraire les colonnes 
##############################################

# ID = première colonne = "Nom_Labo_Échantillons" = Identification de chaques tiges
id <- geno_raw[[1]]

# Colonnes d'allèles (toutes sauf les 2 premières car ID + #population de 1 à 12)
allele_df <- geno_raw[, -(1:2)]

# Conversion en numérique
X <- allele_df %>%
  mutate(across(everything(), as.numeric)) %>%
  as.matrix()

##############################################
# 3. Filtrer les individus complet (donc ceux avec allele completes) pour la fonction dist()
##############################################

keep <- complete.cases(X)
Xc  <- X[keep, , drop = FALSE]
idc <- id[keep]

cat("Individus complets :", length(idc), "\n")

##############################################
# 4. DISTANCE EUCLIDIENNE entre chaque pair d'indiv.
##############################################

#creation
D  <- dist(Xc, method = "euclidean")
Dm <- as.matrix(D)

rownames(Dm) <- idc
colnames(Dm) <- idc

##############################################
# 5. Paires de clones (distance = 0)
##############################################

tol <- 1e-9 #rapporte a 0 si infime
idx <- which(Dm <= tol, arr.ind = TRUE)

# enlever doublons et diagonale pour avoir un ind chaque
idx <- idx[idx[,1] < idx[,2], , drop = FALSE]

#creation de tableau de paire de clones
clone_pairs <- data.frame(
  ind1 = rownames(Dm)[idx[,1]],
  ind2 = colnames(Dm)[idx[,2]],
  distance = Dm[idx],
  row.names = NULL
)
#nombre de pair clonale
cat("Nombre de paires de clones :", nrow(clone_pairs), "\n")
clone_pairs
