############################################################
## POPPR EXCEL (.xlsx) -> GENEPOP (1 file per population)
## For ML-Relate (needs 1 or 2 POP blocks per file)
############################################################

############################################################
### 0. WORKING DIRECTORY ------------------------------------
############################################################
setwd("/Users/zacharygagnon/Desktop/R data workshop")
getwd()
list.files()   # should show poppr.xlsx

############################################################
### 1. PACKAGES ---------------------------------------------
############################################################
pkgs <- c("readxl", "dplyr", "stringr", "tibble")
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install) > 0) install.packages(to_install, quiet = TRUE)

library(readxl)
library(dplyr)
library(stringr)
library(tibble)

############################################################
### 2. INPUT / OUTPUT ---------------------------------------
############################################################
infile   <- "poppr.xlsx"
out_dir  <- "MLRelate_byPop"   # folder where files will be written
dir.create(out_dir, showWarnings = FALSE)

############################################################
### 3. READ EXCEL -------------------------------------------
############################################################
dat <- read_excel(infile, sheet = 1) %>% as.data.frame()

############################################################
### 4. SET THE ID + POP COLUMNS (YOUR FILE) -----------------
############################################################
id_col  <- "Nom_Labo_Échantillons"
pop_col <- "Numéro_Population"

if (!id_col %in% names(dat))  stop(paste0("❌ Can't find ID column: ", id_col))
if (!pop_col %in% names(dat)) stop(paste0("❌ Can't find POP column: ", pop_col))

dat[[id_col]]  <- as.character(dat[[id_col]])
dat[[pop_col]] <- as.character(dat[[pop_col]])

############################################################
### 5. DETECT LOCUS ALLELE PAIRS (ROBUST) -------------------
############################################################
find_locus_pairs <- function(dat, id_col, pop_col) {
  
  cols <- names(dat)
  cols <- setdiff(cols, c(id_col, pop_col))
  cols <- cols[!grepl("^Unnamed", cols, ignore.case = TRUE)]
  
  # allele-1 columns end with _1 or .1
  a1_cols <- cols[grepl("(_1$)|(\\.1$)", cols)]
  
  if (length(a1_cols) == 0) {
    stop("❌ I couldn't find any allele-1 columns ending with _1 or .1.")
  }
  
  out <- list()
  
  for (a1 in a1_cols) {
    
    if (grepl("_1$", a1)) {
      base <- sub("_1$", "", a1)
      a2   <- paste0(base, "_2")
    } else {
      base <- sub("\\.1$", "", a1)
      a2   <- paste0(base, ".2")
    }
    
    if (a2 %in% cols) {
      out[[length(out) + 1]] <- data.frame(
        locus = base,
        a1col = a1,
        a2col = a2,
        stringsAsFactors = FALSE
      )
    }
  }
  
  if (length(out) == 0) {
    stop("❌ No allele pairs detected. Check column names for _1/_2 or .1/.2 pairs.")
  }
  
  dplyr::bind_rows(out) %>%
    dplyr::distinct(locus, .keep_all = TRUE)
}

pairs <- find_locus_pairs(dat, id_col, pop_col)
loci  <- pairs$locus

############################################################
### 6. MAKE 6-DIGIT GENEPOP GENOTYPES -----------------------
############################################################
pad3 <- function(x) {
  x <- trimws(as.character(x))
  x[x == "" | is.na(x)] <- NA
  x <- ifelse(is.na(x), NA, stringr::str_extract(x, "\\d+"))
  x <- suppressWarnings(as.integer(x))
  ifelse(is.na(x), "000", sprintf("%03d", x))
}

geno <- dat[, c(id_col, pop_col), drop = FALSE]

for (i in seq_len(nrow(pairs))) {
  loc <- pairs$locus[i]
  a1  <- pad3(dat[[pairs$a1col[i]]])
  a2  <- pad3(dat[[pairs$a2col[i]]])
  geno[[loc]] <- paste0(a1, a2)  # aaabbb format (e.g., 112/114 -> 112114)
}

geno <- geno %>% arrange(.data[[pop_col]], .data[[id_col]])
pops <- unique(geno[[pop_col]])

############################################################
### 7. WRITE ONE GENEPOP FILE PER POPULATION ----------------
############################################################
safe_name <- function(x) {
  x <- as.character(x)
  x <- str_replace_all(x, "[^A-Za-z0-9._-]+", "_")
  x <- str_replace_all(x, "_+", "_")
  x <- str_replace_all(x, "^_|_$", "")
  ifelse(nchar(x) == 0, "POP", x)
}

write_genepop_onepop <- function(sub, pop_value, loci, out_path) {
  
  con <- file(out_path, open = "wt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  
  # Title line
  writeLines(paste0("Microsatellite genotypes for Fagus grandifolia - Pop ", pop_value), con)
  writeLines("", con)
  
  # Locus names (one per line)
  writeLines(loci, con)
  writeLines("", con)
  
  # Single POP block
  writeLines("POP", con)
  
  for (r in seq_len(nrow(sub))) {
    id <- sub[[id_col]][r]
    g  <- as.character(sub[r, loci, drop = TRUE])
    line <- paste0(id, " , ", paste(g, collapse = " "))
    writeLines(line, con)
  }
  
  invisible(TRUE)
}

written <- character(0)

for (p in pops) {
  sub <- geno[geno[[pop_col]] == p, , drop = FALSE]
  
  out_file <- file.path(out_dir, paste0("MLRelate_pop_", safe_name(p), ".txt"))
  write_genepop_onepop(sub = sub, pop_value = p, loci = loci, out_path = out_file)
  
  written <- c(written, out_file)
}

message("✅ Done! Wrote ", length(written), " files into: ", normalizePath(out_dir))
message("ℹ️ Pops: ", paste(pops, collapse = ", "))
message("ℹ️ # loci: ", length(loci))

############################################################
### 8. OPTIONAL: write a 2-pop file (ML-Relate also allows 2)
############################################################
# Example: make one file containing exactly TWO pops (p1 and p2)
# p1 <- pops[1]
# p2 <- pops[2]
# two <- geno[geno[[pop_col]] %in% c(p1, p2), , drop = FALSE]
#
# out2 <- file.path(out_dir, paste0("MLRelate_2pops_", safe_name(p1), "_", safe_name(p2), ".txt"))
# con <- file(out2, open="wt", encoding="UTF-8"); on.exit(close(con), add=TRUE)
# writeLines(paste0("Fagus grandifolia - Pops ", p1, " & ", p2), con); writeLines("", con)
# writeLines(loci, con); writeLines("", con)
# for (pp in c(p1, p2)) {
#   writeLines("POP", con)
#   sub <- two[two[[pop_col]] == pp, , drop=FALSE]
#   for (r in seq_len(nrow(sub))) {
#     line <- paste0(sub[[id_col]][r], " , ", paste(as.character(sub[r, loci, drop=TRUE]), collapse=" "))
#     writeLines(line, con)
#   }
# }
# message("✅ Wrote 2-pop file: ", out2)
