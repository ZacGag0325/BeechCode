# scripts/_load_objects.R
############################################################
# Shared loader for BeechCode population genetics workflow
# - Finds project root
# - Loads gi, gi_mll, df_ids, meta from outputs/v1/objects
# - Ensures site/pop labels are available
# - Creates standardized output directories
############################################################

suppressPackageStartupMessages({
  library(poppr)
  library(adegenet)
  library(vegan)
  library(hierfstat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# ----------------------------
# 1) Find project root
# ----------------------------
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

PROJECT_ROOT <- find_project_root()
setwd(PROJECT_ROOT)

# ----------------------------
# 2) Standard output folders
# ----------------------------
OUTPUT_DIR   <- file.path(PROJECT_ROOT, "outputs")
TABLES_DIR   <- file.path(OUTPUT_DIR, "tables")
FIGURES_DIR  <- file.path(OUTPUT_DIR, "figures")
MATRICES_DIR <- file.path(OUTPUT_DIR, "matrices")

TABLES_SUPP_DIR  <- file.path(TABLES_DIR, "supplementary")
FIGURES_SUPP_DIR <- file.path(FIGURES_DIR, "supplementary")

for (d in c(OUTPUT_DIR, TABLES_DIR, FIGURES_DIR, MATRICES_DIR, TABLES_SUPP_DIR, FIGURES_SUPP_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ----------------------------
# 3) df_ids schema helpers
# ----------------------------
DF_IDS_ID_CHOICES <- c("ind", "individual", "sample", "sampleid", "id", "ind_id")
DF_IDS_SITE_CHOICES <- c("Site", "site", "pop", "population")

normalize_id <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("\uFEFF", "", x, fixed = TRUE)
  x <- gsub("[[:cntrl:]]", "", x)
  toupper(x)
}

resolve_col_ci <- function(df, choices) {
  nms <- names(df)
  idx <- match(TRUE, tolower(nms) %in% tolower(choices), nomatch = 0)
  if (idx == 0) return(NA_character_)
  nms[idx]
}

resolve_df_ids_columns <- function(df_ids_tbl, context = "[df_ids]", require = TRUE) {
  id_col <- resolve_col_ci(df_ids_tbl, DF_IDS_ID_CHOICES)
  site_col <- resolve_col_ci(df_ids_tbl, DF_IDS_SITE_CHOICES)
  
  if (require && (is.na(id_col) || is.na(site_col))) {
    stop(
      context,
      " df_ids must contain an individual ID column and a Site column.",
      "\nAccepted ID columns: ", paste(DF_IDS_ID_CHOICES, collapse = ", "),
      "\nAccepted Site columns: ", paste(DF_IDS_SITE_CHOICES, collapse = ", ")
    )
  }
  
  list(id_col = id_col, site_col = site_col)
}

# ----------------------------
# 4) Load canonical objects
# ----------------------------
OBJ_DIR <- file.path(PROJECT_ROOT, "outputs", "v1", "objects")
required_files <- c("gi.rds", "gi_mll.rds", "df_ids.rds", "meta.rds")
missing_files <- required_files[!file.exists(file.path(OBJ_DIR, required_files))]

if (length(missing_files) > 0) {
  stop(
    "Missing required object files in ", OBJ_DIR,
    "\nMissing: ", paste(missing_files, collapse = ", "),
    "\nRun 00_master_pipeline.R first."
  )
}

gi     <- readRDS(file.path(OBJ_DIR, "gi.rds"))
gi_mll <- readRDS(file.path(OBJ_DIR, "gi_mll.rds"))
df_ids <- readRDS(file.path(OBJ_DIR, "df_ids.rds"))
meta   <- readRDS(file.path(OBJ_DIR, "meta.rds"))

df_ids_cols <- resolve_df_ids_columns(df_ids, context = "[_load_objects]", require = TRUE)

# ----------------------------
# 5) Site resolver and pop assignment
# ----------------------------
resolve_site_for_genind <- function(gobj, df_ids_tbl, cols = NULL) {
  p <- adegenet::pop(gobj)
  if (!is.null(p) && length(p) == adegenet::nInd(gobj) && all(!is.na(p))) {
    return(as.factor(p))
  }
  
  if (is.null(cols)) {
    cols <- resolve_df_ids_columns(df_ids_tbl, context = "[_load_objects]", require = TRUE)
  }
  
  id_vec <- as.character(df_ids_tbl[[cols$id_col]])
  site_vec <- as.character(df_ids_tbl[[cols$site_col]])
  
  mapper <- setNames(site_vec, normalize_id(id_vec))
  sites <- mapper[normalize_id(adegenet::indNames(gobj))]
  
  if (any(is.na(sites))) {
    missing_ids <- adegenet::indNames(gobj)[is.na(sites)]
    stop(
      "Could not map Site for all individuals in genind object. Missing count: ",
      sum(is.na(sites)),
      "\nExample missing IDs: ",
      paste(head(missing_ids, 10), collapse = ", ")
    )
  }
  
  as.factor(sites)
}

adegenet::pop(gi) <- resolve_site_for_genind(gi, df_ids, cols = df_ids_cols)
adegenet::pop(gi_mll) <- resolve_site_for_genind(gi_mll, df_ids, cols = df_ids_cols)

# ----------------------------
# 6) Console summary
# ----------------------------
message("[_load_objects] Loaded objects from: ", OBJ_DIR)
message("[_load_objects] nInd(gi) = ", adegenet::nInd(gi), ", nLoc(gi) = ", adegenet::nLoc(gi))
message("[_load_objects] nInd(gi_mll) = ", adegenet::nInd(gi_mll), ", nLoc(gi_mll) = ", adegenet::nLoc(gi_mll))
message("[_load_objects] Detected df_ids ID column: ", df_ids_cols$id_col)
message("[_load_objects] Detected df_ids Site column: ", df_ids_cols$site_col)
message("[_load_objects] Output directories ready:")
message("  - ", TABLES_DIR)
message("  - ", FIGURES_DIR)
message("  - ", MATRICES_DIR)
message("  - ", TABLES_SUPP_DIR)
message("  - ", FIGURES_SUPP_DIR)