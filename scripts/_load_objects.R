# scripts/_load_objects.R
############################################################
# Shared loader for BeechCode population genetics workflow
# - Finds project root
# - Loads gi, gi_mll, df_ids, meta from outputs/v1/objects
# - Ensures site/pop labels are available
# - Creates standardized output directories
############################################################

suppressPackageStartupMessages({
  library(adegenet)
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
# 3) Load canonical objects
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

# ----------------------------
# 4) Site resolver and pop assignment
# ----------------------------
resolve_col <- function(df, choices) {
  nms <- names(df)
  nms_low <- tolower(nms)
  idx <- match(TRUE, nms_low %in% tolower(choices), nomatch = 0)
  if (idx == 0) return(NA_character_)
  nms[idx]
}

resolve_site_for_genind <- function(gobj, df_ids_tbl) {
  p <- adegenet::pop(gobj)
  if (!is.null(p) && length(p) == adegenet::nInd(gobj) && all(!is.na(p))) {
    return(as.factor(p))
  }
  
  id_col <- resolve_col(df_ids_tbl, c("ind", "individual", "sample", "sampleid", "id"))
  site_col <- resolve_col(df_ids_tbl, c("site", "population", "pop"))
  
  if (is.na(id_col) || is.na(site_col)) {
    stop("df_ids must include an individual ID column and a Site column.")
  }
  
  mapper <- setNames(as.character(df_ids_tbl[[site_col]]), as.character(df_ids_tbl[[id_col]]))
  sites <- mapper[adegenet::indNames(gobj)]
  
  if (any(is.na(sites))) {
    missing_n <- sum(is.na(sites))
    stop("Could not map Site for all individuals in genind object. Missing count: ", missing_n)
  }
  
  as.factor(sites)
}

adegenet::pop(gi) <- resolve_site_for_genind(gi, df_ids)
adegenet::pop(gi_mll) <- resolve_site_for_genind(gi_mll, df_ids)

# ----------------------------
# 5) Console summary
# ----------------------------
message("[_load_objects] Loaded objects from: ", OBJ_DIR)
message("[_load_objects] nInd(gi) = ", adegenet::nInd(gi), ", nLoc(gi) = ", adegenet::nLoc(gi))
message("[_load_objects] nInd(gi_mll) = ", adegenet::nInd(gi_mll), ", nLoc(gi_mll) = ", adegenet::nLoc(gi_mll))
message("[_load_objects] Output directories ready:")
message("  - ", TABLES_DIR)
message("  - ", FIGURES_DIR)
message("  - ", MATRICES_DIR)
message("  - ", TABLES_SUPP_DIR)
message("  - ", FIGURES_SUPP_DIR)