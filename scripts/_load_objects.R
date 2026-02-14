# scripts/_load_objects.R
############################################################
# scripts/_load_objects.R
# Loads gi/df_ids/meta/gi_mll from outputs/v1/objects
############################################################

suppressPackageStartupMessages(library(adegenet))

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

PROJECT_ROOT <- find_project_root()
setwd(PROJECT_ROOT)

RUN_TAG <- "v1"
RUN_OUT <- file.path(PROJECT_ROOT, "outputs", RUN_TAG)
OBJ_DIR <- file.path(RUN_OUT, "objects")

need <- c("gi.rds", "df_ids.rds", "meta.rds", "gi_mll.rds")
missing <- need[!file.exists(file.path(OBJ_DIR, need))]
if (length(missing) > 0) {
  stop(
    "Missing saved objects in: ", OBJ_DIR,
    "\nMissing: ", paste(missing, collapse = ", "),
    "\nFix: run 00_master_pipeline.R first so outputs/v1/objects is fully created."
  )
}

gi     <- readRDS(file.path(OBJ_DIR, "gi.rds"))
df_ids <- readRDS(file.path(OBJ_DIR, "df_ids.rds"))
meta   <- readRDS(file.path(OBJ_DIR, "meta.rds"))
gi_mll <- readRDS(file.path(OBJ_DIR, "gi_mll.rds"))

cat("Loaded objects:\n")
cat(" - nInd(gi)     =", nInd(gi), "\n")
cat(" - nLoc(gi)     =", nLoc(gi), "\n")
cat(" - df_ids rows  =", nrow(df_ids), "\n")
cat(" - meta rows    =", nrow(meta), "\n")
cat(" - nInd(gi_mll) =", nInd(gi_mll), "\n")
