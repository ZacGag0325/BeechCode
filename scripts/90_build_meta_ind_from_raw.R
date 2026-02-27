############################################################
# scripts/90_build_meta_ind_from_raw.R
# Build individual-level databases required by STRUCTURE mapping
############################################################

options(repos = c(CRAN = "https://cloud.r-project.org"))

pkgs <- c("readxl", "readr", "dplyr", "stringi")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

suppressPackageStartupMessages({
  library(readxl)
  library(readr)
  library(dplyr)
  library(stringi)
})

find_project_root <- function() {
  cur <- normalizePath(getwd(), mustWork = FALSE)
  repeat {
    if (file.exists(file.path(cur, "00_master_pipeline.R"))) return(cur)
    parent <- dirname(cur)
    if (identical(parent, cur)) break
    cur <- parent
  }
  stop("MISSING FROM YOUR SIDE:\n- Cannot find project root (00_master_pipeline.R).")
}

to_utf8 <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  if (!is.character(x)) return(x)
  y <- suppressWarnings(iconv(x, from = "", to = "UTF-8", sub = ""))
  y[is.na(y)] <- ""
  y
}

clean_names <- function(x) {
  x <- to_utf8(as.character(x))
  x <- trimws(x)
  x <- stringi::stri_trans_general(x, "Latin-ASCII")
  x <- gsub("\\s+", "_", x)
  x <- gsub("[^a-zA-Z0-9_]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  tolower(x)
}

read_any_tabular <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("xlsx", "xls")) {
    x <- readxl::read_excel(path)
  } else if (ext == "csv") {
    x <- suppressWarnings(readr::read_csv(path, show_col_types = FALSE, guess_max = 50000))
  } else {
    x <- utils::read.table(path, header = TRUE, fill = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  }
  x <- as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE)
  names(x) <- clean_names(names(x))
  for (j in seq_along(x)) {
    if (is.character(x[[j]]) || is.factor(x[[j]])) x[[j]] <- to_utf8(as.character(x[[j]]))
  }
  x
}

choose_col <- function(nms, candidates) {
  nms <- clean_names(nms)
  cands <- clean_names(candidates)
  i <- which(nms %in% cands)
  if (length(i) == 0) return(NA_character_)
  nms[i[1]]
}

PROJECT_ROOT <- find_project_root()
RAW_DIR <- file.path(PROJECT_ROOT, "data", "raw")
IN_DIR <- file.path(PROJECT_ROOT, "inputs")
OUT_IDS <- file.path(PROJECT_ROOT, "outputs", "v1", "structure_runs")
dir.create(IN_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_IDS, recursive = TRUE, showWarnings = FALSE)

candidates <- c(
  file.path(RAW_DIR, "poppr.xlsx"),
  file.path(RAW_DIR, "poppr.csv"),
  file.path(RAW_DIR, "genotypes.xlsx"),
  file.path(RAW_DIR, "genotypes.csv")
)
source_file <- candidates[file.exists(candidates)][1]
if (is.na(source_file)) {
  stop(
    "MISSING FROM YOUR SIDE:\n",
    "- No genotype raw file found in data/raw/.\n",
    "- Provide one of: poppr.xlsx, poppr.csv, genotypes.xlsx, genotypes.csv."
  )
}

cat("[build_meta_ind] source:", source_file, "\n")
dat <- read_any_tabular(source_file)

id_candidates <- c("ind_id", "id", "ind", "sample", "sample_id", "nom_labo_echantillons", "nom_labo_echantillon")
site_candidates <- c("site", "population", "pop", "numero_population", "num_population", "site_id")

id_col <- choose_col(names(dat), id_candidates)
site_col <- choose_col(names(dat), site_candidates)

if (is.na(id_col) || is.na(site_col)) {
  stop(
    "MISSING FROM YOUR SIDE:\n",
    "- Cannot identify ID/site columns in raw genotype file.\n",
    "- Found columns: ", paste(names(dat), collapse = ", "), "\n",
    "- Expected minimum columns: ind_id + site (or synonyms)."
  )
}

meta_ind <- dat %>%
  transmute(
    ind_id = to_utf8(as.character(.data[[id_col]])),
    site = to_utf8(as.character(.data[[site_col]]))
  ) %>%
  mutate(ind_id = trimws(ind_id), site = trimws(site)) %>%
  filter(ind_id != "", site != "")

if (nrow(meta_ind) == 0) {
  stop("MISSING FROM YOUR SIDE:\n- Extracted meta_ind is empty after cleaning.")
}

meta_ind <- meta_ind %>% distinct(ind_id, .keep_all = TRUE)

meta_ind_path <- file.path(IN_DIR, "meta_ind.csv")
ids_order_path <- file.path(OUT_IDS, "ids_order_from_raw.csv")

readr::write_csv(meta_ind, meta_ind_path)
readr::write_csv(meta_ind %>% transmute(ind_id = ind_id), ids_order_path)

cat("[build_meta_ind] written:\n")
cat(" - ", meta_ind_path, " (n=", nrow(meta_ind), ")\n", sep = "")
cat(" - ", ids_order_path, " (n=", nrow(meta_ind), ")\n", sep = "")
cat("[build_meta_ind] columns used: id_col=", id_col, " | site_col=", site_col, "\n", sep = "")
