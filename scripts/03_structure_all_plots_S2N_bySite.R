############################################################
# 03_structure_all_plots_S2N_bySite.R  (BeechCode)
############################################################

### 0) PROJECT ROOT + PATHS (ROBUST) ----
suppressPackageStartupMessages(library(here))
options(repos = c(CRAN = "https://cloud.r-project.org"))

FALLBACK_ROOT  <- file.path(path.expand("~"), "Desktop", "BeechCode")
candidate_root <- here::here()

if (dir.exists(file.path(candidate_root, "data", "raw"))) {
  PROJECT_ROOT <- candidate_root
} else if (dir.exists(file.path(FALLBACK_ROOT, "data", "raw"))) {
  PROJECT_ROOT <- FALLBACK_ROOT
  setwd(PROJECT_ROOT)
} else {
  stop("Can't find project root.\nFix: open your BeechCode.Rproj OR change FALLBACK_ROOT.")
}

OUTPUTS_DIR <- file.path(PROJECT_ROOT, "outputs")
RUN_TAG <- "v1"
RUN_OUT <- file.path(OUTPUTS_DIR, RUN_TAG)

STRUCTURE_RUNS_DIR <- file.path(RUN_OUT, "structure_runs")
Q_DIR   <- file.path(STRUCTURE_RUNS_DIR, "Q_extracted")
PLOT_DIR <- file.path(RUN_OUT, "structure_plots")

dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)

message("Project root: ", PROJECT_ROOT)
message("RUN_OUT: ", RUN_OUT)
message("Q_DIR: ", Q_DIR)
message("PLOT_DIR: ", PLOT_DIR)

### 1) PACKAGES ----
pkgs <- c("ggplot2", "patchwork", "readxl")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(readxl))

### 2) FIND STRUCTURE RESULTS FOLDER (where *_f files are) ----
cand1 <- file.path(STRUCTURE_RUNS_DIR, "STRUCTURE_ZG_HEG-RESULTS")
cand2 <- file.path(STRUCTURE_RUNS_DIR, "STRUCTURE_ZG_HEG-results")

if (dir.exists(cand1)) {
  RESULTS_DIR <- cand1
} else if (dir.exists(cand2)) {
  RESULTS_DIR <- cand2
} else {
  stop("Could not find STRUCTURE results folder under: ", STRUCTURE_RUNS_DIR)
}
message("RESULTS_DIR: ", RESULTS_DIR)

### 3) LOAD FIELD METADATA from EXCEL DATABASE (genetique sheet) ----
META_FILE  <- file.path(PROJECT_ROOT, "data", "raw", "donnees_modifiees_west_summer2024 copie.xlsx")
META_SHEET <- "genetique"

if (!file.exists(META_FILE)) {
  stop("Excel metadata file not found at:\n  ", META_FILE,
       "\nFix: put the file in data/raw or edit META_FILE.")
}

meta_xl <- readxl::read_excel(META_FILE, sheet = META_SHEET)

req_xl <- c("Nom_Labo_Échantillons", "Site", "Latitude", "Longitude")
missing_xl <- setdiff(req_xl, names(meta_xl))
if (length(missing_xl) > 0) {
  stop("Your Excel sheet '", META_SHEET, "' is missing columns:\n  ",
       paste(missing_xl, collapse = ", "),
       "\nFound columns:\n  ",
       paste(names(meta_xl), collapse = ", "))
}

meta <- data.frame(
  ind       = trimws(as.character(meta_xl[["Nom_Labo_Échantillons"]])),
  Site      = trimws(as.character(meta_xl[["Site"]])),
  Latitude  = as.numeric(meta_xl[["Latitude"]]),
  Longitude = as.numeric(meta_xl[["Longitude"]]),
  stringsAsFactors = FALSE
)

ind_COL  <- "ind"
SITE_COL <- "Site"
LAT_COL  <- "Latitude"
LON_COL  <- "Longitude"

cat("Loaded Excel metadata:\n  ", META_FILE, " (sheet: ", META_SHEET, ")\n", sep = "")
cat("Rows:", nrow(meta), " | Unique inds:", length(unique(meta[[ind_COL]])), "\n")

if (anyNA(meta[[LAT_COL]])) {
  bad <- meta[is.na(meta[[LAT_COL]]), ind_COL]
  stop("Some Latitude values are NA in metadata (first 10 inds): ", paste(head(bad, 10), collapse = ", "))
}
if (anyNA(meta[[LON_COL]])) {
  bad <- meta[is.na(meta[[LON_COL]]), ind_COL]
  stop("Some Longitude values are NA in metadata (first 10 inds): ", paste(head(bad, 10), collapse = ", "))
}
dup_inds <- meta[[ind_COL]][duplicated(meta[[ind_COL]])]
if (length(dup_inds) > 0) {
  stop("Duplicate ind IDs found in metadata (first 10): ",
       paste(head(unique(dup_inds), 10), collapse = ", "))
}

### 4) RECOVER ind ORDER FROM ONE STRUCTURE *_f FILE ----
f_files <- list.files(RESULTS_DIR, recursive = TRUE, full.names = TRUE, pattern = "_f$")
if (length(f_files) == 0) stop("No *_f files found under: ", RESULTS_DIR)

pick_f <- function(ff) {
  k1 <- ff[grepl("K1", ff, ignore.case = TRUE)]
  if (length(k1) > 0) return(k1[1])
  ff[1]
}
f_example <- pick_f(f_files)
message("Using this _f file to recover ind order: ", basename(f_example))

meta_inds <- unique(meta[[ind_COL]])

extract_inds_from_f <- function(f, valid_inds) {
  x <- readLines(f, warn = FALSE)
  
  # find ancestry section header
  indx <- grep("Inferred ancestry of indiv", x, ignore.case = TRUE)
  if (length(indx) == 0) indx <- grep("Inferred ancestry|Membership coefficients|Estimated membership", x, ignore.case = TRUE)
  if (length(indx) == 0) stop("Could not find ancestry section in: ", basename(f))
  
  i0 <- indx[1]
  inds <- character(0)
  
  for (j in (i0 + 1):length(x)) {
    ln <- x[j]
    if (!grepl(":", ln, fixed = TRUE)) next
    
    if (length(inds) > 0 && trimws(ln) == "") break
    
    pos <- gregexpr(":", ln, fixed = TRUE)[[1]]
    if (pos[1] == -1) next
    last_colon <- tail(pos, 1)
    lhs <- trimws(substr(ln, 1, last_colon - 1))
    
    lhs2 <- gsub("\\(.*?\\)", " ", lhs)
    tokens <- strsplit(trimws(lhs2), "\\s+")[[1]]
    
    hit <- tokens[tokens %in% valid_inds]
    if (length(hit) >= 1) {
      inds <- c(inds, hit[1])
    } else {
      if (length(tokens) >= 2 && tokens[2] %in% valid_inds) {
        inds <- c(inds, tokens[2])
      }
    }
  }
  
  inds <- inds[inds %in% valid_inds]
  inds <- inds[!duplicated(inds)]
  inds
}

inds_in_Q_order <- extract_inds_from_f(f_example, meta_inds)
message("Recovered ", length(inds_in_Q_order), " inds in STRUCTURE row order.")

if (length(inds_in_Q_order) != nrow(meta)) {
  stop(paste0(
    "Extracted ", length(inds_in_Q_order), " valid inds from the _f file, but metadata has ",
    nrow(meta), ".\n",
    "This usually means the ancestry block parsing didn't match the file format.\n",
    "Try a different _f file or check the ancestry section in: ", f_example
  ))
}

### 5) BUILD SOUTH TO NORTH ORDER (by SITE mean latitude, then by Latitude) ----
meta2 <- meta[match(inds_in_Q_order, meta[[ind_COL]]), , drop = FALSE]

site_means <- aggregate(meta2[[LAT_COL]], by = list(Site = meta2[[SITE_COL]]), FUN = mean)
names(site_means)[2] <- "MeanLat"
site_order <- site_means$Site[order(site_means$MeanLat)]

message("Population/Site order (South to North):")
print(site_order)

ord <- order(match(meta2[[SITE_COL]], site_order), meta2[[LAT_COL]], meta2[[LON_COL]])

site_sorted <- meta2[[SITE_COL]][ord]
boundary_pos <- which(site_sorted[-1] != site_sorted[-length(site_sorted)]) + 0.5

### 6) READ Q FILES + GROUP BY K ----
qfiles <- list.files(Q_DIR, pattern = "\\.Q$", full.names = TRUE)
if (length(qfiles) == 0) stop("No .Q files found in: ", Q_DIR)

read_Q <- function(f) {
  Q <- as.matrix(read.table(f, header = FALSE))
  storage.mode(Q) <- "numeric"
  if (nrow(Q) != length(inds_in_Q_order)) {
    stop("Row mismatch in ", basename(f), ": Q has ", nrow(Q),
         " rows but extracted ind order has ", length(inds_in_Q_order))
  }
  Q
}

K_map <- list()
for (f in qfiles) {
  Q <- read_Q(f)
  K <- ncol(Q)
  key <- as.character(K)
  if (is.null(K_map[[key]])) K_map[[key]] <- list()
  K_map[[key]][[length(K_map[[key]]) + 1]] <- list(file = f, Q = Q)
}

available_K <- sort(as.integer(names(K_map)))
message("Available K values: ", paste(available_K, collapse = ", "))

### 7) PICK BEST REPLICATE PER K ----
score_Q <- function(Q) mean(apply(Q, 1, max), na.rm = TRUE)
pick_best_run <- function(entries) {
  scores <- vapply(entries, function(e) score_Q(e$Q), numeric(1))
  entries[[which.max(scores)]]
}

### 8) PLOT FUNCTION (ordered S->N + site boundaries) ----
plot_Q <- function(Q, K, show_legend = TRUE, subtitle_txt = NULL) {
  Q <- Q[ord, , drop = FALSE]
  n <- nrow(Q)
  
  df <- data.frame(
    Ind = rep(seq_len(n), times = K),
    Cluster = factor(rep(seq_len(K), each = n), levels = seq_len(K)),
    Q = as.vector(Q)
  )
  
  p <- ggplot(df, aes(x = Ind, y = Q, fill = Cluster)) +
    geom_col(width = 1) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(
      title = paste0("K = ", K),
      subtitle = subtitle_txt,
      x = NULL,
      y = "Ancestry proportion",
      fill = "Cluster"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.title = element_text(size = 11, face = "bold"),
      plot.subtitle = element_text(size = 9),
      legend.title = element_text(size = 9),
      legend.text  = element_text(size = 8)
    )
  
  if (length(boundary_pos) > 0) {
    p <- p + geom_vline(xintercept = boundary_pos, linewidth = 0.3)
  }
  
  if (!show_legend) p <- p + theme(legend.position = "none")
  p
}

### 9) SAVE PER-K + COMBINED + GRIND ----
save_per_k <- function(Ks_to_plot, tag = "Kset") {
  chosen_log <- c(paste0("Chosen best replicate per K for ", tag, " (mean max Q heuristic):"))
  
  for (K in Ks_to_plot) {
    key <- as.character(K)
    if (is.null(K_map[[key]])) {
      message("Skipping K=", K, " (not found).")
      next
    }
    
    best <- pick_best_run(K_map[[key]])
    subtitle_txt <- "Individuals ordered South to North (by Site mean latitude; within Site by Latitude)"
    p <- plot_Q(best$Q, K, show_legend = TRUE, subtitle_txt = subtitle_txt)
    
    pdf_file <- file.path(PLOT_DIR, paste0(tag, "_K", K, "_best_S2N.pdf"))
    png_file <- file.path(PLOT_DIR, paste0(tag, "_K", K, "_best_S2N.png"))
    
    grDevices::pdf(pdf_file, width = 11, height = 4)
    print(p)
    grDevices::dev.off()
    
    grDevices::png(png_file, width = 2200, height = 800, res = 200)
    print(p)
    grDevices::dev.off()
    
    chosen_log <- c(chosen_log, paste0("K=", K, " -> ", basename(best$file)))
    message("Saved ", tag, " K=", K)
  }
  
  log_file <- file.path(PLOT_DIR, paste0(tag, "_chosen_runs_S2N.txt"))
  writeLines(chosen_log, log_file)
  message("Saved log: ", log_file)
}

# Per-K outputs
save_per_k(1:8,  tag = "K1to8")
save_per_k(1:10, tag = "K1to10")

# Combined multi-page PDF for K=1..8
combined_pdf <- file.path(PLOT_DIR, "K1to8_all_pages_S2N.pdf")
grDevices::pdf(combined_pdf, width = 11, height = 4)
for (K in 1:8) {
  key <- as.character(K)
  if (is.null(K_map[[key]])) next
  best <- pick_best_run(K_map[[key]])
  print(plot_Q(best$Q, K, show_legend = TRUE,
               subtitle_txt = "Individuals ordered South to North"))
}
grDevices::dev.off()
message("Saved combined multi-page PDF: ", combined_pdf)

# Grind figure K=1..10 (2 columns)
plots_grind <- list()
chosen_grind <- c("Chosen best replicate per K for GRIND K=1..10:")

for (K in 1:10) {
  key <- as.character(K)
  if (is.null(K_map[[key]])) next
  best <- pick_best_run(K_map[[key]])
  plots_grind[[key]] <- plot_Q(best$Q, K, show_legend = TRUE,
                               subtitle_txt = "South to North ordering")
  chosen_grind <- c(chosen_grind, paste0("K=", K, " -> ", basename(best$file)))
}

grind_plot <- wrap_plots(plots_grind, ncol = 2) +
  plot_annotation(
    title = "STRUCTURE / ADMIXTURE barplots (K = 1 to 10)",
    subtitle = "Individuals ordered South to North using site mean latitude",
    theme = theme(plot.title = element_text(size = 14, face = "bold"))
  )

grind_pdf <- file.path(PLOT_DIR, "STRUCTURE_grind_K1to10_S2N.pdf")
grind_png <- file.path(PLOT_DIR, "STRUCTURE_grind_K1to10_S2N.png")
grind_log <- file.path(PLOT_DIR, "STRUCTURE_grind_K1to10_S2N_chosen_runs.txt")

ggsave(grind_pdf, grind_plot, width = 11, height = 14)
ggsave(grind_png, grind_plot, width = 11, height = 14, dpi = 200)
writeLines(chosen_grind, grind_log)

message("Saved grind PDF: ", grind_pdf)
message("Saved grind PNG: ", grind_png)
message("Saved grind log: ", grind_log)

message("DONE. Outputs in: ", PLOT_DIR)
