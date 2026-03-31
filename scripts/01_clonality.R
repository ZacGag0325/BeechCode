# scripts/01_clonality.R
############################################################
# Clonality summary (full dataset: gi)
# Why gi here?
# - Clonality must be quantified on the full (non-clone-corrected) dataset.
# - Exact MLG identity is reported for threshold-free clone counts.
# - Bruvo-based MLL identity is also reported for microsatellite-aware clone
#   lineages under the configured threshold.
#
# Outputs:
# - outputs/tables/clonality_summary.csv
# - outputs/tables/clonality_individual_assignments.csv
############################################################

suppressPackageStartupMessages({
  library(poppr)
  library(adegenet)
  library(dplyr)
  library(ggplot2)
})

source("scripts/_load_objects.R")

message("[01_clonality] Calculating clonality summaries for both MLG and Bruvo-based MLL on gi...")

DEFAULT_BRUVO_MLL_THRESHOLD <- 0.09
DEFAULT_BRUVO_ALGORITHM <- "farthest_neighbor"

find_synonym_col <- function(df, choices) {
  resolve_col_ci(df, choices)
}

compute_mlg_mll_from_gi <- function(gi, threshold = DEFAULT_BRUVO_MLL_THRESHOLD, algorithm = DEFAULT_BRUVO_ALGORITHM) {
  gc_mlg <- poppr::as.genclone(gi)
  mlg_raw <- tryCatch(poppr::mlg.vector(gc_mlg), error = function(e) as.integer(factor(poppr::mlg(gc_mlg))))
  mlg_labels <- paste0("MLG_", as.integer(factor(mlg_raw)))
  
  replen <- rep(2, adegenet::nLoc(gi))
  names(replen) <- adegenet::locNames(gi)
  
  gc_mll <- gc_mlg
  poppr::mlg.filter(
    gc_mll,
    distance = poppr::bruvo.dist,
    replen = replen,
    algorithm = algorithm
  ) <- threshold
  
  mll_raw <- poppr::mll(gc_mll)
  mll_labels <- paste0("MLL_", as.integer(factor(mll_raw)))
  
  list(
    MLG = mlg_labels,
    MLL = mll_labels,
    Bruvo_MLL_threshold = threshold,
    Bruvo_algorithm = algorithm
  )
}

recover_or_recompute_clonality_columns <- function(df_ids_tbl, gi, gi_mll) {
  out <- df_ids_tbl
  
  if (!all(c("MLG", "MLL") %in% names(out))) {
    mlg_syn <- find_synonym_col(out, c("MLG", "mlg", "MLG_id", "mlg_id", "genotype_id"))
    mll_syn <- find_synonym_col(out, c("MLL", "mll", "clone_id", "lineage", "multilocus_lineage"))
    
    if (!is.na(mlg_syn) && mlg_syn != "MLG") names(out)[names(out) == mlg_syn] <- "MLG"
    if (!is.na(mll_syn) && mll_syn != "MLL") names(out)[names(out) == mll_syn] <- "MLL"
    
    if (all(c("MLG", "MLL") %in% names(out))) {
      message("[01_clonality] Recovered MLG/MLL from existing columns.")
    }
  }
  
  if (!all(c("MLG", "MLL") %in% names(out))) {
    message("[01_clonality] Recomputed MLG/MLL from gi because columns were missing in df_ids.")
    cols <- resolve_df_ids_columns(out, context = "[01_clonality]", require = TRUE)
    recomputed <- compute_mlg_mll_from_gi(gi)
    key <- normalize_id(out[[cols$id_col]])
    gi_key <- normalize_id(adegenet::indNames(gi))
    idx <- match(key, gi_key)
    
    if (any(is.na(idx))) {
      stop("[01_clonality] Failed to align df_ids rows to gi while recomputing MLG/MLL.")
    }
    
    out$MLG <- recomputed$MLG[idx]
    out$MLL <- recomputed$MLL[idx]
    out$Bruvo_MLL_threshold <- recomputed$Bruvo_MLL_threshold
    out$Bruvo_algorithm <- recomputed$Bruvo_algorithm
  }
  
  if (!all(c("MLG", "MLL") %in% names(out))) {
    stop("[01_clonality] Failed to recover or recompute MLG/MLL columns.")
  }
  
  n_mll_df <- length(unique(out$MLL[!is.na(out$MLL)]))
  n_mll_obj <- adegenet::nInd(gi_mll)
  if (!identical(n_mll_df, n_mll_obj)) {
    message("[01_clonality] MLL count in df_ids did not match gi_mll; recomputing MLG/MLL from gi.")
    cols <- resolve_df_ids_columns(out, context = "[01_clonality]", require = TRUE)
    recomputed <- compute_mlg_mll_from_gi(gi)
    key <- normalize_id(out[[cols$id_col]])
    gi_key <- normalize_id(adegenet::indNames(gi))
    idx <- match(key, gi_key)
    out$MLG <- recomputed$MLG[idx]
    out$MLL <- recomputed$MLL[idx]
    out$Bruvo_MLL_threshold <- recomputed$Bruvo_MLL_threshold
    out$Bruvo_algorithm <- recomputed$Bruvo_algorithm
    n_mll_df <- length(unique(out$MLL[!is.na(out$MLL)]))
  }
  
  if (!identical(n_mll_df, n_mll_obj)) {
    stop("[01_clonality] df_ids$MLL is inconsistent with gi_mll after recovery/recomputation.")
  }
  
  out
}

make_clone_group_table <- function(assignments_df, clone_col, individual_col = "Individual", site_col = "Site") {
  clone_label <- rlang::sym(clone_col)
  individual_label <- rlang::sym(individual_col)
  has_site <- site_col %in% names(assignments_df)
  
  repeated_groups <- assignments_df %>%
    filter(!is.na(!!clone_label)) %>%
    group_by(!!clone_label) %>%
    mutate(Group_Size = dplyr::n()) %>%
    ungroup() %>%
    filter(Group_Size > 1) %>%
    arrange(!!clone_label, !!individual_label)
  
  if (!has_site) {
    repeated_groups[[site_col]] <- NA_character_
  }
  
  repeated_groups
}

print_separator <- function(char = "=", width = 72) {
  cat(paste(rep(char, width), collapse = ""), "\n", sep = "")
}

print_clone_summary_block <- function(repeated_groups, clone_col, site_available, title) {
  cat("\n")
  print_separator("-", 72)
  cat(title, "\n")
  print_separator("-", 72)
  
  if (nrow(repeated_groups) == 0) {
    cat("No repeated ", clone_col, " groups detected.\n", sep = "")
    return(invisible(NULL))
  }
  
  split_groups <- split(repeated_groups, repeated_groups[[clone_col]])
  
  for (group_name in names(split_groups)) {
    grp <- split_groups[[group_name]]
    cat("\n")
    cat("* ", clone_col, ": ", group_name, "\n", sep = "")
    cat("  Number of individuals: ", nrow(grp), "\n", sep = "")
    cat("  Individuals:\n")
    
    for (i in seq_len(nrow(grp))) {
      site_value <- grp$Site[i]
      if (site_available && !is.na(site_value) && nzchar(site_value)) {
        cat("    - ", grp$Individual[i], " [Site: ", site_value, "]\n", sep = "")
      } else {
        cat("    - ", grp$Individual[i], "\n", sep = "")
      }
    }
  }
  
  invisible(NULL)
}

print_quick_clone_summary <- function(assignments_df, site_available = TRUE) {
  mlg_repeated <- make_clone_group_table(assignments_df, clone_col = "MLG")
  mll_repeated <- make_clone_group_table(assignments_df, clone_col = "MLL")
  
  n_unique_mlg <- dplyr::n_distinct(assignments_df$MLG, na.rm = TRUE)
  n_unique_mll <- dplyr::n_distinct(assignments_df$MLL, na.rm = TRUE)
  n_repeated_mlg_groups <- dplyr::n_distinct(mlg_repeated$MLG, na.rm = TRUE)
  n_repeated_mll_groups <- dplyr::n_distinct(mll_repeated$MLL, na.rm = TRUE)
  n_repeated_mlg_individuals <- nrow(mlg_repeated)
  n_repeated_mll_individuals <- nrow(mll_repeated)
  
  cat("\n")
  print_separator("=", 72)
  cat("QUICK CLONE SUMMARY\n")
  print_separator("=", 72)
  cat("Dataset: gi (full, non-clone-corrected dataset)\n")
  cat("Individuals analysed: ", nrow(assignments_df), "\n", sep = "")
  if (site_available) {
    cat("Site metadata: available\n")
  } else {
    cat("Site metadata: unavailable for one or more individuals; printing IDs only where needed\n")
  }
  
  cat("\nSUMMARY COUNTS\n")
  print_separator("-", 72)
  cat("MLGs\n")
  cat("  - Number of unique MLGs: ", n_unique_mlg, "\n", sep = "")
  cat("  - Number of repeated MLG groups: ", n_repeated_mlg_groups, "\n", sep = "")
  cat("  - Number of individuals involved in repeated MLGs: ", n_repeated_mlg_individuals, "\n", sep = "")
  cat("MLLs\n")
  cat("  - Number of unique MLLs: ", n_unique_mll, "\n", sep = "")
  cat("  - Number of repeated MLL groups: ", n_repeated_mll_groups, "\n", sep = "")
  cat("  - Number of individuals involved in repeated MLLs: ", n_repeated_mll_individuals, "\n", sep = "")
  
  print_clone_summary_block(
    repeated_groups = mlg_repeated,
    clone_col = "MLG",
    site_available = site_available,
    title = "INDIVIDUALS IN DUPLICATED MLGs"
  )
  
  print_clone_summary_block(
    repeated_groups = mll_repeated,
    clone_col = "MLL",
    site_available = site_available,
    title = "INDIVIDUALS IN DUPLICATED MLLs"
  )
  
  cat("\n")
  print_separator("=", 72)
  cat("END QUICK CLONE SUMMARY\n")
  print_separator("=", 72)
  
  invisible(
    list(
      MLG_repeated = mlg_repeated,
      MLL_repeated = mll_repeated,
      counts = list(
        n_unique_mlg = n_unique_mlg,
        n_repeated_mlg_groups = n_repeated_mlg_groups,
        n_repeated_mlg_individuals = n_repeated_mlg_individuals,
        n_unique_mll = n_unique_mll,
        n_repeated_mll_groups = n_repeated_mll_groups,
        n_repeated_mll_individuals = n_repeated_mll_individuals
      )
    )
  )
}

df_ids <- recover_or_recompute_clonality_columns(df_ids, gi, gi_mll)

df_ids_cols <- resolve_df_ids_columns(df_ids, context = "[01_clonality]", require = TRUE)
id_col <- df_ids_cols$id_col
site_col <- df_ids_cols$site_col
threshold_col <- resolve_col_ci(df_ids, c("Bruvo_MLL_threshold", "bruvo_mll_threshold"))
algorithm_col <- resolve_col_ci(df_ids, c("Bruvo_algorithm", "bruvo_algorithm"))

id_key <- normalize_id(df_ids[[id_col]])
site_map <- setNames(as.character(df_ids[[site_col]]), id_key)
mlg_map <- setNames(as.character(df_ids[["MLG"]]), id_key)
mll_map <- setNames(as.character(df_ids[["MLL"]]), id_key)

ind_key <- normalize_id(adegenet::indNames(gi))
site_labels <- site_map[ind_key]
mlg_labels <- mlg_map[ind_key]
mll_labels <- mll_map[ind_key]

site_available <- !any(is.na(site_labels) | !nzchar(site_labels))
if (!site_available) {
  warning("[01_clonality] Site metadata were unavailable for one or more individuals; quick clone summary will print IDs without site where necessary.")
}
if (any(is.na(mlg_labels))) stop("[01_clonality] Could not map all individuals to MLG.")
if (any(is.na(mll_labels))) stop("[01_clonality] Could not map all individuals to MLL.")

bruvo_threshold <- if (!is.na(threshold_col)) unique(stats::na.omit(df_ids[[threshold_col]])) else numeric(0)
bruvo_algorithm <- if (!is.na(algorithm_col)) unique(stats::na.omit(df_ids[[algorithm_col]])) else character(0)

clonality_df <- data.frame(
  Individual = adegenet::indNames(gi),
  Site = ifelse(is.na(site_labels) | !nzchar(site_labels), NA_character_, site_labels),
  MLG = mlg_labels,
  MLL = mll_labels,
  stringsAsFactors = FALSE
)

calc_R <- function(N, G) ifelse(N > 1, (G - 1) / (N - 1), NA_real_)

add_clonality_metrics <- function(dat) {
  dat %>%
    mutate(
      Clonal_Richness_MLG = calc_R(N_individuals, N_MLG),
      Clonal_Richness_MLL = calc_R(N_individuals, N_MLL),
      Genotypic_Richness_MLG = N_MLG / N_individuals,
      Genotypic_Richness_MLL = N_MLL / N_individuals
    )
}

overall <- clonality_df %>%
  summarise(
    N_individuals = dplyr::n(),
    N_MLG = dplyr::n_distinct(MLG),
    N_MLL = dplyr::n_distinct(MLL)
  ) %>%
  add_clonality_metrics() %>%
  mutate(Level = "overall", Site = "ALL") %>%
  select(Level, Site, everything())

by_site <- clonality_df %>%
  mutate(Site = ifelse(is.na(Site) | !nzchar(Site), "SITE_UNAVAILABLE", Site)) %>%
  group_by(Site) %>%
  summarise(
    N_individuals = dplyr::n(),
    N_MLG = dplyr::n_distinct(MLG),
    N_MLL = dplyr::n_distinct(MLL),
    .groups = "drop"
  ) %>%
  add_clonality_metrics() %>%
  mutate(Level = "site") %>%
  select(Level, Site, everything())

clonality_summary <- bind_rows(overall, by_site) %>%
  mutate(
    N = N_individuals,
    G = N_MLL,
    Clonal_Richness_R = Clonal_Richness_MLL
  )

if (length(bruvo_threshold) > 0) {
  clonality_summary$Bruvo_MLL_threshold <- bruvo_threshold[1]
}
if (length(bruvo_algorithm) > 0) {
  clonality_summary$Bruvo_algorithm <- bruvo_algorithm[1]
}

out_file <- file.path(TABLES_DIR, "clonality_summary.csv")
write.csv(clonality_summary, out_file, row.names = FALSE)

if (length(bruvo_threshold) > 0) {
  clonality_df$Bruvo_MLL_threshold <- bruvo_threshold[1]
}
if (length(bruvo_algorithm) > 0) {
  clonality_df$Bruvo_algorithm <- bruvo_algorithm[1]
}

assign_file <- file.path(TABLES_DIR, "clonality_individual_assignments.csv")
write.csv(clonality_df, assign_file, row.names = FALSE)

print_quick_clone_summary(clonality_df, site_available = site_available)

# -----------------------------------------------------------------------------
# MLG barplots by site (FR + EN) for presentation
# -----------------------------------------------------------------------------
# Objective:
# - x axis: Site
# - y axis: Number of distinct MLGs observed in each site
# - reference line at theoretical maximum = 24 sampled individuals per site
# - robust behavior even when some sites have <24 valid individuals

find_latitude_col <- function(df) {
  find_synonym_col(
    df,
    c("Latitude", "latitude", "lat", "LAT", "Lat", "y", "Y")
  )
}

build_site_mlg_summary <- function(assignments_df, df_ids_tbl = NULL) {
  summary_tbl <- assignments_df %>%
    mutate(Site = ifelse(is.na(Site) | !nzchar(Site), "SITE_UNAVAILABLE", Site)) %>%
    group_by(Site) %>%
    summarise(
      N_MLG = dplyr::n_distinct(MLG, na.rm = TRUE),
      N_individuals = dplyr::n(),
      .groups = "drop"
    )
  
  if (!is.null(df_ids_tbl)) {
    lat_col <- find_latitude_col(df_ids_tbl)
    if (!is.na(lat_col)) {
      df_ids_cols_local <- resolve_df_ids_columns(df_ids_tbl, context = "[01_clonality]", require = TRUE)
      lat_by_site <- df_ids_tbl %>%
        transmute(
          Site = as.character(.data[[df_ids_cols_local$site_col]]),
          Latitude = suppressWarnings(as.numeric(.data[[lat_col]]))
        ) %>%
        filter(!is.na(Site), nzchar(Site)) %>%
        group_by(Site) %>%
        summarise(
          Latitude = if (all(is.na(Latitude))) NA_real_ else mean(Latitude, na.rm = TRUE),
          .groups = "drop"
        )
      
      summary_tbl <- summary_tbl %>%
        left_join(lat_by_site, by = "Site")
    } else {
      summary_tbl$Latitude <- NA_real_
    }
  } else {
    summary_tbl$Latitude <- NA_real_
  }
  
  if ("Latitude" %in% names(summary_tbl) && any(!is.na(summary_tbl$Latitude))) {
    # Sort south -> north when latitude is available
    summary_tbl <- summary_tbl %>%
      arrange(Latitude, Site)
  } else {
    summary_tbl <- summary_tbl %>%
      arrange(Site)
  }
  
  summary_tbl$Site <- factor(summary_tbl$Site, levels = summary_tbl$Site)
  summary_tbl
}

make_mlg_barplot <- function(summary_tbl, lang = c("fr", "en"), max_theoretical = 24) {
  lang <- match.arg(lang)
  
  labels <- switch(
    lang,
    fr = list(
      title = "Nombre de génotypes multilocus (MLG) par site",
      subtitle = "Maximum théorique de 24 individus échantillonnés par site",
      x = "Site",
      y = "Nombre de MLG"
    ),
    en = list(
      title = "Number of multilocus genotypes (MLGs) per site",
      subtitle = "Theoretical maximum of 24 sampled individuals per site",
      x = "Site",
      y = "Number of MLGs"
    )
  )
  
  ggplot(summary_tbl, aes(x = Site, y = N_MLG)) +
    geom_col(fill = "#4C78A8", width = 0.75) +
    geom_hline(
      yintercept = max_theoretical,
      linetype = "dashed",
      linewidth = 0.8,
      color = "#595959"
    ) +
    scale_y_continuous(
      limits = c(0, max(max_theoretical, summary_tbl$N_MLG, na.rm = TRUE) * 1.05),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(
      title = labels$title,
      subtitle = labels$subtitle,
      x = labels$x,
      y = labels$y
    ) +
    theme_minimal(base_size = 16) +
    theme(
      plot.title = element_text(face = "bold", size = 19),
      plot.subtitle = element_text(size = 14),
      axis.title = element_text(size = 15),
      axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1),
      axis.text.y = element_text(size = 12),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank()
    )
}

save_mlg_plot_dual_language <- function(summary_tbl, fig_dir, max_theoretical = 24) {
  plot_fr <- make_mlg_barplot(summary_tbl, lang = "fr", max_theoretical = max_theoretical)
  plot_en <- make_mlg_barplot(summary_tbl, lang = "en", max_theoretical = max_theoretical)
  
  png_fr <- file.path(fig_dir, "mlg_per_site_barplot_fr.png")
  pdf_fr <- file.path(fig_dir, "mlg_per_site_barplot_fr.pdf")
  png_en <- file.path(fig_dir, "mlg_per_site_barplot_en.png")
  pdf_en <- file.path(fig_dir, "mlg_per_site_barplot_en.pdf")
  
  ggsave(filename = png_fr, plot = plot_fr, width = 12, height = 7, dpi = 320)
  ggsave(filename = pdf_fr, plot = plot_fr, width = 12, height = 7)
  ggsave(filename = png_en, plot = plot_en, width = 12, height = 7, dpi = 320)
  ggsave(filename = pdf_en, plot = plot_en, width = 12, height = 7)
  
  list(
    fr_png = png_fr, fr_pdf = pdf_fr,
    en_png = png_en, en_pdf = pdf_en
  )
}

site_mlg_summary <- build_site_mlg_summary(clonality_df, df_ids_tbl = df_ids)

cat("\n[01_clonality] Tableau résumé utilisé pour le barplot MLG par site:\n")
print(site_mlg_summary)

site_mlg_summary_file <- file.path(TABLES_DIR, "mlg_per_site_summary.csv")
write.csv(site_mlg_summary, site_mlg_summary_file, row.names = FALSE)

mlg_plot_files <- save_mlg_plot_dual_language(
  summary_tbl = site_mlg_summary,
  fig_dir = FIGURES_DIR,
  max_theoretical = 24
)

message("[01_clonality] Saved: ", out_file)
message("[01_clonality] Saved: ", assign_file)
message("[01_clonality] Saved: ", site_mlg_summary_file)
message("[01_clonality] Saved: ", mlg_plot_files$fr_png)
message("[01_clonality] Saved: ", mlg_plot_files$fr_pdf)
message("[01_clonality] Saved: ", mlg_plot_files$en_png)
message("[01_clonality] Saved: ", mlg_plot_files$en_pdf)