# scripts/excavation_validation.R
############################################################
# Excavation validation analysis for physically excavated pairs.
#
# This analysis is intentionally separate from the main BeechCode 2024
# population-genetic workflow. Excavation/LG2 samples are not added to the
# main diversity, STRUCTURE, AMOVA, PCA, DAPC, Jost's D, Mantel, or
# population-level clonality analyses.
#
# Important ID policy:
# - Field stem IDs (for example Id_Tige_1 / Id_Tige_2) are not assumed to be
#   laboratory genotype IDs.
# - Genetic comparisons are only made after explicit lab_id_1 / lab_id_2
#   mapping is available either in the genetic_excavation sheet or in a filled
#   outputs/excavation/excavation_field_to_lab_id_mapping_TEMPLATE.csv file.
# - The script does not infer field-to-lab genotype links from row order.
############################################################

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(janitor)
  library(ggplot2)
  library(readr)
  library(poppr)
  library(adegenet)
  library(officer)
  library(flextable)
})

FIELD_FILE <- file.path("data", "raw", "donnees_modifiees_west_summer2024 copie.xlsx")
FIELD_SHEET <- "genetic_excavation"
GENO_FILE <- file.path("data", "raw", "poppr_excavation.xlsx")
GENO_SHEET <- "Excav_Sample"
OUTPUT_DIR <- file.path("outputs", "excavation")
MAPPING_TEMPLATE_FILE <- file.path(OUTPUT_DIR, "excavation_field_to_lab_id_mapping_TEMPLATE.csv")
AVAILABLE_LAB_IDS_FILE <- file.path(OUTPUT_DIR, "excavation_available_lab_ids.csv")
MISSINGNESS_THRESHOLD <- 0.35
BRUVO_MLL_THRESHOLD <- 0.09
BRUVO_ALGORITHM <- "farthest_neighbor"
SAMPLE_ID_COL_RAW <- "Nom_Labo_Échantillons"

# The excavation validation uses the same 15-locus set requested for the
# microsatellite workflow. DZ447_A_0 is intentionally not included; if present
# in a workbook it is ignored because it was excluded from the main analysis.
REQUESTED_LOCI <- c(
  "sfc_0036", "csolfagus_29", "csolfagus_31", "csolfagus_19", "FS1-15F",
  "csolfagus_06", "sfc_1143", "csolfagus_05", "EMILY_A_0", "FCM5",
  "EEU75_A_0", "FG5", "EJV8T_A_0", "DUKCT_A_0", "ERHBI_A_0"
)

LAB_ID_1_CANDIDATES <- c("lab_id_1", "labo_id_1", "id_labo_1", "genotype_id_1", "poppr_id_1")
LAB_ID_2_CANDIDATES <- c("lab_id_2", "labo_id_2", "id_labo_2", "genotype_id_2", "poppr_id_2")
missing_tokens <- c("", "na", "n/a", "null", ".", "-", "?", "0", "-9")

stop_if_missing_file_or_sheet <- function(path, sheet) {
  if (!file.exists(path)) {
    stop("Input file not found: ", path, call. = FALSE)
  }
  sheets <- readxl::excel_sheets(path)
  if (!(sheet %in% sheets)) {
    stop(
      "Sheet '", sheet, "' was not found in ", path, ". Available sheets: ",
      paste(sheets, collapse = ", "),
      call. = FALSE
    )
  }
}

clean_value <- function(x) {
  out <- str_trim(as.character(x))
  out[is.na(x) | str_to_lower(out) %in% c("", "na", "n/a", "null")] <- NA_character_
  out
}

is_missing_allele <- function(x) {
  x_chr <- str_trim(as.character(x))
  is.na(x) | str_to_lower(x_chr) %in% missing_tokens
}

find_column_exact_or_clean <- function(raw_names, raw_target) {
  exact <- which(raw_names == raw_target)
  if (length(exact) > 0) return(raw_names[exact[1]])
  
  cleaned_names <- janitor::make_clean_names(raw_names)
  target_clean <- janitor::make_clean_names(raw_target)
  cleaned_match <- which(cleaned_names == target_clean)
  if (length(cleaned_match) > 0) return(raw_names[cleaned_match[1]])
  
  NA_character_
}

resolve_first_existing_column <- function(nms, candidates) {
  hit <- candidates[candidates %in% nms]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

resolve_allele_col <- function(nms, locus, suffix) {
  candidates <- c(paste0(locus, "_", suffix), paste0(locus, ".", suffix))
  hit <- candidates[candidates %in% nms]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

validate_requested_loci <- function(geno_raw) {
  nms <- names(geno_raw)
  locus_pairs <- tibble(
    locus = REQUESTED_LOCI,
    allele_1_col = vapply(REQUESTED_LOCI, function(loc) resolve_allele_col(nms, loc, "1"), character(1)),
    allele_2_col = vapply(REQUESTED_LOCI, function(loc) resolve_allele_col(nms, loc, "2"), character(1))
  )
  
  missing_cols <- locus_pairs %>% filter(is.na(allele_1_col) | is.na(allele_2_col))
  if (nrow(missing_cols) > 0) {
    stop(
      "The Excav_Sample sheet is missing paired allele columns for requested loci: ",
      paste(missing_cols$locus, collapse = ", "), ". Available columns: ",
      paste(nms, collapse = ", "),
      call. = FALSE
    )
  }
  
  locus_pairs
}

standardize_connection <- function(x) {
  x_chr <- str_to_lower(str_trim(as.character(x)))
  case_when(
    is.na(x) | x_chr %in% c("", "na", "n/a", "unclear", "incertain") ~ "unclear",
    x_chr %in% c("oui", "yes", "y", "o") ~ "yes",
    x_chr %in% c("non", "no", "n") ~ "no",
    TRUE ~ "unclear"
  )
}

read_genotype_data <- function() {
  stop_if_missing_file_or_sheet(GENO_FILE, GENO_SHEET)
  
  geno_raw <- readxl::read_excel(GENO_FILE, sheet = GENO_SHEET, .name_repair = "minimal") %>%
    as.data.frame(check.names = FALSE, stringsAsFactors = FALSE)
  sample_id_col <- find_column_exact_or_clean(names(geno_raw), SAMPLE_ID_COL_RAW)
  if (is.na(sample_id_col)) {
    stop(
      "Required sample ID column '", SAMPLE_ID_COL_RAW, "' was not found in ", GENO_FILE,
      " sheet '", GENO_SHEET, "'. Available columns: ", paste(names(geno_raw), collapse = ", "),
      call. = FALSE
    )
  }
  
  locus_pairs <- validate_requested_loci(geno_raw)
  geno_raw$sample_id_for_analysis <- clean_value(geno_raw[[sample_id_col]])
  geno_raw <- geno_raw %>% filter(!is.na(sample_id_for_analysis), sample_id_for_analysis != "")
  
  list(data = geno_raw, sample_id_col = sample_id_col, locus_pairs = locus_pairs)
}

calculate_sample_missingness <- function(geno_raw, locus_pairs) {
  n_missing_loci <- integer(nrow(geno_raw))
  for (i in seq_len(nrow(locus_pairs))) {
    m1 <- is_missing_allele(geno_raw[[locus_pairs$allele_1_col[i]]])
    m2 <- is_missing_allele(geno_raw[[locus_pairs$allele_2_col[i]]])
    n_missing_loci <- n_missing_loci + as.integer(m1 & m2)
  }
  
  tibble(
    sample_id = geno_raw$sample_id_for_analysis,
    n_loci_total = nrow(locus_pairs),
    n_loci_missing = n_missing_loci,
    missingness_prop = n_loci_missing / n_loci_total,
    missingness_percent = missingness_prop * 100,
    failing_missingness_threshold = missingness_prop > MISSINGNESS_THRESHOLD,
    retained_at_35_percent = !failing_missingness_threshold
  )
}

write_available_lab_ids <- function(sample_missingness) {
  available_lab_ids <- sample_missingness %>%
    transmute(
      lab_id = sample_id,
      missingness_percent,
      retained_at_35_percent,
      n_loci_missing,
      n_loci_total
    ) %>%
    arrange(lab_id)
  
  readr::write_csv(available_lab_ids, AVAILABLE_LAB_IDS_FILE)
  available_lab_ids
}

read_and_clean_field_data <- function() {
  stop_if_missing_file_or_sheet(FIELD_FILE, FIELD_SHEET)
  
  field_raw <- readxl::read_excel(FIELD_FILE, sheet = FIELD_SHEET, .name_repair = "minimal") %>%
    as.data.frame(check.names = FALSE, stringsAsFactors = FALSE)
  field <- janitor::clean_names(field_raw)
  
  expected_clean <- c(
    "site", "id_paire", "id_tige_1", "id_tige_2", "circ_tige_1_cm", "dhp_tige_1_cm",
    "circ_tige_2_cm", "dhp_tige_2_cm", "distance_cm", "connexion_sout", "photo", "notes"
  )
  missing_field_cols <- setdiff(expected_clean, names(field))
  if (length(missing_field_cols) > 0) {
    stop(
      "The field sheet '", FIELD_SHEET, "' is missing required columns after clean_names(): ",
      paste(missing_field_cols, collapse = ", "), ". Available cleaned columns: ",
      paste(names(field), collapse = ", "),
      call. = FALSE
    )
  }
  
  lab_id_1_source <- resolve_first_existing_column(names(field), LAB_ID_1_CANDIDATES)
  lab_id_2_source <- resolve_first_existing_column(names(field), LAB_ID_2_CANDIDATES)
  has_sheet_mapping <- !is.na(lab_id_1_source) && !is.na(lab_id_2_source)
  
  out <- field %>%
    transmute(
      site = as.character(site),
      field_pair_id = as.character(id_paire),
      stem_1_field_id = as.character(id_tige_1),
      stem_2_field_id = as.character(id_tige_2),
      circ_stem_1_cm = suppressWarnings(as.numeric(circ_tige_1_cm)),
      dhp_stem_1_cm = suppressWarnings(as.numeric(dhp_tige_1_cm)),
      circ_stem_2_cm = suppressWarnings(as.numeric(circ_tige_2_cm)),
      dhp_stem_2_cm = suppressWarnings(as.numeric(dhp_tige_2_cm)),
      distance_cm = suppressWarnings(as.numeric(distance_cm)),
      distance_m = distance_cm / 100,
      connection_observed = standardize_connection(connexion_sout),
      photo = as.character(photo),
      notes = as.character(notes)
    ) %>%
    mutate(
      excavation_pair_id = paste0("EXC_", str_pad(row_number(), 2, pad = "0")),
      .before = site
    )
  
  if (has_sheet_mapping) {
    out <- out %>%
      mutate(
        lab_id_1 = clean_value(field[[lab_id_1_source]]),
        lab_id_2 = clean_value(field[[lab_id_2_source]])
      )
    attr(out, "mapping_source") <- paste0("genetic_excavation columns: ", lab_id_1_source, ", ", lab_id_2_source)
  } else {
    out <- attach_mapping_from_template_or_stop(out)
  }
  
  out
}

mapping_template_columns <- function(field) {
  field %>%
    transmute(
      excavation_pair_id,
      site,
      field_pair_id,
      stem_1_field_id,
      stem_2_field_id,
      distance_cm,
      distance_m,
      connection_observed,
      lab_id_1 = "",
      lab_id_2 = "",
      notes
    )
}

attach_mapping_from_template_or_stop <- function(field) {
  template <- mapping_template_columns(field)
  
  if (file.exists(MAPPING_TEMPLATE_FILE)) {
    mapping <- readr::read_csv(MAPPING_TEMPLATE_FILE, show_col_types = FALSE, col_types = cols(.default = col_character()))
    required_cols <- names(template)
    missing_mapping_cols <- setdiff(required_cols, names(mapping))
    if (length(missing_mapping_cols) > 0) {
      stop(
        "Existing mapping template is missing required columns: ",
        paste(missing_mapping_cols, collapse = ", "), ". Recreate or repair ", MAPPING_TEMPLATE_FILE,
        call. = FALSE
      )
    }
    
    mapping <- mapping %>%
      mutate(
        lab_id_1 = clean_value(lab_id_1),
        lab_id_2 = clean_value(lab_id_2)
      )
    
    if (any(!is.na(mapping$lab_id_1) | !is.na(mapping$lab_id_2))) {
      field_key <- field$excavation_pair_id
      idx <- match(field_key, mapping$excavation_pair_id)
      if (any(is.na(idx))) {
        stop(
          "The filled mapping template does not contain all excavation_pair_id values from the current field sheet. Missing: ",
          paste(field_key[is.na(idx)], collapse = ", "),
          call. = FALSE
        )
      }
      field$lab_id_1 <- mapping$lab_id_1[idx]
      field$lab_id_2 <- mapping$lab_id_2[idx]
      attr(field, "mapping_source") <- MAPPING_TEMPLATE_FILE
      return(field)
    }
  }
  
  readr::write_csv(template, MAPPING_TEMPLATE_FILE)
  stop(
    "No lab ID mapping columns were found in genetic_excavation. Fill outputs/excavation/excavation_field_to_lab_id_mapping_TEMPLATE.csv or add lab_id_1 and lab_id_2 columns to the genetic_excavation sheet, then rerun the script.",
    call. = FALSE
  )
}

validate_allele_data_present <- function(geno_raw, locus_pairs) {
  allele_cols <- c(locus_pairs$allele_1_col, locus_pairs$allele_2_col)
  allele_matrix_has_data <- any(vapply(geno_raw[allele_cols], function(x) any(!is_missing_allele(x)), logical(1)))
  if (!allele_matrix_has_data) {
    stop(
      "poppr_excavation.xlsx has sample IDs but no allele data in the Excav_Sample sheet. Add genotype allele values before running excavation validation.",
      call. = FALSE
    )
  }
}

make_id_check <- function(field, genotype_ids) {
  field %>%
    transmute(
      excavation_pair_id, field_pair_id, site, stem_1_field_id, stem_2_field_id,
      lab_id_1, lab_id_2,
      lab_id_1_found = !is.na(lab_id_1) & lab_id_1 %in% genotype_ids,
      lab_id_2_found = !is.na(lab_id_2) & lab_id_2 %in% genotype_ids,
      id_check_status = case_when(
        is.na(lab_id_1) | is.na(lab_id_2) ~ "one_or_both_lab_ids_unmapped",
        lab_id_1_found & lab_id_2_found ~ "both_lab_ids_found",
        lab_id_1_found | lab_id_2_found ~ "one_lab_id_found",
        TRUE ~ "neither_lab_id_found"
      )
    )
}

build_genind <- function(geno_raw, locus_pairs) {
  geno <- data.frame(row.names = geno_raw$sample_id_for_analysis, check.names = FALSE)
  
  for (i in seq_len(nrow(locus_pairs))) {
    locus <- locus_pairs$locus[i]
    locus_safe <- janitor::make_clean_names(locus)
    a1 <- str_trim(as.character(geno_raw[[locus_pairs$allele_1_col[i]]]))
    a2 <- str_trim(as.character(geno_raw[[locus_pairs$allele_2_col[i]]]))
    miss <- is_missing_allele(a1) | is_missing_allele(a2)
    g <- paste(a1, a2, sep = "/")
    g[miss] <- NA_character_
    geno[[locus_safe]] <- g
  }
  
  adegenet::df2genind(
    geno,
    ploidy = 2,
    ind.names = rownames(geno),
    type = "codom",
    sep = "/",
    ncode = 3
  )
}

get_mlg_labels <- function(gi) {
  gc <- poppr::as.genclone(gi)
  raw <- tryCatch(poppr::mlg.vector(gc), error = function(e) as.integer(factor(poppr::mlg(gc))))
  setNames(paste0("MLG_", as.integer(factor(raw))), adegenet::indNames(gi))
}

calculate_bruvo_distance_matrix <- function(gi) {
  # Repeat lengths were not available in the project files supplied for this
  # standalone excavation script, so all loci use repeat length = 2. This is
  # the same conservative fallback used elsewhere in BeechCode when locus-
  # specific repeat lengths are unavailable.
  replen <- rep(2, adegenet::nLoc(gi))
  names(replen) <- adegenet::locNames(gi)
  as.matrix(poppr::bruvo.dist(gi, replen = replen))
}

pairwise_results <- function(field, sample_missingness, mlg_labels, bruvo_matrix) {
  ids_in_geno <- sample_missingness$sample_id
  miss_lookup <- setNames(sample_missingness$failing_missingness_threshold, sample_missingness$sample_id)
  
  field %>%
    rowwise() %>%
    mutate(
      lab_id_1_found = !is.na(lab_id_1) && lab_id_1 %in% ids_in_geno,
      lab_id_2_found = !is.na(lab_id_2) && lab_id_2 %in% ids_in_geno,
      sample_1_fails_missingness = if (lab_id_1_found) isTRUE(miss_lookup[[lab_id_1]]) else NA,
      sample_2_fails_missingness = if (lab_id_2_found) isTRUE(miss_lookup[[lab_id_2]]) else NA,
      both_samples_pass_missingness = lab_id_1_found && lab_id_2_found &&
        identical(sample_1_fails_missingness, FALSE) && identical(sample_2_fails_missingness, FALSE),
      bruvo_distance = if (
        lab_id_1_found && lab_id_2_found && lab_id_1 %in% rownames(bruvo_matrix) && lab_id_2 %in% colnames(bruvo_matrix)
      ) {
        bruvo_matrix[lab_id_1, lab_id_2]
      } else {
        NA_real_
      },
      same_mlg = case_when(
        !both_samples_pass_missingness ~ "uncertain",
        mlg_labels[lab_id_1] == mlg_labels[lab_id_2] ~ "yes",
        TRUE ~ "no"
      ),
      same_mll = case_when(
        !both_samples_pass_missingness | is.na(bruvo_distance) ~ "uncertain",
        bruvo_distance <= BRUVO_MLL_THRESHOLD ~ "yes",
        TRUE ~ "no"
      ),
      clone_call = case_when(
        !both_samples_pass_missingness | is.na(bruvo_distance) ~ "uncertain",
        same_mll == "yes" ~ "clone",
        same_mll == "no" ~ "not_clone",
        TRUE ~ "uncertain"
      )
    ) %>%
    ungroup()
}

write_word_table <- function(results) {
  word_tbl <- results %>%
    select(
      excavation_pair_id, field_pair_id, site, stem_1_field_id, stem_2_field_id,
      lab_id_1, lab_id_2, distance_m, connection_observed,
      bruvo_distance, same_mlg, same_mll, clone_call, notes
    ) %>%
    mutate(
      distance_m = round(distance_m, 3),
      bruvo_distance = round(bruvo_distance, 4)
    )
  
  ft <- flextable::flextable(word_tbl) %>%
    flextable::autofit()
  
  doc <- officer::read_docx() %>%
    officer::body_add_par("Excavation validation table", style = "heading 1") %>%
    flextable::body_add_flextable(ft)
  
  print(doc, target = file.path(OUTPUT_DIR, "excavation_validation_table.docx"))
}

save_plot_both <- function(plot, basename, width = 7, height = 5) {
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
  pdf_file <- file.path(OUTPUT_DIR, paste0(basename, ".pdf"))
  png_file <- file.path(OUTPUT_DIR, paste0(basename, ".png"))
  
  ggplot2::ggsave(pdf_file, plot, width = width, height = height, device = "pdf")
  ggplot2::ggsave(png_file, plot, width = width, height = height, dpi = 300)
  
  c(pdf_file, png_file)
}

make_figures <- function(results) {
  plot_data <- results %>%
    filter(!is.na(bruvo_distance)) %>%
    mutate(
      connection_observed = factor(connection_observed, levels = c("yes", "no", "unclear")),
      clone_call = factor(clone_call, levels = c("clone", "not_clone", "uncertain"))
    )
  
  theme_exc <- theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(), legend.position = "right")
  
  fig1 <- ggplot(plot_data, aes(x = connection_observed, fill = clone_call)) +
    geom_bar(color = "grey20") +
    scale_x_discrete(drop = FALSE) +
    labs(x = "Field-observed underground connection", y = "Number of excavated pairs", fill = "Genetic result") +
    theme_exc
  figure_files <- save_plot_both(fig1, "excavation_connection_vs_clone_call")
  
  fig2 <- ggplot(plot_data, aes(x = connection_observed, y = bruvo_distance)) +
    geom_point(aes(color = clone_call), position = position_jitter(width = 0.08, height = 0), size = 2.5, alpha = 0.85) +
    scale_x_discrete(drop = FALSE) +
    labs(x = "Field-observed underground connection", y = "Bruvo distance", color = "Genetic result") +
    theme_exc
  enough_for_boxplot <- plot_data %>% count(connection_observed) %>% filter(n >= 2) %>% nrow() > 0
  if (enough_for_boxplot) {
    fig2 <- fig2 + geom_boxplot(aes(group = connection_observed), width = 0.45, alpha = 0.25, outlier.shape = NA)
  }
  figure_files <- c(figure_files, save_plot_both(fig2, "excavation_bruvo_distance_by_connection"))
  
  fig3 <- ggplot(plot_data, aes(x = distance_m, y = bruvo_distance, color = connection_observed, shape = connection_observed)) +
    geom_point(size = 2.7, alpha = 0.85) +
    labs(x = "Physical distance between stems (m)", y = "Bruvo distance", color = "Connection", shape = "Connection") +
    theme_exc
  figure_files <- c(figure_files, save_plot_both(fig3, "excavation_bruvo_distance_by_physical_distance"))
  
  fig4 <- ggplot(plot_data, aes(x = clone_call, y = distance_m, color = connection_observed)) +
    geom_point(position = position_jitter(width = 0.08, height = 0), size = 2.7, alpha = 0.85) +
    labs(x = "Clone call", y = "Physical distance between stems (m)", color = "Connection") +
    theme_exc
  figure_files <- c(figure_files, save_plot_both(fig4, "excavation_clone_call_by_distance"))
  
  cat("\n[excavation_validation] Figure files written:\n")
  cat(paste0("- ", figure_files, collapse = "\n"), "\n", sep = "")
  
  invisible(figure_files)
}

write_summaries <- function(results) {
  connection_summary <- results %>%
    count(connection_observed, clone_call, name = "n_pairs") %>%
    arrange(connection_observed, clone_call)
  
  distance_summary <- bind_rows(
    results %>% group_by(clone_call) %>% summarise(summary_type = "mean_distance_m_by_clone_call", group = first(clone_call), mean_distance_m = mean(distance_m, na.rm = TRUE), n_pairs = n(), .groups = "drop") %>% select(summary_type, group, mean_distance_m, n_pairs),
    results %>% group_by(connection_observed) %>% summarise(summary_type = "mean_distance_m_by_connection_observed", group = first(connection_observed), mean_distance_m = mean(distance_m, na.rm = TRUE), n_pairs = n(), .groups = "drop") %>% select(summary_type, group, mean_distance_m, n_pairs)
  )
  
  readr::write_csv(connection_summary, file.path(OUTPUT_DIR, "excavation_connection_vs_genotype_summary.csv"))
  readr::write_csv(distance_summary, file.path(OUTPUT_DIR, "excavation_distance_summary.csv"))
  
  list(connection_summary = connection_summary, distance_summary = distance_summary)
}

print_console_summary <- function(field, id_check, results) {
  connected <- results %>% filter(connection_observed == "yes")
  unconnected <- results %>% filter(connection_observed == "no")
  
  cat("\n[excavation_validation] Summary\n")
  cat("Field-to-lab mapping source: ", attr(field, "mapping_source"), "\n", sep = "")
  cat("Number of excavation field pairs: ", nrow(field), "\n", sep = "")
  cat("Number of pairs with both lab IDs mapped: ", sum(!is.na(field$lab_id_1) & !is.na(field$lab_id_2)), "\n", sep = "")
  cat("Number of pairs with both lab IDs found in genotype file: ", sum(id_check$id_check_status == "both_lab_ids_found"), "\n", sep = "")
  cat("Number of connected pairs: ", nrow(connected), "\n", sep = "")
  cat("Number of connected pairs confirmed as clones: ", sum(connected$clone_call == "clone", na.rm = TRUE), "\n", sep = "")
  cat("Number of connected pairs genetically distinct: ", sum(connected$clone_call == "not_clone", na.rm = TRUE), "\n", sep = "")
  cat("Number of unconnected pairs confirmed as clones: ", sum(unconnected$clone_call == "clone", na.rm = TRUE), "\n", sep = "")
  cat("Number of unconnected pairs genetically distinct: ", sum(unconnected$clone_call == "not_clone", na.rm = TRUE), "\n", sep = "")
  cat("Mean Bruvo distance overall: ", round(mean(results$bruvo_distance, na.rm = TRUE), 4), "\n", sep = "")
  cat("Mean Bruvo distance for connected pairs: ", round(mean(connected$bruvo_distance, na.rm = TRUE), 4), "\n", sep = "")
  cat("Mean Bruvo distance for unconnected pairs: ", round(mean(unconnected$bruvo_distance, na.rm = TRUE), 4), "\n", sep = "")
  cat("Mean physical distance overall: ", round(mean(results$distance_m, na.rm = TRUE), 3), "\n", sep = "")
  cat("Mean physical distance for clone pairs: ", round(mean(results$distance_m[results$clone_call == "clone"], na.rm = TRUE), 3), "\n", sep = "")
  cat("Mean physical distance for non-clone pairs: ", round(mean(results$distance_m[results$clone_call == "not_clone"], na.rm = TRUE), 3), "\n", sep = "")
}

main <- function() {
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
  
  geno_info <- read_genotype_data()
  sample_missingness <- calculate_sample_missingness(geno_info$data, geno_info$locus_pairs)
  write_available_lab_ids(sample_missingness)
  
  field <- read_and_clean_field_data()
  genotype_ids <- geno_info$data$sample_id_for_analysis
  
  id_check <- make_id_check(field, genotype_ids)
  readr::write_csv(id_check, file.path(OUTPUT_DIR, "excavation_id_check.csv"))
  readr::write_csv(sample_missingness, file.path(OUTPUT_DIR, "excavation_sample_missingness.csv"))
  
  validate_allele_data_present(geno_info$data, geno_info$locus_pairs)
  
  gi <- build_genind(geno_info$data, geno_info$locus_pairs)
  mlg_labels <- get_mlg_labels(gi)
  bruvo_matrix <- calculate_bruvo_distance_matrix(gi)
  
  results <- pairwise_results(field, sample_missingness, mlg_labels, bruvo_matrix)
  summaries <- write_summaries(results)
  
  readr::write_csv(field, file.path(OUTPUT_DIR, "excavation_cleaned.csv"))
  readr::write_csv(results, file.path(OUTPUT_DIR, "excavation_pairwise_genetic_results.csv"))
  
  write_word_table(results)
  make_figures(results)
  print_console_summary(field, id_check, results)
  
  invisible(list(field = field, id_check = id_check, sample_missingness = sample_missingness, results = results, summaries = summaries))
}

main()
# end of file