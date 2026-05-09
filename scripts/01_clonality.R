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
# - outputs/tables/clonality_summary_table.csv
# - outputs/tables/clonality_summary_table.xlsx
# - outputs/tables/clonality_summary_table.docx
# - outputs/tables/clonality_summary_table_numeric.csv
############################################################

suppressPackageStartupMessages({
  library(poppr)
  library(adegenet)
  library(dplyr)
  library(ggplot2)
  library(readxl)
  library(openxlsx)
})

source("scripts/_load_objects.R")

message("[01_clonality] Calculating clonality summaries for both MLG and Bruvo-based MLL on gi...")

DEFAULT_BRUVO_MLL_THRESHOLD <- 0.09
DEFAULT_BRUVO_ALGORITHM <- "farthest_neighbor"

SITE_LOOKUP_REQUIRED_COLUMNS <- c("Site", "Site_label", "Region", "Site_order", "Numéro_Population")
SITE_LOOKUP_SHEET <- "site_lookup"
GENETIC_SHEET_CANDIDATES <- c("genetique", "génétique", "genetic", "genetics")
ARBRE_SHEET_CANDIDATES <- c("arbre", "trees")

normalize_lookup_key <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("\uFEFF", "", x, fixed = TRUE)
  x <- gsub("[[:cntrl:]]", "", x)
  x[is.na(x)] <- ""
  x
}

find_site_lookup_workbook <- function(raw_dir = file.path(PROJECT_ROOT, "data", "raw"), sheet = SITE_LOOKUP_SHEET) {
  if (!dir.exists(raw_dir)) {
    message("[01_clonality] site_lookup unavailable because data/raw does not exist: ", raw_dir)
    return(NULL)
  }
  excel_files <- list.files(raw_dir, pattern = "\\.(xlsx|xls)$", full.names = TRUE, ignore.case = TRUE)
  if (length(excel_files) == 0) {
    message("[01_clonality] site_lookup unavailable: no Excel workbook was found in data/raw.")
    return(NULL)
  }
  
  has_lookup <- vapply(
    excel_files,
    function(path) {
      sheet_names <- tryCatch(readxl::excel_sheets(path), error = function(e) character(0))
      sheet %in% sheet_names
    },
    logical(1)
  )
  lookup_files <- excel_files[has_lookup]
  
  if (length(lookup_files) == 0) {
    message(
      "[01_clonality] site_lookup unavailable: none of the Excel workbooks in data/raw contains a '",
      sheet,
      "' sheet. Checked: ",
      paste(basename(excel_files), collapse = ", ")
    )
    return(NULL)
  }
  if (length(lookup_files) > 1) {
    message(
      "[01_clonality] Multiple workbooks contain site_lookup; using the first one returned by list.files(): ",
      basename(lookup_files[1]),
      ". Other matches: ",
      paste(basename(lookup_files[-1]), collapse = ", ")
    )
  } else {
    message("[01_clonality] Reading site_lookup from: ", basename(lookup_files[1]))
  }
  
  lookup_files[1]
}

load_site_lookup <- function() {
  workbook <- find_site_lookup_workbook()
  if (is.null(workbook)) return(NULL)
  lookup <- suppressMessages(readxl::read_excel(workbook, sheet = SITE_LOOKUP_SHEET)) %>%
    as.data.frame(stringsAsFactors = FALSE)
  
  missing_cols <- setdiff(SITE_LOOKUP_REQUIRED_COLUMNS, names(lookup))
  if (length(missing_cols) > 0) {
    stop(
      "[01_clonality] site_lookup is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      "\nRequired columns are: ",
      paste(SITE_LOOKUP_REQUIRED_COLUMNS, collapse = ", ")
    )
  }
  
  lookup <- lookup %>%
    mutate(
      Site = normalize_lookup_key(Site),
      Site_label = normalize_lookup_key(Site_label),
      Region = normalize_lookup_key(Region),
      Site_order = suppressWarnings(as.numeric(Site_order)),
      Numéro_Population = normalize_lookup_key(Numéro_Population)
    )
  
  if (any(is.na(lookup$Site_order))) {
    stop("[01_clonality] site_lookup$Site_order must be numeric and non-missing for every row.")
  }
  if (any(!nzchar(lookup$Site))) {
    stop("[01_clonality] site_lookup$Site contains blank or missing site codes.")
  }
  if (any(!nzchar(lookup$Site_label))) {
    stop("[01_clonality] site_lookup$Site_label contains blank or missing display labels.")
  }
  if (anyDuplicated(lookup$Site)) {
    dup_sites <- unique(lookup$Site[duplicated(lookup$Site)])
    stop("[01_clonality] site_lookup$Site contains duplicate site codes: ", paste(dup_sites, collapse = ", "))
  }
  if (anyDuplicated(lookup$Site_label)) {
    dup_labels <- unique(lookup$Site_label[duplicated(lookup$Site_label)])
    stop("[01_clonality] site_lookup$Site_label contains duplicate display labels: ", paste(dup_labels, collapse = ", "))
  }
  if (anyDuplicated(lookup$Site_order)) {
    dup_orders <- unique(lookup$Site_order[duplicated(lookup$Site_order)])
    stop("[01_clonality] site_lookup$Site_order contains duplicate ordering values: ", paste(dup_orders, collapse = ", "))
  }
  pop_nonblank <- nzchar(lookup$Numéro_Population)
  if (anyDuplicated(lookup$Numéro_Population[pop_nonblank])) {
    dup_pops <- unique(lookup$Numéro_Population[pop_nonblank][duplicated(lookup$Numéro_Population[pop_nonblank])])
    stop("[01_clonality] site_lookup$Numéro_Population contains duplicate non-blank population numbers: ", paste(dup_pops, collapse = ", "))
  }
  
  lookup %>%
    arrange(Site_order, Site_label)
}

validate_clonality_plot_site_labels <- function(summary_tbl, lookup) {
  ordered_lookup <- lookup %>%
    arrange(Site_order, Site_label)
  expected_labels <- as.character(ordered_lookup$Site_label)
  expected_orders <- ordered_lookup$Site_order
  plot_tbl <- summary_tbl %>%
    arrange(Site_order, Site_label)
  plot_labels <- as.character(plot_tbl$Site_label)
  plot_orders <- plot_tbl$Site_order
  
  if (any(is.na(plot_labels) | !nzchar(plot_labels))) {
    stop("[01_clonality] Clonality plot Site_label values contain missing or blank display labels.")
  }
  if (anyDuplicated(plot_labels)) {
    dup_labels <- unique(plot_labels[duplicated(plot_labels)])
    stop("[01_clonality] Clonality plot Site_label values are duplicated: ", paste(dup_labels, collapse = ", "))
  }
  if (any(is.na(plot_orders))) {
    stop("[01_clonality] Clonality plot Site_order values contain missing ordering values.")
  }
  if (anyDuplicated(plot_orders)) {
    dup_orders <- unique(plot_orders[duplicated(plot_orders)])
    stop("[01_clonality] Clonality plot Site_order values are duplicated: ", paste(dup_orders, collapse = ", "))
  }
  if (is.unsorted(plot_orders, strictly = TRUE)) {
    stop("[01_clonality] Clonality plot Site_label values are not ordered by strictly increasing Site_order.")
  }
  
  missing_labels <- setdiff(expected_labels, plot_labels)
  extra_labels <- setdiff(plot_labels, expected_labels)
  if (length(missing_labels) > 0) {
    stop("[01_clonality] Clonality plot is missing Site_label values from site_lookup: ", paste(missing_labels, collapse = ", "))
  }
  if (length(extra_labels) > 0) {
    stop("[01_clonality] Clonality plot contains Site_label values not found in site_lookup: ", paste(extra_labels, collapse = ", "))
  }
  if (!identical(plot_labels, expected_labels) || !identical(as.numeric(plot_orders), as.numeric(expected_orders))) {
    stop(
      "[01_clonality] Clonality plot Site_label order does not match site_lookup$Site_order. Expected: ",
      paste(expected_labels, collapse = ", "),
      "; got: ",
      paste(plot_labels, collapse = ", ")
    )
  }
  
  message("[01_clonality] Final clonality plot site order: ", paste(plot_labels, collapse = ", "))
  invisible(plot_labels)
}

resolve_population_col <- function(df) {
  resolve_col_ci(df, c("Numéro_Population", "Numero_Population", "numero_population", "population_number", "Population_Number"))
}

join_site_lookup <- function(df, lookup, context = "[01_clonality] clonality data") {
  out <- df
  out <- out %>%
    select(-any_of(c("Site_label", "Region", "Site_order", "Numéro_Population_lookup")))
  site_col_local <- resolve_col_ci(out, c("Site", "site"))
  pop_col_local <- resolve_population_col(out)
  
  if (is.null(lookup)) {
    if (is.na(site_col_local)) {
      stop(context, " must contain a Site column when site_lookup is unavailable.")
    }
    out$Site <- normalize_lookup_key(out[[site_col_local]])
    out$Site_label <- out$Site
    out$Region <- NA_character_
    site_order_levels <- unique(out$Site[nzchar(out$Site)])
    out$Site_order <- match(out$Site, site_order_levels)
    return(out)
  }
  
  if (!is.na(site_col_local)) {
    out$.__site_lookup_key <- normalize_lookup_key(out[[site_col_local]])
    missing_keys <- sort(unique(out$.__site_lookup_key[!out$.__site_lookup_key %in% lookup$Site]))
    missing_keys <- missing_keys[nzchar(missing_keys)]
    if (length(missing_keys) > 0) {
      stop(
        context,
        " contains Site values missing from site_lookup$Site: ",
        paste(missing_keys, collapse = ", ")
      )
    }
    
    out <- out %>%
      left_join(
        lookup %>%
          select(Site, Site_label, Region, Site_order, Numéro_Population) %>%
          rename(.__site_lookup_key = Site, Numéro_Population_lookup = Numéro_Population),
        by = ".__site_lookup_key"
      ) %>%
      select(-.__site_lookup_key)
  } else if (!is.na(pop_col_local)) {
    out$.__population_lookup_key <- normalize_lookup_key(out[[pop_col_local]])
    lookup_pop <- lookup %>%
      filter(nzchar(Numéro_Population))
    missing_keys <- sort(unique(out$.__population_lookup_key[!out$.__population_lookup_key %in% lookup_pop$Numéro_Population]))
    missing_keys <- missing_keys[nzchar(missing_keys)]
    if (length(missing_keys) > 0) {
      stop(
        context,
        " contains Numéro_Population values missing from site_lookup$Numéro_Population: ",
        paste(missing_keys, collapse = ", ")
      )
    }
    
    out <- out %>%
      left_join(
        lookup_pop %>%
          select(Site, Site_label, Region, Site_order, Numéro_Population) %>%
          rename(.__population_lookup_key = Numéro_Population),
        by = ".__population_lookup_key"
      ) %>%
      select(-.__population_lookup_key)
  } else {
    stop(context, " must contain either a Site column or a Numéro_Population column for site_lookup joining.")
  }
  
  if (any(is.na(out$Site_label) | !nzchar(as.character(out$Site_label)))) {
    stop(context, " could not be fully annotated with Site_label from site_lookup.")
  }
  if (any(is.na(out$Site_order))) {
    stop(context, " could not be fully annotated with Site_order from site_lookup.")
  }
  
  label_levels <- lookup %>%
    arrange(Site_order) %>%
    pull(Site_label) %>%
    unique()
  out$Site_label <- factor(out$Site_label, levels = label_levels, ordered = TRUE)
  out
}

clonality_plot_theme <- function(base_size = 22) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(size = 26, face = "bold"),
      plot.subtitle = element_text(size = 22),
      axis.title = element_text(size = 26, face = "bold"),
      axis.text = element_text(size = 22),
      axis.text.x = element_text(size = 22, angle = 45, hjust = 1, vjust = 1),
      axis.text.y = element_text(size = 22),
      axis.line = element_line(linewidth = 0.4, color = "black"),
      axis.ticks = element_line(linewidth = 1.0),
      axis.ticks.length = grid::unit(0.28, "cm"),
      panel.border = element_blank(),
      plot.background = element_blank(),
      panel.background = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      legend.title = element_text(size = 22),
      legend.text = element_text(size = 22)
    )
}


normalize_ascii_key <- function(x) {
  x <- iconv(as.character(x), from = "", to = "ASCII//TRANSLIT")
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x[is.na(x)] <- ""
  x
}

compact_match_key <- function(...) {
  parts <- list(...)
  raw <- do.call(paste, c(lapply(parts, as.character), sep = "_"))
  key <- toupper(gsub("[^A-Za-z0-9]+", "", raw))
  key[is.na(key)] <- ""
  key
}

coerce_numeric_clonality <- function(x) {
  if (is.numeric(x)) return(x)
  x_chr <- as.character(x)
  x_chr <- gsub(",", ".", x_chr, fixed = TRUE)
  x_chr <- gsub("[^0-9.\\-]+", "", x_chr)
  suppressWarnings(as.numeric(x_chr))
}

normalize_tree_number <- function(x) {
  x_chr <- trimws(as.character(x))
  x_chr[x_chr %in% c("", "NA", "NaN", "NULL")] <- NA_character_
  x_num <- suppressWarnings(as.numeric(gsub(",", ".", x_chr, fixed = TRUE)))
  x_chr <- ifelse(!is.na(x_num) & abs(x_num - round(x_num)) < .Machine$double.eps^0.5, as.character(as.integer(round(x_num))), x_chr)
  x_chr
}

find_preferred_column_ci_norm <- function(df, candidates) {
  cn <- names(df)
  cn_norm <- normalize_ascii_key(cn)
  candidate_norm <- normalize_ascii_key(candidates)
  idx <- match(candidate_norm, cn_norm, nomatch = 0)
  idx <- idx[idx > 0]
  if (length(idx) == 0) NA_character_ else cn[idx[1]]
}

pick_pattern_column_base <- function(df, patterns, label, required = FALSE, exclude = character()) {
  cn <- names(df)
  cn_norm <- normalize_ascii_key(cn)
  exclude_norm <- normalize_ascii_key(exclude)
  idx <- which(grepl(paste(patterns, collapse = "|"), cn_norm) & !(cn_norm %in% exclude_norm))
  if (length(idx) == 0) {
    if (required) warning("[01_clonality] Could not detect required column for ", label, call. = FALSE)
    return(NA_character_)
  }
  if (length(idx) > 1) {
    message("[01_clonality] Multiple matches for ", label, ": ", paste(cn[idx], collapse = ", "), ". Using: ", cn[idx[1]])
  }
  cn[idx[1]]
}

find_dbh_col <- function(df) {
  find_preferred_column_ci_norm(df, c("dbh", "DBH", "dbh_cm", "DBH_cm", "dhp", "DHP", "dhp_cm", "DHP_cm"))
}

find_individual_id_col <- function(df) {
  find_preferred_column_ci_norm(
    df,
    c(
      DF_IDS_ID_CHOICES,
      "individual_id", "Individual_ID", "individu", "Individu", "id_individu", "ID_Individu",
      "sample_id", "Sample_ID", "sample_name", "Sample_Name", "echantillon", "échantillon",
      "id_echantillon", "id_échantillon", "identifiant", "genetic_id", "Genetic_ID"
    )
  )
}

find_tree_number_col <- function(df) {
  find_preferred_column_ci_norm(
    df,
    c(
      "tree", "Tree", "tree_id", "Tree_ID", "tree_number", "Tree_Number", "tree_no", "Tree_No",
      "arbre", "Arbre", "no_arbre", "No_Arbre", "numero_arbre", "Numéro_Arbre",
      "numero_tree", "num_tree", "no_tree", "id_arbre", "id_tige", "stem", "stem_id",
      "numero_individu", "num_individu", "no_individu", "id_individu",
      "sample_number", "sample_no", "numero", "num", "no"
    )
  )
}

strip_site_prefix_from_individual_id <- function(individual_id, site) {
  out <- trimws(as.character(individual_id))
  site_chr <- trimws(as.character(site))
  if (length(site_chr) == 1L && length(out) != 1L) {
    site_chr <- rep(site_chr, length(out))
  }
  for (i in seq_along(out)) {
    if (is.na(out[i]) || is.na(site_chr[i]) || !nzchar(site_chr[i])) next
    site_nchar <- nchar(site_chr[i], type = "chars", allowNA = FALSE, keepNA = FALSE)
    if (nchar(out[i], type = "chars", allowNA = FALSE, keepNA = FALSE) < site_nchar) next
    possible_prefix <- substr(out[i], 1L, site_nchar)
    if (!identical(toupper(possible_prefix), toupper(site_chr[i]))) next
    out[i] <- substr(out[i], site_nchar + 1L, nchar(out[i], type = "chars", allowNA = FALSE, keepNA = FALSE))
    while (nzchar(out[i]) && substr(out[i], 1L, 1L) %in% c("-", "_", ".", " ")) {
      out[i] <- substr(out[i], 2L, nchar(out[i], type = "chars", allowNA = FALSE, keepNA = FALSE))
    }
  }
  out
}

derive_tree_number_from_individual_id <- function(individual_id, site = NULL) {
  x <- trimws(as.character(individual_id))
  if (!is.null(site)) {
    x <- strip_site_prefix_from_individual_id(x, site)
  }
  out <- sub(".*?([0-9]+)[^0-9]*$", "\\1", x)
  unchanged <- out == x
  unchanged[is.na(unchanged)] <- TRUE
  out[unchanged | is.na(x) | !grepl("[0-9]", x)] <- NA_character_
  normalize_tree_number(out)
}

find_raw_workbook_for_clonality <- function(raw_dir = file.path(PROJECT_ROOT, "data", "raw")) {
  if (!dir.exists(raw_dir)) stop("[01_clonality] Directory not found while searching for DBH workbook: ", raw_dir, call. = FALSE)
  files <- list.files(raw_dir, pattern = "\\.(xlsx|xls)$", full.names = TRUE, ignore.case = TRUE)
  if (length(files) == 0) stop("[01_clonality] No Excel file found in data/raw for DBH matching.", call. = FALSE)
  preferred_idx <- which(grepl("donnee|modifie|west|summer|copie", normalize_ascii_key(basename(files))))
  chosen <- if (length(preferred_idx) > 0) files[preferred_idx[1]] else files[1]
  message("[01_clonality] Excel files found in data/raw for DBH matching:")
  message(paste0(" - ", basename(files), collapse = "\n"))
  message("[01_clonality] Selected workbook for DBH matching: ", basename(chosen))
  chosen
}

match_sheet_clonality <- function(sheets, target) {
  idx <- which(normalize_ascii_key(sheets) == normalize_ascii_key(target))
  if (length(idx) == 0) NA_character_ else sheets[idx[1]]
}

make_retained_identity_table <- function(clonality_assignments, df_ids_tbl) {
  df_ids_cols_local <- resolve_df_ids_columns(df_ids_tbl, context = "[01_clonality]", require = TRUE)
  id_key_local <- normalize_id(df_ids_tbl[[df_ids_cols_local$id_col]])
  tree_col <- find_tree_number_col(df_ids_tbl)
  tree_map <- if (!is.na(tree_col)) setNames(as.character(df_ids_tbl[[tree_col]]), id_key_local) else NULL
  retained_key <- normalize_id(clonality_assignments$Individual)
  retained_tree <- if (!is.null(tree_map)) normalize_tree_number(tree_map[retained_key]) else rep(NA_character_, nrow(clonality_assignments))
  missing_tree <- is.na(retained_tree) | !nzchar(retained_tree)
  if (any(missing_tree)) {
    retained_tree[missing_tree] <- derive_tree_number_from_individual_id(
      clonality_assignments$Individual[missing_tree],
      clonality_assignments$Site[missing_tree]
    )
  }
  data.frame(
    Individual = as.character(clonality_assignments$Individual),
    Site = as.character(clonality_assignments$Site),
    Site_label = as.character(clonality_assignments$Site_label),
    Tree_Number = retained_tree,
    stringsAsFactors = FALSE
  )
}

get_metadata_dbh_for_retained <- function(retained_tbl, metadata_tbl, source_name) {
  dbh_col <- find_dbh_col(metadata_tbl)
  if (is.na(dbh_col)) return(NULL)
  cols <- resolve_df_ids_columns(metadata_tbl, context = paste0("[01_clonality] ", source_name), require = FALSE)
  if (is.na(cols$id_col)) {
    message("[01_clonality] ", source_name, " contains DBH column '", dbh_col, "' but no recognized individual ID column; skipping this DBH source.")
    return(NULL)
  }
  key <- normalize_id(metadata_tbl[[cols$id_col]])
  if (anyDuplicated(key)) {
    dup_ids <- unique(as.character(metadata_tbl[[cols$id_col]])[duplicated(key)])
    stop("[01_clonality] ", source_name, " has duplicate individual IDs while matching DBH. Examples: ", paste(head(dup_ids, 10), collapse = ", "), call. = FALSE)
  }
  idx <- match(normalize_id(retained_tbl$Individual), key)
  if (any(is.na(idx))) {
    missing_ids <- retained_tbl$Individual[is.na(idx)]
    stop("[01_clonality] ", source_name, " has a DBH column but could not be aligned to all retained individuals. Missing examples: ", paste(head(missing_ids, 10), collapse = ", "), call. = FALSE)
  }
  list(
    dbh = coerce_numeric_clonality(metadata_tbl[[dbh_col]][idx]),
    matched = rep(TRUE, nrow(retained_tbl)),
    source = paste0(source_name, "$", dbh_col)
  )
}