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
  usable_key <- nzchar(key)
  if (anyDuplicated(key[usable_key])) {
    dup_ids <- unique(as.character(metadata_tbl[[cols$id_col]])[usable_key][duplicated(key[usable_key])])
    stop("[01_clonality] ", source_name, " has duplicate individual IDs while matching DBH. Examples: ", paste(head(dup_ids, 10), collapse = ", "), call. = FALSE)
  }
  idx <- match(normalize_id(retained_tbl$Individual), key)
  dbh <- rep(NA_real_, nrow(retained_tbl))
  matched <- !is.na(idx)
  dbh[matched] <- coerce_numeric_clonality(metadata_tbl[[dbh_col]][idx[matched]])
  list(
    dbh = dbh,
    matched = matched,
    source = paste0(source_name, "$", dbh_col)
  )
}

make_raw_dbh_source_table <- function(raw_df, sheet_name) {
  dbh_col <- find_dbh_col(raw_df)
  if (is.na(dbh_col)) return(NULL)
  id_col <- find_individual_id_col(raw_df)
  site_col <- resolve_col_ci(raw_df, c("Site", "site", "pop", "population"))
  if (is.na(site_col)) {
    site_col <- pick_pattern_column_base(raw_df, c("^site$", "site_?id", "id_?site", "placette", "parcelle", "plot", "station", "localite", "location"), paste0(sheet_name, " site"))
  }
  tree_col <- find_tree_number_col(raw_df)
  source_id <- if (!is.na(id_col)) as.character(raw_df[[id_col]]) else rep(NA_character_, nrow(raw_df))
  source_site <- if (!is.na(site_col)) normalize_lookup_key(raw_df[[site_col]]) else rep(NA_character_, nrow(raw_df))
  source_tree <- if (!is.na(tree_col)) normalize_tree_number(raw_df[[tree_col]]) else rep(NA_character_, nrow(raw_df))
  missing_tree <- is.na(source_tree) | !nzchar(source_tree)
  if (any(missing_tree) && !is.na(id_col)) {
    source_tree[missing_tree] <- derive_tree_number_from_individual_id(source_id[missing_tree], source_site[missing_tree])
  }
  data.frame(
    Source_ID = source_id,
    Site = source_site,
    Tree_Number = source_tree,
    DBH_cm = coerce_numeric_clonality(raw_df[[dbh_col]]),
    stringsAsFactors = FALSE
  )
}

match_retained_dbh_from_source <- function(retained_tbl, source_tbl, source_name) {
  attempts <- list()
  if (any(!is.na(source_tbl$Source_ID) & nzchar(source_tbl$Source_ID))) {
    attempts[["individual ID"]] <- list(
      retained_key = compact_match_key(retained_tbl$Individual),
      source_key = compact_match_key(source_tbl$Source_ID)
    )
  }
  if (any(!is.na(retained_tbl$Site) & nzchar(retained_tbl$Site)) && any(!is.na(source_tbl$Site) & nzchar(source_tbl$Site)) && any(!is.na(retained_tbl$Tree_Number) & nzchar(retained_tbl$Tree_Number)) && any(!is.na(source_tbl$Tree_Number) & nzchar(source_tbl$Tree_Number))) {
    attempts[["site + tree number"]] <- list(
      retained_key = compact_match_key(retained_tbl$Site, retained_tbl$Tree_Number),
      source_key = compact_match_key(source_tbl$Site, source_tbl$Tree_Number)
    )
  }
  if (length(attempts) == 0) return(NULL)
  best <- NULL
  for (method in names(attempts)) {
    retained_key <- attempts[[method]]$retained_key
    source_key <- attempts[[method]]$source_key
    usable_source <- nzchar(source_key)
    source_key_use <- source_key[usable_source]
    source_tbl_use <- source_tbl[usable_source, , drop = FALSE]
    if (anyDuplicated(source_key_use)) next
    idx <- match(retained_key, source_key_use)
    n_match <- sum(!is.na(idx))
    candidate <- list(idx = idx, source_tbl = source_tbl_use, n_match = n_match, method = method)
    if (is.null(best) || candidate$n_match > best$n_match) best <- candidate
  }
  if (is.null(best) || best$n_match == 0) return(NULL)
  dbh <- rep(NA_real_, nrow(retained_tbl))
  matched <- !is.na(best$idx)
  dbh[matched] <- best$source_tbl$DBH_cm[best$idx[matched]]
  list(dbh = dbh, matched = matched, source = paste0(source_name, " via ", best$method))
}
get_raw_dbh_for_retained <- function(retained_tbl) {
  workbook <- find_raw_workbook_for_clonality()
  sheets <- readxl::excel_sheets(workbook)
  message("[01_clonality] Workbook sheets available for DBH matching: ", paste(sheets, collapse = ", "))
  preferred_sheets <- c(GENETIC_SHEET_CANDIDATES, ARBRE_SHEET_CANDIDATES)
  best <- NULL
  for (target in preferred_sheets) {
    sheet <- match_sheet_clonality(sheets, target)
    if (is.na(sheet)) next
    raw_df <- suppressMessages(readxl::read_excel(workbook, sheet = sheet)) %>% as.data.frame(stringsAsFactors = FALSE)
    source_tbl <- make_raw_dbh_source_table(raw_df, sheet)
    if (is.null(source_tbl)) next
    matched <- match_retained_dbh_from_source(retained_tbl, source_tbl, paste0(basename(workbook), "::", sheet))
    if (!is.null(matched) && normalize_ascii_key(sheet) %in% normalize_ascii_key(GENETIC_SHEET_CANDIDATES)) {
      matched$source <- paste0(matched$source, " [genetic sheet DBH]")
    }
    if (is.null(matched)) next
    if (is.null(best) || sum(matched$matched) > sum(best$matched)) best <- matched
    if (all(matched$matched)) break
  }
  best
}

attach_dbh_to_clonality_assignments <- function(clonality_assignments, df_ids_tbl, meta_tbl = NULL) {
  retained_tbl <- make_retained_identity_table(clonality_assignments, df_ids_tbl)
  message("[01_clonality] Retained sampled stems used in clonality analysis: ", nrow(retained_tbl))
  dbh_match <- get_metadata_dbh_for_retained(retained_tbl, df_ids_tbl, "df_ids")
  if (!is.null(dbh_match) && sum(!is.na(dbh_match$dbh)) == 0) {
    message("[01_clonality] df_ids contains a DBH column, but no numeric DBH values matched retained stems; trying the next DBH source.")
    dbh_match <- NULL
  }
  if (is.null(dbh_match) && !is.null(meta_tbl) && is.data.frame(meta_tbl)) {
    dbh_match <- get_metadata_dbh_for_retained(retained_tbl, meta_tbl, "meta")
    if (!is.null(dbh_match) && sum(!is.na(dbh_match$dbh)) == 0) {
      message("[01_clonality] meta contains a DBH column, but no numeric DBH values matched retained stems; trying the raw Excel workbook.")
      dbh_match <- NULL
    }
  }
  if (is.null(dbh_match)) {
    message("[01_clonality] No usable DBH values detected in genetic metadata; attempting DBH matching from the raw Excel workbook.")
    dbh_match <- get_raw_dbh_for_retained(retained_tbl)
  }
  if (is.null(dbh_match)) {
    stop("[01_clonality] Could not match DBH for retained sampled stems. No usable DBH source was found in df_ids, meta, or the raw Excel genetique/arbre sheets.", call. = FALSE)
  }
  unmatched <- retained_tbl$Individual[!dbh_match$matched]
  n_with_dbh <- sum(!is.na(dbh_match$dbh))
  n_missing_dbh <- sum(is.na(dbh_match$dbh) | !dbh_match$matched)
  message("[01_clonality] DBH source used for retained sampled stems: ", dbh_match$source)
  message("[01_clonality] Retained individuals with DBH values matched: ", n_with_dbh, " / ", nrow(retained_tbl))
  message("[01_clonality] Retained individuals missing DBH values: ", n_missing_dbh)
  if (n_with_dbh == 0) {
    stop("[01_clonality] No DBH values could be matched for the retained sampled stems. Check DBH column names and individual/sample ID matching in df_ids, meta, or the raw Excel genetique/arbre sheets.", call. = FALSE)
  }
  if (n_missing_dbh > 0) {
    warning(
      "[01_clonality] Some retained sampled stems are missing DBH values; continuing and calculating site means with available DBH values only. Missing count: ",
      n_missing_dbh,
      if (length(unmatched) > 0) paste0(". Unmatched examples: ", paste(head(unmatched, 10), collapse = ", ")) else "",
      call. = FALSE
    )
  }
  clonality_assignments$DBH_cm <- dbh_match$dbh
  clonality_assignments
}

validate_article_clonality_summary_table <- function(summary_tbl, site_summary_tbl) {
  if (anyDuplicated(as.character(summary_tbl$Site))) {
    dup_sites <- unique(as.character(summary_tbl$Site)[duplicated(as.character(summary_tbl$Site))])
    stop("[01_clonality] Final article-ready clonality table has duplicated sites: ", paste(dup_sites, collapse = ", "), call. = FALSE)
  }
  n_unique_sites <- dplyr::n_distinct(as.character(site_summary_tbl$Site_label))
  if (nrow(summary_tbl) > n_unique_sites) {
    stop(
      "[01_clonality] Final article-ready clonality table has more rows (", nrow(summary_tbl),
      ") than the number of unique sites (", n_unique_sites, ").",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

site_lookup <- load_site_lookup()

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

make_clone_group_table <- function(assignments_df, clone_col, individual_col = "Individual", site_col = "Site_label") {
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
      site_value <- as.character(grp[["Site_label"]][i])
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
df_ids <- join_site_lookup(df_ids, site_lookup, context = "[01_clonality] df_ids")

df_ids_cols <- resolve_df_ids_columns(df_ids, context = "[01_clonality]", require = TRUE)
id_col <- df_ids_cols$id_col
site_col <- df_ids_cols$site_col
threshold_col <- resolve_col_ci(df_ids, c("Bruvo_MLL_threshold", "bruvo_mll_threshold"))
algorithm_col <- resolve_col_ci(df_ids, c("Bruvo_algorithm", "bruvo_algorithm"))

id_key <- normalize_id(df_ids[[id_col]])
site_map <- setNames(as.character(df_ids[[site_col]]), id_key)
site_label_map <- setNames(as.character(df_ids[["Site_label"]]), id_key)
site_order_map <- setNames(as.numeric(df_ids[["Site_order"]]), id_key)
region_map <- setNames(as.character(df_ids[["Region"]]), id_key)
mlg_map <- setNames(as.character(df_ids[["MLG"]]), id_key)
mll_map <- setNames(as.character(df_ids[["MLL"]]), id_key)

ind_key <- normalize_id(adegenet::indNames(gi))
site_labels <- site_map[ind_key]
site_display_labels <- site_label_map[ind_key]
site_orders <- site_order_map[ind_key]
regions <- region_map[ind_key]
mlg_labels <- mlg_map[ind_key]
mll_labels <- mll_map[ind_key]

site_available <- !any(is.na(site_display_labels) | !nzchar(site_display_labels) | is.na(site_orders))
if (!site_available) {
  warning("[01_clonality] Site_label or Site_order metadata were unavailable for one or more individuals; quick clone summary will print IDs without site where necessary.")
}
if (any(is.na(mlg_labels))) stop("[01_clonality] Could not map all individuals to MLG.")
if (any(is.na(mll_labels))) stop("[01_clonality] Could not map all individuals to MLL.")

bruvo_threshold <- if (!is.na(threshold_col)) unique(stats::na.omit(df_ids[[threshold_col]])) else numeric(0)
bruvo_algorithm <- if (!is.na(algorithm_col)) unique(stats::na.omit(df_ids[[algorithm_col]])) else character(0)

site_label_levels <- if (!is.null(site_lookup)) {
  site_lookup %>%
    arrange(Site_order) %>%
    pull(Site_label) %>%
    unique()
} else {
  unique(site_display_labels[!is.na(site_display_labels) & nzchar(site_display_labels)])
}

clonality_df <- data.frame(
  Individual = adegenet::indNames(gi),
  Site = ifelse(is.na(site_labels) | !nzchar(site_labels), NA_character_, site_labels),
  Site_label = site_display_labels,
  Region = regions,
  Site_order = site_orders,
  MLG = mlg_labels,
  MLL = mll_labels,
  stringsAsFactors = FALSE
) %>%
  mutate(Site_label = factor(Site_label, levels = site_label_levels, ordered = TRUE))

clonality_df <- attach_dbh_to_clonality_assignments(clonality_df, df_ids, meta_tbl = meta)

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

format_three_decimals <- function(x) {
  out <- ifelse(is.na(x), NA_character_, sprintf("%.3f", round(as.numeric(x), 3)))
  out
}

format_one_decimal <- function(x) {
  out <- ifelse(is.na(x), NA_character_, sprintf("%.1f", round(as.numeric(x), 1)))
  out
}
build_article_clonality_summary_table <- function(site_summary_tbl, use_site_lookup_order = TRUE, formatted = TRUE) {
  out <- site_summary_tbl
  if (isTRUE(use_site_lookup_order) && "Site_order" %in% names(out) && any(!is.na(out$Site_order))) {
    out <- out %>% arrange(Site_order, Site_label)
  } else if ("Latitude" %in% names(out) && any(!is.na(out$Latitude))) {
    out <- out %>% arrange(Latitude, Site_label)
  } else if ("Site_order" %in% names(out) && any(!is.na(out$Site_order))) {
    out <- out %>% arrange(Site_order, Site_label)
  }
  
  numeric_tbl <- out %>%
    transmute(
      Site = as.character(Site_label),
      N = as.integer(N_individuals),
      `Mean DBH (cm)` = round(Mean_DBH_cm, 1),
      MLG = as.integer(N_MLG),
      R_MLG = round(Clonal_Richness_MLG, 3),
      MLL = as.integer(N_MLL),
      R_MLL = round(Clonal_Richness_MLL, 3)
    )
  
  if (!isTRUE(formatted)) return(numeric_tbl)
  
  numeric_tbl %>%
    mutate(
      `Mean DBH (cm)` = format_one_decimal(`Mean DBH (cm)`),
      R_MLG = format_three_decimals(R_MLG),
      R_MLL = format_three_decimals(R_MLL)
    )
}

write_article_clonality_summary_xlsx <- function(summary_tbl, path) {
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "Clonality summary")
  openxlsx::writeData(wb, "Clonality summary", summary_tbl)
  header_style <- openxlsx::createStyle(
    textDecoration = "bold",
    halign = "center",
    valign = "center",
    border = "bottom",
    fgFill = "#D9EAD3"
  )
  body_style <- openxlsx::createStyle(halign = "center", valign = "center")
  openxlsx::addStyle(
    wb,
    "Clonality summary",
    header_style,
    rows = 1,
    cols = seq_len(ncol(summary_tbl)),
    gridExpand = TRUE
  )
  if (nrow(summary_tbl) > 0) {
    openxlsx::addStyle(
      wb,
      "Clonality summary",
      body_style,
      rows = 2:(nrow(summary_tbl) + 1),
      cols = seq_len(ncol(summary_tbl)),
      gridExpand = TRUE
    )
  }
  openxlsx::freezePane(wb, "Clonality summary", firstRow = TRUE)
  openxlsx::setColWidths(wb, "Clonality summary", cols = seq_len(ncol(summary_tbl)), widths = "auto")
  openxlsx::saveWorkbook(wb, path, overwrite = TRUE)
}

word_export_package_status <- function() {
  list(
    officer_available = requireNamespace("officer", quietly = TRUE),
    flextable_available = requireNamespace("flextable", quietly = TRUE)
  )
}

xml_escape_word <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub('"', "&quot;", x, fixed = TRUE)
  x <- gsub("'", "&apos;", x, fixed = TRUE)
  x
}

word_cell_xml <- function(value, bold = FALSE) {
  text <- xml_escape_word(value)
  bold_xml <- if (bold) "<w:rPr><w:b/></w:rPr>" else ""
  paste0(
    "<w:tc>",
    "<w:tcPr><w:tcW w:w=\"2400\" w:type=\"dxa\"/></w:tcPr>",
    "<w:p><w:r>", bold_xml, "<w:t>", text, "</w:t></w:r></w:p>",
    "</w:tc>"
  )
}

word_paragraph_xml <- function(value, bold = FALSE) {
  bold_xml <- if (bold) "<w:rPr><w:b/></w:rPr>" else ""
  paste0("<w:p><w:r>", bold_xml, "<w:t>", xml_escape_word(value), "</w:t></w:r></w:p>")
}

write_basic_docx_table <- function(df, path, title, note = NULL) {
  tmp_dir <- tempfile("clonality_docx_")
  dir.create(file.path(tmp_dir, "_rels"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(tmp_dir, "word", "_rels"), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  
  df <- as.data.frame(lapply(df, as.character), stringsAsFactors = FALSE, check.names = FALSE)
  rows <- apply(df, 1, function(row) {
    paste0("<w:tr>", paste(vapply(row, word_cell_xml, character(1)), collapse = ""), "</w:tr>")
  })
  header <- paste0("<w:tr>", paste(vapply(names(df), word_cell_xml, character(1), bold = TRUE), collapse = ""), "</w:tr>")
  note_xml <- if (!is.null(note) && nzchar(note)) word_paragraph_xml(note) else ""
  
  document_xml <- paste0(
    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
    '<w:document xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main" ',
    'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">',
    '<w:body>',
    word_paragraph_xml(title, bold = TRUE),
    '<w:tbl>',
    '<w:tblPr><w:tblStyle w:val="TableGrid"/><w:tblW w:w="0" w:type="auto"/>',
    '<w:tblBorders>',
    '<w:top w:val="single" w:sz="4" w:space="0" w:color="auto"/>',
    '<w:left w:val="single" w:sz="4" w:space="0" w:color="auto"/>',
    '<w:bottom w:val="single" w:sz="4" w:space="0" w:color="auto"/>',
    '<w:right w:val="single" w:sz="4" w:space="0" w:color="auto"/>',
    '<w:insideH w:val="single" w:sz="4" w:space="0" w:color="auto"/>',
    '<w:insideV w:val="single" w:sz="4" w:space="0" w:color="auto"/>',
    '</w:tblBorders></w:tblPr>',
    header,
    paste(rows, collapse = ""),
    '</w:tbl>',
    note_xml,
    '<w:sectPr><w:pgSz w:w="12240" w:h="15840"/>',
    '<w:pgMar w:top="720" w:right="720" w:bottom="720" w:left="720" w:header="360" w:footer="360" w:gutter="0"/></w:sectPr>',
    '</w:body></w:document>'
  )
  
  writeLines(c(
    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
    '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">',
    '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>',
    '<Default Extension="xml" ContentType="application/xml"/>',
    '<Override PartName="/word/document.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.document.main+xml"/>',
    '</Types>'
  ), file.path(tmp_dir, "[Content_Types].xml"), useBytes = TRUE)
  writeLines(c(
    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
    '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">',
    '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="word/document.xml"/>',
    '</Relationships>'
  ), file.path(tmp_dir, "_rels", ".rels"), useBytes = TRUE)
  writeLines(c(
    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
    '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"></Relationships>'
  ), file.path(tmp_dir, "word", "_rels", "document.xml.rels"), useBytes = TRUE)
  writeLines(document_xml, file.path(tmp_dir, "word", "document.xml"), useBytes = TRUE)
  
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(tmp_dir)
  if (file.exists(path)) unlink(path)
  utils::zip(
    zipfile = path,
    files = list.files(tmp_dir, recursive = TRUE, all.files = TRUE, no.. = TRUE),
    flags = "-q"
  )
  if (!file.exists(path)) stop("Failed to create Word document at: ", path, call. = FALSE)
  invisible(path)
}

build_article_clonality_word_table <- function(summary_tbl) {
  summary_tbl %>%
    transmute(
      Site = as.character(Site),
      N = as.integer(N),
      `Mean DBH (cm)` = as.character(`Mean DBH (cm)`),
      MLG = as.integer(MLG),
      R_MLG = as.character(R_MLG),
      MLL = as.integer(MLL),
      R_MLL = as.character(R_MLL)
    )
}

write_article_clonality_summary_docx <- function(summary_tbl, path) {
  word_tbl <- build_article_clonality_word_table(summary_tbl)
  required_word_cols <- c("Site", "N", "Mean DBH (cm)", "MLG", "R_MLG", "MLL", "R_MLL")
  missing_word_cols <- setdiff(required_word_cols, names(word_tbl))
  if (length(missing_word_cols) > 0) {
    stop("[01_clonality] Word clonality table is missing required column(s): ", paste(missing_word_cols, collapse = ", "), call. = FALSE)
  }
  message("[01_clonality] Word clonality summary table columns: ", paste(names(word_tbl), collapse = ", "))
  table_title <- "Table X. Clonal structure summary by site."
  table_note <- "N = number of individuals analyzed; Mean DBH (cm) = mean diameter at breast height of the retained sampled stems from the genetic/clonality dataset; MLG = number of multilocus genotypes; MLL = number of multilocus lineages after Bruvo-distance clustering; R_MLG and R_MLL represent clonal richness, calculated as (G − 1)/(N − 1)."
  pkg_status <- word_export_package_status()
  
  docx_status <- tryCatch(
    {
      if (pkg_status$officer_available && pkg_status$flextable_available) {
        ft <- flextable::flextable(word_tbl) |>
          flextable::theme_booktabs() |>
          flextable::autofit() |>
          flextable::align(align = "center", part = "all") |>
          flextable::align(j = "Site", align = "left", part = "body") |>
          flextable::bold(part = "header")
        
        doc <- officer::read_docx()
        doc <- officer::body_add_par(doc, table_title, style = "heading 1")
        doc <- officer::body_add_flextable(doc, ft)
        doc <- officer::body_add_par(doc, table_note, style = "Normal")
        print(doc, target = path)
        "created with flextable"
      } else if (pkg_status$officer_available) {
        message(
          "[01_clonality] flextable is not available/loadable, so the script is creating the Word table with officer instead. ",
          "If you want flextable styling, restart R and run install.packages(\"flextable\") once before sourcing this script."
        )
        doc <- officer::read_docx()
        doc <- officer::body_add_par(doc, table_title, style = "heading 1")
        doc <- officer::body_add_table(doc, value = as.data.frame(word_tbl), style = "table_template")
        doc <- officer::body_add_par(doc, table_note, style = "Normal")
        print(doc, target = path)
        "created with officer fallback"
      } else {
        message(
          "[01_clonality] Neither flextable nor officer is loadable, so the script is creating a basic Word-compatible .docx file without extra packages. ",
          "For package-based Word export, run install.packages(c(\"officer\", \"flextable\"))."
        )
        write_basic_docx_table(word_tbl, path, title = table_title, note = table_note)
        "created with built-in docx fallback"
      }
    },
    error = function(e) {
      message(
        "[01_clonality] Package-based Word export failed: ",
        conditionMessage(e),
        ". Creating a basic Word-compatible .docx file without extra packages instead."
      )
      write_basic_docx_table(word_tbl, path, title = table_title, note = table_note)
      "created with built-in docx fallback after package export failed"
    }
  )
  
  if (!file.exists(path)) {
    stop(
      "[01_clonality] Word clonality summary table was not created: ",
      path,
      ". If package-based export is needed, install packages with install.packages(c(\"officer\", \"flextable\")).",
      call. = FALSE
    )
  }
  
  message("[01_clonality] Saved Word clonality summary table: ", path, " (", docx_status, ")")
  invisible(path)
}

print_article_clonality_summary_table <- function(summary_tbl) {
  cat("\n")
  print_separator("=", 72)
  cat("ARTICLE/THESIS-READY CLONALITY SUMMARY TABLE\n")
  print_separator("=", 72)
  utils::write.table(summary_tbl, row.names = FALSE, quote = FALSE, sep = "\t", na = "NA")
  cat("\nCopy-paste table:\n")
  utils::write.table(summary_tbl, row.names = FALSE, quote = FALSE, sep = "\t", na = "NA")
  print_separator("=", 72)
  invisible(summary_tbl)
}

verify_article_clonality_summary_exports <- function(summary_tbl, csv_path, xlsx_path) {
  csv_tbl <- utils::read.csv(csv_path, stringsAsFactors = FALSE, check.names = FALSE, colClasses = "character")
  xlsx_tbl <- openxlsx::read.xlsx(xlsx_path, colNames = TRUE, detectDates = FALSE)
  xlsx_tbl <- as.data.frame(lapply(xlsx_tbl, as.character), stringsAsFactors = FALSE, check.names = FALSE)
  names(xlsx_tbl) <- names(summary_tbl)
  
  expected_tbl <- as.data.frame(lapply(summary_tbl, as.character), stringsAsFactors = FALSE, check.names = FALSE)
  names(expected_tbl) <- names(summary_tbl)
  
  if (!identical(csv_tbl, expected_tbl)) {
    stop("[01_clonality] CSV clonality summary table does not match the printed article-ready table.")
  }
  if (!identical(xlsx_tbl, expected_tbl)) {
    stop("[01_clonality] XLSX clonality summary table does not match the printed article-ready table.")
  }
  
  message("[01_clonality] Verified CSV and XLSX clonality summary tables match the printed article-ready table.")
  invisible(TRUE)
}

overall <- clonality_df %>%
  summarise(
    N_individuals = dplyr::n(),
    Mean_DBH_cm = ifelse(all(is.na(DBH_cm)), NA_real_, mean(DBH_cm, na.rm = TRUE)),
    N_MLG = dplyr::n_distinct(MLG),
    N_MLL = dplyr::n_distinct(MLL)
  ) %>%
  add_clonality_metrics() %>%
  mutate(
    Level = "overall",
    Site = "ALL",
    Site_label = "ALL",
    Region = "ALL",
    Site_order = 0
  ) %>%
  select(Level, Site_label, Site_order, Site, Region, everything())

by_site <- clonality_df %>%
  group_by(Site, Site_label, Region, Site_order) %>%
  summarise(
    N_individuals = dplyr::n(),
    Mean_DBH_cm = ifelse(all(is.na(DBH_cm)), NA_real_, mean(DBH_cm, na.rm = TRUE)),
    N_MLG = dplyr::n_distinct(MLG),
    N_MLL = dplyr::n_distinct(MLL),
    .groups = "drop"
  ) %>%
  add_clonality_metrics() %>%
  mutate(Level = "site") %>%
  arrange(Site_order, Site_label) %>%
  select(Level, Site_label, Site_order, Site, Region, everything())

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

dir.create(TABLES_DIR, recursive = TRUE, showWarnings = FALSE)

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
# Clonality percentage barplots by site (FR + EN) for presentation
# -----------------------------------------------------------------------------
# Objective:
# - x axis: Site_label display labels from site_lookup, ordered by Site_order
# - y axis: Clonality (%) based on MLL clonality proportion
# - robust behavior even when some sites have variable sample sizes

find_latitude_col <- function(df) {
  find_synonym_col(
    df,
    c("Latitude", "latitude", "lat", "LAT", "Lat", "y", "Y")
  )
}

build_site_clonality_summary <- function(assignments_df, df_ids_tbl = NULL) {
  summary_tbl <- assignments_df %>%
    group_by(Site, Site_label, Region, Site_order) %>%
    summarise(
      N_individuals = dplyr::n(),
      Mean_DBH_cm = ifelse(all(is.na(DBH_cm)), NA_real_, mean(DBH_cm, na.rm = TRUE)),
      N_MLG = dplyr::n_distinct(MLG, na.rm = TRUE),
      N_MLL = dplyr::n_distinct(MLL, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    add_clonality_metrics() %>%
    mutate(
      Clonality_MLL = 1 - (N_MLL / N_individuals),
      Clonality_MLL_percent = Clonality_MLL * 100
    ) %>%
    arrange(Site_order, Site_label)
  
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
  
  summary_tbl <- summary_tbl %>%
    arrange(Site_order, Site_label)
  summary_tbl$Site_label <- factor(summary_tbl$Site_label, levels = site_label_levels, ordered = TRUE)
  summary_tbl
}

make_clonality_percent_barplot <- function(summary_tbl, lang = c("fr", "en")) {
  lang <- match.arg(lang)
  
  labels <- switch(
    lang,
    fr = list(
      title = "Clonalité par site",
      subtitle = "Pourcentage d'individus appartenant à des lignées multilocus répétées",
      x = "Site",
      y = "Clonalité (%)"
    ),
    en = list(
      title = "Clonality by site",
      subtitle = "Percentage of individuals belonging to repeated multilocus lineages",
      x = "Site",
      y = "Clonality (%)"
    )
  )
  
  ggplot(summary_tbl, aes(x = Site_label, y = Clonality_MLL_percent)) +
    geom_col(fill = "#2E8B57", width = 0.75) +
    scale_y_continuous(
      limits = c(0, 100),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(
      title = labels$title,
      subtitle = labels$subtitle,
      x = labels$x,
      y = labels$y
    ) +
    clonality_plot_theme()
}

save_clonality_percent_plot_dual_language <- function(summary_tbl, fig_dir) {
  plot_fr <- make_clonality_percent_barplot(summary_tbl, lang = "fr")
  plot_en <- make_clonality_percent_barplot(summary_tbl, lang = "en")
  
  png_fr <- file.path(fig_dir, "clonality_percent_per_site_fr.png")
  pdf_fr <- file.path(fig_dir, "clonality_percent_per_site_fr.pdf")
  png_en <- file.path(fig_dir, "clonality_percent_per_site_en.png")
  pdf_en <- file.path(fig_dir, "clonality_percent_per_site_en.pdf")
  
  ggsave(filename = png_fr, plot = plot_fr, width = 12, height = 7, dpi = 320)
  ggsave(filename = pdf_fr, plot = plot_fr, width = 12, height = 7)
  ggsave(filename = png_en, plot = plot_en, width = 12, height = 7, dpi = 320)
  ggsave(filename = pdf_en, plot = plot_en, width = 12, height = 7)
  
  list(
    fr_png = png_fr, fr_pdf = pdf_fr,
    en_png = png_en, en_pdf = pdf_en
  )
}

site_clonality_summary <- build_site_clonality_summary(clonality_df, df_ids_tbl = df_ids)

clonality_summary_table_numeric <- build_article_clonality_summary_table(
  site_clonality_summary,
  use_site_lookup_order = !is.null(site_lookup),
  formatted = FALSE
)
clonality_summary_table <- build_article_clonality_summary_table(
  site_clonality_summary,
  use_site_lookup_order = !is.null(site_lookup),
  formatted = TRUE
)
validate_article_clonality_summary_table(clonality_summary_table, site_clonality_summary)
validate_article_clonality_summary_table(clonality_summary_table_numeric, site_clonality_summary)
clonality_summary_table_csv <- file.path(TABLES_DIR, "clonality_summary_table.csv")
clonality_summary_table_xlsx <- file.path(TABLES_DIR, "clonality_summary_table.xlsx")
clonality_summary_table_docx <- file.path(TABLES_DIR, "clonality_summary_table.docx")
clonality_summary_table_numeric_csv <- file.path(TABLES_DIR, "clonality_summary_table_numeric.csv")
write.csv(clonality_summary_table, clonality_summary_table_csv, row.names = FALSE, na = "NA")
write.csv(clonality_summary_table_numeric, clonality_summary_table_numeric_csv, row.names = FALSE, na = "NA")
write_article_clonality_summary_xlsx(clonality_summary_table, clonality_summary_table_xlsx)
clonality_summary_table_docx_written <- write_article_clonality_summary_docx(
  clonality_summary_table,
  clonality_summary_table_docx
)
verify_article_clonality_summary_exports(clonality_summary_table, clonality_summary_table_csv, clonality_summary_table_xlsx)
print_article_clonality_summary_table(clonality_summary_table)

cat("\n[01_clonality] Tableau résumé utilisé pour le barplot de clonalité (%) par site:\n")
print(site_clonality_summary)

site_clonality_summary_file <- file.path(TABLES_DIR, "clonality_percent_per_site_summary.csv")
write.csv(site_clonality_summary, site_clonality_summary_file, row.names = FALSE)

if (!is.null(site_lookup)) {
  validate_clonality_plot_site_labels(site_clonality_summary, site_lookup)
} else {
  message("[01_clonality] Skipping site_lookup plot-order validation because site_lookup is unavailable.")
}

clonality_plot_files <- list(
  fr_png = NA_character_,
  fr_pdf = NA_character_,
  en_png = NA_character_,
  en_pdf = NA_character_
)
clonality_plot_files <- save_clonality_percent_plot_dual_language(
  summary_tbl = site_clonality_summary,
  fig_dir = FIGURES_DIR
)

saved_output_files <- c(
  out_file,
  clonality_summary_table_csv,
  clonality_summary_table_xlsx,
  clonality_summary_table_docx_written,
  clonality_summary_table_numeric_csv,
  assign_file,
  site_clonality_summary_file,
  unname(unlist(clonality_plot_files, use.names = FALSE))
)
saved_output_files <- saved_output_files[!is.na(saved_output_files) & nzchar(saved_output_files)]
for (saved_output_file in saved_output_files) {
  message("[01_clonality] Saved: ", saved_output_file)
}