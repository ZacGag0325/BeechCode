# scripts/sampled_stem_dbh_distribution.R
############################################################
# Sampled-stem DBH distribution for the clonality/genetic dataset
#
# Purpose:
# - Describe the diameter distribution of the Fagus grandifolia stems used
#   in the full clonality/genetic dataset (gi; not clone-corrected gi_mll).
# - Export a DBH-class summary table, an optional site-level DBH summary,
#   and a clean bar graph for article/thesis supplementary use.
#
# Outputs:
# - outputs/tables/sampled_stem_dbh_distribution.csv
# - outputs/tables/sampled_stem_dbh_distribution.xlsx
# - outputs/tables/sampled_stem_dbh_distribution.docx
# - outputs/tables/sampled_stem_dbh_distribution_by_site.csv
# - outputs/tables/sampled_stem_dbh_distribution_by_site.xlsx
# - outputs/figures/sampled_stem_dbh_distribution.png
# - outputs/figures/sampled_stem_dbh_distribution.pdf
# - outputs/figures/supplementary/sampled_stem_dbh_distribution.png
# - outputs/figures/supplementary/sampled_stem_dbh_distribution.pdf
############################################################

suppressPackageStartupMessages({
  library(poppr)
  library(adegenet)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readxl)
  library(openxlsx)
})

source("scripts/_load_objects.R")

message("[sampled_stem_dbh_distribution] Project root: ", PROJECT_ROOT)
message("[sampled_stem_dbh_distribution] Using full filtered genetic/clonality genind object: gi")
message("[sampled_stem_dbh_distribution] Object directory: ", OBJ_DIR)

SITE_LOOKUP_REQUIRED_COLUMNS <- c("Site", "Site_label", "Region", "Site_order", "Numéro_Population")
SITE_LOOKUP_SHEET <- "site_lookup"
GENETIC_DBH_SHEET <- "genetique"
DBH_CLASS_LEVELS <- c(paste0(1:8, " cm"), ">8 cm")
SITE_LABEL_LEVELS_SOUTH_TO_NORTH <- c(paste0("S", 1:6), paste0("N", 1:6))
DBH_COLUMN_CHOICES <- c(
  "DBH", "dbh", "Dbh", "dbh_cm", "DBH_cm", "Dbh_cm",
  "dhp", "DHP", "Dhp", "dhp_cm", "DHP_cm", "Dhp_cm",
  "dhp_tige", "DHP_tige", "Dhp_tige", "dhp_tige_cm", "DHP_tige_cm", "Dhp_tige_cm"
)

# -----------------------------
# General helpers
# -----------------------------
normalize_lookup_key <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("\uFEFF", "", x, fixed = TRUE)
  x <- gsub("[[:cntrl:]]", "", x)
  x[is.na(x)] <- ""
  x
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

find_preferred_column_ci_norm <- function(df, candidates) {
  cn <- names(df)
  cn_norm <- normalize_ascii_key(cn)
  candidate_norm <- normalize_ascii_key(candidates)
  idx <- match(candidate_norm, cn_norm, nomatch = 0)
  idx <- idx[idx > 0]
  if (length(idx) == 0) NA_character_ else cn[idx[1]]
}

find_dbh_col <- function(df) {
  find_preferred_column_ci_norm(df, DBH_COLUMN_CHOICES)
}

find_individual_id_col <- function(df) {
  find_preferred_column_ci_norm(
    df,
    c(
      DF_IDS_ID_CHOICES,
      "individual_id", "Individual_ID", "individu", "Individu", "id_individu", "ID_Individu",
      "sample_id", "Sample_ID", "sample_name", "Sample_Name", "echantillon", "échantillon",
      "id_echantillon", "id_échantillon", "identifiant", "genetic_id", "Genetic_ID",
      "Nom_Labo_Échantillons", "Nom_Labo_Echantillons", "Nom_Labo_Echantillon",
      "ID", "Id", "id", "Nom", "nom", "Name", "name"
    )
  )
}

find_tree_number_col <- function(df) {
  find_preferred_column_ci_norm(
    df,
    c(
      "tree", "Tree", "tree_id", "Tree_ID", "tree_number", "Tree_Number", "tree_no", "Tree_No",
      "arbre", "Arbre", "no_arbre", "No_Arbre", "numero_arbre", "Numéro_Arbre",
      "numero_tree", "num_tree", "no_tree", "id_arbre", "id_tige", "Id_tige", "stem", "stem_id",
      "numero_individu", "num_individu", "no_individu", "id_individu",
      "sample_number", "sample_no", "numero", "num", "no"
    )
  )
}

find_site_code_col <- function(df) {
  candidates <- c(
    "Site", "site", "Code_site", "code_site", "site_code", "Site_code", "CodeSite",
    "pop", "Pop", "population", "Population", "localite", "Localite", "localité", "Localité",
    "Numéro_Population", "Numero_Population", "numero_population", "Population_Number"
  )
  out <- find_preferred_column_ci_norm(df, candidates)
  if (!is.na(out)) return(out)
  cn_norm <- normalize_ascii_key(names(df))
  idx <- which(grepl("(^site$|site_?code|code_?site|^pop$|population|localite|locality)", cn_norm))
  if (length(idx) == 0) NA_character_ else names(df)[idx[1]]
}

coerce_numeric_dbh <- function(x) {
  if (is.numeric(x)) return(as.numeric(x))
  x_chr <- trimws(as.character(x))
  x_chr[x_chr %in% c("", "NA", "N/A", "NaN", "NULL", "null", "na")] <- NA_character_
  x_chr <- gsub(",", ".", x_chr, fixed = TRUE)
  x_chr <- gsub("[^0-9eE+.-]+", "", x_chr)
  suppressWarnings(as.numeric(x_chr))
}

normalize_exact_genetic_id <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("\uFEFF", "", x, fixed = TRUE)
  x <- gsub("[[:cntrl:]]", "", x)
  x[x %in% c("", "NA", "N/A", "NaN", "NULL", "null", "na")] <- NA_character_
  x
}

normalize_tree_number <- function(x) {
  x_chr <- trimws(as.character(x))
  x_chr[x_chr %in% c("", "NA", "NaN", "NULL")] <- NA_character_
  x_num <- suppressWarnings(as.numeric(gsub(",", ".", x_chr, fixed = TRUE)))
  x_chr <- ifelse(!is.na(x_num) & abs(x_num - round(x_num)) < .Machine$double.eps^0.5, as.character(as.integer(round(x_num))), x_chr)
  x_chr
}

strip_site_prefix_from_individual_id <- function(individual_id, site) {
  out <- trimws(as.character(individual_id))
  site_chr <- trimws(as.character(site))
  if (length(site_chr) == 1L && length(out) != 1L) site_chr <- rep(site_chr, length(out))
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
  if (!is.null(site)) x <- strip_site_prefix_from_individual_id(x, site)
  out <- sub(".*?([0-9]+)[^0-9]*$", "\\1", x)
  unchanged <- out == x
  unchanged[is.na(unchanged)] <- TRUE
  out[unchanged | is.na(x) | !grepl("[0-9]", x)] <- NA_character_
  normalize_tree_number(out)
}

compact_match_key <- function(...) {
  parts <- list(...)
  raw <- do.call(paste, c(lapply(parts, as.character), sep = "_"))
  key <- toupper(gsub("[^A-Za-z0-9]+", "", raw))
  key[is.na(key)] <- ""
  key
}

format_percent <- function(x, digits = 1) {
  ifelse(is.na(x), NA_character_, paste0(sprintf(paste0("%.", digits, "f"), x), "%"))
}

# -----------------------------
# Site lookup helpers
# -----------------------------
find_site_lookup_workbook <- function(raw_dir = file.path(PROJECT_ROOT, "data", "raw"), sheet = SITE_LOOKUP_SHEET) {
  if (!dir.exists(raw_dir)) return(NULL)
  excel_files <- list.files(raw_dir, pattern = "\\.(xlsx|xls)$", full.names = TRUE, ignore.case = TRUE)
  if (length(excel_files) == 0) return(NULL)
  has_lookup <- vapply(
    excel_files,
    function(path) sheet %in% readxl::excel_sheets(path),
    logical(1)
  )
  lookup_files <- excel_files[has_lookup]
  if (length(lookup_files) == 0) return(NULL)
  lookup_files[1]
}

load_site_lookup <- function() {
  workbook <- find_site_lookup_workbook()
  if (is.null(workbook)) {
    message("[sampled_stem_dbh_distribution] site_lookup unavailable; using Site values already stored in df_ids/gi metadata.")
    return(NULL)
  }
  lookup <- suppressMessages(readxl::read_excel(workbook, sheet = SITE_LOOKUP_SHEET)) %>%
    as.data.frame(stringsAsFactors = FALSE)
  missing_cols <- setdiff(SITE_LOOKUP_REQUIRED_COLUMNS, names(lookup))
  if (length(missing_cols) > 0) {
    stop(
      "[sampled_stem_dbh_distribution] site_lookup is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
  lookup <- lookup %>%
    mutate(
      Site = normalize_lookup_key(Site),
      Site_label = normalize_lookup_key(Site_label),
      Region = normalize_lookup_key(Region),
      Site_order = suppressWarnings(as.numeric(Site_order)),
      Numéro_Population = normalize_lookup_key(Numéro_Population)
    ) %>%
    arrange(Site_order, Site_label)
  message("[sampled_stem_dbh_distribution] Reading site_lookup from: ", workbook)
  lookup
}

apply_site_lookup <- function(stem_tbl, lookup) {
  out <- stem_tbl %>% select(-any_of(c("Site_label", "Region", "Site_order")))
  out$Site <- normalize_lookup_key(out$Site)
  if (is.null(lookup)) {
    out$Site_label <- out$Site
    out$Region <- NA_character_
    out$Site_order <- match(out$Site_label, SITE_LABEL_LEVELS_SOUTH_TO_NORTH)
    missing_order <- is.na(out$Site_order)
    if (any(missing_order)) out$Site_order[missing_order] <- length(SITE_LABEL_LEVELS_SOUTH_TO_NORTH) + match(out$Site_label[missing_order], unique(out$Site_label[missing_order]))
    return(out)
  }
  missing_sites <- sort(unique(out$Site[!out$Site %in% lookup$Site & nzchar(out$Site)]))
  if (length(missing_sites) > 0) {
    stop(
      "[sampled_stem_dbh_distribution] Site values are missing from site_lookup$Site: ",
      paste(missing_sites, collapse = ", "),
      call. = FALSE
    )
  }
  out %>%
    left_join(
      lookup %>% select(Site, Site_label, Region, Site_order),
      by = "Site"
    )
}

# -----------------------------
# Source and match DBH values
# -----------------------------
build_stem_table_from_gi <- function(gobj, df_ids_tbl) {
  aligned <- align_df_ids_to_genind(gobj, df_ids_tbl, context = "[sampled_stem_dbh_distribution gi]")
  data.frame(
    Individual = adegenet::indNames(gobj),
    Site = as.character(aligned$Site),
    stringsAsFactors = FALSE
  ) %>%
    bind_cols(aligned %>% select(-ind_id, -Site))
}

use_dbh_from_table_if_available <- function(stem_tbl, candidate_tbl, source_name) {
  dbh_col <- find_dbh_col(candidate_tbl)
  if (is.na(dbh_col)) return(NULL)
  id_col <- find_individual_id_col(candidate_tbl)
  if (is.na(id_col)) {
    if (nrow(candidate_tbl) == nrow(stem_tbl)) {
      out <- stem_tbl
      out$DBH_cm <- coerce_numeric_dbh(candidate_tbl[[dbh_col]])
      return(list(stem_tbl = out, dbh_col = dbh_col, source_label = paste0(source_name, " (row-order match)")))
    }
    return(NULL)
  }
  source_ids <- normalize_exact_genetic_id(candidate_tbl[[id_col]])
  retained_ids <- normalize_exact_genetic_id(stem_tbl$Individual)
  source_use <- data.frame(
    Individual_key = source_ids,
    DBH_cm = coerce_numeric_dbh(candidate_tbl[[dbh_col]]),
    stringsAsFactors = FALSE
  ) %>%
    filter(!is.na(Individual_key), nzchar(Individual_key))
  source_retained <- source_use[source_use$Individual_key %in% retained_ids, , drop = FALSE]
  if (nrow(source_retained) == 0) return(NULL)
  if (anyDuplicated(source_retained$Individual_key)) {
    dup_ids <- unique(source_retained$Individual_key[duplicated(source_retained$Individual_key)])
    stop(
      "[sampled_stem_dbh_distribution] ", source_name, " has duplicate retained individual IDs in column '", id_col,
      "'. Example duplicate(s): ", paste(head(dup_ids, 10), collapse = ", "),
      call. = FALSE
    )
  }
  idx <- match(retained_ids, source_use$Individual_key)
  match_rate <- sum(!is.na(idx)) / nrow(stem_tbl)
  if (match_rate < 0.95) return(NULL)
  out <- stem_tbl
  out$DBH_cm <- ifelse(!is.na(idx), source_use$DBH_cm[idx], NA_real_)
  list(
    stem_tbl = out,
    dbh_col = dbh_col,
    source_label = paste0(source_name, "$", dbh_col, " matched by ", source_name, "$", id_col)
  )
}

find_raw_dbh_workbook <- function(raw_dir = file.path(PROJECT_ROOT, "data", "raw")) {
  if (!dir.exists(raw_dir)) return(NULL)
  excel_files <- list.files(raw_dir, pattern = "\\.(xlsx|xls)$", full.names = TRUE, ignore.case = TRUE)
  if (length(excel_files) == 0) return(NULL)
  candidates <- list()
  for (path in excel_files) {
    sheets <- readxl::excel_sheets(path)
    sheet_order <- c(sheets[normalize_ascii_key(sheets) == normalize_ascii_key(GENETIC_DBH_SHEET)], sheets)
    sheet_order <- unique(sheet_order)
    for (sheet in sheet_order) {
      raw_df <- suppressMessages(readxl::read_excel(path, sheet = sheet)) %>% as.data.frame(stringsAsFactors = FALSE)
      dbh_col <- find_dbh_col(raw_df)
      if (!is.na(dbh_col)) {
        candidates[[length(candidates) + 1]] <- list(path = path, sheet = sheet, raw_df = raw_df, dbh_col = dbh_col)
      }
    }
  }
  if (length(candidates) == 0) return(NULL)
  scores <- vapply(candidates, function(x) {
    score <- 0
    if (normalize_ascii_key(x$sheet) == normalize_ascii_key(GENETIC_DBH_SHEET)) score <- score + 100
    if (normalize_ascii_key(x$dbh_col) %in% normalize_ascii_key(c("Dhp_tige", "dhp_tige"))) score <- score + 50
    if (grepl("donnee|modifie|west|summer|copie", normalize_ascii_key(basename(x$path)))) score <- score + 10
    score
  }, numeric(1))
  candidates[[which.max(scores)]]
}

match_raw_dbh_to_stems <- function(stem_tbl) {
  source <- find_raw_dbh_workbook()
  if (is.null(source)) return(NULL)
  raw_df <- source$raw_df
  dbh_col <- source$dbh_col
  raw_dbh <- coerce_numeric_dbh(raw_df[[dbh_col]])
  retained_ids <- normalize_exact_genetic_id(stem_tbl$Individual)
  candidate_id_cols <- setdiff(names(raw_df), dbh_col)
  id_report <- lapply(candidate_id_cols, function(col) {
    values <- normalize_exact_genetic_id(raw_df[[col]])
    matched_values <- values[!is.na(values) & values %in% retained_ids]
    data.frame(
      column = col,
      n_exact_matches = sum(retained_ids %in% values, na.rm = TRUE),
      match_rate = sum(retained_ids %in% values, na.rm = TRUE) / nrow(stem_tbl),
      n_duplicate_matched_ids = sum(duplicated(matched_values)),
      stringsAsFactors = FALSE
    )
  }) %>%
    bind_rows() %>%
    arrange(desc(n_exact_matches), column)
  best_id <- if (nrow(id_report) > 0) id_report[1, , drop = FALSE] else NULL
  if (!is.null(best_id) && is.finite(best_id$match_rate) && best_id$match_rate >= 0.95) {
    if (best_id$n_duplicate_matched_ids > 0) {
      values <- normalize_exact_genetic_id(raw_df[[best_id$column]])
      dup_ids <- unique(values[!is.na(values) & values %in% retained_ids & duplicated(values[!is.na(values) & values %in% retained_ids])])
      stop(
        "[sampled_stem_dbh_distribution] Raw DBH source has duplicate retained IDs in column '", best_id$column,
        "'. Example duplicate(s): ", paste(head(dup_ids, 10), collapse = ", "),
        call. = FALSE
      )
    }
    values <- normalize_exact_genetic_id(raw_df[[best_id$column]])
    idx <- match(retained_ids, values)
    out <- stem_tbl
    out$DBH_cm <- ifelse(!is.na(idx), raw_dbh[idx], NA_real_)
    return(list(
      stem_tbl = out,
      dbh_col = dbh_col,
      source_label = paste0(source$path, " | sheet '", source$sheet, "' | column '", dbh_col, "' matched by exact retained ID column '", best_id$column, "'")
    ))
  }
  
  raw_site_col <- find_site_code_col(raw_df)
  raw_tree_col <- find_tree_number_col(raw_df)
  stem_tree_col <- find_tree_number_col(stem_tbl)
  if (!is.na(raw_site_col) && !is.na(raw_tree_col)) {
    stem_tree <- if (!is.na(stem_tree_col)) normalize_tree_number(stem_tbl[[stem_tree_col]]) else derive_tree_number_from_individual_id(stem_tbl$Individual, stem_tbl$Site)
    raw_key <- compact_match_key(normalize_lookup_key(raw_df[[raw_site_col]]), normalize_tree_number(raw_df[[raw_tree_col]]))
    stem_key <- compact_match_key(normalize_lookup_key(stem_tbl$Site), stem_tree)
    source_use <- data.frame(match_key = raw_key, DBH_cm = raw_dbh, stringsAsFactors = FALSE) %>%
      filter(nzchar(match_key))
    source_retained <- source_use[source_use$match_key %in% stem_key, , drop = FALSE]
    if (nrow(source_retained) > 0 && !anyDuplicated(source_retained$match_key)) {
      idx <- match(stem_key, source_use$match_key)
      match_rate <- sum(!is.na(idx)) / nrow(stem_tbl)
      if (match_rate >= 0.95) {
        out <- stem_tbl
        out$DBH_cm <- ifelse(!is.na(idx), source_use$DBH_cm[idx], NA_real_)
        return(list(
          stem_tbl = out,
          dbh_col = dbh_col,
          source_label = paste0(source$path, " | sheet '", source$sheet, "' | column '", dbh_col, "' matched by site + tree number")
        ))
      }
    }
  }
  
  top_report <- if (!is.null(best_id)) paste(utils::capture.output(print(head(id_report, 10), row.names = FALSE)), collapse = "\n") else "No ID-like columns available."
  stop(
    "[sampled_stem_dbh_distribution] A DBH column was found in a raw workbook, but DBH values could not be matched to at least 95% of genetic/clonality stems.\n",
    "Raw source: ", source$path, " | sheet '", source$sheet, "' | DBH column '", dbh_col, "'\n",
    "Top exact-ID matching report:\n", top_report,
    call. = FALSE
  )
}

resolve_sampled_stem_dbh <- function() {
  stem_tbl <- build_stem_table_from_gi(gi, df_ids)
  
  from_df_ids <- use_dbh_from_table_if_available(stem_tbl, stem_tbl, "df_ids")
  if (!is.null(from_df_ids)) return(from_df_ids)
  
  if (exists("meta", inherits = TRUE) && is.data.frame(meta)) {
    from_meta <- use_dbh_from_table_if_available(stem_tbl, meta, "meta")
    if (!is.null(from_meta)) return(from_meta)
  }
  
  raw_match <- match_raw_dbh_to_stems(stem_tbl)
  if (!is.null(raw_match)) return(raw_match)
  
  stop(
    "[sampled_stem_dbh_distribution] No DBH column could be found in df_ids, meta, or raw Excel workbooks. ",
    "Accepted DBH column names include: ", paste(DBH_COLUMN_CHOICES, collapse = ", "),
    call. = FALSE
  )
}

# -----------------------------
# Summaries, plotting, and exports
# -----------------------------
bin_dbh_class <- function(dbh_cm) {
  dbh_int <- round(dbh_cm)
  dbh_int[!is.na(dbh_int) & dbh_int < 1] <- 1
  out <- ifelse(is.na(dbh_int), NA_character_, ifelse(dbh_int > 8, ">8 cm", paste0(dbh_int, " cm")))
  factor(out, levels = DBH_CLASS_LEVELS, ordered = TRUE)
}

summarise_dbh_classes <- function(stem_tbl) {
  n_with_dbh <- sum(!is.na(stem_tbl$DBH_cm))
  stem_tbl %>%
    filter(!is.na(DBH_class)) %>%
    count(DBH_class, name = "Number of stems") %>%
    complete(DBH_class = factor(DBH_CLASS_LEVELS, levels = DBH_CLASS_LEVELS, ordered = TRUE), fill = list(`Number of stems` = 0L)) %>%
    arrange(DBH_class) %>%
    mutate(
      `DBH class` = as.character(DBH_class),
      `Percent of stems` = ifelse(n_with_dbh > 0, 100 * `Number of stems` / n_with_dbh, NA_real_),
      `Cumulative number` = cumsum(`Number of stems`),
      `Cumulative percent` = ifelse(n_with_dbh > 0, 100 * `Cumulative number` / n_with_dbh, NA_real_)
    ) %>%
    select(`DBH class`, `Number of stems`, `Percent of stems`, `Cumulative number`, `Cumulative percent`)
}

summarise_dbh_by_site <- function(stem_tbl) {
  stem_tbl %>%
    mutate(
      Site = as.character(Site_label),
      Site = factor(Site, levels = unique(as.character(stem_tbl$Site_label[order(stem_tbl$Site_order, stem_tbl$Site_label)])), ordered = TRUE)
    ) %>%
    group_by(Site) %>%
    summarise(
      N = dplyr::n(),
      `mean DBH` = ifelse(all(is.na(DBH_cm)), NA_real_, mean(DBH_cm, na.rm = TRUE)),
      `median DBH` = ifelse(all(is.na(DBH_cm)), NA_real_, stats::median(DBH_cm, na.rm = TRUE)),
      `minimum DBH` = ifelse(all(is.na(DBH_cm)), NA_real_, min(DBH_cm, na.rm = TRUE)),
      `maximum DBH` = ifelse(all(is.na(DBH_cm)), NA_real_, max(DBH_cm, na.rm = TRUE)),
      `percent 1 cm` = ifelse(sum(!is.na(DBH_cm)) > 0, 100 * sum(DBH_class == "1 cm", na.rm = TRUE) / sum(!is.na(DBH_cm)), NA_real_),
      `percent ≤2 cm` = ifelse(sum(!is.na(DBH_cm)) > 0, 100 * sum(!is.na(DBH_cm) & DBH_cm_rounded <= 2, na.rm = TRUE) / sum(!is.na(DBH_cm)), NA_real_),
      `percent ≤5 cm` = ifelse(sum(!is.na(DBH_cm)) > 0, 100 * sum(!is.na(DBH_cm) & DBH_cm_rounded <= 5, na.rm = TRUE) / sum(!is.na(DBH_cm)), NA_real_),
      .groups = "drop"
    ) %>%
    mutate(Site = as.character(Site))
}

write_summary_xlsx <- function(class_summary, site_summary, class_path, site_path) {
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "DBH classes")
  openxlsx::writeData(wb, "DBH classes", class_summary)
  header_style <- openxlsx::createStyle(textDecoration = "bold", halign = "center", valign = "center", border = "bottom", fgFill = "#D9EAD3")
  body_style <- openxlsx::createStyle(halign = "center", valign = "center")
  openxlsx::addStyle(wb, "DBH classes", header_style, rows = 1, cols = seq_len(ncol(class_summary)), gridExpand = TRUE)
  openxlsx::addStyle(wb, "DBH classes", body_style, rows = 2:(nrow(class_summary) + 1), cols = seq_len(ncol(class_summary)), gridExpand = TRUE)
  openxlsx::freezePane(wb, "DBH classes", firstRow = TRUE)
  openxlsx::setColWidths(wb, "DBH classes", cols = seq_len(ncol(class_summary)), widths = "auto")
  openxlsx::saveWorkbook(wb, class_path, overwrite = TRUE)
  
  wb_site <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb_site, "DBH by site")
  openxlsx::writeData(wb_site, "DBH by site", site_summary)
  openxlsx::addStyle(wb_site, "DBH by site", header_style, rows = 1, cols = seq_len(ncol(site_summary)), gridExpand = TRUE)
  if (nrow(site_summary) > 0) {
    openxlsx::addStyle(wb_site, "DBH by site", body_style, rows = 2:(nrow(site_summary) + 1), cols = seq_len(ncol(site_summary)), gridExpand = TRUE)
  }
  openxlsx::freezePane(wb_site, "DBH by site", firstRow = TRUE)
  openxlsx::setColWidths(wb_site, "DBH by site", cols = seq_len(ncol(site_summary)), widths = "auto")
  openxlsx::saveWorkbook(wb_site, site_path, overwrite = TRUE)
}

sampled_dbh_plot_theme <- function(base_size = 22) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 4, margin = margin(b = 6)),
      plot.subtitle = element_text(size = base_size, margin = margin(b = 12)),
      axis.title = element_text(size = base_size + 2),
      axis.text = element_text(size = base_size),
      axis.line = element_line(linewidth = 0.8, colour = "black"),
      axis.ticks = element_line(linewidth = 0.8, colour = "black"),
      axis.ticks.length = grid::unit(0.25, "cm"),
      panel.border = element_blank(),
      plot.background = element_blank(),
      panel.background = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(linewidth = 0.25, colour = "grey88")
    )
}

make_dbh_distribution_plot <- function(class_summary, total_n) {
  plot_tbl <- class_summary %>%
    mutate(
      `DBH class` = factor(`DBH class`, levels = DBH_CLASS_LEVELS, ordered = TRUE),
      label = ifelse(`Number of stems` > 0, paste0(`Number of stems`, "\n(", sprintf("%.1f", `Percent of stems`), "%)"), "")
    )
  y_max <- max(plot_tbl$`Number of stems`, na.rm = TRUE)
  ggplot(plot_tbl, aes(x = `DBH class`, y = `Number of stems`)) +
    geom_col(fill = "#2E8B57", width = 0.75) +
    geom_text(aes(label = label), vjust = -0.25, size = 5.0, lineheight = 0.9) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.14)), limits = c(0, y_max * 1.18)) +
    labs(
      title = "DBH distribution of sampled Fagus grandifolia stems",
      subtitle = paste0("Total N = ", total_n, " sampled stems with DBH values"),
      x = "DBH class",
      y = "Number of sampled stems"
    ) +
    sampled_dbh_plot_theme()
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
  paste0("<w:tc><w:tcPr><w:tcW w:w=\"2400\" w:type=\"dxa\"/></w:tcPr><w:p><w:r>", bold_xml, "<w:t>", text, "</w:t></w:r></w:p></w:tc>")
}

word_paragraph_xml <- function(value, bold = FALSE) {
  bold_xml <- if (bold) "<w:rPr><w:b/></w:rPr>" else ""
  paste0("<w:p><w:r>", bold_xml, "<w:t>", xml_escape_word(value), "</w:t></w:r></w:p>")
}

write_basic_docx_table <- function(df, path, title, note = NULL) {
  tmp_dir <- tempfile("sampled_dbh_docx_")
  dir.create(file.path(tmp_dir, "_rels"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(tmp_dir, "word", "_rels"), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  df <- as.data.frame(lapply(df, as.character), stringsAsFactors = FALSE, check.names = FALSE)
  rows <- apply(df, 1, function(row) paste0("<w:tr>", paste(vapply(row, word_cell_xml, character(1)), collapse = ""), "</w:tr>"))
  header <- paste0("<w:tr>", paste(vapply(names(df), word_cell_xml, character(1), bold = TRUE), collapse = ""), "</w:tr>")
  note_xml <- if (!is.null(note) && nzchar(note)) word_paragraph_xml(note) else ""
  document_xml <- paste0(
    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
    '<w:document xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><w:body>',
    word_paragraph_xml(title, bold = TRUE),
    note_xml,
    '<w:tbl><w:tblPr><w:tblStyle w:val="TableGrid"/><w:tblW w:w="0" w:type="auto"/><w:tblBorders>',
    '<w:top w:val="single" w:sz="4" w:space="0" w:color="auto"/><w:left w:val="single" w:sz="4" w:space="0" w:color="auto"/>',
    '<w:bottom w:val="single" w:sz="4" w:space="0" w:color="auto"/><w:right w:val="single" w:sz="4" w:space="0" w:color="auto"/>',
    '<w:insideH w:val="single" w:sz="4" w:space="0" w:color="auto"/><w:insideV w:val="single" w:sz="4" w:space="0" w:color="auto"/>',
    '</w:tblBorders></w:tblPr>', header, paste(rows, collapse = ""), '</w:tbl>',
    '<w:sectPr><w:pgSz w:w="12240" w:h="15840"/><w:pgMar w:top="720" w:right="720" w:bottom="720" w:left="720" w:header="360" w:footer="360" w:gutter="0"/></w:sectPr>',
    '</w:body></w:document>'
  )
  writeLines(c('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>', '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">', '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>', '<Default Extension="xml" ContentType="application/xml"/>', '<Override PartName="/word/document.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.document.main+xml"/>', '</Types>'), file.path(tmp_dir, "[Content_Types].xml"), useBytes = TRUE)
  writeLines(c('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>', '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">', '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="word/document.xml"/>', '</Relationships>'), file.path(tmp_dir, "_rels", ".rels"), useBytes = TRUE)
  writeLines(c('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>', '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"></Relationships>'), file.path(tmp_dir, "word", "_rels", "document.xml.rels"), useBytes = TRUE)
  writeLines(document_xml, file.path(tmp_dir, "word", "document.xml"), useBytes = TRUE)
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(tmp_dir)
  if (file.exists(path)) unlink(path)
  utils::zip(zipfile = path, files = list.files(tmp_dir, recursive = TRUE, all.files = TRUE, no.. = TRUE), flags = "-q")
  if (!file.exists(path)) stop("Failed to create Word document at: ", path, call. = FALSE)
  invisible(path)
}

add_flextable_to_docx <- function(doc, ft) {
  if (requireNamespace("flextable", quietly = TRUE) && "body_add_flextable" %in% getNamespaceExports("flextable")) {
    return(flextable::body_add_flextable(doc, value = ft))
  }
  if (requireNamespace("officer", quietly = TRUE) && "body_add_flextable" %in% getNamespaceExports("officer")) {
    return(officer::body_add_flextable(doc, value = ft))
  }
  stop(
    "No compatible body_add_flextable() function is exported by the installed officer/flextable packages.",
    call. = FALSE
  )
}

write_dbh_docx <- function(class_summary, path, figure_path) {
  title <- "Sampled-stem DBH distribution"
  note <- "This table summarizes DBH classes for the sampled Fagus grandifolia stems included in the full genetic/clonality dataset. DBH values are rounded to the nearest centimetre for class assignment; stems rounded above 8 cm are grouped as >8 cm."
  word_tbl <- class_summary %>%
    mutate(
      `Percent of stems` = sprintf("%.1f", `Percent of stems`),
      `Cumulative percent` = sprintf("%.1f", `Cumulative percent`)
    )
  
  status <- tryCatch(
    {
      if (requireNamespace("officer", quietly = TRUE) && requireNamespace("flextable", quietly = TRUE)) {
        ft <- flextable::flextable(word_tbl) %>%
          flextable::theme_booktabs() %>%
          flextable::autofit() %>%
          flextable::align(align = "center", part = "all") %>%
          flextable::bold(part = "header")
        doc <- officer::read_docx()
        doc <- officer::body_add_par(doc, title, style = "heading 1")
        doc <- officer::body_add_par(doc, note, style = "Normal")
        doc <- add_flextable_to_docx(doc, ft)
        if (file.exists(figure_path)) {
          doc <- officer::body_add_par(doc, "Figure. DBH distribution of sampled Fagus grandifolia stems.", style = "Normal")
          doc <- officer::body_add_img(doc, src = figure_path, width = 6.5, height = 4.5)
        }
        print(doc, target = path)
        "created with compatible officer/flextable export"
      } else if (requireNamespace("officer", quietly = TRUE)) {
        doc <- officer::read_docx()
        doc <- officer::body_add_par(doc, title, style = "heading 1")
        doc <- officer::body_add_par(doc, note, style = "Normal")
        doc <- officer::body_add_table(doc, value = as.data.frame(word_tbl), style = "table_template")
        if (file.exists(figure_path)) {
          doc <- officer::body_add_par(doc, "Figure. DBH distribution of sampled Fagus grandifolia stems.", style = "Normal")
          doc <- officer::body_add_img(doc, src = figure_path, width = 6.5, height = 4.5)
        }
        print(doc, target = path)
        "created with officer fallback"
      } else {
        write_basic_docx_table(word_tbl, path, title = title, note = note)
        "created with built-in docx fallback"
      }
    },
    error = function(e) {
      message("[sampled_stem_dbh_distribution] Package-based Word export failed: ", conditionMessage(e), ". Creating a basic Word-compatible .docx fallback.")
      write_basic_docx_table(word_tbl, path, title = title, note = note)
      "created with built-in docx fallback after package export failed"
    }
  )
  message("[sampled_stem_dbh_distribution] Saved Word table: ", path, " (", status, ")")
  invisible(path)
}

# -----------------------------
# Main analysis
# -----------------------------
resolved <- resolve_sampled_stem_dbh()
stem_dbh <- resolved$stem_tbl
stem_dbh <- apply_site_lookup(stem_dbh, load_site_lookup())
stem_dbh <- stem_dbh %>%
  mutate(
    DBH_cm = coerce_numeric_dbh(DBH_cm),
    DBH_cm_rounded = round(DBH_cm),
    DBH_cm_rounded = ifelse(!is.na(DBH_cm_rounded) & DBH_cm_rounded < 1, 1, DBH_cm_rounded),
    DBH_class = bin_dbh_class(DBH_cm)
  )

if (sum(!is.na(stem_dbh$DBH_cm)) == 0) {
  stop("[sampled_stem_dbh_distribution] DBH values were matched, but all matched values are missing/non-numeric.", call. = FALSE)
}

n_total <- nrow(stem_dbh)
n_missing <- sum(is.na(stem_dbh$DBH_cm))
n_with_dbh <- sum(!is.na(stem_dbh$DBH_cm))
missing_pct <- 100 * n_missing / n_total
if (n_missing > 0) {
  warning(
    "[sampled_stem_dbh_distribution] Missing DBH values: ", n_missing, " / ", n_total,
    " (", sprintf("%.1f", missing_pct), "%). Summaries and percentages use stems with non-missing DBH values as the denominator.",
    call. = FALSE
  )
}
if (missing_pct >= 10) {
  warning(
    "[sampled_stem_dbh_distribution] Many DBH values are missing (>=10%). Please verify DBH matching before interpreting the distribution.",
    call. = FALSE
  )
}

class_summary <- summarise_dbh_classes(stem_dbh)
site_summary <- summarise_dbh_by_site(stem_dbh)

n_1cm <- sum(stem_dbh$DBH_class == "1 cm", na.rm = TRUE)
pct_1cm <- 100 * n_1cm / n_with_dbh
n_le2 <- sum(!is.na(stem_dbh$DBH_cm_rounded) & stem_dbh$DBH_cm_rounded <= 2)
pct_le2 <- 100 * n_le2 / n_with_dbh
n_le5 <- sum(!is.na(stem_dbh$DBH_cm_rounded) & stem_dbh$DBH_cm_rounded <= 5)
pct_le5 <- 100 * n_le5 / n_with_dbh

csv_path <- file.path(TABLES_DIR, "sampled_stem_dbh_distribution.csv")
xlsx_path <- file.path(TABLES_DIR, "sampled_stem_dbh_distribution.xlsx")
docx_path <- file.path(TABLES_DIR, "sampled_stem_dbh_distribution.docx")
site_csv_path <- file.path(TABLES_DIR, "sampled_stem_dbh_distribution_by_site.csv")
site_xlsx_path <- file.path(TABLES_DIR, "sampled_stem_dbh_distribution_by_site.xlsx")
fig_png_path <- file.path(FIGURES_DIR, "sampled_stem_dbh_distribution.png")
fig_pdf_path <- file.path(FIGURES_DIR, "sampled_stem_dbh_distribution.pdf")
fig_supp_png_path <- file.path(FIGURES_SUPP_DIR, "sampled_stem_dbh_distribution.png")
fig_supp_pdf_path <- file.path(FIGURES_SUPP_DIR, "sampled_stem_dbh_distribution.pdf")

write.csv(class_summary, csv_path, row.names = FALSE, na = "NA")
write.csv(site_summary, site_csv_path, row.names = FALSE, na = "NA")
write_summary_xlsx(class_summary, site_summary, xlsx_path, site_xlsx_path)

p <- make_dbh_distribution_plot(class_summary, total_n = n_with_dbh)
ggsave(filename = fig_png_path, plot = p, width = 11, height = 7, dpi = 320)
ggsave(filename = fig_pdf_path, plot = p, width = 11, height = 7)
ggsave(filename = fig_supp_png_path, plot = p, width = 11, height = 7, dpi = 320)
ggsave(filename = fig_supp_pdf_path, plot = p, width = 11, height = 7)
write_dbh_docx(class_summary, docx_path, fig_png_path)

message("[sampled_stem_dbh_distribution] Object/file used: ", resolved$source_label)
message("[sampled_stem_dbh_distribution] Number of stems included from gi: ", n_total)
message("[sampled_stem_dbh_distribution] Number of stems with DBH values: ", n_with_dbh)
message("[sampled_stem_dbh_distribution] Number of missing DBH values: ", n_missing)
message("[sampled_stem_dbh_distribution] Detected DBH column: ", resolved$dbh_col)
message(
  "[sampled_stem_dbh_distribution] DBH range: ",
  sprintf("%.2f", min(stem_dbh$DBH_cm, na.rm = TRUE)), " to ",
  sprintf("%.2f", max(stem_dbh$DBH_cm, na.rm = TRUE)), " cm"
)
message("[sampled_stem_dbh_distribution] DBH binning: DBH values are rounded to the nearest whole centimetre; rounded values 1-8 are labelled 1 cm through 8 cm, and rounded values >8 are grouped as >8 cm. Values rounded below 1 are assigned to 1 cm.")
message("[sampled_stem_dbh_distribution] 1 cm stems: ", n_1cm, " / ", n_with_dbh, " (", sprintf("%.1f", pct_1cm), "%)")
message("[sampled_stem_dbh_distribution] Stems <=2 cm: ", n_le2, " / ", n_with_dbh, " (", sprintf("%.1f", pct_le2), "%)")
message("[sampled_stem_dbh_distribution] Stems <=5 cm: ", n_le5, " / ", n_with_dbh, " (", sprintf("%.1f", pct_le5), "%)")
message("[sampled_stem_dbh_distribution] Output files:")
for (path in c(csv_path, xlsx_path, docx_path, site_csv_path, site_xlsx_path, fig_png_path, fig_pdf_path, fig_supp_png_path, fig_supp_pdf_path)) {
  message("  - ", path)
}