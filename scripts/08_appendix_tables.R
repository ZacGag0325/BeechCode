# scripts/08_appendix_tables.R
############################################################
# Appendix-ready HWE and allele-frequency tables
#
# Creates Daphnee-Bernier-style appendix tables:
#   1) Hardy-Weinberg equilibrium by locus, overall and by site
#      - full detailed verification table
#      - compact thesis-display Word table like Appendix Table C.1
#   2) Allele frequencies by locus/allele, by site plus weighted and
#      non-weighted overall frequencies
#
# Primary dataset policy:
# - Uses gi_mll, the final MLL clone-corrected genind object loaded by
#   scripts/_load_objects.R, because this is the dataset used by the
#   genetic-diversity workflow and avoids pseudo-replication by clonemates.
# - Uses updated site labels from the site_lookup sheet when available.
#
# Output folder:
#   outputs/appendix_tables/
############################################################

suppressPackageStartupMessages({
  library(adegenet)
  library(pegas)
  library(dplyr)
  library(tidyr)
  library(readxl)
  library(openxlsx)
})

# ----------------------------
# Project-root helper
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
      if (file.exists(file.path(cur, "scripts", "08_appendix_tables.R"))) return(cur)
      parent <- dirname(cur)
      if (identical(parent, cur)) break
      cur <- parent
    }
  }
  
  stop("[08_appendix_tables] Cannot find project root containing scripts/08_appendix_tables.R", call. = FALSE)
}

PROJECT_ROOT <- find_project_root()
setwd(PROJECT_ROOT)

# ----------------------------
# User-adjustable settings
# ----------------------------
APPENDIX_DIR <- file.path(PROJECT_ROOT, "outputs", "appendix_tables")
dir.create(APPENDIX_DIR, recursive = TRUE, showWarnings = FALSE)

EXPECTED_SITE_LABELS <- c(paste0("S", 1:6), paste0("N", 1:6))
HWE_MONTE_CARLO_REPS <- 9999L
HWE_MIN_NON_MISSING_N <- 8L
HWE_MIN_UNIQUE_GENOTYPES <- 2L
ALPHA <- 0.05
FREQ_DIGITS <- 3L

# ----------------------------
# Console/file helpers
# ----------------------------
message_appendix <- function(...) {
  message("[08_appendix_tables] ", ...)
}

write_csv_out <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.csv(df, path, row.names = FALSE, na = "")
  message_appendix("Saved CSV: ", path)
  invisible(path)
}

write_xlsx_out <- function(sheets, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (is.data.frame(sheets)) sheets <- list(Table = sheets)
  openxlsx::write.xlsx(sheets, file = path, overwrite = TRUE, asTable = TRUE)
  message_appendix("Saved Excel: ", path)
  invisible(path)
}

# ----------------------------
# Search first: existing scripts and outputs
# ----------------------------
search_existing_appendix_tables <- function() {
  script_files <- list.files(file.path(PROJECT_ROOT, "scripts"), pattern = "\\.R$", full.names = TRUE)
  output_roots <- c(
    file.path(PROJECT_ROOT, "outputs"),
    file.path(PROJECT_ROOT, "results"),
    file.path(PROJECT_ROOT, "data", "derived")
  )
  output_files <- unlist(lapply(output_roots[dir.exists(output_roots)], function(root) {
    list.files(root, recursive = TRUE, full.names = TRUE)
  }), use.names = FALSE)
  
  script_hits <- character(0)
  for (path in script_files) {
    txt <- readLines(path, warn = FALSE)
    if (any(grepl("Hardy|HWE|Weinberg|allele.frequency|allelic.frequency|Appendix|appendix", txt, ignore.case = TRUE))) {
      script_hits <- c(script_hits, normalizePath(path, winslash = "/", mustWork = FALSE))
    }
  }
  
  output_hits <- output_files[grepl("hwe|hardy|weinberg|allele|frequency|appendix|supplement", basename(output_files), ignore.case = TRUE)]
  output_hits <- normalizePath(output_hits, winslash = "/", mustWork = FALSE)
  
  message_appendix("Existing script search found ", length(script_hits), " relevant script(s):")
  if (length(script_hits) == 0) {
    message_appendix("  - none")
  } else {
    for (x in script_hits) message_appendix("  - ", x)
  }
  
  message_appendix("Existing output search found ", length(output_hits), " potentially relevant output file(s):")
  if (length(output_hits) == 0) {
    message_appendix("  - none")
  } else {
    for (x in output_hits) message_appendix("  - ", x)
  }
  
  invisible(list(scripts = script_hits, outputs = output_hits))
}

# ----------------------------
# Site lookup helpers
# ----------------------------
normalize_name <- function(x) {
  x <- iconv(as.character(x), from = "", to = "ASCII//TRANSLIT")
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  gsub("^_|_$", "", x)
}

normalize_key <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("\uFEFF", "", x, fixed = TRUE)
  x <- gsub("[[:cntrl:]]", "", x)
  x[is.na(x)] <- ""
  x
}

find_site_lookup <- function() {
  source_dirs <- c(file.path(PROJECT_ROOT, "data", "raw"), file.path(PROJECT_ROOT, "inputs"))
  excel_files <- unlist(lapply(source_dirs[dir.exists(source_dirs)], function(root) {
    list.files(root, pattern = "\\.(xlsx|xls)$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  }), use.names = FALSE)
  
  for (path in excel_files) {
    sheets <- readxl::excel_sheets(path)
    hit <- which(normalize_name(sheets) == "site_lookup")
    if (length(hit) == 0) next
    
    lookup <- as.data.frame(readxl::read_excel(path, sheet = sheets[hit[1]], .name_repair = "minimal"), stringsAsFactors = FALSE)
    names_norm <- normalize_name(names(lookup))
    site_col <- names(lookup)[match("site", names_norm, nomatch = 0)]
    label_col <- names(lookup)[match("site_label", names_norm, nomatch = 0)]
    order_col <- names(lookup)[match("site_order", names_norm, nomatch = 0)]
    
    if (length(site_col) == 0 || length(label_col) == 0) next
    
    out <- data.frame(
      Site_original = normalize_key(lookup[[site_col]]),
      Site_label = normalize_key(lookup[[label_col]]),
      Site_order = if (length(order_col) > 0) suppressWarnings(as.numeric(lookup[[order_col]])) else NA_real_,
      stringsAsFactors = FALSE
    ) %>%
      filter(nzchar(Site_original), nzchar(Site_label)) %>%
      distinct(Site_original, .keep_all = TRUE)
    
    attr(out, "source_file") <- normalizePath(path, winslash = "/", mustWork = FALSE)
    attr(out, "source_sheet") <- sheets[hit[1]]
    return(out)
  }
  
  NULL
}

build_site_label_map <- function(df_ids_mll_tbl) {
  lookup <- find_site_lookup()
  
  if (!is.null(lookup) && nrow(lookup) > 0) {
    message_appendix(
      "Using site labels from site_lookup sheet: ", attr(lookup, "source_file"),
      " [", attr(lookup, "source_sheet"), "]"
    )
    return(lookup)
  }
  
  if ("Site_label" %in% names(df_ids_mll_tbl)) {
    message_appendix("No site_lookup sheet found; using Site_label already present in df_ids_mll.")
    return(
      df_ids_mll_tbl %>%
        transmute(Site_original = normalize_key(Site), Site_label = normalize_key(Site_label), Site_order = NA_real_) %>%
        filter(nzchar(Site_original), nzchar(Site_label)) %>%
        distinct(Site_original, .keep_all = TRUE)
    )
  }
  
  message_appendix("No site_lookup sheet found and no Site_label column found; using existing Site values.")
  data.frame(
    Site_original = sort(unique(normalize_key(df_ids_mll_tbl$Site))),
    Site_label = sort(unique(normalize_key(df_ids_mll_tbl$Site))),
    Site_order = NA_real_,
    stringsAsFactors = FALSE
  )
}

apply_site_labels <- function(site_vec, lookup) {
  site_chr <- normalize_key(site_vec)
  idx <- match(site_chr, lookup$Site_original)
  out <- lookup$Site_label[idx]
  out[is.na(out) | !nzchar(out)] <- site_chr[is.na(out) | !nzchar(out)]
  out
}

site_sort_levels <- function(site_labels, lookup) {
  site_labels <- unique(as.character(site_labels))
  expected_present <- EXPECTED_SITE_LABELS[EXPECTED_SITE_LABELS %in% site_labels]
  remaining <- setdiff(site_labels, expected_present)
  
  if (length(remaining) > 0 && !is.null(lookup) && "Site_order" %in% names(lookup)) {
    order_tbl <- lookup %>%
      filter(Site_label %in% remaining) %>%
      arrange(Site_order, Site_label)
    remaining_ordered <- unique(order_tbl$Site_label)
    remaining <- c(remaining_ordered, sort(setdiff(remaining, remaining_ordered)))
  } else {
    remaining <- sort(remaining)
  }
  
  c(expected_present, remaining)
}

# ----------------------------
# Genotype parsing helpers
# ----------------------------
split_genotype <- function(x) {
  x <- trimws(as.character(x))
  if (is.na(x) || x == "" || x %in% c("NA", "0", "0/0", "NA/NA", "-")) {
    return(c(NA_character_, NA_character_))
  }
  
  parts <- strsplit(x, "/", fixed = TRUE)[[1]]
  if (length(parts) != 2) return(c(NA_character_, NA_character_))
  parts <- trimws(parts)
  if (any(!nzchar(parts))) return(c(NA_character_, NA_character_))
  sort(parts)
}

parse_locus_genotypes <- function(geno_vec) {
  parsed <- t(vapply(geno_vec, split_genotype, character(2)))
  colnames(parsed) <- c("a1", "a2")
  parsed
}

# ----------------------------
# HWE helpers
# ----------------------------
calculate_hwe_chisq <- function(geno_vec) {
  parsed <- parse_locus_genotypes(geno_vec)
  keep <- !(is.na(parsed[, "a1"]) | is.na(parsed[, "a2"]))
  n_non_missing <- sum(keep)
  
  if (n_non_missing == 0) {
    return(list(chisq = NA_real_, df = NA_integer_, note = "chi_square_not_available_no_non_missing_genotypes"))
  }
  
  a1 <- parsed[keep, "a1"]
  a2 <- parsed[keep, "a2"]
  alleles <- sort(unique(c(a1, a2)))
  k <- length(alleles)
  
  if (k < 2) {
    return(list(chisq = NA_real_, df = NA_integer_, note = "chi_square_not_available_invariant_locus"))
  }
  
  allele_counts <- table(factor(c(a1, a2), levels = alleles))
  allele_freq <- as.numeric(allele_counts) / (2 * n_non_missing)
  names(allele_freq) <- alleles
  
  genotype_levels <- character(0)
  expected <- numeric(0)
  for (i in seq_along(alleles)) {
    for (j in i:k) {
      genotype_levels <- c(genotype_levels, paste(alleles[i], alleles[j], sep = "/"))
      exp_freq <- if (i == j) allele_freq[i]^2 else 2 * allele_freq[i] * allele_freq[j]
      expected <- c(expected, n_non_missing * exp_freq)
    }
  }
  
  observed_strings <- paste(a1, a2, sep = "/")
  observed <- as.numeric(table(factor(observed_strings, levels = genotype_levels)))
  use <- expected > 0
  chisq <- sum((observed[use] - expected[use])^2 / expected[use])
  df <- k * (k - 1) / 2
  
  list(chisq = as.numeric(chisq), df = as.integer(df), note = "chi_square_approximation_from_observed_expected_genotypes")
}

run_single_locus_hwe <- function(geno_vec,
                                 B = HWE_MONTE_CARLO_REPS,
                                 min_n = HWE_MIN_NON_MISSING_N,
                                 min_unique_genotypes = HWE_MIN_UNIQUE_GENOTYPES) {
  parsed <- parse_locus_genotypes(geno_vec)
  keep <- !(is.na(parsed[, "a1"]) | is.na(parsed[, "a2"]))
  n_non_missing <- sum(keep)
  missing_fraction <- 1 - (n_non_missing / length(geno_vec))
  
  base_row <- list(
    N = as.integer(n_non_missing),
    missing_fraction = as.numeric(missing_fraction),
    n_alleles = NA_integer_,
    n_unique_genotypes = NA_integer_,
    chi_square = NA_real_,
    df = NA_integer_,
    exact_p_value = NA_real_,
    status = "not_tested",
    note = ""
  )
  
  chi <- calculate_hwe_chisq(geno_vec)
  base_row$chi_square <- chi$chisq
  base_row$df <- chi$df
  
  if (length(geno_vec) == 0) {
    base_row$status <- "failed"
    base_row$note <- "empty_genotype_vector"
    return(base_row)
  }
  
  if (n_non_missing < min_n) {
    base_row$status <- "skipped"
    base_row$note <- sprintf("too_few_individuals_non_missing (N=%d < %d); %s", n_non_missing, min_n, chi$note)
    return(base_row)
  }
  
  a1 <- parsed[keep, "a1"]
  a2 <- parsed[keep, "a2"]
  genotype_strings <- paste(a1, a2, sep = "/")
  n_genotypes <- dplyr::n_distinct(genotype_strings)
  base_row$n_unique_genotypes <- as.integer(n_genotypes)
  
  if (n_genotypes < min_unique_genotypes) {
    base_row$status <- "skipped"
    base_row$note <- sprintf("insufficient_genotype_diversity (unique_genotypes=%d); %s", n_genotypes, chi$note)
    return(base_row)
  }
  
  n_alleles <- dplyr::n_distinct(c(a1, a2))
  base_row$n_alleles <- as.integer(n_alleles)
  
  if (n_alleles < 2) {
    base_row$status <- "skipped"
    base_row$note <- paste("invariant_locus", chi$note, sep = "; ")
    return(base_row)
  }
  
  error_msg <- NULL
  pval <- tryCatch({
    loc_df <- data.frame(Locus = factor(genotype_strings), stringsAsFactors = TRUE)
    loci_obj <- pegas::as.loci(loc_df)
    hw_fit <- pegas::hw.test(loci_obj, B = B)
    
    if (is.list(hw_fit) && !is.null(hw_fit$p.value)) {
      as.numeric(hw_fit$p.value[1])
    } else if (is.matrix(hw_fit) || is.data.frame(hw_fit)) {
      hw_mat <- as.matrix(hw_fit)
      cn <- tolower(colnames(hw_mat))
      pick <- which(cn %in% c("p.value", "pvalue", "p", "pr(>chi)", "pr(prob)") | grepl("p", cn))
      if (length(pick) == 0) suppressWarnings(as.numeric(hw_mat[1, ncol(hw_mat)])) else suppressWarnings(as.numeric(hw_mat[1, pick[1]]))
    } else {
      suppressWarnings(as.numeric(hw_fit[1]))
    }
  }, error = function(e) {
    error_msg <<- conditionMessage(e)
    NA_real_
  })
  
  if (!is.finite(pval) || is.na(pval)) {
    base_row$status <- "failed"
    base_row$note <- if (!is.null(error_msg) && nzchar(error_msg)) {
      paste0("exact_hwe_test_failed: ", error_msg, "; ", chi$note)
    } else {
      paste0("exact_hwe_test_failed; ", chi$note)
    }
    return(base_row)
  }
  
  base_row$exact_p_value <- pval
  base_row$status <- "ok"
  base_row$note <- paste("exact_or_monte_carlo_p_value_from_pegas_hw.test", chi$note, sep = "; ")
  base_row
}

add_bonferroni <- function(df, p_col = "exact_p_value", alpha = ALPHA) {
  p <- suppressWarnings(as.numeric(df[[p_col]]))
  ok <- is.finite(p) & !is.na(p)
  p_adj <- rep(NA_real_, length(p))
  if (sum(ok) > 0) p_adj[ok] <- p.adjust(p[ok], method = "bonferroni")
  df$bonferroni_p_value <- p_adj
  df$significant_bonferroni <- !is.na(p_adj) & p_adj < alpha
  df
}

run_hwe_scope <- function(gdf, loci, scope_name, site_label) {
  rows <- vector("list", length(loci))
  for (i in seq_along(loci)) {
    loc <- loci[i]
    test <- run_single_locus_hwe(gdf[[loc]])
    rows[[i]] <- data.frame(
      Scope = scope_name,
      Site = site_label,
      Locus = loc,
      N = as.integer(test$N),
      missing_fraction = as.numeric(test$missing_fraction),
      n_alleles = as.integer(test$n_alleles),
      n_unique_genotypes = as.integer(test$n_unique_genotypes),
      chi_square = as.numeric(test$chi_square),
      df = as.integer(test$df),
      exact_p_value = as.numeric(test$exact_p_value),
      status = as.character(test$status),
      note = as.character(test$note),
      stringsAsFactors = FALSE
    )
  }
  dplyr::bind_rows(rows)
}

format_hwe_wide <- function(hwe_long, site_levels) {
  hwe_for_wide <- hwe_long %>%
    mutate(
      Appendix_scope = ifelse(Scope == "Overall", "All", Site),
      chi_square = round(chi_square, 3),
      exact_p_value = round(exact_p_value, 6),
      bonferroni_p_value = round(bonferroni_p_value, 6)
    ) %>%
    select(
      Locus,
      Appendix_scope,
      N,
      chi_square,
      df,
      exact_p_value,
      bonferroni_p_value,
      significant_bonferroni,
      status
    )
  
  wide <- hwe_for_wide %>%
    pivot_wider(
      names_from = Appendix_scope,
      values_from = c(N, chi_square, df, exact_p_value, bonferroni_p_value, significant_bonferroni, status),
      names_glue = "{Appendix_scope}_{.metric}"
    )
  
  scope_levels <- c("All", site_levels)
  metric_levels <- c("N", "chi_square", "df", "exact_p_value", "bonferroni_p_value", "significant_bonferroni", "status")
  wanted <- c("Locus", as.vector(t(outer(scope_levels, metric_levels, paste, sep = "_"))))
  wanted <- wanted[wanted %in% names(wide)]
  
  wide %>%
    select(all_of(wanted)) %>%
    arrange(Locus)
}


format_hwe_p_value <- function(p) {
  ifelse(
    is.na(p) | !is.finite(p),
    "--",
    ifelse(p < 0.001, "<0.001", format(round(p, 3), nsmall = 3, trim = TRUE, scientific = FALSE))
  )
}

format_hwe_stat_cell <- function(chi_square, df, exact_p_value) {
  chi_txt <- ifelse(is.na(chi_square) | !is.finite(chi_square), "--", format(round(chi_square, 2), nsmall = 2, trim = TRUE, scientific = FALSE))
  df_txt <- ifelse(is.na(df), "--", as.character(df))
  p_txt <- format_hwe_p_value(exact_p_value)
  paste0("X2=", chi_txt, "; df=", df_txt, "; p=", p_txt)
}

build_hwe_display_table <- function(hwe_long, site_levels) {
  scope_levels <- c("All", site_levels)
  
  display_long <- hwe_long %>%
    mutate(
      Appendix_scope = ifelse(Scope == "Overall", "All", Site),
      Display = format_hwe_stat_cell(chi_square, df, exact_p_value),
      significant_bonferroni = !is.na(significant_bonferroni) & significant_bonferroni
    ) %>%
    filter(Appendix_scope %in% scope_levels) %>%
    select(Locus, Appendix_scope, Display, significant_bonferroni)
  
  display_table <- display_long %>%
    select(Locus, Appendix_scope, Display) %>%
    pivot_wider(names_from = Appendix_scope, values_from = Display) %>%
    select(any_of(c("Locus", scope_levels))) %>%
    arrange(Locus)
  
  significant_wide <- display_long %>%
    select(Locus, Appendix_scope, significant_bonferroni) %>%
    pivot_wider(names_from = Appendix_scope, values_from = significant_bonferroni, values_fill = FALSE) %>%
    right_join(display_table %>% select(Locus), by = "Locus") %>%
    select(any_of(c("Locus", scope_levels))) %>%
    arrange(match(Locus, display_table$Locus))
  
  significant_matrix <- as.matrix(significant_wide[, names(display_table), drop = FALSE])
  significant_matrix[is.na(significant_matrix)] <- FALSE
  significant_matrix[, "Locus"] <- FALSE
  mode(significant_matrix) <- "logical"
  
  list(table = display_table, significant = significant_matrix)
}

# ----------------------------
# Allele-frequency helpers
# ----------------------------
allele_counts_for_locus <- function(geno_vec, locus, site_vec, site_levels) {
  parsed <- parse_locus_genotypes(geno_vec)
  keep <- !(is.na(parsed[, "a1"]) | is.na(parsed[, "a2"]))
  
  if (!any(keep)) {
    return(data.frame(
      Locus = character(0), Allele = character(0), Site = character(0), allele_count = integer(0), allele_total = integer(0),
      stringsAsFactors = FALSE
    ))
  }
  
  long <- data.frame(
    Site = rep(site_vec[keep], each = 2),
    Allele = c(rbind(parsed[keep, "a1"], parsed[keep, "a2"])),
    stringsAsFactors = FALSE
  )
  
  alleles <- sort(unique(long$Allele))
  grid <- expand.grid(Site = site_levels, Allele = alleles, stringsAsFactors = FALSE)
  counts <- long %>% count(Site, Allele, name = "allele_count")
  totals <- long %>% count(Site, name = "allele_total")
  
  grid %>%
    left_join(counts, by = c("Site", "Allele")) %>%
    left_join(totals, by = "Site") %>%
    mutate(
      Locus = locus,
      allele_count = ifelse(is.na(allele_count), 0L, as.integer(allele_count)),
      allele_total = ifelse(is.na(allele_total), 0L, as.integer(allele_total))
    ) %>%
    select(Locus, Allele, Site, allele_count, allele_total)
}

build_allele_frequency_table <- function(gdf, loci, site_vec, site_levels) {
  counts <- dplyr::bind_rows(lapply(loci, function(loc) {
    allele_counts_for_locus(gdf[[loc]], loc, site_vec, site_levels)
  }))
  
  if (nrow(counts) == 0) {
    stop("[08_appendix_tables] No non-missing allele data were available to compute allele frequencies.", call. = FALSE)
  }
  
  freq_by_site <- counts %>%
    mutate(Frequency = ifelse(allele_total > 0, allele_count / allele_total, NA_real_))
  
  overall <- freq_by_site %>%
    group_by(Locus, Allele) %>%
    summarise(
      All_W = ifelse(sum(allele_total, na.rm = TRUE) > 0, sum(allele_count, na.rm = TRUE) / sum(allele_total, na.rm = TRUE), NA_real_),
      All_NW = mean(Frequency, na.rm = TRUE),
      .groups = "drop"
    )
  
  wide <- freq_by_site %>%
    mutate(Frequency = round(Frequency, FREQ_DIGITS)) %>%
    select(Locus, Allele, Site, Frequency) %>%
    pivot_wider(names_from = Site, values_from = Frequency) %>%
    left_join(overall, by = c("Locus", "Allele")) %>%
    mutate(
      All_W = round(All_W, FREQ_DIGITS),
      All_NW = round(All_NW, FREQ_DIGITS)
    )
  
  wanted <- c("Locus", "Allele", site_levels, "All_W", "All_NW")
  wide %>%
    select(any_of(wanted)) %>%
    arrange(Locus, suppressWarnings(as.numeric(Allele)), Allele)
}

# ----------------------------
# Basic DOCX writer (no officer/flextable dependency required)
# ----------------------------
xml_escape <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub('"', "&quot;", x, fixed = TRUE)
  x <- gsub("'", "&apos;", x, fixed = TRUE)
  x
}

word_cell_xml <- function(value, bold = FALSE, color = NULL, width = 1800L, font_size = 18L, shade = NULL) {
  run_props <- paste0(
    if (bold) "<w:b/>" else "",
    if (!is.null(color) && nzchar(color)) paste0("<w:color w:val=\"", color, "\"/>") else "",
    if (!is.null(font_size)) paste0("<w:sz w:val=\"", as.integer(font_size), "\"/>") else ""
  )
  run_props <- if (nzchar(run_props)) paste0("<w:rPr>", run_props, "</w:rPr>") else ""
  shade_xml <- if (!is.null(shade) && nzchar(shade)) paste0("<w:shd w:fill=\"", shade, "\"/>") else ""
  
  paste0(
    "<w:tc><w:tcPr><w:tcW w:w=\"", as.integer(width), "\" w:type=\"dxa\"/>", shade_xml,
    "<w:tcMar><w:top w:w=\"40\" w:type=\"dxa\"/><w:left w:w=\"40\" w:type=\"dxa\"/><w:bottom w:w=\"40\" w:type=\"dxa\"/><w:right w:w=\"40\" w:type=\"dxa\"/></w:tcMar>",
    "</w:tcPr><w:p><w:pPr><w:spacing w:before=\"0\" w:after=\"0\"/></w:pPr><w:r>",
    run_props, "<w:t>", xml_escape(value), "</w:t></w:r></w:p></w:tc>"
  )
}

word_paragraph_xml <- function(value, bold = FALSE, font_size = 20L) {
  run_props <- paste0(
    if (bold) "<w:b/>" else "",
    if (!is.null(font_size)) paste0("<w:sz w:val=\"", as.integer(font_size), "\"/>") else ""
  )
  run_props <- if (nzchar(run_props)) paste0("<w:rPr>", run_props, "</w:rPr>") else ""
  paste0("<w:p><w:pPr><w:spacing w:before=\"0\" w:after=\"80\"/></w:pPr><w:r>", run_props, "<w:t>", xml_escape(value), "</w:t></w:r></w:p>")
}

format_docx_df <- function(df) {
  out <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
  for (nm in names(out)) {
    if (is.numeric(out[[nm]])) {
      out[[nm]] <- ifelse(is.na(out[[nm]]), "", format(out[[nm]], trim = TRUE, scientific = FALSE))
    }
    if (is.logical(out[[nm]])) out[[nm]] <- ifelse(out[[nm]], "Yes", "No")
  }
  out
}

write_basic_docx <- function(df,
                             path,
                             title,
                             note,
                             bold_matrix = NULL,
                             red_matrix = NULL,
                             cell_width = 1800L,
                             font_size = 18L,
                             header_font_size = 18L,
                             margins = 720L) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  df <- format_docx_df(df)
  df[] <- lapply(df, as.character)
  
  if (is.null(bold_matrix)) bold_matrix <- matrix(FALSE, nrow = nrow(df), ncol = ncol(df), dimnames = list(NULL, names(df)))
  if (is.null(red_matrix)) red_matrix <- matrix(FALSE, nrow = nrow(df), ncol = ncol(df), dimnames = list(NULL, names(df)))
  bold_matrix <- as.matrix(bold_matrix[, names(df), drop = FALSE])
  red_matrix <- as.matrix(red_matrix[, names(df), drop = FALSE])
  
  header <- paste0(
    "<w:tr>",
    paste(vapply(names(df), word_cell_xml, character(1), bold = TRUE, width = cell_width, font_size = header_font_size, shade = "D9EAF7"), collapse = ""),
    "</w:tr>"
  )
  
  rows <- vapply(seq_len(nrow(df)), function(i) {
    cells <- vapply(seq_along(df), function(j) {
      word_cell_xml(
        df[[j]][i],
        bold = isTRUE(bold_matrix[i, j]),
        color = if (isTRUE(red_matrix[i, j])) "C00000" else NULL,
        width = cell_width,
        font_size = font_size
      )
    }, character(1))
    paste0("<w:tr>", paste(cells, collapse = ""), "</w:tr>")
  }, character(1))
  
  tmp_dir <- tempfile("appendix_docx_")
  dir.create(file.path(tmp_dir, "_rels"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(tmp_dir, "word", "_rels"), recursive = TRUE, showWarnings = FALSE)
  
  document_xml <- paste0(
    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
    '<w:document xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><w:body>',
    word_paragraph_xml(title, bold = TRUE, font_size = 22L),
    word_paragraph_xml(note, font_size = 18L),
    '<w:tbl><w:tblPr><w:tblStyle w:val="TableGrid"/><w:tblW w:w="0" w:type="auto"/><w:tblLayout w:type="fixed"/><w:tblBorders>',
    '<w:top w:val="single" w:sz="4" w:space="0" w:color="auto"/><w:left w:val="single" w:sz="4" w:space="0" w:color="auto"/>',
    '<w:bottom w:val="single" w:sz="4" w:space="0" w:color="auto"/><w:right w:val="single" w:sz="4" w:space="0" w:color="auto"/>',
    '<w:insideH w:val="single" w:sz="4" w:space="0" w:color="auto"/><w:insideV w:val="single" w:sz="4" w:space="0" w:color="auto"/>',
    '</w:tblBorders><w:tblCellMar><w:top w:w="40" w:type="dxa"/><w:left w:w="40" w:type="dxa"/><w:bottom w:w="40" w:type="dxa"/><w:right w:w="40" w:type="dxa"/></w:tblCellMar></w:tblPr>',
    header, paste(rows, collapse = ""), '</w:tbl>',
    '<w:sectPr><w:pgSz w:w="15840" w:h="12240" w:orient="landscape"/><w:pgMar w:top="', margins, '" w:right="', margins, '" w:bottom="', margins, '" w:left="', margins, '" w:header="360" w:footer="360" w:gutter="0"/></w:sectPr>',
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
  if (!file.exists(path)) stop("Failed to create DOCX: ", path, call. = FALSE)
  message_appendix("Saved Word DOCX: ", path)
  invisible(path)
}

# ----------------------------
# Main workflow
# ----------------------------
search_existing_appendix_tables()

source("scripts/_load_objects.R")

if (!inherits(gi_mll, "genind")) {
  stop("[08_appendix_tables] gi_mll must be a genind object. Run scripts/00_master_pipeline.R first.", call. = FALSE)
}

validate_columns(df_ids_mll, c("ind_id", "Site"), df_name = "[08_appendix_tables] df_ids_mll")
if (!all(adegenet::indNames(gi_mll) == df_ids_mll$ind_id)) {
  stop("[08_appendix_tables] gi_mll and df_ids_mll are not aligned.", call. = FALSE)
}

site_lookup <- build_site_label_map(df_ids_mll)
site_labels <- apply_site_labels(df_ids_mll$Site, site_lookup)
site_levels <- site_sort_levels(site_labels, site_lookup)

adegenet::pop(gi_mll) <- factor(site_labels, levels = site_levels)
df_ids_mll_appendix <- df_ids_mll %>%
  mutate(Site_original = as.character(Site), Site = site_labels)

gdf <- adegenet::genind2df(gi_mll, sep = "/")
loci <- setdiff(names(gdf), "pop")
if (length(loci) == 0) stop("[08_appendix_tables] No loci detected in gi_mll.", call. = FALSE)

message_appendix("Dataset used: gi_mll (MLL clone-corrected genind object).")
message_appendix("Individuals included: ", adegenet::nInd(gi_mll))
message_appendix("Loci included: ", adegenet::nLoc(gi_mll))
message_appendix("Sampling sites included: ", paste(site_levels, collapse = ", "))
message_appendix("Output directory: ", APPENDIX_DIR)

# HWE table: overall plus every site, one row per locus in wide appendix format.
hwe_overall <- run_hwe_scope(gdf = gdf, loci = loci, scope_name = "Overall", site_label = "All")
hwe_by_site <- dplyr::bind_rows(lapply(site_levels, function(site) {
  idx <- which(site_labels == site)
  run_hwe_scope(gdf = gdf[idx, , drop = FALSE], loci = loci, scope_name = "Site", site_label = site)
}))

hwe_long <- dplyr::bind_rows(hwe_overall, hwe_by_site) %>%
  add_bonferroni(p_col = "exact_p_value", alpha = ALPHA) %>%
  arrange(factor(Scope, levels = c("Overall", "Site")), factor(Site, levels = c("All", site_levels)), Locus)

hwe_wide <- format_hwe_wide(hwe_long, site_levels)
hwe_display <- build_hwe_display_table(hwe_long, site_levels)
hwe_display_table <- hwe_display$table
hwe_display_significant <- hwe_display$significant

hwe_note <- paste(
  "HWE exact p-values are from pegas::hw.test with", HWE_MONTE_CARLO_REPS,
  "Monte Carlo replicates. Chi-square values are approximate Pearson statistics computed from observed and expected genotype counts under HWE; exact/permutation p-values, not chi-square p-values, should be used for inference. Bonferroni significance is evaluated across all valid overall and site-by-locus tests at alpha =", ALPHA, "."
)

hwe_display_note <- paste(
  "Compact thesis-display HWE table: each dataset/site cell reports X2, df, and the exact/Monte Carlo p-value from pegas::hw.test.",
  "Exact/permutation p-values, not chi-square p-values, are used for inference.",
  "Bold red cells are significant after Bonferroni correction across all valid overall and site-by-locus tests at alpha =",
  ALPHA,
  "."
)

hwe_csv <- file.path(APPENDIX_DIR, "appendix_hwe_wide.csv")
hwe_xlsx <- file.path(APPENDIX_DIR, "appendix_hwe_wide.xlsx")
hwe_docx <- file.path(APPENDIX_DIR, "appendix_hwe_wide.docx")
hwe_long_csv <- file.path(APPENDIX_DIR, "appendix_hwe_long_source.csv")
hwe_display_csv <- file.path(APPENDIX_DIR, "appendix_hwe_table_C1_display.csv")
hwe_display_docx <- file.path(APPENDIX_DIR, "appendix_hwe_table_C1_display.docx")

write_csv_out(hwe_wide, hwe_csv)
write_csv_out(hwe_long, hwe_long_csv)
write_csv_out(hwe_display_table, hwe_display_csv)
write_xlsx_out(list(HWE_appendix_wide = hwe_wide, HWE_long_source = hwe_long), hwe_xlsx)
write_basic_docx(hwe_wide, hwe_docx, "Appendix C. Hardy-Weinberg equilibrium tests by locus and sampling site", hwe_note)
write_basic_docx(
  hwe_display_table,
  hwe_display_docx,
  "Table C.1. Hardy-Weinberg equilibrium tests by locus and sampling site",
  hwe_display_note,
  bold_matrix = hwe_display_significant,
  red_matrix = hwe_display_significant,
  cell_width = 900L,
  font_size = 14L,
  header_font_size = 14L,
  margins = 360L
)

# Allele-frequency table: clone-corrected site frequencies plus All_W and All_NW.
allele_freq_table <- build_allele_frequency_table(gdf = gdf, loci = loci, site_vec = site_labels, site_levels = site_levels)

allele_note <- paste(
  "Allele frequencies were calculated from the MLL clone-corrected gi_mll dataset.",
  "All_W is the weighted overall frequency across all clone-corrected individuals.",
  "All_NW is the non-weighted mean of site-level frequencies across sampling sites.",
  "Frequencies are rounded to", FREQ_DIGITS, "decimals."
)

allele_csv <- file.path(APPENDIX_DIR, "appendix_allele_frequencies.csv")
allele_xlsx <- file.path(APPENDIX_DIR, "appendix_allele_frequencies.xlsx")
allele_docx <- file.path(APPENDIX_DIR, "appendix_allele_frequencies.docx")

write_csv_out(allele_freq_table, allele_csv)
write_xlsx_out(allele_freq_table, allele_xlsx)
write_basic_docx(allele_freq_table, allele_docx, "Appendix D. Microsatellite allele frequencies by locus and sampling site", allele_note)

# Provenance/readme table for reproducibility.
provenance <- data.frame(
  item = c(
    "dataset",
    "clone_correction",
    "individuals",
    "loci",
    "site_labels",
    "hwe_test",
    "hwe_chi_square_note",
    "bonferroni_scope",
    "allele_frequency_weighted_overall",
    "allele_frequency_nonweighted_overall",
    "output_directory"
  ),
  value = c(
    "gi_mll",
    "MLL clone-corrected dataset loaded from outputs/v1/objects/gi_mll.rds",
    as.character(adegenet::nInd(gi_mll)),
    as.character(adegenet::nLoc(gi_mll)),
    paste(site_levels, collapse = ", "),
    paste0("pegas::hw.test exact/Monte Carlo test, B = ", HWE_MONTE_CARLO_REPS),
    "Chi-square columns are approximate Pearson statistics from observed/expected genotype counts; exact p-values are used for significance.",
    "Across all valid overall and site-by-locus HWE tests in this appendix table.",
    "All_W = allele copies for an allele across all individuals divided by total allele copies at that locus.",
    "All_NW = arithmetic mean of site-level allele frequencies across sites.",
    APPENDIX_DIR
  ),
  stringsAsFactors = FALSE
)

provenance_csv <- file.path(APPENDIX_DIR, "appendix_tables_provenance.csv")
write_csv_out(provenance, provenance_csv)

message_appendix("Completed appendix table generation.")
message_appendix("HWE detailed outputs: ", hwe_csv, "; ", hwe_xlsx, "; ", hwe_docx)
message_appendix("HWE compact display outputs: ", hwe_display_csv, "; ", hwe_display_docx)
message_appendix("Allele-frequency outputs: ", allele_csv, "; ", allele_xlsx, "; ", allele_docx)