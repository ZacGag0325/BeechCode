#!/usr/bin/env Rscript
# Appendix M: diagnostics for selection of K in final HEG_ZG_run2 STRUCTURE runs.

# ---- Configuration -----------------------------------------------------------
PROJECT_ROOT <- path.expand("~/Desktop/BeechCode")
VERBOSE_FILE_LISTING <- FALSE
VERBOSE_PARSING <- FALSE

EXPECTED_K <- 1:12
EXPECTED_REPLICATES <- 20L

SITE_LEVELS <- c(
  "S1", "S2", "S3", "S4", "S5", "S6",
  "N1", "N2", "N3", "N4", "N5", "N6"
)

SITE_KEY <- c(
  AMC = "S1", ALB = "S2", IKJ = "S3", IKO = "S4",
  LGG = "S5", LGR = "S6", ML1 = "N1", ML2 = "N2",
  ML3 = "N3", CPF = "N4", PLI = "N5", LDF = "N6"
)

required <- c(
  "dplyr", "tidyr", "readr", "stringr", "purrr", "tibble",
  "ggplot2", "fs", "officer", "flextable"
)

missing <- required[
  !vapply(required, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing)) {
  stop(
    "[appendix_M] Missing required R package(s): ",
    paste(missing, collapse = ", "),
    ". Install them before running this script.",
    call. = FALSE
  )
}

invisible(lapply(required, library, character.only = TRUE))

if (!dir.exists(PROJECT_ROOT)) {
  stop(
    "[appendix_M] Project root does not exist: ",
    PROJECT_ROOT,
    call. = FALSE
  )
}

abort_m <- function(...) {
  stop("[appendix_M] ", paste0(..., collapse = ""), call. = FALSE)
}

out_tables <- fs::path(PROJECT_ROOT, "outputs/tables/supplementary")
out_figures <- fs::path(PROJECT_ROOT, "outputs/figures/supplementary")
out_derived <- fs::path(PROJECT_ROOT, "outputs/derived")

fs::dir_create(
  c(out_tables, out_figures, out_derived),
  recurse = TRUE
)

num_re <- "[-+]?(?:[0-9]*\\.?[0-9]+|[0-9]+\\.?)(?:[eE][-+]?[0-9]+)?"

normalise_name <- function(x) {
  stringr::str_replace_all(tolower(x), "[^a-z0-9]", "")
}

number_after <- function(x, label) {
  z <- stringr::str_match(
    x,
    paste0(
      "(?i)",
      label,
      ".{0,80}?(?:[:=]|\\s)\\s*(",
      num_re,
      ")"
    )
  )[, 2]
  
  z <- z[!is.na(z)]
  
  if (length(z)) as.numeric(z[[1]]) else NA_real_
}

integer_after <- function(x, label) {
  z <- number_after(x, label)
  ifelse(is.finite(z), as.integer(z), NA_integer_)
}

save_plot <- function(p, stem, width = 6.5, height = 4.2) {
  ggplot2::ggsave(
    fs::path(out_figures, paste0(stem, ".png")),
    p,
    width = width,
    height = height,
    units = "in",
    dpi = 600
  )
  
  ggplot2::ggsave(
    fs::path(out_figures, paste0(stem, ".pdf")),
    p,
    width = width,
    height = height,
    units = "in"
  )
}

theme_m <- ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major.x = ggplot2::element_blank()
  )

# ---- Find only final-run inputs ---------------------------------------------
structure_dir <- fs::path(PROJECT_ROOT, "data", "structure")
final_dir <- fs::path(structure_dir, "HEG_ZG_run2")
medk_file <- fs::path(structure_dir, "1769700564.MedK.0.5.tsv")
sum_file <- fs::path(structure_dir, "1769700564.sum.tsv")
delta_file <- fs::path(structure_dir, "1769700564.DeltaK.tsv")
zip_file <- fs::path(structure_dir, "HEG_ZG_run2_logs.zip")

search_roots <- character()

if (dir.exists(final_dir)) {
  search_roots <- c(search_roots, final_dir)
}

if (file.exists(zip_file)) {
  extracted <- fs::path(
    out_derived,
    "HEG_ZG_run2_logs_extracted"
  )
  
  fs::dir_create(extracted, recurse = TRUE)
  
  utils::unzip(
    zip_file,
    exdir = extracted,
    overwrite = FALSE
  )
  
  search_roots <- c(search_roots, extracted)
  
  if (VERBOSE_PARSING) {
    message("[appendix_M] Extracted read-only archive to: ", extracted)
  }
}

if (!length(search_roots)) {
  # Deliberately exclude old first-run, Q_Files, and literature-review material.
  all_dirs <- fs::dir_ls(
    structure_dir,
    recurse = TRUE,
    type = "directory",
    fail = FALSE
  )
  
  search_roots <- all_dirs[
    stringr::str_detect(
      all_dirs,
      "(?i)HEG[_-]?ZG[_-]?run2"
    )
  ]
}

if (!length(search_roots)) {
  abort_m("No final HEG_ZG_run2 STRUCTURE outputs were found.")
}

# ---- Parse and validate final STRUCTURE outputs -----------------------------
parse_ids <- function(path) {
  s <- paste(
    fs::path_file(path),
    fs::path_rel(path, start = PROJECT_ROOT),
    sep = " "
  )
  
  # Do not interpret the final-analysis name HEG_ZG_run2 as replicate 2.
  s_rep <- stringr::str_replace_all(
    s,
    "(?i)HEG[_ .-]?ZG[_ .-]?run2",
    ""
  )
  
  k_matches <- stringr::str_match(
    s,
    "(?i)(?:^|[_ .-])K[_ .-]?([0-9]{1,2})(?:[_ .-]|$)|K=([0-9]{1,2})"
  )
  
  k_values <- k_matches[!is.na(k_matches)]
  k <- if (length(k_values) >= 2) k_values[[2]] else NA_character_
  
  r <- stringr::str_match(
    s_rep,
    "(?i)rep(?:licate)?[_ .-]?([0-9]+)"
  )[, 2]
  
  if (is.na(r)) {
    r <- stringr::str_match(
      s_rep,
      "(?i)(?:^|[_ .-])run[_ .-]?([0-9]+)(?:[_ .-]|$)"
    )[, 2]
  }
  
  if (is.na(r)) {
    r <- stringr::str_match(
      s_rep,
      "(?i)K[_ .-]?[0-9]+[_ .-]+([0-9]+)(?:[_ .-]|$)"
    )[, 2]
  }
  
  list(
    K = suppressWarnings(as.integer(k)),
    replicate = r
  )
}

parse_structure <- function(path) {
  x <- tryCatch(
    readLines(path, warn = FALSE, encoding = "UTF-8"),
    error = function(e) character()
  )
  
  if (!length(x)) {
    return(NULL)
  }
  
  id <- parse_ids(path)
  
  lnp <- number_after(
    x,
    paste0(
      "(?:",
      "estimated\\s+ln\\s+prob(?:ability)?\\s+of\\s+data|",
      "lnp\\s*\\(\\s*d\\s*\\)|",
      "log\\s+probability\\s+of\\s+data",
      ")"
    )
  )
  
  if (!is.finite(lnp)) {
    return(NULL)
  }
  
  tibble::tibble(
    K = id$K,
    replicate = id$replicate,
    LnP = lnp,
    mean_ln_likelihood = number_after(
      x,
      "mean\\s+value\\s+of\\s+ln\\s+likelihood"
    ),
    variance_ln_likelihood = number_after(
      x,
      "variance\\s+of\\s+ln\\s+likelihood"
    ),
    n_individuals = integer_after(
      x,
      "(?:number\\s+of\\s+individuals|numind|individuals)"
    ),
    n_loci = integer_after(
      x,
      "(?:number\\s+of\\s+loci|numloci|loci)"
    ),
    burnin = integer_after(
      x,
      "(?:burnin|burn-in)"
    ),
    mcmc_iterations = integer_after(
      x,
      "(?:numreps|mcmc(?:\\s+iterations)?)"
    ),
    path = path
  )
}

output_files <- unique(unlist(
  lapply(search_roots, function(x) {
    fs::dir_ls(
      x,
      recurse = TRUE,
      type = "file",
      fail = FALSE,
      regexp = "(?i)(?:_f$|\\.out$|\\.txt$|\\.log$|result|structure)"
    )
  }),
  use.names = FALSE
))

if (VERBOSE_FILE_LISTING) {
  for (x in output_files) {
    message("[appendix_M] Candidate: ", x)
  }
}

parsed <- purrr::compact(
  purrr::map(output_files, parse_structure)
)

replicates <- dplyr::bind_rows(parsed) |>
  dplyr::filter(
    !is.na(K),
    K %in% EXPECTED_K,
    !is.na(replicate)
  )

if (!nrow(replicates)) {
  abort_m(
    "No parseable K/replicate/LnP(D) records were found in final HEG_ZG_run2 outputs."
  )
}

dupes <- replicates |>
  dplyr::count(K, replicate) |>
  dplyr::filter(n > 1)

if (nrow(dupes)) {
  abort_m(
    "Duplicate K-replicate combination(s): ",
    paste(
      paste0(
        "K=", dupes$K,
        ", replicate=", dupes$replicate
      ),
      collapse = "; "
    )
  )
}

counts <- replicates |>
  dplyr::count(K, name = "n") |>
  tidyr::complete(
    K = EXPECTED_K,
    fill = list(n = 0L)
  )

bad <- dplyr::filter(
  counts,
  n != EXPECTED_REPLICATES
)

if (nrow(bad)) {
  for (k in bad$K) {
    message(
      "[appendix_M] K=", k,
      " files: ",
      paste(
        replicates$path[replicates$K == k],
        collapse = "; "
      )
    )
  }
  
  abort_m(
    "K = ",
    bad$K[[1]],
    " contains only ",
    bad$n[[1]],
    " valid replicates; expected 20."
  )
}

if (any(!is.finite(replicates$LnP))) {
  abort_m("Non-finite LnP(D) values were detected.")
}

n_individuals <- unique(
  stats::na.omit(replicates$n_individuals)
)

if (length(n_individuals) != 1L) {
  abort_m(
    "Individual counts are absent or inconsistent across final STRUCTURE files."
  )
}

n_individuals <- n_individuals[[1]]

# ---- Evanno diagnostics ------------------------------------------------------
evanno <- replicates |>
  dplyr::group_by(K) |>
  dplyr::summarise(
    n_replicates = dplyr::n(),
    mean_LnP_K = mean(LnP),
    sd_LnP_K = stats::sd(LnP),
    se_LnP_K = sd_LnP_K / sqrt(n_replicates),
    min_LnP_K = min(LnP),
    max_LnP_K = max(LnP),
    .groups = "drop"
  ) |>
  dplyr::arrange(K) |>
  dplyr::mutate(
    L_prime_K = dplyr::lead(mean_LnP_K) - mean_LnP_K,
    L_double_prime_K = abs(
      dplyr::lead(mean_LnP_K) -
        2 * mean_LnP_K +
        dplyr::lag(mean_LnP_K)
    ),
    DeltaK = dplyr::if_else(
      is.finite(L_double_prime_K) & sd_LnP_K > 0,
      L_double_prime_K / sd_LnP_K,
      NA_real_
    )
  )

max_delta <- evanno |>
  dplyr::filter(is.finite(DeltaK)) |>
  dplyr::slice_max(
    DeltaK,
    n = 1,
    with_ties = FALSE
  )

if (!nrow(max_delta)) {
  abort_m("No finite Delta K value could be calculated.")
}

evanno_reps_path <- fs::path(
  out_tables,
  "evanno_replicates.csv"
)

evanno_path <- fs::path(
  out_tables,
  "evanno_summary.csv"
)

readr::write_csv(replicates, evanno_reps_path)
readr::write_csv(evanno, evanno_path)

p_lnp <- ggplot2::ggplot(
  evanno,
  ggplot2::aes(K, mean_LnP_K)
) +
  ggplot2::geom_errorbar(
    ggplot2::aes(
      ymin = mean_LnP_K - sd_LnP_K,
      ymax = mean_LnP_K + sd_LnP_K
    ),
    width = 0.16
  ) +
  ggplot2::geom_line(linewidth = 0.45) +
  ggplot2::geom_point(size = 2) +
  ggplot2::scale_x_continuous(breaks = EXPECTED_K) +
  ggplot2::labs(
    x = "Number of genetic clusters (K)",
    y = "Mean LnP(K)"
  ) +
  theme_m

p_delta <- ggplot2::ggplot(
  evanno,
  ggplot2::aes(K, DeltaK)
) +
  ggplot2::geom_line(
    linewidth = 0.45,
    na.rm = TRUE
  ) +
  ggplot2::geom_point(
    size = 2,
    na.rm = TRUE
  ) +
  ggplot2::geom_point(
    data = max_delta,
    shape = 21,
    fill = "white",
    colour = "#9A4D3E",
    size = 3.5,
    stroke = 0.8
  ) +
  ggplot2::scale_x_continuous(breaks = EXPECTED_K) +
  ggplot2::labs(
    x = "Number of genetic clusters (K)",
    y = expression(Delta * K)
  ) +
  theme_m

save_plot(p_lnp, "appendix_M_mean_LnPK")
save_plot(p_delta, "appendix_M_evanno_deltaK")

# ---- StructureSelector Puechmaille parsing ----------------------------------
clean_table <- function(x) {
  x <- as.data.frame(
    x,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  x <- x[
    apply(
      x,
      1,
      function(z) {
        any(
          !is.na(z) &
            trimws(as.character(z)) != ""
        )
      }
    ),
    ,
    drop = FALSE
  ]
  
  if (nrow(x) > 1) {
    x <- x[
      !apply(
        x,
        1,
        function(z) {
          identical(
            normalise_name(z),
            normalise_name(names(x))
          )
        }
      ),
      ,
      drop = FALSE
    ]
  }
  
  tibble::as_tibble(
    x,
    .name_repair = "minimal"
  )
}

show_attempt <- function(label, x) {
  if (is.null(x)) {
    message("[appendix_M] ", label, ": failed")
    return()
  }
  
  message(
    "[appendix_M] ",
    label,
    ": ",
    nrow(x),
    " rows x ",
    ncol(x),
    " columns; names: ",
    paste(names(x), collapse = " | ")
  )
  
  print(utils::head(x, 10))
}

inspect_tsv <- function(path) {
  raw <- readLines(
    path,
    warn = FALSE,
    encoding = "UTF-8"
  )
  
  if (VERBOSE_PARSING) {
    message("[appendix_M] Inspecting StructureSelector TSV: ", path)
    
    message(
      paste(
        sprintf(
          "%02d: %s",
          seq_len(min(30, length(raw))),
          raw[seq_len(min(30, length(raw)))]
        ),
        collapse = "\n"
      )
    )
  }
  
  # Do not submit StructureSelector comments/metadata to readr. The only
  # rectangular data block in a MedK file starts at its K/MedMed header.
  header_i <- which(
    stringr::str_detect(
      raw,
      "^K\\tMedMed\\tMedMean\\tMaxMed\\tMaxMean\\tReps\\s*$"
    )
  )
  
  summary_i <- which(
    stringr::str_detect(
      raw,
      "^\\s*(?:MedMedK|MedMed)\\t"
    )
  )
  
  attempts <- list()
  
  if (length(header_i)) {
    end_i <- summary_i[
      summary_i > header_i[[1]]
    ][1]
    
    if (is.na(end_i)) {
      end_i <- length(raw) + 1L
    }
    
    table_text <- paste(
      raw[header_i[[1]]:(end_i - 1L)],
      collapse = "\n"
    )
    
    attempts$medk_table <- readr::read_tsv(
      I(table_text),
      show_col_types = FALSE,
      name_repair = "minimal",
      progress = FALSE
    )
  } else if (VERBOSE_PARSING) {
    message(
      "[appendix_M] No valid tabular MedK block identified in: ",
      path
    )
  }
  
  if (VERBOSE_PARSING) {
    purrr::iwalk(
      attempts,
      function(x, label) show_attempt(label, x)
    )
  }
  
  list(
    raw = raw,
    attempts = attempts
  )
}

method_name <- function(x) {
  z <- normalise_name(x)
  
  dplyr::case_when(
    z %in% c("medmed", "medmedk") ~ "MedMedK",
    z %in% c("medmean", "medmeank", "medmeak") ~ "MedMeanK",
    z %in% c("maxmed", "maxmedk") ~ "MaxMedK",
    z %in% c("maxmean", "maxmeank", "maxmeak") ~ "MaxMeanK",
    TRUE ~ NA_character_
  )
}

numeric_k <- function(x) {
  y <- suppressWarnings(
    as.numeric(
      stringr::str_extract(
        as.character(x),
        "[0-9]+(?:\\.[0-9]+)?"
      )
    )
  )
  
  as.integer(y)
}

parse_medk_summary_rows <- function(raw, path, threshold) {
  # StructureSelector MedK files contain a final two-line summary, e.g.
  # "\tMedMedK\tMedMeaK\tMaxMedK\tMaxMeaK" followed by "ALL\t7\t7\t8\t8".
  header_i <- which(
    stringr::str_detect(
      raw,
      "^\\s*(?:MedMedK|MedMed)\\t"
    )
  )
  
  all_i <- which(
    stringr::str_detect(
      raw,
      "^ALL\\t"
    )
  )
  
  if (!length(header_i) || !length(all_i)) {
    return(NULL)
  }
  
  header_i <- header_i[[length(header_i)]]
  
  all_i <- all_i[
    all_i > header_i
  ][1]
  
  if (is.na(all_i)) {
    return(NULL)
  }
  
  headers <- strsplit(
    sub("^\\t", "", raw[[header_i]]),
    "\\t"
  )[[1]]
  
  values <- strsplit(
    raw[[all_i]],
    "\\t"
  )[[1]][-1]
  
  n <- min(length(headers), length(values))
  
  ans <- tibble::tibble(
    Method = method_name(headers[seq_len(n)]),
    Selected_K = numeric_k(values[seq_len(n)])
  ) |>
    dplyr::filter(
      !is.na(Method),
      is.finite(Selected_K),
      Selected_K %in% EXPECTED_K
    ) |>
    dplyr::distinct(Method, .keep_all = TRUE)
  
  if (!nrow(ans)) {
    return(NULL)
  }
  
  list(
    results = dplyr::mutate(
      ans,
      Threshold = threshold,
      Source = paste0(
        "StructureSelector: ",
        fs::path_file(path)
      ),
      Input_file = path,
      Parsing_layout = "StructureSelector ALL summary row"
    ),
    raw = tibble::tibble(raw_line = raw)
  )
}

parse_puech_table <- function(inspected, path, threshold) {
  # Prioritize StructureSelector's explicit final ALL row.
  direct <- parse_medk_summary_rows(
    inspected$raw,
    path,
    threshold
  )
  
  if (!is.null(direct)) {
    return(direct)
  }
  
  for (tab in inspected$attempts) {
    if (is.null(tab) || !ncol(tab)) {
      next
    }
    
    tab <- clean_table(tab)
    names(tab) <- make.unique(names(tab), sep = "_")
    
    norm <- normalise_name(names(tab))
    
    raw_out <- tibble::as_tibble(
      tab,
      .name_repair = "minimal"
    )
    
    # Wide format: estimator columns contain selected K values.
    meth_cols <- method_name(names(tab))
    
    if (any(!is.na(meth_cols))) {
      rows <- purrr::map_dfr(
        which(!is.na(meth_cols)),
        function(j) {
          tibble::tibble(
            Method = meth_cols[[j]],
            Selected_K = numeric_k(tab[[j]])
          )
        }
      ) |>
        dplyr::filter(
          is.finite(Selected_K),
          Selected_K %in% EXPECTED_K
        ) |>
        dplyr::distinct(Method, .keep_all = TRUE)
      
      if (nrow(rows)) {
        return(list(
          results = dplyr::mutate(
            rows,
            Threshold = threshold,
            Source = paste0(
              "StructureSelector: ",
              fs::path_file(path)
            ),
            Input_file = path,
            Parsing_layout = "wide estimator columns"
          ),
          raw = raw_out
        ))
      }
    }
    
    # Long format: one estimator per row with a selected K/value column.
    est_i <- which(
      norm %in% c(
        "estimator",
        "method",
        "criterion",
        "measure"
      )
    )
    
    k_i <- which(
      norm %in% c(
        "k",
        "selectedk",
        "bestk",
        "value",
        "selectedclusters",
        "clusters"
      )
    )
    
    if (length(est_i) && length(k_i)) {
      rows <- tibble::tibble(
        Method = method_name(tab[[est_i[[1]]]]),
        Selected_K = numeric_k(tab[[k_i[[1]]]])
      ) |>
        dplyr::filter(
          !is.na(Method),
          is.finite(Selected_K),
          Selected_K %in% EXPECTED_K
        ) |>
        dplyr::distinct(Method, .keep_all = TRUE)
      
      if (nrow(rows)) {
        return(list(
          results = dplyr::mutate(
            rows,
            Threshold = threshold,
            Source = paste0(
              "StructureSelector: ",
              fs::path_file(path)
            ),
            Input_file = path,
            Parsing_layout = "long estimator/selected-K rows"
          ),
          raw = raw_out
        ))
      }
    }
    
    # Transposed/summary format: estimator names occur in one column and
    # selected K values occur in another column.
    for (j in seq_len(ncol(tab))) {
      for (i in seq_len(ncol(tab))) {
        if (i == j) {
          next
        }
        
        rows <- tibble::tibble(
          Method = method_name(tab[[j]]),
          Selected_K = numeric_k(tab[[i]])
        ) |>
          dplyr::filter(
            !is.na(Method),
            is.finite(Selected_K),
            Selected_K %in% EXPECTED_K
          ) |>
          dplyr::distinct(Method, .keep_all = TRUE)
        
        if (nrow(rows) >= 2) {
          return(list(
            results = dplyr::mutate(
              rows,
              Threshold = threshold,
              Source = paste0(
                "StructureSelector: ",
                fs::path_file(path)
              ),
              Input_file = path,
              Parsing_layout = "summary rows/columns"
            ),
            raw = raw_out
          ))
        }
      }
    }
  }
  
  NULL
}

threshold_from <- function(path, raw, default = 0.5) {
  x <- paste(
    c(fs::path_file(path), raw),
    collapse = " "
  )
  
  z <- number_after(
    x,
    "(?:threshold|membership\\s+cutoff|medk)"
  )
  
  if (is.finite(z) && z >= 0 && z <= 1) {
    z
  } else if (
    stringr::str_detect(
      fs::path_file(path),
      "(?i)MedK[._-]0[._-]5"
    )
  ) {
    0.5
  } else {
    default
  }
}

puech_inputs <- c(
  medk_file,
  sum_file
)[
  file.exists(c(medk_file, sum_file))
]

inspections <- lapply(
  puech_inputs,
  inspect_tsv
)

parsed_puech <- purrr::map2(
  inspections,
  puech_inputs,
  function(x, p) {
    parse_puech_table(
      x,
      p,
      threshold_from(p, x$raw)
    )
  }
)

good <- which(
  !vapply(parsed_puech, is.null, logical(1))
)

if (!length(good)) {
  abort_m(
    "Puechmaille estimators could not be parsed from existing StructureSelector results."
  )
}

chosen <- parsed_puech[[good[[1]]]]

if (length(good) > 1) {
  other <- parsed_puech[[good[[2]]]]$results
  
  shared <- intersect(
    chosen$results$Method,
    other$Method
  )
  
  if (
    length(shared) &&
    any(
      chosen$results$Selected_K[
        match(shared, chosen$results$Method)
      ] != other$Selected_K[
        match(shared, other$Method)
      ]
    )
  ) {
    warning(
      "[appendix_M] StructureSelector TSV files disagree for one or more shared estimators; ",
      "MedK.0.5.tsv is preferred."
    )
  }
}

puech <- chosen$results

raw_path <- fs::path(
  out_tables,
  "puechmaille_parsed_raw_table.csv"
)

puech_path <- fs::path(
  out_tables,
  "puechmaille_summary.csv"
)

readr::write_csv(chosen$raw, raw_path)
readr::write_csv(puech, puech_path)

message(
  "[appendix_M] Parsed existing StructureSelector Puechmaille results from: ",
  puech$Input_file[[1]]
)

message(
  "[appendix_M] Puechmaille threshold: ",
  puech$Threshold[[1]]
)

for (m in c("MedMedK", "MedMeanK", "MaxMedK", "MaxMeanK")) {
  message(
    "[appendix_M] ",
    m,
    " selected K: ",
    if (m %in% puech$Method) {
      puech$Selected_K[match(m, puech$Method)]
    } else {
      "not available"
    }
  )
}

# Existing TSVs report selected K per estimator, so a horizontal point plot is
# scientifically preferable to falsely representing selected values as curves.
figure_type <- "horizontal selected-K point plot"

puech$Method_display <- factor(
  puech$Method,
  levels = rev(c(
    "MedMedK",
    "MedMeanK",
    "MaxMedK",
    "MaxMeanK"
  ))
)

p_puech <- ggplot2::ggplot(
  puech,
  ggplot2::aes(
    Selected_K,
    Method_display,
    colour = Method
  )
) +
  ggplot2::geom_segment(
    ggplot2::aes(
      x = 0,
      xend = Selected_K,
      y = Method_display,
      yend = Method_display
    ),
    colour = "grey70"
  ) +
  ggplot2::geom_point(size = 3) +
  ggplot2::scale_x_continuous(
    breaks = EXPECTED_K,
    limits = c(0, 12)
  ) +
  ggplot2::scale_colour_manual(
    values = c(
      MedMedK = "#4C78A8",
      MedMeanK = "#72A5A1",
      MaxMedK = "#9A7D4F",
      MaxMeanK = "#805D93"
    ),
    drop = FALSE
  ) +
  ggplot2::labs(
    x = "Selected number of genetic clusters (K)",
    y = "Estimator",
    colour = "Estimator"
  ) +
  theme_m +
  ggplot2::theme(legend.position = "none")

save_plot(
  p_puech,
  "appendix_M_puechmaille_estimators",
  width = 6.5,
  height = 3.8
)

# ---- Combined summary and captions ------------------------------------------
methods <- c(
  "MedMedK",
  "MedMeanK",
  "MaxMedK",
  "MaxMeanK"
)

puech_rows <- tibble::tibble(
  Method = methods
) |>
  dplyr::left_join(
    puech |>
      dplyr::select(
        Method,
        Selected_K,
        Threshold,
        Source
      ),
    by = "Method"
  ) |>
  dplyr::transmute(
    Method,
    `Selected K` = Selected_K,
    Criterion = Method,
    Threshold,
    `Replicates per K` = EXPECTED_REPLICATES,
    Source = ifelse(
      is.na(Source),
      "Not available in parsed StructureSelector TSV",
      Source
    )
  )

combined <- dplyr::bind_rows(
  tibble::tibble(
    Method = "Evanno Delta K",
    `Selected K` = max_delta$K,
    Criterion = "Maximum finite Delta K",
    Threshold = NA_real_,
    `Replicates per K` = EXPECTED_REPLICATES,
    Source = "HEG_ZG_run2 STRUCTURE outputs"
  ),
  puech_rows
)

combined_csv <- fs::path(
  out_tables,
  "structure_K_selection_summary.csv"
)

combined_docx <- fs::path(
  out_tables,
  "structure_K_selection_summary.docx"
)

readr::write_csv(
  combined,
  combined_csv,
  na = ""
)

ft <- flextable::flextable(combined) |>
  flextable::bold(part = "header") |>
  flextable::fontsize(size = 9, part = "all") |>
  flextable::theme_booktabs() |>
  flextable::autofit()

doc <- officer::read_docx() |>
  officer::body_add_par(
    "Table M.1. Summary of the criteria used to identify the number of genetic clusters in STRUCTURE analyses.",
    style = "Normal"
  ) |>
  flextable::body_add_flextable(value = ft) |>
  officer::body_add_par(
    paste0(
      "Note. The Evanno method is based on the second-order rate of change ",
      "in STRUCTURE log probabilities across successive values of K. ",
      "Puechmaille estimators use population-level assignment patterns and ",
      "the ancestry-membership threshold reported in the table."
    ),
    style = "Normal"
  )

print(doc, target = combined_docx)

caption_M1 <- paste0(
  "Figure M.1. Mean STRUCTURE log probability of the data across tested values of K. ",
  "Points represent mean LnP(K) calculated from 20 independent runs, ",
  "and error bars represent one standard deviation."
)

caption_M2 <- paste0(
  "Figure M.2. Selection of the number of genetic clusters using the Evanno method. ",
  "ΔK was calculated from the second-order rate of change in mean STRUCTURE ",
  "log probabilities across successive values of K. The highest ΔK value indicates ",
  "the solution most strongly supported by the Evanno criterion."
)

caption_M3 <- paste0(
  "Figure M.3. Selection of the number of genetic clusters using the Puechmaille estimators. ",
  "MedMedK, MedMeanK, MaxMedK, and MaxMeanK were calculated from STRUCTURE ancestry ",
  "coefficients using the membership threshold reported in Table M.1."
)

if (VERBOSE_PARSING) {
  cat(
    "\n",
    caption_M1,
    "\n",
    caption_M2,
    "\n",
    caption_M3,
    "\n",
    sep = ""
  )
  
  message(
    "[appendix_M] Puechmaille figure type: ",
    figure_type
  )
}

message(
  "[appendix_M] Evanno selected K: ",
  max_delta$K
)

message(
  "[appendix_M] Output paths: ",
  paste(
    c(
      evanno_reps_path,
      evanno_path,
      puech_path,
      raw_path,
      combined_csv,
      combined_docx,
      fs::path(
        out_figures,
        c(
          "appendix_M_mean_LnPK.png",
          "appendix_M_mean_LnPK.pdf",
          "appendix_M_evanno_deltaK.png",
          "appendix_M_evanno_deltaK.pdf",
          "appendix_M_puechmaille_estimators.png",
          "appendix_M_puechmaille_estimators.pdf"
        )
      )
    ),
    collapse = "; "
  )
)

if (VERBOSE_PARSING) {
  print(evanno)
  print(puech)
}

message("[appendix_M] All Appendix M outputs written successfully.")