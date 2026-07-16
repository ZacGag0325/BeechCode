#!/usr/bin/env Rscript
# Appendix M: STRUCTURE K-selection diagnostics (Evanno and Puechmaille).
# This script only writes derived diagnostic outputs under outputs/.

# ---- Configuration and package checks ---------------------------------------
project_root <- path.expand("~/Desktop/BeechCode")
if (!dir.exists(project_root)) {
  stop("[appendix_M] Project root does not exist: ", project_root, call. = FALSE)
}

required_packages <- c(
  "dplyr", "tidyr", "readr", "stringr", "purrr", "tibble",
  "ggplot2", "fs", "officer", "flextable"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages)) {
  stop(
    "[appendix_M] Missing required R package(s): ",
    paste(missing_packages, collapse = ", "),
    ". Install them before running this script.",
    call. = FALSE
  )
}

invisible(lapply(required_packages, library, character.only = TRUE))

EXPECTED_K <- 1:12
EXPECTED_REPLICATES <- 20L

SITE_LEVELS <- c(
  "S1", "S2", "S3", "S4", "S5", "S6",
  "N1", "N2", "N3", "N4", "N5", "N6"
)

site_key <- c(
  AMC = "S1", ALB = "S2", IKJ = "S3", IKO = "S4",
  LGG = "S5", LGR = "S6", ML1 = "N1", ML2 = "N2",
  ML3 = "N3", CPF = "N4", PLI = "N5", LDF = "N6"
)

out_table_dir <- fs::path(
  project_root, "outputs", "tables", "supplementary"
)

out_figure_dir <- fs::path(
  project_root, "outputs", "figures", "supplementary"
)

out_derived_dir <- fs::path(project_root, "outputs", "derived")

fs::dir_create(
  c(out_table_dir, out_figure_dir, out_derived_dir),
  recurse = TRUE
)

abort_m <- function(...) {
  stop("[appendix_M] ", paste0(..., collapse = ""), call. = FALSE)
}

num_pattern <- "[-+]?(?:[0-9]*\\.?[0-9]+|[0-9]+\\.?)(?:[eE][-+]?[0-9]+)?"

extract_number <- function(x, label_regex) {
  hit <- stringr::str_match(
    x,
    paste0(
      "(?i)",
      label_regex,
      ".{0,80}?(?:[:=]|\\s)\\s*(",
      num_pattern,
      ")"
    )
  )[, 2]
  
  hit <- hit[!is.na(hit)]
  
  if (!length(hit)) {
    NA_real_
  } else {
    suppressWarnings(as.numeric(hit[[1]]))
  }
}

extract_integer <- function(x, regex) {
  hit <- stringr::str_match(
    x,
    paste0(
      "(?i)",
      regex,
      ".{0,60}?(?:[:=]|\\s)\\s*([0-9]+)"
    )
  )[, 2]
  
  hit <- hit[!is.na(hit)]
  
  if (!length(hit)) {
    NA_integer_
  } else {
    as.integer(hit[[1]])
  }
}

# ---- Input discovery ---------------------------------------------------------
all_files <- fs::dir_ls(
  project_root,
  recurse = TRUE,
  type = "file",
  fail = FALSE
)

all_dirs <- fs::dir_ls(
  project_root,
  recurse = TRUE,
  type = "directory",
  fail = FALSE
)

rel_files <- fs::path_rel(all_files, start = project_root)

candidate_re <- paste0(
  "(?i)(",
  "HEG[_-]?ZG[_-]?run2|",
  "structure|",
  "clumpak|",
  "structure.?harvester|",
  "structureselector|",
  "lnp\\(?d\\)?|",
  "medmedk|",
  "medmeank|",
  "maxmedk|",
  "maxmeank|",
  "_f$|",
  "\\.q$|",
  "qmatrix",
  ")"
)

candidates <- unique(c(
  all_files[stringr::str_detect(rel_files, candidate_re)],
  all_dirs[
    stringr::str_detect(
      fs::path_rel(all_dirs, start = project_root),
      candidate_re
    )
  ]
))

message("[appendix_M] Candidate paths found (all are read-only inputs):")

if (length(candidates)) {
  for (p in candidates) {
    message("  ", p)
  }
} else {
  message("  (none)")
}

run2_paths <- candidates[
  stringr::str_detect(
    fs::path_file(candidates),
    "(?i)HEG[_-]?ZG[_-]?run2"
  ) |
    stringr::str_detect(
      candidates,
      "(?i)HEG[_-]?ZG[_-]?run2"
    )
]

zip_paths <- run2_paths[
  stringr::str_detect(run2_paths, "(?i)\\.zip$")
]

search_roots <- unique(c(
  project_root,
  run2_paths[fs::dir_exists(run2_paths)]
))

if (length(zip_paths)) {
  extract_dir <- fs::path(
    out_derived_dir,
    "HEG_ZG_run2_logs_extracted"
  )
  
  fs::dir_create(extract_dir, recurse = TRUE)
  
  for (z in zip_paths) {
    message(
      "[appendix_M] Extracting read-only archive to: ",
      extract_dir
    )
    
    utils::unzip(
      z,
      exdir = extract_dir,
      overwrite = FALSE
    )
  }
  
  search_roots <- unique(c(search_roots, extract_dir))
}

# ---- STRUCTURE output parsing ------------------------------------------------
parse_k_replicate <- function(path) {
  text <- paste(
    fs::path_file(path),
    fs::path_rel(path, start = project_root),
    sep = " "
  )
  
  k <- stringr::str_match(
    text,
    "(?i)(?:^|[_ .-])K[_ .-]?([0-9]{1,2})(?:[_ .-]|$)"
  )[, 2]
  
  if (is.na(k)) {
    k <- stringr::str_match(
      text,
      "(?i)K=([0-9]{1,2})"
    )[, 2]
  }
  
  r <- stringr::str_match(
    text,
    "(?i)(?:rep(?:licate)?|run)[_ .-]?([0-9]+)"
  )[, 2]
  
  if (is.na(r)) {
    r <- stringr::str_match(
      text,
      "(?i)K[_ .-]?[0-9]+[_ .-]+([0-9]+)(?:[_ .-]|$)"
    )[, 2]
  }
  
  list(
    K = suppressWarnings(as.integer(k)),
    replicate = ifelse(is.na(r), NA_character_, r)
  )
}

parse_structure_file <- function(path) {
  x <- tryCatch(
    readLines(path, warn = FALSE, encoding = "UTF-8"),
    error = function(e) character()
  )
  
  if (!length(x)) {
    return(NULL)
  }
  
  kr <- parse_k_replicate(path)
  
  lnp <- extract_number(
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
    K = kr$K,
    replicate = kr$replicate,
    LnP = lnp,
    mean_ln_likelihood = extract_number(
      x,
      "mean\\s+value\\s+of\\s+ln\\s+likelihood"
    ),
    variance_ln_likelihood = extract_number(
      x,
      "variance\\s+of\\s+ln\\s+likelihood"
    ),
    n_individuals = extract_integer(
      x,
      "(?:number\\s+of\\s+individuals|numind|individuals)"
    ),
    n_loci = extract_integer(
      x,
      "(?:number\\s+of\\s+loci|numloci|loci)"
    ),
    burnin = extract_integer(
      x,
      "(?:burnin|burn-in)"
    ),
    mcmc_iterations = extract_integer(
      x,
      "(?:numreps|mcmc(?:\\s+iterations)?)"
    ),
    path = path
  )
}

possible_outputs <- unique(unlist(
  lapply(search_roots, function(root) {
    fs::dir_ls(
      root,
      recurse = TRUE,
      type = "file",
      fail = FALSE,
      regexp = "(?i)(?:_f$|\\.out$|\\.txt$|\\.log$|result|structure)"
    )
  }),
  use.names = FALSE
))

possible_outputs <- unique(c(
  run2_paths[
    !fs::dir_exists(run2_paths) &
      !stringr::str_detect(run2_paths, "(?i)\\.zip$")
  ],
  possible_outputs
))

if (length(zip_paths)) {
  possible_outputs <- possible_outputs[
    stringr::str_detect(
      possible_outputs,
      "(?i)HEG[_-]?ZG[_-]?run2"
    ) |
      fs::path_has_parent(possible_outputs, extract_dir)
  ]
} else {
  possible_outputs <- possible_outputs[
    stringr::str_detect(
      possible_outputs,
      "(?i)HEG[_-]?ZG[_-]?run2"
    )
  ]
}

parsed <- purrr::map(
  possible_outputs,
  parse_structure_file
) |>
  purrr::compact()

if (!length(parsed)) {
  abort_m("No final HEG_ZG_run2 STRUCTURE outputs were found.")
}

replicates <- dplyr::bind_rows(parsed) |>
  dplyr::filter(
    !is.na(K),
    K %in% EXPECTED_K,
    !is.na(replicate)
  )

if (!nrow(replicates)) {
  abort_m(
    "No parseable K/replicate/LnP(D) records were found ",
    "in final HEG_ZG_run2 outputs."
  )
}

validate_evanno <- function(x) {
  counts <- x |>
    dplyr::count(K, name = "n") |>
    tidyr::complete(K = EXPECTED_K, fill = list(n = 0L))
  
  bad <- dplyr::filter(counts, n != EXPECTED_REPLICATES)
  
  if (nrow(bad)) {
    for (k in bad$K) {
      message("[appendix_M] K = ", k, " candidate files:")
      message(
        paste0("  ", x$path[x$K == k]),
        collapse = "\n"
      )
    }
    
    abort_m(
      "K = ", bad$K[[1]],
      " contains only ", bad$n[[1]],
      " valid replicates; expected ",
      EXPECTED_REPLICATES,
      "."
    )
  }
  
  duplicate <- x |>
    dplyr::count(K, replicate) |>
    dplyr::filter(n > 1)
  
  if (nrow(duplicate)) {
    abort_m(
      "Duplicate K-replicate combination(s): ",
      paste(
        paste0(
          "K=", duplicate$K,
          ", replicate=", duplicate$replicate
        ),
        collapse = "; "
      )
    )
  }
  
  if (any(!is.finite(x$LnP))) {
    abort_m("Non-finite LnP(D) values were detected.")
  }
  
  known_n <- unique(stats::na.omit(x$n_individuals))
  
  if (!length(known_n)) {
    abort_m(
      "The number of individuals could not be determined ",
      "from final STRUCTURE files."
    )
  }
  
  if (length(known_n) != 1L) {
    abort_m(
      "Individual counts are inconsistent across final STRUCTURE outputs: ",
      paste(known_n, collapse = ", ")
    )
  }
  
  invisible(counts)
}

replicate_counts <- validate_evanno(replicates)

n_individuals <- unique(
  stats::na.omit(replicates$n_individuals)
)[[1]]

# ---- Evanno calculations, tables, and figures --------------------------------
evanno_summary <- replicates |>
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
      is.finite(L_double_prime_K) &
        is.finite(sd_LnP_K) &
        sd_LnP_K > 0,
      L_double_prime_K / sd_LnP_K,
      NA_real_
    )
  )

max_delta <- evanno_summary |>
  dplyr::filter(is.finite(DeltaK)) |>
  dplyr::slice_max(DeltaK, n = 1, with_ties = FALSE)

if (!nrow(max_delta)) {
  abort_m("No finite Delta K value could be calculated.")
}

evanno_replicates_path <- fs::path(
  out_table_dir,
  "evanno_replicates.csv"
)

evanno_summary_path <- fs::path(
  out_table_dir,
  "evanno_summary.csv"
)

readr::write_csv(replicates, evanno_replicates_path)
readr::write_csv(evanno_summary, evanno_summary_path)

theme_m <- ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major.x = ggplot2::element_blank()
  )

p_lnp <- ggplot2::ggplot(
  evanno_summary,
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
  evanno_summary,
  ggplot2::aes(K, DeltaK)
) +
  ggplot2::geom_line(linewidth = 0.45, na.rm = TRUE) +
  ggplot2::geom_point(size = 2, na.rm = TRUE) +
  ggplot2::geom_point(
    data = max_delta,
    colour = "#9A4D3E",
    fill = "white",
    shape = 21,
    stroke = 0.8,
    size = 3.5
  ) +
  ggplot2::scale_x_continuous(breaks = EXPECTED_K) +
  ggplot2::labs(
    x = "Number of genetic clusters (K)",
    y = "ΔK"
  ) +
  theme_m

save_plot <- function(plot, stem) {
  ggplot2::ggsave(
    fs::path(out_figure_dir, paste0(stem, ".png")),
    plot,
    width = 6.5,
    height = 4.2,
    units = "in",
    dpi = 600
  )
  
  ggplot2::ggsave(
    fs::path(out_figure_dir, paste0(stem, ".pdf")),
    plot,
    width = 6.5,
    height = 4.2,
    units = "in"
  )
}

save_plot(p_lnp, "appendix_M_mean_LnPK")
save_plot(p_delta, "appendix_M_evanno_deltaK")

# ---- Puechmaille: parse validated existing results only ---------------------
# Puechmaille (2016) must use consistently aligned Q matrices and confirmed
# population membership. This script deliberately refuses to substitute a
# heuristic if no validated project/StructureSelector result is available.

parse_existing_puechmaille <- function(path) {
  dat <- tryCatch(
    readr::read_delim(
      path,
      delim = readr::guess_delim(
        readLines(path, n = 3, warn = FALSE)
      ),
      show_col_types = FALSE,
      name_repair = "minimal"
    ),
    error = function(e) NULL
  )
  
  if (is.null(dat)) {
    return(NULL)
  }
  
  nm <- names(dat)
  
  key <- c(
    "MedMedK",
    "MedMeanK",
    "MaxMedK",
    "MaxMeanK"
  )
  
  if (!all(key %in% nm)) {
    return(NULL)
  }
  
  k_col <- nm[
    stringr::str_detect(
      nm,
      "(?i)^k$|k tested|clusters"
    )
  ][1]
  
  if (is.na(k_col)) {
    return(NULL)
  }
  
  out <- dplyr::transmute(
    dat,
    K = suppressWarnings(as.integer(.data[[k_col]])),
    MedMedK = suppressWarnings(as.numeric(MedMedK)),
    MedMeanK = suppressWarnings(as.numeric(MedMeanK)),
    MaxMedK = suppressWarnings(as.numeric(MaxMedK)),
    MaxMeanK = suppressWarnings(as.numeric(MaxMeanK))
  ) |>
    dplyr::filter(K %in% EXPECTED_K) |>
    dplyr::arrange(K)
  
  if (
    !all(EXPECTED_K %in% out$K) ||
    any(!is.finite(as.matrix(out[, -1])))
  ) {
    return(NULL)
  }
  
  out
}

puech_candidates <- candidates[
  stringr::str_detect(
    candidates,
    "(?i)(puechmaille|structureselector|medmedk|medmeank|maxmedk|maxmeank)"
  )
]

puech_found <- purrr::map(
  puech_candidates,
  parse_existing_puechmaille
)

valid_idx <- which(
  vapply(puech_found, Negate(is.null), logical(1))
)

if (!length(valid_idx)) {
  message(
    "[appendix_M] Evanno outputs were written successfully. ",
    "No complete, validated existing Puechmaille result was found."
  )
  
  abort_m(
    "Puechmaille estimators could not be calculated from the available files. ",
    "A complete validated StructureSelector/project Puechmaille table, ",
    "or reproducibly aligned Q matrices plus exact individual-order site metadata, ",
    "is required; no unverified approximation was made."
  )
}

puech_source <- puech_candidates[[valid_idx[[1]]]]
puech <- puech_found[[valid_idx[[1]]]]

threshold <- NA_real_

source_text <- paste(
  readLines(puech_source, warn = FALSE),
  collapse = " "
)

threshold <- extract_number(
  source_text,
  "(?:threshold|membership\\s+cutoff)"
)

if (!is.finite(threshold)) {
  threshold <- 0.5
  
  warning(
    "[appendix_M] Puechmaille threshold was not recorded; ",
    "using explicit fallback 0.5."
  )
}

puech_summary <- dplyr::mutate(
  puech,
  threshold = threshold,
  n_replicates = EXPECTED_REPLICATES,
  n_individuals = n_individuals,
  source_or_method = paste(
    "Read validated existing result:",
    puech_source
  )
)

puech_path <- fs::path(
  out_table_dir,
  "puechmaille_summary.csv"
)

readr::write_csv(puech_summary, puech_path)

# These estimators identify selected K values, rather than continuous
# quantities; plot the selected K for each estimator as a compact,
# scientifically interpretable bar chart.

selected <- tibble::tibble(
  Method = c(
    "MedMedK",
    "MedMeanK",
    "MaxMedK",
    "MaxMeanK"
  ),
  Selected_K = c(
    puech$K[which.max(puech$MedMedK)],
    puech$K[which.max(puech$MedMeanK)],
    puech$K[which.max(puech$MaxMedK)],
    puech$K[which.max(puech$MaxMeanK)]
  )
)

message(
  "[appendix_M] Puechmaille estimators are shown as selected-K bars ",
  "because their primary interpretation is the K selected by each criterion."
)

p_puech <- ggplot2::ggplot(
  selected,
  ggplot2::aes(Method, Selected_K, fill = Method)
) +
  ggplot2::geom_col(width = 0.62, show.legend = TRUE) +
  ggplot2::scale_y_continuous(
    breaks = EXPECTED_K,
    limits = c(0, max(EXPECTED_K))
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      "#4C78A8",
      "#72A5A1",
      "#9A7D4F",
      "#805D93"
    )
  ) +
  ggplot2::labs(
    x = NULL,
    y = "Selected number of genetic clusters (K)",
    fill = "Estimator"
  ) +
  theme_m +
  ggplot2::theme(legend.position = "bottom")

save_plot(p_puech, "appendix_M_puechmaille_estimators")

# ---- Combined summary and Word table ----------------------------------------
combined <- dplyr::bind_rows(
  tibble::tibble(
    Method = "Evanno ΔK",
    `Selected K` = max_delta$K,
    Criterion = "Maximum finite ΔK",
    Threshold = NA_real_,
    `Replicates per K` = EXPECTED_REPLICATES,
    Source = "HEG_ZG_run2 STRUCTURE outputs"
  ),
  dplyr::transmute(
    selected,
    Method,
    `Selected K` = Selected_K,
    Criterion = Method,
    Threshold = threshold,
    `Replicates per K` = EXPECTED_REPLICATES,
    Source = paste(
      "Read validated existing result:",
      puech_source
    )
  )
)

combined_csv <- fs::path(
  out_table_dir,
  "structure_K_selection_summary.csv"
)

combined_docx <- fs::path(
  out_table_dir,
  "structure_K_selection_summary.docx"
)

readr::write_csv(combined, combined_csv, na = "")

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

# ---- Captions and final report -----------------------------------------------
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

message("\n[appendix_M] Final report")
message("Project root: ", project_root)

message(
  "Final STRUCTURE source directory/archive: ",
  paste(unique(c(search_roots, zip_paths)), collapse = "; ")
)

message("Valid STRUCTURE output files: ", nrow(replicates))

message(
  "Replicates per K: ",
  paste0(
    replicate_counts$K,
    "=",
    replicate_counts$n,
    collapse = ", "
  )
)

message("Individuals detected: ", n_individuals)
message("Maximum ΔK: K = ", max_delta$K)

message(
  "Puechmaille selected K: ",
  paste0(
    selected$Method,
    "=",
    selected$Selected_K,
    collapse = ", "
  )
)

message(
  "Puechmaille threshold: ",
  threshold,
  "; values were read from: ",
  puech_source
)

message(
  "Outputs: ",
  paste(
    c(
      evanno_replicates_path,
      evanno_summary_path,
      puech_path,
      combined_csv,
      combined_docx,
      fs::path(
        out_figure_dir,
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

print(evanno_summary)
print(puech_summary)