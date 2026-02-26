# filename: scripts/09_heterozygosity_fis.R
############################################################
# scripts/09_heterozygosity_fis.R
# Ho / He(Hs) / FIS summaries from clone-corrected genind
############################################################

options(repos = c(CRAN = "https://cloud.r-project.org"))
pkgs <- c("adegenet", "hierfstat", "dplyr", "tidyr", "readr", "ggplot2")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)

suppressPackageStartupMessages({
  library(adegenet)
  library(hierfstat)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
})

if (!exists("gi_mll", inherits = TRUE) || !inherits(get("gi_mll", inherits = TRUE), "genind")) {
  stop("Clone-corrected genind object gi_mll not found; run clone correction first.")
}

gi_mll <- get("gi_mll", inherits = TRUE)

if (nInd(gi_mll) < 1 || nLoc(gi_mll) < 1) {
  stop("gi_mll exists but has no individuals or loci; cannot compute heterozygosity metrics.")
}

project_root <- if (exists("PROJECT_ROOT", inherits = TRUE)) {
  get("PROJECT_ROOT", inherits = TRUE)
} else {
  getwd()
}

run_tag <- if (exists("RUN_TAG", inherits = TRUE)) get("RUN_TAG", inherits = TRUE) else "v1"
out_dir <- if (exists("RUN_OUT", inherits = TRUE)) {
  get("RUN_OUT", inherits = TRUE)
} else {
  file.path(project_root, "outputs", run_tag)
}

diversity_dir <- file.path(out_dir, "diversity")
dir.create(diversity_dir, showWarnings = FALSE, recursive = TRUE)

pop_vec <- as.character(pop(gi_mll))
loci_vec <- as.character(locNames(gi_mll))
if (length(loci_vec) == 0) stop("No locus names found in gi_mll.")

pop_counts <- table(pop_vec)
pop_counts_df <- tibble(
  pop = names(pop_counts),
  n_ind_pop = as.integer(pop_counts)
)

hf <- genind2hierfstat(gi_mll)
bs <- basic.stats(hf)

as_pop_locus <- function(metric_obj, metric_name, pops, loci) {
  base_grid <- tidyr::expand_grid(pop = unique(pops), locus = loci)
  
  if (is.null(metric_obj)) {
    base_grid[[metric_name]] <- NA_real_
    return(base_grid)
  }
  
  if (is.matrix(metric_obj) || is.data.frame(metric_obj)) {
    mm <- as.matrix(metric_obj)
    rr <- rownames(mm)
    cc <- colnames(mm)
    
    if (is.null(rr) || length(rr) == 0) rr <- unique(pops)[seq_len(nrow(mm))]
    if (is.null(cc) || length(cc) == 0) cc <- loci[seq_len(ncol(mm))]
    
    out <- as.data.frame(mm, stringsAsFactors = FALSE)
    out$pop <- rr
    out <- out %>%
      pivot_longer(cols = -pop, names_to = "locus", values_to = metric_name) %>%
      mutate(
        pop = as.character(pop),
        locus = as.character(locus),
        !!metric_name := as.numeric(.data[[metric_name]])
      )
    
    return(base_grid %>% left_join(out, by = c("pop", "locus")))
  }
  
  vv <- as.numeric(metric_obj)
  
  if (length(vv) == length(loci)) {
    if (length(unique(pops)) == 1) {
      out <- tibble(pop = unique(pops), locus = loci, val = vv)
      names(out)[names(out) == "val"] <- metric_name
      return(base_grid %>% left_join(out, by = c("pop", "locus")))
    }
    
    warning(metric_name, " appears to be locus-level only; pop-specific values set to NA.")
    base_grid[[metric_name]] <- NA_real_
    return(base_grid)
  }
  
  base_grid[[metric_name]] <- NA_real_
  warning("Could not safely reshape ", metric_name, "; values set to NA.")
  base_grid
}

ho_long <- as_pop_locus(bs$Ho, "Ho", pop_vec, loci_vec)
he_long <- as_pop_locus(bs$Hs, "He", pop_vec, loci_vec)
fis_long <- as_pop_locus(bs$Fis, "FIS", pop_vec, loci_vec)

hetero_by_pop_locus <- ho_long %>%
  left_join(he_long, by = c("pop", "locus")) %>%
  left_join(fis_long, by = c("pop", "locus")) %>%
  left_join(pop_counts_df, by = "pop") %>%
  select(pop, locus, n_ind_pop, Ho, He, FIS) %>%
  arrange(pop, locus)

hetero_by_pop <- hetero_by_pop_locus %>%
  group_by(pop, n_ind_pop) %>%
  summarise(
    mean_Ho = mean(Ho, na.rm = TRUE),
    sd_Ho = sd(Ho, na.rm = TRUE),
    mean_He = mean(He, na.rm = TRUE),
    sd_He = sd(He, na.rm = TRUE),
    mean_FIS = mean(FIS, na.rm = TRUE),
    sd_FIS = sd(FIS, na.rm = TRUE),
    n_loci_non_na = sum(!is.na(Ho) | !is.na(He) | !is.na(FIS)),
    .groups = "drop"
  ) %>%
  mutate(
    across(c(mean_Ho, sd_Ho, mean_He, sd_He, mean_FIS, sd_FIS),
           ~ ifelse(is.nan(.x), NA_real_, .x))
  ) %>%
  arrange(desc(n_ind_pop), pop)

hetero_by_pop_nge10 <- hetero_by_pop %>%
  filter(n_ind_pop >= 10)

hetero_by_locus <- hetero_by_pop_locus %>%
  group_by(locus) %>%
  summarise(
    mean_Ho = mean(Ho, na.rm = TRUE),
    sd_Ho = sd(Ho, na.rm = TRUE),
    mean_He = mean(He, na.rm = TRUE),
    sd_He = sd(He, na.rm = TRUE),
    mean_FIS = mean(FIS, na.rm = TRUE),
    sd_FIS = sd(FIS, na.rm = TRUE),
    n_pops_non_na = sum(!is.na(Ho) | !is.na(He) | !is.na(FIS)),
    .groups = "drop"
  ) %>%
  mutate(
    across(c(mean_Ho, sd_Ho, mean_He, sd_He, mean_FIS, sd_FIS),
           ~ ifelse(is.nan(.x), NA_real_, .x))
  ) %>%
  arrange(locus)

get_overall_value <- function(overall_obj, key) {
  if (is.null(overall_obj)) return(NA_real_)
  nm <- names(overall_obj)
  if (is.null(nm)) return(NA_real_)
  idx <- which(tolower(nm) == tolower(key))
  if (length(idx) == 0) return(NA_real_)
  as.numeric(overall_obj[idx[1]])
}

overall_ho <- get_overall_value(bs$overall, "Ho")
overall_he <- get_overall_value(bs$overall, "Hs")
overall_fis <- get_overall_value(bs$overall, "Fis")

# Fallbacks if basic.stats() overall keys are unavailable in this hierfstat version.
if (is.na(overall_ho)) overall_ho <- mean(hetero_by_pop_locus$Ho, na.rm = TRUE)
if (is.na(overall_he)) overall_he <- mean(hetero_by_pop_locus$He, na.rm = TRUE)
if (is.na(overall_fis)) {
  if (!is.na(overall_ho) && !is.na(overall_he) && overall_he != 0) {
    overall_fis <- 1 - (overall_ho / overall_he)
  } else {
    overall_fis <- mean(hetero_by_pop_locus$FIS, na.rm = TRUE)
  }
}

overall_df <- tibble(
  n_ind_total = nInd(gi_mll),
  n_loci_total = nLoc(gi_mll),
  n_pops_total = nPop(gi_mll),
  Ho = ifelse(is.nan(overall_ho), NA_real_, overall_ho),
  He = ifelse(is.nan(overall_he), NA_real_, overall_he),
  FIS = ifelse(is.nan(overall_fis), NA_real_, overall_fis)
)

readr::write_csv(hetero_by_pop_locus, file.path(diversity_dir, "heterozygosity_by_pop_locus.csv"), na = "")
readr::write_csv(hetero_by_pop, file.path(diversity_dir, "heterozygosity_by_pop.csv"), na = "")
readr::write_csv(hetero_by_pop_nge10, file.path(diversity_dir, "heterozygosity_by_pop_nge10.csv"), na = "")
readr::write_csv(hetero_by_locus, file.path(diversity_dir, "heterozygosity_by_locus.csv"), na = "")
readr::write_csv(overall_df, file.path(diversity_dir, "heterozygosity_overall.csv"), na = "")

p_ho <- ggplot(hetero_by_pop_locus, aes(x = pop, y = Ho)) +
  geom_boxplot(na.rm = TRUE) +
  labs(title = "Observed heterozygosity (Ho) by population", x = "Population", y = "Ho") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  filename = file.path(diversity_dir, "heterozygosity_Ho_by_population.png"),
  plot = p_ho, width = 10, height = 6, dpi = 150
)

p_he <- ggplot(hetero_by_pop_locus, aes(x = pop, y = He)) +
  geom_boxplot(na.rm = TRUE) +
  labs(title = "Expected heterozygosity (He) by population", x = "Population", y = "He") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  filename = file.path(diversity_dir, "heterozygosity_He_by_population.png"),
  plot = p_he, width = 10, height = 6, dpi = 150
)

p_fis <- ggplot(hetero_by_pop, aes(x = pop, y = mean_FIS)) +
  geom_col(na.rm = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "FIS by population", x = "Population", y = "Mean FIS across loci") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  filename = file.path(diversity_dir, "heterozygosity_FIS_by_population.png"),
  plot = p_fis, width = 10, height = 6, dpi = 150
)

cat("[INFO] Ho/He/FIS outputs saved in: ", diversity_dir, "\n", sep = "")
cat("[INFO] Files: heterozygosity_by_pop_locus.csv, heterozygosity_by_pop.csv, ",
    "heterozygosity_by_pop_nge10.csv, heterozygosity_by_locus.csv, heterozygosity_overall.csv\n", sep = "")
