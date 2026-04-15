# Vegetation / stand structure descriptive pipeline
# Run with: source("R/vegetation_pipeline.R"); run_vegetation_pipeline()

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(readxl)
  library(tidyr)
  library(purrr)
  library(stringr)
})

source("R/utils_vegetation.R")

run_vegetation_pipeline <- function(
    workbook_path = NULL,
    output_dir = file.path("outputs", "vegetation")
) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_dir, "cleaned_data"), recursive = TRUE, showWarnings = FALSE)
  
  summary_lines <- c(
    "Vegetation ecological-context pipeline summary",
    sprintf("Generated on: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    ""
  )
  
  if (is.null(workbook_path)) {
    workbook_path <- find_workbook_path()
  }
  
  if (is.na(workbook_path) || !file.exists(workbook_path)) {
    stop("Workbook not found. Please place 'donnees_modifiees_west_summer2024 copie.xlsx' in project root or data/.")
  }
  
  message(sprintf("Using workbook: %s", workbook_path))
  
  sheets <- readxl::excel_sheets(workbook_path)
  message("Sheets found: ", paste(sheets, collapse = ", "))
  
  summary_lines <- c(summary_lines, "Sheets found:", paste0("- ", sheets), "")
  
  sheet_data <- purrr::map(sheets, ~ safe_read_sheet(workbook_path, .x))
  names(sheet_data) <- sheets
  
  target_sheet_idx <- str_detect(str_to_lower(names(sheet_data)), "arbre|gaule|regen|veget")
  selected <- sheet_data[target_sheet_idx]
  
  if (length(selected) == 0) {
    stop("No relevant sheets detected (expected names containing arbre/gaule/regeneration/vegetation).")
  }
  
  cleaned_list <- list()
  site_tables <- list()
  diversity_tables <- list()
  top_species_regen <- NULL
  basal_area_table <- NULL
  
  for (sheet_name in names(selected)) {
    dat <- selected[[sheet_name]]
    if (is.null(dat) || nrow(dat) == 0) {
      summary_lines <- c(summary_lines, sprintf("- %s: skipped (empty or unreadable)", sheet_name))
      next
    }
    
    std <- coerce_site_species(dat)
    dat2 <- std$data
    stratum <- sheet_stratum(sheet_name)
    
    cleaned_path <- file.path(output_dir, "cleaned_data", paste0("cleaned_", stratum, ".csv"))
    readr::write_csv(dat2, cleaned_path)
    
    cleaned_list[[stratum]] <- dat2
    
    missing_by_col <- colSums(is.na(dat2))
    miss_tbl <- tibble(column = names(missing_by_col), n_missing = as.integer(missing_by_col)) %>%
      arrange(desc(n_missing))
    readr::write_csv(miss_tbl, file.path(output_dir, "tables", paste0("missingness_", stratum, ".csv")))
    
    site_summary <- summarise_site_beech(dat2, stratum)
    if (!is.null(site_summary)) {
      site_tables[[stratum]] <- site_summary
      readr::write_csv(site_summary, file.path(output_dir, "tables", paste0("site_summary_", stratum, ".csv")))
    }
    
    diversity <- calc_diversity(dat2, stratum)
    if (!is.null(diversity)) {
      diversity_tables[[stratum]] <- diversity
      readr::write_csv(diversity, file.path(output_dir, "tables", paste0("diversity_", stratum, ".csv")))
    }
    
    if (stratum == "regeneration" && "species_std" %in% names(dat2) && "site_std" %in% names(dat2)) {
      count_col <- infer_count_col(dat2)
      regen_work <- if (is.null(count_col)) {
        dat2 %>% mutate(.count = 1)
      } else {
        dat2 %>% mutate(.count = coalesce(suppressWarnings(as.numeric(.data[[count_col]])), 0))
      }
      
      top_species_regen <- regen_work %>%
        filter(!is.na(species_std), species_std != "") %>%
        group_by(site_std, species_std) %>%
        summarise(abundance = sum(.count, na.rm = TRUE), .groups = "drop") %>%
        group_by(site_std) %>%
        slice_max(order_by = abundance, n = 5, with_ties = FALSE) %>%
        ungroup()
      
      readr::write_csv(top_species_regen, file.path(output_dir, "tables", "top_species_regeneration_by_site.csv"))
    }
  }
  
  if ("arbre" %in% names(cleaned_list)) {
    basal_area_table <- calc_basal_area(cleaned_list$arbre)
    if (!is.null(basal_area_table)) {
      readr::write_csv(basal_area_table, file.path(output_dir, "tables", "basal_area_by_site.csv"))
    }
  }
  
  if (length(site_tables) > 0) {
    combined_site <- bind_rows(site_tables)
    readr::write_csv(combined_site, file.path(output_dir, "tables", "site_summary_all_strata.csv"))
    
    if (all(c("site_std", "stratum", "beech_abundance") %in% names(combined_site))) {
      fig1 <- combined_site %>%
        filter(!is.na(beech_abundance)) %>%
        ggplot(aes(x = site_std, y = beech_abundance, fill = stratum)) +
        geom_col(position = position_dodge(width = 0.8), width = 0.7) +
        labs(
          x = "Site",
          y = "Beech abundance",
          fill = "Stratum",
          title = "American beech abundance by site and stratum"
        ) +
        theme_minimal(base_size = 12) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      
      ggsave(
        filename = file.path(output_dir, "figures", "beech_abundance_by_site_stratum.png"),
        plot = fig1,
        width = 10,
        height = 6,
        dpi = 300
      )
    }
  }
  
  if (!is.null(top_species_regen) && nrow(top_species_regen) > 0) {
    fig2 <- ggplot(top_species_regen, aes(x = site_std, y = abundance, fill = species_std)) +
      geom_col(position = "fill", width = 0.75) +
      scale_y_continuous(labels = scales::percent_format()) +
      labs(
        x = "Site",
        y = "Relative composition",
        fill = "Species",
        title = "Top regeneration species composition by site"
      ) +
      theme_minimal(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggsave(
      filename = file.path(output_dir, "figures", "regeneration_species_composition_top5.png"),
      plot = fig2,
      width = 11,
      height = 6,
      dpi = 300
    )
  }
  
  if (length(diversity_tables) > 0) {
    diversity_all <- bind_rows(diversity_tables)
    readr::write_csv(diversity_all, file.path(output_dir, "tables", "diversity_all_strata.csv"))
    
    fig3 <- diversity_all %>%
      ggplot(aes(x = site_std, y = species_richness, fill = stratum)) +
      geom_col(position = position_dodge(width = 0.8), width = 0.7) +
      labs(
        x = "Site",
        y = "Species richness",
        fill = "Stratum",
        title = "Species richness by site"
      ) +
      theme_minimal(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggsave(
      filename = file.path(output_dir, "figures", "species_richness_by_site.png"),
      plot = fig3,
      width = 10,
      height = 6,
      dpi = 300
    )
  }
  
  if (!is.null(basal_area_table) && nrow(basal_area_table) > 0) {
    basal_long <- basal_area_table %>%
      select(site_std, total_basal_area_m2, beech_basal_area_m2) %>%
      pivot_longer(
        cols = c(total_basal_area_m2, beech_basal_area_m2),
        names_to = "metric",
        values_to = "value"
      ) %>%
      mutate(
        metric = recode(
          metric,
          total_basal_area_m2 = "Total basal area",
          beech_basal_area_m2 = "Beech basal area"
        )
      )
    
    fig4 <- ggplot(basal_long, aes(x = site_std, y = value, fill = metric)) +
      geom_col(position = position_dodge(width = 0.8), width = 0.7) +
      labs(
        x = "Site",
        y = expression("Basal area (m"^2*")"),
        fill = "Metric",
        title = "Total versus beech basal area by site"
      ) +
      theme_minimal(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggsave(
      filename = file.path(output_dir, "figures", "basal_area_total_vs_beech_by_site.png"),
      plot = fig4,
      width = 10,
      height = 6,
      dpi = 300
    )
  }
  
  compact_final <- NULL
  if (length(site_tables) > 0) {
    combined_site <- bind_rows(site_tables)
    compact_final <- combined_site %>%
      filter(stratum %in% c("arbre", "gaule", "regeneration")) %>%
      select(site_std, stratum, beech_abundance) %>%
      mutate(beech_abundance = coalesce(beech_abundance, 0)) %>%
      pivot_wider(
        names_from = stratum,
        values_from = beech_abundance,
        names_prefix = "beech_abundance_"
      )
    
    if (length(diversity_tables) > 0) {
      rich_total <- bind_rows(diversity_tables) %>%
        group_by(site_std) %>%
        summarise(total_richness_across_strata = sum(species_richness, na.rm = TRUE), .groups = "drop")
      compact_final <- left_join(compact_final, rich_total, by = "site_std")
    }
    
    if (!is.null(basal_area_table)) {
      compact_final <- left_join(
        compact_final,
        basal_area_table %>% select(site_std, total_basal_area_m2),
        by = "site_std"
      )
    }
    
    readr::write_csv(compact_final, file.path(output_dir, "tables", "site_compact_context_table.csv"))
  }
  
  regen_site_summary <- NULL
  if (length(site_tables) > 0) {
    st <- bind_rows(site_tables)
    regen_site_summary <- st %>%
      filter(stratum %in% c("regeneration", "gaule")) %>%
      select(site_std, stratum, total_abundance, beech_abundance, beech_prop) %>%
      arrange(site_std, stratum)
    
    if (nrow(regen_site_summary) > 0) {
      readr::write_csv(regen_site_summary, file.path(output_dir, "tables", "regeneration_gaule_summary_by_site.csv"))
    }
  }
  
  summary_lines <- c(summary_lines, "Analyses run:")
  summary_lines <- c(summary_lines, sprintf("- Cleaned sheets exported: %d", length(cleaned_list)))
  summary_lines <- c(summary_lines, sprintf("- Site-level tables produced: %d", length(site_tables)))
  summary_lines <- c(summary_lines, sprintf("- Diversity tables produced: %d", length(diversity_tables)))
  summary_lines <- c(summary_lines, sprintf("- Basal area computed: %s", ifelse(is.null(basal_area_table), "no", "yes")))
  summary_lines <- c(summary_lines, sprintf("- Regeneration top-species table: %s", ifelse(is.null(top_species_regen), "no", "yes")))
  
  summary_lines <- c(summary_lines, "", "Variable availability by cleaned sheet:")
  for (nm in names(cleaned_list)) {
    vars <- names(cleaned_list[[nm]])
    summary_lines <- c(summary_lines, sprintf("- %s: %s", nm, paste(vars, collapse = ", ")))
  }
  
  if (!is.null(regen_site_summary) && nrow(regen_site_summary) > 0) {
    high_beech_lower <- regen_site_summary %>%
      group_by(stratum) %>%
      summarise(
        mean_beech_prop = mean(beech_prop, na.rm = TRUE),
        median_beech_prop = median(beech_prop, na.rm = TRUE),
        .groups = "drop"
      )
    
    summary_lines <- c(summary_lines, "", "Lower-strata descriptive findings (non-causal):")
    for (i in seq_len(nrow(high_beech_lower))) {
      summary_lines <- c(
        summary_lines,
        sprintf(
          "- %s: mean beech proportion = %.2f; median = %.2f",
          high_beech_lower$stratum[[i]],
          high_beech_lower$mean_beech_prop[[i]],
          high_beech_lower$median_beech_prop[[i]]
        )
      )
    }
    summary_lines <- c(
      summary_lines,
      "- Interpretation note: high small-stem abundance is descriptive and does not alone demonstrate clonality."
    )
  } else {
    summary_lines <- c(
      summary_lines,
      "",
      "Lower-strata descriptive findings: unavailable (required columns not detected for gaule/regeneration)."
    )
  }
  
  writeLines(summary_lines, con = file.path(output_dir, "vegetation_summary.txt"))
  
  message("Vegetation pipeline completed. Outputs written to: ", output_dir)
  
  invisible(
    list(
      workbook = workbook_path,
      sheets = sheets,
      cleaned = cleaned_list,
      site_tables = site_tables,
      diversity_tables = diversity_tables,
      basal_area = basal_area_table,
      compact_table = compact_final
    )
  )
}