#!/usr/bin/env Rscript

# Alias wrapper for US spelling.
# For full script and settings, edit/run:
#   scripts/nearest_neighbour_adjacent_numbers_analysis.R

args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
script_dir <- if (length(file_arg) > 0) {
  dirname(normalizePath(file_arg[1], winslash = "/", mustWork = FALSE))
} else {
  getwd()
}

target <- file.path(script_dir, "nearest_neighbour_adjacent_numbers_analysis.R")
if (!file.exists(target)) {
  target <- "scripts/nearest_neighbour_adjacent_numbers_analysis.R"
}

source(target)