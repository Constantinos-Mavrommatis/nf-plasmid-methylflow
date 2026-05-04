#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(arrow)
  library(dplyr)
  library(fs)
})

option_list <- list(
  make_option(
    "--outdir",
    type = "character",
    help = "Base output directory for combined results"
  ),
  make_option(
    "--pileup-list",
    type = "character",
    default = NULL,
    help = "Text file containing one pileup parquet path per line"
  ),
  make_option(
    "--extract-list",
    type = "character",
    default = NULL,
    help = "Text file containing one collapsed extract parquet path per line"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$outdir)) {
  stop("You must provide --outdir", call. = FALSE)
}

out_dir <- opt$outdir

if (!dir_exists(out_dir)) {
  dir_create(out_dir, recurse = TRUE)
}

read_manifest <- function(path) {
  if (is.null(path) || !file.exists(path)) {
    return(character())
  }

  entries <- trimws(readLines(path, warn = FALSE))
  entries[nzchar(entries)]
}

pileup_files  <- read_manifest(opt$`pileup-list`)
extract_files <- read_manifest(opt$`extract-list`)

if (length(pileup_files) == 0L && length(extract_files) == 0L) {
  stop("No input parquet files were provided to combine_run.R", call. = FALSE)
}

if (length(pileup_files) > 0L) {
  message("Combining ", length(pileup_files), " pileup parquet file(s)")
  pileup_df <- arrow::open_dataset(pileup_files, format = "parquet") %>%
    dplyr::collect()

  needed_cols <- c("treatment", "prob_threshold")
  if (!all(needed_cols %in% colnames(pileup_df))) {
    stop(
      "Pileup parquet files must contain columns: ",
      paste(needed_cols, collapse = ", "),
      call. = FALSE
    )
  }

  if (any(is.na(pileup_df$treatment)) || any(is.na(pileup_df$prob_threshold))) {
    stop("Pileup parquet files contain missing treatment or prob_threshold values", call. = FALSE)
  }

  pileup_targets <- pileup_df %>%
    dplyr::distinct(treatment, prob_threshold) %>%
    dplyr::arrange(treatment, prob_threshold)

  for (i in seq_len(nrow(pileup_targets))) {
    treatment_name <- pileup_targets$treatment[[i]]
    threshold_value <- pileup_targets$prob_threshold[[i]]
    threshold_label <- paste0(as.character(threshold_value), "_modification_threshold")
    combined_dir <- fs::path(out_dir, treatment_name, "05_combined")
    dir_create(combined_dir, recurse = TRUE)

    out_pileup <- fs::path(
      combined_dir,
      paste0(treatment_name, "_", threshold_label, "_pileup_combined.parquet")
    )

    combined_pileup <- pileup_df %>%
      dplyr::filter(treatment == !!treatment_name, prob_threshold == !!threshold_value)

    arrow::write_parquet(
      combined_pileup,
      out_pileup,
      compression = "zstd"
    )

    cat("Wrote combined pileup:", out_pileup, "\n")
  }
}

if (length(extract_files) > 0L) {
  message("Combining ", length(extract_files), " collapsed extract parquet file(s)")
  extract_df <- arrow::open_dataset(extract_files, format = "parquet") %>%
    dplyr::collect()

  if (!"treatment" %in% colnames(extract_df)) {
    stop("Collapsed extract parquet files must contain a treatment column", call. = FALSE)
  }

  if (any(is.na(extract_df$treatment))) {
    stop("Collapsed extract parquet files contain missing treatment values", call. = FALSE)
  }

  treatments <- extract_df %>%
    dplyr::distinct(treatment) %>%
    dplyr::arrange(treatment) %>%
    dplyr::pull(treatment)

  for (treatment_name in treatments) {
    combined_dir <- fs::path(out_dir, treatment_name, "05_combined")
    dir_create(combined_dir, recurse = TRUE)

    out_extract <- fs::path(
      combined_dir,
      paste0(treatment_name, "_extract_combined.parquet")
    )

    combined_extract <- extract_df %>%
      dplyr::filter(treatment == !!treatment_name)

    arrow::write_parquet(
      combined_extract,
      out_extract,
      compression = "zstd"
    )

    cat("Wrote combined extract:", out_extract, "\n")
  }
}
