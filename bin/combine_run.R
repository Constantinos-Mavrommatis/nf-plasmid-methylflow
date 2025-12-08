#!/home/constantinos/miniconda3/envs/nf-env/bin/Rscript
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(arrow)
  library(dplyr)
  library(fs)
})

option_list <- list(
  make_option(
    "--base-outdir",
    type = "character",
    help = "Base output directory of the run (e.g. results or results_run1)"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$`base-outdir`)) {
  stop("You must provide --base-outdir", call. = FALSE)
}

base_dir <- opt$`base-outdir`

if (!dir_exists(base_dir)) {
  stop("Base outdir does not exist: ", base_dir, call. = FALSE)
}

# Find treatment subdirectories: e.g. results/A, results/5mc, ...
treatment_dirs <- fs::dir_ls(
  path = base_dir,
  type = "directory"
)

if (length(treatment_dirs) == 0) {
  stop("No treatment subdirectories found under ", base_dir, call. = FALSE)
}

for (t_dir in treatment_dirs) {
  treatment <- fs::path_file(t_dir)
  message("Processing treatment: ", treatment)

  pileup_dir   <- fs::path(t_dir, "04_pileup")
  extract_dir  <- fs::path(t_dir, "03_extract")
  combined_dir <- fs::path(t_dir, "05_combined")

  dir.create(combined_dir, recursive = TRUE, showWarnings = FALSE)

  ## 1) Combine pileup for this treatment
  if (dir_exists(pileup_dir)) {
    pileup_files <- fs::dir_ls(pileup_dir, glob = "*.parquet")

    if (length(pileup_files) > 0) {
      message("  Combining pileup files from: ", pileup_dir)
      pileup_ds <- arrow::open_dataset(pileup_files, format = "parquet")
      combined_pileup <- pileup_ds %>% dplyr::collect()

      out_pileup <- fs::path(combined_dir, paste0(treatment, "_pileup_combined.parquet"))

      arrow::write_parquet(
        combined_pileup,
        out_pileup,
        compression = "zstd"
      )

      cat("  Wrote combined pileup:", out_pileup, "\n")
    } else {
      message("  No pileup .parquet files found in ", pileup_dir)
    }
  } else {
    message("  Pileup dir does not exist for treatment ", treatment, ": ", pileup_dir)
  }

  ## 2) Combine extract_collapsed for this treatment
  if (dir_exists(extract_dir)) {
    extract_files <- fs::dir_ls(extract_dir, glob = "*.parquet")

    if (length(extract_files) > 0) {
      message("  Combining extract_collapsed files from: ", extract_dir)
      extract_ds <- arrow::open_dataset(extract_files, format = "parquet")
      combined_extract <- extract_ds %>% dplyr::collect()

      out_extract <- fs::path(combined_dir, paste0(treatment, "_extract_combined.parquet"))

      arrow::write_parquet(
        combined_extract,
        out_extract,
        compression = "zstd"
      )

      cat("  Wrote combined extract:", out_extract, "\n")
    } else {
      message("  No extract .parquet files found in ", extract_dir)
    }
  } else {
    message("  Extract dir does not exist for treatment ", treatment, ": ", extract_dir)
  }
}
