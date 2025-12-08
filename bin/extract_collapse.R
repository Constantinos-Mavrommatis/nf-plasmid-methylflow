#!/usr/bin/Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(arrow)
  library(dplyr)
  library(tidyr)
})

option_list <- list(
  make_option("--in-parquet",  type = "character", help = "Input Parquet file"),
  make_option("--out-parquet", type = "character", help = "Output Parquet file"),
  make_option("--sample-id",   type = "character", help = "Sample ID",   default = NULL),
  make_option("--plasmid-id",  type = "character", help = "Plasmid ID",  default = NULL),
  make_option("--treatment",   type = "character", help = "Treatment",   default = NULL),
  make_option("--enzyme-conc", type = "character", help = "Enzyme conc", default = NULL),
  make_option("--replicate",   type = "character", help = "Replicate",   default = NULL),
  make_option("--run-id",      type = "character", help = "Run ID",      default = NULL)
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$`in-parquet`) || is.null(opt$`out-parquet`)) {
  stop("Both --in-parquet and --out-parquet must be provided", call. = FALSE)
}

in_pq     <- opt$`in-parquet`
out_pq    <- opt$`out-parquet`
sample_id   <- opt$`sample-id`
plasmid_id  <- opt$`plasmid-id`
treatment   <- opt$treatment
enzyme_conc <- opt$`enzyme-conc`
replicate   <- opt$replicate
run_id      <- opt$`run-id`

# Ensure output directory exists
out_dir <- dirname(out_pq)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

# Read the per-read Parquet
df <- arrow::read_parquet(in_pq) %>%
  dplyr::group_by(ref_position, read_id) %>%
  # probability mass for canonical = 1 - sum(mod_qual) over all mods at that site
  dplyr::mutate(Prob_canonical = 1 - sum(mod_qual)) %>%
  dplyr::ungroup() %>%
  # normalise placeholder for canonical
  dplyr::mutate(mod_code = ifelse(mod_code == "canonical_placeholder", "-", mod_code)) %>%
  # now pivot canonical vs modification probability into long format
  tidyr::pivot_longer(
    cols      = c(Prob_canonical, mod_qual),
    names_to  = ".tmp",
    values_to = "prob"
  ) %>%
  dplyr::mutate(
    # when the row comes from Prob_canonical, treat mod_code as '-'
    mod_code = ifelse(.tmp == "Prob_canonical", "-", mod_code)
  ) %>%
  dplyr::select(-.tmp) %>%
  dplyr::distinct() %>%
  # choose the highest probability outcome (canonical vs modified) per read×position
  dplyr::group_by(read_id, ref_position) %>%
  dplyr::slice_max(order_by = prob, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup()

# Attach metadata
df <- df %>%
  dplyr::mutate(
    sample_id   = if (!is.null(sample_id))   sample_id   else basename(in_pq),
    plasmid_id  = if (!is.null(plasmid_id))  plasmid_id  else NA_character_,
    treatment   = if (!is.null(treatment))   treatment   else NA_character_,
    enzyme_conc = if (!is.null(enzyme_conc)) enzyme_conc else NA_character_,
    replicate   = if (!is.null(replicate))   replicate   else NA_character_,
    run_id      = if (!is.null(run_id))      run_id      else NA_character_
  )

# Write collapsed Parquet
arrow::write_parquet(df, out_pq, compression = "zstd")

cat("Done writing", out_pq, "\n")
