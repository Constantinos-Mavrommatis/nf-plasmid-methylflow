#!/usr/bin/Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(arrow)
  library(vroom)
})

option_list <- list(
  make_option(
    "--in-tsv",
    type = "character",
    help = "Input TSV file from modkit extract full"
  ),
  make_option(
    "--out-parquet",
    type = "character",
    help = "Output Parquet file path"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$`in-tsv`) || is.null(opt$`out-parquet`)) {
  stop("Both --in-tsv and --out-parquet must be provided", call. = FALSE)
}

in_tsv  <- opt$`in-tsv`
out_pq  <- opt$`out-parquet`

# Ensure output directory exists
out_dir <- dirname(out_pq)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

# Read TSV (only needed columns)
tbl <- vroom::vroom(
  file          = in_tsv,
  delim         = "\t",
  col_select    = c(
    read_id,
    ref_position,
    ref_strand,
    read_length,
    mod_qual,
    mod_code,
    base_qual,
    ref_kmer
  ),
  progress      = FALSE,
  show_col_types = FALSE
)

# Write Parquet with compression
arrow::write_parquet(tbl, out_pq, compression = "zstd")

cat("Done writing", out_pq, "\n")
