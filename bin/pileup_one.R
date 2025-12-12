#!/home/constantinos/miniconda3/envs/nf-env/bin/Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(arrow)
  library(dplyr)
  library(tidyr)
})

option_list <- list(
  make_option("--mod-prob",type = "double", default = 0.8, help = "Probability threshold for modification calls"),
  make_option("--mod-code",type = "character",default = "m",help = "Modification code of interest (e.g. 'm' for 5mC)"),
  make_option("--in-parquet",  type = "character", help = "Input Parquet file"),
  make_option("--out-parquet", type = "character", help = "Output Parquet file"),
  make_option("--sample-id",   type = "character", help = "Sample ID",   default = NULL),
  make_option("--plasmid-id",  type = "character", help = "Plasmid ID",  default = NULL),
  make_option("--treatment",   type = "character", help = "Treatment",   default = NULL),
  make_option("--enzyme-conc", type = "character", help = "Enzyme conc", default = NULL),
  make_option("--replicate",   type = "character", help = "Replicate",   default = NULL),
  make_option("--run-id",      type = "character", help = "Run ID",      default = NULL),
  make_option("--mod-base",    type = "character", help = "Modified base for this sample (A/C)", default = NA)
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
mod_prob  <- opt$`mod-prob`
mod_code  <- opt$`mod-code`
mod_base_raw <- opt$`mod-base`

# Ensure output directory exists
out_dir <- dirname(out_pq)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

if (!is.null(mod_base_raw) && !is.na(mod_base_raw) && nzchar(mod_base_raw)) {
  mod_base_val <- toupper(mod_base_raw)
  if (!mod_base_val %in% c("A","C","G","T")) {
    stop("mod-base must be one of: A, C, G, T (got: ", mod_base_raw, ")", call. = FALSE)
  }
} else {
  mod_base_val <- NA_character_
}

# Core logic: collapse per read×position, then aggregate per position
df <- arrow::read_parquet(in_pq) %>%
  dplyr::group_by(ref_position, read_id) %>%
  dplyr::mutate(Prob_canonical = 1 - sum(mod_qual)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(mod_code = ifelse(mod_code == "canonical_placeholder", "-", mod_code)) %>%
  tidyr::pivot_longer(
    cols      = c(Prob_canonical, mod_qual),
    names_to  = ".tmp",
    values_to = "prob"
  ) %>%
  dplyr::mutate(
    mod_code = ifelse(.tmp == "Prob_canonical", "-", mod_code)
  ) %>%
  dplyr::select(-.tmp) %>%
  dplyr::distinct() %>%
  dplyr::group_by(read_id, ref_position) %>%
  dplyr::slice_max(order_by = prob, n = 1, with_ties = FALSE) %>%
  # Same logic as your original: only enforce threshold for the mod code "m"
  # (we could generalise to mod_code, but we keep semantics identical)
  dplyr::filter(mod_code != "m" | prob >= mod_prob) %>%
  # Now aggregate to site level
  dplyr::group_by(ref_position, ref_strand, ref_kmer) %>%
  dplyr::summarise(
    Nmod              = sum(mod_code == mod_code),
    Ncanonical        = sum(mod_code == "-"),
    Nother_mod        = sum(mod_code != mod_code & mod_code != "-"),
    Nvalid_cov        = Nmod + Ncanonical + Nother_mod,
    fraction_modified = Nmod / Nvalid_cov,
    .groups           = "drop"
  )

# Attach metadata
df <- df %>%
  dplyr::mutate(
    sample_id   = if (!is.null(sample_id))   sample_id   else basename(in_pq),
    plasmid_id  = if (!is.null(plasmid_id))  plasmid_id  else NA_character_,
    treatment   = if (!is.null(treatment))   treatment   else NA_character_,
    enzyme_conc = if (!is.null(enzyme_conc)) enzyme_conc else NA_character_,
    replicate   = if (!is.null(replicate))   replicate   else NA_character_,
    run_id      = if (!is.null(run_id))      run_id      else NA_character_,
    mod_base    = mod_base_val,
    prob_threshold = mod_prob

  )


# Write pileup Parquet
arrow::write_parquet(df, out_pq, compression = "zstd")

cat("Done writing", out_pq, "\n")
