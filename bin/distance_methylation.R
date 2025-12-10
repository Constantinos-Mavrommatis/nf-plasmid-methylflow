#!/home/constantinos/miniconda3/envs/nf-env/bin/Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(arrow)
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(readr)
})

## ---- CLI options ----
option_list <- list(
  make_option(
    "--in-parquet",
    type = "character",
    help = "Input combined pileup Parquet file (per treatment)"
  ),
  make_option(
    "--enzyme",
    type = "character",
    help = "Treatment / enzyme name (e.g. Dam, CpG, EcoGII)"
  ),
  make_option(
    "--out-tsv",
    type = "character",
    help = "Output TSV: distance vs methylation summary"
  ),
  make_option(
    "--out-plot",
    type = "character",
    help = "Output plot file (PDF or PNG)"
  ),
  make_option(
    "--call-thr",
    type = "double",
    default = 0.70,
    help = "Methylation calling threshold on fraction_modified [default %default]"
  ),
  make_option(
    "--max-d",
    type = "integer",
    default = 25L,
    help = "Maximum nearest-motif distance to consider [default %default]"
  ),
  make_option(
    "--min-cov",
    type = "integer",
    default = 250L,
    help = "Minimum Nvalid_cov per site [default %default]"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$`in-parquet`) ||
    is.null(opt$enzyme)       ||
    is.null(opt$`out-tsv`)    ||
    is.null(opt$`out-plot`)) {
  stop("You must provide --in-parquet, --enzyme, --out-tsv and --out-plot", call. = FALSE)
}

in_pq    <- opt$`in-parquet`
enzyme   <- opt$enzyme
out_tsv  <- opt$`out-tsv`
out_plot <- opt$`out-plot`
call_thr <- opt$`call-thr`
max_d    <- opt$`max-d`
min_cov  <- opt$`min-cov`

enzyme <- trimws(enzyme)

## ---- 1. Load pileup and filter by coverage ----
message("Reading pileup from: ", in_pq)

df <- arrow::read_parquet(in_pq) %>%
  dplyr::filter(Nvalid_cov >= min_cov)

needed_cols <- c("ref_kmer", "ref_position", "fraction_modified", "sample_id", "Nvalid_cov")
if (!all(needed_cols %in% colnames(df))) {
  stop(
    "Input parquet must contain columns: ",
    paste(needed_cols, collapse = ", "),
    call. = FALSE
  )
}

## ---- 2. Determine k-mer length and center index ----
k <- unique(nchar(df$ref_kmer))
if (length(k) != 1L || k %% 2L != 1L) {
  stop("ref_kmer must have a single odd length across all rows (e.g. k=5)", call. = FALSE)
}
k <- k[[1]]
mid_idx <- (k + 1L) %/% 2L

## ---- 3. Mark motif sites using ref_kmer ----
center  <- substr(df$ref_kmer, mid_idx, mid_idx)

is_motif <- FALSE

if (enzyme == "Dam") {
  # motif GATC, mod on A at center; 5-mer should be xGATC
  motif_kmer <- "GATC"
  is_motif <- substr(df$ref_kmer, mid_idx-1, mid_idx+2) == motif_kmer

} else if (enzyme == "CpG") {
  # motif CG, mod on C at center; 2-mer starting at center should be CG
  motif_kmer <- "CG"
  is_motif <- substr(df$ref_kmer, mid_idx, mid_idx+1) == motif_kmer

} else if (enzyme == "EcoGII") {
  # motif is just A; any A at center
  motif_kmer <- "A"
  is_motif <- center == "A"

} else if (enzyme %in% c("BamHI", "Dcm")) {
  stop("Enzyme ", enzyme,
       " not yet implemented with ref_kmer-only logic; treat as a limitation for now.",
       call. = FALSE)

} else {
  stop("Unknown enzyme '", enzyme,
       "'. Supported so far: Dam, CpG, EcoGII (BamHI/Dcm to add later).",
       call. = FALSE)
}

df_motif <- df %>%
  dplyr::filter(is_motif)

if (nrow(df_motif) == 0L) {
  stop("No motif-like sites found in ref_kmer for enzyme ", enzyme,
       " after coverage filter. Check k-mer definition or enzyme name.",
       call. = FALSE)
}

## ---- 4. Total possible motif sites (within data) ----
# We treat "total possible" as all distinct positions with the motif pattern
# that appear in the parquet after coverage filtering.
motif_positions <- df_motif %>%
  dplyr::distinct(ref_position)

n_motif_total <- nrow(motif_positions)
message("Found ", n_motif_total,
        " motif-like positions for enzyme ", enzyme,
        " in the pileup data (k-mer based).")

## ---- 5. Per-sample coverage & methylation at motif sites ----
df_motif <- df_motif %>%
  dplyr::mutate(
    current_call = fraction_modified >= call_thr
  )

coverage_summary <- df_motif %>%
  dplyr::group_by(sample_id) %>%
  dplyr::summarise(
    n_motif_covered    = dplyr::n_distinct(ref_position),
    n_motif_methylated = dplyr::n_distinct(ref_position[current_call]),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    n_motif_total = n_motif_total
  )

## ---- 6. Nearest-motif distance per site (per sample) ----
df_nn <- df_motif %>%
  dplyr::group_by(sample_id) %>%
  dplyr::arrange(ref_position, .by_group = TRUE) %>%
  dplyr::mutate(
    dist_fwd  = dplyr::lead(ref_position) - ref_position,
    dist_back = ref_position - dplyr::lag(ref_position),
    use_fwd   = dplyr::case_when(
      is.na(dist_back)      ~ TRUE,
      is.na(dist_fwd)       ~ FALSE,
      dist_fwd <= dist_back ~ TRUE,
      TRUE                  ~ FALSE
    ),
    nearest_d = dplyr::if_else(use_fwd, dist_fwd, dist_back)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::filter(!is.na(nearest_d), nearest_d >= 1) %>%
  dplyr::mutate(nearest_d = as.integer(nearest_d)) %>%
  dplyr::filter(nearest_d <= max_d)

if (nrow(df_nn) == 0L) {
  stop(
    "No motif sites left after computing nearest_d and applying max_d filter. ",
    "Try lowering max_d or min_cov.",
    call. = FALSE
  )
}

## ---- 7. Per-sample, per-distance stats ----
by_sample <- df_nn %>%
  dplyr::group_by(sample_id, nearest_d) %>%
  dplyr::summarise(
    n_sites = dplyr::n(),
    frac_current_methylated = mean(current_call, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::left_join(coverage_summary, by = "sample_id") %>%
  dplyr::mutate(
    enzyme     = enzyme,
    call_thr   = call_thr,
    min_cov    = min_cov,
    max_d_used = max_d,
    kmer_len   = k
  )

## ---- 8. Save summary as TSV ----
dir.create(dirname(out_tsv), recursive = TRUE, showWarnings = FALSE)
readr::write_tsv(by_sample, out_tsv)

cat("Wrote k-mer-based motif distance vs methylation summary TSV:", out_tsv, "\n")

## ---- 9. Make and save plot ----
dir.create(dirname(out_plot), recursive = TRUE, showWarnings = FALSE)

p <- ggplot(by_sample, aes(nearest_d, frac_current_methylated, group = sample_id)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ sample_id) +
  scale_x_continuous(breaks = seq_len(max_d)) +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = paste0(
      "Methylated motif-like sites (≥ ", call_thr, ") vs nearest-motif distance\n",
      "Enzyme: ", enzyme,
      " | total motif-like sites (in data): ", n_motif_total
    ),
    x = "Distance to nearest motif-like site (bp)",
    y = "Fraction of motif-like sites methylated"
  )

ggsave(out_plot, plot = p, width = 8, height = 6)
cat("Wrote k-mer-based motif distance vs methylation plot:", out_plot, "\n")


