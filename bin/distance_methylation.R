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
    help = "Input combined pileup Parquet file"
  ),
  make_option(
    "--out-tsv",
    type = "character",
    help = "Output TSV file with distance vs methylation summary"
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
    help = "Maximum nearest-A distance to consider [default %default]"
  ),
  make_option(
    "--min-cov",
    type = "integer",
    default = 250L,
    help = "Minimum Nvalid_cov per site [default %default]"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$`in-parquet`) || is.null(opt$`out-tsv`) || is.null(opt$`out-plot`)) {
  stop("You must provide --in-parquet, --out-tsv and --out-plot", call. = FALSE)
}

in_pq      <- opt$`in-parquet`
out_tsv    <- opt$`out-tsv`
out_plot   <- opt$`out-plot`
call_thr   <- opt$`call-thr`
max_d      <- opt$`max-d`
min_cov    <- opt$`min-cov`

## ---- Load & filter pileup data ----
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

## ---- K-mer sanity check & center index ----
k <- unique(nchar(df$ref_kmer))
stopifnot(length(k) == 1L, k %% 2L == 1L)  # one odd k-mer length

mid_idx <- (k + 1L) %/% 2L

## ---- Restrict to sites where center base is 'A' ----
df_A <- df %>%
  dplyr::filter(substr(ref_kmer, mid_idx, mid_idx) == "A") %>%
  dplyr::group_by(sample_id, ref_position) %>%  # ensure uniqueness per site
  dplyr::summarise(
    fraction_modified = mean(fraction_modified, na.rm = TRUE),
    .groups = "drop"
  )

if (nrow(df_A) == 0L) {
  stop("No A-centered positions after filtering. Check your ref_kmer or filters.", call. = FALSE)
}

## ---- Compute nearest-A distance per site ----
df_nn <- df_A %>%
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
    nearest_d    = dplyr::if_else(use_fwd, dist_fwd, dist_back),
    current_call = fraction_modified >= call_thr
  ) %>%
  dplyr::ungroup() %>%
  dplyr::filter(!is.na(nearest_d), nearest_d >= 1) %>%
  dplyr::mutate(nearest_d = as.integer(nearest_d)) %>%
  dplyr::filter(nearest_d <= max_d)

if (nrow(df_nn) == 0L) {
  stop("No sites left after computing nearest_d and applying max_d filter. ",
       "Try lowering max_d or min_cov.", call. = FALSE)
}

## ---- Per-sample, per-distance stats ----
by_sample <- df_nn %>%
  dplyr::group_by(sample_id, nearest_d) %>%
  dplyr::summarise(
    n_sites = dplyr::n(),
    frac_current_methylated = mean(current_call, na.rm = TRUE),
    .groups = "drop"
  )

## ---- Save summary as TSV ----
dir.create(dirname(out_tsv), recursive = TRUE, showWarnings = FALSE)

readr::write_tsv(by_sample, out_tsv)
cat("Wrote distance vs methylation summary TSV:", out_tsv, "\n")

## ---- Make and save plot ----
dir.create(dirname(out_plot), recursive = TRUE, showWarnings = FALSE)

p <- ggplot(by_sample, aes(nearest_d, frac_current_methylated, group = sample_id)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ sample_id) +
  scale_x_continuous(breaks = seq_len(max_d)) +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = paste0("Methylated positions (≥ ", call_thr, ") vs nearest-A distance"),
    x = "Distance to nearest A (bp)",
    y = "Fraction of sites methylated"
  )

ggsave(out_plot, plot = p, width = 8, height = 6)
cat("Wrote distance vs methylation plot:", out_plot, "\n")

