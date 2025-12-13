#!/home/constantinos/miniconda3/envs/nf-env/bin/Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(arrow)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(scales)
  library(ggrepel)
})

option_list <- list(
  make_option(
    "--in-dir",
    type = "character",
    help = "Directory containing combined pileup Parquet files for all thresholds"
  ),
  make_option(
    "--positions-tsv",
    type = "character",
    help = "TSV file with a column 'ref_position' listing validated methylated positions"
  ),
  make_option(
    "--out-tsv",
    type = "character",
    default = "positional_metrics_multi_threshold.tsv",
    help = "Output TSV with per-sample per-threshold metrics [default %default]"
  ),
  make_option(
    "--out-plot",
    type = "character",
    default = "positional_metrics_multi_threshold.pdf",
    help = "Output PDF plot [default %default]"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$`in-dir`) || is.null(opt$`positions-tsv`)) {
  stop("You must provide --in-dir and --positions-tsv", call. = FALSE)
}

in_dir       <- opt$`in-dir`
positions_tsv <- opt$`positions-tsv`
out_tsv      <- opt$`out-tsv`
out_plot     <- opt$`out-plot`

## 1. Read validated positions
pos_df <- readr::read_csv(positions_tsv, show_col_types = FALSE)
if (!all(c("ref_position", "truth_status") %in% colnames(pos_df))) {
  stop(
    "positions-tsv must contain columns 'ref_position' and 'truth_status'",
    call. = FALSE
  )
}

pos <- unique(pos_df$ref_position)

## 2. Read all combined Parquet files (possibly multiple thresholds)
cat("Reading combined Parquet files from:", in_dir, "\n")
ds <- arrow::open_dataset(in_dir)

needed_cols <- c("sample_id", "ref_position", "Nmod", "Ncanonical", "Nother_mod", "prob_threshold")
if (!all(needed_cols %in% colnames(ds))) {
  stop(
    "Input dataset must contain columns: ",
    paste(needed_cols, collapse = ", "),
    call. = FALSE
  )
}

df <- ds %>%
  filter(ref_position %in% pos) %>%   # optional but efficient
  collect() %>%
  inner_join(pos_df, by = "ref_position")

if (nrow(df) == 0L) {
  stop("No rows found at the validated positions in the provided dataset.", call. = FALSE)
}

## 3. Compute TP, FN, total_reads per sample & prob_threshold
metrics_meth <- df %>%
  filter(truth_status == "Methylated") %>%
  group_by(sample_id, prob_threshold, ref_position) %>%
  summarise(
    TP          = sum(Nmod),                               # called mod at mod sites
    FN          = sum(Ncanonical + Nother_mod),            # anything not mod
    total_reads = sum(Nmod + Ncanonical + Nother_mod),
    truth_status = "Methylated",
    .groups     = "drop"
  )

metrics_can <- df %>%
  filter(truth_status == "Canonical") %>%
  group_by(sample_id, prob_threshold, ref_position) %>%
  summarise(
    TP          = sum(Ncanonical),                         # called canonical at canonical sites
    FN          = sum(Nmod + Nother_mod),                  # called modified at canonical sites
    total_reads = sum(Nmod + Ncanonical + Nother_mod),
    truth_status = "Canonical",
    .groups     = "drop"
  )

## 4. Compute per-sample baseline: max total_reads across thresholds

metrics <- bind_rows(metrics_meth, metrics_can)

# Baseline per sample + truth_status (so methylated and canonical each get their own baseline)
baseline <- metrics %>%
  group_by(sample_id, truth_status) %>%
  summarise(
    baseline_reads = max(total_reads, na.rm = TRUE),
    .groups        = "drop"
  )

metrics2 <- metrics %>%
  left_join(baseline, by = c("sample_id", "truth_status")) %>%
  mutate(
    recall = ifelse(TP + FN > 0, TP / (TP + FN), NA_real_),
    retained_fraction = ifelse(
      baseline_reads > 0,
      total_reads / baseline_reads,
      NA_real_
    )
  )

## 5. Save TSV
readr::write_tsv(metrics2, out_tsv)
cat("Wrote metrics to:", out_tsv, "\n")

## 6. Plot recall vs retained_fraction
p <- ggplot(
  metrics2,
  aes(
    x = retained_fraction,
    y = recall,
    color = factor(prob_threshold),
    group = interaction(ref_position, truth_status),
    label = ref_position
  )
) +
  geom_point(alpha = 0.7) +
  geom_text_repel(size = 2.5, max.overlaps = 25) +
  facet_grid(truth_status ~ sample_id) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title = "Per-Position Recall vs Retained Fraction at Validated Positions",
    x = "Retained Fraction of Reads",
    y = "Recall",
    color = "Prob Threshold"
  ) +
  theme_bw()

ggsave(out_plot, p, width = 8, height = 6)
cat("Wrote plot to:", out_plot, "\n")