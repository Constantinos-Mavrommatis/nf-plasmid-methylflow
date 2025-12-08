# plasmid-methylation-qc

Nextflow-based pipeline for **QC of plasmid methylation** using Oxford Nanopore **modBAMs** and `modkit`.

The pipeline wraps:

- `modkit extract full` → per-read modification TSV
- TSV → Parquet conversion
- Per-read collapse to **one call per read×position**
- Site-level pileup with filtering by modification probability
- Per-treatment combining across samples
- Optional distance-based QC (e.g. methylation vs spacing between modified bases or motifs)

## Features

- Supports multiple samples via a **samplesheet** (CSV)
- Groups output by **treatment** (e.g. BamHI, CpG, Dam, EcoGII, Dcm)
- Stores metadata per call:
  - `sample_id`, `plasmid_id`, `treatment`, `enzyme_conc`,
  - `replicate`, `run_id`, `mod_base`
- Optional per-run combine step:
  - combines per-sample Parquets into per-treatment Parquets
- Extra QC scripts:
  - `distance_methylation.R` – distance vs methylation for chosen base
  - `distance_motif_methylation.R` – **motif-based** distance vs methylation
    (e.g. GGATCC, CG, CCWGG…)

## Installation

Requirements:

- Unix-like environment (Linux or WSL2)
- [Nextflow](https://www.nextflow.io/) ≥ 25
- Python 3 (for `modkit_extract_full.py`)
- R ≥ 4.3 with packages:
  - `optparse, arrow, dplyr, tidyr, fs, ggplot2, scales, readr, Biostrings`

Optionally, create a conda env:

```bash
conda env create -f envs/nf-env.yml
conda activate nf-env

## Samplesheet

Input samples are defined via a CSV:

sample_id,plasmid_id,treatment,mod_base,enzyme_conc,replicate,run_id,modbam,ref_fasta
barcode03,Plasmid1,EcoGII,A,4,1,toy_run,data/barcode03_aln_sorted.bam,data/EGF_met_1.fa
barcode04,Plasmid1,Dam,A,4,1,toy_run,data/barcode04_aln_sorted.bam,data/EGF_met_1.fa

treatment – descriptive name (EcoGII, Dam, BamHI, CpG, Dcm, etc.)

mod_base – modified base for this assay (A or C typically)

modbam – path to modBAM file

ref_fasta – path to reference FASTA for that plasmid


## Quickstart

nextflow run main.nf \
  --samplesheet samplesheet.csv \
  --outdir results \
  -profile local

Outputs (--outdir results) are grouped by treatment:

results/
  EcoGII/
    01_extract_full/         # modkit extract full TSVs
    02_parquet/              # per-read Parquets
    03_extract/              # collapsed per-read calls
    04_pileup/               # site-level pileups
    05_combined/             # combined per-treatment files
  Dam/
    ...

## Reproducibility

All processes run via Nextflow (main.nf).

Tools and R packages are pinned in envs/nf-env.yml (or containers, if you add them).

