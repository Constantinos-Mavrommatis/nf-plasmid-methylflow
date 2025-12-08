# nf-plasmid-methylflow

`plasmid-methylation-qc` is a **Nextflow pipeline** for analysing and QC-ing plasmid methylation from Oxford Nanopore **modBAM** files using `modkit`.

It is designed for:

- **Quality control (QC)** of plasmid methylation in routine workflows  
- **Research analysis**, e.g. testing different methylases or treatments  
- **Scalable settings**, such as **biofoundries**, where many plasmids and treatments must be processed in a consistent, automated way

The pipeline is primarily optimised for **plasmid-scale references** (small, well-defined sequences) but many components can be adapted to other small constructs.

---

## Why Nextflow?

This project started as a collection of Python and R scripts. Nextflow turns that into a **reproducible, scalable workflow**:

- ✅ **Reproducible**  
  Every step (modkit, TSV → Parquet, pileup, combining, QC) is versioned and part of a single, traceable pipeline.

- ✅ **Robust & resumable**  
  If something fails (e.g. a cluster/job glitch), you can restart with `-resume` and only the failed parts are rerun.

- ✅ **Scalable**  
  Handles **dozens to hundreds of samples** across multiple treatments in one go – ideal for a **biofoundry** setting where many plasmids and methylation conditions are tested in parallel.

- ✅ **Flexible execution**  
  You can run the same workflow:
  - on a laptop (via WSL2 or Linux),
  - on a lab server,
  - or on an HPC/cluster with a different Nextflow profile.

- ✅ **Separation of concerns**  
  The pipeline logic (in `main.nf`) is separate from:
  - the analysis scripts (`bin/`),
  - the environment (`envs/`),
  - and configuration (`nextflow.config` / profiles).

In short: **Nextflow makes this usable as a day-to-day QC tool and as a research pipeline that can grow with your experiments**.

---

## What the pipeline does

The pipeline wraps a series of steps around `modkit` + R analysis:

1. **Per-sample modkit extract**  
   - `modkit extract full` → per-read modification TSV

2. **TSV → Parquet**  
   - Faster loading and more efficient storage for downstream analysis.

3. **Per-read collapse**  
   - Collapse per-read calls to **one call per (read × position)**.  
   - Attach metadata: `sample_id`, `plasmid_id`, `treatment`, `enzyme_conc`, `replicate`, `run_id`, `mod_base`.

4. **Site-level pileup**  
   - Aggregate across reads per site to get:
     - `Nmod`, `Ncanonical`, `Nother_mod`, `Nvalid_cov`,
     - `fraction_modified` (with a configurable probability threshold).

5. **Per-treatment combining**  
   - Combine all samples of the same `treatment` (e.g. EcoGII, Dam, BamHI) into:
     - a combined **pileup** file, and  
     - a combined **collapsed-per-read** file.

6. **Optional downstream QC**  
   - `distance_methylation.R`  
     - distance vs methylation for a chosen modified base (A/C/G/T), using `mod_base` metadata.
   - `distance_motif_methylation.R`  
     - **motif-based** distance vs methylation (e.g. GGATCC, CG, CCWGG), counting:
       - total motif sites on the plasmid,
       - how many are covered,
       - how many are methylated.

This makes it useful both as:

- a **QC tool** (are my methylases doing what I expect on this plasmid?), and  
- a **research platform** (compare conditions, enzymes, constructs, etc.).

---

## Features

- **Plasmid-focused**  
  Best suited for small, well-defined plasmid constructs where motifs and positions are known and interpretable.

- **Multi-sample, multi-treatment support** via a **samplesheet** (CSV).
- Automatic grouping of output by **treatment** (e.g. BamHI, CpG, Dam, EcoGII, Dcm).
- Per-call metadata:
  - `sample_id`, `plasmid_id`, `treatment`, `mod_base`,
  - `enzyme_conc`, `replicate`, `run_id`.
- Optional **per-run combine** step:
  - per-sample Parquets → combined per-treatment Parquets.
- Extra QC / analysis scripts:
  - `distance_methylation.R`: base-level spacing vs methylation.
  - `distance_motif_methylation.R`: **motif-level** spacing vs methylation (e.g. out of 10 possible sites, 8 methylated).

---

## Installation

### Requirements

- Unix-like environment:
  - Linux, or
  - Windows with **WSL2**
- [Nextflow](https://www.nextflow.io/) ≥ 25
- Python 3 (for `modkit_extract_full.py`)
- R ≥ 4.3 with packages:
  - `optparse`
  - `arrow`
  - `dplyr`
  - `tidyr`
  - `fs`
  - `ggplot2`
  - `scales`
  - `readr`
  - `Biostrings`

You will also need `modkit` and a basecaller that produces **modBAMs** (e.g. Dorado).

### Conda environment (recommended)

```bash
conda env create -f envs/nf-env.yml
conda activate nf-env
```

This can be adapted to your lab/server setup.

---

## Samplesheet

Inputs are defined via a CSV with one row per sample, for example:

```csv
sample_id,plasmid_id,treatment,mod_base,enzyme_conc,replicate,run_id,modbam,ref_fasta
barcode03,Plasmid1,EcoGII,A,4,1,toy_run,data/barcode03_aln_sorted.bam,data/EGF_met_1.fa
barcode04,Plasmid1,Dam,A,4,1,toy_run,data/barcode04_aln_sorted.bam,data/EGF_met_1.fa
```

- `treatment`  
  Descriptive name (e.g. EcoGII, Dam, BamHI, CpG, Dcm). Often corresponds to the methylase or experimental condition.

- `mod_base`  
  Modified base for this assay (typically `A` or `C`), used by downstream QC scripts to know which base to analyse.

- `modbam`  
  Path to the **modBAM** file for that sample.

- `ref_fasta`  
  Path to the **reference FASTA** for the plasmid used in that sample.

---

## Quickstart

Run the pipeline:

```bash
nextflow run main.nf   --samplesheet samplesheet.csv   --outdir results   -profile local
```

Outputs (`--outdir results`) are grouped by `treatment`, for example:

```text
results/
  EcoGII/
    01_extract_full/         # modkit extract full TSVs
    02_parquet/              # per-read Parquets
    03_extract/              # collapsed per-read calls
    04_pileup/               # site-level pileups
    05_combined/             # combined per-treatment files
  Dam/
    01_extract_full/
    02_parquet/
    03_extract/
    04_pileup/
    05_combined/
  ...
```

This layout is friendly for:

- routine **QC** (check each treatment/plasmid),
- and later **aggregation** across treatments or runs.

---

## Example motif-based QC

For a given treatment, you can run motif-based distance QC.  
Example for EcoGII (motif `A`, modified base index `1`):

```bash
bin/distance_motif_methylation.R   --in-parquet results/EcoGII/05_combined/EcoGII_pileup_combined.parquet   --ref-fasta  data/EGF_met_1.fa   --motif      A   --mod-index  1   --out-tsv    results/EcoGII/06_distance/EcoGII_motif_distance.tsv   --out-plot   results/EcoGII/06_distance/EcoGII_motif_distance.pdf   --call-thr   0.70   --max-d      25   --min-cov    250
```

For other methylases:

- **BamHI** – motif `GGATCC`, modified A at index 3  
- **CpG** – motif `CG`, modified C at index 1  
- **Dam** – motif `GATC`, modified A at index 2  
- **Dcm** – motif `CCWGG`, modified second C (`mod-index = 2`)

The output tells you:

- how many motif sites exist on the plasmid (total possible),  
- how many are covered,  
- how many are methylated,  
- and how methylation depends on motif spacing.

---

## Use cases

- **Biofoundry / high-throughput environments**
  - Many plasmids and treatments in parallel.
  - Need reproducible, scripted QC and summary plots.
  - Easy to run the whole panel overnight and inspect per-treatment readouts in the morning.

- **Method development**
  - Compare different methylases, enzyme concentrations, or reaction conditions.
  - Investigate motif spacing effects, partial methylation, and off-target effects.

- **Routine plasmid QC**
  - Check that a plasmid has the expected methylation pattern before sending it downstream (e.g. into an expression screen).

---

## Reproducibility

- All logic is encoded in `main.nf` and `bin/` scripts.
- Runs are resumable (`-resume`) and traceable.
- Environments can be captured with:
  - conda (`envs/nf-env.yml`), or

---

## Citation

MIT
