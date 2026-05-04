nextflow.enable.dsl = 2

// Define parameter defaults without clobbering explicit false/zero values
if( params.samplesheet == null ) params.samplesheet = 'samplesheet.csv'
if( params.outdir == null )      params.outdir      = 'results'
if( params.mod_prob == null )    params.mod_prob    = 0.5
if( params.mod_code == null )    params.mod_code    = 'a'
if( params.combine_run == null ) params.combine_run = true

// Optional distance analysis
if( params.do_distance == null ) params.do_distance = false
if( params.call_thr == null )    params.call_thr    = 0.7
if( params.max_d == null )       params.max_d       = 25
if( params.min_cov == null )     params.min_cov     = 250

include { DISTANCE_KMER_MOTIF } from './modules/distance_kmer_motif.nf'

def outdir_abs = new File(params.outdir).getAbsolutePath()

workflow {

    /*
     * 1. Read samplesheet.csv and turn each row into:
     *    (meta, bam, ref)
     */
    Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row ->
            def meta = [
                sample_id   : row.sample_id,
                plasmid_id  : row.plasmid_id,
                treatment   : row.treatment,
                mod_base    : row.mod_base,
                enzyme_conc : row.enzyme_conc,
                replicate   : row.replicate,
                run_id      : row.run_id
            ]

            def bam_path = row.modbam
            def ref_path = row.ref_fasta

            tuple(
                meta,
                file(bam_path),
                file(ref_path)
            )
        }
        .set { sample_ch }

    /*
     * 2. Run modkit_extract_full.py for each sample → extract TSV
     */
    extract_ch = EXTRACT_FULL(sample_ch)

    /*
     * 3. Convert TSV → Parquet for each sample
     */
    parquet_ch = TSV_TO_PARQUET(extract_ch)

    /* 
     * 4b. Per-read Parquet -> collapsed per-read calls
     */
    collapsed_ch = EXTRACT_COLLAPSED(parquet_ch)

    /*
     * 4a. Convert TSV → Parquet for each sample
     */
    pileup_ch = PILEUP_R(parquet_ch)

    /*
     * 6. Optional per-run combine step
     *    We "collect" the channels so COMBINE_RUN only runs once,
     *    after all per-sample tasks are done.
     */
    if( params.combine_run ) {
        def pileup_files_ch  = pileup_ch
            .map { meta, pileup, ref -> pileup }
            .collect()
        def extract_files_ch = collapsed_ch
            .map { meta, parquet, ref -> parquet }
            .collect()

        COMBINE_RUN( pileup_files_ch, extract_files_ch )
    }

    // 7. OPTIONAL: run distance QC per sample, if requested
    if( params.do_distance ) {
        distance_ch = DISTANCE_KMER_MOTIF(pileup_ch)
    }

}

/*
 * Process: EXTRACT_FULL
 * One sample in → one extract-full TSV out
 */
process EXTRACT_FULL {

    tag "${meta.sample_id}"
    publishDir "${params.outdir}/${meta.treatment}/01_extract_full", mode: 'copy'

    input:
    tuple val(meta), path(modbam), path(ref)

    output:
    tuple val(meta), path("${meta.sample_id}.extract_full.tsv"), path(ref)

    script:
    """
    modkit_extract_full.py \
      --in-bam ${modbam} \
      --ref ${ref} \
      --out-tsv ${meta.sample_id}.extract_full.tsv
    """
}

process TSV_TO_PARQUET {

    tag "${meta.sample_id}"
    publishDir "${params.outdir}/${meta.treatment}/02_parquet", mode: 'copy'

    input:
    tuple val(meta), path(tsv), path(ref)

    output:
    // carry meta + ref forward, but now with a Parquet file
    tuple val(meta), path("${meta.sample_id}.parquet"), path(ref)

    script:
    """
    convert_one.R \
      --in-tsv ${tsv} \
      --out-parquet ${meta.sample_id}.parquet
    """
}

process EXTRACT_COLLAPSED {

    tag "${meta.sample_id}"
    publishDir "${params.outdir}/${meta.treatment}/03_extract", mode: 'copy'

    input:
    // matches the output of TSV_TO_PARQUET: (meta, <parquet>, ref)
    tuple val(meta), path(parquet), path(ref)

    output:
    // carry meta + ref forward, but now with collapsed Parquet
    tuple val(meta), path("${meta.sample_id}.extract_collapsed.parquet"), path(ref)

    script:
    """
    extract_collapse.R \
    --in-parquet ${parquet} \
    --out-parquet ${meta.sample_id}.extract_collapsed.parquet \
    --sample-id ${meta.sample_id} \
    --plasmid-id ${meta.plasmid_id} \
    --treatment ${meta.treatment} \
    --enzyme-conc ${meta.enzyme_conc} \
    --replicate ${meta.replicate} \
    --run-id ${meta.run_id} \
    --mod-base ${meta.mod_base}
    """
}

process PILEUP_R {

    tag "${meta.sample_id}"
    publishDir "${params.outdir}/${meta.treatment}/04_pileup/${params.mod_prob}_modification_threshold/", mode: 'copy'

    input:
    // matches the output of TSV_TO_PARQUET: (meta, <per-read.parquet>, ref)
    tuple val(meta), path(parquet), path(ref)

    output:
    // carry meta + ref forward with a pileup Parquet
    tuple val(meta), path("${meta.sample_id}_${params.mod_prob}.pileup.parquet"), path(ref)

    script:
    """
    pileup_one.R \
    --in-parquet ${parquet} \
    --out-parquet ${meta.sample_id}_${params.mod_prob}.pileup.parquet \
    --mod-prob ${params.mod_prob} \
    --mod-code ${params.mod_code} \
    --sample-id ${meta.sample_id} \
    --plasmid-id ${meta.plasmid_id} \
    --treatment ${meta.treatment} \
    --enzyme-conc ${meta.enzyme_conc} \
    --replicate ${meta.replicate} \
    --run-id ${meta.run_id} \
    --mod-base ${meta.mod_base}
    """
}

/*
 * Process: COMBINE_RUN
 * Runs ONCE at the end of the workflow.
 * For each treatment folder under ${params.outdir}:
 *   - combine 03_extract/*.parquet → 05_combined/<treatment>_extract_combined.parquet
 *   - combine 04_pileup/*.parquet → 05_combined/<treatment>_pileup_combined.parquet
 */
process COMBINE_RUN {

    tag "combine_run"

    input:
    // Use the actual upstream parquet outputs so combine depends on real workflow data
    val pileup_files
    val extract_files

    output:
    // a marker file so Nextflow knows this step ran
    path "combine_run.done"

    script:
    def pileup_manifest = pileup_files.collect { it.toString() }.join('\n')
    def extract_manifest = extract_files.collect { it.toString() }.join('\n')
    """
    printf '%s\n' "${pileup_manifest}" > pileup_files.txt
    printf '%s\n' "${extract_manifest}" > extract_files.txt

    combine_run.R \
      --outdir ${outdir_abs} \
      --pileup-list pileup_files.txt \
      --extract-list extract_files.txt

    touch combine_run.done
    """
}
