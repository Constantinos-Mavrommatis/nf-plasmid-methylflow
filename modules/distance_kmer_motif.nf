nextflow.enable.dsl = 2

process DISTANCE_KMER_MOTIF {

    tag "${meta.sample_id}"
    publishDir "${params.outdir}/${meta.treatment}/06_distance", mode: 'copy'

    input:
    // Re-use the same tuple type as PILEUP_R output:
    // (meta, pileup_parquet, ref)
    tuple val(meta), path(pileup), path(ref)

    output:
    // Keep meta so we can chain later if we want
    tuple val(meta), path("${meta.sample_id}.distance.tsv"), path("${meta.sample_id}.distance.pdf")

    script:
    """
    distance_methylation.R \
      --in-parquet ${pileup} \
      --enzyme ${meta.treatment} \
      --out-tsv  ${meta.sample_id}.distance.tsv \
      --out-plot ${meta.sample_id}.distance.pdf \
      --call-thr ${params.call_thr} \
      --max-d    ${params.max_d} \
      --min-cov  ${params.min_cov}
    """
}

