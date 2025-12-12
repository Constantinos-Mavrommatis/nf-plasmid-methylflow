nextflow.enable.dsl = 2

process DISTANCE_KMER_MOTIF {

    tag "${meta.sample_id}"
    publishDir "${params.outdir}/${meta.treatment}/06_distance/${meta.mod_prob}", mode: 'copy'

    input:
    // Re-use the same tuple type as PILEUP_R output:
    // (meta, pileup_parquet, ref)
    tuple val(meta), path(pileup), path(ref)

    output:
    // Keep meta so we can chain later if we want
    tuple val(meta), path("${meta.sample_id}_${meta.mod_prob}.distance.tsv"), path("${meta.sample_id}_${meta.mod_prob}.distance.pdf")

    script:
    """
    distance_methylation.R \
      --in-parquet ${pileup} \
      --enzyme ${meta.treatment} \
      --out-tsv  ${meta.sample_id}_${meta.mod_prob}.distance.tsv \
      --out-plot ${meta.sample_id}_${meta.mod_prob}_.distance.pdf \
      --call-thr ${params.call_thr} \
      --max-d    ${params.max_d} \
      --min-cov  ${params.min_cov}
    """
}

