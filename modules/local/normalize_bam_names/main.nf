process NORMALIZE_BAM_NAMES {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.id}.bam"), emit: bam

    script:
    """
    ln -s ${bam} ${meta.id}.bam
    """
}
