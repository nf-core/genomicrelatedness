process MAKE_SCAFFOLD_BED {
    tag "$meta.id"
    label 'process_single'

    input:
    val(scaffolds)
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path("${meta.id}.bed"), emit: bed

    script:
    def regex = scaffolds.collect().join("|")
    """
    {
        echo "# BED file of mitochondrial scaffolds"
        grep -P "^(?:${regex})\\t" ${bed}
    } > "${meta.id}.bed"
    """
}
