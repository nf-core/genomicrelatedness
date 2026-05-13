process BCFTOOLS_MPILEUP {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0b/0b4d52ca9a56d07be3f78a12af654e5116f5112908dba277e6796fd9dfb83fe5/data'
        : 'community.wave.seqera.io/library/bcftools_htslib:1.23.1--9f08ec665533d64a'}"

    input:
    tuple val(meta),  path(intervals), path(bam)
    tuple val(meta2), path(fasta)
    val save_mpileup   // boolean

    output:
    tuple val(meta), path("${prefix}.${extension}")    , emit: vcf, optional: true
    tuple val(meta), path("${prefix}.${extension}.tbi"), emit: tbi, optional: true
    tuple val(meta), path("${prefix}.${extension}.csi"), emit: csi, optional: true
    tuple val(meta), path("${meta.id}.mpileup.gz")     , emit: raw_mpileup, optional: true
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | sed '1!d; s/^.*bcftools //'"), topic: versions, emit: versions_bcftools


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    def intervalsOpt = intervals ? "-T ${intervals}" : ""
    //def bams = bam.findAll { file -> !(file instanceof List) }.collect()

    // Optional saving of the raw textual mpileup
    def raw_mpileup_cmd = save_mpileup ? "| tee ${prefix}.mpileup" : ""
    def compress_raw_mpileup = save_mpileup ? "bgzip -f ${prefix}.mpileup" : ""

    extension = args.contains("--output-type b") || args.contains("-Ob") || args.contains("-O b") ? "bcf.gz" :
                args.contains("--output-type u") || args.contains("-Ou") || args.contains("-O u") ? "bcf" :
                args.contains("--output-type z") || args.contains("-Oz") || args.contains("-O z") ? "vcf.gz" :
                args.contains("--output-type v") || args.contains("-Ov") || args.contains("-O v") ? "vcf" :
                "vcf"

    """
    bcftools mpileup \
        --fasta-re ${fasta} \
        -o ${prefix}.${extension} \
        ${intervalsOpt} \
        ${args} \
        ${bam} \
        ${raw_mpileup_cmd}

    ${compress_raw_mpileup}
    """
}
