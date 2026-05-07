process BCFTOOLS_MPILEUP {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/47/474a5ea8dc03366b04df884d89aeacc4f8e6d1ad92266888e7a8e7958d07cde8/data':
        'community.wave.seqera.io/library/bcftools_htslib:0a3fa2654b52006f' }"

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
