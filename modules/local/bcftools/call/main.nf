process BCFTOOLS_CALL {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/47/474a5ea8dc03366b04df884d89aeacc4f8e6d1ad92266888e7a8e7958d07cde8/data':
        'community.wave.seqera.io/library/bcftools_htslib:0a3fa2654b52006f' }"

    input:
    tuple val(meta), path(vcf), path(intervals)

    output:
    tuple val(meta), path("${prefix}.${extension}")    , emit: vcf, optional: true
    tuple val(meta), path("${prefix}.${extension}.tbi"), emit: tbi, optional: true
    tuple val(meta), path("${prefix}.${extension}.csi"), emit: csi, optional: true
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | sed '1!d; s/^.*bcftools //'"), topic: versions, emit: versions_bcftools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    def intervalsOpt = intervals ? "-T ${intervals}" : ""

    extension = args.contains("--output-type b") || args.contains("-Ob") || args.contains("-O b") ? "bcf.gz" :
                args.contains("--output-type u") || args.contains("-Ou") || args.contains("-O u") ? "bcf" :
                args.contains("--output-type z") || args.contains("-Oz") || args.contains("-O z") ? "vcf.gz" :
                args.contains("--output-type v") || args.contains("-Ov") || args.contains("-O v") ? "vcf" :
                "vcf"

    """
    bcftools call \
        -o ${prefix}.${extension} \
        ${intervalsOpt} \
        ${args} \
        ${vcf}
    """
}
