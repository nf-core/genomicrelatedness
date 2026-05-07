process ANGSD_NGSRELATE {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/ngsrelate:2.0--hea85c65_0':
          'biocontainers/ngsrelate:2.0--hea85c65_0' }"

    input:
    tuple val(meta), path(vcf), path(sample_names)

    output:
    tuple val(meta), path("*.res"), emit: result
    tuple val("${task.process}"), val('ngsrelate'), eval('echo "2.0"'), topic: versions, emit: versions_ngsrelate

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix   = task.ext.prefix ?: "${meta.id}"
    def args_list = args.tokenize()
    def forbidden_args = ["-h", "-O", "-p", "-z"].intersect(args_list)

    if (forbidden_args) {
      error "ANGSD_NGSRELATE: Reserved arguments found in task.ext.args (${forbidden_args.join(', ')}). The module sets -h, -O, -p, and -z automatically."
    }

    def arg_sample_names = sample_names ? "-z ${sample_names}" : ""

    """
    ngsRelate \\
      -p ${task.cpus} \\
      -h ${vcf} \\
      -O ${prefix}.res \\
      ${arg_sample_names} \\
      ${args}
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix   = task.ext.prefix ?: "${meta.id}"
    def args_list = args.tokenize()
    def forbidden_args = ["-h", "-O", "-p", "-z"].intersect(args_list)

    if (forbidden_args) {
      error "ANGSD_NGSRELATE: Reserved arguments found in task.ext.args (${forbidden_args.join(', ')}). The module sets -h, -O, -p, and -z automatically."
    }

    """
    touch ${prefix}.res
    """
}
