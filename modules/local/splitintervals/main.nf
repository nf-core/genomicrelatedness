process SPLIT_INTERVALS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(intervals)

    output:
    tuple val(meta), path("*.bed"), emit: bed

    script:
    """
    split_intervals.py \
        --bed ${intervals} \
        --target-number-files ${params.target_number_of_interval_files} \
        --max-intervals ${params.max_number_of_intervals_per_file} \
        --out-prefix interval
    """
}
