/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow COMBINE_CRAM_CRAI_INTERVALS {
    take:
    intervals   // channel: [ meta, intervals, number_of_intervals]
    cram        // channel: [ meta, cram]
    crai        // channel: [ meta, crai]

    main:
    // Combine CRAM with intervals
    cram.join(crai)
    .combine(intervals)
    .map { cram_meta, cram_file, crai_file, interval_meta, interval_file, num_intervals ->
        // Construct new ID: sampleID_intervalName
        def new_id = "${cram_meta.id}_${interval_meta.interval_name}"

        // Merge metadata and overwrite id
        def meta = cram_meta + interval_meta + [ id: new_id ] + [ num_intervals: num_intervals ]

        tuple(meta, cram_file, crai_file, interval_file)
    }.set { cram_crai_intervals }

    emit:
    cram_crai_intervals
}
