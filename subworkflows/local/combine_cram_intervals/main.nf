/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow COMBINE_CRAM_INTERVALS {
    take:
    intervals   // channel: [ meta, intervals, number_of_intervals]
    cram        // channel: [ meta, cram]

    main:

    // Combine CRAM with intervals
    cram
    .combine(intervals)
    .map { cram_meta, cram_file, interval_meta, interval_file, num_intervals ->
        // Construct new ID: sampleID_intervalName
        def new_id = "${cram_meta.id}_${interval_meta.interval_name}"

        // Merge metadata and overwrite id
        def meta = cram_meta + interval_meta + [ id: new_id ] + [ num_intervals: num_intervals ]

        tuple(meta, cram_file, interval_file)
    }.set { cram_intervals }

    emit:
    cram_intervals
}
