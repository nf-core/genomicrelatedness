/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { GAWK as BUILD_INTERVALS } from '../../../modules/nf-core/gawk'
include { SPLIT_INTERVALS         } from '../../../modules/local/splitintervals'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow PREPARE_INTERVALS {
    take:
    fai // channel: [ meta, fai]

    main:

    // Build intervals from FASTA index
    BUILD_INTERVALS(fai, [], false)

    intervals_combined = BUILD_INTERVALS.out.output
    .map { meta, intervals ->
        def num_intervals = intervals.readLines().size()
        tuple(meta + [ reference_fasta: meta.id ], intervals, num_intervals)
    }

    intervals_combined_branched = BUILD_INTERVALS.out.output
        .branch { _tuple ->
            load_from_file:  params.intervals
            do_split:       !params.intervals
        }

    if (params.intervals) {
        // Load intervals from directory
        split_with_meta = fai
            .map { meta, _fai ->
                def beds = file("${params.intervals}/interval_*.bed").sort()
                tuple(meta, beds)
            }
    } else {
        // Split intervals into separate files
        SPLIT_INTERVALS(intervals_combined_branched.do_split)
        split_with_meta = SPLIT_INTERVALS.out.bed
    }

    intervals_split = split_with_meta
        .flatMap { meta, beds ->
            def list = beds instanceof List ? beds : [beds]
            def count = list.size()
            list.collect { bed ->
                // Extract numeric index from filename: intervals_001.bed -> 001
                def idx_str = (bed.baseName =~ /interval_(\d+)/)[0][1]
                def idx_int = idx_str as int
                def interval_name = "I${idx_str}"

                def new_meta = meta + [
                    id              : interval_name,
                    interval_name   : interval_name,
                    interval_idx    : idx_int,
                    reference_fasta : meta.id
                ]

                tuple(new_meta, bed, count)
            }
        }

    emit:
    intervals_combined  // [[id:'reference_fasta', reference_fasta: 'reference_fasta']], interval.bed, number_of_intervals]
    intervals_split     // [[id:'interval_name', interval_name:'I001', interval_idx: 001,reference_fasta: 'reference_fasta'], interval.bed, number_of_intervals]
}
