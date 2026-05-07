/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { GATK4_BASERECALIBRATOR  } from '../../../modules/nf-core/gatk4/baserecalibrator'
include { GATK4_GATHERBQSRREPORTS } from '../../../modules/nf-core/gatk4/gatherbqsrreports/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow CRAM_BASERECALIBRATOR {
    take:
    fasta   // channel: [ meta, fasta]
    fai     // channel: [ meta, fai]
    dict    // channel: [ meta, dict]
    cram    // channel: [ meta, cram, crai, intervals ]
    vcf     // channel: [ meta, vcf]
    tbi     // channel: [ meta, tbi]

    main:

    // Run BaseRecalibrator
    GATK4_BASERECALIBRATOR(
        cram,
        fasta,
        fai,
        dict,
        vcf.map { _meta, files -> [[id:'known_sites'], files] },
        tbi.map { _meta, files -> [[id:'known_sites'], files] }
    )

    // Figuring out if there is one or more table(s) from the same sample
    ch_table_to_merge = GATK4_BASERECALIBRATOR.out.table
        .map{ meta, table ->
            // Use sample name and bootstrapping stage as key, ensure num_intervals is available
            def sample_name = meta.sample ?: meta.id.split('_')[0]
            def new_id = "${sample_name}" + (meta.bootstrapping_round ? "_${meta.bootstrapping_round}" : "") + (meta.pass == 2 ? "_second_pass" : "")

            // Remove interval_name from meta in order to group by sample only
            def new_meta = meta + [id: new_id] - meta.subMap('interval_name', 'interval_idx', 'reference_fasta')
            tuple(new_meta, table)
        }
        .groupTuple()
        .branch{ tuple ->
            // Use meta.num_intervals to asses number of intervals
            single:   tuple[0].num_intervals <= 1
            multiple: tuple[0].num_intervals > 1
        }

    // Only when using intervals
    GATK4_GATHERBQSRREPORTS(ch_table_to_merge.multiple)

    // Mix intervals and no_intervals channels together
    table_bqsr = GATK4_GATHERBQSRREPORTS.out.table
        .mix(ch_table_to_merge.single
            .map{ meta, table ->
                [ meta, table[0] ]
            }
        )
        // Remove no longer necessary field: num_intervals
        .map{ meta, table -> [ meta - meta.subMap('num_intervals'), table ] }

    emit:
    table_bqsr
}
