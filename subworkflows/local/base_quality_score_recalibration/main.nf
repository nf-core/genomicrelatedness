/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { GATK4_ANALYZECOVARIATES } from '../../../modules/nf-core/gatk4/analyzecovariates'
include { GATK4_APPLYBQSR         } from '../../../modules/nf-core/gatk4/applybqsr'
include { SAMTOOLS_INDEX          } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_MERGE          } from '../../../modules/nf-core/samtools/merge/main'

include { COMBINE_CRAM_CRAI_INTERVALS                                            } from '../combine_cram_crai_intervals'
include { COMBINE_CRAM_CRAI_INTERVALS as COMBINE_CRAM_CRAI_INTERVALS_SECOND_PASS } from '../combine_cram_crai_intervals'
include { CRAM_BASERECALIBRATOR                                                  } from '../cram_baserecalibrator'
include { CRAM_BASERECALIBRATOR as CRAM_BASERECALIBRATOR_SECOND_PASS             } from '../cram_baserecalibrator'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow BASE_QUALITY_SCORE_RECALIBRATION {
    take:
    fasta       // channel: [ meta, fasta]
    fai         // channel: [ meta, fai]
    dict        // channel: [ meta, dict]
    intervals   // channel: [ meta, intervals, number_of_intervals]
    cram        // channel: [ meta, cram]
    crai        // channel: [ meta, crai]
    vcf         // channel: [ meta, vcf]
    tbi         // channel: [ meta, tbi]

    main:
    versions = channel.empty()
    multiqc_files = channel.empty()

    // Combine CRAM with intervals
    COMBINE_CRAM_CRAI_INTERVALS(intervals, cram, crai)
    combined_cram_crai_intervals = COMBINE_CRAM_CRAI_INTERVALS.out.cram_crai_intervals
        .map { meta, cram_file, crai_file, interval_file ->
            def new_id = meta.id + (meta.bootstrapping_round ? "_${meta.bootstrapping_round}" : "")
            def new_meta = meta + [ id: new_id ]
            tuple(new_meta, cram_file, crai_file, interval_file)
        }

    // Run BaseRecalibrator
    CRAM_BASERECALIBRATOR(fasta, fai, dict, combined_cram_crai_intervals, vcf, tbi)

    // Combine CRAM with BQSR table
    ch_cram_with_table = combined_cram_crai_intervals
        .combine(CRAM_BASERECALIBRATOR.out.table_bqsr)
        .filter { meta_cc, _cram_file, _crai_file, _interval_file, meta_tab, _table ->
            // only keep pairs where sample IDs match
            meta_cc.sample == meta_tab.sample
        }
        .map { meta_cram, cram_file, crai_file, interval_file, _meta_table, table ->
            tuple(meta_cram, cram_file, crai_file, table, interval_file)
        }

    // Combine fasta, fai, dict
    ch_reference = fasta
        .combine(fai)
        .combine(dict)
        .map {
            meta_fasta, fasta_file,
            _meta_fai,  fai_file,
            _meta_dict, dict_file ->

            tuple(meta_fasta, fasta_file, fai_file, dict_file)
        }
        .first() // Only one fasta, fai, dict, respectively, so we can take the first

    // Run ApplyBQSR
    GATK4_APPLYBQSR(
        ch_cram_with_table,
        ch_reference,
        "cram"
    )

    // Merge recalibrated CRAMs if needed
    ch_cram_branch = GATK4_APPLYBQSR.out.cram
        .map{ meta, cram_file ->
            def new_id = (meta.sample ?: meta.id.split('_')[0]) + (meta.bootstrapping_round ? "_${meta.bootstrapping_round}" : "") + "_recalibrated"
            def new_meta = meta + [ id: new_id ] - meta.subMap('interval_name', 'interval_idx')
            tuple(new_meta, cram_file)
        }
        .groupTuple()
        .branch { tuple ->
            single:   tuple[0].num_intervals <= 1
            multiple: tuple[0].num_intervals > 1
        }

    // Build reference tuple for SAMTOOLS_MERGE input signature
    ch_merge_reference = fasta.join(fai)
        .map { meta, fasta_file, fai_file -> tuple(meta, fasta_file, fai_file, []) }
        .collect()

    // Build input for samtools/merge
    merge_input = ch_cram_branch.multiple
        .map { meta, input_files -> tuple(meta, input_files, [])}

    // Merge CRAMs if multiple intervals
    SAMTOOLS_MERGE(
        merge_input,
        ch_merge_reference
    )

    // Mix intervals and no_intervals channels together
    ch_recalibrated_cram = SAMTOOLS_MERGE.out.cram
        .mix(ch_cram_branch.single)
        .map{ meta, cram_file ->
            tuple(meta - meta.subMap('interval_name', 'num_intervals'), cram_file)
        }

    // Index CRAM
    SAMTOOLS_INDEX(ch_recalibrated_cram)

    // Remove 'recalibrated' from ID
    ch_recalibrated_cram = ch_recalibrated_cram
        .map { meta, cram_file ->
            def new_id = (meta.sample ?: meta.id.split('_')[0]) + (meta.bootstrapping_round ? "_${meta.bootstrapping_round}" : "")
            tuple(meta + [id: new_id], cram_file)
        }

    // Remove 'recalibrated' from ID
    ch_recalibrated_crai = SAMTOOLS_INDEX.out.index
        .map { meta, crai_file ->
            def new_id = (meta.sample ?: meta.id.split('_')[0]) + (meta.bootstrapping_round ? "_${meta.bootstrapping_round}" : "")
            tuple(meta + [id: new_id], crai_file)
        }

    // Combine CRAM with intervals for (second pass, for quality control)
    COMBINE_CRAM_CRAI_INTERVALS_SECOND_PASS(intervals, ch_recalibrated_cram, ch_recalibrated_crai)

    combined_cram_crai_intervals_second_pass = COMBINE_CRAM_CRAI_INTERVALS_SECOND_PASS.out.cram_crai_intervals
        .map { meta, cram_file, crai_file, interval_file ->
            def new_meta = meta + [id: meta.id + "_second_pass"] + [ pass: 2 ]
            tuple(new_meta, cram_file, crai_file, interval_file)
        }

    // Run BaseRecalibrator (second pass, for quality control)
    CRAM_BASERECALIBRATOR_SECOND_PASS(
        fasta,
        fai,
        dict,
        combined_cram_crai_intervals_second_pass,
        vcf,
        tbi
    )

    ch_bqsr_first = CRAM_BASERECALIBRATOR.out.table_bqsr
        .map { meta, table ->
            def new_id = (meta.sample ?: meta.id.split('_')[0]) + (meta.bootstrapping_round ? "_${meta.bootstrapping_round}" : "")
            def new_meta = meta - meta.subMap('sample', 'single_end') + [id: new_id] + [ pass: 2 ]
            tuple(new_meta, table)
        }
    ch_bqsr_second = CRAM_BASERECALIBRATOR_SECOND_PASS.out.table_bqsr
        .map { meta, table ->
            def new_id = (meta.sample ?: meta.id.split('_')[0]) + (meta.bootstrapping_round ? "_${meta.bootstrapping_round}" : "")
            def new_meta = meta - meta.subMap('sample', 'single_end') + [id: new_id]
            tuple(new_meta, table)
        }
    ch_bqsr_tables = ch_bqsr_first
        .join(ch_bqsr_second)
        .map { meta, before_table, after_table ->
            tuple(meta, before_table, after_table, [])
        }

    // Run AnalyzeCovariates
    GATK4_ANALYZECOVARIATES(ch_bqsr_tables)

    emit:
    recalibrated_cram = ch_recalibrated_cram
    recalibrated_crai = ch_recalibrated_crai
    multiqc_files
    versions
}
