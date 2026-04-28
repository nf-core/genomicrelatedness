/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BCFTOOLS_CALL                      } from '../../../modules/local/bcftools/call/main'
include { BCFTOOLS_CONCAT                    } from '../../../modules/nf-core/bcftools/concat/main'
include { BCFTOOLS_MPILEUP                   } from '../../../modules/local/bcftools/mpileup/main'
include { SAMTOOLS_CONVERT                   } from '../../../modules/nf-core/samtools/convert/main'

include { COMBINE_CRAM_INTERVALS } from '../combine_cram_intervals'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow CALL_VARIANTS_BCFTOOLS {
    take:
    fasta       // channel: [ meta, fasta]
    fai         // channel: [ meta, fai]
    dict        // channel: [ meta, dict]
    intervals   // channel: [ meta, intervals, number_of_intervals]
    cram        // channel: [ meta, cram]
    crai        // channel: [ meta, crai]

    main:
    versions = channel.empty()
    multiqc_files = channel.empty()

    // Add variant caller metadata to channels
    fasta.map { meta, fasta_file ->
        tuple( meta + [variantcaller: 'bcftools'], fasta_file ) }
        .set { fasta }
    fai.map { meta, fai_file ->
        tuple( meta + [variantcaller: 'bcftools'], fai_file ) }
        .set { fai }
    dict.map { meta, dict_file ->
        tuple( meta + [variantcaller: 'bcftools'], dict_file ) }
        .set { dict }
    intervals.map { meta, interval_file, num_intervals ->
        tuple( meta + [variantcaller: 'bcftools'], interval_file, num_intervals ) }
        .set { intervals }
    cram.map { meta, cram_file ->
        tuple( meta + [variantcaller: 'bcftools'], cram_file ) }
        .set { cram }
    crai.map { meta, crai_file ->
        tuple( meta + [variantcaller: 'bcftools'], crai_file ) }
        .set { crai }

    // Build reference tuple for SAMTOOLS_CONVERT input signature
    ch_reference = fasta.join(fai).collect()

    // Convert CRAM/CAI to BAM/BAI
    ch_cram_crai_to_convert = cram.join(crai)
    SAMTOOLS_CONVERT(ch_cram_crai_to_convert, ch_reference)

    // Collect a list of all BAM files
    ch_bams = SAMTOOLS_CONVERT.out.bam
        .map { _meta, bam_file -> bam_file }
        .collect()
        .map { bam_list -> [bam_list]}

    ch_intervals = intervals
        .map { meta, interval_file, _num_intervals ->
            def new_id = meta.id + ".mpileup"
            def new_meta = meta + [id: new_id]
            tuple(new_meta, interval_file)
        }

    ch_mpileup_input = ch_intervals
        .combine(ch_bams)

    // Run Bcftools mpileup
    keep_bcftools_mpileup = false
    BCFTOOLS_MPILEUP(ch_mpileup_input, fasta, keep_bcftools_mpileup)

    ch_call_input = BCFTOOLS_MPILEUP.out.vcf
        .map { meta, vcf_file ->
            def new_meta = meta + [id: meta.interval_name]
            tuple(new_meta, vcf_file)
        }
        .join(intervals)
        .map { meta, vcf_file, interval_file, _num_intervals ->
            def new_id = meta.interval_name + ".called"
            def new_meta = meta + [id: new_id]
            tuple(new_meta, vcf_file, interval_file)
        }

    // Run Bcftools call
    BCFTOOLS_CALL(ch_call_input)

    vcf_list = BCFTOOLS_CALL.out.vcf
        .map { meta, vcf_file ->
            def new_id = "called_variants" + ".${meta.variantcaller}"
            def new_meta = meta + [id: new_id] - meta.subMap('interval_name', 'interval_idx')
            tuple(new_meta, vcf_file)
        }
        .groupTuple()

    index_list = BCFTOOLS_CALL.out.tbi
        .map { meta, vcf_file ->
            def new_id = "called_variants" + ".${meta.variantcaller}"
            def new_meta = meta + [id: new_id] - meta.subMap('interval_name', 'interval_idx')
            tuple(new_meta, vcf_file)
        }
        .groupTuple()

    ch_concat_input = vcf_list
        .join(index_list)

    // Run Bcftools concat
    BCFTOOLS_CONCAT(ch_concat_input)

    emit:
    vcf = BCFTOOLS_CONCAT.out.vcf
    tbi = BCFTOOLS_CONCAT.out.tbi
    multiqc_files
    versions
}
