/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BASE_QUALITY_SCORE_RECALIBRATION as BQSR_BOOTSTRAP } from '../base_quality_score_recalibration'
include { CALL_VARIANTS_GATK as CALL_VARIANTS_GATK_BOOTSTRAP } from '../call_variants_gatk'
include { FILTER_VARIANTS as FILTER_VARIANTS_BOOTSTRAP       } from '../filter_variants'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow BOOTSTRAP_VARIANT_SET {
    take:
    fasta       // channel: [ meta, fasta]
    fai         // channel: [ meta, fai]
    dict        // channel: [ meta, dict]
    intervals   // channel: [ meta, intervals, number_of_intervals]
    cram        // channel: [ meta, cram]
    crai        // channel: [ meta, crai]
    round       // channel: integer bootstrap round number

    main:
    versions = channel.empty()
    multiqc_files = channel.empty()

    // Add bootstrapping metadata to channels
    fasta.map { meta, fasta_file ->
        tuple( meta + [bootstrapping_round: round], fasta_file ) }
        .set { fasta }
    fai.map { meta, fai_file ->
        tuple( meta + [bootstrapping_round: round], fai_file ) }
        .set { fai }
    dict.map { meta, dict_file ->
        tuple( meta + [bootstrapping_round: round], dict_file ) }
        .set { dict }
    intervals.map { meta, interval_file, num_intervals ->
        tuple( meta + [bootstrapping_round: round], interval_file, num_intervals ) }
        .set { intervals }
    cram.map { meta, cram_file ->
        tuple( meta + [bootstrapping_round: round], cram_file ) }
        .set { cram }
    crai.map { meta, crai_file ->
        tuple( meta + [bootstrapping_round: round], crai_file ) }
        .set { crai }

    //
    // SUBWORKFLOW: CALL_VARIANTS_GATK_BOOTSTRAP
    //
    CALL_VARIANTS_GATK_BOOTSTRAP(
        fasta,
        fai,
        dict,
        intervals,
        cram,
        crai
    )
    versions = versions.mix(CALL_VARIANTS_GATK_BOOTSTRAP.out.versions)
    multiqc_files = multiqc_files.mix(CALL_VARIANTS_GATK_BOOTSTRAP.out.multiqc_files)

    //
    // SUBWORKFLOW: FILTER_VARIANTS
    //
    FILTER_VARIANTS_BOOTSTRAP(
        fasta,
        fai,
        dict,
        CALL_VARIANTS_GATK_BOOTSTRAP.out.vcf,
        CALL_VARIANTS_GATK_BOOTSTRAP.out.tbi
    )

    //
    // SUBWORKFLOW: BQSR_BOOTSTRAP
    //
    BQSR_BOOTSTRAP(
        fasta,
        fai,
        dict,
        intervals,
        cram,
        crai,
        FILTER_VARIANTS_BOOTSTRAP.out.vcf.collect(),
        FILTER_VARIANTS_BOOTSTRAP.out.tbi.collect()
    )
    versions = versions.mix(BQSR_BOOTSTRAP.out.versions)
    multiqc_files = multiqc_files.mix(BQSR_BOOTSTRAP.out.multiqc_files)

    // Remove bootstrapping metadata from CRAM channel
    BQSR_BOOTSTRAP.out.recalibrated_cram.map { meta, cram_file ->
        tuple( meta - meta.subMap('bootstrapping_round'), cram_file ) }
        .set { ch_cram_output }
    BQSR_BOOTSTRAP.out.recalibrated_crai.map { meta, crai_file ->
        tuple( meta - meta.subMap('bootstrapping_round'), crai_file ) }
        .set { ch_crai_output }
    FILTER_VARIANTS_BOOTSTRAP.out.vcf.map { meta, vcf_file ->
        tuple( meta - meta.subMap('bootstrapping_round'), vcf_file ) }
        .set { ch_vcf_output }
    FILTER_VARIANTS_BOOTSTRAP.out.tbi.map { meta, tbi_file ->
        tuple( meta - meta.subMap('bootstrapping_round'), tbi_file ) }
        .set { ch_tbi_output }

    emit:
    cram = ch_cram_output
    crai = ch_crai_output
    vcf  = ch_vcf_output.collect()
    tbi  = ch_tbi_output.collect()
    multiqc_files
    versions
}
