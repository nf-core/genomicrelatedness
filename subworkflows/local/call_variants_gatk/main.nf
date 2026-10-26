/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BCFTOOLS_SORT          } from '../../../modules/nf-core/bcftools/sort/main'
include { BCFTOOLS_STATS         } from '../../../modules/nf-core/bcftools/stats'
include { GATK4_GENOMICSDBIMPORT } from '../../../modules/nf-core/gatk4/genomicsdbimport'
include { GATK4_GENOTYPEGVCFS    } from '../../../modules/nf-core/gatk4/genotypegvcfs'
include { GATK4_HAPLOTYPECALLER  } from '../../../modules/nf-core/gatk4/haplotypecaller'
include { GATK4_MERGEVCFS        } from '../../../modules/nf-core/gatk4/mergevcfs'

include { COMBINE_CRAM_CRAI_INTERVALS } from '../combine_cram_crai_intervals'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow CALL_VARIANTS_GATK {
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
        tuple( meta + [variantcaller: 'gatk'], fasta_file ) }
        .set { fasta }
    fai.map { meta, fai_file ->
        tuple( meta + [variantcaller: 'gatk'], fai_file ) }
        .set { fai }
    dict.map { meta, dict_file ->
        tuple( meta + [variantcaller: 'gatk'], dict_file ) }
        .set { dict }
    intervals.map { meta, interval_file, num_intervals ->
        tuple( meta + [variantcaller: 'gatk'], interval_file, num_intervals ) }
        .set { intervals }
    cram.map { meta, cram_file ->
        tuple( meta + [variantcaller: 'gatk'], cram_file ) }
        .set { cram }
    crai.map { meta, crai_file ->
        tuple( meta + [variantcaller: 'gatk'], crai_file ) }
        .set { crai }

    // Combine CRAM with CRAI and intervals
    COMBINE_CRAM_CRAI_INTERVALS(intervals, cram, crai)

    // Prepare HaplotypeCaller input
    ch_haplotypecaller_input = COMBINE_CRAM_CRAI_INTERVALS.out.cram_crai_intervals
        .map { meta, cram_file, crai_file, interval_file ->
            def new_id = meta.id + (meta.bootstrapping_round ? "_${meta.bootstrapping_round}" : "") + "_${meta.variantcaller}"
            def new_meta = meta + [ id: new_id ]
            tuple(new_meta, cram_file, crai_file, interval_file, [])
        }

    // Run GATK HaplotypeCaller
    GATK4_HAPLOTYPECALLER(ch_haplotypecaller_input, fasta, fai, dict, [[id: 'no_dbsnp'], []], [[id: 'no_dbsnp_tbi'], []])

    // Prepare for GenomicsDBImport
    ch_gvcfs = GATK4_HAPLOTYPECALLER.out.vcf
        .join(GATK4_HAPLOTYPECALLER.out.tbi)
        .map { meta, vcf, tbi ->
            def key = meta + [ id: meta.interval_name ] - meta.subMap('sample', 'single_end', 'num_intervals')
            tuple(key, vcf, tbi)
        }

    // Prepare GenomicsDBImport input by grouping GVCFs by interval_name
    ch_gdb_input = ch_gvcfs
        .groupTuple()
        .join(intervals)
        .map { meta, vcfs, tbis, bed, _num_intervals ->
            def new_meta = meta + [id: meta.interval_name + (meta.bootstrapping_round ? "_${meta.bootstrapping_round}" : "") + "_joint"]
            tuple(
                new_meta,
                vcfs,
                tbis,
                bed,
                [],
                file('.')
            )
        }

    // Run GATK GenomicsDBImport
    GATK4_GENOMICSDBIMPORT(ch_gdb_input, false, false, false)

    // Run GATK GenotypeGVCFs
    ch_gtp_input = GATK4_GENOMICSDBIMPORT.out.genomicsdb.map { meta, genomicsdb -> tuple(meta, genomicsdb, [], [], []) }
    GATK4_GENOTYPEGVCFS(ch_gtp_input, fasta, fai, dict, [[id: 'no_dbsnp'], []], [[id: 'no_dbsnp_tbi'], []])

    // Sort each interval VCF before merging
    ch_vcfs = GATK4_GENOTYPEGVCFS.out.vcf
        .map { meta, vcf ->
            def new_meta = meta + [ id: "${meta.id}.sorted" ]
            tuple(new_meta, vcf)
        }

    BCFTOOLS_SORT(ch_vcfs)

    ch_merge_vcfs = BCFTOOLS_SORT.out.vcf
            .toSortedList { a, b -> a[0].interval_idx <=> b[0].interval_idx }
            .flatMap { list ->
                if (!list || list.isEmpty())
                    return []

                def metas = list.collect { tuple -> tuple[0] }
                def vcfs  = list.collect { tuple -> tuple[1] }
                def base  = metas[0]

                def new_meta = base + [
                    id: "called_variants" +
                        (base.bootstrapping_round ? "_${base.bootstrapping_round}" : "") +
                        ".${base.variantcaller}"
                ] - base.subMap('interval_name', 'interval_idx')

                return [ tuple(new_meta, vcfs) ]
            }

    // Merge all intervals into one VCF
    GATK4_MERGEVCFS(ch_merge_vcfs, dict)

    // Run BCFtools stats on merged VCF
    ch_merged_vcf_tbi = GATK4_MERGEVCFS.out.vcf
        .join(GATK4_MERGEVCFS.out.tbi)
        .map { meta, vcf, tbi -> tuple(meta, vcf, tbi) }

    BCFTOOLS_STATS(
        ch_merged_vcf_tbi,
        [[id: 'no_regions'], []],
        [[id: 'no_targets'], []],
        [[id: 'no_samples'], []],
        [[id: 'no_exons'], []],
        fasta
    )
    multiqc_files = multiqc_files.mix(BCFTOOLS_STATS.out.stats.map { tuple -> tuple[1] })

    emit:
    vcf = GATK4_MERGEVCFS.out.vcf
    tbi = GATK4_MERGEVCFS.out.tbi
    multiqc_files
    versions
}
