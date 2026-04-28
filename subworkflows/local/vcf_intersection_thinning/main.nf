/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BCFTOOLS_ISEC                                  } from '../../../modules/nf-core/bcftools/isec'
include { MAKE_SCAFFOLD_BED as MAKE_SCAFFOLD_EXCLUDE_BED } from '../../../modules/local/make_scaffold_bed/'
include { MAKE_SCAFFOLD_BED as MAKE_SCAFFOLD_INCLUDE_BED } from '../../../modules/local/make_scaffold_bed/'
include { VCFTOOLS as VCFTOOLS_EXCLUDE                   } from '../../../modules/nf-core/vcftools/'
include { VCFTOOLS as VCFTOOLS_THIN                      } from '../../../modules/nf-core/vcftools/'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow VCF_INTERSECTION_THINNING {
    take:
    vcf_tool1 // channel: [ meta, vcf]
    tbi_tool1 // channel: [ meta, tbi]
    vcf_tool2 // channel: [ meta, vcf]
    tbi_tool2 // channel: [ meta, tbi]
    intervals // channel: [ meta, bed, number_of_intervals]

    main:
    versions = channel.empty()

    vcf_tool1_prepared = vcf_tool1
        .map { meta, vcf ->
            tuple(meta + [id: 'called_variants'] - meta.subMap('variantcaller'), vcf)
        }
    vcf_tool2_prepared = vcf_tool2
        .map { meta, vcf ->
            tuple(meta + [id: 'called_variants'] - meta.subMap('variantcaller'), vcf)
        }
    tbi_tool1_prepared = tbi_tool1
        .map { meta, tbi ->
            tuple(meta + [id: 'called_variants'] - meta.subMap('variantcaller'), tbi)
        }
    tbi_tool2_prepared = tbi_tool2
        .map { meta, tbi ->
            tuple(meta + [id: 'called_variants'] - meta.subMap('variantcaller'), tbi)
        }

    ch_vcfs = vcf_tool1_prepared
        .mix(vcf_tool2_prepared)
        .groupTuple()
        .map { meta, vcfs -> tuple(meta, vcfs) }

    ch_tbis = tbi_tool1_prepared
        .mix(tbi_tool2_prepared)
        .groupTuple()
        .map { meta, tbis -> tuple(meta, tbis) }

    ch_isec_input = ch_vcfs
        .join(ch_tbis)
        .map { meta, vcfs, tbis -> tuple(meta, vcfs, tbis, [], [], []) }

    // Run BCFTOOLS_ISEC
    BCFTOOLS_ISEC(ch_isec_input)

    // Collect intersection output
    def has_include = params.include_scaffolds
    def has_exclude = params.exclude_scaffolds
    def need_scaffold_filter = has_include || has_exclude
    intersection = BCFTOOLS_ISEC.out.results
        .map { meta, dir ->
            def file_common = file("${dir}/0002.vcf")
            tuple(meta, file_common)
        }
        .branch { _tuple ->
            filter:      need_scaffold_filter
            passthrough: !need_scaffold_filter
        }

    def include_scaffolds = normalize_scaffold_param(params.include_scaffolds)
    def exclude_scaffolds = normalize_scaffold_param(params.exclude_scaffolds)
    include_ch = include_scaffolds ? channel.value(include_scaffolds) : channel.empty()
    exclude_ch = exclude_scaffolds ? channel.value(exclude_scaffolds) : channel.empty()

    // Make scaffold BED for inclusion
    bed_include = intervals
        .map { meta, bed_file, _number_of_intervals ->
            tuple(meta + [id: "include_scaffolds"], bed_file)
        }
    MAKE_SCAFFOLD_INCLUDE_BED(
        include_ch,
        bed_include
    )

    // Make scaffold BED for exclusion
    bed_exclude = intervals
        .map { meta, bed_file, _number_of_intervals ->
            tuple(meta + [id: "exclude_scaffolds"], bed_file)
        }
    MAKE_SCAFFOLD_EXCLUDE_BED(
        exclude_ch,
        bed_exclude
    )

    bed = MAKE_SCAFFOLD_INCLUDE_BED.out.bed
        .mix(MAKE_SCAFFOLD_EXCLUDE_BED.out.bed)
        .map { _meta, bed_file -> bed_file }.collect()

    vcftools_exclude_input = intersection.filter
        .map { meta, vcf_file ->
            tuple(meta + [id: meta.id + "_scaffolds_excluded"], vcf_file)
        }

    VCFTOOLS_EXCLUDE(
        vcftools_exclude_input,
        bed,
        [] // diff_variant_file: unused
    )
    versions = versions.mix(VCFTOOLS_EXCLUDE.out.versions)

    vcf_cleaned = VCFTOOLS_EXCLUDE.out.vcf
        .mix(intersection.passthrough)

    vcftools_thin_input = vcf_cleaned
        .map { meta, vcf_file ->
            tuple(meta + [id: meta.id + "_thinned"], vcf_file)
        }

    VCFTOOLS_THIN(
        vcftools_thin_input,
        [], // bed: unused
        []  // diff_variant_file: unused
    )
    versions = versions.mix(VCFTOOLS_THIN.out.versions)

    emit:
    intersection = VCFTOOLS_THIN.out.vcf
    versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def normalize_scaffold_param(param) {
    if (!param)
        return []

    // CASE 1: Array of strings
    if (param instanceof List)
        return param*.trim()

    // CASE 2: Single string
    if (param instanceof String) {
        // CASE 2a: a single string
        def f = file(param)
        if (f.exists()) {
            return f.readLines()
                .findAll { tuple -> tuple && !tuple.startsWith("#") }
                .collect { tuple -> tuple.trim().tokenize()[0] }
        }

        // CASE 2b: a single string
        return [param.trim()]
    }

    return []
}
