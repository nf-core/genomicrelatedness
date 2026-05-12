/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BWAMEM2_MEM                    } from '../../../modules/nf-core/bwamem2/mem'
include { FASTP                          } from '../../../modules/nf-core/fastp'
include { GATK4_ADDORREPLACEREADGROUPS   } from '../../../modules/nf-core/gatk4/addorreplacereadgroups'
include { GATK4_MARKDUPLICATES           } from '../../../modules/nf-core/gatk4/markduplicates'
include { MOSDEPTH                       } from '../../../modules/nf-core/mosdepth'
include { NORMALIZE_BAM_NAMES            } from '../../../modules/local/normalize_bam_names'
include { PRESEQ_CCURVE                  } from '../../../modules/nf-core/preseq/ccurve'
include { PRESEQ_LCEXTRAP                } from '../../../modules/nf-core/preseq/lcextrap'
include { SAMTOOLS_INDEX                 } from '../../../modules/nf-core/samtools/index'
include { SAMTOOLS_MERGE                 } from '../../../modules/nf-core/samtools/merge/main'
include { SAMTOOLS_STATS                 } from '../../../modules/nf-core/samtools/stats'
include { SPRING_DECOMPRESS              } from '../../../modules/nf-core/spring/decompress'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow PREPROCESS {
    take:
    samplesheet // channel: [ meta, list(fastq) ]
    fasta       // channel: [ meta, fasta]
    fai         // channel: [ meta, fai]
    bwamem2     // channel: [ meta, bwamem2 ]

    main:
    versions = channel.empty()
    multiqc_files = channel.empty()

    // Split by file type (spring vs fastq)
    samplesheet.branch { row ->
        spring: row[1].every { file -> file.getName().endsWith('.spring') }
        fastq : row[1].every { file -> file.getName().endsWith('.fastq') || file.getName().endsWith('.fastq.gz') || file.getName().endsWith('.fq.gz') }
        bam   : row[1].every { file -> file.getName().endsWith('.bam') }
        cram  : row[1].every { file -> file.getName().endsWith('.cram') }
    }.set{ ch_input_branches }

    // Decompress SPRING to FASTQ pairs
    SPRING_DECOMPRESS(ch_input_branches.spring, false)

    // Merge with normal FASTQs into one unified channel
    merged_fastqs = ch_input_branches.fastq.mix(SPRING_DECOMPRESS.out.fastq)

    // Trim and QC with FASTP
    ch_fastp_input = merged_fastqs
        .map { meta, reads ->
            def new_id = meta.RGSM + "_RGID${meta.RGID}"
            def new_meta = meta + [id: new_id]
            tuple(new_meta, reads, [])
        }
    FASTP(ch_fastp_input, false, false, false)
    multiqc_files = multiqc_files.mix(FASTP.out.json.collect{ _meta, json -> json })
    multiqc_files = multiqc_files.mix(FASTP.out.html.collect{ _meta, html -> html })

    // Map to reference
    BWAMEM2_MEM(FASTP.out.reads, bwamem2, fasta, true)
    bam = BWAMEM2_MEM.out.bam.mix(ch_input_branches.bam)

    // Add read groups
    GATK4_ADDORREPLACEREADGROUPS(bam, fasta, fai)

    grouped_bams = GATK4_ADDORREPLACEREADGROUPS.out.bam
                        .map { meta, bam_file -> tuple(meta.sample, meta, bam_file) }
                        .groupTuple()

    single_bams = grouped_bams
        .filter { _sample, _metas, bams -> bams.size() == 1 }
        .map { sample, metas, bams ->
            def m = metas[0]

            def meta = [
                id        : sample,
                sample    : sample,
                single_end: m.single_end
            ]

            tuple(meta, bams[0])
        }

    multi_bams = grouped_bams
        .filter { _sample, _metas, bams -> bams.size() > 1 }
        .map { sample, metas, bams ->
            def m = metas[0]

            def meta = [
                id        : sample,
                sample    : sample,
                single_end: m.single_end
            ]

            tuple(meta, bams, [])
        }

    // Build reference tuple for SAMTOOLS_MERGE input signature
    ch_merge_reference = fasta.join(fai)
        .map { meta, fasta_file, fai_file -> tuple(meta, fasta_file, fai_file, []) }
        .collect()

    // Merge BAMs per-sample
    SAMTOOLS_MERGE(
        multi_bams,
        ch_merge_reference
    )

    NORMALIZE_BAM_NAMES(single_bams)

    merged_bams = NORMALIZE_BAM_NAMES.out.bam.mix(SAMTOOLS_MERGE.out.bam)

    // Mark duplicates
    GATK4_MARKDUPLICATES(merged_bams, fasta.map { tuple -> tuple[1] }, fai.map{ tuple -> tuple[1] })
    multiqc_files = multiqc_files.mix(GATK4_MARKDUPLICATES.out.metrics.map { tuple -> tuple[1] })
    cram = GATK4_MARKDUPLICATES.out.cram
        .mix(ch_input_branches.cram)

    // Compute index
    SAMTOOLS_INDEX(cram)
    crai = SAMTOOLS_INDEX.out.index

    // Preseq analyses
    PRESEQ_CCURVE(cram)
    versions = versions.mix(PRESEQ_CCURVE.out.versions)
    multiqc_files = multiqc_files.mix(PRESEQ_CCURVE.out.c_curve.map { _meta, file -> file }).mix(PRESEQ_CCURVE.out.log.map{ _meta, file -> file })

    PRESEQ_LCEXTRAP(cram)
    multiqc_files = multiqc_files.mix(PRESEQ_LCEXTRAP.out.lc_extrap.map { _meta, file -> file }).mix(PRESEQ_LCEXTRAP.out.log.map{ _meta, file -> file })

    // Samtools stats on final CRAMs
    samstats_input = cram
        .join(crai)
        .map { meta, cram_file, crai_file ->
            tuple(meta, cram_file, crai_file)
        }
    samtools_reference = fasta
        .join(fai)
        .map { meta, fasta_file, fai_file ->
            tuple(meta, fasta_file, fai_file)
        }
        .collect()
    SAMTOOLS_STATS(samstats_input, samtools_reference)
    multiqc_files = multiqc_files.mix(SAMTOOLS_STATS.out.stats.map { tuple -> tuple[1] })

    // Coverage calculation with mosdepth
    mosdepth_input = cram.join(crai).map { meta, cram_file, crai_file -> tuple(meta, cram_file, crai_file, []) }
    MOSDEPTH(mosdepth_input, fasta, [])
    multiqc_files = multiqc_files.mix(MOSDEPTH.out.global_txt.map { _meta, file -> file }).mix(MOSDEPTH.out.summary_txt.map { _meta, file -> file })

    emit:
    cram
    crai
    multiqc_files
    versions
}
