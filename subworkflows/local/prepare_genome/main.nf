/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BWAMEM2_INDEX                  } from '../../../modules/nf-core/bwamem2/index'
include { GATK4_CREATESEQUENCEDICTIONARY } from '../../../modules/nf-core/gatk4/createsequencedictionary'
include { SAMTOOLS_FAIDX                 } from '../../../modules/nf-core/samtools/faidx'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow PREPARE_GENOME {
    take:
    fasta   // channel: [ meta, fasta]

    main:

    // Build the BWA index from the provided FASTA
    BWAMEM2_INDEX(fasta)

    // Build the sequence dictionary
    GATK4_CREATESEQUENCEDICTIONARY(fasta)

    // Build reference tuple for SAMTOOLS_FAIDX input signature
    ch_reference = fasta
        .map { meta, fasta_file -> tuple(meta, fasta_file, []) }

    // Build the FASTA index (fai)
    SAMTOOLS_FAIDX(ch_reference, false)

    emit:
    bwamem2_index   = BWAMEM2_INDEX.out.index.collect()
    dict            = GATK4_CREATESEQUENCEDICTIONARY.out.dict.collect()
    fasta_fai       = SAMTOOLS_FAIDX.out.fai.collect()
}
