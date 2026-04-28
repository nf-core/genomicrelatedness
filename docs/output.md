# nf-core/genomicrelatedness: Output

## Introduction

This document describes the output produced by the nfcore/genomicrelatedness pipeline. The output primarily consists of the following main components: output files (e.g. txt, CRAM, BAM or VCF files), and summary statistics of the whole run presented in a [`MultiQC`](https://multiqc.info) report. Intermediate files and module-specific statistics files are also retained.
The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory. The results directory of the pipeline needs to be specified with the `--outdir` flag when running the pipeline. During the run, intermediate files will be written to the work directory, which can be specified with the `-work-dir`parameter.

```text
{outdir}
├── bootstrapping
│   ├── round_1
│   │   ├── bqsr
│   │   |   ├── cram
|   |   |   |   └── merged
│   │   |   └── qc
|   |   ├── stats
│   │   └── variants
│   │       ├── db
│   │       ├── filtered
│   │       └── merged
|   |
│   ├── round_2
│   │   ├── …
|   ⋮    ⋮
|   ⋮
│   └── round_3
│       ├── …
|       ⋮
|
├── intervals
|
├── multiqc
|
├── pipeline_info
|
├── preprocessing
│   ├── alignment
|   |    ├── bam
|   |    |   └── bwamem2
│   |    └── cram
|   |
│   ├── coverage
│   ├── fastp
│   ├── preseq
│   └── stats
|
├── variant_calling
│   ├── bcftools
|   |    ├── merged
│   |    └── bam
|   |
│   └── gatk
|        ├── I<interval>_joint
|        ⋮
|        ├── merged
│        └── stats
|
└── relatedness_estimation
    ├── exclude
    |
    ├── intersection
    |
    ├── thinned
    |
    └── angsd_ngsrelate
work/
.nextflow.log
```

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [Preprocessing Section](#preprocessing-section)
  - [Prepare Reference Genome](#prepare-reference-genome)
  - [Prepare Intervals](#prepare-intervals)
  - [Prepare Input Files](#prepare-input-files)
  - [Map to Reference](#map-to-reference)
  - [Mark Duplicates](#mark-duplicates)
  - [Preprocessing Statistics](#preprocessing-statistics)
- [Bootstrapping Section](#bootstrapping-section)
  - [Call Variants](#call-variants)
  - [Hard Filter Variants](#hard-filter-variants)
  - [Base Quality Score Recalibration](#base-quality-score-recalibration)
- [Variant Calling Section](#variant-calling-section)
- [Relatedness Estimation Section](#relatedness-estimatio-section)
- [MultiQC](#multiqc)
- [Pipeline Information](#pipeline-information)

## Preprocessing section

The first section prepares the input files for downstream analyses.

### Prepare Reference Genome

The subworkflow `prepare_genome` indexes the reference genome with [`BWA-mem2 -index`](https://bio-bwa.sourceforge.net/bwa.shtml), creates a sequencing dictionary with the [`GATK4 -createsequencedictionary`](https://gatk.broadinstitute.org/hc/en-us/articles/360036729911-CreateSequenceDictionary-Picard), and generates a fasta index file with [`samtools -faidx`](https://www.htslib.org/doc/samtools-faidx.html), providing the basis for mapping of the sequencing reads.

### Prepare Intervals

The subworkflow `build_intervals` builds intervals with [`gawk -build_intervals`] and splits them with [`split_intervals`], providing a computational efficiente basis for the Bootstrapping and Variant Calling Sections.

### Prepare Input Files

The raw sequencing reads are prepared for downstream analysis. The input files are the sample table, raw sequencing reads FASTQ files, and the reference genome in FASTA format. The sequencing reads undergo quality control, trimming, filtering, and merging of paired reads with [`fastp`](https://github.com/OpenGene/fastp).

### Map to Reference

The processed sequencing reads are mapped to the reference genome with the [`BWA-mem2 -mem`](https://bio-bwa.sourceforge.net/bwa.shtml). Sorts the reads with [`samtools -sort`](https://www.htslib.org/doc/samtools-sort.html) and add read groups with [`GATK4 -addorreplacereadgroups`](https://janis.readthedocs.io/en/latest/tools/bioinformatics/gatk4/gatk4createsequencedictionary.html).

### Mark Duplicates

Reads stemming from the same sample are merged with [`samtools -merge`](https://www.htslib.org/doc/samtools-merge.html), duplicates removed with [`GATK4 -markduplicates`](https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard), and indexed with [`samtools -index`](https://www.htslib.org/doc/samtools-index.html), providing the cram file for further analysis.

### Preprocessing statistics

Library complexity is assessed with [`preseq -c_curve`](https://preseq.readthedocs.io/en/latest/) and [`preseq -lcextrap`](https://preseq.readthedocs.io/en/latest/) (errorStrategy = ‘ignore’; maxRetries = 1). [`mosdepth`](https://github.com/brentp/mosdepth) and [`samtools -stats`](https://www.htslib.org/doc/samtools-stats.html) are used to summarize coverage and mapping statistics.

<details markdown="1">
<summary>Output files</summary>

- `preprocessing/`
  - `alignment/`: directory containing the bam files for each sample_RGID (in the subdirectory `bam/bwamem2/`) and the deduplicated cram and crai files for each individual with corresponding deduplication metrics (in the subdirectory `cram/`).
  - `coverage/`: directory containing text files with coverage statistics output from mosdepth (global and summarized for each contig) for each individual.
  - `fastp/`: directory containing for each individual log, html and json files of the fastp preprocessing statistics, as well as trimmed and filterd fastq files.
  - `preseq/`: directory containing the estimates on library complexity for each individual provided by c-Curve and lc_extrap
  - `stats/`: directory containing the sequencing read statistics for each individuals.
  </details>

## Bootstrapping Section

The second section, **bootstrapping**, corrects systematic errors introduced during sequencing by recalibrating base quality scores. If an external VCF file with a reference variant set is provided this section directly starts with the BQSR subworkflow.

### Call Variants

In the absence of a known variant set, the workflow first performs an internal bootstrapping procedure to create a temporary high-confidence variant resource. The cram files from the Preprocessing section, the bed file with intervals from the `prepare_intervals` subworkflow, and the fasta reference genome are combined and used for an initial round of variant calling with the subworkflow `GATK_variant_calling`. [`GATK4 -HaplotypeCaller`](https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller) produces gVCFs for each sample and scaffold (-ERC GVCF), which are then combined with [`GATK4 -GenomicsDBImport`](https://gatk.broadinstitute.org/hc/en-us/articles/360036883491-GenomicsDBImport) and jointly genotyped with [`GATK4 -GenotypeGVCFs`](https://gatk.broadinstitute.org/hc/en-us/articles/360037057852-GenotypeGVCFs) and [`GATK4 -mergevcfs`](https://gatk.broadinstitute.org/hc/en-us/articles/360036713331-MergeVcfs-Picard).

### Hard Filter Variants

The subworkflow `variant_filtering` then uses hard filter criteria to create a reference variant set in vcf format with et (.vcf) with [`GATK4 -VariantFiltration`](https://gatk.broadinstitute.org/hc/en-us/articles/360037434691-VariantFiltration) (Qual >=100, QD < 2.0; MQ < 35.0; FS >60; HaplotypeScore > 13.0; MQRankSum < -12.5; ReadPosRankSum < -8.0) and [`GATK4 -SelectVariants`](https://gatk.broadinstitute.org/hc/en-us/articles/360037055952-SelectVariants).

### Base Quality Score Recalibration

In the subworkflow `BQSR` [`GATK4 -BaseRecalibrator`](https://gatk.broadinstitute.org/hc/en-us/articles/360036898312-BaseRecalibrator) takes as input the output from the preprocessing section and the vcf file with the reference variant set (either the one produced previously or an existing one) to compute recalibration tables. These are compiled with [`GATK4 -GatherBQSRReport`](https://gatk.broadinstitute.org/hc/en-us/articles/360037433771-GatherBQSRReports) and applied with [`GATK4 -ApplyBQSR`](https://gatk.broadinstitute.org/hc/en-us/articles/360037268511-ApplyBQSR) followed by [`samtools -merge`](https://www.htslib.org/doc/samtools-merge.html) and [`samtools -index`](https://www.htslib.org/doc/samtools-index.html) to produce recalibrated CRAMs. The BQSR subworkflow is designed to be iterated, using the recalibrated CRAM files instead of the ones from the preprocessing section. For each round, ith [`GATK4 -AnalyzeCovariates`](https://gatk.broadinstitute.org/hc/en-us/articles/360037066912-AnalyzeCovariates) generates diagnostic plots to evaluate the effectiveness of recalibration, and[`bcftools -stats`](https://samtools.github.io/bcftools/bcftools.html#stats) produces variant quality summaries.

<details markdown="1">
<summary>Output files</summary>

- `bootstrapping/round_X`
  - `bqsr/`: directory containing the cram files for each individual and interval (in the subdirectory `cram/`) as well as the recalibrated cram with corresponding index .crai files merged for each individual (in the subdirectory `cram/merged/`), and the recalibration diagnostics as pdf and csv files for each sample (in the subdirectory `qc/`).
  - `stats/`: directory containing the text file with variant statistics.
  - `variants/`: directory containing the called variants for each individuals as vcf and tbi files for each interval, the merged vcf and tbi file (in the subdirectory `merged/`) and the hard filtered variant set (in the subdirectory `filtered/`).
  </details>

## Variant Calling Section

In the third section, **variant calling**, the recalibrated CRAMs are processed to generate genotype likelihoods and multi-sample VCFs. Two variant calling approaches, bcftools and GATK4, are chosen to mitigate caller-specific biases. The subworkflow _call_variants_gatk_ is executed like in the previous section with [`GATK4 -HaplotypeCaller`](https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller) , [`GATK4 -GenomicsDBImport`](https://gatk.broadinstitute.org/hc/en-us/articles/360036883491-GenomicsDBImport), [`GATK4 -GenotypeGVCFs`](https://gatk.broadinstitute.org/hc/en-us/articles/360037057852-GenotypeGVCFs), [`GATK4 -mergevcfs`](https://gatk.broadinstitute.org/hc/en-us/articles/360036713331-MergeVcfs-Picard), the subworkflow _call_variants_bcftools_ converts the recalibrated CRAM with samtools `convert`, employs bcftools [`bcftools -mpileup`](https://samtools.github.io/bcftools/bcftools.html#mpileup) (--output-type z -d 100) and [`bcftools -call`](https://samtools.github.io/bcftools/bcftools.html#call) (--output-type z -m -v –write-index=tbi) per scaffold, and concatenates results into a full cohort VCF with [`bcftools -concat`](https://samtools.github.io/bcftools/bcftools.html#concat) (--output-type z –write-index=tbi). The callsets from both subworkflows are accompanied by variant-level quality summaries from [`bcftools stats`](https://samtools.github.io/bcftools/bcftools.html#stats). This stage outputs two harmonized, multi-sample VCFs — one from GATK and one from bcftools. Both subworkflows run in parallel. Additionally, summary statistics are produced such as transition/transversion ratios that can be used to judge the quality of/improvement in the called variants and whether subsequent rounds of BQSR could be beneficial. Finally, in the subworkflow `intersect_variants`, the two callsets are intersected using [`bcftools -isec`](https://samtools.github.io/bcftools/bcftools.html#isec) to retain only sites called by both subworkflows and filtered using [`vcftools -exclude`](https://vcftools.sourceforge.net/man_latest.html) (using parameter -include_scaffolds or -exclude_scaffolds), e.g. to exclude mitochondrial or gonosomal scaffolds or only include autosomal scaffolds, and [`vcftools -thin`](https://vcftools.sourceforge.net/man_latest.html) (--remove-filtered-all –remove-indels –maf 0.025 –recode –recode-INFO-all –max-missing 0.75) producing the final variant set.

<details markdown="1">
<summary>Output files</summary>

- `variant_calling/`
  - `bcftools/`: directory containing the called variants as vcf and tbi files for each interval, the individual bam files per sample (in the subdirectory `bam/`) and the merged vcf and tbi file (in the subdirectory `merged/`).
  - `gatk/`: directory containing the called variants for each individuals as vcf and tbi files for each interval, the merged vcf and tbi file (in the subdirectory `merged/`) and a text file with variant statistics (in the subfolder `stats/`).

</details>

## Relatedness Estimation Section

The fourth section, **relatedness estimation**, uses the final variant set to produce robust estimates of pairwise relatedness suitable for low-coverage whole-genome sequencing data. It infers relatedness and inbreeding from genotype likelihood data using ANGSD's [`NgsRelatev2`](https://github.com/ANGSD/NgsRelate). The final output is a pairwise relatedness matrix.

<details markdown="1">
<summary>Output files</summary>

- `relatedness_estimation/`
  - `angsd_ngsrelate/`: directory containing the tab-delimited file with the pairwise relatedness matrix of the analyzed samples.
  - `excluded/`: directory containing the vcf file that was filtered to only included variants of the scaffolds defined to be included.
  - `include_contigs/`: directory containing bed file with the scaffolds to define which variants to use for relatedness estimation.
  - `intersection/`: directory containing the vcf files produced by bcftools `isec`, a README.txt to explain the content of the four different vcf files, and a text file with the sites.
  - `thinned/`: a directory containing the vcf file with the filtered variant data used as input for the relatedness estimation.

</details>

### MultiQC

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results from the preprocessing and bootstrapping section are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastP. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.
An excellent overview of MultiQC outputs and how to interpret the various plots and values is provided in the [nf-core/eager](https://nf-co.re/eager/2.5.3/docs/output/#multiqc-report) pipeline.

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file, which includes the most important overview statistics and can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

### Pipeline information

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>
