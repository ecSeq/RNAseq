#!/usr/bin/env nextflow
// This file defines individual processes (separated for portability)


// index the genome
process "STAR_index" {

    label "low"
    label "finish"

    maxForks "${params.fork}".toInteger()

    input:
    path fasta

    output:
    path "genome/STAR"

    script:
    """
    mkdir genome genome/STAR
    STAR --runMode genomeGenerate \\
    --genomeDir genome/STAR/ \\
    --outFileNamePrefix genome/STAR/${fasta.baseName}. \\
    --genomeFastaFiles ${fasta} --runThreadN ${task.cpus}
    """
}



// align the reads to the genome
process "STAR" {

    label "low"
    label "finish"
    tag "$sample"

    maxForks "${params.fork}".toInteger()

    publishDir "${params.output}", pattern: "mapping/${sample}/*", mode: 'copy'

    input:
    tuple val(sample), path(reads)
    // eg. [sample, [/path/to/sample_1.fastq, /path/to/sample_2.fastq]]
    // eg. [sample, /path/to/sample.fastq]
    path index

    output:
    tuple val(sample), path("mapping/${sample}/${sample}.Aligned.sortedByCoord.out.bam")
    // eg. [sample, /path/to/sample.Aligned.sortedByCoord.out.bam]
    path "mapping/${sample}/*"

    script:
    """
    mkdir mapping mapping/${sample}
    STAR --genomeDir ${index}/ --readFilesIn ${reads} \\
    --outFileNamePrefix mapping/${sample}/${sample}. --outSAMattributes All \\
    --runThreadN ${task.cpus} --outSAMtype BAM SortedByCoordinate \\
    --readFilesCommand zcat
    """
}


// perform QC of alignments
process "bamQC" {

    label "low"
    label "ignore"
    tag "$sample"

    maxForks "${params.fork}".toInteger()

    publishDir "${params.output}", pattern: "mapping/${sample}/*.{txt,pdf,log}", mode: 'move'

    input:
    tuple val(sample), path(bam)
    // eg. [sample, /path/to/sample.bam]

    output:
    path "mapping/${sample}/*.{txt,pdf,log}"

    when:
    params.bamQC

    script:
    """
    mkdir mapping mapping/${sample}
    qualimap bamqc -bam ${bam} -outdir mapping/${sample} -outformat pdf > mapping/${sample}/bamqc.log 2>&1
    """
}


// generate feature counts based on gene annotations
process "featureCounts" {

    label "low"
    label "finish"

    input:
    path bams
    // eg. [/path/to/sample1.bam, ... /path/to/sampleN.bam]
    path gtf

    output:
    path "featureCounts/counts.tsv"

    script:
    """
    mkdir featureCounts
    featureCounts -M --fraction -g gene_id -T ${task.cpus} -s 0 \\
    -a ${gtf} -o featureCounts/table.featureCounts ${bams} || exit \$?

    cut -f 1,7- featureCounts/table.featureCounts | 
    sed 's/.Aligned.sortedByCoord.out.bam//g' > featureCounts/counts.tsv
    """

}


// perform differential expression analysis
process "DESeq2" {

    label "low"
    label "finish"

    publishDir "${params.output}", pattern: "output.tsv", mode: 'move'
    publishDir "${params.output}", pattern: "DESeq2.log", mode: 'move'

    input:
    path counts
    path "samples.tsv"

    output:
    path "output.tsv"
    path "DESeq2.log"

    script:
    """
    Rscript ${baseDir}/bin/DESeq2.R ${counts} samples.tsv ${params.control} 2> DESeq2.log
    """
}