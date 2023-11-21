#!/usr/bin/env nextflow

// INCLUDES # here you must give the relevant processes from the modules/process directory 
include { STAR_index;STAR;bamQC;featureCounts;DESeq2 } from "${projectDir}/modules/process"

// SUB-WORKFLOWS
workflow 'RNAseq' {

    // take the initial Channels and paths
    take:
        READS
        fasta
        gtf
        samples

    // here we define the structure of our workflow i.e. how the different processes lead into each other
    // eg. process(input1, input2, etc.)
    // eg. process.out[0], process.out[1], etc.
    // index numbers [0],[1],etc. refer to different outputs defined for processes in process.nf
    // ALWAYS PAY ATTENTION TO CARDINALITY!!
    main:
        // first we will perform indexing
        STAR_index(fasta)
        // then we can begin the mapping 
        STAR(READS,STAR_index.out)

        // now we have a series of alignments from STAR e.g. [sample, /path/to/sample.bam]
        // this is exactly what we need for bamQC
        bamQC(STAR.out[0])

        // however, featureCounts needs a collection of all bams and has no use anymore for sample names
        featureCounts_input = STAR.out[0].map{ it[1] }.collect()
        // eg. [/path/to/sample.bam, ... /path/to/sampleN.bam]

        // now we can provide the new channel to featureCounts alongside the gtf file
        featureCounts(featureCounts_input,gtf)

        // finally we are able to run DESeq2 on the featureCount table alongside the samplesheet
        DESeq2(featureCounts.out,samples)

}