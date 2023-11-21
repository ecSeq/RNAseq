#!/usr/bin/env nextflow

// PRINT HELP AND EXIT
if(params.help){
    include { printHelp } from './lib/functions.nf'
    printHelp()
}

// PRINT VERSION AND EXIT
if(params.version){
    include { printVersion } from './lib/functions.nf'
    printVersion()
}

// DEFINE PATHS # these are strings which are used to define input Channels,
// but they are specified here as they may be referenced in LOGGING
fasta = file("${params.reference}", checkIfExists: true, glob: false)
gtf = file("${params.features}", checkIfExists: true, glob: false)
samples = file("${params.samples}", checkIfExists: true, glob: false)
reads_path = params.SE ? "${params.input}/*.fastq.gz" : "${params.input}/*_{1,2}.fastq.gz"

// We should write an error-handler to check if "control" has been specified
if(!params.control){exit 1, "ERROR: please specify a control group with --control"}
// We should also write an error-handler to check if "control" is present in the "samples" file
count = 0
samples
    .readLines()
    .each { line ->
        group = line.toString().tokenize('\t').get(1)
        if(group == params.control) { count++ }
    }

// now test that control group was present
if(count == 0) {exit 1, "ERROR: --control group does not match with groups specified in: ${params.samples}"}


// PRINT STANDARD LOGGING INFO
include { printLogging } from './lib/functions.nf'
printLogging()


////////////////////
// STAGE CHANNELS //
////////////////////

/*
 *   Channels are where you define the input for the different
 *    processes which make up a pipeline. Channels indicate
 *    the flow of data, i.e. the "route" that a file will take.
 */

// STAGE BAM FILES FROM TEST PROFILE # this establishes the test data to use with -profile test
if ( workflow.profile.tokenize(",").contains("test") ){

        include { check_test_data } from './lib/functions.nf' params(readPaths: params.readPaths, singleEnd: params.SE)
        READS = check_test_data(params.readPaths, params.SE)

} else {

    // STAGE READS CHANNELS # this defines the normal input when test profile is not in use
    READS = Channel
        .fromFilePairs(reads_path, size: params.SE ? 1 : 2)
        .ifEmpty{ exit 1, """Cannot find valid read files in dir: ${params.input}
        The pipeline will expect PE reads in compressed *_{1,2}.fastq.gz format
        unless you have specified the --SE parameter in which case it is *.fastq.gz"""}
        .map{ tuple(it[0], it[1]) }
        .take(params.take.toInteger())

}



////////////////////
// BEGIN PIPELINE //
////////////////////

/*
 *   Workflows are where you define how different processes link together. They
 *    may be modularised into "sub-workflows" which must be named eg. 'RNAseq'
 *    and there must always be one MAIN workflow to link them together, which
 *    is always unnamed.
 */

include { RNAseq } from "${projectDir}/modules/workflow"

// MAIN WORKFLOW 
workflow {

    // call sub-workflows eg. WORKFLOW(Channel1, Channel2, Channel3, etc.)
    main:
        RNAseq(READS, fasta, gtf, samples)

}


//////////////////
// END PIPELINE //
//////////////////

// WORKFLOW TRACING # what to display when the pipeline finishes
// eg. with errors
workflow.onError {
    log.info "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}

// eg. in general
include { printSummary } from './lib/functions.nf'
workflow.onComplete {

    printSummary()

    // run a small clean-up script to remove "work" directory after successful completion 
    if (!params.debug && workflow.success) {
        ["bash", "${baseDir}/bin/clean.sh", "${workflow.sessionId}"].execute() }
}
