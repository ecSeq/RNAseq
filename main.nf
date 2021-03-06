#!/usr/bin/env nextflow

// ENABLE DSL2
nextflow.enable.dsl=2

// PRINT HELP AND EXIT
if(params.help){
    println """\

         ===========================================
          E C S E Q - R N A s e q   P I P E L I N E
         ===========================================
         ~ version ${workflow.manifest.version}

         Usage: 
              nextflow run ecseq/rnaseq [OPTIONS]...

         Options: GENERAL
              --input [path/to/input/dir]     [REQUIRED] Provide the directory containing fastq file(s) in "*_{1,2}.fastq.gz" format

              --reference [path/to/ref.fa]    [REQUIRED] Provide the path to the reference genome file in fasta format

              --features [path/to/ref.gtf]    [REQUIRED] Provide the path to the feature annotation file in gtf format

              --samples [path/to/samples.tsv] [REQUIRED] Provide the path to the samplesheet file required for DESeq2 in tsv format

              --output [STR]                  A string that can be given to name the output directory. [default: "."]


         Options: MODIFIERS
              --control [STR]                 [REQUIRED] Specify a group from the samplesheet to act as the control group

              --SE                            Indicate to the pipeline whether fastq files are SE reads in "*.fastq.gz" format. [default: off]

              --bamQC                         Generate bamQC report of alignments. [default: off]


         Options: ADDITIONAL
              --help                          Display this help information and exit
              --version                       Display the current pipeline version and exit
              --debug                         Run the pipeline in debug mode    


         Example: 
              nextflow run ecseq/rnaseq \
              --input /path/to/input/dir --reference /path/to/genome.fa \
              --samples /path/to/samples.tsv --features /path/to/features.gtf \
              --bamQC --keepReads

    """
    ["bash", "${baseDir}/bin/clean.sh", "${workflow.sessionId}"].execute()
    exit 0
}

// PRINT VERSION AND EXIT
if(params.version){
    println """\
         ===========================================
          E C S E Q - R N A s e q   P I P E L I N E
         ===========================================
         ~ version ${workflow.manifest.version}
    """
    ["bash", "${baseDir}/bin/clean.sh", "${workflow.sessionId}"].execute()
    exit 0
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
log.info ""
log.info "         ==========================================="
log.info "          E C S E Q - R N A s e q   P I P E L I N E"
if(params.debug){
log.info "         (debug mode enabled)"
log.info "         ===========================================" }
else {
log.info "         ===========================================" }
log.info "         ~ version ${workflow.manifest.version}"
log.info ""
log.info "         input dir    : ${workflow.profile.tokenize(",").contains("test") ? "-" : "${reads_path}"}"
log.info "         reference    : ${params.reference}"
log.info "         features     : ${params.features}"
log.info "         samples      : ${params.samples}"
log.info "         control      : ${params.control}"
log.info "         output dir   : ${params.output}"
log.info "         mode         : ${params.SE ? "single-end" : "paired-end"}"
log.info "         QC options   : ${params.bamQC ? "bamQC" : "-"}"
log.info ""
log.info "         ================================================"
log.info "         RUN NAME: ${workflow.runName}"
log.info ""



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
        .ifEmpty{ exit 1, "ERROR: cannot find valid read files in dir: ${params.input}\n \
        The pipeline will expect PE reads in compressed *_{1,2}.fastq.gz format\n \
        unless you have specified the --SE parameter in which case it is *.fastq.gz"}
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

// INCLUDES # here you must give the relevant process files from the lib directory 
include {STAR_index;STAR;bamQC;featureCounts;DESeq2} from './lib/process.nf' params(params)

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
workflow.onComplete {

    log.info ""
    log.info "         Pipeline execution summary"
    log.info "         ---------------------------"
    log.info "         Name         : ${workflow.runName}${workflow.resume ? " (resumed)" : ""}"
    log.info "         Profile      : ${workflow.profile}"
    log.info "         Launch dir   : ${workflow.launchDir}"    
    log.info "         Work dir     : ${workflow.workDir} ${!params.debug && workflow.success ? "(cleared)" : "" }"
    log.info "         Status       : ${workflow.success ? "success" : "failed"}"
    log.info "         Error report : ${workflow.errorReport ?: "-"}"
    log.info ""

    // run a small clean-up script to remove "work" directory after successful completion 
    if (!params.debug && workflow.success) {
        ["bash", "${baseDir}/bin/clean.sh", "${workflow.sessionId}"].execute() }
}
