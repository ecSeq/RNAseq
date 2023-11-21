#!/usr/bin/env nextflow
// This file for loading custom functions into the main.nf script (separated for portability)

// print --help information
def printHelp() {
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

// print --version information
def printVersion() {
    println """\
         ===========================================
          E C S E Q - R N A s e q   P I P E L I N E
         ===========================================
         ~ version ${workflow.manifest.version}
    """
    ["bash", "${baseDir}/bin/clean.sh", "${workflow.sessionId}"].execute()
    exit 0
}

// print pipeline initiation logging
def printLogging() {
    log.info """
         ===========================================
          E C S E Q - R N A s e q   P I P E L I N E """
    if(params.debug){ log.info "         (debug mode enabled)" }
    log.info """         ===========================================
         ~ version ${workflow.manifest.version}

         input dir    : ${workflow.profile.tokenize(',').contains('test') ? '-' : reads_path}
         reference    : ${params.reference}
         features     : ${params.features}
         samples      : ${params.samples}
         control      : ${params.control}
         output dir   : ${params.output}
         mode         : ${params.SE ? 'single-end' : 'paired-end'}
         QC options   : ${params.bamQC ? 'bamQC' : ''}

         ===========================================
         RUN NAME: ${workflow.runName}

    """
}

// print pipeline execution summary
def printSummary() {

    log.info """
         Pipeline execution summary
         ---------------------------
         Name         : ${workflow.runName}${workflow.resume ? ' (resumed)' : ''}
         Profile      : ${workflow.profile}
         Launch dir   : ${workflow.launchDir}
         Work dir     : ${workflow.workDir} ${params.debug || !workflow.success ? '' : '(cleared)' }
         Status       : ${workflow.success ? 'success' : 'failed'}
         Error report : ${workflow.errorReport ?: '-'}
    """
}

// FUNCTION TO LOAD DATASETS IN TEST PROFILE
def check_test_data(readPaths, singleEnd) {

    // Set READS testdata
    if( singleEnd ) {
        READS = Channel.from(readPaths)
            .map { row -> [ row[0], file(row[1][0]) ] }
            .ifEmpty { exit 1, "test profile readPaths was empty - no input files supplied" }
    } else {
        READS = Channel.from(readPaths)
            .map { row -> [ row[0], [file(row[1][0]),file(row[1][1])] ] }
            .ifEmpty { exit 1, "test profile readPaths was empty - no input files supplied" }
    }
    // Return READS channel
    return READS
}
