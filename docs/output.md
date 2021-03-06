# ecSeq-RNAseq Output
This document describes the output produced by the pipeline.

## Pipeline overview
The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

* [Read Alignment](#read-alignment) - mapping RNAseq reads with STAR
* [Quality Control](#quality-control) - generating bamQC reports
* [Differential Expression](#differential-expression) - performing differential expression analysis
* [Pipeline Info](#pipeline-info) - reports from nextflow about the pipeline run


## Read Alignment
RNAseq reads will be aligned with STAR.

**Output directory: `./mapping`**


## Quality Control
Following alignment, the pipeline will optionally generate bamQC reports for each BAM file.

**Output directory: `./mapping/<sample>`**


## Differential Expression
Feature counts will be estimated with Subread featureCounts, prior to differential expression analysis with DESeq2. The final table of results `output.tsv` will be located directly in the user-specified output directory.

**Output directory: `./`**


## Pipeline Info
Nextflow has several built-in reporting tools that give information about the pipeline run.

**Output directory: `./`**

* `dag.svg`
  * DAG graph giving a diagrammatic view of the pipeline run.
  * NB: If [Graphviz](http://www.graphviz.org/) was not installed when running the pipeline, this file will be in [DOT format](http://www.graphviz.org/content/dot-language) instead of SVG.
* `report.html`
  * Nextflow report describing parameters, computational resource usage and task bash commands used.
* `timeline.html`
  * A waterfall timeline plot showing the running times of the workflow tasks.
* `trace.txt`
  * A text file with machine-readable statistics about every task executed in the pipeline.