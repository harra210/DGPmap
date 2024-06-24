# DGPmap
NHGRI Dog Genome Project Whole Genome Sequencing Alignment Pipeline

#### Last updated 24 June 2024

This is a working pipeline for processing Whole Genome Sequencing data for the NHGRI Dog Genome Project working group. This approach is not meant to be definitive and is a work in progress.

## Software

```
bwa-mem 2 version 2.2.1
samtools version ≥ 1.9
GATK version 4.4.0.0
sratoolkit verion ≥ 3.0.0
tabix and bgzip from htslib version ≥ 1.9
```
## Conceptual overview of pipeline

The main script is designed to be executed in a unix environment using our cluster environment using global shared disk file system. This pipeline is designed in shell and is designed to work via subshells in order to iterate through subdirectories to excute commands and write a properly formatted swarm file which can then be submitted to the cluster at the very end after user confirmation that the swarmfile syntax looks correct.

To invoke the main DGPmap script:

```
DGPmap.sh [OPTION...]
	-i, -input                       Parent input directory of the files initially looking to process.
	-f, --fastq                      Set input fastq style format (default: NISC).
	-s, --swarm-name                 Set the base swarm name for the pipeline, job specific names will be added on job submission.
	-h, --help                       Print this message and exit.
```
Only the fastq flag is optional and the output of the files are to their respective directories. This comes into play later once QC is performed.

**Note:** If you are using a system with faste storage options, or a Google Spark Cluster then using GATK's spark set of tools would be more beneficial to you. 

All paths given are used as placeholders and should be modiefied for personal use.

The companion SRA download package is a script that also integrates with our cluster environment but will intake a tab delimited file or manually inputted SRA Run numbers, download, unpackage the fastQ files into a format ready for the main script to be able process.

