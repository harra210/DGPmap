# DGPmap
NHGRI Dog Genome Project Whole Genome Sequencing Alignment Pipeline

#### Last updated 24 June 2024

This is a working pipeline for processing Whole Genome Sequencing data for the NHGRI Dog Genome Project working group.  It is heavily inspiredby the pipeline created by [Jeff Kidd for Dog10K](https://github.com/jmkidd/dogmap)

## Required files

Associated files for alignment can be found at https://kiddlabshare.med.umich.edu/public-data/UU_Cfam_GSD_1.0-Y/

```
total 2.4G
  87K Feb 26 14:19 canFam4.chromAlias.txt
  144 Feb 26 14:20 chrY.chromAlias.txt
 898M Mar  3 13:57 SRZ189891_722g.simp.GSD1.0.vcf.gz
 1.6M Mar  3 14:01 SRZ189891_722g.simp.GSD1.0.vcf.gz.tbi
  11M Mar  5 16:00 SRZ189891_722g.simp.header.Axiom_K9_HD.names.GSD_1.0.filter.vcf.gz
 1.2M Mar  5 16:00 SRZ189891_722g.simp.header.Axiom_K9_HD.names.GSD_1.0.filter.vcf.gz.tbi
  11M Mar  5 16:00 SRZ189891_722g.simp.header.CanineHDandAxiom_K9_HD.GSD_1.0.vcf.gz
 1.2M Mar  5 16:00 SRZ189891_722g.simp.header.CanineHDandAxiom_K9_HD.GSD_1.0.vcf.gz.tbi
 2.5M Mar  5 16:00 SRZ189891_722g.simp.header.CanineHD.names.GSD_1.0.filter.vcf.gz
 798K Mar  5 16:00 SRZ189891_722g.simp.header.CanineHD.names.GSD_1.0.filter.vcf.gz.tbi
 680M Mar 23 09:56 UU_Cfam_GSD_1.0.BQSR.DB.bed.gz
 2.0M Mar 23 09:55 UU_Cfam_GSD_1.0.BQSR.DB.bed.gz.tbi
 343K Feb 26 14:20 UU_Cfam_GSD_1.0_ROSY.dict
 104K Feb 26 14:20 UU_Cfam_GSD_1.0_ROSY.fa.fai
 750M Feb 26 14:20 UU_Cfam_GSD_1.0_ROSY.fa.gz
```
Specifics regarding the use of these files can be found [here](https://github.com/jmkidd/dogmap/blob/main/README.md)

## Software

```
bwa-mem 2 version 2.2.1
samtools version ≥ 1.9
GATK version 4.4.0.0
sratoolkit verion ≥ 3.0.0
tabix and bgzip from htslib version ≥ 1.9
```
## Conceptual overview of pipeline

The main script is designed to be executed in a unix environment using our [cluster environment](https://hpc.nih.gov/) using global shared disk file system. This pipeline is designed in shell and is designed to work via subshells in order to iterate through subdirectories to excute commands and write a properly formatted swarm file which can then be submitted to the cluster at the very end after user confirmation that the swarmfile syntax looks correct.

To invoke the main DGPmap script:

```
DGPmap.sh [OPTION...]
	-i, --input                      Parent input directory of the files initially looking to process.
	-f, --fastq                      Set input fastq style format (default: NISC).
	-s, --swarm-name                 Set the base swarm name for the pipeline, job specific names will be added on job submission.
	-h, --help                       Print this message and exit.
```
Only the fastq flag is optional and the output of the files are to their respective directories. This comes into play later once QC is performed.

**Note:** If you are using a system with fast storage options, or a Google Spark Cluster then using GATK's spark set of tools would be more beneficial to you. 

All paths given are used as placeholders and should be modiefied for personal use.

The companion SRA download package is a script that also integrates with our cluster environment but will intake a tab delimited file or manually inputted SRA Run numbers, download, unpackage the fastQ files into a format ready for the main script to be able process.

## Step 1 - Table Generation

The first step in the process is to create a reference table for the programs to pull data from in order to be able to place into their respective programs. This is mainly important for the alignment step with BWAMEM2 but becomes useful downstream. There are two main methods in which tables would be generated, based on the initial input flag for fastQ's (NISC or SRA)

For fastQ's where individual reads were Illumina fastq's typically lab processed for samples sequenced by our lab at NISC:
```
Parent directory
	|
	 -> Folder containing fastQ reads
```

For fastQ's downloaded from the SRA database:
```
Parent directory
	|
	 -> Biosample Folder name
		|
		 -> Experiment Run Folders
```
For runs created by the SRA database the table generation should be able to handle both single run samples as well as multi-run samples under the auspice of using the BioSample number as the means of naming the Sample Name.
