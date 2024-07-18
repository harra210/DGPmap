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

The first step in the process is to create a reference table for the programs to pull data from in order to be able to place into their respective programs. This is mainly important for the alignment step with BWAMEM2 but becomes useful downstream. There are two main methods in which tables would be generated, based on the initial input flag for fastQ's (Illumina or SRA).

For FastQ's where individual reads are Illumina FastQ's the following is the expected hierarchy:
```
Parent directory
	|
	 -> Folder containing fastQ reads
```

The expected FastQ file name structure for the pipeline to work on is as follows:
```
*_S1_L001_R1_001.fastq.gz
```

Illumina caveat: The table-generation script does rely on 1 pair of FastQ's per subfolder in the parent directory so that it can iterate correctly
```bash
#Line 47 of table-generation.sh
find $PWD -maxdepth 2 -name "*L004_R1_001.fastq.gz" -printf '%h\n' &> "$tmpdir"/single_dir.tmp # Note: This will need to be modified to garner a unique result
```
For FastQ's downloaded from the SRA database:
```
Parent directory
	|
	 -> Biosample Folder name
		|
		 -> Experiment Run Folders
```
For runs created by the SRA database the table generation should be able to handle both single run samples as well as multi-run samples under the auspice of using the BioSample number as the means of naming the Sample Name.

## Step 2 - Alignment with BWA-MEM2

Example command:
```
bwa-mem2 mem -K 100000000  -t NUM_THREADS -Y \
-R READGROUPINFO \
PATH_TO_UU_Cfam_GSD_1.0_ROSY.fa_with_index \
fq1.gz fq2.gz | samtools view -h | samtools sort -@ $SLURM_CPUS_PER_TASK -T tmp/dir/foo -o "$FQ_DIR"/bwa_bam/sort_foo.bar.bam 
```

This command is inspired by [Reiger et al. Functional equivalence of genome sequencing analysis](https://www.nature.com/articles/s41467-018-06159-4) and by 
discussions with colleagues. The -K option removes non-determinism in the fragment length distributions identified using different numbers of threads. -Y uses
soft clipping for supplementary alignments, which may aid down stream structural variation analyses. 

## Step 3 - Merging and Mark Duplicates

One may wonder why in Step 2, samples were processed through Samtools' sort tool, this was done as in this step Samtools merge tool requires coordinate sorted bam files as inputs. Merge will combine all of the supplied bam files and then sort them prior to output. From there the pipeline will use GATK MarkDuplicates to mark duplicate reads, but not remove them in the event the end user would like to be able to see these duplicates.

Example command:
```
samtools merge /tmp/sort_SAMPLE.bam [list of input bam files] && \
gatk MarkDuplicates \
I=/tmp/sort_SAMPLE.bam \
O=../dedup_bam/SAMPLE.sort.md.bam \
M=../txt_dir/SAMPLE.sort.md.metrics.txt \
REMOVE_DUPLICATES=false \
ASSUME_SORTED=true \
TMP_DIR=/tmp
```

## Step 4 - BQSR Step 1: BaseRecalibration

BQSR is run in parallel with 40 seperate jobs for chr1-chr38, chrX, and all other sequences together. The resulting recalibration data tables are then combined using GatherBQSRReports.

Example command:
```
gatk --java-options "-Xmx4G" BaseRecalibrator \
--tmp-dir /tmp \
-I ../dedup_bam/SAMPLE.sort.md.bam \
-R PATH_TO_UU_Cfam_GSD_1.0_ROSY.fa \
--intervals INTERVAL_FILE_TO_PROCESS \
--known-sites PATH_TO_UU_Cfam_GSD_1.0.BQSR.DB.bed.gz \
-O ../BQSR/tables/SAMPLE_INTERVAL_recal.table
```
Then when completed the recalibration data is combined:

```
gatk --java-options "-Xmx6G" GatherBQSRReports \
--input ../BQSR/tables/SAMPLE_INTERVAL-1_recal.table \
--input ../BQSR/tables/SAMPLE_INTERVAL-2_recal.table \
--input ../BQSR/tables/SAMPLE_INTERVAL-3_recal.table \
...
--output ../BQSR/tables/SAMPLE.bqsrgathered.reports.list
```

## Step 5 - BQSR Step 2: ApplyBQSR

The base calibration is then performed with ApplyBQSR based on the combined reports. This is done in parallel with
41 separate jobs for chr1-chr38, chrX, all other sequences together, and unmapped reads.
The recalibrated qualities scores are discretized into 3 bins with low quality values unchanged.

Example command:
```
gatk --java-options "-Xmx4G" ApplyBQSR \
--tmp-dir /tmp \
-I ../dedup_bam/SAMPLE.sort.md.bam \
-R PATH_TO_UU_Cfam_GSD_1.0_ROSY.fa \
-O ../BQSR/BQSR_chr/SAMPLE.INTERVAL.bqsr.bam \
--intervals INTERVAL_FILE_TO_PROCESS \
--bqsr-recal-file ../BQSR/tables/SAMPLE.bqsrgathered.reports.list \
--preserve-qscores-less-than 6 \
--static-quantized-quals 10 \
--static-quantized-quals 20 \
--static-quantized-quals 30
```
The 41 recalibrated BAMs are then combined using GatherBamFiles

```
gatk --java-options "-Xmx6G" GatherBamFiles \
--CREATE_INDEX true \
-I ../BQSR/BQSR_chr/SAMPLE.INTERVAL-1.bqsr.bam \
-I ../BQSR/BQSR_chr/SAMPLE.INTERVAL-2.bqsr.bam \
-I ../BQSR/BQSR_chr/SAMPLE.INTERVAL-3.bqsr.bam \
...
-O ../BQSR/SAMPLE.BQSR.bam
```

## Step 6 - Create gVCF with HaplotypeCaller

HaplotypeCaller is then run to create a per-sample GVCF file for subsequent cohort
short variant identification.  This is run in parallel across 39 separate jobs 
for chr1-chr38, chrX, and all other sequences together.  The resulting GVCF files
are then combined using GatherVcfs.

Example command:
```
gatk --java-options "-Xmx4G" HaplotypeCaller \
--tmp-dir /tmp \
-R PATH_TO_UU_Cfam_GSD_1.0_ROSY.fa \
-I ../BQSR/SAMPLE.BQSR.bam \
--intervals INTERVAL_FILE_TO_PROCESS \
-O ../gVCF/HC/SAMPLE_INTERVAL.g.vcf.gz \
-ERC GVCF 
```
The GVCFs are combined using:

```
gatk --java-options "-Xmx6G" GatherVcfs \
-I ../gVCF/HC/SAMPLE_INTERVAL-1.g.vcf.gz \
-I ../gVCF/HC/SAMPLE_INTERVAL-2.g.vcf.gz \
-I ../gVCF/HC/SAMPLE_INTERVAL-3.g.vcf.gz \
...
-O ../gVCF/SAMPLE_g.vcf.gz
 ```

## Step 7 - DGPmap-postprocessing - Conversion, Recompression, Metrics

Once the main computation portion is now complete, postprocessing can now begin. This involves conversion of the BAM files to CRAM, the recompression of gVCF files and the collection of various metrics required to generate stats per sample.

# Flags
```
        -d, --directory; Parent directory of pipeline processed files. (Req'd)
        -f, --file; Space-delimited file containing Original name followed by new name.
        -s, --swarm-name; Set the base swarm name, any job specific names will be added on job submission. (Req'd, assuming swarm environment)
        -h, --help; Print this message and exit.
```

This script allows for the user to input a two column space or tab-delimited file with column 1 detailing the original name, which should match the sample label in the runs table (ex: **SAMPLE**.runs.table) and column 2 containing the respective changed name. This file does not have to line match, be sorted, nor contain the exact number of samples being renamed. It is very tolerant of number of samples. The naming conversion uses an awk indexing array to search for the value of the converting sample within the input conversion file, match it and output the correct conversion syntax for GATK PrintReads.

In the event no file is supplied, the user is then prompted to supply a new name, to wit the user can opt to input the old name.

# Conversion of BAM to CRAM
```
gatk --java-options "-Xmx6G" PrintReads \
-R PATH_TO_UU_Cfam_GSD_1.0_ROSY.fa \
-I ../BQSR/SAMPLE.BQSR.bam \
-O /CHANGEDNAME.BQSR.cram \
-OBM true \
-OBI false \
&& samtools index -@ threads -c /CHANGEDNAME.BQSR.cram \
&& cp /txt_dir/SAMPLE.sort.md.metrics.txt /txt_dir/CHANGEDNAME.sort.md.metrics.txt
```

This outputs an md5 file to be generated, and a rename of the markduplicates metrics file which will be used later generating stats.

# gVCF Recompression
```
zcat ../gVCF/SAMPLE_g.vcf.gz | bgzip -@6 -l 9 -c > /CHANGEDNAME.g.vcf.gz && \
tabix -p vcf -f /CHANGEDNAME.g.vcf.gz && \
md5sum /CHANGEDNAME.g.vcf.gz > /CHANGEDNAME.g.vcf.gz.md5
```
