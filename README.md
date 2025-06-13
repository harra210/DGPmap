# NHGRI Dog Genome Project WGS mapping Pipeline

This repository contains the pipeline toolkit which WGS samples are processed at NHGRI's Dog Genome Project. This pipeline is heavily inspired by the pipeline created by [Jeff Kidd for Dog10K](https://github.com/jmkidd/dogmap). It consists of three main sections: Alignment, Jointcalling, and cohort evaluation. This repository is designed to be a stand alone package that can be modified to a user's particular situation. This pipeline is designed to be executed in a unix environment using our [cluster environment](https://hpc.nih.gov/) using a global shared disk file system. This pipeline is designed in shell and is designed to work via subshells in order to iterate through subdirectories to excute commands and write a properly formatted swarm file which can then be submitted to the cluster at the very end after user confirmation that the swarmfile syntax looks correct.

## Required files to download

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

## Software required

* bwa-mem 2 version 2.2.1
* samtools version ≥ 1.9
* sratoolkit verion ≥ 3.0.0
* tabix and bgzip from htslib version ≥ 1.9
* bcftools ≥ 1.16
* vcflib 1.0.3
* Python ≥ 3.7
* R ≥ 4.3.2
* GATK version ≤ 4.2.0.0 or ≥ 4.6.0.0
>[!CAUTION]
> Running this pipeline using a version of GATK between 4.2.0.0 and 4.6.0.0 will result with a high likelihood of corrupted CRAMs and missing calls in joint calls notated in a non-VCF standard format of 0/0 DP=0 instead of ./.

## Further Reading
Each section will have respective details on how the pipeline will work and what inputs are required from the user in order to keep this page brief.

## Discussions
Report issues or enhancement ideas using the issue tracker with additional discussions happening within pull requests. Please note that only one person is coding this pipeline with other responsiblities so responses may be delayed.
