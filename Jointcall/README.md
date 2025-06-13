## DGPmap joint call pipeline

### Requirements

#### Files
```
 104K Feb 26 14:20 UU_Cfam_GSD_1.0_ROSY.fa.fai
 750M Feb 26 14:20 UU_Cfam_GSD_1.0_ROSY.fa.gz
  11M Mar  5 16:00 SRZ189891_722g.simp.header.CanineHDandAxiom_K9_HD.GSD_1.0.vcf.gz
 1.2M Mar  5 16:00 SRZ189891_722g.simp.header.CanineHDandAxiom_K9_HD.GSD_1.0.vcf.gz.tbi
```
#### Tools
* bcftools ≥ 1.16
* GATK version ≤ 4.2.0.0 or ≥ 4.6.0.0
>[!CAUTION]
> Running this pipeline using a version of GATK between 4.2.0.0 and 4.6.0.0 will cause missing calls in joint calls to be notated in a non-VCF standard format of 0/0 DP=0 instead of ./.
