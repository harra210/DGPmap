#!/bin/bash
#
## This is an internal script for the post-processing of the DGPmap pipeline that will rename and recompress gVCF files maximally.
#
## Variable definition section
#
FQ_DIR=$(pwd)
#
cd txt_files
TXT_DIR=$(pwd)
cd "$FQ_DIR"
#
cd gVCF
gVCF_DIR=$(pwd)
cd "$FQ_DIR"
#
## Section to define the original name. Since the rename string has already been defined by the BQSR section no if/then statement required to set the variable FN.
FN=( $(awk 'NR==1{print $1; exit}' *.runs.table ) )
#
cd "$TXT_DIR"
IFS=,$'\n' read -d '' -r -a changedname < rename.txt # import the previously created rename file into an array for use as a variable
#
cd "$gVCF_DIR"
gVCF_Gb=( $(stat --printf="%s" "${FN[0]}"_g.vcf.gz | awk '{print $1/1024/1024/1024; }' ) )
#
echo "$gVCF_Gb" > "$TXT_DIR"/gVCFsize_old.txt
#
echo "cd "$gVCF_DIR"; cp "${FN[0]}"_g.vcf.gz /lscratch/\$SLURM_JOB_ID/"${FN[0]}"_g.vcf.gz && zcat /lscratch/\$SLURM_JOB_ID/"${FN[0]}"_g.vcf.gz | bgzip -@6 -l 9 -c > "$FQ_DIR"/"${changedname[0]}".g.vcf.gz && tabix -p vcf -f "$FQ_DIR"/"${changedname[0]}".g.vcf.gz && md5sum "$FQ_DIR"/"${changedname[0]}".g.vcf.gz > "$FQ_DIR"/"${changedname[0]}".g.vcf.gz.md5" >> "$homedir"/postprocess-gVCFRecompress.swarm
