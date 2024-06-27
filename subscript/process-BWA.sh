#!/bin/bash
#
# This an internal script from the parent DGPmap.sh. This aligns fastq files using
# BWAMEM2 based from the tables generated previously.
#
## Variable Definition Section
FQ_DIR=$(pwd)
#
cd txt_files
TXT_DIR=$(pwd)
#
cd $FQ_DIR
#
# Here we read the table into its component variables.
# Table is broken down into SM LB ID PU R1 R2
# Sample (SM) will be the external sample name and will be renamed with the internal ID names at the very end.
# 
if [ "$FQ_IN" = "NISC" ];
then
# Capture the first row, first column value of the table which will give the array the file name and allow dynamic usage.
FN=( $(awk 'NR==1{print $1; exit}' "$TXT_DIR"/final_header.txt ) )
#
# Set the variables for the table
while read SM LB ID PU R1 R2 ; do
echo "cd "$FQ_DIR"; bwa-mem2 mem -K 100000000 -t \$SLURM_CPUS_PER_TASK -Y -R '@RG\\tID:"${ID}"\\tSO:coordinate\\tLB:"${LB}"\\tPL:ILLUMINA\\tSM:"${SM}"\\tPU:"${PU}"' "$CanFam4_Ref" "${R1}" "${R2}" | samtools view -h | samtools sort -@ \$SLURM_CPUS_PER_TASK -T /lscratch/\$SLURM_JOB_ID/"${LB}" -o "$FQ_DIR"/bwa_bam/sort_"${SM}"."${PU}".bam" >> "$homedir"/bwa_to_picard.swarm
done < "${FN[0]}".runs.table
#
elif [ "$FQ_IN" = "SRA" ];
then
FN=( $(awk 'NR==1{print $1; exit}' *.runs.table ) )
#
# Set the variables for the SRA table
while read SM LB ID PU R1 R2 ; do
echo "cd "$FQ_DIR"; bwa-mem2 mem -K 100000000 -t \$SLURM_CPUS_PER_TASK -Y -R '@RG\\tID:"${ID}"\\tSO:coordinate\\tLB:"${LB}"\\tPL:ILLUMINA\\tSM:"${SM}"\\tPU:"${PU}"' "$CanFam4_Ref" "${R1}" "${R2}" | samtools view -h | samtools sort -@ \$SLURM_CPUS_PER_TASK -T /lscratch/\$SLURM_JOB_ID/"${LB}" -o "$FQ_DIR"/bwa_bam/sort_"${SM}"."${PU}".bam" >> "$homedir"/bwa_to_picard.swarm
done < "${FN[0]}".runs.table
#
else
  exit 1
fi
