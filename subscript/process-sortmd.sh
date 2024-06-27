#!/bin/bash
#
## Variable Definition Section
FQ_DIR=$(pwd)
#
cd bwa_bam
BWA_DIR=$(pwd)
cd $FQ_DIR
#
cd txt_files
TXT_DIR=$(pwd)
cd $FQ_DIR
#
cd dedup_bam
MD_DIR=$(pwd)
cd $FQ_DIR
#
if [ "$FQ_IN" = "NISC" ];
then
# Capture the first row, first column value of the table which will give the array the file name and allow dynamic usage
FN=( $(awk 'NR==1{print $1; exit}' "$TXT_DIR"/final_header.txt ) )
#
# Use the table to generate the sorted file names and push that to the samples text files folder
#
awk '{print "sort_"$1"."$4".bam"}' "${FN[0]}".runs.table | paste -sd " " - &> "$TXT_DIR"/merge_list.tmp
IFS=,$'\n' read -d '' -r -a lsfiles < "$TXT_DIR"/merge_list.tmp
#
#
echo "cd "$BWA_DIR"; samtools merge /lscratch/\$SLURM_JOB_ID/"${FN[0]}".bam ${lsfiles[0]} && gatk MarkDuplicates I=/lscratch/\$SLURM_JOB_ID/"${FN[0]}".bam O="$MD_DIR"/"${FN[0]}".sort.md.bam M="$TXT_DIR"/"${FN[0]}".sort.md.metrics.txt REMOVE_DUPLICATES=false ASSUME_SORTED=true TMP_DIR=/lscratch/\$SLURM_JOB_ID && samtools index "$MD_DIR"/"${FN[0]}".sort.md.bam" >> "$homedir"/mergeDedup.swarm
#
# Now if the original inputs were SRA data
elif [ "$FQ_IN" = "SRA" ];
then
FN=( $(awk 'NR==1{print $1; exit}' *.runs.table ) )
#
# Generate the sorted file names and push that to the samples text folder. This allows for singular SRR runs or multiple runs
awk '{print "sort_"$1"."$4".bam"}' "${FN[0]}".runs.table | paste -sd " " - &> "$TXT_DIR"/merge_list.tmp
IFS=,$'\n' read -d '' -r -a lsfiles < "$TXT_DIR"/merge_list.tmp
#
echo "cd "$BWA_DIR"; samtools merge /lscratch/\$SLURM_JOB_ID/"${FN[0]}".bam ${lsfiles[0]} && gatk MarkDuplicates I=/lscratch/\$SLURM_JOB_ID/"${FN[0]}".bam O="$MD_DIR"/"${FN[0]}".sort.md.bam M="$TXT_DIR"/"${FN[0]}".sort.md.metrics.txt REMOVE_DUPLICATES=false ASSUME_SORTED=true TMP_DIR=/lscratch/\$SLURM_JOB_ID && samtools index "$MD_DIR"/"${FN[0]}".sort.md.bam" >> "$homedir"/mergeDedup.swarm
#
else
  exit 1
fi
