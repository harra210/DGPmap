#!/bin/bash
#
# This is an internal script for DGPmap to apply the BQSR tables and output per chromsome.
#
FQ_DIR=$(pwd)
#
## Variable Definition Section
#
cd txt_files
TXT_DIR=$(pwd)
cd $FQ_DIR
#
cd dedup_bam
MD_DIR=$(pwd)
cd $FQ_DIR
#
cd BQSR/tables
TBL_DIR=$(pwd)
cd $FQ_DIR
#
cd BQSR_chr
BQSRchr_DIR=$(pwd)
cd $FQ_DIR
#
##
# 
FN=( $(awk 'NR==1{print$1; exit}' "$TXT_DIR"/final_header.txt ) )
#
# Write the Y chromosome first to the swarm file as the swarm will be bunched so to take advantage of the swarm functionality, use the non-loop command to share deduplicated bam file across all bunched jobs.
#
echo "cd "$MD_DIR"; cp "${FN[0]}".sort.md.bam /lscratch/\$SLURM_JOB_ID; cp "${FN[0]}".sort.md.bam.bai /lscratch/\$SLURM_JOB_ID && gatk --java-options \"-Xmx4G\" ApplyBQSR -R "$CanFam4_Ref" --tmp-dir /lscratch/\$SLURM_JOB_ID -I lscratch/\$SLURM_JOB_ID/"${FN[0]}".sort.md.bam -XL "$interval_list" -O "$BQSR_chr"/"${FN[0]}".chrY.BQSR.bam -bqsr-recal-file "$TBL_DIR"/"${FN[0]}".bqsrgathered.reports.list --preserve-qscores-less-than 6 --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30" >> "$homedir"/bqsr_ApplyBQSR.swarm
#
# Generate a while loop for the intervals to create the per chromosome BAM files (possible to use CRAM here?)
while read g
do
	echo "cd "$BQSRchr_DIR"; gatk --java-options \"-Xmx4G\" ApplyBQSR -R "$CanFam4_Ref" --tmp-dir /lscratch/\$SLURM_JOB_ID -I /lscratch/\$SLURM_JOB_ID/"${FN[0]}".sort.md.bam -L "$g" -O "$BQSRchr_DIR"/"${FN[0]}"."$g".BQSR.bam --bqsr-recal-file "$TBL_DIR"/"${FN[0]}".bqsrgathered.reports.list --preserve-qscores-less-than 6 --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30" >> "$homedir"/bqsr_ApplyBQSR.swarm
done < "$interval_list"
#
# This completes the ApplyBQSR script
