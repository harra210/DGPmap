#!/bin/bash
#
# This is an internal script for DGPmap to generate BQSR tables per chromosome.
#
## Variable Definition Section
FQ_DIR=$(pwd)
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
## End Variable Definitions
#
# Capture the first row, first column value of this table to give the array the file name and allow dynamic usage
FN=( $(awk 'NR==1{print $1; exit}' "$TXT_DIR"/final_header.txt) )
#
# Now, to generate the swarm file using a while loop to create files chr1 - chrX
#
while read g
do
	echo "cd "$MD_DIR"; gatk --java-options \"-Xmx4G\" BaseRecalibrator -R "$CanFam4_Ref" --tmp-dir /lscratch/\$SLURM_JOB_ID -I "${FN[0]}".sort.md.bam --known-sites "$knownsite" -L "$g" -O "$TBL_DIR"/"${FN[0]}"_"$g"_recal.table" >> "$homedir"/bqsr_BaseRecalibrator.swarm
done < "$interval_list"
#
echo "cd "$MD_DIR"; gatk --java-options \"-Xmx4G\" BaseRecalibrator -R "$CanFam4_Ref" --tmp-dir /lscratch/\$SLURM_JOB_ID -I "${FN[0]}".sort.md.bam --known-sites "$knownsite" -XL "$interval_list" -O "$TBL_DIR"/"${FN[0]}"_chrY_recal.table" >> "$homedir"/bqsr_BaseRecalibrator.swarm
#
## This completes the Recalibration Table script
