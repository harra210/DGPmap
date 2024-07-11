#!/bin/bash
#
## This is an internal script for the post-processing of the DGPmap pipeline that will convert BQSR bam files to CRAM files.
## This script will also rename the MarkDuplicate metrics file to reflect the newly renamed CRAM file for later when stats are collected.
#
## Variable definition section
#
FQ_DIR=$(pwd)
#
cd BQSR
BQSR_DIR=$(pwd)
cd "$FQ_DIR"
#
cd txt_files
TXT_DIR=$(pwd)
cd "$FQ_DIR"
#
FN=( $(awk 'NR==1{print $1; exit}' "$FQ_DIR"/*.runs.table ) )
#
## Section to define the original name to be renamed and for usage down the line.
#
## If file flag was originally passed
if [[ ! -z "$FILE_IN" ]]
then
	awk 'NR==FNR{a[$1]; next} {for (i in a) if (index($0, i)) print $2}' "$FQ_DIR"/*.runs.table "$FILE_IN" &> "$TXT_DIR"/rename.txt
elif [[ -z "$FILE_IN" ]]
then
	read -ep "What is the Ostrander Lab Name for this sample: "${FN[0]}"? " rename;
	echo ""$rename"" &> "$TXT_DIR"/rename.txt
#
else
	echo "Script error!"; exit 1
fi
#
IFS=,$'\n' read -d '' -r -a changedname < "$TXT_DIR"/rename.txt
## User prompt to create a name for the sample that will be overwritten through each iteration of the for loop.
#echo "What is the Ostrander Lab Name for this sample: "${FN[0]}" ?"
#$read -ep "New Sample Name: " changedname
#
#cd "$TXT_DIR"
#echo "$changedname" &> rename.txt
#
#cd "$FILE_DIR"
#
echo "cd "$BQSR_DIR"; cp "${FN[0]}"_BQSR.bam /lscratch/\$SLURM_JOB_ID/"${FN[0]}"_BQSR.bam && gatk --java-options \"-Xmx6G\" PrintReads -R "$CanFam4_Ref" -I /lscratch/\$SLURM_JOB_ID/"${FN[0]}"_BQSR.bam -O "$FILE_DIR""${changedname[0]}".BQSR.cram -OBM true -OBI false && samtools index -@ \$SLURM_CPUS_PER_TASK -c "$FILE_DIR""${changedname[0]}".BQSR.cram && cp "$TXT_DIR"/"${FN[0]}".sort.md.metrics.txt "$TXT_DIR"/"${changedname[0]}".sort.md.metrics.txt" >> "$homedir"/BQSR2CRAM_Rename.swarm
