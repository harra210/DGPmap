#!/bin/bash
#
# This is an internal script from the parent DGPmap.sh to gather bqsr recalibration tables.
#
## Variable Definition Section
FQ_DIR=$(pwd)
#
cd txt_files
TXT_DIR=$(pwd)
> bqsrreports.list
cd $FQ_DIR
#
cd BQSR/tables
TBL_DIR=$(pwd)
cd $FQ_FIR
#
PREFIX="--input "
#
##
# Capture a specific value to place into an array to complete the filename properly
FN=( $(awk 'NR==1{print$1; exit}' "$TXT_DIR"/final_header.txt) )
#
# Generate list of file names from chr 1 to chr X to later import into an array
#
while read g
do
	echo ""${FN[0]}"_"$g"_recal.table" >> "$TXT_DIR"/bqsrreports.list
done < "$interval_list"
#
# Now the same for chrY and the unmapped reads
echo ""${FN[0]}"_chrY_recal.table" >> "$TXT_DIR"/bqsrreports.list
#
# Import the now complete bqsrreports list file into an array for the swarm.
#
IFS=,$'\n' read -d '' -r -a samplereports < "$TXT_DIR"/bqsrreports.list
declare -a samplereports
#
# Now to fully generate the swarmfile
#
echo "cd "$TBL_DIR"; gatk --java-options \"-Xmx6G\" GatherBQSRReports "${samplereports[*]/#/$PREFIX}" --output "$TBL_DIR"/"${FN[0]}".bqsrgathered.reports.list" >> "$homedir"/bqsr_gatherBQSRReports.swarm
#
## This completes the gathering of BQSR Reports into one cohesive report
