#!/bin/bash
#
# This is an internal script for DGPmap to gather bam files previously created per chromosome
#
FQ_DIR=$(pwd)
#
## Variable Definition Section
#
cd txt_files
TXT_DIR=$(pwd)
> bqsr_gatherfiles.tmp
cd $FQ_DIR
#
cd BQSR/BQSR_chr
BQSRchr_DIR=$(pwd)
cd $FQ_DIR
#
cd BQSR
BQSR_DIR=$(pwd)
cd $FQ_DIR
#
PREFIX="-I "
#
## For NISC-style FastQ inputs
#
if [ "$FQ_IN" = "NISC" ];
	then
	FN=( $(awk 'NR==1{print$1; exit}' "$TXT_DIR"/final_header.txt ) )
#
# Generate BQSR sample name list into a tmp file in the text files folder
#
	while read g
		do
			echo ""${FN[0]}"."$g".BQSR.bam" >> "$TXT_DIR"/bqsr_gatherfiles.tmp
	done < "$interval_list"
#
# Now for chr Y and the unmapped reads
	echo ""${FN[0]}".chrY.BQSR.bam" >> "$TXT_DIR"/bqsr_gatherfiles.tmp
#
# Import the chromosomal bam names into an array
#
	IFS=,$'\n' read -d '' -r -a bqsrbamfiles < "$TXT_DIR"/bqsr_gatherfiles.tmp
	declare -a bqsrbamfiles
#
# Now to fully generate the swarmfile
#
	echo "cd "$BQSRchr_DIR"; gatk --java-options \"-Xmx6G\" GatherBamFiles --CREATE_INDEX true "${bqsrbamfiles[*]/#/$PREFIX}" -O "$BQSR_DIR"/"${FN[0]}".BQSR.bam" >> "$homedir"/bqsr_gatherBQSRBams.swarm
#
## For SRA-style FastQ inputs
#
elif [ "$FQ_IN" = "SRA" ];
	then
 	FN=( $(awk 'NR==1{print$1; exit}' "$FQ_DIR"/*.runs.table ) )
  #
  	while read g
		do
			echo ""${FN[0]}"."$g".BQSR.bam" >> "$TXT_DIR"/bqsr_gatherfiles.tmp
	done < "$interval_list"
#
# Now for chr Y and the unmapped reads
	echo ""${FN[0]}".chrY.BQSR.bam" >> "$TXT_DIR"/bqsr_gatherfiles.tmp
#
# Import the chromosomal bam names into an array
#
	IFS=,$'\n' read -d '' -r -a bqsrbamfiles < "$TXT_DIR"/bqsr_gatherfiles.tmp
	declare -a bqsrbamfiles
#
# Now to fully generate the swarmfile
#
	echo "cd "$BQSRchr_DIR"; gatk --java-options \"-Xmx6G\" GatherBamFiles --CREATE_INDEX true "${bqsrbamfiles[*]/#/$PREFIX}" -O "$BQSR_DIR"/"${FN[0]}".BQSR.bam" >> "$homedir"/bqsr_gatherBQSRBams.swarm
else
	exit 1
fi
