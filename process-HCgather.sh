#!/bin/bash
#
# This is an internal script for DGPmap to gather haplotypecaller shards into a singular gVCF. This is the last script before performing internal name ID changes and conversions to CRAM and removal of extraneous files.
#
## Variable Definition Section
FQ_DIR=$(pwd)
#
cd txt_files
TXT_DIR=$(pwd)
> HC_samplenames.tmp
cd $FQ_DIR
#
cd gVCF/HC
HC_DIR=$(pwd)
cd $FQ_DIR
#
cd gVCF
GVCF_DIR=$(pwd)
cd $FQ_DIR
#
PREFIX="-I "
#
## End Variable Definitions
#
## For NISC-style inputs
if [ "$FQ_IN" = "NISC" ];
	then
# Capture samplename from table using awk
	FN=( $(awk 'NR==1{print$1; exit}' "$TXT_DIR"/final_header.txt ) )
#
# Generate list of file names from chr1 to chrX using a while loop to later import into an array
	while read g
		do
			echo ""${FN[0]}"_"$g".g.vcf.gz" >> "$TXT_DIR"/HC_samplenames.tmp
		done < "$interval_list"
#
# Now echo a line for chrY and the unmapped reads
	echo ""${FN[0]}"_chrY.g.vcf.gz" >> "$TXT_DIR"/HC_samplenames.tmp
#
# Import the completed samplename file into an array
#
	IFS=,$'\n' read -d '' -r -a samplename < "$TXT_DIR"/HC_samplenames.tmp
	declare -a samplename # Redundancy measure to ensure that variable samplename is declared an array
#
# Now to generate the swarmfile
#
	echo "cd "$HC_DIR"; gatk --java-options \"-Xmx6G\" GatherVcfs "${samplename[*]/#/$PREFIX}" -O "$GVCF_DIR"/"${FN[0]}"_g.vcf.gz && tabix -p vcf -f "$GVCF_DIR"/"${FN[0]}"_g.vcf.gz" >> "$homedir"/hc_gathergvcfs.swarm
#
## For SRA-style input
elif [ "$FQ_IN" = "SRA" ];
	then
 	FN=( $(awk 'NR==1{print$1; exit}' "$FQ_DIR"/*.runs.table ) )
#
# Generate list of file names from chr1 to chrX using a while loop to later import into an array
	while read g
		do
			echo ""${FN[0]}"_"$g".g.vcf.gz" >> "$TXT_DIR"/HC_samplenames.tmp
		done < "$interval_list"
#
# Now echo a line for chrY and the unmapped reads
	echo ""${FN[0]}"_chrY.g.vcf.gz" >> "$TXT_DIR"/HC_samplenames.tmp
#
# Import the completed samplename file into an array
#
	IFS=,$'\n' read -d '' -r -a samplename < "$TXT_DIR"/HC_samplenames.tmp
	declare -a samplename # Redundancy measure to ensure that variable samplename is declared an array
#
# Now to generate the swarmfile
#
	echo "cd "$HC_DIR"; gatk --java-options \"-Xmx6G\" GatherVcfs "${samplename[*]/#/$PREFIX}" -O "$GVCF_DIR"/"${FN[0]}"_g.vcf.gz && tabix -p vcf -f "$GVCF_DIR"/"${FN[0]}"_g.vcf.gz" >> "$homedir"/hc_gathergvcfs.swarm
#
else
	exit 1
fi
#
## This completes the main computation portion of the Pipeline.
