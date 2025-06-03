#!/bin/bash
#
# This is an internal script for DGPmap to generate gVCFs per chromosome
#
## Variable Definition Section
#
FQ_DIR=$(pwd)
#
cd txt_files
TXT_DIR=$(pwd)
cd $FQ_DIR
#
cd BQSR
BQSR_DIR=$(pwd)
cd $FQ_DIR
#
cd gVCF/HC
HC_DIR=$(pwd)
cd $FQ_DIR
#
## End Variable Definitions
#
## For Illumina-style inputs
if [ "$FQ_IN" = "Illumina" ];
	then
# Generate the sample name via awk from the original sample table
	FN=( $(awk 'NR==1{print $1; exit}' "$TXT_DIR"/final_header.txt ) )
#
# Now to generate the swarm file. For chr1 - chrX use a while loop and the for the remainder a simple echo line.
#
	while read g
		do
		echo "cd "$BQSR_DIR"; gatk --java-options \"-Xmx4G\" HaplotypeCaller -R "$CanFam4_Ref" --tmp-dir /lscratch/\$SLURM_JOB_ID -I "${FN[0]}".BQSR.bam -L "$g" -O "$HC_DIR"/"${FN[0]}"_"$g".g.vcf.gz -ERC GVCF" >> "$homedir"/haplotypecaller.swarm
	done < "$interval_list"
#
	echo "cd "$BQSR_DIR"; gatk --java-options \"-Xmx4G\" HaplotypeCaller -R "$CanFam4_Ref" --tmp-dir /lscratch/\$SLURM_JOB_ID -I "${FN[0]}".BQSR.bam -XL "$interval_list" -O "$HC_DIR"/"${FN[0]}"_chrY.g.vcf.gz -ERC GVCF" >> "$homedir"/haplotypecaller.swarm
#
## For SRA-style inputs
elif [ "$FQ_IN" = "SRA" ];
	then
 	FN=( $(awk 'NR==1{print$1; exit}' *runs.table ) )
#
# Generate list of file names from chr 1 to chr X to later import into an array
#
	while read g
		do
echo "cd "$BQSR_DIR"; gatk --java-options \"-Xmx4G\" HaplotypeCaller -R "$CanFam4_Ref" --tmp-dir /lscratch/\$SLURM_JOB_ID -I "${FN[0]}".BQSR.bam -L "$g" -O "$HC_DIR"/"${FN[0]}"_"$g".g.vcf.gz -ERC GVCF" >> "$homedir"/haplotypecaller.swarm
	done < "$interval_list"
#
	echo "cd "$BQSR_DIR"; gatk --java-options \"-Xmx4G\" HaplotypeCaller -R "$CanFam4_Ref" --tmp-dir /lscratch/\$SLURM_JOB_ID -I "${FN[0]}".BQSR.bam -XL "$interval_list" -O "$HC_DIR"/"${FN[0]}"_chrY.g.vcf.gz -ERC GVCF" >> "$homedir"/haplotypecaller.swarm
#
else
	exit 1
fi
