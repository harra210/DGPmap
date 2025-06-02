#!/bin/bash
#
# This begins after the creation of sample tables and will use the table to parse out all of the
# swarmfiles to submit to the cluster. This script will also call additional scripts as needed and be more in line
# with the original Dog10K pipeline.
#
# Flag Processing Section to parse flags from the invokation of pipeline and export flags variables to their respective scripts.
#
# Currently only have one flag of use to be used in generating proper tables but this function will allow for modification later if different alignments are desired to be used. If that becomes the case, then that will necessitate the modification of subscripts.
#
# Variable defining defaults
FILE_IN=
OUT_DIR=
NAME=
SWARM_NAME=
YEAR=$(date +%Y)

#
# getopts string
opts=":mo:n:s:h"
#
# Gets the command name without path
cmd(){ echo `basename $0`; }
#
# Help command output
usage(){
        echo "\
        `cmd` [OPTION...]
                -m, --mapfile; Optional: Supplied tab-delimited file containing two columns [samplename] [location].
                -o, --output; Req'd: Base Output directory to place joint call files.
                -n, --name; Req'd: Name of the joint-called VCF
                -s, --swarm-name; Set the base swarm name for the pipeline, job specific names will be added on job submission.
                -h, --help; Print this message and exit.
        " | column -t -s ";"
}
#
# Error message
error(){
        echo "`cmd`: invalid option -- '$1'";
        echo "Try '`cmd` -h' or '`cmd`' --help for more information.";
        exit 1;
}
#
# Two passes for the flag parsing. The first pass handles the long options and any short option that is in canonical form.
# The second pass uses `getopt` to canonicalize any erroneous short options and handle them.
# This second pass will allow for future expansion of the script.
#
for pass in 1 2; do
        while [ -n "$1" ]; do
                case $1 in
                        --) shift; break;;
                        -*) case $1 in
                                -m|--mapfile)     export FILE_IN=$2; shift;;
                                -o|--output)     export OUT_DIR=$2; shift;;
                                -n|--name)      export NAME=$2; shift;;
                                -s|--swarm-name)        export SWARM_NAME=$2; shift;;
                                -h|--help) usage; exit 1;;
                                --*)    error $1;;
                                -*)     if [ $pass -eq 1 ]; then ARGS="$ARGS $1";
                                        else error $1; fi;;
                                esac;;
                        *) if [ $pass -eq 1 ]; then ARGS="$ARGS $1";
                                else error $1; fi;;
                esac
                shift
        done
        if [ $pass -eq 1 ]; then ARGS=`getopt $opts $ARGS`
                if [ $? != 0 ]; then usage; exit 2; fi; set -- $ARGS
        fi
done
#
# Handle positional arguments
if [ -n "$*" ]; then
        echo "`cmd`: Extra arguments -- $*"
        echo "Try '`cmd` -h' for more information."
        exit 1
fi
#
# Check to see if required flags are passed via command line
if [[ -z $OUT_DIR ]] && [[ -z $NAME ]] && [[ -z $SWARM_NAME ]]
then
        echo "ERROR: Missing required flags!"
        usage
        exit 1
elif [[ -z $OUT_DIR ]]
then
        echo "ERROR: Missing output directory!"
        echo "Try '`cmd` -h' for more information."
        exit 1
elif [[ -z $NAME ]]
then
        echo "ERROR: Missing flag for VCF name!"
        echo "Try '`cmd` -h' for more information."
        exit 1
elif [[ -z $SWARM_NAME ]]
then
        echo "ERROR: Missing swarm name flag!"
        echo "Try '`cmd` -h' for more information."
        exit 1
elif [[ -z $FILE_IN ]]
then
        while true; do
                read -p "WARNING: No mapfile was provided. Do you wish to continue? (Yes or No) " promptA # The reason for this prompt is a sanity check for the user as performing a full joint call is extremely resource-intensive and data heavy. The purpose for this script is so that users can input their own mapfile and create custom joint calls without the use of trimming.
                case $promptA in
                        [YyEeSs]* ) break;;
                        [NnOo]* ) echo "Try '`cmd` -h' for more information."; exit 1;;
                        * ) echo "Re-enter answer: Yes or No";;
                esac
        done
else
        shift
fi
# Handle if mapfile is set to null to set an additional variable to be used in generating the samplefile.
if [[ -z "${FILE_IN-x}" ]]
then
        GVCF_DIR="/data/Ostrander/CanFam4_GSD/gVCF/"
        export GVCF_DIR
else
        GVCF_DIR=""
        export GVCF_DIR
fi
#
## Variables to be set and/or exported to subshells
#
CF4_Ref="/data/Ostrander/Resources/CanFam4_GSD/BWAMEM2/UU_Cfam_GSD_1.0_ROSY.fa" # CanFam 4.0 GSD Reference. Tolerant to add references in this section later by adding an additional flag if there is need for it.
export CF4_Ref
#
CF4_VQSR="/data/Ostrander/Resources/CanFam4_GSD/BWAMEM2/SRZ189891_722g.simp.header.CanineHDandAxiom_K9_HD.GSD_1.0.vcf.gz" # CanFam 4.0 GSD VQSR Training Set
export CF4_VQSR
#
homedir=$(pwd)
export homedir
#
mkdir -p swarmfiles
cd swarmfiles
swarmdir=$(pwd)
> GenomicsDBImport_GenotypeGVCFs.swarm
> GenomicsDBImport_GenotypeGVCFs_XNonPAR.swarm
> GATK_Gathervcfs.swarm
> DGPjointcall-concatAutoXPAR.swarm
> DGPjointcall-concatXNonPAR.swarm
> VQSR_SelectVariants.swarm
> VQSR_VariantFiltration.swarm
> VQSR_VariantRecalibrator.swarm
export swarmdir
cd $homedir
#
cd ../tmp
tmpdir=$(pwd)
export tmpdir
#
mkdir -p bcftools_intervalmap # make this directory in the tmp directory if it does not exist. This folder sets the bcftools concatenation intervals per chromosome for parallelization.
cd bcftools_intervalmap
bcfintervalmap=$(pwd)
export bcfintervalmap
#
cd "$tmpdir"
mkdir -p gatk_chr_samplemap
cd gatk_chr_samplemap
gatkchrmap=$(pwd)
export gatkchrmap
#
cd "$homedir"
#
## Now we traverse to the desired output directory to generate all of the folders required to complete the joint call.
#
cd $OUT_DIR
mkdir -p RAW_AutoXPAR
cd RAW_AutoXPAR/
#
RAW_Auto=$(pwd)
#
cd $OUT_DIR
#
mkdir -p RAW_AutoXPAR_chr
cd RAW_AutoXPAR_chr
#
RAW_AutoCHR=$(pwd)
#
cd $OUT_DIR
#
mkdir -p RAW_XNonPAR
cd RAW_XNonPAR/
RAW_NONPAR=$(pwd)
cd $tmpdir
#
cd $OUT_DIR
#
mkdir -p RAW_VQSR
cd RAW_VQSR/
VQSR_DIR=$(pwd)
mkdir -p Final_VQSR
#
cd Final_VQSR
FILTER_DIR=$(pwd)
#
## Important BASH note for functions, within functions do not use the declare tag to formalize array variables as seen in previous scripts. Ex. declare -a sample ; what this does within the function is blank out the contents of the variable and will cause issues when referencing the variable in other areas of the script.
#
## Now we create an if/else statement to generate a samplemap file if the flag is null, and if one is supplied the script will pause to prompt the user if the file supplied is the correct file
#
mapmake(){
        cd "$tmpdir"
        > GenomicsDB_samplenames.txt
        > GenomicsDB_directories.txt
        > gvcf_final.txt
        cd "$GVCF_DIR"
        find . -maxdepth 1 -name '*.g.vcf.gz' -printf '%f\n' | sed 's/.g.vcf.gz//' | sort -n &> "$tmpdir"/GenomicsDB_samplenames.txt
        find $PWD -maxdepth 1 -name '*.g.vcf.gz' -printf '%h\n' &> "$tmpdir"/GenomicsDB_directories.txt
        cd "$tmpdir"
        IFS=,$'\n' read -d '' -r -a samplename < GenomicsDB_samplenames.txt
        for sample in "${samplename[@]}";
        do
                echo ""$GVCF_DIR""$sample".g.vcf.gz" >> gvcf_final.txt;
        done
        paste "$tmpdir"/GenomicsDB_samplenames.txt "$tmpdir"/gvcf_final.txt > "$swarmdir"/GenomicsDB_samplemap.txt
	tail "$swarmdir"/GenomicsDB_samplemap.txt
	read -p "This is the last 10 lines of that file. Is this the correct file and is the syntax correct? (Yes or No) " promptA
                while true; do
                        case "$promptA" in
                                [YyEeSs]* ) break ;;
                                [NnOo]* ) echo "Verify that the inputs are correct and try again"; exit;;
                                *) echo "Enter yes or no" ;;
                esac
        done
        
}
#
samplemap(){
	if [[ ! -z "$GVCF_DIR" ]]
	then
        	mapmake
	else
        	echo "You have selected "$FILE_IN" as your samplemap."
	        echo ""
        	tail "$FILE_IN"
	        echo ""
	        read -p "This is the last 10 lines of that file. Is this the correct file and is the syntax correct? (Yes or No) " promptA
        	while true; do
                	case "$promptA" in
                        	[YyEeSs]* ) cp "$FILE_IN" "$swarmdir"/GenomicsDB_samplemap.txt; break ;;
                        	[NnOo]* ) echo "Verify that the inputs are correct and try again"; exit;;
                        	*) echo "Enter yes or no" ;;
        	esac
	done
	fi
}
#
## Now that the sample maps are created now we can get to the meat of the joint call prep.
#
## At this point we are going to create all of the functions that will be used.
#
## First are creating the swarm files for the GenomicsDBImport and GenotypeGVCFs step for both the AutoXPAR and nonXPAR regions.
#
AutoXPAR_gDB(){
        cd "$swarmdir" # Make sure we are in the swarm directory so that we can properly redirect commands
#        > GenomicsDBImport_GenotypeGVCFs.swarm
        #
        while read i
        do
                echo "ulimit -u 16384 && gatk --java-options \"-Xmx7g -Xms7g\" GenomicsDBImport --tmp-dir /lscratch/\$SLURM_JOB_ID/ --L $i --sample-name-map GenomicsDB_samplemap.txt --batch-size 50 --genomicsdb-workspace-path /lscratch/\$SLURM_JOB_ID/"$i" --genomicsdb-shared-posixfs-optimizations && gatk --java-options \"-Xmx6g -Xms6g\" GenotypeGVCFs --tmp-dir /lscratch/\$SLURM_JOB_ID/ -R $CF4_Ref -O "$RAW_Auto"/"$NAME"."$i".AutoXPAR.RAW.vcf.gz -V gendb:///lscratch/\$SLURM_JOB_ID/"$i"" >> GenomicsDBImport_GenotypeGVCFs.swarm
        done < /data/Ostrander/Resources/CanFam4_GSD/Intervals/Full_intervals.intervals
	#
	cp GenomicsDB_samplemap.txt "$OUT_DIR"
}
#
XNonPAR_gDB(){
        cd "$swarmdir" # Make sure we are in the home directory. This is included even in this function needs to be called apart from the AutoXPAR function.
#        > GenomicsDBImport_GenotypeGVCFs_XNonPAR.swarm
        while read i
        do
                echo "ulimit -u 16384 && gatk --java-options \"-Xmx7g -Xms7g\" GenomicsDBImport --tmp-dir /lscratch/\$SLURM_JOB_ID/ --L $i --sample-name-map GenomicsDB_samplemap.txt --batch-size 50 --genomicsdb-workspace-path /lscratch/\$SLURM_JOB_ID/"$i" --genomicsdb-shared-posixfs-optimizations && gatk --java-options \"-Xmx6g -Xms6g\" GenotypeGVCFs --tmp-dir /lscratch/\$SLURM_JOB_ID/ -R $CF4_Ref -O "$RAW_NONPAR"/"$NAME"."$i".chrX.NONPAR.RAW.vcf.gz -V gendb:///lscratch/\$SLURM_JOB_ID/"$i"" >> GenomicsDBImport_GenotypeGVCFs_XNonPAR.swarm
        done < /data/Ostrander/Resources/CanFam4_GSD/Intervals/XNonPAR_intervals.intervals
}
#
#
## That completes the functions of generating the GenomicsDBImports and GenotypeGVCFs steps and now we have to concatenate all of the files together to make sure they all fit.
## First we will change into that directory, perform a test to make sure that the directory is empty of all files, and if not then remove all files.
#
## Check to see if the AutoXPAR interval map directory is empty. If not, delete all files within the directory so that interval files can be populated
## within the directory without duplicates.
#
bcfAutoXPAR_intervals(){
	cd "$bcfintervalmap"
        #
	for f in {1..38};
        do
                while read i
                do
                        echo ""$RAW_Auto"/"$NAME"."$i".AutoXPAR.RAW.vcf.gz" >> "$bcfintervalmap"/chr"$f"_AutoXPARsamples.txt
                done < /data/Ostrander/Resources/CanFam4_GSD/Intervals/chr/chr"$f"/chr"$f".intervals
        done
        #
        while read i
        do
                echo ""$RAW_Auto"/"$NAME"."$i".AutoXPAR.RAW.vcf.gz" >> "$bcfintervalmap"/chrX_AutoXPARsamples.txt
        done < /data/Ostrander/Resources/CanFam4_GSD/Intervals/chr/chrX/chrX.intervals
}
#
## Issue to figure out for this function is why when rerunning the echo repeats until complete rather than just one echo'ed line.
#
bcfAutoXPAR_rerun(){
        cd "$bcfintervalmap"
	#
	for file in "$bcfintervalmap"/*_AutoXPARsamples.txt
        do
                if [ -s "$file" ]
                then
                        echo "Blanking old AutoXPAR interval files..." && rm -R *_AutoXPARsamples.txt && printf "done" && bcfAutoXPAR_intervals
                else
                        bcfAutoXPAR_intervals
                fi
        done
}
#
## Now we do the same for chromosome X NONPAR regions
#
bcfXNonPAR_intervals(){
        cd "$bcfintervalmap"
	#
	while read i
        do
                echo ""$RAW_NONPAR"/"$NAME"."$i".chrX.NONPAR.RAW.vcf.gz" >> "$bcfintervalmap"/XNonPARsamples.txt
        done < /data/Ostrander/Resources/CanFam4_GSD/Intervals/XNonPAR_intervals.intervals
	#
        IFS=,$'\n' read -d '' -r -a NonPARvcf < "$bcfintervalmap"/XNonPARsamples.txt
        sortedNonPARvcf=( $(printf "%s\n" ${NonPARvcf[*]} | sort -V ) )
        #
        cd $swarmdir
        > DGPjointcall-XNonPARsamplemap.txt
        #
        printf "%s\n" "${sortedNonPARvcf[@]}" >> "$swarmdir"/DGPjointcall-XNonPARsamplemap.txt
}
#
bcfXNonPAR_rerun(){
        cd "$bcfintervalmap"
        #
	for file in "$bcfintervalmap"/XNonPARsamples.txt
        do
                if [ -e "$file" ]
                then
                        echo "Blanking old XNonPAR interval file..." && rm -R XNonPARsamples.txt && printf "done" && bcfXNonPAR_intervals
                else
                        bcfXNonPAR_intervals
                fi
        done
}
#
## Now, confirmatory checks need to be done to concat everything to chromosome level as this check is performed at the GenomicsDBImport steps. The AutoXPAR will have 2 commands written to the swarm file, chr 1-38 and then chrX.
#
bcfconcat_AutoXPAR(){
	cd "$swarmdir"
	#
	for i in {1..38};
	do
        echo "cd "$swarmdir"; bcftools concat -a -D --threads \$SLURM_CPUS_PER_TASK -f "$bcfintervalmap"/chr"$i"_AutoXPARsamples.txt -O z -o "$RAW_AutoCHR"/"$NAME".chr"$i".AutoXPAR.vcf.gz && gatk IndexFeatureFile -I "$RAW_AutoCHR"/"$NAME".chr"$i".AutoXPAR.vcf.gz --tmp-dir /lscratch/\$SLURM_JOB_ID" >> DGPjointcall-concatAutoXPAR.swarm
	done
#
	echo "cd $swarmdir; bcftools concat -a -D --threads \$SLURM_CPUS_PER_TASK -f "$bcfintervalmap"/chrX_AutoXPARsamples.txt -O z -o "$RAW_AutoCHR"/"$NAME".chrX.AutoXPAR.vcf.gz && gatk IndexFeatureFile -I "$RAW_AutoCHR"/"$NAME".chrX.AutoXPAR.vcf.gz --tmp-dir /lscratch/\$SLURM_JOB_ID" >> DGPjointcall-concatAutoXPAR.swarm
}
#
## Here is the bcfconcat for the XNonPAR areas.
#
bcfconcat_XNonPAR(){
	cd "$swarmdir"
	#
        echo "cd $swarmdir; bcftools concat -a -D --threads \$SLURM_CPUS_PER_TASK -f DGPjointcall-XNonPARsamplemap.txt -O z -o "$VQSR_DIR"/"$NAME".chrX.NONPAR.vcf.gz && gatk IndexFeatureFile -I "$VQSR_DIR"/"$NAME".chrX.NONPAR.vcf.gz --tmp-dir /lscratch/\$SLURM_JOB_ID" >> DGPjointcall-concatXNonPAR.swarm
}
#
GATK_gathervcfs(){
        cd "$gatkchrmap"
        #
        > chr_AutoXPARsamplemap.txt
        #
        for i in {1..38};
        do
                echo ""$NAME".chr"$i".AutoXPAR.vcf.gz" >> chr_AutoXPARsamplemap.txt
        done
        #
        echo ""$NAME".chrX.AutoXPAR.vcf.gz" >> chr_AutoXPARsamplemap.txt
        #
        cd "$swarmdir"
	#echo ""$gatkchrmap"/chr_AutoXPARsamplemap.txt" # This is a debug line to verify the contents of the text file.
	#sleep 30 # Debug line to pause the script; Ctrl+C to exit the script or wait 30 seconds to continue.
        #
        IFS=,$'\n' read -d '' -r -a vcf < "$gatkchrmap"/chr_AutoXPARsamplemap.txt
        #echo "${vcf[*]}" # This is a debug line to verify $vcf is not null and its contents is there.
	#sleep 30 # Debug line to pause the script; Ctrl+C to exit the script or wait 30 seconds to continue.
	export vcf
	sortedvcf=( $(printf "%s\n" ${vcf[*]} | sort -V ) )
        export sortedvcf
	#echo "${sortedvcf[*]}" # This is a debug line to verify $sortedvcf is not null and its contents is there.
#	sleep 30 # Debug line to pause the script; Ctrl+C to exit the script or wait 30 seconds to continue.
	#
        PREFIX="-I "
	export PREFIX
        #
        echo "cd "$RAW_AutoCHR"; gatk GatherVcfsCloud "${sortedvcf[*]/#/$PREFIX}" -O "$VQSR_DIR"/"$NAME".AutoXPAR.vcf.gz --tmp-dir /lscratch/\$SLURM_JOB_ID && gatk IndexFeatureFile -I "$VQSR_DIR"/"$NAME".AutoXPAR.vcf.gz --tmp-dir /lscratch/\$SLURM_JOB_ID" >> GATK_Gathervcfs.swarm
        #
}
#
## This section details functions for the AutoXPAR and NONPAR Variants selecting for SNPs
#
GATK_SelectVariants(){
        cd "$swarmdir"
        #
        echo "cd "$VQSR_DIR"; gatk --java-options \"-Xmx4g\" SelectVariants -V "$NAME".AutoXPAR.vcf.gz -O "$NAME".AutoXPAR.SNPs.vcf.gz -select-type SNP && gatk --java-options \"-Xmx4g\" SelectVariants -V "$NAME".AutoXPAR.vcf.gz -O "$NAME".AutoXPAR.nonSNPs.vcf.gz -xl-select-type SNP && gatk --java-options \"-Xmx4g\" SelectVariants -V "$NAME".chrX.NONPAR.vcf.gz -O "$NAME".chrX.NONPAR.SNPs.vcf.gz -select-type SNP && gatk --java-options \"-Xmx4g\" SelectVariants -V "$NAME".chrX.NONPAR.vcf.gz -O "$NAME".chrX.NONPAR.nonSNPs.vcf.gz -xl-select-type SNP" > VQSR_SelectVariants.swarm
        #
}
#
## This section details functions for GATK Variant Filtration function for Indels.
#
GATK_VariantFiltration(){
        cd "$swarmdir"
        #
        echo "cd "$VQSR_DIR"; gatk --java-options \"-Xmx4g\" VariantFiltration -V "$NAME".AutoXPAR.nonSNPs.vcf.gz -O "$FILTER_DIR"/"$NAME".AutoXPAR.nonSNPs.filtered.vcf.gz --verbosity ERROR -filter \"QD < 2.0\" --filter-name \"QD2\" -filter \"FS > 200.0\" --filter-name \"FS200\" -filter \"ReadPosRankSum < -2.0\" --filter-name \"ReadPosRankSum-2\" -filter \"SOR > 10.0\" --filter-name \"SOR-10\" && gatk --java-options \"-Xmx4g\" VariantFiltration -V "$NAME".chrX.NONPAR.nonSNPs.vcf.gz -O "$FILTER_DIR"/"$NAME".chrX.NONPAR.nonSNPs.filtered.vcf.gz --verbosity ERROR -filter \"QD < 2.0\" --filter-name \"QD2\" -filter \"FS > 200.0\" --filter-name \"FS200\" -filter \"ReadPosRankSum < -2.0\" --filter-name \"ReadPosRankSum-2\" -filter \"SOR > 10.0\" --filter-name \"SOR-10\"" > VQSR_VariantFiltration.swarm
        #
}
#
GATK_VariantRecalibrator(){
	cd "$swarmdir"
	#
	echo "cd "$VQSR_DIR"; gatk --java-options \"-Xmx59g\" VariantRecalibrator -R "$CF4_Ref" -V "$NAME".AutoXPAR.SNPs.vcf.gz -O "$NAME".AutoXPAR.SNPs.recal -resource:array,known=false,training=true,truth=true,prior=12.0 "$CF4_VQSR" --use-annotation QD --use-annotation MQ --use-annotation MQRankSum --use-annotation ReadPosRankSum --use-annotation FS --use-annotation SOR --use-annotation DP --trust-all-polymorphic true -mode SNP --rscript-file "$NAME".AutoXPAR.SNPs.plots.R --tranches-file "$NAME".AutoXPAR.SNPs.tranches -tranche 100.00 -tranche 99.9 -tranche 99.8 -tranche 99.7 -tranche 99.5 -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 92.5 -tranche 90.0 --tmp-dir /lscratch/\$SLURM_JOB_ID/ && gatk --java-options \"-Xmx59g\" ApplyVQSR -R "$CF4_Ref" -V "$NAME".AutoXPAR.SNPs.vcf.gz -O "$FILTER_DIR"/"$NAME".AutoXPAR.SNPs.VQSR.vcf.gz --tmp-dir /lscratch/\$SLURM_JOB_ID/ --truth-sensitivity-filter-level 99.0 --tranches-file "$NAME".AutoXPAR.SNPs.tranches --recal-file "$NAME".AutoXPAR.SNPs.recal -mode SNP && gatk --java-options \"-Xmx59g\" VariantRecalibrator -R "$CF4_Ref" -V "$NAME".chrX.NONPAR.SNPs.vcf.gz -O "$NAME".chrX.NONPAR.SNPs.recal -resource:array,known=false,training=true,truth=true,prior=12.0 "$CF4_VQSR" --use-annotation QD --use-annotation MQ --use-annotation MQRankSum --use-annotation ReadPosRankSum --use-annotation FS --use-annotation SOR --use-annotation DP --trust-all-polymorphic true -mode SNP --max-gaussians 4 --rscript-file "$NAME".X.output.plots.R --tranches-file "$NAME".X.SNP.output.tranches -tranche 100.00 -tranche 99.9 -tranche 99.5 -tranche 99.3 -tranche 99.0 -tranche 98.5 -tranche 98.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 92.5 -tranche 90.0 --tmp-dir /lscratch/\$SLURM_JOB_ID/ && gatk --java-options \"-Xmx59g\" ApplyVQSR -R "$CF4_Ref" -V "$NAME".chrX.NONPAR.SNPs.vcf.gz -O "$FILTER_DIR"/"$NAME".chrX.NONPAR.SNPs.VQSR.vcf.gz --tmp-dir /lscratch/\$SLURM_JOB_ID/ --truth-sensitivity-filter-level 99.0 --tranches-file "$NAME".X.SNP.output.tranches --recal-file "$NAME".chrX.NONPAR.SNPs.recal -mode SNP" > VQSR_VariantRecalibrator.swarm
	#
}
####
#### This section will contain functions that perform checks for each sections which will then be used to dictate how to submit the pipeline to the cluster
####
#
## First section is for the GenomicsDBImport of the AutoXPAR and NONPAR regions
#
AutoXPARfilecheck(){
        for file in "$RAW_Auto"/*.AutoXPAR.RAW.vcf.gz.tbi
        do
                if [ -e "$file" ]
                then
                        wcRAW=$(find "$RAW_Auto" -name "*.AutoXPAR.RAW.vcf.gz.tbi" | wc -l);
                else
                        wcRAW=0
                fi
		export wcRAW
        done
}
#
XNonPARfilecheck(){
        for file in "$RAW_NONPAR"/*.chrX.NONPAR.RAW.vcf.gz
        do
                if [ -e "$file" ]
                then

                        wcXNON=$(find "$RAW_NONPAR" -name "*.chrX.NONPAR.RAW.vcf.gz.tbi" | wc -l);
                else
                        wcXNON=0
                fi
		export wcXNON
        done
}
#
bcfconcat_AutoXPAR_check(){
        for file in "$RAW_AutoCHR"/*.AutoXPAR.vcf.gz.tbi
        do
                if [ -e "$file" ]
                then
                        bcAuto=$(find "$RAW_AutoCHR" -name "*.AutoXPAR.vcf.gz.tbi" | wc -l);
                else
                        bcAuto=0
                fi
		export bcAuto
        done
}
#
bcfconcat_XNonPAR_check(){
	for file in "$VQSR_DIR"/*.chrX.NONPAR.vcf.gz.tbi
	do
		if [ -e "$file" ]
		then
			bcXNPAR=$(find "$VQSR_DIR" -name "*.chrX.NONPAR.vcf.gz.tbi" | wc -l)
		else
			bcXNPAR=0
		fi
		export bcXNPAR
	done
}
#
GATK_gathervcfs_check(){
	for file in "$VQSR_DIR"/*.AutoXPAR.vcf.gz.tbi
	do
		if [ -e "$file" ]
		then
			gatkGV=$(find "$VQSR_DIR" -name "*.AutoXPAR.vcf.gz.tbi" | wc -l)
		else
			gatkGV=0
		fi
		export gatkGV
	done
}
#
GATK_SelectVariants_check(){
        for file in "$VQSR_DIR"/*.vcf.gz.tbi
        do
                if [ -e "$file" ]
                then
                        vcfSV=$(find "$VQSR_DIR" -name "*.vcf.gz.tbi" | wc -l)
                else
                        vcfSV=0
                fi
		export vcfSV
        done
}
#
GATK_VariantFiltration_check(){
        for file in "$VQSR_DIR"/*.nonSNPs.vcf.gz.tbi
        do
                if [ -e "$file" ]
                then
                        vcfVF=$(find "$VQSR_DIR" -name "*.nonSNPs.vcf.gz.tbi" | wc -l)
                else
                        vcfVF=0
                fi
		export vcfVF
	done
}
#
GATK_VariantRecalibrator_check(){
	for file in "$FILTER_DIR"/*.SNPs.VQSR.vcf.gz
	do
		if [ -e "$file" ]
		then
			vcfVR=$(find "$FILTER_DIR" -name "*.SNPs.VQSR.vcf.gz.tbi" | wc -l)
		else
			vcfVR=0
		fi
		export vcfVR
	done
}
#
####
#### This section contains all of the retry functions required to make the pipeline work
####
#
AutoXPARtruncatedretry(){
        cd "$RAW_Auto"
        find . -name "*.AutoXPAR.RAW.vcf.gz.tbi" -printf "%f\n" | sed 's/.tbi//' | awk -F '.' '{print $2}' &> "$tmpdir"/DGPjointcall-truncatedretry-tbishard.tmp
        #
        ## Now output the difference between the two files and output that difference to a new file
        cd "$tmpdir"
        comm -3 <(sort DGPjointcall-truncatedretry-tbishard.tmp) <(sort /data/Ostrander/Resources/CanFam4_GSD/Intervals/Full_intervals.intervals) | tr -d '\t' > "$tmpdir"/truncatedinterval_retry.tmp
}
#
## If AutoXPARtruncatedretry is invoked, then invoke the following function to recreate the swarmfile with just the truncated shards to resubmit as GDBI.
#
AutoXPARtruncatedGdbI(){
        cd "$swarmdir" # Make sure we are in the swarm directory so that we can properly redirect commands
        > GenomicsDBImport_GenotypeGVCFs.swarm
        #
        while read i
        do
                echo "ulimit -u 16384 && gatk --java-options \"-Xmx7g -Xms7g\" GenomicsDBImport --tmp-dir /lscratch/\$SLURM_JOB_ID/ --L $i --sample-name-map GenomicsDB_samplemap.txt --batch-size 50 --genomicsdb-workspace-path /lscratch/\$SLURM_JOB_ID/"$i" --genomicsdb-shared-posixfs-optimizations && gatk --java-options \"-Xmx6g -Xms6g\" GenotypeGVCFs --tmp-dir /lscratch/\$SLURM_JOB_ID/ -R $CF4_Ref -O "$RAW_Auto"/"$NAME"."$i".AutoXPAR.RAW.vcf.gz -V gendb:///lscratch/\$SLURM_JOB_ID/"$i"" >> GenomicsDBImport_GenotypeGVCFs.swarm
        done < "$tmpdir"/truncatedinterval_retry.tmp
}
#
XNonPARtruncatedretry(){
        cd $RAW_NONPAR
        find . -name "*.chrX.NONPAR.RAW.vcf.gz.tbi" -printf "%f\n" | sed 's/.tbi//' | awk -F '.' '{print $2}' &> "$tmpdir"/DGPjointcall-truncatedretry-NONPARtbishard.tmp
        #
        cd "$tmpdir"
        comm -3 <(sort DGPjointcall-truncatedretry-NONPARtbishard.tmp) <(sort /data/Ostrander/Resources/CanFam4_GSD/Intervals/XNonPAR_intervals.intervals) | tr -d '\t' > "$tmpdir"/NONPARtruncatedinterval_retry.tmp
}
#
XNonPARtruncatedGdbI(){
        cd "$swarmdir" # Make sure we are in the home directory. This is included even in this function needs to be called apart from the AutoXPAR function.
        > GenomicsDBImport_GenotypeGVCFs_XNonPAR.swarm
        while read i
        do
                echo "ulimit -u 16384 && gatk --java-options \"-Xmx7g -Xms7g\" GenomicsDBImport --tmp-dir /lscratch/\$SLURM_JOB_ID/ --L $i --sample-name-map GenomicsDB_samplemap.txt --batch-size 50 --genomicsdb-workspace-path /lscratch/\$SLURM_JOB_ID/"$i" --genomicsdb-shared-posixfs-optimizations && gatk --java-options \"-Xmx6g -Xms6g\" GenotypeGVCFs --tmp-dir /lscratch/\$SLURM_JOB_ID/ -R $CF4_Ref -O "$RAW_NONPAR"/"$NAME"."$i".chrX.NONPAR.RAW.vcf.gz -V gendb:///lscratch/\$SLURM_JOB_ID/"$i"" >> GenomicsDBImport_GenotypeGVCFs_XNonPAR.swarm
        done < "$tmpdir"/NONPARtruncatedinterval_retry.tmp
}
#
####
#### This section now that we've generated variable counts that are global, we'll utilize them to call the other functions and generate a usable pipeline that can be modal
####
#
## First we call all of the functions to get values for each variable check
#
AutoXPARfilecheck
XNonPARfilecheck
bcfconcat_AutoXPAR_check
bcfconcat_XNonPAR_check
GATK_gathervcfs_check
GATK_SelectVariants_check
GATK_VariantFiltration_check
GATK_VariantRecalibrator_check
#
##########
## DEBUG AREA ##
##########
#
# Here we just echo the number the files considered as each variable. These variable numbers are numbers to correleate with the conditional branching.
#
#echo "wcRAW $wcRAW"
#echo "wcXNON $wcXNON"
#echo "bcAuto $bcAuto"
#echo "bcXNPAR $bcXNPAR"
#echo "gatkGV $gatkGV"
#echo "vcfSV $vcfSV"
#echo "vcfVF $vcfVF"
#echo "vcfVR $vcfVR"

#sleep 30 ## AS OF 5/21/25 THIS PART IS VALIDATED.
#########
## END DEBUG AREA ##
#########
#
## Here we create all of the functions with all of the possibilities, later we'll create the conditional branching
#
fullsubmit(){
	samplemap
	AutoXPAR_gDB
	XNonPAR_gDB
	bcfAutoXPAR_intervals
	bcfXNonPAR_intervals
	bcfconcat_AutoXPAR
	bcfconcat_XNonPAR
	GATK_gathervcfs
	GATK_SelectVariants
	GATK_VariantFiltration
	GATK_VariantRecalibrator
}
#
full_rerun(){
	samplemap
	AutoXPARtruncatedretry
	AutoXPARtruncatedGdbI
	XNonPARtruncatedretry
	XNonPARtruncatedGdbI
	bcfAutoXPAR_rerun
        bcfXNonPAR_rerun
        bcfconcat_AutoXPAR
        bcfconcat_XNonPAR
        GATK_gathervcfs
        GATK_SelectVariants
        GATK_VariantFiltration
        GATK_VariantRecalibrator
}
#
AutoXPARGdbI_rerun(){
	samplemap
	AutoXPARtruncatedretry
	AutoXPARtruncatedGdbI
	bcfAutoXPAR_rerun
	bcfXNonPAR_rerun
	bcfconcat_AutoXPAR
	bcfconcat_XNonPAR
	GATK_gathervcfs
	GATK_SelectVariants
	GATK_VariantFiltration
	GATK_VariantRecalibrator
}
#
XNonPARGdbI_rerun(){
	samplemap
	XNonPARtruncatedretry
	XNonPARtruncatedGdbI
	bcfAutoXPAR_rerun
	bcfXNonPAR_rerun
	bcfconcat_AutoXPAR
	bcfconcat_XNonPAR
	GATK_gathervcfs
	GATK_SelectVariants
	GATK_VariantFiltration
	GATK_VariantRecalibrator
}
#
bcfconcatFULL_rerun(){
	bcfconcat_AutoXPAR
	bcfconcat_XNonPAR
	GATK_gathervcfs
	GATK_SelectVariants
	GATK_VariantFiltration
	GATK_VariantRecalibrator
}
#
bcfconcatAutoXPAR_rerun(){
	bcfconcat_AutoXPAR
	GATK_gathervcfs
	GATK_SelectVariants
	GATK_VariantFiltration
	GATK_VariantRecalibrator
}
#
bcfconcatXNonPAR_rerun(){
	bcfconcat_XNonPAR
	GATK_gathervcfs
	GATK_SelectVariants
	GATK_VariantFiltration
	GATK_VariantRecalibrator
}
#
GATK_gathervcfs_rerun(){
	GATK_gathervcfs
        GATK_SelectVariants
        GATK_VariantFiltration
        GATK_VariantRecalibrator
}
GATK_SelectVariants_rerun(){
	GATK_SelectVariants
        GATK_VariantFiltration
        GATK_VariantRecalibrator
}
GATK_VariantFiltration_rerun(){
	GATK_VariantFiltration
        GATK_VariantRecalibrator
}
GATK_VariantRecalibrator_rerun(){
	GATK_VariantRecalibrator
}
#
## This section contains the conditional branching to populate the required swarm files. This mainly helps for when there are node failures that cause conditional swarm failures. As the swarm is submitted to chained only if the previous job exits with a exit code of 0, if one subjob fails the entire pipeline will stop.
#
if [ $wcRAW == 0 ] && [ $wcXNON == 0 ] # No AutoXPAR or X NonPAR vcf shards = full pipeline submit
then
	fullsubmit
elif [ $wcRAW -ge 0 ] && ! [ $wcRAW == 2255 ] && [ $wcXNON -ge 0 ] && ! [ $wcXNON == 119 ] #If AutoXPAR AND XNonPAR shards are missing, will perform a full pipeline rerun.
then
	full_rerun
elif  [ $wcRAW -ge 0 ] && ! [ $wcRAW == 2255 ]  && [ $wcXNON == 119 ] #If missing AutoXPAR shards, will only run AutoXPAR GdbI. XPAR will not run.
then
	AutoXPARGdbI_rerun
elif [ $wcRAW == 2255 ] && [ $wcXNON -ge 0 ] && ! [ $wcXNON == 119 ] #If AutoXPAR shards are complete but missing XNonPAR shards detected, will rerun XNonPAR.
then
	XNonPARGdbI_rerun
elif [ $wcRAW == 2255 ] && [ $wcXNON == 119 ] && ! [ $bcAuto == 39 ] #GdbI completed successfully but missing 39 concatenated vcf files then BCFtools concat will rerun
then
	bcfconcatAutoXPAR_rerun
elif [ $wcRAW == 2255 ] && [ $wcXNON == 119 ] && [ $bcAuto == 39 ] && [ $bcXNPAR == 0 ] && [ $gatkGV == 1 ] #GdbI completed successfully and AutoXPAR concat completed successfully but XNonPAR failed and needs redoing. Note, can add an additional branch since gathervcfs is separate to xnonpar.
then
	bcfconcatXNonPAR_rerun
elif [ $wcRAW == 2255 ] && [ $wcXNON == 119 ] && [ $bcAuto == 39 ] && [ $bcXNPAR == 1 ] && ! [ $gatkGV == 1 ] #If AutoXPAR GATK Gathervcfs failed and needs to be rerun from here.
then
	GATK_gathervcfs_rerun
elif [ $wcRAW == 2255 ] && [ $wcXNON == 119 ] && [ $bcAuto == 39 ] && [ $bcXNPAR == 1 ] && [ $gatkGV == 1 ] && ! [ $vcfSV == 6 ] #This checks to see if there are all of the proper files in the VQSR directory
then
	GATK_SelectVariants_rerun
elif [ $wcRAW == 2255 ] && [ $wcXNON == 119 ] && [ $bcAuto == 39 ] && [ $bcXNPAR == 1 ] && [ $gatkGV == 1 ] && [ $vcfSV == 6 ] && ! [ $vcfVF == 2 ] # Checks if nonSNPs are good or need to be rerun
then
	GATK_VariantFiltration_rerun
elif [ $wcRAW == 2255 ] && [ $wcXNON == 119 ] && [ $bcAuto == 39 ] && [ $bcXNPAR == 1 ] && [ $gatkGV == 1 ] && [ $vcfSV == 6 ] && [ $vcfVF == 2 ] && ! [ $vcfVR == 2 ] #Checks to make sure that all of the VQSR vcfs are ready.
then
	GATK_VariantRecalibrator_rerun
else
	echo "Error detected! Troubleshooting required!"
fi
#
## Here we create the functions to submit swarm files to the cluster
#
fullsubmit_swarm(){
	local jobid1=$(swarm -f GenomicsDBImport_GenotypeGVCFs.swarm -g 36 -p 2 --time 36:00:00 --module GATK/4.6.0.0 --gres=lscratch:100 --logdir /data/"$USER"/job_outputs/"$YEAR"/GenomicsDBImport/"$SWARM_NAME" --sbatch "--mail-type=ALL,TIME_LIMIT_80 --job-name "$SWARM_NAME"_GDBI --requeue")
	echo "GenomicsDBImport AutoXPAR Swarm ID: "$jobid1""
	#
	local jobid2=$(swarm -f GenomicsDBImport_GenotypeGVCFs_XNonPAR.swarm -g 36 -p 2 --time 24:00:00 --module GATK/4.6.0.0 --gres=lscratch:100 --logdir /data/"$USER"/job_outputs/"$YEAR"/GenomicsDBImport/"$SWARM_NAME"_XNonPAR --sbatch "--mail-type=ALL,TIME_LIMIT_80 --job-name "$SWARM_NAME"_GDBI_XNonPAR --requeue")
	echo "GenomicsDBImport XNonPAR Swarm ID: "$jobid2""
	#
	local jobid3=$(swarm -f DGPjointcall-concatAutoXPAR.swarm -g 12 -t 12 --time 2-0 --gres=lscratch:150 --module bcftools/1.16,GATK/4.6.0.0 --logdir ~/job_outputs/BCFtools/Concat/"$SWARM_NAME"_AutoXPAR --sbatch "--mail-type=ALL,TIME_LIMIT_80 --job-name "$SWARM_NAME"_AutoXPAR_Concat --dependency=afterok:"$jobid1","$jobid2"")
	echo "BCFTools Chr Concat Swarm ID: "$jobid3""
	#
	local jobid4=$(swarm -f DGPjointcall-concatXNonPAR.swarm -g 32 -t 8 --time 24:00:00 --gres=lscratch:150 --module bcftools/1.16,GATK/4.6.0.0 --logdir ~/job_outputs/BCFtools/Concat/"$SWARM_NAME"_XNonPAR --sbatch "--mail-type=ALL,TIME_LIMIT_80 --job-name "$SWARM_NAME"_XNonPAR_Concat --dependency=afterok:"$jobid2"")
	echo "BCFTools XNonPar Concat Swarm ID: "$jobid4""
	#
	local jobid5=$(swarm -f GATK_Gathervcfs.swarm -g 5 -t 4 --time 16:00:00 --gres=lscratch:125 --module GATK/4.6.0.0 --logdir ~/job_outputs/gatk/gathervcfs/"$SWARM_NAME"_gathervcfs --sbatch "--mail-type=ALL,TIME_LIMIT_80 --job-name "$SWARM_NAME"_gathervcfs --dependency=afterok:"$jobid3"")
	echo "Gathervcfs Swarm ID: "$jobid5""
	#
	local jobid6=$(swarm -f VQSR_SelectVariants.swarm -g 8 -t 10 --time 24:00:00 --module GATK/4.6.0.0 --logdir ~/job_outputs/gatk/SelectVariants/"$SWARM_NAME" --sbatch "--mail-type=ALL,TIME_LIMIT_80 --job-name "$SWARM_NAME"_SelectVariants --dependency=afterok:"$jobid5","$jobid4"")
	echo "SelectVariants Swarm ID: "$jobid6""
	#
	local jobid7=$(swarm -f VQSR_VariantFiltration.swarm -g 8 -t 10 --time 12:00:00 --module GATK/4.6.0.0 --logdir ~/job_outputs/gatk/VariantFiltration/"$SWARM_NAME" --sbatch "--mail-type=ALL,TIME_LIMIT_80 --job-name "$SWARM_NAME"_VariantFiltration --dependency=afterok:"$jobid6"")
	echo "VariantFiltration Swarm ID: "$jobid7""
	#
	local jobid8=$(swarm -f VQSR_VariantRecalibrator.swarm -g 72 -t 12 --time 3-0 --gres=lscratch:250 --module GATK/4.6.0.0 --logdir ~/job_outputs/gatk/VariantRecalibrator/"$SWARM_NAME" --sbatch "--mail-type=ALL,TIME_LIMIT_80 --job-name "$SWARM_NAME"_VariantRecal --dependency=afterok:"$jobid7"")
	echo "VariantRecalibrator Swarm ID: "$jobid8""
}	# Fullsubmit_swarm is to be used for also the full retry submit
#
AutoXPAR_rerun_swarm(){
	local jobid1=$(swarm -f GenomicsDBImport_GenotypeGVCFs.swarm -g 36 -p 2 --time 36:00:00 --module GATK/4.6.0.0 --gres=lscratch:100 --logdir /data/"$USER"/job_outputs/"$YEAR"/GenomicsDBImport/"$SWARM_NAME" --sbatch "--mail-type=ALL,TIME_LIMIT_80 --job-name "$SWARM_NAME"_GDBI --requeue")
        echo "GenomicsDBImport AutoXPAR Swarm ID: "$jobid1""
        #
	local jobid2=$(swarm -f DGPjointcall-concatAutoXPAR.swarm -g 12 -t 12 --time 2-0 --gres=lscratch:150 --module bcftools/1.16,GATK/4.6.0.0 --logdir ~/job_outputs/BCFtools/Concat/"$SWARM_NAME"_AutoXPAR --sbatch "--mail-type=ALL,TIME_LIMIT_80 --job-name "$SWARM_NAME"_AutoXPAR_Concat --dependency=afterok:"$jobid1"")
        echo "BCFTools Chr Concat Swarm ID: "$jobid2""
	#
	local jobid3=$(swarm -f GATK_Gathervcfs.swarm -g 5 -t 4 --time 16:00:00 --gres=lscratch:125 --module GATK/4.6.0.0 --logdir ~/job_outputs/gatk/gathervcfs/"$SWARM_NAME"_gathervcfs --sbatch "--mail-type=ALL,TIME_LIMIT_80 --job-name "$SWARM_NAME"_gathervcfs --dependency=afterok:"$jobid2"")
        echo "Gathervcfs Swarm ID: "$jobid3""
	#
	local jobid4=$(swarm -f VQSR_SelectVariants.swarm -g 8 -t 10 --time 24:00:00 --module GATK/4.6.0.0 --logdir ~/job_outputs/gatk/SelectVariants/"$SWARM_NAME" --sbatch "--mail-type=ALL,TIME_LIMIT_80 --job-name "$SWARM_NAME"_SelectVariants --dependency=afterok:"$jobid3"")
        echo "SelectVariants Swarm ID: "$jobid4""
	#
	local jobid5=$(swarm -f VQSR_VariantFiltration.swarm -g 8 -t 10 --time 12:00:00 --module GATK/4.6.0.0 --logdir ~/job_outputs/gatk/VariantFiltration/"$SWARM_NAME" --sbatch "--mail-type=ALL,TIME_LIMIT_80 --job-name "$SWARM_NAME"_VariantFiltration --dependency=afterok:"$jobid4"")
        echo "VariantFiltration Swarm ID: "$jobid5""
	#
	local jobid6=$(swarm -f VQSR_VariantRecalibrator.swarm -g 72 -t 12 --time 3-0 --gres=lscratch:250 --module GATK/4.6.0.0 --logdir ~/job_outputs/gatk/VariantRecalibrator/"$SWARM_NAME" --sbatch "--mail-type=ALL,TIME_LIMIT_80 --job-name "$SWARM_NAME"_VariantRecal --dependency=afterok:"$jobid5"")
        echo "VariantRecalibrator Swarm ID: "$jobid6""
}
#
XNonPAR_rerun_swarm(){
	local jobid1=$(swarm -f GenomicsDBImport_GenotypeGVCFs_XNonPAR.swarm -g 36 -p 2 --time 24:00:00 --module GATK/4.6.0.0 --gres=lscratch:100 --logdir /data/"$USER"/job_outputs/"$YEAR"/GenomicsDBImport/"$SWARM_NAME"_XNonPAR --sbatch "--mail-type=ALL,TIME_LIMIT_80 --job-name "$SWARM_NAME"_GDBI_XNonPAR --requeue")
        echo "GenomicsDBImport XNonPAR Swarm ID: "$jobid1""
	#
	local jobid2=$(swarm -f DGPjointcall-concatXNonPAR.swarm -g 32 -t 8 --time 24:00:00 --gres=lscratch:150 --module bcftools/1.16,GATK/4.6.0.0 --logdir ~/job_outputs/BCFtools/Concat/"$SWARM_NAME"_XNonPAR --sbatch "--mail-type=ALL,TIME_LIMIT_80 --job-name "$SWARM_NAME"_XNonPAR_Concat --dependency=afterok:"$jobid1"")
        echo "BCFTools XNonPar Concat Swarm ID: "$jobid2""
	#
	local jobid3=$(swarm -f GATK_Gathervcfs.swarm -g 5 -t 4 --time 16:00:00 --gres=lscratch:125 --module GATK/4.6.0.0 --logdir ~/job_outputs/gatk/gathervcfs/"$SWARM_NAME"_gathervcfs --sbatch "--mail-type=ALL,TIME_LIMIT_80 --job-name "$SWARM_NAME"_gathervcfs --dependency=afterok:"$jobid2"")
        echo "Gathervcfs Swarm ID: "$jobid3""
        #
        local jobid4=$(swarm -f VQSR_SelectVariants.swarm -g 8 -t 10 --time 24:00:00 --module GATK/4.6.0.0 --logdir ~/job_outputs/gatk/SelectVariants/"$SWARM_NAME" --sbatch "--mail-type=ALL,TIME_LIMIT_80 --job-name "$SWARM_NAME"_SelectVariants --dependency=afterok:"$jobid3"")
        echo "SelectVariants Swarm ID: "$jobid4""
        #
        local jobid5=$(swarm -f VQSR_VariantFiltration.swarm -g 8 -t 10 --time 12:00:00 --module GATK/4.6.0.0 --logdir ~/job_outputs/gatk/VariantFiltration/"$SWARM_NAME" --sbatch "--mail-type=ALL,TIME_LIMIT_80 --job-name "$SWARM_NAME"_VariantFiltration --dependency=afterok:"$jobid4"")
        echo "VariantFiltration Swarm ID: "$jobid5""
        #
        local jobid6=$(swarm -f VQSR_VariantRecalibrator.swarm -g 72 -t 12 --time 3-0 --gres=lscratch:250 --module GATK/4.6.0.0 --logdir ~/job_outputs/gatk/VariantRecalibrator/"$SWARM_NAME" --sbatch "--mail-type=ALL,TIME_LIMIT_80 --job-name "$SWARM_NAME"_VariantRecal --dependency=afterok:"$jobid5"")
        echo "VariantRecalibrator Swarm ID: "$jobid6""
}
#
bcfconcat_rerun_swarm(){
	local jobid1=$(swarm -f DGPjointcall-concatAutoXPAR.swarm -g 12 -t 12 --time 2-0 --gres=lscratch:150 --module bcftools/1.16,GATK/4.6.0.0 --logdir ~/job_outputs/BCFtools/Concat/"$SWARM_NAME"_AutoXPAR --sbatch "--mail-type=ALL,TIME_LIMIT_80 --job-name "$SWARM_NAME"_AutoXPAR_Concat")
        echo "BCFTools Chr Concat Swarm ID: "$jobid1""
	#
	local jobid2=$(swarm -f DGPjointcall-concatXNonPAR.swarm -g 32 -t 8 --time	24:00:00 --gres=lscratch:150 --module bcftools/1.16,GATK/4.6.0.0 --logdir ~/job_outputs/BCFtools/Concat/"$SWARM_NAME"_XNonPAR --sbatch "--mail-type=ALL,TIME_LIMIT_80 --job-name "$SWARM_NAME"_XNonPAR_Concat")
        echo "BCFTools XNonPar Concat Swarm ID: "$jobid2""
	#
	local jobid3=$(swarm -f GATK_Gathervcfs.swarm -g 5 -t 4 --time 16:00:00 --gres=lscratch:125 --module GATK/4.6.0.0 --logdir ~/job_outputs/gatk/gathervcfs/"$SWARM_NAME"_gathervcfs --sbatch "--mail-type=ALL,TIME_LIMIT_80 --job-name "$SWARM_NAME"_gathervcfs --dependency=afterok:"$jobid1"")
        echo "Gathervcfs Swarm ID: "$jobid3""
        #
        local jobid4=$(swarm -f VQSR_SelectVariants.swarm -g 8 -t 10 --time 24:00:00 --module GATK/4.6.0.0 --logdir ~/job_outputs/gatk/SelectVariants/"$SWARM_NAME" --sbatch "--mail-type=ALL,TIME_LIMIT_80 --job-name "$SWARM_NAME"_SelectVariants --dependency=afterok:"$jobid2","$jobid3"")
        echo "SelectVariants Swarm ID: "$jobid4""
        #
        local jobid5=$(swarm -f VQSR_VariantFiltration.swarm -g 8 -t 10 --time 12:00:00 --module GATK/4.6.0.0 --logdir ~/job_outputs/gatk/VariantFiltration/"$SWARM_NAME" --sbatch "--mail-type=ALL,TIME_LIMIT_80 --job-name "$SWARM_NAME"_VariantFiltration --dependency=afterok:"$jobid4"")
        echo "VariantFiltration Swarm ID: "$jobid5""
        #
        local jobid6=$(swarm -f VQSR_VariantRecalibrator.swarm -g 72 -t 12 --time 3-0 --gres=lscratch:250 --module GATK/4.6.0.0 --logdir ~/job_outputs/gatk/VariantRecalibrator/"$SWARM_NAME" --sbatch "--mail-type=ALL,TIME_LIMIT_80 --job-name "$SWARM_NAME"_VariantRecal --dependency=afterok:"$jobid5"")
        echo "VariantRecalibrator Swarm ID: "$jobid6""
}
#
bcfconcatXNonPAR_rerun_swarm(){
	local jobid2=$(swarm -f DGPjointcall-concatXNonPAR.swarm -g 32 -t 8 --time 24:00:00 --gres=lscratch:150 --module bcftools/1.16,GATK/4.6.0.0 --logdir ~/job_outputs/BCFtools/Concat/"$SWARM_NAME"_XNonPAR --sbatch "--mail-type=ALL,TIME_LIMIT_80 --job-name "$SWARM_NAME"_XNonPAR_Concat")
        echo "BCFTools XNonPar Concat Swarm ID: "$jobid2""
        #
        local jobid4=$(swarm -f VQSR_SelectVariants.swarm -g 8 -t 10 --time 24:00:00 --module GATK/4.6.0.0 --logdir ~/job_outputs/gatk/SelectVariants/"$SWARM_NAME" --sbatch "--mail-type=ALL,TIME_LIMIT_80 --job-name "$SWARM_NAME"_SelectVariants --dependency=afterok:"$jobid2"")
        echo "SelectVariants Swarm ID: "$jobid4""
        #
        local jobid5=$(swarm -f VQSR_VariantFiltration.swarm -g 8 -t 10 --time 12:00:00 --module GATK/4.6.0.0 --logdir ~/job_outputs/gatk/VariantFiltration/"$SWARM_NAME" --sbatch "--mail-type=ALL,TIME_LIMIT_80 --job-name "$SWARM_NAME"_VariantFiltration --dependency=afterok:"$jobid4"")
        echo "VariantFiltration Swarm ID: "$jobid5""
        #
        local jobid6=$(swarm -f VQSR_VariantRecalibrator.swarm -g 72 -t 12 --time 3-0 --gres=lscratch:250 --module GATK/4.6.0.0 --logdir ~/job_outputs/gatk/VariantRecalibrator/"$SWARM_NAME" --sbatch "--mail-type=ALL,TIME_LIMIT_80 --job-name "$SWARM_NAME"_VariantRecal --dependency=afterok:"$jobid5"")
        echo "VariantRecalibrator Swarm ID: "$jobid6""
}
#
GATKgathervcfs_rerun_swarm(){
	local jobid1=$(swarm -f GATK_Gathervcfs.swarm -g 5 -t 4 --time 16:00:00 --gres=lscratch:125 --module GATK/4.6.0.0 --logdir ~/job_outputs/gatk/gathervcfs/"$SWARM_NAME"_gathervcfs --sbatch "--mail-type=ALL,TIME_LIMIT_80 --job-name "$SWARM_NAME"_gathervcfs")
        echo "Gathervcfs Swarm ID: "$jobid1""
        #
        local jobid2=$(swarm -f VQSR_SelectVariants.swarm -g 8 -t 10 --time 24:00:00 --module GATK/4.6.0.0 --logdir ~/job_outputs/gatk/SelectVariants/"$SWARM_NAME" --sbatch "--mail-type=ALL,TIME_LIMIT_80 --job-name "$SWARM_NAME"_SelectVariants --dependency=afterok:"$jobid1"")
        echo "SelectVariants Swarm ID: "$jobid2""
        #
        local jobid3=$(swarm -f VQSR_VariantFiltration.swarm -g 8 -t 10 --time 12:00:00 --module GATK/4.6.0.0 --logdir ~/job_outputs/gatk/VariantFiltration/"$SWARM_NAME" --sbatch "--mail-type=ALL,TIME_LIMIT_80 --job-name "$SWARM_NAME"_VariantFiltration --dependency=afterok:"$jobid2"")
        echo "VariantFiltration Swarm ID: "$jobid3""
        #
        local jobid4=$(swarm -f VQSR_VariantRecalibrator.swarm -g 72 -t 12 --time 3-0 --gres=lscratch:250 --module GATK/4.6.0.0 --logdir ~/job_outputs/gatk/VariantRecalibrator/"$SWARM_NAME" --sbatch "--mail-type=ALL,TIME_LIMIT_80 --job-name "$SWARM_NAME"_VariantRecal --dependency=afterok:"$jobid3"")
        echo "VariantRecalibrator Swarm ID: "$jobid4""
}
#
GATKselectvariants_rerun_swarm(){
	local jobid1=$(swarm -f VQSR_SelectVariants.swarm -g 8 -t 10 --time 24:00:00 --module GATK/4.6.0.0 --logdir ~/job_outputs/gatk/SelectVariants/"$SWARM_NAME" --sbatch "--mail-type=ALL,TIME_LIMIT_80 --job-name "$SWARM_NAME"_SelectVariants")
        echo "SelectVariants Swarm ID: "$jobid1""
        #
        local jobid2=$(swarm -f VQSR_VariantFiltration.swarm -g 8 -t 10 --time 12:00:00 --module GATK/4.6.0.0 --logdir ~/job_outputs/gatk/VariantFiltration/"$SWARM_NAME" --sbatch "--mail-type=ALL,TIME_LIMIT_80 --job-name "$SWARM_NAME"_VariantFiltration --dependency=afterok:"$jobid1"")
        echo "VariantFiltration Swarm ID: "$jobid2""
        #
        local jobid3=$(swarm -f VQSR_VariantRecalibrator.swarm -g 72 -t 12 --time 3-0 --gres=lscratch:250 --module GATK/4.6.0.0 --logdir ~/job_outputs/gatk/VariantRecalibrator/"$SWARM_NAME" --sbatch "--mail-type=ALL,TIME_LIMIT_80 --job-name "$SWARM_NAME"_VariantRecal --dependency=afterok:"$jobid2"")
        echo "VariantRecalibrator Swarm ID: "$jobid3""
}
#
GATKvariantfiltration_rerun_swarm(){
	local jobid1=$(swarm -f VQSR_VariantFiltration.swarm -g 8 -t 10 --time 12:00:00 --module GATK/4.6.0.0 --logdir ~/job_outputs/gatk/VariantFiltration/"$SWARM_NAME" --sbatch "--mail-type=ALL,TIME_LIMIT_80 --job-name "$SWARM_NAME"_VariantFiltration")
        echo "VariantFiltration Swarm ID: "$jobid1""
        #
        local jobid2=$(swarm -f VQSR_VariantRecalibrator.swarm -g 72 -t 12 --time 3-0 --gres=lscratch:250 --module GATK/4.6.0.0 --logdir ~/job_outputs/gatk/VariantRecalibrator/"$SWARM_NAME" --sbatch "--mail-type=ALL,TIME_LIMIT_80 --job-name "$SWARM_NAME"_VariantRecal --dependency=afterok:"$jobid1"")
        echo "VariantRecalibrator Swarm ID: "$jobid2""
}
#
GATKvariantrecal_rerun_swarm(){
	local jobid1=$(swarm -f VQSR_VariantRecalibrator.swarm -g 72 -t 12 --time 3-0 --gres=lscratch:250 --module GATK/4.6.0.0 --logdir ~/job_outputs/gatk/VariantRecalibrator/"$SWARM_NAME" --sbatch "--mail-type=ALL,TIME_LIMIT_80 --job-name "$SWARM_NAME"_VariantRecal")
        echo "VariantRecalibrator Swarm ID: "$jobid1""
}
#
## Here we gather all of the modification times of the swarmfiles now that the script has completed all of the swarm generation stages since epoch.
#
#cd "$swarmdir" # This will ensure that the script is now being executed in the swarm directory
#
#current_time=$(date +%s) # This gets the current time in seconds since epoch
#AutoGdbI_mod=$((current_time - $(stat -c %Y GenomicsDBImport_GenotypeGVCFs.swarm)))
#XNonPAR_mod=$((current_time - $(stat -c %Y GenomicsDBImport_GenotypeGVCFs_XNonPAR.swarm)))
#bcfAutoconcat_mod=$((current_time - $(stat -c %Y bcftools_concatAutoXPAR.swarm)))
#bcfXNonPARconcat_mod=$((current_time - $(stat -c %Y DGPjointcall-concatXNonPAR.swarm)))
#gathervcfs_mod=$((current_time - $(stat -c %Y GATK_Gathervcfs.swarm)))
#selectvariants_mod=$((current_time - $(stat -c %Y VQSR_SelectVariants.swarm)))
#variantfiltration_mod=$((current_time - $(stat -c %Y VQSR_VariantFiltration.swarm)))
#variantrecal_mod=$((current_time - $(stat -c %Y VQSR_VariantRecalibrator.swarm)))
#
## Now we create an inverse conditional branching tree for the submission based if swarmfiles are empty (0 bytes). Needs to start at VariantRecal and work its way up to the GDBI section as the script will read from top down.
#
#
## First is to verify swarm syntax before submission.
#
cd $swarmdir
if [ -s GATK_Gathervcfs.swarm ] && { [ -s GenomicsDBImport_GenotypeGVCFs.swarm ] || [ -s GenomicsDBImport_GenotypeGVCFs_XNonPAR.swarm ]; }
then
	if [ -s GenomicsDBImport_GenotypeGVCFs.swarm ] && [ -s GenomicsDBImport_GenotypeGVCFs_XNonPAR.swarm ]
	then
		echo "GenomicsDBImport GenotypeGVCFs AutoXPAR swarm:"
		head -n 1 GenomicsDBImport_GenotypeGVCFs.swarm
		echo ""
		echo "GenomicsDBImport GenotypeGVCFs XNonPAR swarm:"
		head -n 1 GenomicsDBImport_GenotypeGVCFs_XNonPAR.swarm
		echo ""
		echo "BCFtools AutoXPAR concat swarm:"
		head -n 1 DGPjointcall-concatAutoXPAR.swarm
		echo ""
		echo "BCFtools XNonPAR concat swarm:"
		head -n 1 DGPjointcall-concatXNonPAR.swarm
		echo ""
		echo "GATK Gathervcfs swarm:"
		head -n 1 GATK_Gathervcfs.swarm
		echo ""
		echo "GATK SelectVariants swarm:"
		head -n 1 VQSR_SelectVariants.swarm
		echo ""
		echo "GATK VariantFiltration.swarm"
		head -n 1 VQSR_VariantFiltration.swarm
		echo ""
		echo "GATK VariantRecalibrator swarm"
		head -n 1 VQSR_VariantRecalibrator.swarm
		echo ""
	elif [ -s GenomicsDBImport_GenotypeGVCFs.swarm ] && [ ! -s GenomicsDBImport_GenotypeGVCFs_XNonPAR.swarm ]
	then
		echo "GenomicsDBImport GenotypeGVCFs AutoXPAR swarm:"
                head -n 1 GenomicsDBImport_GenotypeGVCFs.swarm
                echo ""
		echo "BCFtools AutoXPAR concat swarm:"
                head -n 1 DGPjointcall-concatAutoXPAR.swarm
                echo ""
		echo "GATK Gathervcfs swarm:"
                head -n 1 GATK_Gathervcfs.swarm
                echo ""
                echo "GATK SelectVariants swarm:"
                head -n 1 VQSR_SelectVariants.swarm
                echo ""
                echo "GATK VariantFiltration.swarm"
                head -n 1 VQSR_VariantFiltration.swarm
                echo ""
                echo "GATK VariantRecalibrator swarm"
                head -n 1 VQSR_VariantRecalibrator.swarm
                echo ""
	else
		echo "GenomicsDBImport GenotypeGVCFs XNonPAR swarm:"
                head -n 1 GenomicsDBImport_GenotypeGVCFs_XNonPAR.swarm
                echo ""
		echo "BCFtools XNonPAR concat swarm:"
                head -n 1 DGPjointcall-concatXNonPAR.swarm
                echo ""
		echo "GATK Gathervcfs swarm:"
                head -n 1 GATK_Gathervcfs.swarm
                echo ""
                echo "GATK SelectVariants swarm:"
                head -n 1 VQSR_SelectVariants.swarm
                echo ""
                echo "GATK VariantFiltration.swarm"
                head -n 1 VQSR_VariantFiltration.swarm
                echo ""
                echo "GATK VariantRecalibrator swarm"
                head -n 1 VQSR_VariantRecalibrator.swarm
                echo ""
	fi
elif [ -s GATK_Gathervcfs.swarm ] && { [ -s DGPjointcall-concatAutoXPAR.swarm ] || [ -s DGPjointcall-concatXNonPAR.swarm ]; }
then
	if [ -s DGPjointcall-concatAutoXPAR.swarm ] && [ -s DGPjointcall-concatXNonPAR.swarm ]
	then
		echo "BCFtools AutoXPAR concat swarm:"
		head -n 1 DGPjointcall-concatAutoXPAR.swarm
		echo ""
        	echo "BCFtools XNonPAR concat swarm:"
		head -n 1 DGPjointcall-concatXNonPAR.swarm
        	echo ""
        	echo "GATK Gathervcfs swarm:"
        	head -n 1 GATK_Gathervcfs.swarm
        	echo ""
		echo "GATK SelectVariants swarm:"
		head -n 1 VQSR_SelectVariants.swarm
		echo ""
		echo "GATK VariantFiltration.swarm"
        	head -n 1 VQSR_VariantFiltration.swarm
        	echo ""
        	echo "GATK VariantRecalibrator swarm"
        	head -n 1 VQSR_VariantRecalibrator.swarm
        	echo ""
	elif [ -s DGPjointcall-concatAutoXPAR.swarm ] && [ ! -s DGPjointcall-concatXNonPAR.swarm ]
	then
		echo "BCFtools AutoXPAR concat swarm:"
                head -n 1 DGPjointcall-concatAutoXPAR.swarm
                echo ""
                echo "GATK Gathervcfs swarm:"
                head -n 1 GATK_Gathervcfs.swarm
                echo ""
                echo "GATK SelectVariants swarm:"
                head -n 1 VQSR_SelectVariants.swarm
                echo ""
                echo "GATK VariantFiltration.swarm"
                head -n 1 VQSR_VariantFiltration.swarm
                echo ""
                echo "GATK VariantRecalibrator swarm"
                head -n 1 VQSR_VariantRecalibrator.swarm
                echo ""
	else
                echo "BCFtools XNonPAR concat swarm:"
                head -n 1 DGPjointcall-concatXNonPAR.swarm
                echo ""
                echo "GATK Gathervcfs swarm:"
                head -n 1 GATK_Gathervcfs.swarm
                echo ""
                echo "GATK SelectVariants swarm:"
                head -n 1 VQSR_SelectVariants.swarm
                echo ""
                echo "GATK VariantFiltration.swarm"
                head -n 1 VQSR_VariantFiltration.swarm
                echo ""
                echo "GATK VariantRecalibrator swarm"
                head -n 1 VQSR_VariantRecalibrator.swarm
                echo ""
	fi
elif [ ! -s GATK_Gathervcfs.swarm ] && [ -s VQSR_SelectVariants.swarm ] && [ -s VQSR_VariantFiltration.swarm ] && [ -s VQSR_VariantRecalibrator.swarm ]
then
	echo "GATK SelectVariants swarm:"
        head -n 1 VQSR_SelectVariants.swarm
        echo ""
        echo "GATK VariantFiltration.swarm"
        head -n 1 VQSR_VariantFiltration.swarm
        echo ""
        echo "GATK VariantRecalibrator swarm"
        head -n 1 VQSR_VariantRecalibrator.swarm
        echo ""
elif [ ! -s VQSR_SelectVariants.swarm ] && [ -s VQSR_VariantFiltration.swarm ] && [ -s VQSR_VariantRecalibrator.swarm ]
then
	echo "GATK VariantFiltration.swarm"
        head -n 1 VQSR_VariantFiltration.swarm
        echo ""
        echo "GATK VariantRecalibrator swarm"
        head -n 1 VQSR_VariantRecalibrator.swarm
        echo ""
elif [ ! -s VQSR_VariantFiltration.swarm ] && [ -s VQSR_VariantRecalibrator.swarm ]
then
	echo "GATK VariantRecalibrator swarm"
        head -n 1 VQSR_VariantRecalibrator.swarm
        echo ""
else
	echo "Error! No swarmfiles detected to submit!"; exit 1
fi
##
while true; do
read -p "Do the swarmfiles above have proper syntax? (Yes or No) " promptA
case $promptA in
        [YyEeSs]* ) break;;
        [NnOo]* ) echo "Troubleshoot scripts"; exit 1;;
        * ) echo "Re-enter answer: Yes or No";;
esac
done
##
if [ -s GATK_Gathervcfs.swarm ] && { [ -s GenomicsDBImport_GenotypeGVCFs.swarm ] || [ -s GenomicsDBImport_GenotypeGVCFs_XNonPAR.swarm ]; }
then
	if [ -s GenomicsDBImport_GenotypeGVCFs.swarm ] && [ -s GenomicsDBImport_GenotypeGVCFs_XNonPAR.swarm ]
	then
		fullsubmit_swarm
	elif [ -s GenomicsDBImport_GenotypeGVCFs.swarm ] && [ ! -s GenomicsDBImport_GenotypeGVCFs_XNonPAR.swarm ]
	then
		AutoXPAR_rerun_swarm
	else
		XNonPAR_rerun_swarm
	fi
elif [ -s GATK_Gathervcfs.swarm ] && { [ -s DGPjointcall-concatAutoXPAR.swarm ] || [ -s DGPjointcall-concatXNonPAR.swarm ]; }
then
	if [ -s DGPjointcall-concatAutoXPAR.swarm ] && [ -s DGPjointcall-concatXNonPAR.swarm ]
	then
		bcfconcat_rerun_swarm
	elif [ -s DGPjointcall-concatAutoXPAR.swarm ] && [ ! -s DGPjointcall-concatXNonPAR.swarm ]
	then
		bcfconcat_rerun_swarm
	else
		bcfconcatXNonPAR_rerun_swarm
	fi
elif [ ! -s DGPjointcall-concatAutoXPAR.swarm ] && [ ! -s DGPjointcall-concatXNonPAR.swarm ] && [ -s GATK_Gathervcfs.swarm ] && [ -s VQSR_SelectVariants.swarm ] && [ -s VQSR_VariantFiltration.swarm ] && [ -s VQSR_VariantRecalibrator.swarm ]
then
	GATK_gathervcfs_rerun_swarm
elif [ ! -s GATK_Gathervcfs.swarm ] && [ -s VQSR_SelectVariants.swarm ] && [ -s VQSR_VariantFiltration.swarm ] && [ -s VQSR_VariantRecalibrator.swarm ]
then
	GATKselectvariants_rerun_swarm
elif [ -s VQSR_VariantFiltration.swarm ] && [ -s VQSR_VariantRecalibrator.swarm ]
then
	GATKvariantfiltration_rerun_swarm
elif [ ! -s VQSR_VariantFiltration.swarm ] && [ -s VQSR_VariantRecalibrator.swarm ]
then
	GATKvariantrecal_rerun_swarm
else
        echo "Error! No swarmfiles detected to submit!"; exit 1
fi
