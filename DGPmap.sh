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
# Fastq defaults
FQ_IN="Illumina"
FILE_DIR=
SWARM_NAME=
#
# getopts string
opts="i:f:s:h"
#
# Gets the command name without path
cmd(){ echo `basename $0`; }
#
# Help command output
usage(){
	echo "\
		`cmd` [OPTION...]
	-i, --input; Parent input directory of the files initially looking to process.
	-f, --fastq; Set input fastq style format (default:"$FQ_IN").
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
				-i|--input)	export FILE_DIR=$2; shift;;
				-f|--fastq)	export FQ_IN=$2; shift;;
				-s|--swarm-name)	export SWARM_NAME=$2; shift;;
				-h|--help) usage; exit 1;;
				--*)	error $1;;
				-*)	if [ $pass -eq 1 ]; then ARGS="$ARGS $1";
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
#echo "$FQ_IN"
#
# Variable Definition Sections that are able to be exported to subscripts. Note there are additional variables in subscripts but those are unable to be exported due to dynamics.
#
CanFam4_Ref="/data/Ostrander/Resources/CanFam4_GSD/BWAMEM2/UU_Cfam_GSD_1.0_ROSY.fa"
export CanFam4_Ref
#
knownsite="/data/Ostrander/Resources/CanFam4_GSD/BWAMEM2/UU_Cfam_GSD_1.0.BQSR.DB.bed"
export knownsite
#
interval_list="/data/Ostrander/Resources/CanFam4_GSD/Intervals/Intervals.list"
export interval_list
#
homedir=$(pwd)
export homedir
#
mkdir -p ../tmp # Make tmp directory to place generated temp files if it doesn't already exist.
cd ../tmp
tmpdir=$(pwd)
export tmpdir
#
cd $homedir
#
# Here, we blank all of the swarmfiles to be used from alignment to completion
> bwa_to_picard.swarm
> mergeDedup.swarm
> bqsr_BaseRecalibrator.swarm
> bqsr_gatherBQSRReports.swarm
> bqsr_ApplyBQSR.swarm
> bqsr_gatherBQSRBams.swarm
> haplotypecaller.swarm
> hc_gathergvcfs.swarm
#
#
cd scripts/
scriptdir=$(pwd)
#
cd $tmpdir
> DatasetDirectories.tmp
#
cd $FILE_DIR
#
find $PWD -mindepth 1 -maxdepth 1 -type d &> "$tmpdir"/DatasetDirectories.tmp
#
cd $tmpdir
#
IFS=,$'\n' read -d '' basedir < DatasetDirectories.tmp
#
cd $FILE_DIR # Change into the parent directory
#
## Begin the loop section to work on generating swarmfiles. The idea here is to create each section into a function so that the script will go in to the directories and check to see if certain files exist and if so, skip the creation of the swarm file. If it is only half done then it *should* only create those that are needed.
# 
#
tables(){ # This function is to generate the initial sample tables for the rest of the pipeline
		for dir in ${basedir[@]};
	do
		 cd "$dir";
		(  for file in "$dir"/*.runs.table
	 	do
		if ! [ -e "$file" ]
		then
			echo "Building table for "$dir""; bash "$scriptdir"/table-generation.sh
		else	
			echo "Skipping "$dir""; exit 0 
		fi 
	done )
	done
}
#
aln(){
	for dir in ${basedir[@]};
	do
		cd "$dir";
		( for file in "$dir"/bwa_bam/sort_*.bam
		do
			wcf=$(find "$dir" -name "*.fastq.gz" | wc -l);
			#echo $wcf
				if [ -e "$file" ]
				then
					wcb=$(find "$dir"/bwa_bam/ -name "sort_*.bam" | wc -l); # echo "$wcb" for debug
				else
					wcb=0; # echo $wcb for debug
				fi
			if ! [ $((wcf/2)) == "$wcb" ]
			then
				bash "$scriptdir"/process-BWA.sh; echo "BWA swarm made for ""$dir""";
			else
				echo "Skipping "$dir""; exit 0
			fi
		done )
	done
}
#
# This part builds the samtools and GATK MarkDuplicates swarmfile
#
markdup(){
        for dir in ${basedir[@]};
        do
                cd "$dir";
                ( for file in "$dir"/dedup_bam/*.sort.md.bam
                do
			wcs=$(find "$dir"/bwa_bam/ -name "*.bam" | wc -l);
			wcf=$(find "$dir" -name "*fastq.gz" | wc -l);
			#echo $wcs
                                if [ -e "$file" ]
                                then
					wcd=$(find "$dir"/dedup_bam/ -name "*.sort.md.bam" | wc -l); # echo $wcd for debug
                                else
                                        wcd=0; # echo $wcd for debug
                                fi
				if [ -e $file ] && [ $((wcf/2)) == "$wcs" ]
                        then
                                echo "Skipping "$dir""; exit 0
			else
                                bash "$scriptdir"/process-sortmd.sh; echo "Markduplicate swarm made for ""$dir""";
                        fi
                done )
        done
}
#
# Now build the BQSR Recal table swarmfile
#
bqsrtables(){
	for dir in ${basedir[@]};
	do
		cd "$dir";
		( for file in "$dir"/dedup_bam/*.sort.md.bam
			do
				wcr=$(find "$dir"/BQSR/tables -name "*_recal.table" | wc -l);
				#echo $wcr # debug line
					if [ -e $file ] && [ $wcr == 40 ]
					then
						echo "Skipping ""$dir"""; exit 0
					else
						bash "$scriptdir"/process-BQSR.sh; echo "BQSR swarm made for ""$dir""";
					fi
		done )
	done
}
#
# Gather BQSR Reports
#
gatherbqsr(){
	for dir in ${basedir[@]};
	do
		cd "$dir";
		( for file in "$dir"/BQSR/tables/*.bqsrgathered.reports.list
			do
				wcr=$(find "$dir"/BQSR/tables -name "*_recal.table" | wc -l);
				#echo $wcr # debug line
				if [ -e $file ] && [ $wcr == 40 ]
				then
					echo "Skipping ""$dir"""; exit 0
				else
					echo "Gathering BQSR tables for ""$dir"""; bash "$scriptdir"/process-BQSRreportsGather.sh
				fi
		done )
	done
}
# Apply BaseRecalibration by chromosome
#
applybqsr(){
	for dir in ${basedir[@]};
	do
		cd "$dir";
		( for file in "$dir"/BQSR/tables/*.bqsrgathered.reports.list
			do
				wcb=$(find "$dir"/BQSR/BQSR_chr/ -name "*.BQSR.bam" | wc -l);
				#echo $wcb # debug line
				if [ -e $file ] && [ $wcb == 40 ]
				then
					echo "Skipping ""$dir"""; exit 0
				else
					echo "Creating ApplyBQSR swarm for ""$dir"""; bash "$scriptdir"/process-applybqsr.sh
				fi
		done )
	done
}
#
# Gather BQSR bam files into a singular BQSR bam file
#
bqsrbamgather(){
	for dir in ${basedir[@]};
	do
		cd "$dir";
		( for file in "$dir"/BQSR/*.BQSR.bam
			do
				wcb=$(find "$dir"/BQSR/BQSR_chr/ -name "*.BQSR.bam" | wc -l);
				#echo $wcb # debug line
				if [ -e $file ] && [ $wcb == 40 ]
				then
					echo "Skipping ""$dir"""; exit 0
				else
					echo "Creating BQSR bam consolidation swarm for ""$dir"""; bash "$scriptdir"/process-gatherBQSRbams.sh
				fi
		done )
	done
}
#
# Perform HaplotypeCaller per chromosome
#
haplotypecaller(){
	for dir in ${basedir[@]};
	do
		cd "$dir";
		( for file in "$dir"/BQSR/*.BQSR.bam
			do
				wch=$(find "$dir"/gVCF/HC/ -name "*.g.vcf.gz" | wc -l);
				#echo $wch # debug line
				if [ -e $file ] && [ $wch == 40 ]
				then
					echo "Skipping ""$dir"""; exit 0
				else
					echo "Creating swarm for ""$dir""'s haplotypecaller step"; bash "$scriptdir"/process-haplotypecaller.sh
				fi
		done )
	done
}
#
# Gather HaplotypeCaller gVCFs into a singular gVCF file
#
hcgather(){
	for dir in ${basedir[@]};
	do
		cd "$dir";
		( for file in "$dir"/HC/*_g.vcf.gz
			do
				wch=$(find "$dir"/gVCF/HC/ -name "*.g.vcf.gz" | wc -l);
				#echo $wch # debug line
				if [ -e $file ] && [ $wch == 40 ]
				then
					echo "Skipping ""$dir"""; exit 0
				else
					echo "Creating gVCF consolidation swarm for ""$dir"""; bash "$scriptdir"/process-HCgather.sh
				fi
		done )
	done
}
#
## Functions have now been created up until the haplotypecaller step. Here the pipeline can be forked from here to the next step of conversion from BAM to CRAM, recompression of gVCFs, cleanup and stats generation
#
## Now we call the aforementioned functions to create the swarmfiles
#
## First step is to create tables.
echo "Generating sample tables"
sleep 1
echo ""
tables || { echo "Error generating tables"; exit 0; }
echo ""
#
## Create BWA swarmfile.
echo "Generating BWAMEM swarmfile"
sleep 1
echo ""
aln || { echo "Error generating BWA swarmfile"; exit 0; }
echo ""
#
## Create markduplicates swarmfile.
echo "Generating merging and marking duplicates swarmfile"
sleep 1
echo ""
markdup || { echo "Error generating Markduplicates swarmfile"; exit 0; }
echo ""
#
## Perform BQSR
echo "Generating BQSR swarmfile"
sleep 1
echo ""
bqsrtables || { echo "Error generating BQSR swarmfile"; exit 0; }
echo ""
#
## Gather BQSR tables
echo "Generating swarmfile to gather BQSR tables"
sleep 1
echo ""
gatherbqsr || { echo "Error generating GatherBQSR swarmfile"; exit 0; }
echo ""
#
## ApplyBQSR using gathered BQSR tables
#
echo "Generating ApplyBQSR swarmfile"
sleep 1
echo ""
applybqsr || { echo "Error generating ApplyBQSR swarmfile"; exit 0; }
echo ""
#
## GatherBQSRbam files
#
echo "Generating swarmfile to gather BQSR bam files"
sleep 1
echo ""
bqsrbamgather || { echo "Error generating swarmfile to gather BQSR bam files"; exit 0; }
echo ""
#
## Generate gVCF files via HaplotypeCaller
#
echo "Generating swarmfile to create gVCF files"
sleep 1
echo ""
haplotypecaller || { echo "Error generating HaplotypeCaller swarm"; exit 0; }
echo ""
#
## Consolidate gVCF files into a singular gVCF file
#
echo "Generating swarmfile to consolidate gVCF files"
sleep 1
echo ""
hcgather || { echo "Error generating gathergVCF file swarm"; exit 0; }
echo ""
#
## Now all swarmfiles are hopefully generated successfull and now the functions to define whether to submit fully to the cluster or partially submit files to the cluster.
#
pipelinesubmit() {
jobid1=$(swarm -f bwa_to_picard.swarm -g 72 -t 32 --gres=lscratch:300 --time 2-0 --module bwa-mem2,samtools --logdir ~/job_outputs/bwa_to_picard/"$SWARM_NAME"_FQ --sbatch "--mail-type=ALL,TIME_LIMIT_80 --job-name "$SWARM_NAME"_FQ")
echo "BWAMEM Alignment swarm ID: "$jobid1""
#
jobid2=$(swarm -f mergeDedup.swarm -g 36 -t 20 --gres=lscratch:350 --time 2-0 --module samtools,GATK/4.4.0.0 --logdir ~/job_outputs/samtools/"$SWARM_NAME"_merge --sbatch "--mail-type=ALL,TIME_LIMIT_90 --dependency=afterok:"$jobid1" --job-name "$SWARM_NAME"_Merge")
echo "Merge and Markduplicates Swarm ID: "$jobid2""
#
jobid3=$(swarm -f bqsr_BaseRecalibrator.swarm -g 6 -t 6 -b 10 --gres=lscratch:75 --time 4:00:00 --module GATK/4.4.0.0 --logdir ~/job_outputs/gatk/BaseRecalibrator/"$SWARM_NAME"_BR --sbatch "--mail-type=ALL,TIME_LIMIT_90 --dependency=afterok:"$jobid2" --job-name "$SWARM_NAME"_BR")
echo "BaseRecalibrator Swarm ID: "$jobid3""
#
jobid4=$(swarm -f bqsr_gatherBQSRReports.swarm -g 8 -t 4 -b 4 --gres=lscratch:75 --time 30:00 --partition=quick --module GATK/4.4.0.0 --logdir ~/job_outputs/gatk/GatherBQSRReports/"$SWARM_NAME"_GatherReports --sbatch "--mail-type=ALL,TIME_LIMIT_90 --dependency=afterok:"$jobid3" --job-name "$SWARM_NAME"_GatherBQSRReports")
echo "GatherBQSRReports Swarm ID: "$jobid4""
#
jobid5=$(swarm -f bqsr_ApplyBQSR.swarm -g 8 -t 6 -b 40 --gres=lscratch:350 --time 6:00:00 --module GATK/4.4.0.0 --logdir ~/job_outputs/gatk/BaseRecalibrator/"$SWARM_NAME"_ApplyBQSR --sbatch "--mail-type=ALL,TIME_LIMIT_80 --dependency=afterok:"$jobid4" --job-name "$SWARM_NAME"_ApplyBQSR")
echo "ApplyBQSR Swarm ID: "$jobid5""
#
jobid6=$(swarm -f bqsr_gatherBQSRBams.swarm -g 4 -t 4 --gres=lscratch:75 --time 4:00:00 --partition=quick --module GATK/4.4.0.0,samtools --logdir ~/job_outputs/gatk/GatherBams/"$SWARM_NAME"_GatherBams --sbatch "--mail-type=ALL,TIME_LIMIT_90 --dependency=afterok:"$jobid5" --job-name "$SWARM_NAME"_GatherBams")
echo "Gather BQSR Bams Swarm ID: "$jobid6""
#
jobid7=$(swarm -f haplotypecaller.swarm -g 8 -t 6 -b 5 --gres=lscratch:120 --time 48:00:00 --module GATK/4.4.0.0 --logdir ~/job_outputs/gatk/HaplotypeCaller/"$SWARM_NAME"_HC --sbatch "--mail-type=ALL,TIME_LIMIT_80 --dependency=afterok:"$jobid6" --job-name "$SWARM_NAME"_HC")
echo "HaplotypeCaller Swarm ID: "$jobid7""
#
jobid8=$(swarm -f hc_gathergvcfs.swarm -g 3 -t 4 -b 2 --gres=lscratch:75 --time 60:00 --partition=quick --module GATK/4.4.0.0,samtools --logdir ~/job_outputs/gatk/gathervcfs/"$SWARM_NAME"_GatherHC --sbatch "--mail-type=ALL,TIME_LIMIT_90 --dependency=afterok:"$jobid7" --job-name "$SWARM_NAME"_GatherHC")
echo "Gather gVCFs Swarm ID: "$jobid8""
}
#
#
## First is if mergeDedup fails, note if the BWA portion fails, its easier to simply just resubmit the entire pipeline so no function req'd for just BWA onwards.
#
mergeresubmit(){
jobid2=$(swarm -f mergeDedup.swarm -g 36 -t 20 --gres=lscratch:350 --time 2-0 --module samtools,GATK/4.4.0.0 --logdir ~/job_outputs/samtools/"$SWARM_NAME"_merge --sbatch "--mail-type=ALL,TIME_LIMIT_80 --job-name "$SWARM_NAME"_Merge")
echo "Merge and Markduplicates Swarm ID: "$jobid2""
#
jobid3=$(swarm -f bqsr_BaseRecalibrator.swarm -g 6 -t 6 -b 10 --gres=lscratch:75 --time 4:00:00 --module GATK/4.4.0.0 --logdir ~/job_outputs/gatk/BaseRecalibrator/"$SWARM_NAME"_BR --sbatch "--mail-type=ALL,TIME_LIMIT_90 --dependency=afterok:"$jobid2" --job-name "$SWARM_NAME"_BR")
echo "BaseRecalibrator Swarm ID: "$jobid3""
#
jobid4=$(swarm -f bqsr_gatherBQSRReports.swarm -g 8 -t 4 -b 4 --gres=lscratch:75 --time 30:00 --partition=quick --module GATK/4.4.0.0 --logdir ~/job_outputs/gatk/GatherBQSRReports/"$SWARM_NAME"_GatherReports --sbatch "--mail-type=ALL,TIME_LIMIT_90 --dependency=afterok:"$jobid3" --job-name "$SWARM_NAME"_GatherBQSRReports")
echo "GatherBQSRReports Swarm ID: "$jobid4""
#
jobid5=$(swarm -f bqsr_ApplyBQSR.swarm -g 8 -t 6 -b 40 --gres=lscratch:350 --time 6:00:00 --module GATK/4.4.0.0 --logdir ~/job_outputs/gatk/BaseRecalibrator/"$SWARM_NAME"_ApplyBQSR --sbatch "--mail-type=ALL,TIME_LIMIT_80 --dependency=afterok:"$jobid4" --job-name "$SWARM_NAME"_ApplyBQSR")
echo "ApplyBQSR Swarm ID: "$jobid5""
#
jobid6=$(swarm -f bqsr_gatherBQSRBams.swarm -g 4 -t 4 --gres=lscratch:75 --time 4:00:00 --partition=quick --module GATK/4.4.0.0,samtools --logdir ~/job_outputs/gatk/GatherBams/"$SWARM_NAME"_GatherBams --sbatch "--mail-type=ALL,TIME_LIMIT_90 --dependency=afterok:"$jobid5" --job-name "$SWARM_NAME"_GatherBams")
echo "Gather BQSR Bams Swarm ID: "$jobid6""
#
jobid7=$(swarm -f haplotypecaller.swarm -g 8 -t 6 -b 5 --gres=lscratch:120 --time 48:00:00 --module GATK/4.4.0.0 --logdir ~/job_outputs/gatk/HaplotypeCaller/"$SWARM_NAME"_HC --sbatch "--mail-type=ALL,TIME_LIMIT_80 --dependency=afterok:"$jobid6" --job-name "$SWARM_NAME"_HC")
echo "HaplotypeCaller Swarm ID: "$jobid7""
#
jobid8=$(swarm -f hc_gathergvcfs.swarm -g 3 -t 4 -b 2 --gres=lscratch:75 --time 60:00 --partition=quick --module GATK/4.4.0.0,samtools --logdir ~/job_outputs/gatk/gathervcfs/"$SWARM_NAME"_GatherHC --sbatch "--mail-type=ALL,TIME_LIMIT_90 --dependency=afterok:"$jobid7" --job-name "$SWARM_NAME"_GatherHC")
echo "Gather gVCFs Swarm ID: "$jobid8""
}
#
BQSRresubmit(){
jobid3=$(swarm -f bqsr_BaseRecalibrator.swarm -g 6 -t 6 -b 10 --gres=lscratch:75 --time 4:00:00 --module GATK/4.4.0.0 --logdir ~/job_outputs/gatk/BaseRecalibrator/"$SWARM_NAME"_BR --sbatch "--mail-type=ALL,TIME_LIMIT_80 --job-name "$SWARM_NAME"_BR")
echo "BaseRecalibrator Swarm ID: "$jobid3""
#
jobid4=$(swarm -f bqsr_gatherBQSRReports.swarm -g 8 -t 4 -b 4 --gres=lscratch:75 --time 30:00 --partition=quick --module GATK/4.4.0.0 --logdir ~/job_outputs/gatk/GatherBQSRReports/"$SWARM_NAME"_GatherReports --sbatch "--mail-type=ALL,TIME_LIMIT_90 --dependency=afterok:"$jobid3" --job-name "$SWARM_NAME"_GatherBQSRReports")
echo "GatherBQSRReports Swarm ID: "$jobid4""
#
jobid5=$(swarm -f bqsr_ApplyBQSR.swarm -g 8 -t 6 -b 40 --gres=lscratch:350 --time 6:00:00 --module GATK/4.4.0.0 --logdir ~/job_outputs/gatk/BaseRecalibrator/"$SWARM_NAME"_ApplyBQSR --sbatch "--mail-type=ALL,TIME_LIMIT_80 --dependency=afterok:"$jobid4" --job-name "$SWARM_NAME"_ApplyBQSR")
echo "ApplyBQSR Swarm ID: "$jobid5""
#
jobid6=$(swarm -f bqsr_gatherBQSRBams.swarm -g 4 -t 4 --gres=lscratch:75 --time 4:00:00 --partition=quick --module GATK/4.4.0.0,samtools --logdir ~/job_outputs/gatk/GatherBams/"$SWARM_NAME"_GatherBams --sbatch "--mail-type=ALL,TIME_LIMIT_80 --dependency=afterok:"$jobid5" --job-name "$SWARM_NAME"_GatherBams")
echo "Gather BQSR Bams Swarm ID: "$jobid6""
#
jobid7=$(swarm -f haplotypecaller.swarm -g 8 -t 6 -b 5 --gres=lscratch:120 --time 48:00:00 --module GATK/4.4.0.0 --logdir ~/job_outputs/gatk/HaplotypeCaller/"$SWARM_NAME"_HC --sbatch "--mail-type=ALL,TIME_LIMIT_80 --dependency=afterok:"$jobid6" --job-name "$SWARM_NAME"_HC")
echo "HaplotypeCaller Swarm ID: "$jobid7""
#
jobid8=$(swarm -f hc_gathergvcfs.swarm -g 3 -t 4 -b 2 --gres=lscratch:75 --time 60:00 --partition=quick --module GATK/4.4.0.0,samtools --logdir ~/job_outputs/gatk/gathervcfs/"$SWARM_NAME"_GatherHC --sbatch "--mail-type=ALL,TIME_LIMIT_90 --dependency=afterok:"$jobid7" --job-name "$SWARM_NAME"_GatherHC")
echo "Gather gVCFs Swarm ID: "$jobid8""
}
#
ReportsGatherresubmit(){
jobid4=$(swarm -f bqsr_gatherBQSRReports.swarm -g 8 -t 4 -b 4 --gres=lscratch:75 --time 30:00 --partition=quick --module GATK/4.4.0.0 --logdir ~/job_outputs/gatk/GatherBQSRReports/"$SWARM_NAME"_GatherReports --sbatch "--mail-type=ALL,TIME_LIMIT_90 --job-name "$SWARM_NAME"_GatherBQSRReports")
echo "GatherBQSRReports Swarm ID: "$jobid4""
#
jobid5=$(swarm -f bqsr_ApplyBQSR.swarm -g 8 -t 6 -b 40 --gres=lscratch:350 --time 6:00:00 --module GATK/4.4.0.0 --logdir ~/job_outputs/gatk/BaseRecalibrator/"$SWARM_NAME"_ApplyBQSR --sbatch "--mail-type=ALL,TIME_LIMIT_80 --dependency=afterok:"$jobid4" --job-name "$SWARM_NAME"_ApplyBQSR")
echo "ApplyBQSR Swarm ID: "$jobid5""
#
jobid6=$(swarm -f bqsr_gatherBQSRBams.swarm -g 4 -t 4 --gres=lscratch:75 --time 4:00:00 --partition=quick --module GATK/4.4.0.0,samtools --logdir ~/job_outputs/gatk/GatherBams/"$SWARM_NAME"_GatherBams --sbatch "--mail-type=ALL,TIME_LIMIT_80 --dependency=afterok:"$jobid5" --job-name "$SWARM_NAME"_GatherBams")
echo "Gather BQSR Bams Swarm ID: "$jobid6""
#
jobid7=$(swarm -f haplotypecaller.swarm -g 8 -t 6 -b 5 --gres=lscratch:120 --time 48:00:00 --module GATK/4.4.0.0 --logdir ~/job_outputs/gatk/HaplotypeCaller/"$SWARM_NAME"_HC --sbatch "--mail-type=ALL,TIME_LIMIT_80 --dependency=afterok:"$jobid6" --job-name "$SWARM_NAME"_HC")
echo "HaplotypeCaller Swarm ID: "$jobid7""
#
jobid8=$(swarm -f hc_gathergvcfs.swarm -g 3 -t 4 -b 2 --gres=lscratch:75 --time 60:00 --partition=quick --module GATK/4.4.0.0,samtools --logdir ~/job_outputs/gatk/gathervcfs/"$SWARM_NAME"_GatherHC --sbatch "--mail-type=ALL,TIME_LIMIT_90 --dependency=afterok:"$jobid7" --job-name "$SWARM_NAME"_GatherHC")
echo "Gather gVCFs Swarm ID: "$jobid8""
}
#
ApplyBQSRresubmit(){
jobid5=$(swarm -f bqsr_ApplyBQSR.swarm -g 8 -t 6 -b 40 --gres=lscratch:350 --time 6:00:00 --module GATK/4.4.0.0 --logdir ~/job_outputs/gatk/BaseRecalibrator/"$SWARM_NAME"_ApplyBQSR --sbatch "--mail-type=ALL,TIME_LIMIT_80 --job-name "$SWARM_NAME"_ApplyBQSR")
echo "ApplyBQSR Swarm ID: "$jobid5""
#
jobid6=$(swarm -f bqsr_gatherBQSRBams.swarm -g 4 -t 4 --gres=lscratch:75 --time 4:00:00 --partition=quick --module GATK/4.4.0.0,samtools --logdir ~/job_outputs/gatk/GatherBams/"$SWARM_NAME"_GatherBams --sbatch "--mail-type=ALL,TIME_LIMIT_80 --dependency=afterok:"$jobid5" --job-name "$SWARM_NAME"_GatherBams")
echo "Gather BQSR Bams Swarm ID: "$jobid6""
#
jobid7=$(swarm -f haplotypecaller.swarm -g 8 -t 6 -b 5 --gres=lscratch:120 --time 48:00:00 --module GATK/4.4.0.0 --logdir ~/job_outputs/gatk/HaplotypeCaller/"$SWARM_NAME"_HC --sbatch "--mail-type=ALL,TIME_LIMIT_80 --dependency=afterok:"$jobid6" --job-name "$SWARM_NAME"_HC")
echo "HaplotypeCaller Swarm ID: "$jobid7""
#
jobid8=$(swarm -f hc_gathergvcfs.swarm -g 3 -t 4 -b 2 --gres=lscratch:75 --time 60:00 --partition=quick --module GATK/4.4.0.0,samtools --logdir ~/job_outputs/gatk/gathervcfs/"$SWARM_NAME"_GatherHC --sbatch "--mail-type=ALL,TIME_LIMIT_90 --dependency=afterok:"$jobid7" --job-name "$SWARM_NAME"_GatherHC")
echo "Gather gVCFs Swarm ID: "$jobid8""
}
#
BQSRbamsGatherresubmit(){
jobid6=$(swarm -f bqsr_gatherBQSRBams.swarm -g 4 -t 4 --gres=lscratch:75 --time 4:00:00 --partition=quick --module GATK/4.4.0.0,samtools --logdir ~/job_outputs/gatk/GatherBams/"$SWARM_NAME"_GatherBams --sbatch "--mail-type=ALL,TIME_LIMIT_80 --job-name "$SWARM_NAME"_GatherBams")
echo "Gather BQSR Bams Swarm ID: "$jobid6""
#
jobid7=$(swarm -f haplotypecaller.swarm -g 8 -t 6 -b 5 --gres=lscratch:120 --time 48:00:00 --module GATK/4.4.0.0 --logdir ~/job_outputs/gatk/HaplotypeCaller/"$SWARM_NAME"_HC --sbatch "--mail-type=ALL,TIME_LIMIT_80 --dependency=afterok:"$jobid6" --job-name "$SWARM_NAME"_HC")
echo "HaplotypeCaller Swarm ID: "$jobid7""
#
jobid8=$(swarm -f hc_gathergvcfs.swarm -g 3 -t 4 -b 2 --gres=lscratch:75 --time 60:00 --partition=quick --module GATK/4.4.0.0,samtools --logdir ~/job_outputs/gatk/gathervcfs/"$SWARM_NAME"_GatherHC --sbatch "--mail-type=ALL,TIME_LIMIT_90 --dependency=afterok:"$jobid7" --job-name "$SWARM_NAME"_GatherHC")
echo "Gather gVCFs Swarm ID: "$jobid8""
}
#
HCresubmit(){
jobid7=$(swarm -f haplotypecaller.swarm -g 8 -t 6 -b 5 --gres=lscratch:120 --time 48:00:00 --module GATK/4.4.0.0 --logdir ~/job_outputs/gatk/HaplotypeCaller/"$SWARM_NAME"_HC --sbatch "--mail-type=ALL,TIME_LIMIT_80 --job-name "$SWARM_NAME"_HC")
echo "HaplotypeCaller Swarm ID: "$jobid7""
#
jobid8=$(swarm -f hc_gathergvcfs.swarm -g 3 -t 4 -b 2 --gres=lscratch:75 --time 60:00 --partition=quick --module GATK/4.4.0.0,samtools --logdir ~/job_outputs/gatk/gathervcfs/"$SWARM_NAME"_GatherHC --sbatch "--mail-type=ALL,TIME_LIMIT_90 --dependency=afterok:"$jobid7" --job-name "$SWARM_NAME"_GatherHC")
echo "Gather gVCFs Swarm ID: "$jobid8""
}
#
HCgatherresubmit(){
jobid1=$(swarm -f hc_gathergvcfs.swarm -g 3 -t 4 -b 2 --gres=lscratch:75 --time 60:00 --partition=quick --module GATK/4.4.0.0,samtools --logdir ~/job_outputs/gatk/gathervcfs/"$SWARM_NAME"_GatherHC --sbatch "--mail-type=ALL,TIME_LIMIT_90 --job-name "$SWARM_NAME"_GatherHC")
echo "Gather gVCFs Swarm ID: "$jobid1""
}
#
#
cd "$homedir"
#
## From here we create a long if elif else statement using the basis of whether a swarmfile is 0 bytes in size (aka, nothing in it) to determine at what stage to submitto the cluster
## However first, to determine syntax we'll perform a shortened version of that if elif else statement prior to full submission as a user check.
echo "Verify that the syntax of the following swarmfiles are correct:"
echo ""
if [ -s bwa_to_picard.swarm ]
then
	echo "BWA swarm:"
	head -n 1 bwa_to_picard.swarm
	echo ""
	echo "Merge and Markduplicates swarm:"
	head -n 1 mergeDedup.swarm
	echo ""
	echo "GATK BaseRecalibrator swarm:"
	head -n 1 bqsr_BaseRecalibrator.swarm
	echo ""
	echo "GatherBQSRReports swarm:"
	head -n 1 bqsr_gatherBQSRReports.swarm
	echo ""
	echo "ApplyBQSR swarm"
	head -n 1 bqsr_ApplyBQSR.swarm
	echo ""
	echo "GatherBQSR Bams swarm:"
	head -n 1 bqsr_gatherBQSRBams.swarm
	echo ""
	echo "HaplotypeCaller swarm:"
	head -n 1 haplotypecaller.swarm
	echo ""
	echo "Gather gVCFs swarm:"
	head -n 1 hc_gathergvcfs.swarm
	echo ""
elif [ -s mergeDedup.swarm ]
then
	echo "Merge and Markduplicates swarm:"
	head -n 1 mergeDedup.swarm
	echo ""
	echo "GATK BaseRecalibrator swarm:"
	head -n 1 bqsr_BaseRecalibrator.swarm
	echo ""
	echo "GatherBQSRReports swarm:"
	head -n 1 bqsr_gatherBQSRReports.swarm
	echo ""
	echo "ApplyBQSR swarm"
	head -n 1 bqsr_ApplyBQSR.swarm
	echo ""
	echo "GatherBQSR Bams swarm:"
	head -n 1 bqsr_gatherBQSRBams.swarm
	echo ""
	echo "HaplotypeCaller swarm:"
	head -n 1 haplotypecaller.swarm
	echo ""
	echo "Gather gVCFs swarm:"
	head -n 1 hc_gathergvcfs.swarm
	echo ""
elif [ -s bqsr_BaseRecalibrator.swarm ]
then
	echo "GATK BaseRecalibrator swarm:"
        head -n 1 bqsr_BaseRecalibrator.swarm
        echo ""
        echo "GatherBQSRReports swarm:"
        head -n 1 bqsr_gatherBQSRReports.swarm
        echo ""
        echo "ApplyBQSR swarm"
        head -n 1 bqsr_ApplyBQSR.swarm
        echo ""
        echo "GatherBQSR Bams swarm:"
        head -n 1 bqsr_gatherBQSRBams.swarm
        echo ""
        echo "HaplotypeCaller swarm:"
        head -n 1 haplotypecaller.swarm
        echo ""
        echo "Gather gVCFs swarm:"
        head -n 1 hc_gathergvcfs.swarm
        echo ""
elif [ -s bqsr_gatherBQSRReports.swarm ]
then
	echo "GatherBQSRReports swarm:"
        head -n 1 bqsr_gatherBQSRReports.swarm
        echo ""
        echo "ApplyBQSR swarm"
        head -n 1 bqsr_ApplyBQSR.swarm
        echo ""
        echo "GatherBQSR Bams swarm:"
        head -n 1 bqsr_gatherBQSRBams.swarm
        echo ""
        echo "HaplotypeCaller swarm:"
        head -n 1 haplotypecaller.swarm
        echo ""
        echo "Gather gVCFs swarm:"
        head -n 1 hc_gathergvcfs.swarm
        echo ""
elif [ -s bqsr_ApplyBQSR.swarm ]
then
	echo "ApplyBQSR swarm"
        head -n 1 bqsr_ApplyBQSR.swarm
        echo ""
        echo "GatherBQSR Bams swarm:"
        head -n 1 bqsr_gatherBQSRBams.swarm
        echo ""
        echo "HaplotypeCaller swarm:"
        head -n 1 haplotypecaller.swarm
        echo ""
        echo "Gather gVCFs swarm:"
        head -n 1 hc_gathergvcfs.swarm
        echo ""
elif [ -s bqsr_gatherBQSRBams.swarm ]
then
	echo "GatherBQSR Bams swarm:"
        head -n 1 bqsr_gatherBQSRBams.swarm
        echo ""
        echo "HaplotypeCaller swarm:"
        head -n 1 haplotypecaller.swarm
        echo ""
        echo "Gather gVCFs swarm:"
        head -n 1 hc_gathergvcfs.swarm
        echo ""
elif [ -s haplotypecaller.swarm ]
then
	echo "HaplotypeCaller swarm:"
        head -n 1 haplotypecaller.swarm
        echo ""
        echo "Gather gVCFs swarm:"
        head -n 1 hc_gathergvcfs.swarm
        echo ""
elif [ -s hc_gathergvcfs.swarm ]
then
	echo "Gather gVCFs swarm:"
        head -n 1 hc_gathergvcfs.swarm
        echo ""
else
	echo "Did you intend to resubmit samples to the pipeline?"
	exit 1
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
#
## If the user confirms that the syntax is correct we now run a if elif else loop to submit the right function to resubmit to the cluster. Currently how the script is written if a step fails, the rest of the pipeline fails. 
## With exception to a node failure, usually a failed step usually means a syntax error or missing files in the swarm file that failed so the user must go back to check and verify.
## Upon checking they may find some samples completed further than others which makes this loop feature useful. Saves on compute time.
#
if [ -s bwa_to_picard.swarm ]
then
	pipelinesubmit
elif [ -s mergeDedup.swarm ]
then
	mergeresubmit
elif [ -s bqsr_BaseRecalibrator.swarm ]
then
	BQSRresubmit
elif [ -s bqsr_gather_BQSRReports.swarm ]
then
	ReportsGatherresubmit
elif [ -s bqsr_ApplyBQSR.swarm ]
then
	ApplyBQSRresubmit
elif [ -s bqsr_gatherBQSRBams.swarm ]
then
	BQSRbamsGatherresubmit
elif [ -s haplotypecaller.swarm ]
then
	HCresubmit
elif [ -s hc_gathergvcfs.swarm ]
then
	HCgatherresubmit
else
	echo "Error! No swarmfiles to submit!"
	exit 1
fi
