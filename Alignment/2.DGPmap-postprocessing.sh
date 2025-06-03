#!/bin/bash
#
## This script is to set Ostrander Lab internal names for samples in order to convert the BQSR bams to CRAM files
#
## Initial Variable Set
#
FILE_DIR=
FILE_IN=
SWARM_NAME=
#
## getopts string
#
opts="d:fs:h"
#
# Gets the command name without path
cmd(){ echo `basename $0`; }
#
# Help command output
usage(){
	echo "\
	`cmd` [OPTION...]
	-d, --directory; Parent directory of pipeline processed files.
	-f, --file; Space-delimited file containing Original name followed by new name.
	-s, --swarm-name; Set the base swarm name, any job specific names will be added on job submission.
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
				-d|--directory)	export FILE_DIR=$2; shift;;
				-f|--file)	export FILE_IN=$2; shift;;
				-s|--swarm-name)	export SWARM_NAME=$2; shift;;
				-h|--help)	usage; exit 1;;
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
#
## Verify that required flags were passed via command line
if [[ -z $FILE_DIR ]] && [[ -z $SWARM_NAME ]]
then
	echo "Missing required flags!"
	usage
	exit 1
elif [[ -z $FILE_DIR ]]
then
	echo "Missing input directory flag!"
	echo "Try '`cmd` -h or `cmd` --help' for more information."
	exit 1
elif [[ -z $SWARM_NAME ]]
then
	echo "Missing swarm-name flag!"
	echo "Try '`cmd` -h or `cmd` --help' for more information."
        exit 1
else
	shift
fi
#
## Running script variable definitions to be exported to subscript
#
CanFam4_Ref="/data/Ostrander/Resources/CanFam4_GSD/BWAMEM2/UU_Cfam_GSD_1.0_ROSY.fa"
export CanFam4_Ref
#
knownsites="/data/Ostrander/Resources/CanFam4_GSD/BWAMEM2/SRZ189891_722g.simp.header.CanineHDandAxiom_K9_HD.GSD_1.0.vcf.gz"
export knownsites
#
homedir=$(pwd)
export homedir
#
cd ../tmp
tmpdir=$(pwd)
export tmpdir
#
cd $homedir
#
cd swarmfiles
swarmdir=$(pwd)
export swarmdir
#
## Blank swarmfiles required to complete script
> postprocess-BQSR2CRAM.swarm
> postprocess-gVCFRecompress.swarm
> postprocess-knownsites.swarm
> postprocess-insertmetrics.swarm
> postprocess-alignmentmetrics.swarm
> postprocess-flagstat.swarm
#
cd $homedir
#
cd subscripts/
scriptdir=$(pwd)
#
cd "$tmpdir"
> DatasetDirectories.tmp
#
cd "$FILE_DIR"
#
find $PWD -mindepth 1 -maxdepth 1 -type d &> "$tmpdir"/DatasetDirectories.tmp
#
cd "$tmpdir"
#
IFS=,$'\n' read -d '' basedir < DatasetDirectories.tmp
#
cd "$FILE_DIR" # Change into the parent directory
#
bqsr2cram(){
	for dir in ${basedir[@]};
	do
		cd "$dir";
		( for file in "$dir"/*.BQSR.cram
		do
			if [ -e "$file" ]
			then
				echo "Skipping "$dir""; exit 0
			else
				echo "Building cram conversion swarmfile for "$dir""; bash "$scriptdir"/postprocess-BQSR2CRAM.sh
			fi
		done )
	done
}
#
gVCFcompress(){
	for dir in ${basedir[@]};
	do
		cd "$dir";
		( for file in "$dir"/*.g.vcf.gz
		do
			if [ -e "$file" ]
			then
				echo "Skipping "$dir""; exit 0
			else
				echo "Building gVCF recompression swarm for "$dir""; bash "$scriptdir"/postprocess-gvcfRecompress.sh
			fi
		done )
	done
}
#
knownsites(){
	for dir in ${basedir[@]};
        do
                cd "$dir";
                ( for file in "$dir"/txt_files/*.knownsites.vcf.gz.tbi
                do
                        if [ -e "$file" ]
                        then
				echo "Skipping "$dir""; exit 0
			else
                                echo "Building known sites swarm for "$dir""; bash "$scriptdir"/postprocess-knownsites.sh
                        fi
                done )
        done
}
#
flagstat(){
	for dir in ${basedir[@]};
	do
		cd "$dir";
		( for file in "$dir"/txt_files/*.flagstat
		do
			if [ -e "$file" ]
			then
				echo "Skipping "$dir""; exit 0
			else
				echo "Building flagstat swarm for "$dir""; bash "$scriptdir"/postprocess-flagstat.sh
			fi
		done )
	done
}
#
insertmetrics(){
	for dir in ${basedir[@]};
	do
		cd "$dir";
		( for file in "$dir"/txt_files/*.insert_size_metrics.txt
		do
			if [ -e "$file" ]
			then
				echo "Skipping "$dir""; exit 0
			else
				echo "Building insert metrics swarm for "$dir""; bash "$scriptdir"/postprocess-insertmetrics.sh
			fi
		done )
	done
}
#
alignmentmetrics(){
	for dir in ${basedir[@]};
	do
		cd "$dir";
		( for file in "$dir"/txt_files/*.alignment_summary_metrics
		do
			if [ -e "$file" ]
			then
				echo "Skipping "$dir""; exit 0
			else
				echo "Building alignment metrics swarm for "$dir""; bash "$scriptdir"/postprocess-alignmentmetrics.sh
			fi
		done )
	done
}
#
## Now that functions have been created, now we call them pipeline style with restart capabilities.
#
## Create the CRAM conversion swarmfile
#
echo "Generating CRAM conversion swarmfiles"
sleep 1
echo ""
bqsr2cram || { echo "Error generating conversion swarmfile"; exit 0; }
echo ""
#
## Create the gVCF Recompression swarmfile
#
echo "Generating gVCF Recompression swarmfile"
sleep 1
echo ""
gVCFcompress || { echo "Error generating recompression swarmfile"; exit 0; }
echo ""
#
## Create the knownsites swarmfile
#
echo "Generating knownsites VCF files"
sleep 1
echo ""
knownsites || { echo "Error generating knownsites swarmfile"; exit 0; }
echo ""
#
## Create the flagstat swarmfile
#
echo "Generating flagstat swarmfile"
sleep 1
echo ""
flagstat || { echo "Error generating flagstat swarmfile"; exit 0; }
echo ""
#
## Create the insert size metrics swarmfile
#
echo "Generating insert size metrics swarmfile"
sleep 1
echo ""
insertmetrics || { echo "Error generating insert size metrics swarmfile"; exit 0; }
echo ""
#
## Create the alignment metrics swarmfile
#
echo "Generating the alignment metrics swarmfile"
sleep 1
echo ""
alignmentmetrics || { echo "Error generating alignment metrics swarmfile"; exit 0; }
echo ""
#
## Now all swarmfiles are generated and below are all of the possible functions of how to submit to the cluster, fully or partially.
#
fullsubmit(){
jobid1=$(swarm -f postprocess-BQSR2CRAM.swarm -g 8 -t 6 --gres=lscratch:175 --time 24:00:00 --module GATK/4.6.0.0,samtools --logdir ~/job_outputs/gatk/PrintReads/"$SWARM_NAME"_PR --sbatch "--mail-type=ALL,TIME_LIMIT_80 --job-name "$SWARM_NAME"_PR")
echo "CRAM conversion Swarm ID: " $jobid1
#
jobid2=$(swarm -f postprocess-gVCFRecompress.swarm -g 8 -t 8 --gres=lscratch:20 --time 24:00:00 --module samtools,bcftools --logdir ~/job_outputs/Recompress/"$SWARM_NAME" --sbatch "--mail-type=ALL,TIME_LIMIT_80 --job-name "$SWARM_NAME"_Recompress")
echo "gVCF Recompress Swarm ID: " $jobid2
#
jobid3=$(swarm -f postprocess-knownsites.swarm -g 8 -t 8 --gres=lscratch:15 --time 24:00:00 --module GATK/4.6.0.0 --logdir ~/job_outputs/gatk/Knownsites/"$SWARM_NAME"_KS --sbatch "--mail-type=ALL,TIME_LIMIT_80 --dependency=afterok:"$jobid2" --job-name "$SWARM_NAME"_knownsites")
echo "Knownsites Swarm ID: " $jobid3
#
jobid4=$(swarm -f postprocess-flagstat.swarm -b 10 --time 1:00:00 --module samtools --logdir ~/job_outputs/samtools/flagstat/"$SWARM_NAME"_flagstat --sbatch "--mail-type=ALL,TIME_LIMIT_80 --job-name "$SWARM_NAME"_flagstat --dependency=afterok:"$jobid1"")
echo "Flagstat Swarm ID: " $jobid4
#
jobid5=$(swarm -f postprocess-insertmetrics.swarm -g 8 -t 4 --time 5:00:00 --module GATK/4.6.0.0,R/4.4 --logdir ~/job_outputs/gatk/InsertMetrics/"$SWARM_NAME"_InsertMet --sbatch "--mail-type=ALL,TIME_LIMIT_80 --job-name "$SWARM_NAME"_InsertMet --dependency=afterok:"$jobid1"")
echo "Insert Size Metrics Swarm ID: " $jobid5
#
jobid6=$(swarm -f postprocess-alignmentmetrics.swarm -g 8 -t 4 --time 5:00:00 --module GATK/4.6.0.0 --logdir ~/job_outputs/gatk/AlignmentMetrics/"$SWARM_NAME"_AlignMet --sbatch "--mail-type=ALL,TIME_LIMIT_80 --job-name "$SWARM_NAME"_AlignMet --dependency=afterok:"$jobid1"")
echo "Alignment Size Metrics Swarm ID: " $jobid6
}
#
cramresub(){
jobid1=$(swarm -f postprocess-BQSR2CRAM.swarm -g 8 -t 6 --gres=lscratch:175 --time 24:00:00 --module GATK/4.6.0.0,samtools --logdir ~/job_outputs/gatk/PrintReads/"$SWARM_NAME"_PR --sbatch "--mail-type=ALL,TIME_LIMIT_80 --job-name "$SWARM_NAME"_PR")
echo "CRAM conversion Swarm ID: " $jobid1
#
jobid4=$(swarm -f postprocess-flagstat.swarm -b 10 --time 1:00:00 --module samtools --logdir ~/job_outputs/samtools/flagstat/"$SWARM_NAME"_flagstat --sbatch "--mail-type=ALL,TIME_LIMIT_80 --job-name "$SWARM_NAME"_flagstat --dependency=afterok:"$jobid1"")
echo "Flagstat Swarm ID: " $jobid4
#
jobid5=$(swarm -f postprocess-insertmetrics.swarm -g 8 -t 4 --time 5:00:00 --module GATK/4.6.0.0,R/4.4 --logdir ~/job_outputs/gatk/InsertMetrics/"$SWARM_NAME"_InsertMet --sbatch "--mail-type=ALL,TIME_LIMIT_80 --job-name "$SWARM_NAME"_InsertMet --dependency=afterok:"$jobid1"")
echo "Insert Size Metrics Swarm ID: " $jobid5
#
jobid6=$(swarm -f postprocess-alignmentmetrics.swarm -g 8 -t 4 --time 5:00:00 --module GATK/4.6.0.0 --logdir ~/job_outputs/gatk/AlignmentMetrics/"$SWARM_NAME"_AlignMet --sbatch "--mail-type=ALL,TIME_LIMIT_80 --job-name "$SWARM_NAME"_AlignMet --dependency=afterok:"$jobid1"")
echo "Alignment Size Metrics Swarm ID: " $jobid6
}
#
recompressonlyresub(){
jobid2=$(swarm -f postprocess-gVCFRecompress.swarm -g 8 -t 8 --gres=lscratch:20 --time 24:00:00 --module samtools,bcftools --logdir ~/job_outputs/Recompress/"$SWARM_NAME" --sbatch "--mail-type=ALL,TIME_LIMIT_80 --job-name "$SWARM_NAME"_Recompress")
echo "gVCF Recompress Swarm ID: " $jobid2
}
#
knownsitesonlyresub(){
jobid3=$(swarm -f postprocess-knownsites.swarm -g 8 -t 8 --gres=lscratch:15 --time 24:00:00 --module GATK/4.6.0.0 --logdir ~/job_outputs/gatk/Knownsites/"$SWARM_NAME"_KS --sbatch "--mail-type=ALL,TIME_LIMIT_80 --job-name "$SWARM_NAME"_knownsites")
echo "Knownsites Swarm ID: " $jobid3
}
gvcfandKSresub(){
jobid2=$(swarm -f postprocess-gVCFRecompress.swarm -g 8 -t 8 --gres=lscratch:20 --time 24:00:00 --module samtools,bcftools --logdir ~/job_outputs/Recompress/"$SWARM_NAME" --sbatch "--mail-type=ALL,TIME_LIMIT_80 --job-name "$SWARM_NAME"_Recompress")
echo "gVCF Recompress Swarm ID: " $jobid2
#
jobid3=$(swarm -f postprocess-knownsites.swarm -g 8 -t 8 --gres=lscratch:15 --time 24:00:00 --module GATK/4.6.0.0 --logdir ~/job_outputs/gatk/Knownsites/"$SWARM_NAME"_KS --sbatch "--mail-type=ALL,TIME_LIMIT_80 --dependency=afterok:"$jobid2" --job-name "$SWARM_NAME"_knownsites")
echo "Knownsites Swarm ID: " $jobid3
}
#
fullmetricsresub(){
jobid4=$(swarm -f postprocess-flagstat.swarm -b 10 --time 1:00:00 --module samtools --logdir ~/job_outputs/samtools/flagstat/"$SWARM_NAME"_flagstat --sbatch "--mail-type=ALL,TIME_LIMIT_80 --job-name "$SWARM_NAME"_flagstat")
echo "Flagstat Swarm ID: " $jobid4
#
jobid5=$(swarm -f postprocess-insertmetrics.swarm -g 8 -t 4 --time 5:00:00 --module GATK/4.6.0.0,R/4.4 --logdir ~/job_outputs/gatk/InsertMetrics/"$SWARM_NAME"_InsertMet --sbatch "--mail-type=ALL,TIME_LIMIT_80 --job-name "$SWARM_NAME"_InsertMet")
echo "Insert Size Metrics Swarm ID: " $jobid5
#
jobid6=$(swarm -f postprocess-alignmentmetrics.swarm -g 8 -t 4 --time 5:00:00 --module GATK/4.6.0.0 --logdir ~/job_outputs/gatk/AlignmentMetrics/"$SWARM_NAME"_AlignMet --sbatch "--mail-type=ALL,TIME_LIMIT_80 --job-name "$SWARM_NAME"_AlignMet")
echo "Alignment Size Metrics Swarm ID: " $jobid6
}
#
metricsresub(){
jobid5=$(swarm -f postprocess-insertmetrics.swarm -g 8 -t 4 --time 12:00:00 --module GATK/4.6.0.0,R/4.4 --logdir ~/job_outputs/gatk/InsertMetrics/"$SWARM_NAME"_InsertMet --sbatch "--mail-type=ALL,TIME_LIMIT_80 --job-name "$SWARM_NAME"_InsertMet")
echo "Insert Size Metrics Swarm ID: " $jobid5
#
jobid6=$(swarm -f postprocess-alignmentmetrics.swarm -g 8 -t 4 --time 12:00:00 --module GATK/4.6.0.0 --logdir ~/job_outputs/gatk/AlignmentMetrics/"$SWARM_NAME"_AlignMet --sbatch "--mail-type=ALL,TIME_LIMIT_80 --job-name "$SWARM_NAME"_AlignMet")
echo "Alignment Size Metrics Swarm ID: " $jobid6
}
#
insertresub(){
jobid5=$(swarm -f postprocess-insertmetrics.swarm -g 8 -t 4 --time 12:00:00 --module GATK/4.6.0.0,R/4.4 --logdir ~/job_outputs/gatk/InsertMetrics/"$SWARM_NAME"_InsertMet --sbatch "--mail-type=ALL,TIME_LIMIT_80 --job-name "$SWARM_NAME"_InsertMet")
echo "Insert Size Metrics Swarm ID: " $jobid5
}
#
alignmentresub(){
jobid6=$(swarm -f postprocess-alignmentmetrics.swarm -g 8 -t 4 --time 24:00:00 --module GATK/4.6.0.0 --logdir ~/job_outputs/gatk/AlignmentMetrics/"$SWARM_NAME"_AlignMet --sbatch "--mail-type=ALL,TIME_LIMIT_80 --job-name "$SWARM_NAME"_AlignMet")
echo "Alignment Size Metrics Swarm ID: " $jobid6
}
#
## Now we make the if elif else statement to determine if a swarmfile has 0 bytes in size to determine what to submit to the cluster.
## Though first we make a shortened version of that statement to show the swarmfile to the user to verify correct syntax
#
cd "$homedir"
#
echo "Verify that the syntax of the following swarmfiles are correct:"
echo ""
if [ -s postprocess-BQSR2CRAM.swarm ]
then
	echo "CRAM conversion swarm:"
	head -n 1 postprocess-BQSR2CRAM.swarm
	echo ""
	echo "gVCF Recompress swarm:"
	head -n 1 postprocess-gVCFRecompress.swarm
	echo ""
	echo "Knownsites swarm:"
	head -n 1 postprocess-knownsites.swarm
	echo ""
	echo "Flagstat swarm:"
	head -n 1 postprocess-flagstat.swarm
	echo ""
	echo "Insert Metrics swarm:"
	head -n 1 postprocess-insertmetrics.swarm
	echo ""
	echo "Alignment Metrics swarm:"
	head -n 1 postprocess-alignmentmetrics.swarm
	echo ""
elif [ -s postprocess-BQSR2CRAM.swarm ] && [ ! -s postprocess-gVCFRecompress.swarm ] && [ ! -s postprocess-knownsites.swarm ]
then
	echo "CRAM conversion swarm:"
        head -n 1 postprocess-BQSR2CRAM.swarm
        echo ""
	echo "Flagstat swarm:"
        head -n 1 postprocess-flagstat.swarm
        echo ""
        echo "Insert Metrics swarm:"
        head -n 1 postprocess-insertmetrics.swarm
        echo ""
        echo "Alignment Metrics swarm:"
        head -n 1 postprocess-alignmentmetrics.swarm
        echo ""
elif [ ! -s postprocess-BQSR2CRAM.swarm ] && { [ -s postprocess-gVCFRecompress.swarm ] || [ -s postprocess-knownsites.swarm ]; }
then
	if [ -s postprocess-gVCFRecompress.swarm ] && [ -s postprocess-knownsites.swarm ]
	then
		echo "gVCF Recompress swarm:"
                head -n 1 postprocess-gVCFRecompress.swarm
                echo ""
                echo "Knownsites swarm:"
                head -n 1 postprocess-knownsites.swarm
                echo ""
	elif [ ! -s postprocess-gVCFRecompress.swarm ] && [ -s postprocess-knownsites.swarm ]
	then
		echo "Knownsites swarm:"
                head -n 1 postprocess-knownsites.swarm
                echo ""
	else
		echo "gVCF Recompress swarm:"
                head -n 1 postprocess-gVCFRecompress.swarm
                echo ""
	fi
elif [ -s postprocess-flagstat.swarm ] && [ -s postprocess-insertmetrics.swarm ] && [ -s postprocess-alignmentmetrics.swarm ]
then
	echo "Flagstat swarm:"
	head -n 1 postprocess-flagstat.swarm
	echo ""
	echo "Insert Metrics swarm:"
	head -n 1 postprocess-insertmetrics.swarm
	echo ""
	echo "Alignment Metrics swarm:"
	head -n 1 postprocess-alignmentmetrics.swarm
	echo ""
elif [ ! -s postprocess-flagstat.swarm ] && { [ -s postprocess-insertmetrics.swarm ] || [ -s postprocess-alignmentmetrics.swarm ]; }
then
	if [ -s postprocess-insertmetrics.swarm ] && [ -s postprocess-alignmentmetrics.swarm ]
	then
		echo "Insert Metrics swarm:"
		head -n 1 postprocess-insertmetrics.swarm
		echo ""
		echo "Alignment Metrics swarm:"
		head -n 1 postprocess-alignmentmetrics.swarm
		echo ""
	elif [ ! -s postprocess-alignmentmetrics.swarm ]
	then
		echo "Insert Metrics swarm:"
		head -n 1 postprocess-insertmetrics.swarm
		echo ""
	else
		echo "Alignment Metrics swarm:"
		head -n 1 postprocess-alignmentmetrics.swarm
		echo""
	fi
else
	echo "Error! No swarmfiles to submit?"; exit 1
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
if [ -s postprocess-BQSR2CRAM.swarm ]
then
	fullsubmit
elif [ -s postprocess-BQSR2CRAM.swarm ] && [ ! -s postprocess-gVCFRecompress.swarm ] && [ ! -s postprocess-knownsites.swarm ]
then
	CRAMresub
elif [ ! -s postprocess-BQSR2CRAM.swarm ] && { [ -s postprocess-gVCFRecompress.swarm ] || [ -s postprocess-knownsites.swarm ]; }
then
	if [ -s postprocess-gVCFRecompress.swarm ] && [ -s postprocess-knownsites.swarm ]
	then
		gvcfandKSresub
	elif [ ! -s postprocess-gVCFRecompress.swarm ] && [ -s postprocess-knownsites.swarm ]
	then
		knownsitesonlyresub
	else
		recompressonlyresub
	fi
elif [ -s postprocess-flagstat.swarm ] && [ -s postprocess-insertmetrics.swarm ] && [ -s postprocess-alignmentmetrics.swarm ]
then
	fullmetricsresub
elif [ ! -s postprocess-flagstat.swarm ] && { [ -s postprocess-insertmetrics.swarm ] || [ -s postprocess-alignmentmetrics.swarm ]; }
then
	if [ -s postprocess-insertmetrics.swarm ] && [ -s postprocess-alignmentmetrics.swarm ]
	then
		metricsresub
	elif [ ! -s postprocess-alignmentmetrics.swarm ]
	then
		insertresub
	else
		alignmentresub
	fi
else
	echo "Error! No swarmfiles submitted!"; exit 1
fi
