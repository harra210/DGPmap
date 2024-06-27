#!/bin/bash
#
## This script is intended to utilize the prefetch and fastq-dump tools of the SRA Toolkit package released by NCBI to download uploaded sequence data locally 
# and sort files by SAMN numbers, if given, or parse out files with a folder hierarchy compatible with the DGPmap pipeline.
#
## Flag Processing section to parse flags from the script invokation. Required flags are Output location (o) and Swarm Name (s). Optional flag is an input file (i)
# if the dataset to be downloaded is known to have multiple reads per sample. If there is only 1 SRA run per sample then an if statement will create folders properly.
#
# getopts string
opts="io:s:" # input (optional), output and swarm_name (req'd)
#
# Variables to set:
INPUT=
OUT_DIR=
SWARM_NAME=
#
# Gets the command name without path
cmd(){ echo `basename $0`; }
#
# Help command output
usage(){
	echo "\
		`cmd` [OPTION...]
	-i, --input; 2 column tab-delimited file containing SRA Sample numbers (SAMN) and their associated SRA Run accession number (SRR) [OPTIONAL]
	-o, --output; Parent output directory to download and expand files [Req'd]
	-s, --swarm-name; Swarm name to give your swarm, use underscores in place of spaces [Req'd]
	-h, --help; Print this message and exit.
	" | column -t -s ";"
}
#
# Error message
error(){
	echo "`cmd`: invalid option -- '$1'";
	echo "Try '`cmd` -h or `cmd` --help' for more information.";
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
				-i|--input)	export INPUT=$2; shift;;
				-o|--output)	export OUT_DIR=$2; shift;;
				-s|--swarm-name)	export SWARM_NAME=$2; shift;;
				-h|--help)	usage; exit 1;;
				--*)	error $1;;
				 -*)     if [ $pass -eq 1 ]; then ARGS="$ARGS $1";
                                        else error $1; fi;;
                                esac;;
                        *) if [ $pass -eq 1 ]; then ARGS="$ARGS $1";
                                else error $1; fi;;
                esac
                shift
        done
        if [ $pass -eq 1 ]; then ARGS=`getopt $opts $ARGS`
                if [ $? != 0 ]; then error; exit 2; fi; set -- $ARGS
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
## Variable defining sections that are unrelated to the flag invokations
#
homedir=$(pwd) # Make the current directory the script is running out of the "Home" directory.
# Blank the swarmfiles
> SRAtoolkit-prefetch.swarm
> SRAtoolkit-fastqdump.swarm
#
cd ../tmp
tmpdir=$(pwd)
> SRA_sampledirectories.tmp
> SRA_rundirectories.tmp
cd "$homedir"
#
#mkdir -p "$OUT_DIR"
#
## From here we can perform an if/then statement based on the status of input flag; whether or not a user-supplied file is given to parse.
#
## This first section is to input a user-supplied file
if [ ! -z "$INPUT" ]; then # Test to see if variable INPUT is not null
	awk '{print $1}' "$INPUT" &> "$tmpdir"/SRA_sampledirectories.tmp # Place the first column of the the user supplied file (containing SAMN numbers) into a tmp file for import into an array
	awk '{print $2}' "$INPUT" &> "$tmpdir"/SRA_rundirectories.tmp
	IFS=,$'\n' read -d '' -r -a SAMNdir < "$tmpdir"/SRA_sampledirectories.tmp
	IFS=,$'\n' read -d '' -r -a SRR < "$tmpdir"/SRA_rundirectories.tmp
	declare -a SAMNdir # SAMN accession directory
	declare -a SRR # SRA Run accession directory
	unset IFS
	#
	cd $OUT_DIR
	#
	for dir in ${SAMNdir[@]}
	do
		mkdir -p "$dir" # This creates the directories based on SAMN number with the -p flag allowing for both multi-accession vs single-accession
	done
	#
	echo "Successfully created directories for the following SAMN accessions:"
	printf "%s\n" "${SAMNdir[@]}" | sort -u # This lists the accessions uniquely so multiple values from the supplied file aren't reported
	echo ""
	#
	## Here we link SAMNdir and SRR together in a for loop across the previously created arrays. This is why no sorting or anything else was done to the original file.
	#
	for ((i = 0; i < ${#SAMNdir[@]}; i++))
	do
		echo "prefetch -X 75g -O "$OUT_DIR""${SAMNdir[$i]}" "${SRR[$i]}"" >> "$homedir"/SRAtoolkit-prefetch.swarm
	done
	#
else
	echo "Input SRA runs (SRRxxx) and separate each run by a space. Hit enter and type done when finished."
	while read terminputs
	do
		[ "$terminputs" == "done" ] && break
		man_input=("${array[@]}" $terminputs)
	done
	echo "These are the SRA runs requested to prefetch FastQ files of:"
	echo "${man_input[@]}"
	#
	## Here we prompt the user to verify that the runs they want to download are correct
	runverify(){
		read -r -p "Are these the SRA accessions you wish to prefetch? [Y/N] " i
		case $i in
			[yY])
				echo "";;
			[nN])
				echo -e "Exit and retry the script"
				exit 1;;
			*)
				echo "Invalid Option"
				runverify
				;;
		esac
	}
	runverify
	#
	### For later: make directories out of the manual inputs as sratoolkit will not automatically make directories.
	cd $OUT_DIR
	#
	for dir in ${man_input[@]}
	do
		mkdir -p "$dir"
	done
	#
	echo "Successfully created directories for the following SRA Run accessions: "
	printf "%s\n" "${man_input[@]}" | sort -u 
	echo ""
	#
	## Now we generate the swarm file for the manual version of the for loop.
	#
	for name in "${man_input[@]}"
	do
		echo "prefetch -X 75g -O "$OUT_DIR""$name" "$name"" >> "$homedir"/SRAtoolkit-prefetch.swarm
	done
	#
fi
#
## Now we prompt the user to verify that the prefetch commands are formatted correctly.
more "$homedir"/SRAtoolkit-prefetch.swarm
echo ""
prefetchverify(){
	read -r -p "Is the previous swarmfile formatted correctly? [Y/N] " j
	case $j in
		[y/Y])
			echo "Continuing to fastq-dump";;
		[n/N])
			echo -e "Exit and retry the script"
			exit 1;;
		*)
			echo "Invalid Option, enter (Y)es or (N)o"
			prefetchverify
			;;
	esac
}
prefetchverify
#
## The prefetch section of the script is complete, now we transition into the fastq-dump portion of the script
#
if [ ! -z "$INPUT" ]; then
	for ((i = 0; i < ${#SAMNdir[@]}; i++))
	do
		echo "fastq-dump --split-files --gzip --dumpbase -O "$OUT_DIR""${SAMNdir[$i]}"/"${SRR[$i]}" "$OUT_DIR""${SAMNdir[$i]}"/"${SRR[$i]}"/"${SRR[$i]}".sra && mv "$OUT_DIR""${SAMNdir[$i]}"/"${SRR[$i]}"/"${SRR[$i]}"_1.fastq.gz "$OUT_DIR""${SAMNdir[$i]}"/ && mv "$OUT_DIR""${SAMNdir[$i]}"/"${SRR[$i]}"/"${SRR[$i]}"_2.fastq.gz "$OUT_DIR""${SAMNdir[$i]}"/" >> "$homedir"/SRAtoolkit-fastqdump.swarm
	done
else
	for name in "${man_input[@]}"
	do
		echo "fastq-dump --split-files --gzip --dumpbase -O "$OUT_DIR""$name"/"$name" "$OUT_DIR""$name"/"$name"/"$name".sra && mv "$OUT_DIR""$name"/"$name"/"$name"_1.fastq.gz "$OUT_DIR""$name"/ && mv "$OUT_DIR""$name"/"$name"/"$name"_2.fastq.gz "$OUT_DIR""$name"" >> "$homedir"/SRAtoolkit-fastqdump.swarm
	done
fi
#
more "$homedir"/SRAtoolkit-fastqdump.swarm
echo ""
fqverify(){
	read -r -p "Is the previous swarmfile formatted correctly? [Y/N] " k
	case $k in
		[y/Y])
			echo "Continuing to submitting pipeline to cluster";;
		[n/N])
			echo -e "Exit and retry script"
			exit 1;;
		*)
			echo "Invalid Option, enter (Y)es or (N)o"
			fqverify
			;;
	esac
}
fqverify
#
jobid1=$(swarm -f SRAtoolkit-prefetch.swarm -g 4 -t 4 --time 24:00:00 --gres=lscratch:75 --module sratoolkit --logdir ~/job_outputs/SRA_Toolkit/Prefetch/"$SWARM_NAME" --sbatch "--mail-type=ALL,TIME_LIMIT_80 --job-name "$SWARM_NAME"_Prefetch")
echo "SRAtoolkit Prefetch Swarm ID: " $jobid1
#
jobid2=$(swarm -f SRAtoolkit-fastqdump.swarm -g 4 -t 4 --time 36:00:00 --gres=lscratch:75 --module sratoolkit --logdir ~/job_outputs/SRA_Toolkit/fastq-dump/"$SWARM_NAME" --sbatch "--mail-type=ALL,TIME_LIMIT_80 --dependency=afterok:"$jobid1" --job-name "$SWARM_Name"_Fqdump")
echo "Sratoolkit Fastq-dump Swarm ID: " $jobid2
#
## End SRAtoolkit pipeline
