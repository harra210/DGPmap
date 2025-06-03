#!/bin/bash
#
## This is the parent script to collect stats for samples pushed through the pipeline.
## It is intended to be the 2nd to last script run in the WGS alignment pipeline.
## This script is designed to run in a non-cluster environment but NOT on the cluster login node. Open an interactive node or run this script on the file system nodes.
#
homedir=$(pwd) # This is the main script directory
export homedir
#
cd ../tmp/
tmpdir=$(pwd)
cd "$homedir"
#
## Here we set up the variables for the script.
#
#
PARENT_DIR=
#
## getopts string
#
opts="d:h" # This string tells the getopt command that the 'd' flag is a required flag to input when invoking this script. The -d flag is the parent directory containing all of the samples needing stats collated.
#
# Gets the command name without path
cmd(){ echo `basename $0`; }
#
# Help command output
usage(){
        echo "\
                `cmd` [OPTION...]
        -d, --directory; Parent directory of pipeline processed files.
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
                                -d|--directory) export PARENT_DIR=$2; shift;;
                                -h|--help)      usage; exit 1;;
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
## Check to verify that the directory flag was passed on the command line
if [[ -z $PARENT_DIR ]]
then
	echo "Missing directory flag!"
	echo "Try '`cmd` -h or `cmd` --help' for more information."
	exit 1
else
	shift
fi
#
cd scripts/
scriptdir=$(pwd)
#
cd "$tmpdir"
> DatasetDirectories.tmp
#
cd "$PARENT_DIR"
#
find $PWD -mindepth 1 -maxdepth 1 -type d &> "$tmpdir"/DatasetDirectories.tmp
#
cd "$tmpdir"
#
IFS=,$'\n' read -d '' basedir < DatasetDirectories.tmp
#
cd "$PARENT_DIR" # Change into the parent directory
#
## Unlike the previous scripts that use functions, this script will not require the use of functions as the subscript has checks and all of the files required are already created.
#
confirm(){
	echo "WARNING! You are about to delete many files!"
	echo "This entails deleting:"
	echo "All BAM files"
	echo "Non-concatenated BQSR tables"
	echo "Non-final gVCF files"
	echo "Non-final statistics files"
	echo ""
	while true; do
	read -p "Are you absolutely sure you want to do this? [Y/N] " j
	case $j in
                [y/Y])
			cleanup; break;;
                [n/N])
                        echo -e "Exiting without deleting any files";
                        exit 1;;
                *)
                        echo "Invalid Option, enter (Y)es or (N)o"
                        confirm
                        ;;
        esac
done
}
#
finalconfirm(){
	echo "ATTENTION!!"
	echo "Final stats file not found! This may mean some files are missing!"
	while true; do
	read -r -p "Are you absolutely sure you want to force cleanup? [Y/N] " l
	case $l in
		[y/Y])
			forcecleanup; break;;
		[n/N])
			echo -e "Exiting this folder without deleting files. Check for missing files."
			exit 1;;
		*)
			echo "Invalid option, enter (Y)es or (N)o"
			finalconfirm
			;;
	esac
done
}	
#
forcecleanup(){
#echo "Successful forcecleanup!"
		printf "Performing cleanup for "$dir"..."; rm -r "$dir"/bwa_bam; rm -r "$dir"/dedup_bam; mv "$dir"/BQSR/tables/*.reports.list "$dir"/txt_files; rm -r "$dir"/BQSR; rm -r "$dir"/gVCF; rm "$dir"/txt_files/*.tmp; rm "$dir"/txt_files/bqsrreports.list; rm "$dir"/txt_files/*.knownsites.vcf.gz*; mv "$dir"/*.runs.table "$dir"/txt_files/ && printf "done"
}
#
cleanup(){
	for dir in ${basedir[@]}
	do
		cd "$dir";
		( for file in "$dir"/txt_files/*.stats.txt
		do
			if [ -e "$file" ]
			then
				printf "Performing cleanup for "$dir"..."; rm -r "$dir"/bwa_bam && rm -r "$dir"/dedup_bam && mv "$dir"/BQSR/tables/*.reports.list "$dir"/txt_files && rm -r "$dir"/BQSR && rm -r "$dir"/gVCF && rm "$dir"/txt_files/*.tmp && rm "$dir"/txt_files/bqsrreports.list && rm "$dir"/txt_files/*.knownsites.vcf.gz* && mv "$dir"/*.runs.table "$dir"/txt_files/ && printf "done\n"
			else
				finalconfirm
			fi
		done )
	done
}
#
confirm
