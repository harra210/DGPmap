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
cd scripts/
scriptdir=$(pwd)
export scriptdir
combinepy=""$scriptdir"/postprocess-combinestats.py"
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
#done
#
# Handle positional arguments
if [ -n "$*" ]; then
        echo "`cmd`: Extra arguments -- $*"
        echo "Try '`cmd` -h' for more information."
        exit 1
fi
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
for dir in ${basedir[@]}
do
	cd "$dir"
	IFS=,$'\n' read -d '' -r -a samplename < "$dir"/txt_files/rename.txt; echo "Collecting stats for "${samplename[0]}""
	export file_dir=$(pwd)
	cd txt_files/
	export txt_dir=$(pwd)
	cd $file_dir
	cd BQSR/	
	export bqsr_dir=$(pwd)
	cd $file_dir
	cd gVCF
	export gvcf_dir=$(pwd)
	python "$combinepy"	
done
