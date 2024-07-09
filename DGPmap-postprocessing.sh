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
## Running script variable definitions to be exported to subscript
#
CanFam4_Ref="/data/Ostrander/Resources/CanFam4_GSD/BWAMEM2/UU_Cfam_GSD_1.0_ROSY.fa"
export CanFam4_Ref
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
## Blank swarmfiles required to complete script
> BQSR2CRAM_Rename.swarm
> gVCFRecompress.swarm
#
cd scripts/
scriptdir=$(pwd)
#
cd "$tmpdir"
> DatasetDirectories.tmp
#
cd "$FILE_DIR"
#
find $PWD -mindepth 1 -maxdepth 1 -type d &> "$tmpdir"/DatasetDirectories.tmp
#
cd $tmpdir
#
IFS=,$'\n' read -d '' basedir < DatasetDirectories.tmp
#
cd $FILE_DIR # Change into the parent directory
#
bqsr2cram(){
	for dir in ${basedir[@]};
	do
		cd "$dir";
		( for file in "$dir"/*.runs.table
		do
			if [ -e "$file" ]
			then
				echo "Building cram conversion swarmfile for "$dir""; bash "$scriptdir"/postprocess-BQSR2CRAM.sh
			else
				echo "Skipping "$dir""; exit 0
			fi
		done )
	done
}
#
gVCFcompress(){
	for dir in ${basedir[@]};
	do
		cd "$dir";
		( for file in "$dir"/gVCF/*_g.vcf.gz
		do
			if [ -e "$file" ]
			then
				echo "Building gVCF recompression swarm for "$dir""; bash "$scriptdir"/postprocess-gVCFrecompress.sh
			else
				echo "Skipping "$dir""; exit 0
			fi
		done )
	done
}
#
