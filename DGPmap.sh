#!/bin/bash
#
# This begins after the creation of sample tables and will use the table to parse out all of the
# swarmfiles to submit to the cluster. This script will also call additional scripts as needed and be more in line
# with the original Dog10K pipeline.
#
# Script will have section which can have flags added in to allow for the ability to use different alignments.
# Would be done at a later date.
#
#set -e
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
# Here, we blank all of the swarmfiles to be used from alignment to completion
> bwa_to_picard.swarm
#> mergeDedup.swarm
#> bqsr_BaseRecalibrator.swarm
#> bqsr_gatherBQSRReports.swarm
#> bqsr_ApplyBQSR.swarm
#> bqsr_gatherBQSRBams.swarm
#> haplotypecaller.swarm
#
#
cd scripts/
scriptdir=$(pwd)
#
cd $tmpdir
> DatasetDirectories.tmp
#
echo "What parent directory are your fastq files that you want to align?"
if ! read -e -t 60 FILE_DIR || (($? > 128)); then
        echo "User timeout" >&2
        exit 1
fi
#
echo "What do you want to name your base swarm?"
if ! read -e -t 60 SWARM_NAME || (($? > 128)); then
        echo "User timeout" >&2
        exit 1
fi
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
# Begin the loop section to work on generating swarmfiles.
# Firstly, we generate the alignment sample tables
for dir in ${basedir[@]};
	do
		( [ -d "$dir" ] && cd $dir && echo "Entering into $dir and generating alignment tables" && bash "$scriptdir"/table-generation.sh )
done
#
# This part builds the BWAMEM2 alignment swarmfile
for dir in ${basedir[@]};
	do
		( [ -d "$dir" ] && cd $dir && echo "Entering into $dir and generating alignment swarm" && bash "$scriptdir"/process-BWA.sh )
done
#
# This part builds the samtools and GATK MarkDuplicates swarmfile
for dir in ${basedir[@]};
        do
                ( [ -d "$dir" ] && cd $dir && echo "Entering into $dir and generating alignment swarm" && bash "$scriptdir"/process-sortmd.sh )
done
#
