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
FQ_IN="NISC"
#
# getops string
opts="f:"
#
# Gets the command name without path
cmd(){ echo `basename $0`; }
#
# Help command output
usage(){
        echo "\
                `cmd` [OPTION...]
        -f, --fastq; Set input fastq style format (default:$OPT)
        " | column -t -s ";"
}
#
# Error message
error(){
        echo "`cmd`: invalid option -- '$1'";
        echo "Try '`cmd` -h' for more information.";
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
                                -f|--fastq)     export FQ_IN=$2; shift;;
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
echo "$FQ_IN"
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
cd ../tmp
tmpdir=$(pwd)
export tmpdir
#
cd $homedir
#
# Here, we blank all of the swarmfiles to be used from alignment to completion
#> bwa_to_picard.swarm
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
#
for dir in ${basedir[@]};
        do
                ( [ -d "$dir" ] && cd $dir && echo "Entering into $dir and generating alignment tables" && bash "$scriptdir"/table-generation.sh )
done
#
# This part builds the BWAMEM2 alignment swarmfile
#
for dir in ${basedir[@]};
        do
                ( [ -d "$dir" ] && cd $dir && echo "Entering into $dir and generating alignment swarm" && bash "$scriptdir"/process-BWA.sh )
done
#
# This part builds the samtools and GATK MarkDuplicates swarmfile
#
for dir in ${basedir[@]};
        do
                ( [ -d "$dir" ] && cd $dir && echo "Entering into $dir and generating sorted MarkDuplicates swarm" && bash "$scriptdir"/process-sortmd.sh )
done
#
# Now build the BQSR Recal table swarmfile
#
for dir in ${basedir[@]};
        do
                ( [ -d "$dir" ] && cd $dir && echo "Entering into $dir and generating BQSR tables swarm" && bash "$scriptdir"/process-BQSR.sh )
done
#
# Gather BQSR Reports
#
for dir in ${basedir[@]};
        do
                ( [ -d "$dir" ] && cd $dir && echo "Entering into $dir and gathering BQSR tables to write to swarm" && bash "$scriptdir"/process-BQSRreportsGather.sh )
done
#
# Apply BaseRecalibration by chromosome
#
for dir in ${basedir[@]};
        do
                ( [ -d "$dir" ] && cd $dir && echo "Entering into $dir and applying BQSR tables to write to swarm" && bash "$scriptdir"/process-applybqsr.sh )
done
#
# Gather BQSR bam files into a singular BQSR bam file
#
for dir in ${basedir[@]};
        do
                ( [ -d "$dir" ] && cd $dir && echo "Entering into $dir and gathering BQSR bams to merge" && bash "$scriptdir"/process-gatherBQSRbams.sh )
done
#
# Perform HaplotypeCaller per chromosome
for dir in ${basedir[@]};
        do
                ( [ -d "$dir" ] && cd $dir && echo "Entering into $dir and performing HaplotypeCaller per chromsome to write to swarm" && bash "$scriptdir"/process-haplotypecaller.sh )
done
#
## Pre-Release v2.0-alpha complete
