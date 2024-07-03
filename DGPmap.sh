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
                                -i|--input)     export FILE_DIR=$2; shift;;
                                -f|--fastq)     export FQ_IN=$2; shift;;
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
                        echo "Skipping "$dir""; exit 1
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
                        echo $wcf
                                if [ -e "$file" ]
                                then
                                        wcb=$(find "$dir"/bwa_bam/ -name "sort_*.bam" | wc -l); echo $wcb;
                                else
                                        wcb=0; echo $wcb
                                fi
                        if ! [ $((wcf/2)) == "$wcb" ]
                        then
                                echo "Making alignment swarm for "$dir""; bash "$scriptdir"/process-BWA.sh
                        else
                                echo "Skipping "$dir""; exit 1
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
                        echo $wcs
                                if [ -e "$file" ]
                                then
                                        wcd=$(find "$dir"/dedup_bam/ -name "*.sort.md.bam" | wc -l); echo $wcd;
                                else
                                        wcd=0; echo $wcd
                                fi
                                if ! [ -e $file ] && [ $((wcf/2)) == "$wcs" ]
                        then
                                echo "Making sorted mark duplicates swarm for "$dir""; bash "$scriptdir"/process-sortmd.sh
                        else
                                echo "Skipping "$dir""; exit 1
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
                                echo $wcr # debug line
                                        if [ -e $file ] && [ $wcr == 40 ]
                                        then
                                                echo "Skipping ""$dir"""; exit 1
                                        else
                                                echo "Making BQSR tables for ""$dir"""; bash "$scriptdir"/process-BQSR.sh
                                        fi
                done )
        done
}
#
# Gather BQSR Reports
#
bqsrgather(){
        for dir in ${basedir[@]};
        do
                cd "$dir";
                ( for file in "$dir"/BQSR/tables/*.bqsrgathered.reports.list
                        do
                                wcr=$(find "$dir"/BQSR/tables -name "*_recal.table" | wc -l);
                                echo $wcr # debug line
                                if [ -e $file ] && [ $wcr == 40 ]
                                then
                                        echo "Skipping ""$dir"""; exit 1
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
                                echo $wcb # debug line
                                if [ -e $file ] && [ $wcb == 40 ]
                                then
                                        echo "Skipping ""$dir"""; exit 1
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
                                echo $wcb # debug line
                                if [ -e $file ] && [ $wcb == 40 ]
                                then
                                        echo "Skipping ""$dir"""; exit 1
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
                                echo $wch # debug line
                                if [ -e $file ] && [ $wch == 40 ]
                                then
                                        echo "Skipping ""$dir"""; exit 1
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
                                echo $wch # debug line
                                if [ -e $file ] && [ $wch == 40 ]
                                then
                                        echo "Skipping ""$dir"""; exit 1
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
