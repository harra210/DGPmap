#!/bin/bash
#
REGIONS="UU_Cfam_GSD_1.0_5MB_Auto_regions.txt"
VCF_IN=
OUTNAME=
WORKDIR=
#
# getopts string
opts="i:o:h"
#
# Gets the command name without path
cmd(){ echo `basename $0`; }
#
# Help command output
usage(){
        echo "\
        `cmd` [OPTION...]
        -i, --input; input VCF to QC.
        -w, --workdir; Work directory of the QC files.
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
                                -i|--input) export VCF_IN=$2; shift;;
                                -w|--workdir)     export WORKDIR=$2; shift;;
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
## Loop section to verify that required flags have been passed.
if [[ -z $VCF_IN ]] && [[ -z $WORKDIR ]]
then
        echo "Missing multiple flags!"
        echo "Try '`cmd` -h or `cmd` --help' for more information."
        exit 1
elif [[ -z $VCF_IN ]]
then
        echo "Missing input flag!"
        echo "Try '`cmd` -h or `cmd` --help' for more information."
        exit 1
elif [[ -z $WORKDIR ]]
then
        echo "Missing output directory flag!"
        echo "Try '`cmd` -h or `cmd` --help' for more information."
        exit 1
else
        shift
fi
#
OUTNAME=$(echo $VCF_IN | awk -F '/' '{print $NF}' | sed 's/.vcf.gz//')
#
# DEV AREA
#echo "DEV MODE"
#echo ""
#echo "REGIONS: " $REGIONS
#echo "WORKDIR: " $WORKDIR
#echo "VCF_IN: " $VCF_IN
#echo "OUTNAME: " $OUTNAME
#echo ""
#echo "DEV MODE END"
#exit
#
#
echo -e "Metric.all.txt loop running...\n"
#
for METRIC in $(echo -e "AB_het\nAB_ref\nAB_alt\nGQ_het\nGQ_ref\nGQ_alt\nDP"); do
        if [[ -f $WORKDIR/$OUTNAME.$METRIC.all.txt ]]; then
                rm $WORKDIR/$OUTNAME.$METRIC.all.txt
        fi
        for i in $(cat $WORKDIR/$OUTNAME.count_tbl_AB_GQ_DP.txt); do echo "$i 0" >> $WORKDIR/$OUTNAME.$METRIC.all.txt ; done
        for REG in $(cat $REGIONS); do
                if [[ -f $WORKDIR/$OUTNAME.$REG.${METRIC}_counts_sum.txt ]]; then
                        paste -d " " $WORKDIR/$OUTNAME.$REG.${METRIC}_counts_sum.txt $WORKDIR/$OUTNAME.$METRIC.all.txt | awk '{print $1" "$2+$4}' > $WORKDIR/$OUTNAME.$METRIC.all.txt0
                        mv $WORKDIR/$OUTNAME.$METRIC.all.txt0 $WORKDIR/$OUTNAME.$METRIC.all.txt
                else
                        echo $WORKDIR/$OUTNAME.$REG.${METRIC}_counts_sum.txt not found
                fi
        done
done
#
echo -e "Loop complete!\n"
echo -e "GT loop running...\n"
#
METRIC=GT
for i in $(cat $WORKDIR/$OUTNAME.count_tbl_GT.txt); do echo "$i 0" >> $WORKDIR/$OUTNAME.$METRIC.all.txt ; done
for REG in $(cat $REGIONS); do
        if [[ -f $WORKDIR/$OUTNAME.$REG.${METRIC}_counts_sum.txt ]]; then
                paste -d " " $WORKDIR/$OUTNAME.$REG.${METRIC}_counts_sum.txt $WORKDIR/$OUTNAME.$METRIC.all.txt | awk '{print $1" "$2+$4}' > $WORKDIR/$OUTNAME.$METRIC.all.txt0
                mv $WORKDIR/$OUTNAME.$METRIC.all.txt0 $WORKDIR/$OUTNAME.$METRIC.all.txt
        else
                echo $WORKDIR/$OUTNAME.$REG.${METRIC}_counts_sum.txt not found
        fi
done
#
echo -e "Loop complete!\n"
echo -e "Regions loop running...\n"
#
for REG in $(cat $REGIONS); do
        cat $WORKDIR/$OUTNAME.$REG.GT_counts_sum.txt | awk -v reg=$REG '$3 = reg' | sed "s/,/\t/" | sed "s/ /\t/g" >> $WORKDIR/$OUTNAME.GT.sum.all.txt
done

for i in $(find $WORKDIR/*all*); do sed -i "s/,/\t/" $i ; sed -i "s/ /\t/" $i ; done
echo "Script complete!"
