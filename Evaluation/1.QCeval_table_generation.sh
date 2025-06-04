#!/bin/bash
REGIONS="UU_Cfam_GSD_1.0_5MB_Auto_regions.txt"
VCF_IN=
WORKDIR=
homedir=$(pwd)
#
# This script generates the base data that will be the foundation to perform downstream cohort analyses on generated datasets. It is designed to function within a slurm/HPC environment. Minor modifications to this script can be made to run this locally.
# Word of warning, as datasets get larger the intermediate QC files generated will in turn get larger and require a larger footprint than the dataset. Approximate dataset size + 30% for space required to run this pipeline. Once the full pipeline is complete and final files are completed then the intermediate files can be deleted to save space.
#
# getopts string
opts="i:o:s:h"
#
# Gets the command name without path
cmd(){ echo `basename $0`; }
#
# Help command output
usage(){
        echo "\
        `cmd` [OPTION...]
        -i, --input; input VCF to QC.
        -o, --output; output directory of the QC files.
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
                                -i|--input) export VCF_IN=$2; shift;;
                                -o|--output)     export WORKDIR=$2; shift;;
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
#
## Loop section to verify that required flags have been passed.
if [[ -z $VCF_IN ]] && [[ -z $WORKDIR ]] && [[ -z $SWARM_NAME ]]
then
        echo "Missing required flags!"
        usage
        exit 1
elif [[ -z $VCF_IN ]] && [[ -z $WORKDIR ]]
then
        echo "Missing multiple flags!"
        echo "Try '`cmd` -h or `cmd` --help' for more information."
        exit 1
elif [[ -z $VCF_IN ]] && [[ -z $SWARM_NAME ]]
then
        echo "Missing multiple flags!"
        echo "Try '`cmd` -h or `cmd` --help' for more information."
        exit 1
elif [[ -z $WORKDIR ]] && [[ -z $SWARM_NAME ]]
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
elif [[ -z $SWARM_NAME ]]
then
        echo "Missing swarm name flag!"
        echo "Try '`cmd` -h or `cmd` --help' for more information."
        exit 1
else
        shift
fi
#
OUTNAME=$(echo $VCF_IN | awk -F '/' '{print $NF}' | sed 's/.vcf.gz//')
#
### DEBUG AREA ###
#
#echo "BUG CHECK"
#echo ""
#echo "REGIONS= " $REGIONS
#echo "VCF_IN= " $VCF_IN
#echo "WORKDIR= " $WORKDIR
#echo "SWARM_NAME= " $SWARM_NAME
#echo "OUTNAME= "$OUTNAME
#exit
###
#
cd swarmfiles/
swarmdir=$(pwd)
> collect_dist_GT.swarm
#
cd "$homedir"
#
#
module load bcftools
module load vcflib
#
echo "modules loaded:"
module list
#
read -sp "`echo -e 'If modules are properly loaded, press any key or Ctrl+C to abort \n\b'`" -n1 key
echo "Script running..."
#
bcftools query -l "$VCF_IN" > "$WORKDIR"/"$OUTNAME".samples.txt
for i in $(cat $WORKDIR/$OUTNAME.samples.txt); do for j in $(seq 0 1 99); do echo "$i,$j" >> $WORKDIR/$OUTNAME.count_tbl_AB_GQ_DP.txt ;done ; done
for i in $(cat $WORKDIR/$OUTNAME.samples.txt); do for j in $(echo -e "0/0\n0/1\n1/1\n./."); do echo "$i,$j" >> $WORKDIR/$OUTNAME.count_tbl_GT.txt ;done ; done
printf "done"
#
for REG in $(cat $REGIONS); do

echo -e "bcftools view -f PASS -r $REG $VCF_IN | vcf2tsv -n NA -g /dev/stdin > /lscratch/\$SLURM_JOB_ID/$OUTNAME.$REG.GT.txt ; \
SMCOL=\$(head -n1 /lscratch/\$SLURM_JOB_ID/$OUTNAME.$REG.GT.txt | tr \"\\\t\" \"\\\n\" | grep -n \"^SAMPLE\$\" | cut -f1 -d\":\") ; \
GTCOL=\$(head -n1 /lscratch/\$SLURM_JOB_ID/$OUTNAME.$REG.GT.txt | tr \"\\\t\" \"\\\n\" | grep -n \"^GT\$\" | cut -f1 -d\":\") ; \
ADCOL=\$(head -n1 /lscratch/\$SLURM_JOB_ID/$OUTNAME.$REG.GT.txt | tr \"\\\t\" \"\\\n\" | grep -n \"^AD\$\" | cut -f1 -d\":\") ; \
DPCOL=\$(head -n1 /lscratch/\$SLURM_JOB_ID/$OUTNAME.$REG.GT.txt | tr \"\\\t\" \"\\\n\" | grep -n \"^DP\$\" | cut -f1 -d\":\" | tail -n1) ; \
GQCOL=\$(head -n1 /lscratch/\$SLURM_JOB_ID/$OUTNAME.$REG.GT.txt | tr \"\\\t\" \"\\\n\" | grep -n \"^GQ\$\" | cut -f1 -d\":\") ; \
cat /lscratch/\$SLURM_JOB_ID/$OUTNAME.$REG.GT.txt | awk -v smcol=\$SMCOL -v gtcol=\$GTCOL -v adcol=\$ADCOL -v dpcol=\$DPCOL -v gqcol=\$GQCOL '{print \$smcol\"\\\t\"\$gtcol\"\\\t\"\$adcol\"\\\t\"\$dpcol\"\\\t\"\$gqcol}' > $WORKDIR/$OUTNAME.$REG.GT_qual_data.txt ; \
cat $WORKDIR/$OUTNAME.$REG.GT_qual_data.txt | awk '\$2==\"0/1\"' | sed \"s/,/\\\t/\" | awk '{\$7 = int(\$3 / (\$5 +.001) * 100)}1' | cut -f1,7 -d\" \" | sed \"s/ /,/\" > $WORKDIR/$OUTNAME.$REG.AB_het_counts.txt ; \
awk 'NR==FNR{a[\$0]=0;next}\$0 in a{a[\$0]+=1}END{for(i in a)print i,a[i]}' $WORKDIR/$OUTNAME.count_tbl_AB_GQ_DP.txt $WORKDIR/$OUTNAME.$REG.AB_het_counts.txt | sort -V > $WORKDIR/$OUTNAME.$REG.AB_het_counts_sum.txt ; \
rm $WORKDIR/$OUTNAME.$REG.AB_het_counts.txt ; \
cat $WORKDIR/$OUTNAME.$REG.GT_qual_data.txt | awk '\$2==\"0/0\"' | sed \"s/,/\\\t/\" | awk '{\$7 = int(\$3 / (\$5 +.001) * 100)}1' | cut -f1,7 -d\" \" | sed \"s/ /,/\" > $WORKDIR/$OUTNAME.$REG.AB_ref_counts.txt ; \
awk 'NR==FNR{a[\$0]=0;next}\$0 in a{a[\$0]+=1}END{for(i in a)print i,a[i]}' $WORKDIR/$OUTNAME.count_tbl_AB_GQ_DP.txt $WORKDIR/$OUTNAME.$REG.AB_ref_counts.txt | sort -V > $WORKDIR/$OUTNAME.$REG.AB_ref_counts_sum.txt ; \
rm $WORKDIR/$OUTNAME.$REG.AB_ref_counts.txt ; \
cat $WORKDIR/$OUTNAME.$REG.GT_qual_data.txt | awk '\$2==\"1/1\"' | sed \"s/,/\\\t/\" | awk '{\$7 = int(\$3 / (\$5 +.001) * 100)}1' | cut -f1,7 -d\" \" | sed \"s/ /,/\" > $WORKDIR/$OUTNAME.$REG.AB_alt_counts.txt ; \
awk 'NR==FNR{a[\$0]=0;next}\$0 in a{a[\$0]+=1}END{for(i in a)print i,a[i]}' $WORKDIR/$OUTNAME.count_tbl_AB_GQ_DP.txt $WORKDIR/$OUTNAME.$REG.AB_alt_counts.txt | sort -V > $WORKDIR/$OUTNAME.$REG.AB_alt_counts_sum.txt ; \
rm $WORKDIR/$OUTNAME.$REG.AB_alt_counts.txt ; \
cat $WORKDIR/$OUTNAME.$REG.GT_qual_data.txt | awk '\$2==\"0/1\"' | cut -f1,5 | sed \"s/\\\t/,/\" > $WORKDIR/$OUTNAME.$REG.GQ_het_counts.txt ; \
awk 'NR==FNR{a[\$0]=0;next}\$0 in a{a[\$0]+=1}END{for(i in a)print i,a[i]}' $WORKDIR/$OUTNAME.count_tbl_AB_GQ_DP.txt $WORKDIR/$OUTNAME.$REG.GQ_het_counts.txt | sort -V > $WORKDIR/$OUTNAME.$REG.GQ_het_counts_sum.txt ; \
rm $WORKDIR/$OUTNAME.$REG.GQ_het_counts.txt ; \
cat $WORKDIR/$OUTNAME.$REG.GT_qual_data.txt | awk '\$2==\"0/0\"' | cut -f1,5 | sed \"s/\\\t/,/\" > $WORKDIR/$OUTNAME.$REG.GQ_ref_counts.txt ; \
awk 'NR==FNR{a[\$0]=0;next}\$0 in a{a[\$0]+=1}END{for(i in a)print i,a[i]}' $WORKDIR/$OUTNAME.count_tbl_AB_GQ_DP.txt $WORKDIR/$OUTNAME.$REG.GQ_ref_counts.txt | sort -V > $WORKDIR/$OUTNAME.$REG.GQ_ref_counts_sum.txt ; \
rm $WORKDIR/$OUTNAME.$REG.GQ_ref_counts.txt ; \
cat $WORKDIR/$OUTNAME.$REG.GT_qual_data.txt | awk '\$2==\"1/1\"' | cut -f1,5 | sed \"s/\\\t/,/\" > $WORKDIR/$OUTNAME.$REG.GQ_alt_counts.txt ; \
awk 'NR==FNR{a[\$0]=0;next}\$0 in a{a[\$0]+=1}END{for(i in a)print i,a[i]}' $WORKDIR/$OUTNAME.count_tbl_AB_GQ_DP.txt $WORKDIR/$OUTNAME.$REG.GQ_alt_counts.txt | sort -V > $WORKDIR/$OUTNAME.$REG.GQ_alt_counts_sum.txt ; \
rm $WORKDIR/$OUTNAME.$REG.GQ_alt_counts.txt ; \
cat $WORKDIR/$OUTNAME.$REG.GT_qual_data.txt | cut -f1,4 | awk -v n=2 -v maxval=99 -v repval=99 '{\$n = \$n > maxval ? repval : \$n; print}' | sed \"s/ /,/\" > $WORKDIR/$OUTNAME.$REG.DP_counts.txt ; \
awk 'NR==FNR{a[\$0]=0;next}\$0 in a{a[\$0]+=1}END{for(i in a)print i,a[i]}' $WORKDIR/$OUTNAME.count_tbl_AB_GQ_DP.txt $WORKDIR/$OUTNAME.$REG.DP_counts.txt | sort -V > $WORKDIR/$OUTNAME.$REG.DP_counts_sum.txt ; \
rm $WORKDIR/$OUTNAME.$REG.DP_counts.txt ; \
cat $WORKDIR/$OUTNAME.$REG.GT_qual_data.txt | cut -f1,2 | sed \"s/\\\t/,/\" > $WORKDIR/$OUTNAME.$REG.GT_counts.txt ; \
awk 'NR==FNR{a[\$0]=0;next}\$0 in a{a[\$0]+=1}END{for(i in a)print i,a[i]}' $WORKDIR/$OUTNAME.count_tbl_GT.txt  $WORKDIR/$OUTNAME.$REG.GT_counts.txt | sort -V > $WORKDIR/$OUTNAME.$REG.GT_counts_sum.txt ; \
rm $WORKDIR/$OUTNAME.$REG.GT_counts.txt" >> "$swarmdir"/collect_dist_GT.swarm

done
#
head -n 1 collect_dist_GT.swarm
#
while true; do
read -p "Does the swarmfile above have proper syntax? (Yes or No) " promptA
case $promptA in
        [YyEeSs]* ) break;;
        [NnOo]* ) echo "Troubleshoot scripts"; exit 1;;
        * ) echo "Re-enter answer: Yes or No";;
esac
done
#
jobid1=$(swarm -f collect_dist_GT.swarm -g 4 -t 4 --module bcftools,vcflib --time=4:00:00 --gres=lscratch:50 --logdir ~/job_outputs/QC/"$SWARM_NAME" --sbatch "--mail-type=ALL,TIME_LIMIT_80 --job-name "$SWARM_NAME"")
echo "QC Regions SwarmID: " $jobid1
