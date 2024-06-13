#!/bin/bash
#
# Function to create tables to input per folder prior to running the WGS alignment pipeline after being passed from parent script.
#
FQ_DIR=$(pwd)
#
cd $tmpdir
#
# Blanking tmp files in the tmp folder
> raw_files.tmp
> directories.tmp
> LB.tmp
> samplenames.tmp
> single_dir.tmp
> IDs.tmp
> R1.tmp
> R2.tmp
> PU.tmp
> SM.tmp
#
cd $FQ_DIR
#
mkdir -p txt_files
mkdir -p bwa_bam
mkdir -p dedup_bam
mkdir -p BQSR
cd BQSR/
mkdir -p tables
mkdir -p BQSR_chr
cd ../
mkdir -p gVCF
cd gVCF/
mkdir -p HC
cd ../
#
if [ -f *.txt ] ; then
	rm *.txt
fi
#
# The idea is that these tmp files will be pasted together to create a table
#
find $PWD -maxdepth 2 -name "*L00?_R1_001.fastq.gz" -printf '%h\n' &> "$tmpdir"/directories.tmp
find $PWD -maxdepth 2 -name "*L004_R1_001.fastq.gz" -printf '%h\n' &> "$tmpdir"/single_dir.tmp # Note: This will need to be modified to garner a unique result
#
cd $tmpdir
#
IFS=,$'\n' read -d '' -r -a singledir < single_dir.tmp
IFS=,$'\n' read -d '' -r -a multidir < directories.tmp
#
for ((i = 0; i < ${#singledir[@]}; i++))
do
       cd ${singledir[$i]}; find . -maxdepth 2 -name "*L00?_R1_001*" -printf '%f\n' | sort -V >> "$tmpdir"/raw_files.tmp
done
#
for ((i = 0; i < ${#singledir[@]}; i++))
do
	cd ${singledir[$i]}; find $PWD -maxdepth 2 -name "*L00?_R1_001*" | sort -V >> "$tmpdir"/R1.tmp
done
#
for ((i = 0; i < ${#singledir[@]}; i++))
do
	cd ${singledir[$i]}; find $PWD -maxdepth 2 -name "*L00?_R2_001*" | sort -V >> "$tmpdir"/R2.tmp
done
#
awk -F '_' '{print $2}' "$tmpdir"/raw_files.tmp > "$tmpdir"/LB.tmp
#
awk -F '_' '{print $2}' "$tmpdir"/raw_files.tmp > "$tmpdir"/IDs.tmp
#
cd $tmpdir
IFS=,$'\n' read -d '' -r -a R1 < R1.tmp
IFS=,$'\n' read -d '' -r -a R2 < R2.tmp
IFS=,$'\n' read -d '' -r -a LB < LB.tmp
IFS=,$'\n' read -d '' -r -a rawname < raw_files.tmp
#
for ((i = 0; i < ${#singledir[@]}; i++))
do
#Multiple SRA VERSION
#       cd ${singledir[$i]}; pwd | awk -F "/" '{print $8}' &> samplename.txt #IMPORTANT: THE AWK PRINT COLUMN NUMBER IS SUSCEPTIBLE TO CHANGING. VERIFY BEFORE CONTINUING
#NISC, uBAM, Single SRA VERSION
        cd ${singledir[$i]}; basename "`pwd`" &> samplename.txt
done
# Iterate through each directory with fastqs and head the first line of the raw R1 fastq and pipe that into file named raw_header.txt that is placed in the samples folder
#
for ((i = 0; i < ${#multidir[@]}; i++))
do
        cd "${multidir[$i]}"; zcat "${R1[$i]}" | head -n 1 >> "${multidir[$i]}"/raw_header.txt
done
#
# Iterates through each directory again, and then out of that print out just the Flowcell and Sample Tag information and label that as the trimmed_header
#
for ((i = 0; i < ${#singledir[@]}; i++))
do
	cd ${singledir[$i]}; awk -F ':' '{print $3,$4,$10}' raw_header.txt > trimmed_header.txt
done
#
for ((i = 0; i < ${#singledir[@]}; i++))
do
        cd ${singledir[$i]}; paste -d ' ' "$tmpdir"/IDs.tmp trimmed_header.txt > final_header.txt
done
#
for ((i = 0; i < ${#singledir[@]}; i++))
do
#Comment for Old Illumina
        cd ${singledir[$i]}; awk '{print$1"."$2"."$3"."$4}' final_header.txt > Final_ID.txt
done
#
for ((i = 0; i < ${#singledir[@]}; i++))
do
#Comment for Old Illumina
#For NISC samples
        cd ${singledir[$i]}; awk '{print $1"."$2}' trimmed_header.txt > final_PU.txt
done
#
for ((i = 0; i < ${#singledir[@]}; i++))
do
	cd ${singledir[$i]}; paste -d ' ' "$tmpdir"/LB.tmp "$tmpdir"/LB.tmp Final_ID.txt final_PU.txt "$tmpdir"/R1.tmp "$tmpdir"/R2.tmp >> "${LB[0]}".runs.table
done	
#
# Here the table is generated and now we perform cleanup
mv *.txt txt_files