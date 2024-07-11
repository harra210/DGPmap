#!/bin/bash
#
## This is an internal script for the post-processing of the DGPmap pipeline that will perform the flagstat statistics for later collection.
#
## Variable definition section
#
FQ_DIR=$(pwd)
#
cd BQSR
BQSR_DIR=$(pwd)
cd "$FQ_DIR"
#
cd txt_files
TXT_DIR=$(pwd)
cd "$FQ_DIR"
#
IFS=,$'\n' read -d '' -r -a samplename < "$TXT_DIR"/rename.txt # Sets the samplename
#
echo "cd "$FQ_DIR"; samtools flagstat "${samplename[0]}".BQSR.cram > "$TXT_DIR"/"${samplename[0]}".flagstat" >> "$homedir"/postprocess-flagstat.swarm
