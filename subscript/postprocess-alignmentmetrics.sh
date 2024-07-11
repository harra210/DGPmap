#!/bin/bash
#
## This is an internal script for the post-processing of the DGPmap pipeline that will perform the alignment statistics for later collection.
#
## Variable definition section
#
FQ_DIR=$(pwd)
#
cd txt_files
TXT_DIR=$(pwd)
cd "$FQ_DIR"
#
IFS=,$'\n' read -d '' -r -a samplename < "$TXT_DIR"/rename.txt # Sets the samplename
#
echo "cd "$FQ_DIR"; gatk CollectAlignmentSummaryMetrics -R "$CanFam4_Ref" -I "${samplename[0]}".BQSR.cram -O "$TXT_DIR"/"${samplename[0]}".alignment_summary_metrics" >> "$homedir"/postprocess-alignmentmetrics.swarm
