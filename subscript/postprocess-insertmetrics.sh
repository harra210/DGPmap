#!/bin/bash
#
## This is an internal script for the post-processing of the DGPmap pipeline that will perform the insert size statistics for later collection.
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
echo "cd "$FQ_DIR"; gatk CollectInsertSizeMetrics -R "$CanFam4_Ref" -I "${samplename[0]}".BQSR.cram -O "$TXT_DIR"/"${samplename[0]}".insert_size_metrics.txt -H "$TXT_DIR"/"${samplename[0]}".insert_size_metrics.hist.pdf" >> "$homedir"/postprocess-insertmetrics.swarm
