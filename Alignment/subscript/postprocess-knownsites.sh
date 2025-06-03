#!/bin/bash
#
## This is an internal script for the post-processing of the DGPmap pipeline that will genotype gVCFs at known sites to be used to generate depth statistics later.
#
## Variable definition section
#
FQ_DIR=$(pwd)
#
cd txt_files
TXT_DIR=$(pwd)
cd "$FQ_DIR"
#
cd gVCF
gVCF_DIR=$(pwd)
cd "$FQ_DIR"
#
## Here we need to get the changed internal IDs.
#
cd "$TXT_DIR"
#
IFS=,$'\n' read -d '' -r -a samplename < rename.txt
#
echo "cd "$FQ_DIR"; gatk --java-options \"-Xmx6g -Xms6g\" GenotypeGVCFs -R "$CanFam4_Ref" -V "${samplename[0]}".g.vcf.gz -O "$TXT_DIR"/"${samplename[0]}".knownsites.vcf.gz --include-non-variant-sites --intervals "$knownsites"" >> "$homedir"/postprocess-knownsites.swarm
