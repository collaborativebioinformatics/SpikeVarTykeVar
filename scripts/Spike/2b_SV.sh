#!/bin/bash

# Input vairables are: READ_LENGTH, VAF, MERGED_VCF, MOD_BAM and OUTPUT_DIR
# Optional input variables for short_paragraph mode:
# 1. REFERENCE : reference for MOD_BAM in fasta format
# 2. OUTPUT_PFX : output prefix for paragraph
# READ_LENGTH can be either SHORT or LONG or short_paragraph

READ_LENGTH=$1
VAF=$2
#$VCF_A=$3
#VCF_B=$4
MERGED_VCF=$3
MOD_BAM=$4
OUTPUT_VCF=$5
REFERENCE=$6
OUTPUT_PFX=$7

echo "MERGED_VCF: $MERGED_VCF"
echo "MOD_BAM: $MOD_BAM"

if [ "$READ_LENGTH" = "short" ]; then

    # https://github.com/Illumina/paragraph
    bin/idxdepth -b ${MOD_BAM} -r ${REFERENCE} -o ${OUTPUT_PFX}.depth.json
    printf  "id\tpath\tidxdepth\nmxsample\t${MOD_BAM}\t${OUTPUT_PFX}.depth.json\n" > ${OUTPUT_PFX}.sample_manifest.txt
    python3 bin/multigrmpy.py -i  $MERGED_VCF -m   ${OUTPUT_PFX}.sample_manifest.txt  -r $REFERENCE  -o ${OUTPUT_PFX}_out

elif [ "$READ_LENGTH" = "long" ]; then
    echo "Input is 'long' starting Sniffles2"
    sniffles --input $MOD_BAM --genotype-vcf $MERGED_VCF --sample-id mosaicSim --vcf "${OUTPUT_VCF}"

else
    echo "Error: Invalid input"
fi

