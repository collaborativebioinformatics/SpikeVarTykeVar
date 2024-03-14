#!/bin/bash

# Input vairables are: READ_LENGTH, VAF, MERGED_VCF, MOD_BAM and OUTPUT_DIR
# Optional input variables for short_paragraph mode:
# 1. REFERENCE : reference for MOD_BAM in fasta format
# 2. OUTPUT_PFX : output prefix for paragraph
# READ_LENGTH can be either SHORT or LONG or short_paragraph

READ_LENGTH=$1
VAF=$2
VCF_A=$3
VCF_B=$4
MOD_BAM=$5
OUTPUT_VCF=$6
REFERENCE=$7
OUTPUT_PFX=$8

echo "VCF_A: $VCF_A"
echo "VCF_B: $VCF_B"
echo "MOD_BAM: $MOD_BAM"

if [ "$READ_LENGTH" = "short" ]; then

    # https://github.com/Illumina/paragraph
    bin/idxdepth -b ${MOD_BAM} -r ${REFERENCE} -o ${OUTPUT_PFX}.depth.json
    printf  "id\tpath\tidxdepth\nmxsample\t${MOD_BAM}\t${OUTPUT_PFX}.depth.json\n" > ${OUTPUT_PFX}.sample_manifest.txt
    python3 bin/multigrmpy.py -i  $MERGED_VCF -m   ${OUTPUT_PFX}.sample_manifest.txt  -r $REFERENCE  -o ${OUTPUT_PFX}_out

elif [ "$READ_LENGTH" = "long" ]; then
    echo "Input is 'long' starting Sniffles2"
    sniffles --input $MOD_BAM --genotype-vcf $VCF_A --sample-id mosaicSim --vcf "${OUTPUT_VCF}.vcf1.vcf"
    sniffles --input $MOD_BAM --genotype-vcf $VCF_B --sample-id mosaicSim --vcf "${OUTPUT_VCF}.vcf2.vcf"

    # Merge the two VCF files with bcftools
    bcftools merge --force-samples -O z --write-index -o $OUTPUT_VCF $VCF_A $VCF_B

else
    echo "Error: Invalid input"
fi

