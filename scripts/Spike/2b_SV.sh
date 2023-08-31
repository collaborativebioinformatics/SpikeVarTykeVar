#!/bin/bash

# Input vairables are: READ_LENGTH, VAF, MERGED_VCF, MOD_BAM and OUTPUT_DIR
# READ_LENGTH can be either SHORT or LONG

READ_LENGTH=$1
VAF=$2
MERGED_VCF=$3
MOD_BAM=$4
OUTPUT_VCF=$5

# Descide which READ_LENGTH was selected

if [ "$READ_LENGTH" = "short" ]; then
    echo "Input is 'short' starting svtyper"
    # change to conda env which uses python 2.7 called mosaicSim27
    conda activate mosaicSim27
    svtyper -B $MOD_BAM -i $MERGED_VCF -o $OUTPUT_VCF
elif [ "$READ_LENGTH" = "long" ]; then
    echo "Input is 'long' starting Sniffles2"
    # change to conda env which uses python 3.6 called mosaicSim
    conda activate mosaicSim
    sniffles --input sample.bam --genotype-vcf $MERGED_VCF --sample-id mosaicSim --vcf $OUTPUT_VCF
else
    echo "Error: Invalid input"
fi

