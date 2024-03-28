#!/bin/bash

# Input variables are: VARIANT, VAF, VCF_A, VCF_B, MOD_BAM, OUTPUT_DIR, READ_LENGTH
# VARIANT can be either SNV or SV
# READ_LENGTH can be either SHORT or LONG for SVs

# Example command:
# For SNV:
# ./2b_re-genotyping_main.sh SNV 0.1 /path/to/vcf1.vcf /path/to/vcf2.vcf /path/to/mod.bam /path/to/output_dir
# or 
# For SV:
# ./2b_re-genotyping_main.sh SV 0.1 /path/to/vcf1.vcf /path/to/vcf2.vcf /path/to/mod.bam /path/to/output_dir SHORT
# ./2b_re-genotyping_main.sh SV 0.1 /path/to/vcf1.vcf /path/to/vcf2.vcf /path/to/mod.bam /path/to/output_dir LONG

VARIANT=$1
VAF=$2
VCF_A=$3
VCF_B=$4
MOD_BAM=$5
OUTPUT_DIR=$6
READ_LENGTH=$7
MERGED_VCF="${OUTPUT_DIR}/merged.vcf.gz"
COLLAPSED_MERGED_VCF="${OUTPUT_DIR}/collapsed_merged.vcf.gz"
COLLAPSED_VARIANTS_VCF="${OUTPUT_DIR}/collapsed_variants.vcf.gz"
OUTPUT_VCF="${OUTPUT_DIR}/output_regenotyped.vcf.gz"
OUTPUT_VCF_FILTERED="${OUTPUT_DIR}/output_genotypes_filtered.vcf"
REFERENCE=${8}

mkdir -p $OUTPUT_DIR

# Merge the two VCF files with bcf tools
bcftools merge --force-samples -O z --write-index -o $MERGED_VCF $VCF_A $VCF_B
# Collapse variants with truvari
truvari collapse -i $MERGED_VCF -o $COLLAPSED_MERGED_VCF -c $COLLAPSED_VARIANTS_VCF

if [ "$VARIANT" = "SNV" ]; then
    echo "Input is 'SNV'"
    ./2b_SNV.sh $READ_LENGTH $VAF $COLLAPSED_MERGED_VCF $MOD_BAM $OUTPUT_VCF $REFERENCE
elif [ "$VARIANT" = "SV" ]; then
    echo "Input is 'SV'"
    ./2b_SV2.sh $READ_LENGTH $VAF $COLLAPSED_MERGED_VCF $OUTPUT_VCF "${OUTPUT_DIR}/output_genotypes" $REFERENCE "Spike_out"
else
    echo "Error: Invalid input"
fi

# Filter the output VCF file according to the VAF
python vaf_seperator.py $OUTPUT_VCF $VAF $OUTPUT_VCF_FILTERED
