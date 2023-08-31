#!/usr/bin/env python
# -*- coding: utf-8 -*-
import pysam
import argparse
from vcf_line_parser import VCFLineSV
import math
import random

def argument_parser():
    parser = argparse.ArgumentParser(
        description="SV caller",
        prog="SV_caller",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "-v",
        "--vcf_file",
        dest="vcf_file",
        type=str,
        help="vcf_file",
        required=True,
    )
    parser.add_argument(
        "-b",
        "--bam_file",
        dest="bam_file",
        type=str,
        help="bam_file",
    )

    return parser.parse_args()


def vcf_number_variants(input_vcf_file,input_bam_file):
    """return the number of each variant seprately"""
    bamfile = pysam.AlignmentFile('test.bam', 'rb')
    #print(variant_reads)
    with open(input_vcf_file, "r") as f:
        lines = f.readlines()
        for line in lines[:]:
            if line[0] != "#":
                obj = VCFLineSV(line)
                if obj.ERROR:
                    continue
                chromosome = '1'
                start_position = int(obj.POS)
                end_position = int(obj.POS) + int(abs(obj.SVLEN))
                #print(start_position,end_position)
                variant_reads = []
                for read in bamfile.fetch(chromosome, start_position,end_position):
                    if read.query_sequence is not None:
                        variant_reads.append(read.query_sequence[:10])
                print(f"DR={obj.DR},AF={obj.AF},SVTYPE={obj.SVTYPE},CHROMOE={obj.CHROM},POS={obj.POS},SV_LENGTH={obj.SVLEN}")
                no_reads=max(1,math.ceil(obj.DR*obj.AF))
                if len(variant_reads) > no_reads:
                    selected_elements = random.sample(variant_reads, no_reads)
                # print(f"{obj.DR}AF={obj.AF},SVTYPE={obj.SVTYPE},CHROMOE={obj.CHROM},POS={obj.POS},SV_LENGTH={obj.SVLEN},SV={obj.ALT}")
                print(selected_elements)
    return
    bamfile.close()


def main():
    args=argument_parser()
    print(args.bam_file)
    vcf_number_variants('test.vcf','test.bam')
if __name__ == "__main__":
	main()
