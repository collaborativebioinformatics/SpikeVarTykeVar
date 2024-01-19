#!/usr/bin/env python
# -*- coding: utf-8 -*-
import pysam
import argparse
from vcf_line_parser import VCFLineSV
import math
import random

from extract_read import edit_read, read_ref, write_fastx_record

def argument_parser():
    parser = argparse.ArgumentParser(
        description="Tyke var simulator",
        prog="Tyke var simulator",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "-v",
        "--vcf_file",
        dest="vcf_file",
        type=str,
        help="vcf file from vcf generator",
        required=True,
    )
    parser.add_argument(
        "-b",
        "--bam_file",
        dest="bam_file",
        type=str,
        help="bam file",
    )
    parser.add_argument(
        "-r",
        "--ref_file",
        type=str,
        help="Reference genome file for BAM",
    )
    parser.add_argument(
        "-o",
        "--out_file",
        type=str,
        help="Output fastq file",
    )

    return parser.parse_args()

def vcf_number_variants(input_vcf_file, input_bam_file, refs, outfile):
    processed_read_ids = set()
    """return the number of each variant seprately"""
    bamfile = pysam.AlignmentFile(input_bam_file, 'rb')
    vcffile = pysam.VariantFile(input_vcf_file, 'r') # SV or SNV vcf
    
    with open(outfile, 'w') as out_fh:
        for v in vcffile:
            var_type = v.info.get('SVTYPE') if 'SVTYPE' in v.info else 'SNV'
            ## does pysam return errors? maybe try & except ##
            # if obj.ERROR:
            #     continue
            chromosome = v.chrom
            start_position = v.start
            variant_seq = v.alts[0]
            variant_len = len(variant_seq)
            end_position = start_position + variant_len
            #print(start_position,end_position)
            #svtype = "INS"
            #variant_len = 100
            #variant_seq = 'T' * variant_len
            #print(f"{chromosome}:{start_position} {variant_len} {variant_seq}")
            variant_reads = []
            for read in bamfile.fetch(chromosome, start_position, end_position):
                if read.query_sequence is not None and read.flag == 0 and read.query_name not in processed_read_ids:
                    variant_reads.append(read)
            no_reads = max(1, math.ceil(len(variant_reads) * float(v.info.get('AF'))))
            if no_reads > len(variant_reads):
                #print(f"Coverage {len(variant_reads)} is below required minimum reads for variant {chromosome}:{start_position} AF {obj.AF}")
                continue
            random.shuffle(variant_reads)
            #print(f"DR={obj.DR},AF={obj.AF},SVTYPE={obj.SVTYPE},CHROMOE={obj.CHROM},POS={obj.POS},SV_LENGTH={obj.SVLEN}")
            #print(f"{obj.DR}AF={obj.AF},SVTYPE={obj.SVTYPE},CHROMOE={obj.CHROM},POS={obj.POS},SV_LENGTH={obj.SVLEN},SV={obj.ALT}")
            ref_seq = refs[read.reference_name]
            num_reads_edited = 0
            for read in variant_reads:
                read_edited = False
                if num_reads_edited < no_reads:
                    new_seq, new_qual = edit_read(ref_seq, read.reference_start, read.query_sequence, read.query_qualities, 
                                                  read.cigartuples, start_position, variant_len, var_type, variant_seq)
                    if new_seq:
                        #write_fastx_record(out_fh, read.query_name, read.query_sequence, read.query_qualities)
                        #query_name = f"{read.query_name}|{read.reference_name}|{start_position}|{svtype}|{variant_len}"
                        query_name = read.query_name
                        write_fastx_record(out_fh, query_name, new_seq, new_qual)
                        num_reads_edited += 1
                        read_edited = True
                        processed_read_ids.add(read.query_name)
                #if not read_edited:
                #    write_fastx_record(out_fh, read.query_name, read.query_sequence, read.query_qualities)
            if num_reads_edited != no_reads:
                print("Could not fulfill required minimum edited reads for variant at {chromosome}:{start_position}:{end_position} type {svtype}")

    bamfile.close()
    vcffile.close()
    return


def main():
    args=argument_parser()
    print(args.bam_file)
    refs = read_ref(args.ref_file)
    vcf_number_variants(args.vcf_file, args.bam_file, refs, args.out_file)
    out_file = args.out_file
if __name__ == "__main__":
	main()
