#!/usr/bin/env python

import pysam
import os
from argparse import ArgumentParser
from argparse import ArgumentDefaultsHelpFormatter


def filter_and_merge_bam(bam_file, modified_bam, out_dir):
    print('start filtering...')
    modified_reads = fetch_modified_reads(modified_bam, out_dir)

    bam = pysam.AlignmentFile(bam_file, 'rb')
    filtered_bam = pysam.AlignmentFile(f'{out_dir}/filtered.bam', 'wb', template=bam)
    for aln in bam:
        if throw_away_aln(aln):
            continue
        if aln.query_name not in modified_reads:
            filtered_bam.write(aln)
    bam.close()
    filtered_bam.close()

    print('merge and sort...')
    bams_to_merge = [f'{out_dir}/filtered.bam', f'{out_dir}/modified_primary.bam']
    pysam.merge("-f", f'{out_dir}/merged.bam', *bams_to_merge)

    final_out = f'{out_dir}/tyke.bam'
    pysam.sort("-o", final_out, f'{out_dir}/merged.bam', "-@", "6")
    pysam.index(final_out)

    os.remove(f'{out_dir}/merged.bam')
    os.remove(f'{out_dir}/filtered.bam')
    os.remove(f'{out_dir}/modified_primary.bam')

def fetch_modified_reads(modified_bam, out_dir):
    read_ids = set()
    tmp_out = f'{out_dir}/modified_primary.bam'
    with pysam.AlignmentFile(modified_bam, 'rb') as bam:
        with pysam.AlignmentFile(tmp_out, 'wb', template=bam) as out_bam:
            for aln in bam:
                read_ids.add(aln.query_name)
                if not throw_away_aln(aln):
                    out_bam.write(aln)
    return read_ids

def throw_away_aln(aln):
    # can add more filters
    if aln.is_secondary or aln.is_supplementary:
        return True
    return False

def main(args):
    os.makedirs(args.out_dir, exist_ok=True)
    filter_and_merge_bam(args.bam, 
                         args.modified_bam, 
                         args.out_dir)

if __name__ == '__main__':
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter,
        add_help=False
    )
    parser.add_argument("-b", "--bam", help="unmodified bam file")
    parser.add_argument("-m", "--modified_bam", help="modified bam file")
    parser.add_argument("-o", "--out_dir", help="output directory", default=os.getcwd())
    args = parser.parse_args()
    main(args)

# run:
# python ~/scripts/filter_merge_bam.py -b chr22.HG002_hs37d5_ONT-UL_GIAB_20200122.phased.bam \
#   -m modified.chr22.bam 
