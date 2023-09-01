#!/usr/bin/env python

import pysam
import os
from argparse import ArgumentParser
from argparse import ArgumentDefaultsHelpFormatter


def filter_and_merge_bam(bam_file, mod_bam, out_dir, out_prefix, primary=False):
    mod_read_ids = fetch_mod_read_ids(mod_bam, primary, out_dir)

    print('start filtering...')
    bam = pysam.AlignmentFile(bam_file, 'rb')
    filt_bam = pysam.AlignmentFile(f'{out_dir}/filtered.bam', 'wb', template=bam)
    for aln in bam:
        if primary and throw_aln(aln):
            continue
        if aln.query_name not in mod_read_ids:
            filt_bam.write(aln)
    bam.close()
    filt_bam.close()

    print('merge and sort...')
    mod_bam = f'{out_dir}/mod_primary.bam' if primary else mod_bam
    bams2merge = [f'{out_dir}/filtered.bam', mod_bam]
    pysam.merge("-f", f'{out_dir}/merged.bam', *bams2merge)

    final_out = f'{out_dir}/{out_prefix}.sorted.merged.bam'
    pysam.sort("-o", final_out, f'{out_dir}/merged.bam', "-@", "10")
    pysam.index(final_out)

    os.remove(f'{out_dir}/merged.bam')
    os.remove(f'{out_dir}/filtered.bam')

def fetch_mod_read_ids(mod_bam, primary, out_dir):
    read_ids = set()
    bam = pysam.AlignmentFile(mod_bam, 'rb')
    if primary:
        assert out_dir is not None, 'specify out_dir'
        out = f'{out_dir}/mod_primary.bam'
        out_bam = pysam.AlignmentFile(out, 'wb', template=bam)

    for aln in bam:
        read_ids.add(aln.query_name)
        if primary and not throw_aln(aln):
            out_bam.write(aln)

    bam.close()
    if primary: out_bam.close()
    return read_ids

def throw_aln(aln):
    # can add more filters later on
    if aln.is_secondary or aln.is_supplementary:
        return True
    return False

def main(args):
    os.makedirs(args.out_dir, exist_ok=True)
    filter_and_merge_bam(args.bam, 
                         args.mod_bam, 
                         args.out_dir, 
                         args.prefix,
                         args.primary)

if __name__ == '__main__':
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter,
        add_help=False
    )
    parser.add_argument("-b", "--bam", help="unmodified bam file")
    parser.add_argument("-m", "--mod_bam", help="modified fastq file")
    parser.add_argument("--primary", action="store_true", help="only keep primary alignments for all bams")
    parser.add_argument("-o", "--out_dir", help="output directory")
    parser.add_argument("--prefix", help="prefix for output bam")
    args = parser.parse_args()
    main(args)

# run:
# python ~/scripts/filter_merge_bam.py -b chr22.HG002_hs37d5_ONT-UL_GIAB_20200122.phased.bam \
#   -m mod.chr22.bam -o . --prefix mod_chr22 --primary
