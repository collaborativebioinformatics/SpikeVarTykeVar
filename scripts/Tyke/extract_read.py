import pysam
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse

DEBUG = False

def edit_read(ref_seq, ref_start, read, cigartuples, variant_pos, variant_length, variant_op, insert_seq):
    #if DEBUG:
    if True:
        print(f"----ref start {ref_start} read len {len(read)}")
        print(f"----variant pos {variant_pos} var length {variant_length} var op {variant_op} var seq {insert_seq}")
        print(f"----cigar tuples {cigartuples}")
    ref_pos = ref_start # ref pos, variant pos are 1 indexed
    last_ref_pos = ref_pos
    read_var_pos = -1
    new_seq = ""
    read_pos = 0 # read pos is 0 indexed
    last_read_pos = 0
    for (c, l) in cigartuples:
        if c in (0, 7, 8):
            if DEBUG:
                new_seq += read[read_pos:read_pos + l]
            read_pos += l
            ref_pos += l
        elif c == 1:
            if DEBUG:
                new_seq += read[read_pos: read_pos + l]
            read_pos += l
        elif c == 2:
            ref_pos += l
        elif c == 3:
            ref_pos += l
        elif c == 4:
            if DEBUG:
                new_seq += read[read_pos: read_pos + l]
            read_pos += l
        elif c == 5 or c == 6:
            pass
        else:
            raise RuntimeError("Unknown cigr ops found " + c)
        if ref_pos > variant_pos and read_var_pos == -1:
            if DEBUG:
                print(f"Current ref pos {ref_pos} last ref pos {last_ref_pos} variant pos {variant_pos}")
            if c == 2 or c == 3:
                # if there's an del variant here, then we remove it from the following seq
                # if there's an insertion or snv do nothing
                if variant_op == "DEL":
                    if variant_length > (ref_pos - variant_pos):
                        variant_length -= (ref_pos - variant_pos)
                        variant_pos = ref_pos
            else:
                read_var_pos = last_read_pos + (variant_pos - last_ref_pos) - 1
        last_ref_pos = ref_pos
        last_read_pos = read_pos
    if DEBUG and new_seq != read:
        print(f"Something went wrong! Constructed read {new_seq} doesn't match query {read}")
    if read_var_pos != -1:
        if variant_op == "INS":
            new_seq = read[:read_var_pos + 1] + insert_seq + read[read_var_pos + 1:]
        elif variant_op == "SNV":
            new_seq = read[:read_var_pos] + insert_seq + read[read_var_pos+1:]
        else:
            new_seq = read[:read_var_pos] + ( read[read_var_pos + variant_length:] if (read_var_pos + variant_length) < len(read) else "")
        if DEBUG:
            print(f"Insert variant at {read_var_pos}")
            print(f"Bef - {read}")
            print(f"Aft - {new_seq}")
            print(f"Len diff {len(new_seq) - len(read)}")
        rec_before = SeqRecord(Seq(read), id="before|", description="before edit")
        rec_after = SeqRecord(Seq(new_seq), id="after|", description="after edit")
        recs = [rec_before, rec_after]
        SeqIO.write(recs, "my_example.fa", "fasta")
    else:
        if DEBUG:
            print("No variant location found. Does the last cigar entry contain insertions to the ref or soft clips?")

def read_ref(fq_file):
    refs = {}
    for seq_record in SeqIO.parse(fq_file, "fasta"):
        refs[seq_record.id] = repr(seq_record.seq)
    return refs

def load_bam_queries(bam_file, refs):
    import random
    bam = pysam.AlignmentFile(bam_file, "r")
    for i, row in enumerate(bam):
        if i == 15:
            if not row.query_sequence:
                print("No query sequence for ", row.query_name)
                continue
            if (row.reference_name and row.query_sequence):
                variant_pos = row.reference_start + random.randint(0, len(row.query_sequence) - 1)
                ops = ["INS", "DEL", "SNV"]
                variant_op = ops[random.randint(0, len(ops) - 1)]
                variant_length = 1 if variant_op == "SNV" else random.randint(1, 10000)
                bases = ['A', 'T', 'C', 'G']
                var_seq = "".join([bases[random.randint(0, 3)] for _ in range(variant_length)])
                edit_read(refs[row.reference_name], row.reference_start, row.query_sequence, row.cigartuples, variant_pos, variant_length, variant_op, var_seq)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
                    prog='Read editor',
                    description='edit reeads based on a variant and a BAM file')
    parser.add_argument('--bam', help='bam file')
    parser.add_argument('--ref', help='ref file')
    args = parser.parse_args()
    refs = read_ref(args.ref)
    bam1_queries = load_bam_queries(args.bam, refs)
