import pytest

from scripts.Tyke import extract_read

def test_read_ref(tmpdir):
    read_id = "read_id"
    seq = "AAACCCTTTGGG"
    qual = [20] * len(seq)

    tmp_fq = tmpdir.join('tmp.fq')

    with open(tmp_fq, 'w') as fh:
        extract_read.write_fastx_record(fh, read_id, seq, qual)

    parsed_fq = extract_read.read_ref(tmp_fq, "fastq")
    parsed_seq = parsed_fq[read_id]
    assert(parsed_seq == seq)
