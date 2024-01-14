import pytest
from unittest.mock import Mock
from scripts.Tyke import extract_read
from scripts.Tyke.extract_read import edit_read, load_bam_queries



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


def test_edit_reads_shoud_retun_none_if_cigartuples_is_empty():
    # Arrange & Act
    return_value = edit_read(
        ref_seq="",
        ref_start=0,
        read="",
        quals="", 
        cigartuples=[], 
        variant_pos=0,
        variant_length=0, 
        variant_op="",
        insert_seq="",
    )

    # Assert
    assert return_value == (None, None)

def test_load_bam_queries(monkeypatch, tmp_path):
    # Arrange
    row_mock = Mock(
        query_sequence="ATTTCGGG",
        query_qualities=[60],
        query_name="test", 
        reference_name="test",
        reference_start=0,
    )
    monkeypatch.setattr("pysam.AlignmentFile", Mock(return_value=[row_mock]))
    monkeypatch.setattr("scripts.Tyke.extract_read.edit_read",Mock(return_value=("ATGC", [60])))
    output_file_path = tmp_path / "output_file"

    # Act
    load_bam_queries(Mock(), {"test":"test"}, output_file_path)
    actual_output_file_content = open(output_file_path,'r').readlines()
    expected_output_file_content = [
        "@test|before\n",
        "ATTTCGGG\n",
        "+\n",
        "]\n",
        "@test|after|INS|1242|1\n",
        "ATGC\n",
        "+\n",
        "]\n",
     ]
    
    # Assert 
    assert actual_output_file_content[0:4] == expected_output_file_content[0:4]