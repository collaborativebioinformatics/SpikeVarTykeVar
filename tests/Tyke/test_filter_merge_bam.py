import pytest
from scripts.Tyke.filter_merge_bam import *
import pysam

def create_test_alignment(is_secondary=False, is_supplementary=False):
    aln = pysam.AlignedSegment()
    aln.is_secondary = is_secondary
    aln.is_supplementary = is_supplementary
    return aln

def test_throw_away_aln():
    assert throw_away_aln(create_test_alignment()) == False
    assert throw_away_aln(create_test_alignment(is_secondary=True)) == True
    assert throw_away_aln(create_test_alignment(is_supplementary=True)) == True
