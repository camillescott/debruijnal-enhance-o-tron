import pytest

from debruijnal_enhance_o_tron.fixtures.sequence import using_ksize
from debruijnal_enhance_o_tron.sequence import (mutate_base,
                                                mutate_sequence,
                                                mutate_position,
                                                get_random_sequence,
                                                reads,
                                                kmers)


@using_ksize([5, 7])
def test_kmers(ksize):
    S = 'A' * ksize + 'T'
    res = list(kmers(S, ksize))
    assert res[0] == 'A' * ksize
    assert res[-1] == ('A' * (ksize - 1)) + 'T'


def test_mutate_sequence():
    for _ in range(100):
        assert 'A' not in mutate_sequence('A' * 10, 10)
        assert 'T' not in mutate_sequence('T' * 10, 10)
        assert 'C' not in mutate_sequence('C' * 10, 10)
        assert 'G' not in mutate_sequence('G' * 10, 10)


def test_mutate_position():
    assert mutate_position('AAAA', 2) in ['AACA', 'AAGA']
    assert mutate_position('TTTT', 2) in ['TTCT', 'TTGT']
    assert mutate_position('CCCC', 2) in ['CCAC', 'CCTC']
    assert mutate_position('GGGG', 2) in ['GGAG', 'GGTG']


def test_reads(ksize):
    contig = get_random_sequence(1000, ksize)

    for read in reads(contig, ksize):
        assert read in contig

    for read in reads(contig, ksize):
        assert mutate_sequence(read) not in contig


