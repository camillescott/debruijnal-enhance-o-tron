import pytest

from debruijnal_enhance_o_tron.sequence import (mutate_base,
                                                 mutate_sequence,
                                                 mutate_position,
                                                 reads,
                                                 kmers)

from debruijnal_enhance_o_tron.fixtures.sequence import (using_ksize,
                                                         using_length)

def test_ksize_default(ksize):
    assert ksize == 21


@using_ksize(31)
def test_ksize_override(ksize):
    assert ksize == 31


@using_ksize([25,31])
def test_ksize_override_seq(ksize):
    assert ksize in [25,31]
    if ksize != 25:
        assert ksize == 31
    elif ksize != 31:
        assert ksize == 25
    else:
        assert False


def test_length_default(length):
    assert length in [500, 1000]
    if length != 500:
        assert length == 1000
    elif length != 1000:
        assert length == 500
    else:
        assert False


@using_length(777)
def test_length_override(length):
    assert length == 777


@using_ksize(25)
@using_length(100)
def test_using_compose(length, ksize):
    assert length == 100
    assert ksize == 25


@using_ksize(5)
@using_length(20)
def test_random_sequence(random_sequence, ksize, length):
    seq1 = random_sequence()
    seq2 = random_sequence(exclude=seq1)

    assert len(seq1) == 20
    assert len(seq2) == 20

    for kmer in kmers(seq1, ksize):
        assert kmer not in seq2


def test_graph_adapter_noimpl(graph, ksize):
    with pytest.raises(NotImplementedError):
        graph.get(None)

    with pytest.raises(NotImplementedError):
        graph.add(None)

    with pytest.raises(NotImplementedError):
        graph.left_degree(None)

    with pytest.raises(NotImplementedError):
        graph.right_degree(None)

    with pytest.raises(NotImplementedError):
        graph.degree(None)


def test_linear_path_noimpl(linear_path, graph, ksize):

    with pytest.raises(NotImplementedError):
        # calling the inner function should trip this
        # on graph.add
        _graph, sequence = linear_path()



