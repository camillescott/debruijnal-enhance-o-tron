import pytest

from debruijnal_enhance_o_tron.sequence import (mutate_base,
                                                 mutate_sequence,
                                                 mutate_position,
                                                 reads,
                                                 kmers)

from debruijnal_enhance_o_tron.fixtures.sequence import (using_ksize,
                                                         using_length)

def test_ksize_default(ksize):
    assert ksize in [21,51]
    if ksize != 21:
        assert ksize == 51
    elif ksize != 51:
        assert ksize == 21
    else:
        assert False


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


def test_consumer_noimpl(consumer):
    graph, consume = consumer
    with pytest.raises(NotImplementedError):
        graph.get(None)
    with pytest.raises(NotImplementedError):
        consume(('AAAAAAA'))


def test_linear_path(linear_path, ksize, length):
    sequence = linear_path()
    assert len(sequence) == length


def test_right_tip(right_tip, ksize, length):
    (core, tip), S = right_tip()
    assert len(core) == length
    assert len(tip) == ksize

    assert tip[:-1] == core[S+1:S+ksize]


def test_right_fork(right_fork, ksize, length):
    (core, branch), S = right_fork()
    
    assert branch[:ksize-1] == core[S+1:S+ksize]


def test_right_triple_fork(right_triple_fork, ksize, length):
    (core, top, bottom), S = right_triple_fork()

    assert top[:ksize-1] == core[S+1:S+ksize]
    assert bottom[:ksize-1] == core[S+1:S+ksize]


def test_snp_bubble(snp_bubble, ksize, length):
    (wildtype, snp), L, R = snp_bubble()
    
    assert wildtype[:L+ksize] == snp[:L+ksize]
    assert wildtype[R:] == snp[R:]
    assert wildtype[R-1] != snp[R-1]

    for kmer in kmers(wildtype[L+1:R], ksize):
        assert kmer not in snp


def test_tandem_repeat_lt_ksize(tandem_repeats_lt_ksize, ksize):
    (repeat, tandem_repeats), n_repeats = tandem_repeats_lt_ksize()

    assert tandem_repeats.count(repeat) == n_repeats


def test_tandem_repeat_gt_ksize(tandem_repeats_gt_ksize, ksize):
    (repeat, tandem_repeats), n_repeats = tandem_repeats_gt_ksize()

    assert tandem_repeats.count(repeat) == n_repeats


def test_circular(circular, ksize, length):
    sequence = circular()
    assert len(sequence) == 2 * length

