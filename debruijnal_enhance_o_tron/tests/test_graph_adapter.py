import pytest

from debruijnal_enhance_o_tron.fixtures import subgraphs
from debruijnal_enhance_o_tron.fixtures.sequence import (using_ksize,
                                                         using_length)
from debruijnal_enhance_o_tron.sequence import kmers


class BaseGraph(subgraphs.GraphAdapter):
    '''Super basic Graph implementation following the
    provided `GraphAdapater` interface.
    '''

    def __init__(self, ksize, *args, **kwargs):
        self.store = set()
        self.ksize = ksize
        super().__init__(*args, **kwargs)

    def get(self, item):
        return item in self.store

    def add(self, item):
        if len(item) < self.ksize:
            raise ValueError()
        elif len(item) == self.ksize:
            self.store.add(item)
        else:
            for kmer in kmers(item, self.ksize):
                self.store.add(kmer)

    def left_degree(self, item):
        return sum((self.get(b + item[:-1]) for b in 'ACGT'))

    def right_degree(self, item):
        return sum((self.get(item[1:] + b) for b in 'ACGT'))



@pytest.fixture
def graph(ksize):
    '''Test override of conftest-injected graph fixture.
    '''

    return BaseGraph(ksize)


@using_length(100)
def test_linear_path_basegraph(linear_path, graph, ksize, length):
    _, (sequence,) = linear_path()
    assert len(sequence) == length

    for kmer in kmers(sequence, ksize):
        assert graph.get(kmer)


@using_length(100)
def test_right_tip_basegraph(right_tip, graph, ksize, length):
    _, (sequence, tip), S = right_tip()
    assert len(sequence) == length

    for kmer in kmers(sequence, ksize):
        assert graph.get(kmer)


def test_right_fork_basegraph(right_fork, ksize, length):
    graph, (sequence, branch), S = right_fork()
