import pytest

from debruijnal_enhance_o_tron.fixtures import (using_ksize,
                                                using_length,
                                                GraphAdapter)
from debruijnal_enhance_o_tron.sequences import (kmers)

class BaseGraph(GraphAdapter):
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
def test_linear_structure_basegraph(linear_structure, graph, ksize, length):
    _, sequence = linear_structure()
    assert len(sequence) == length

    for kmer in kmers(sequence, ksize):
        assert graph.get(kmer)


@using_length(100)
def test_right_tip_basegraph(right_tip_structure, graph, ksize, length):
    _, sequence, tip, S = right_tip_structure()
    assert len(sequence) == length

    for kmer in kmers(sequence, ksize):
        assert graph.get(kmer)
