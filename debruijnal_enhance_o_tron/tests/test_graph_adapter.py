import pytest

from debruijnal_enhance_o_tron.fixtures import (using_ksize,
                                                using_length,
                                                GraphAdapter)
from debruijnal_enhance_o_tron.sequences import (kmers)

class BaseGraph(GraphAdapter):
    '''Super basic Graph implementation following the
    provided GraphAdapater interface.
    '''

    def __init__(self, ksize, *args, **kwargs):
        self.store = set()
        self.ksize = ksize

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

    def degree(self, item):
        d = 0
        for b in 'ACGT':
            if store.get(item[1:] + b):
                d += 1
            if store.get(b + item[:-1]):
                d += 1
        return d


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
