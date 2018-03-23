.. image:: https://travis-ci.org/camillescott/debruijnal-enhance-o-tron.svg?branch=master
    :target: https://travis-ci.org/camillescott/debruijnal-enhance-o-tron

the Debruijnal Enhance-O-Tron!
=========================

pytest fixtures for testing de Bruin Graphs. WIP.


Basic Usage
-----------

The only required setup is to set up an adapter for your de Bruijn Graph
implementation. This is done by creating a `graph` fixture that overrides
the enhance-o-tron's basic one. It should conform to (either implement
the interface or subclass) the provided `GraphAdapter` class.

Here's a minimal example::

    import pytest
    from debruijnal_enhance_o_tron.fixtures import GraphAdapter


    class BaseGraph(GraphAdapter):
        '''Super basic Graph implementation following the
        provided GraphAdapater interface.

        All your graph adapter needs (for now) is `get()`, `add()`, 
        and `degree()`.
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

The `graph` fixture will be injected into the enhance-o-trons subgraph
generator fixtures. Your tests can then be structured like...::

    def test_something_with_unitig(linear_structure, graph, ksize, length):
        _, sequence = liner_structure()
        shiny_assembler = Assembler(graph)

        assert shiny_assembler.assemble(sequence[:ksize]) == sequence


Both these snippets assume you've injected the fixtures into your
test session via your `conftest.py`, in which you could just stick::

    from debruijnal_enhance_o_tron.fixtures import *
