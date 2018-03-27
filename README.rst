.. image:: https://travis-ci.org/camillescott/debruijnal-enhance-o-tron.svg?branch=master
    :target: https://travis-ci.org/camillescott/debruijnal-enhance-o-tron

the Debruijnal Enhance-O-Tron!
=========================

pytest fixtures for testing de Bruin Graphs. WIP.


Basic Setup
-----------

The only required setup is to set up an adapter for your de Bruijn Graph
implementation. This is done by creating a `graph` fixture that overrides
the enhance-o-tron's basic one. It should conform to (either implement
the interface or subclass) the provided `GraphAdapter` class.

Here's a minimal example:

.. code:: python

    import pytest
    from debruijnal_enhance_o_tron.fixtures import GraphAdapter


    class BaseGraph(GraphAdapter):
        '''Super basic Graph implementation following the
        provided `GraphAdapater` interface. The minimum requirements
        are to expose `add`, `get`, `left_degree`, and `right_degree`.
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

The `graph` fixture will be injected into the enhance-o-trons subgraph
generator fixtures. Your tests can then be structured like...:

.. code:: python

    def test_something_with_unitig(linear_path, graph, ksize, length):
        sequence = linear_path()
        shiny_assembler = ShinyAssembler(graph)

        assert shiny_assembler.assemble(sequence[:ksize]) == sequence


Both these snippets assume you've injected the fixtures into your
test session via your `conftest.py`, in which you could just stick:

.. code:: python

    from debruijnal_enhance_o_tron.fixtures import *


Using the Fixtures
------------------

You'll want to look at the documentation for each fixture in
`debruijnal_enhance_o_tron/fixtures/`, but for the most part
they return the constituent sequences and the position of
important structural decision nodes. For example, the `right_fork`
fixtures gives the core sequence, the branch sequence, and the left
index of the decision node within the core sequence:

.. code:: python

    def test_fork(right_fork, length, ksize):
        (core_sequence, fork_sequence), position = right_fork()
        #stuff

Note that the subgraph fixture return a *function* than can be called repeatedly
to generate more subgraphs.

By default, the sequences are not added to any dBG; however, if you've
overridden the `graph` fixture and provided an adapter as described above,
you can bring in the `consumer` fixture which will add add the sequences:

.. code:: python

    def test_consume_fork(right_fork, ksize, graph, consumer):
        (core_sequence, fork_sequence), position = right_fork()

        assert graph.get(core_sequence[:ksize])
