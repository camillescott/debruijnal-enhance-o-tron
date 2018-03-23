from functools import wraps

import pytest

from debruijnal_enhance_o_tron.sequences import (mutate_base,
                                                 mutate_sequence,
                                                 mutate_position,
                                                 get_random_sequence,
                                                 reads,
                                                 kmers,
                                                 revcomp)

@pytest.fixture(params=[21], ids=['K=21'])
def ksize(request):
    '''
    Core fixture to set the K parameter.
    '''
    return request.param


def using_ksize(K):
    '''
    Convenience wrapper to reduce verbosity of indirect parametrization
    '''
    def wrapped(fixture_func):
        if isinstance(K, int):
            ksize = [K]
        else:
            ksize = list(K)
        return pytest.mark.parametrize('ksize', 
                                       ksize,
                                       indirect=['ksize'],
                                       ids=lambda k: 'K={0}'.format(k))(fixture_func)
    return wrapped


@pytest.fixture(params=[500, 1000], ids=lambda l: 'L={0}'.format(l))
def length(request):
    '''
    Basic length fixture.
    '''
    return request.param


def using_length(L):
    '''
    Convenience wrapper for basic length fixure.
    '''

    def wrapped(fixture_func):
        if isinstance(L, int):
            length = [L]
        else:
            length = list(L)
        return pytest.mark.parametrize('length', 
                                       length,
                                       indirect=['length'],
                                       ids=lambda l: 'L={0}'.format(l))(fixture_func)
    return wrapped


@pytest.fixture
def random_sequence(request, ksize, length):
    global_seen = set()

    def get(exclude=None):

        try:
            sequence = get_random_sequence(length, 
                                           ksize,
                                           exclude=exclude,
                                           seen=global_seen)
        except ValueError:
            request.applymarker(pytest.mark.xfail)

        for i in range(len(sequence)-ksize):
            global_seen.add(sequence[i:i+ksize-1])
            global_seen.add(revcomp(sequence[i:i+ksize-1]))

        return sequence

    return get


class GraphAdapter(object):

    def __init__(self, *args, **kwargs):
        pass

    def get(self, item):
        raise NotImplementedError()

    def add(self, item):
        raise NotImplementedError()

    def left_degree(self, item):
        raise NotImplementedError()

    def right_degree(self, item):
        raise NotImplementedError()

    def degree(self, item):
        return self.right_degree(item) + self.left_degree(item)


def count_decision_nodes(sequence, graph, ksize):
    '''Get the degree distribution of nodes with degree more than 2.
    '''

    dnodes = {}
    for kmer in kmers(sequence, ksize):
        d = (graph.left_degree(kmer), graph.right_degree(kmer))
        ld, rd = d
        if ld > 1 or rd > 1:
            dnodes[d] = dnodes.get(d, 0) + 1

    return dnodes


@pytest.fixture
def graph(ksize):
    '''Main point of client customization. Clients
    must override this fixture locally to expose their 
    dBG implementation. The implementation should, at minimum,
    support the methods from in the `GraphAdapter` interface.
    '''

    return GraphAdapter()


@pytest.fixture
def linear_structure(request, graph, ksize, random_sequence):
    '''Sets up a simple linear path graph structure.

    sequence
    [0]→o→o~~o→o→[-1]
    '''
    def get():
        sequence = random_sequence()
        graph.add(sequence)

        # Check for false positive neighbors in our graph
        # Mark as an expected failure if any are found
        if count_decision_nodes(sequence, graph, ksize):
            request.applymarker(pytest.mark.xfail)

        return graph, sequence

    return get
