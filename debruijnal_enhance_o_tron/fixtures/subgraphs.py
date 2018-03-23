import pytest

from debruijnal_enhance_o_tron.sequence import (mutate_base,
                                                mutate_sequence,
                                                mutate_position,
                                                get_random_sequence,
                                                reads,
                                                kmers,
                                                revcomp)


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
def linear_path(request, graph, ksize, random_sequence):
    '''Simple linear path graph structure.

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

        return graph, (sequence,)

    return get


@pytest.fixture
def right_tip(request, graph, ksize, random_sequence):
    '''
    Sets up a graph structure like so:
                                 ([S+1:S+K]+B tip)
    sequence                   ↗
    [0]→o→o~~o→(L)→([S:S+K] HDN)→(R)→o→o→o~~o→[-1]

    Where S is the start position of the high degreen node (HDN).
    That is, it has a single branch at the Sth K-mer.

    HDN: S : S+K
    L:   S-1 : S-1+K
    R:   S+1 : S+1+K
    '''
    def get():
        sequence = random_sequence()
        S = len(sequence) // 2
        # right of the HDN
        R = S + 1
        # the branch kmer
        tip = mutate_position(sequence[R:R+ksize], -1)

        graph.add(sequence)
        graph.add(tip)

        # Check for false positive neighbors and mark as expected failure if found
        if count_decision_nodes(sequence, graph, ksize) != {(1,2): 1}:
            request.applymarker(pytest.mark.xfail)

        return graph, (sequence, tip), S

    return get


@pytest.fixture
def right_fork(request, ksize, right_tip, random_sequence):
    '''
    Sets up a graph structure like so:
                                               branch
                                 ([S+1:S+K]+B)→o~~o→o
    core_sequence               ↗
    [0]→o→o~~o→(L)→([S:S+K] HDN)→(R)→o→o→o~~o→[-1]

    Where S is the start position of the high degreen node (HDN)
    and B is the mutated base starting the branch.

    This is a tip with a longer fork length.
    '''

    def get():
        graph, (core_sequence, tip), S = right_tip()
        print('\nCore Len:', len(core_sequence))
        branch_sequence = random_sequence()
        print('Branch len:', len(branch_sequence))

        branch_sequence = tip + random_sequence()

        graph.add(core_sequence)
        graph.add(branch_sequence)

        # Check for false positive neighbors
        core_decision_nodes = count_decision_nodes(core_sequence, graph, ksize)

        # the core sequence should conain a decision node 
        # with ldegree of 1 and rdegre of 2
        if core_decision_nodes != {(1,2): 1}:
            request.applymarker(pytest.mark.xfail)

        return graph, (core_sequence, branch_sequence), S

    return get

