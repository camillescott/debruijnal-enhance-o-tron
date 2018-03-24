import random
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
def consumer(request, graph):
    def consume(sequences):
        for sequence in sequences:
            graph.add(sequence)
    return graph, consume


def do_consume(request, *args):
    if 'consumer' in request.fixturenames:
        graph, consume = request.getfixturevalue('consumer')
        consume(args)
        return graph
    return False


@pytest.fixture
def linear_path(request, ksize, random_sequence):
    '''Simple linear path graph structure.

    sequence
    [0]→o→o~~o→o→[-1]
    '''
    def get():
        sequence = random_sequence()

        # Check for false positive neighbors in our graph
        # Mark as an expected failure if any are found
        graph = do_consume(request, sequence)
        if graph and count_decision_nodes(sequence, graph, ksize):
            request.applymarker(pytest.mark.xfail)

        return sequence

    return get


@pytest.fixture
def right_tip(request, ksize, random_sequence):
    '''
    Sets up a graph structure like so:
                                 ([S+1:S+K]+B tip)
    sequence                   ↗
    [0]→o→o~~o→(L)→([S:S+K] HDN)→(R)→o→o→o~~o→[-1]

    Where S is the start position of the high degreen node (HDN).
    That is, it has a single branch at the Sth K-mer.

    HDN: S:S+K
    L:   S-1:S-1+K
    R:   S+1:S+1+K

    The mutated base itself is at S+K
    '''
    def get():
        sequence = random_sequence()
        S = len(sequence) // 2
        # right of the HDN
        R = S + 1
        # the branch kmer
        tip = mutate_position(sequence[R:R+ksize], -1)

        # Check for false positive neighbors and mark as expected failure if found
        graph = do_consume(request, sequence, tip)
        if graph and count_decision_nodes(sequence,
                                          graph,
                                          ksize) != {(1,2): 1}:
            request.applymarker(pytest.mark.xfail)

        return (sequence, tip), S

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
        (core_sequence, tip), S = right_tip()
        print('\nCore Len:', len(core_sequence))
        branch_sequence = random_sequence()
        print('Branch len:', len(branch_sequence))

        branch_sequence = tip + random_sequence()

        graph = do_consume(request, core_sequence, branch_sequence)
        if graph:
            # Check for false positive neighbors
            core_decision_nodes = count_decision_nodes(core_sequence, graph, ksize)

            # the core sequence should conain a decision node 
            # with ldegree of 1 and rdegre of 2
            if core_decision_nodes != {(1,2): 1}:
                request.applymarker(pytest.mark.xfail)

        return (core_sequence, branch_sequence), S

    return get


@pytest.fixture
def right_triple_fork(request, ksize, right_fork, random_sequence):
    '''
    Sets up a graph structure like so:

                                       top_branch
                                ([:S+1]+B)→o~~o→o
    core_sequence              ↗
    [0]→o→o~~o→(L)→([S:S+K] HDN)→(R)→o→o→o~~o→[-1]
                               ↘
                                ([:S+1]+B)→o~~o→o
                                     bottom_branch

    Where S is the start position of the high degreen node (HDN).
    '''
    
    def get():
        (core_sequence, top_branch), S = right_fork()
        bottom_branch = random_sequence()
        print(len(core_sequence), len(top_branch), len(bottom_branch))

        # the branch sequence, mutated at position S+1
        # choose a base not already represented at that position
        bases = {'A', 'C', 'G', 'T'}
        used = {core_sequence[S+ksize], top_branch[ksize-1]}
        mutated = random.choice(list(bases - used))

        bottom_branch = top_branch[:ksize - 1] + mutated + bottom_branch

        
        graph = do_consume(request, core_sequence, bottom_branch, top_branch)
        if graph:
            core_decision_nodes = count_decision_nodes(core_sequence,
                                                       graph,
                                                       ksize)

            if not (core_decision_nodes == {(1,3): 1}):
                request.applymarker(pytest.mark.xfail)

        return (core_sequence, top_branch, bottom_branch), S 

    return get



@pytest.fixture
def snp_bubble(request, ksize, linear_path):
    '''
    Sets up a graph structure resulting from a SNP (Single Nucleotide
    Polymorphism).

                (HDN_L[1:]+SNP)→o~~o→(SNP+)
              ↗                           ↘
    o~~(HDN_L)                             (HDN_R)~~o
              ↘                           ↗
                (HDN_L[1:]+W)→o~~o~~o→(W+)


    HDN_L: S:S+K
    HDN_R: S+K+1:S+2K+1

    Where S is the start position of HDN directly left of the SNP (HDN_L),
    SNP is the mutated base, and W is the wildtype (original) base.
    Of course, W and SNP could be interchanged here, we don't actually
    know which is which ;)

    Note our parametrization: we need a bit more room from the ends,
    so we bring the rightmost SNP a tad left.
    '''

    def get():
        wildtype_sequence = linear_path()
        HDN_L = len(wildtype_sequence) // 2
        HDN_R = HDN_L + ksize + 1

        snp_sequence = mutate_position(wildtype_sequence, HDN_L + ksize)

        graph = do_consume(request, wildtype_sequence, snp_sequence)
        if graph:
            wildtype_decision_nodes = count_decision_nodes(wildtype_sequence,
                                                           graph,
                                                           ksize)
            snp_decision_nodes = count_decision_nodes(snp_sequence,
                                                      graph,
                                                      ksize)
            if not (snp_decision_nodes == \
                    wildtype_decision_nodes == \
                    {(1,2): 1, (2,1):1}):
                request.applymarker(pytest.mark.xfail)

        return (wildtype_sequence, snp_sequence), HDN_L, HDN_R

    return get

