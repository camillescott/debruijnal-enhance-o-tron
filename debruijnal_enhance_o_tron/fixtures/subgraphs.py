import random
import pytest

from debruijnal_enhance_o_tron.sequence import (mutate_base,
                                                mutate_sequence,
                                                mutate_position,
                                                _get_random_sequence,
                                                get_random_sequence,
                                                reads,
                                                kmers,
                                                left_kmers,
                                                right_kmers,
                                                revcomp)
from debruijnal_enhance_o_tron.fixtures.collectors import (consume_collector, check_fp_collector)


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

    def reset(self):
        raise NotImplementedError()

    def shallow_clone(self):
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
def linear_path(request, ksize, random_sequence, consume_collector, check_fp_collector):
    '''Simple linear path graph structure.

    sequence
    [0]→o→o~~o→o→[-1]
    '''
    def _linear_path():
        sequence = random_sequence()

        consume_collector(sequence)
        check_fp_collector((lambda G: count_decision_nodes(sequence, G, ksize), {}))

        return sequence

    return _linear_path


@pytest.fixture
def right_sea(request, ksize, random_sequence, consume_collector, check_fp_collector):
    '''Sets up a C shaped graph structure:
                ([1:1+K])→o~~o→[-1] top
              ↗
    ([:K] HDN)→ ([1:1+K])→o~~o→[-1] bottom

    That is, HDN has in-degree of 0 and out-degree of 2:
    the minimal-degree decision node.
    '''

    def _right_sea():
        top = random_sequence()
        bottom = random_sequence()
        hdn = random_sequence(length=ksize)
        top = hdn + top
        bottom = list(hdn + bottom)
        # make sure the HDN really is the HDN...
        bottom[ksize] = mutate_base(top[ksize])
        bottom = ''.join(bottom)

        consume_collector(top, bottom)
        check_fp_collector((lambda G : count_decision_nodes(core,
                                                             G,
                                                             ksize),
                             {(0,2): 1}),
                            (lambda G: count_decision_nodes(top,
                                                            G,
                                                            ksize),
                             {}),
                            (lambda G: count_decision_nodes(bottom,
                                                            G,
                                                            ksize),
                             {}))

        return top, bottom

    return _right_sea


@pytest.fixture
def right_tip(request, ksize, random_sequence, consume_collector, check_fp_collector):
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
    def _right_tip():
        sequence = random_sequence()
        S = (len(sequence) // 2) - (ksize // 2)
        # right of the HDN
        R = S + 1

        if S < 1:
            raise ValueError("ksize too large for length")

        # the branch kmer
        tip = mutate_position(sequence[R:R+ksize], -1)

        consume_collector(sequence, tip)
        check_fp_collector((lambda G: count_decision_nodes(sequence, G, ksize), {(1,2): 1}),
                            (lambda G: count_decision_nodes(tip, G, ksize), {}))

        return (sequence, tip), S

    return _right_tip


@pytest.fixture
def right_fork(request, ksize, length, right_tip, random_sequence,
               consume_collector, check_fp_collector):
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

    def _right_fork():
        (core_sequence, tip), S = right_tip()
        branch_sequence = random_sequence()

        branch_sequence = tip + random_sequence()[:length-S-ksize]

        consume_collector(core_sequence, branch_sequence)
        check_fp_collector((lambda G: count_decision_nodes(core_sequence,
                                                            G,
                                                            ksize),
                             {(1,2): 1}),
                            (lambda G: count_decision_nodes(branch_sequence,
                                                            G,
                                                            ksize),
                             {}))

        return (core_sequence, branch_sequence), S

    return _right_fork


@pytest.fixture
def left_fork(request, right_fork):

    def _left_fork():
        (core_sequence, branch), pos = right_fork()
        
        


@pytest.fixture
def right_triple_fork(request, ksize, length, right_fork, random_sequence,
                      consume_collector, check_fp_collector):
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
    
    def _right_triple_fork():
        (core_sequence, top_branch), S = right_fork()
        bottom_branch = random_sequence()[:length-S-ksize]

        # the branch sequence, mutated at position S+1
        # choose a base not already represented at that position
        bases = {'A', 'C', 'G', 'T'}
        used = {core_sequence[S+ksize], top_branch[ksize-1]}
        mutated = random.choice(list(bases - used))

        bottom_branch = top_branch[:ksize - 1] + mutated + bottom_branch

        consume_collector(core_sequence, bottom_branch, top_branch)
        check_fp_collector((lambda G: count_decision_nodes(core_sequence, G, ksize),
                             {(1,3): 1}),
                            (lambda G: count_decision_nodes(bottom_branch, G, ksize),
                             {}))
        
        return (core_sequence, top_branch, bottom_branch), S 

    return _right_triple_fork


@pytest.fixture
def snp_bubble(request, ksize, linear_path, consume_collector, check_fp_collector):
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

    def _snp_bubble():
        wildtype_sequence = linear_path()
        HDN_L = (len(wildtype_sequence) // 2) - ksize
        HDN_R = HDN_L + ksize + 1

        if HDN_L < 1:
            raise ValueError("ksize too long for length")

        snp_sequence = mutate_position(wildtype_sequence, HDN_L + ksize)

        consume_collector(wildtype_sequence, snp_sequence)
        check_fp_collector((lambda G: count_decision_nodes(wildtype_sequence, G, ksize),
                             {(1,2): 1, (2,1):1}),
                            (lambda G: count_decision_nodes(snp_sequence, G, ksize),
                             {(1,2): 1, (2,1):1}))

        return (wildtype_sequence, snp_sequence), HDN_L, HDN_R

    return _snp_bubble


@pytest.fixture
def tandem_quad_forks(request, ksize, length, linear_path, random_sequence,
                      consume_collector, check_fp_collector):

    def _tandem_quad_forks():
        core = linear_path()
        S_l = (len(core) // 2) - ksize
        S_r = S_l + 1

        left_branches = [kmer + random_sequence(exclude=kmer) \
                         for kmer in right_kmers(core[S_l:S_l+ksize]) \
                         if kmer not in core]
        right_branches = [kmer + random_sequence(exclude=kmer) \
                          for kmer in right_kmers(core[S_r:S_r+ksize]) \
                          if kmer not in core]

        consume_collector(core, *left_branches, *right_branches)
        check_fp_collector((lambda G: count_decision_nodes(core, G, ksize), {(1,4): 2}),
                            *[(lambda G: count_decision_nodes(branch, G, ksize), {}) \
                              for branch in left_branches + right_branches])

        return (core, left_branches, right_branches), S_l, S_r
    
    return _tandem_quad_forks


@pytest.fixture(params=[2,6,10], ids=lambda r: 'repeats={0}'.format(r))
def tandem_repeats_lt_ksize(request, ksize, consume_collector, check_fp_collector):
    if ksize < 4:
        raise ValueError('Must use ksize >= 4')

    def _tandem_repeats_lt_ksize():
        repeat = _get_random_sequence(ksize - 2)
        tandem_repeats = repeat * request.param

        consume_collector(tandem_repeats)
        check_fp_collector((lambda G: count_decision_nodes(tandem_repeats, G, ksize), {}))

        return (repeat, tandem_repeats), request.param

    return _tandem_repeats_lt_ksize


@pytest.fixture(params=[2,6,10], ids=lambda r: 'repeats={0}'.format(r))
def tandem_repeats_gt_ksize(request, ksize, consume_collector, check_fp_collector):

    def _tandem_repeats_gt_ksize():
        repeat = _get_random_sequence(ksize * 2)
        tandem_repeats = repeat * request.param

        consume_collector(tandem_repeats)
        check_fp_collector((lambda G: count_decision_nodes(tandem_repeats, G, ksize), {}))

        return (repeat, tandem_repeats), request.param

    return _tandem_repeats_gt_ksize


@pytest.fixture
def circular(request, linear_path, consume_collector, check_fp_collector):

    def _circular():
        sequence = linear_path()
        sequence += sequence

        consume_collector(sequence)
        check_fp_collector((lambda G: count_decision_nodes(sequence, G, ksize), {}))

        return sequence

    return _circular
