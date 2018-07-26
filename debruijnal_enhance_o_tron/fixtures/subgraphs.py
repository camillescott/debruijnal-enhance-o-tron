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
def consumer(request, graph):
    def consume(sequences):
        for sequence in sequences:
            graph.add(sequence)
    return graph, consume


def conditional_consume(request, *args):
    '''Check if the consumer fixture is active in this request,
    and if so, use it to consume the sequences into the graph.

    Returns False if consumer is not active OR the check_fp
    marker is not active on the test.
    '''

    if 'consumer' in request.fixturenames:
        graph, consume = request.getfixturevalue('consumer')
        consume(args)


def conditional_check_fp(request, *args):

    if request.node.get_marker('check_fp') is not None \
       and 'graph' in request.fixturenames:

        graph = request.getfixturevalue('graph')
        check_graph = graph.shallow_clone()

        for sequence in args:
            for kmer in kmers(sequence, request.getfixturevalue('ksize')):
                assert graph.get(kmer) == False
        for sequence in args:
            check_graph.add(sequence)
        return check_graph
    else:
        return False


def check_fp_xfail(request):
    pytest.xfail('False-positive check FAILED ({}).'.format(request.fixturename))


def check_fp_pass(request):
    print('False-positive check PASS ({}).'.format(request.fixturename))


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
        conditional_consume(request, sequence)
        graph = conditional_check_fp(request, sequence)
        if graph:
            if count_decision_nodes(sequence, graph, ksize):
                check_fp_xfail(request)
            else:
                check_fp_pass(request)

        return sequence

    return get


@pytest.fixture
def right_sea(request, ksize, random_sequence):
    '''Sets up a C shaped graph structure:
                ([1:1+K])→o~~o→[-1] top
              ↗
    ([:K] HDN)→ ([1:1+K])→o~~o→[-1] bottom

    That is, HDN has in-degree of 0 and out-degree of 2:
    the minimal-degree decision node.
    '''

    def get():
        top = random_sequence()
        bottom = random_sequence()
        hdn = random_sequence(length=ksize)
        top = hdn + top
        bottom = list(hdn + bottom)
        # make sure the HDN really is the HDN...
        bottom[ksize] = mutate_base(top[ksize])
        bottom = ''.join(bottom)

        conditional_consume(request, top, bottom)
        graph = conditional_check_fp(request, top, bottom)
        if graph:
            if count_decision_nodes(core,
                                    graph,
                                    ksize) != {(0,2): 1}:
                check_fp_xfail(request)
            else:
                check_fp_pass(request)

        return top, bottom

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
        S = (len(sequence) // 2) - (ksize // 2)
        # right of the HDN
        R = S + 1

        if S < 1:
            raise ValueError("ksize too large for length")

        # the branch kmer
        tip = mutate_position(sequence[R:R+ksize], -1)

        # Check for false positive neighbors and mark as expected failure if found
        conditional_consume(request, sequence, tip)
        graph = conditional_check_fp(request, sequence, tip)
        if graph:
            if count_decision_nodes(sequence,
                                    graph,
                                    ksize) != {(1,2): 1}:
                check_fp_xfail(request)
            else:
                check_fp_pass(request)

        return (sequence, tip), S

    return get


@pytest.fixture
def right_fork(request, ksize, length, right_tip, random_sequence):
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
        branch_sequence = random_sequence()

        branch_sequence = tip + random_sequence()[:length-S-ksize]

        conditional_consume(request, core_sequence, branch_sequence)
        graph = conditional_check_fp(request, core_sequence, branch_sequence)
        if graph:
            # Check for false positive neighbors
            core_decision_nodes = count_decision_nodes(core_sequence, graph, ksize)

            # the core sequence should conain a decision node 
            # with ldegree of 1 and rdegre of 2
            if core_decision_nodes != {(1,2): 1}:
                check_fp_xfail(request)
            else:
                check_fp_pass(request)

        return (core_sequence, branch_sequence), S

    return get


@pytest.fixture
def right_triple_fork(request, ksize, length, right_fork, random_sequence):
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
        bottom_branch = random_sequence()[:length-S-ksize]

        # the branch sequence, mutated at position S+1
        # choose a base not already represented at that position
        bases = {'A', 'C', 'G', 'T'}
        used = {core_sequence[S+ksize], top_branch[ksize-1]}
        mutated = random.choice(list(bases - used))

        bottom_branch = top_branch[:ksize - 1] + mutated + bottom_branch

        
        conditional_consume(request, core_sequence, bottom_branch, top_branch)
        graph = conditional_check_fp(request, core_sequence, bottom_branch, top_branch)
        if graph:
            core_decision_nodes = count_decision_nodes(core_sequence,
                                                       graph,
                                                       ksize)

            if not (core_decision_nodes == {(1,3): 1}):
                check_fp_xfail(request)
            else:
                check_fp_pass(request)

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
        HDN_L = (len(wildtype_sequence) // 2) - ksize
        HDN_R = HDN_L + ksize + 1

        if HDN_L < 1:
            raise ValueError("ksize too long for length")

        snp_sequence = mutate_position(wildtype_sequence, HDN_L + ksize)

        conditional_consume(request, wildtype_sequence, snp_sequence)
        graph = conditional_check_fp(request, wildtype_sequence, snp_sequence)
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
                check_fp_xfail(request)
            else:
                check_fp_pass(request)

        return (wildtype_sequence, snp_sequence), HDN_L, HDN_R

    return get


@pytest.fixture
def tandem_quad_forks(request, ksize, length, linear_path, random_sequence):

    def get():
        core = linear_path()
        S_l = (len(core) // 2) - ksize
        S_r = S_l + 1

        left_branches = [kmer + random_sequence(exclude=kmer) \
                         for kmer in right_kmers(core[S_l:S_l+ksize]) \
                         if kmer not in core]
        right_branches = [kmer + random_sequence(exclude=kmer) \
                          for kmer in right_kmers(core[S_r:S_r+ksize]) \
                          if kmer not in core]

        conditional_consume(request, core, *left_branches, *right_branches)
        graph = conditional_check_fp(request, core, *left_branches, *right_branches)
        if graph:
            decision_nodes = count_decision_nodes(core, graph, ksize)
            if not decision_nodes == {(1,4): 2}:
                check_fp_xfail(request)
            else:
                check_fp_pass(request)

        return (core, left_branches, right_branches), S_l, S_r
    
    return get


@pytest.fixture(params=[2,6,10], ids=lambda r: 'repeats={0}'.format(r))
def tandem_repeats_lt_ksize(request, ksize):
    if ksize < 4:
        raise ValueError('Must use ksize >= 4')

    def get():
        repeat = _get_random_sequence(ksize - 2)
        tandem_repeats = repeat * request.param

        conditional_consume(request, tandem_repeats)
        graph = conditional_check_fp(request, tandem_repeats)
        if graph:
            if count_decision_nodes(tandem_repeats, graph, ksize):
                check_fp_xfail(request)
            else:
                check_fp_pass(request)


        return (repeat, tandem_repeats), request.param

    return get


@pytest.fixture(params=[2,6,10], ids=lambda r: 'repeats={0}'.format(r))
def tandem_repeats_gt_ksize(request, ksize):

    def get():
        repeat = _get_random_sequence(ksize * 2)
        tandem_repeats = repeat * request.param

        conditional_consume(request, tandem_repeats)
        graph = conditional_check_fp(request, tandem_repeats)
        if graph:
            if count_decision_nodes(tandem_repeats, graph, ksize):
                check_fp_xfail(request)
            else:
                check_fp_pass(request)

        return (repeat, tandem_repeats), request.param

    return get


@pytest.fixture
def circular(request, linear_path):

    def get():
        sequence = linear_path()
        sequence += sequence

        conditional_consume(request, sequence)
        graph = conditional_check_fp(request, sequence)
        if graph:
            if count_decision_nodes(sequence, graph, ksize):
                check_fp_xfail(request)
            else:
                check_fp_pass(request)

        return sequence

    return get
