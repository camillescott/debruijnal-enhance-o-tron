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
from debruijnal_enhance_o_tron.fixtures.sequence import using_pivot


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
    for i, kmer in enumerate(kmers(sequence, ksize)):
        d = (graph.left_degree(kmer), graph.right_degree(kmer))
        ld, rd = d
        if ld > 1 or rd > 1:
            dnodes[d] = dnodes.get(d, 0) + 1
            print(i, d)

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
        check_fp_collector((lambda G : count_decision_nodes(top,
                                                            G,
                                                            ksize),
                             {(0,2): 1}),
                            (lambda G: count_decision_nodes(bottom,
                                                            G,
                                                            ksize),
                             {(0,2): 1}))

        return top, bottom

    return _right_sea


@pytest.fixture
def left_sea(request, ksize, random_sequence, consume_collector, check_fp_collector):

    def _left_sea():
        top = random_sequence()
        bottom = random_sequence()
        hdn = random_sequence(length=ksize)
        top = top + hdn
        bottom = list(bottom + hdn)
        # make sure the HDN really is the HDN...
        bottom[-(ksize + 1)] = mutate_base(top[-(ksize+1)])
        bottom = ''.join(bottom)

        consume_collector(top, bottom)
        check_fp_collector((lambda G : count_decision_nodes(top,
                                                            G,
                                                            ksize),
                             {(2,0): 1}),
                            (lambda G: count_decision_nodes(bottom,
                                                            G,
                                                            ksize),
                             {(2,0): 1}))

        return top, bottom

    return _left_sea


@pytest.fixture
def right_tip(request, ksize, internal_pivot, random_sequence,
              consume_collector, check_fp_collector):
    '''
    Sets up a graph structure like so:

                                    ([pivot+1:] tip)
        sequence                   ↗
        [0]→o→o~~o→(L)→([pivot:pivot+K] decision)→(R)→o→o→o~~o→[-1]

    Where pivot is the start position of the decision k-mer:
    pivot:pivot+K has right-degree of two.

    The mutated base is at S+K
    '''
    def _right_tip():
        sequence = random_sequence()
        pivot = internal_pivot

        if pivot < 1:
            raise ValueError("ksize too large for length")

        # the branch kmer
        tip = mutate_position(sequence[pivot+1:pivot+1+ksize], -1)

        consume_collector(sequence, tip)
        check_fp_collector((lambda G: count_decision_nodes(sequence, G, ksize), {(1,2): 1}),
                            (lambda G: count_decision_nodes(tip, G, ksize), {}))

        return (sequence, tip), pivot

    return _right_tip


@pytest.fixture
def right_fork(request, ksize, length, 
               right_tip, random_sequence,
               consume_collector, check_fp_collector):
    '''
    Sets up a graph structure like so:

                                                   branch
                                     ([pivot+1:])→o~~o→o
        core_sequence               ↗
        [0]→o→o~~o→(L)→([pivot:pivot+K] decision)→(R)→o→o→o~~o→[-1]

    Where pivot is the start position of the decision k-mer:
    pivot:pivot+K has right-degree of two.

    This is a tip with a longer fork length.
    '''

    def _right_fork():
        (core_sequence, tip), pivot = right_tip()
        branch_sequence = random_sequence()

        branch_sequence = tip + random_sequence()[:length-pivot-ksize]

        consume_collector(core_sequence, branch_sequence)
        check_fp_collector((lambda G: count_decision_nodes(core_sequence,
                                                           G,
                                                           ksize),
                            {(1,2): 1}),
                           (lambda G: count_decision_nodes(branch_sequence,
                                                           G,
                                                           ksize),
                            {}))

        return (core_sequence, branch_sequence), pivot

    return _right_fork


@pytest.fixture
def left_fork(request, ksize, length, right_fork,
              consume_collector, check_fp_collector):
    '''
    The opposite of a right-fork: pivot:pivot+K has left-degree of two;
    the mutated base in the core sequence is at pivot-1.

    This is a tip with a longer fork length.
    '''

    def _left_fork():
        (core_sequence, branch), pivot = right_fork()
        core_sequence = revcomp(core_sequence)
        branch = revcomp(branch)
        pivot = length - pivot - ksize
        
        _collector = consume_collector()
        _collector.pop() # remove previous two fixtures that compose right_fork
        _collector.pop()
        consume_collector(core_sequence, branch)
        _collector = check_fp_collector()
        _collector.pop()
        _collector.pop()
        check_fp_collector((lambda G: count_decision_nodes(core_sequence,
                                                           G,
                                                           ksize),
                            {(2,1): 1}),
                           (lambda G: count_decision_nodes(branch, G, ksize), {}))

        return (core_sequence, branch), pivot

    return _left_fork
        


@pytest.fixture
def right_triple_fork(request, ksize, length, right_fork, random_sequence,
                      consume_collector, check_fp_collector):
    '''
    Sets up a graph structure like so:

                                           top_branch
                                    ([pivot+1:])→o~~o→o
        core_sequence              ↗
        [0]→o→o~~o→(L)→([pivot:pivot+K] decision)→(R)→o→o→o~~o→[-1]
                                   ↘
                                    ([pivot+1:]+)→o~~o→o
                                         bottom_branch

    Where S is the start position of the high degreen node (HDN).
    '''
    
    def _right_triple_fork():
        (core_sequence, top_branch), pivot = right_fork()
        bottom_branch = random_sequence()[:length-pivot-ksize]

        # the branch sequence, mutated at position S+1
        # choose a base not already represented at that position
        bases = {'A', 'C', 'G', 'T'}
        used = {core_sequence[pivot+ksize], top_branch[ksize-1]}
        mutated = random.choice(list(bases - used))

        bottom_branch = top_branch[:ksize - 1] + mutated + bottom_branch

        consume_collector(core_sequence, bottom_branch, top_branch)
        check_fp_collector((lambda G: count_decision_nodes(core_sequence, G, ksize),
                             {(1,3): 1}),
                            (lambda G: count_decision_nodes(bottom_branch, G, ksize),
                             {}))
        
        return (core_sequence, top_branch, bottom_branch), pivot 

    return _right_triple_fork


@pytest.fixture
def snp_bubble(request, ksize, linear_path, middle_pivot,
               consume_collector, check_fp_collector):
    '''
    Sets up a graph structure resulting from a SNP (Single Nucleotide
    Polymorphism).

                        (SNP_L[1:]+SNP)→o~~o→(SNP+)
                      ↗                            ↘
        o~~(decision_L)                             (decision_R)~~o
                      ↘                            ↗
                        (SNP_L[1:]+W)→o~~o~~o→(W+)

    Returns two equal-length sequences, with a mutation at position
    pivot+K in the second sequence.

    decision_L is pivot:pivot+K in both sequences.
    decision_R is pivot+K+1:pivot+2K+1 in both sequences

    SNP is the mutated base, and W is the wildtype (original) base.
    Of course, W and SNP could be interchanged here, we don't actually
    know which is which ;)

    '''

    def _snp_bubble():
        wildtype_sequence = linear_path()
        decision_L = middle_pivot
        decision_R = decision_L + ksize + 1

        if decision_L < 1:
            raise ValueError("ksize too long for length")

        snp_sequence = mutate_position(wildtype_sequence, decision_L + ksize)

        consume_collector(wildtype_sequence, snp_sequence)
        check_fp_collector((lambda G: count_decision_nodes(wildtype_sequence, G, ksize),
                             {(1,2): 1, (2,1): 1}),
                            (lambda G: count_decision_nodes(snp_sequence, G, ksize),
                             {(1,2): 1, (2,1): 1}))

        return (wildtype_sequence, snp_sequence), decision_L, decision_R

    return _snp_bubble


@pytest.fixture
def hourglass_tangle(request, ksize, length, middle_pivot,
                     linear_path, random_sequence,
                     consume_collector, check_fp_collector):
    '''
    The hourglass tangle consists of two sequences, the "top" and "bottom,"
    and has four decision k-mers: two in the top, neighboring each other,
    and two in the bottom, neighboring each other as well, with the left decision
    of the bottom sharing the right decision of the top as a neighbor,
    and the right decision of the bottom having the left decision of the top as
    a neighbor.

        top                
        [0]~~([pivot:pivot+K] L_top)→([pivot+1:pivot+1+K] R_top)~~o
                                   ↘ ↗
                                   ↗ ↘
        [0]~~([pivot:pivot+K] L_bot)→([pivot+1:pivot+1+K] R_bot)~~o
        bottom

    L_top and L_bot have left-degree 1 and right-degree 2,
    while R_top and R_bot have left-degree 2 and right-degree 1.
    '''

    def _hourglass_tangle():
        top = linear_path()
        L = middle_pivot
        decision_segment = top[L:L+ksize+1]
        decision_segment = mutate_position(decision_segment, 0)
        decision_segment = mutate_position(decision_segment, -1)

        bottom = random_sequence(exclude=decision_segment)[:L] \
                 + decision_segment
        bottom += random_sequence(exclude=bottom)[:length-L-len(decision_segment)]

        consume_collector(top, bottom)
        check_fp_collector((lambda G: count_decision_nodes(top, G, ksize),
                            {(1,2): 1, (2,1): 1}),
                           (lambda G: count_decision_nodes(bottom, G, ksize),
                            {(1,2): 1, (2,1): 1}))

        return (top, bottom), L

    return _hourglass_tangle


@pytest.fixture
def triple_chain_tangle(request, ksize, length, middle_pivot,
                     linear_path, random_sequence,
                     consume_collector, check_fp_collector):
    '''
    The triple chain tangle consists of two sequences, the "top" and "bottom,"
    and has six decision k-mers: three in the top, neighboring each other,
    and three in the bottom, neighboring each other as well, with the left decision
    of the bottom sharing the right decision of the top as a neighbor,
    and the right decision of the bottom having the left decision of the top as
    a neighbor.

        top                
        [0]~~([pivot:pivot+K] L_top)→([pivot+1:pivot+1+K] R_top)~~o
                                   ↘ ↗
                                   ↗ ↘
        [0]~~([pivot:pivot+K] L_bot)→([pivot+1:pivot+1+K] R_bot)~~o
        bottom

    L_top and L_bot have left-degree 1 and right-degree 2,
    while R_top and R_bot have left-degree 2 and right-degree 1.
    '''

    def _triple_chain_tangle():
        top = linear_path()
        L = middle_pivot
        decision_segment = top[L:L+ksize+1]
        decision_segment = mutate_position(decision_segment, 0)
        decision_segment = mutate_position(decision_segment, -1)
        third_dsegment = top[L+3:L+3+ksize]
        mutate_position(third_dsegment, -1)
        decision_segment += third_dsegment

        bottom = random_sequence(exclude=decision_segment)[:L] \
                 + decision_segment
        bottom += random_sequence(exclude=bottom)[:length-L-len(decision_segment)]

        consume_collector(top, bottom)
        check_fp_collector((lambda G: count_decision_nodes(top, G, ksize),
                            {(1,2): 2, (2,1): 1}),
                           (lambda G: count_decision_nodes(bottom, G, ksize),
                            {(1,2): 1, (2,1): 2}))

        return (top, bottom), L

    return _triple_chain_tangle


@pytest.fixture
def bowtie_tangle(request, ksize, length, middle_pivot,
                  linear_path, random_sequence,
                  consume_collector, check_fp_collector):

    def _bowtie_tangle():
        top = linear_path()
        decision_segment = top[middle_pivot:middle_pivot+ksize+2]
        decision_segment = mutate_position(decision_segment, 0)
        decision_segment = mutate_position(decision_segment, -1)

        bottom = random_sequence(exclude=decision_segment)[:middle_pivot] \
                 + decision_segment
        bottom += random_sequence(exclude=bottom)[:length-middle_pivot-len(decision_segment)]

        consume_collector(top, bottom)
        check_fp_collector((lambda G: count_decision_nodes(top, G, ksize),
                            {(2,2): 1}),
                           (lambda G: count_decision_nodes(bottom, G, ksize),
                            {(2,2): 1}))

        return (top, bottom), middle_pivot

    return _bowtie_tangle


@pytest.fixture
def left_hairpin(request, ksize, linear_path, consume_collector, check_fp_collector):
    '''
    Sets up a left hairpin graph structure with HDN at pos.
    '''

    def _left_hairpin():
        core = linear_path()
        pos = len(core) // 2
        if core[pos - 1] == core[-1]:
            core = mutate_position(core, -1)
        hdn = core[pos:pos+ksize]
        result = core + hdn

        _collector = consume_collector()
        _collector.pop()
        consume_collector(result)
        _check_fp_collector = check_fp_collector()
        _check_fp_collector.pop()
        check_fp_collector((lambda G: count_decision_nodes(result, G, ksize), {(2,1) : 2}))

        return result, pos

    return _left_hairpin


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
def circular(request, ksize, linear_path, consume_collector, check_fp_collector):

    def _circular():
        sequence = linear_path()
        sequence += sequence

        consume_collector(sequence)
        check_fp_collector((lambda G: count_decision_nodes(sequence, G, ksize), {}))

        return sequence

    return _circular

@pytest.fixture
def suffix_circular(request, ksize, linear_path, consume_collector, check_fp_collector):

    def _suffix_circular():
        sequence = linear_path()
        sequence += sequence[:ksize-1]

        consume_collector(sequence)
        check_fp_collector((lambda G: count_decision_nodes(sequence, G, ksize), {}))

        return sequence

    return _suffix_circular


@pytest.fixture
def suffix_circular_tangle(request, ksize, length, suffix_circular, middle_pivot,
                           random_sequence, consume_collector, check_fp_collector):

    def _suffix_circular_tangle():
        base = suffix_circular()
        
        L = middle_pivot
        decision_segment = base[L:L+ksize+1]
        decision_segment = mutate_position(decision_segment, 0)
        decision_segment = mutate_position(decision_segment, -1)

        inducer = random_sequence(exclude=decision_segment)[:L] \
                 + decision_segment
        inducer += random_sequence(exclude=inducer)[:length-L-len(decision_segment)]

        consume_collector(base, inducer)
        check_fp_collector((lambda G: count_decision_nodes(base, G, ksize),
                            {(1,2): 1, (2,1): 1}),
                           (lambda G: count_decision_nodes(inducer, G, ksize),
                            {(1,2): 1, (2,1): 1}))

        return (base, inducer), L

    return _suffix_circular_tangle

@pytest.fixture
def circular_key(request, ksize, length, pivot,
                 linear_path, consume_collector, check_fp_collector):

    def _circular_key():
        loop = linear_path()
        loop = loop + loop[:ksize-1]
        if pivot in (length - ksize, length-ksize-1):
            _pivot = pivot + ksize - 1
        else:
            _pivot = pivot

        loop_kmers = list(kmers(loop, ksize))
        tail_kmers = [loop_kmers[i % len(loop_kmers)] for i in range(_pivot+1, _pivot+1+ksize)]
        tail = ''.join((kmer[0] for kmer in tail_kmers))
        tail = mutate_position(tail, -1)
        
        consume_collector(loop, tail)
        check_fp_collector((lambda G: count_decision_nodes(loop, G, ksize), {(1,2) : 1}))

        return (loop, tail), _pivot

    return _circular_key
