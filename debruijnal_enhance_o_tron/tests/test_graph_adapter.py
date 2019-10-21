#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2019
# File   : test_graph_adapter.py
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 20.05.2019
from itertools import chain

import pytest

from debruijnal_enhance_o_tron.fixtures import subgraphs
from debruijnal_enhance_o_tron.fixtures.sequence import (using_ksize,
                                                         using_length)
from debruijnal_enhance_o_tron.fixtures.collectors import (consume, check_fp)
from debruijnal_enhance_o_tron.graph import GraphAdapter
from debruijnal_enhance_o_tron.sequence import kmers, get_random_sequence



@pytest.fixture
def graph(ksize):
    '''Test override of conftest-injected graph fixture.
    '''

    return GraphAdapter(ksize)


def test_linear_path_noconsume(linear_path, graph, ksize, length):
    sequence = linear_path()

    for kmer in kmers(sequence, ksize):
        assert not graph.get(kmer)


def test_linear_path_consume(linear_path, graph, ksize, length, consume):
    sequence = linear_path()
    consume()

    assert subgraphs.count_decision_nodes(sequence, graph, ksize) == {}

    for kmer in kmers(sequence, ksize):
        assert graph.get(kmer)


def test_right_comb_consume(right_comb, graph, ksize, tip_length, n_branches, consume):
    seqs = right_comb()
    consume()

    assert subgraphs.count_decision_nodes(seqs[0], graph, ksize) == {(0, n_branches): 1}

    for kmer in chain(*(kmers(s, ksize) for s in seqs)):
        print(kmer)
        assert graph.get(kmer)


def test_left_comb_consume(left_comb, graph, ksize, tip_length, n_branches, consume):
    seqs = left_comb()
    consume()

    assert subgraphs.count_decision_nodes(seqs[0], graph, ksize) == {(n_branches, 0): 1}

    for kmer in chain(*(kmers(s, ksize) for s in seqs)):
        assert graph.get(kmer)


def test_right_tip_noconsume(right_tip, graph, ksize, length):
    (sequence, tip), S = right_tip()

    for kmer in chain(kmers(sequence, ksize),
                      kmers(tip, ksize)):
        assert not graph.get(kmer)


def test_right_tip_consume(right_tip, graph, ksize, length, consume):
    (sequence, tip), S = right_tip()
    consume()

    assert subgraphs.count_decision_nodes(sequence, graph, ksize) == {(1,2): 1}

    for kmer in chain(kmers(sequence, ksize),
                      kmers(tip, ksize)):
        assert graph.get(kmer)


def test_right_fork_noconsume(right_fork, graph, ksize, length):
    (sequence, branch), S = right_fork()
    assert len(sequence) == length

    for kmer in chain(kmers(sequence, ksize),
                      kmers(branch, ksize)):
        assert not graph.get(kmer)


def test_right_fork_consume(right_fork, graph, ksize, length, consume):
    (sequence, branch), S = right_fork()
    consume()

    assert subgraphs.count_decision_nodes(sequence, graph, ksize) == {(1,2): 1}

    for kmer in chain(kmers(sequence, ksize),
                      kmers(branch, ksize)):
        assert graph.get(kmer)


def test_left_fork_consume(left_fork, graph, ksize, length, consume):
    (sequence, branch), pos = left_fork()
    consume()

    assert subgraphs.count_decision_nodes(sequence, graph, ksize) == {(2,1): 1}
    assert graph.left_degree(sequence[pos:pos+ksize]) == 2
    assert graph.right_degree(sequence[pos:pos+ksize]) == 1

    for kmer in chain(kmers(sequence, ksize),
                      kmers(branch, ksize)):
        assert graph.get(kmer)


def test_right_triple_fork_noconsume(right_triple_fork, graph, ksize, length):
    (core, top, bottom), S = right_triple_fork()

    for kmer in chain(kmers(core, ksize),
                      kmers(top, ksize),
                      kmers(bottom, ksize)):
        assert not graph.get(kmer)


def test_right_triple_fork_consume(right_triple_fork, graph, 
                                   ksize, length, consume):
    (core, top, bottom), S = right_triple_fork()
    consume()

    assert subgraphs.count_decision_nodes(core, graph, ksize) == {(1,3): 1}

    for kmer in chain(kmers(core, ksize),
                      kmers(top, ksize),
                      kmers(bottom, ksize)):
        assert graph.get(kmer)


def test_snp_bubble_noconsume(snp_bubble, graph, ksize, length):
    (wildtype, snp), L, R = snp_bubble()

    for kmer in chain(kmers(wildtype, ksize),
                      kmers(snp, ksize)):
        assert not graph.get(kmer)


def test_snp_bubble_consume(snp_bubble, graph, ksize, length, consume, check_fp):
    (wildtype, snp), L, R = snp_bubble()
    consume()
    check_fp()

    assert subgraphs.count_decision_nodes(wildtype, graph, ksize) \
            == {(1,2): 1, (2,1):1}
    assert subgraphs.count_decision_nodes(snp, graph, ksize) \
            == {(1,2): 1, (2,1):1}

    for kmer in chain(kmers(wildtype, ksize),
                      kmers(snp, ksize)):
        assert graph.get(kmer)


def test_left_hairpin_consume(left_hairpin, graph, ksize, length, consume):
    sequence, pos = left_hairpin()
    consume()

    for kmer in kmers(sequence, ksize):
        assert graph.get(kmer)

    print(sequence, pos)
    assert subgraphs.count_decision_nodes(sequence, graph, ksize) == {(2,1): 2}


def test_tandem_quad_forks(tandem_quad_forks, graph, ksize, length, consume):
    (core, left_branches, right_branches), S_l, S_r = tandem_quad_forks()
    consume()

    assert subgraphs.count_decision_nodes(core, graph, ksize) == {(1,4): 2}
    
    assert graph.left_degree(core[S_l:S_l+ksize]) == 1
    assert graph.right_degree(core[S_l:S_l+ksize]) == 4

    assert graph.left_degree(core[S_r:S_r+ksize]) == 1
    assert graph.right_degree(core[S_r:S_r+ksize]) == 4


def test_hourglass_tangle(hourglass_tangle, graph, ksize, length, consume):
    (top, bottom), L = hourglass_tangle()
    assert len(bottom) == len(top)
    consume()

    for kmer in chain(kmers(top, ksize), kmers(bottom, ksize)):
        assert graph.get(kmer)

    assert subgraphs.count_decision_nodes(top, graph, ksize) == {(1,2): 1, (2,1): 1}
    assert subgraphs.count_decision_nodes(bottom, graph, ksize) == {(1,2): 1, (2,1): 1}
    assert graph.left_degree(top[L:L+ksize]) == 1
    assert graph.right_degree(top[L:L+ksize]) == 2
    assert graph.left_degree(bottom[L:L+ksize]) == 1
    assert graph.right_degree(bottom[L:L+ksize]) == 2


'''
@using_ksize(5)
@using_length(25)
def test_triple_chain_tangle(triple_chain_tangle, graph, ksize, length, consume):
    (top, bottom), L = triple_chain_tangle()
    print(L, top, L * ' ' + '|', bottom, sep='\n')
    assert len(bottom) == len(top)
    consume()

    for kmer in chain(kmers(top, ksize), kmers(bottom, ksize)):
        assert graph.get(kmer)

    assert subgraphs.count_decision_nodes(top, graph, ksize) == {(1,2): 2, (2,1): 1}
    assert subgraphs.count_decision_nodes(bottom, graph, ksize) == {(1,2): 1, (2,1): 2}
    assert graph.left_degree(top[L:L+ksize]) == 1
    assert graph.right_degree(top[L:L+ksize]) == 2
    assert graph.left_degree(bottom[L:L+ksize]) == 1
    assert graph.right_degree(bottom[L:L+ksize]) == 2
'''


def test_bowtie_tangle(bowtie_tangle, graph, ksize, length, consume):
    (top, bottom), L = bowtie_tangle()
    assert len(bottom) == len(top)
    consume()

    for kmer in chain(kmers(top, ksize), kmers(bottom, ksize)):
        assert graph.get(kmer)

    assert subgraphs.count_decision_nodes(top, graph, ksize) == {(2,2): 1}
    assert subgraphs.count_decision_nodes(bottom, graph, ksize) == {(2,2): 1}
    assert graph.left_degree(top[L+1:L+1+ksize]) == 2
    assert graph.right_degree(top[L+1:L+1+ksize]) == 2


def test_circular_key(circular_key, graph, ksize, length, consume, check_fp):
    (loop, tail), pos = circular_key()
    consume()

    assert graph.left_degree(loop[pos:pos+ksize]) == 1
    assert graph.right_degree(loop[pos:pos+ksize]) == 2


def test_tandem_repeat_lt_ksize_noconsume(tandem_repeats_lt_ksize,
                                          ksize,
                                          graph):
    (repeat, tandem_repeats), n_repeats = tandem_repeats_lt_ksize()

    for kmer in kmers(tandem_repeats, ksize):
        assert not graph.get(kmer)


def test_tandem_repeat_lt_ksize_consume(tandem_repeats_lt_ksize,
                                          ksize,
                                          graph,
                                          consume):
    (repeat, tandem_repeats), n_repeats = tandem_repeats_lt_ksize()
    consume()

    assert not subgraphs.count_decision_nodes(tandem_repeats, graph, ksize)

    for kmer in kmers(tandem_repeats, ksize):
        assert graph.get(kmer)


def test_tandem_repeat_gt_ksize_noconsume(tandem_repeats_gt_ksize,
                                          ksize,
                                          graph):
    (repeat, tandem_repeats), n_repeats = tandem_repeats_gt_ksize()

    for kmer in kmers(tandem_repeats, ksize):
        assert not graph.get(kmer)


def test_tandem_repeat_gt_ksize_consume(tandem_repeats_gt_ksize,
                                          ksize,
                                          graph,
                                          consume):
    (repeat, tandem_repeats), n_repeats = tandem_repeats_gt_ksize()
    consume()

    assert not subgraphs.count_decision_nodes(tandem_repeats, graph, ksize)

    for kmer in kmers(tandem_repeats, ksize):
        assert graph.get(kmer)


def test_circular_noconsume(circular, graph, length, ksize):
    sequence = circular()
    
    for kmer in kmers(sequence, ksize):
        assert not graph.get(kmer)


def test_circular_consume(circular, graph, length, ksize, consume):
    sequence = circular()
    consume()

    assert not subgraphs.count_decision_nodes(sequence, graph, ksize)
    assert graph.left_degree(sequence[:ksize]) == 1

    for kmer in kmers(sequence, ksize):
        assert graph.get(kmer)


@using_ksize([11,13,15])
@using_length(500)
def test_sequence_generator_degree(graph, length, ksize, random_sequence):
    seqs = [random_sequence() for _ in range(10)]
    for seq in seqs:
        graph.add(seq)
    for seq in seqs:
        seq_kmers = list(kmers(seq, ksize))
        assert graph.left_degree(seq_kmers[0]) == 0
        assert graph.right_degree(seq_kmers[0]) == 1
        for kmer in seq_kmers[1:-1]:
            assert graph.left_degree(kmer) == 1
            assert graph.right_degree(kmer) == 1
        assert graph.left_degree(seq_kmers[-1]) == 1
        assert graph.right_degree(seq_kmers[-1]) == 0


@using_ksize([7,9,11])
@using_length(500)
def test_sequence_generator_degree_with_fp(graph, length, ksize, linear_path, check_fp):
    seqs = [linear_path() for _ in range(10)]
    check_fp()

    for seq in seqs:
        graph.add(seq)
    for seq in seqs:
        seq_kmers = list(kmers(seq, ksize))
        assert graph.left_degree(seq_kmers[0]) == 0
        assert graph.right_degree(seq_kmers[0]) == 1
        for kmer in seq_kmers[1:-1]:
            assert graph.left_degree(kmer) == 1
            assert graph.right_degree(kmer) == 1
        assert graph.left_degree(seq_kmers[-1]) == 1
        assert graph.right_degree(seq_kmers[-1]) == 0


