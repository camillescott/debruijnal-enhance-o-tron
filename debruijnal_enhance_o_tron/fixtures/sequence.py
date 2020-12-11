#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2019
# File   : sequence.py
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 20.05.2019
import pytest

from debruijnal_enhance_o_tron.sequence import (mutate_base,
                                                mutate_sequence,
                                                mutate_position,
                                                get_random_sequence,
                                                reads,
                                                kmers,
                                                revcomp)
from debruijnal_enhance_o_tron.generators.generator import SequenceGenerator

@pytest.fixture(params=[21], ids=lambda k: 'K={0}'.format(k))
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


@pytest.fixture(params=[1, None], ids=lambda l: 'L=K*2' if l is None else 'L=1')
def tip_length(request, ksize):
    if request.param is None:
        return ksize * 2
    else:
        return request.param


@pytest.fixture(params=[2,3,4], ids=lambda l: 'branches={0}'.format(l))
def n_branches(request):
    return request.param


def using_n_branches(n):
    def wrapped(fixture_func):
        if isinstance(n, int):
            N = [n]
        else:
            N = list(n)
        return pytest.mark.parametrize('n_branches', 
                                       N,
                                       indirect=['n_branches'],
                                       ids=lambda l: 'branches={0}'.format(l))(fixture_func)
    return wrapped


def get_internal_pivot(pivot_var, length, ksize):
    if pivot_var == 'flank_left':
        return 1
    elif pivot_var == 'middle':
        return length // 2
    elif pivot_var == 'flank_right':
        return length - ksize - 1
    else:
        raise ValueError("Invalid pivot")


def get_end_pivot(pivot_var, length, ksize):
    if pivot_var == 'left':
        return 0
    elif pivot_var == 'right':
        return length - ksize
    else:
        raise ValueError("Invalid pivot")


@pytest.fixture(params=['left', 'right'], ids=lambda l: 'pivot={0}'.format(l))
def end_pivot(request, length, ksize):
    try:
        return get_end_pivot(request.param, length, ksize)
    except ValueError:
        return request.param


@pytest.fixture(params=['flank_left', 'middle', 'flank_right'],
                ids=lambda l: 'pivot={0}'.format(l))
def internal_pivot(request, length, ksize):
    try:
        return get_internal_pivot(request.param, length, ksize)
    except ValueError:
        return request.param


@pytest.fixture(params=['middle'], ids=lambda l: 'pivot={0}'.format(l))
def middle_pivot(request, length, ksize):
    try:
        return get_internal_pivot(request.param, length, ksize)
    except ValueError:
        return request.param


@pytest.fixture(params=['left', 'flank_left', 'middle', 'flank_right', 'right'],
                ids=lambda l: 'pivot={0}'.format(l))
def pivot(request, length, ksize):
    try:
        return get_internal_pivot(request.param, length, ksize)
    except ValueError:
        pass

    try:
        return get_end_pivot(request.param, length, ksize)
    except ValueError:
        return request.param


def using_pivot(P):
    '''
    Convenience wrapper for basic length fixure.
    '''

    def wrapped(fixture_func):
        try:
            _ = iter(P)
        except:
            piv = [P]
        else:
            piv = P

        return pytest.mark.parametrize('pivot', 
                                       piv,
                                       indirect=['pivot'],
                                       ids=lambda l: 'pivot={0}'.format(l))(fixture_func)
    return wrapped


@pytest.fixture
def sequence_generator(request, ksize):
    return SequenceGenerator(ksize)
    
    

@pytest.fixture
def random_sequence(request, ksize, length):
    global_seen = set()

    def get(exclude=None, length=length):

        try:
            sequence = get_random_sequence(length, 
                                           ksize,
                                           exclude=exclude,
                                           seen=global_seen)
        except ValueError:
            pytest.xfail('Failed to generate sequence of given length and k-size')

        for kmer in kmers(sequence, ksize-1):
            global_seen.add(kmer)

        return sequence

    return get

