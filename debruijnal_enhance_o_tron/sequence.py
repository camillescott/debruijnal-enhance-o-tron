#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2019
# File   : sequence.py
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 18.07.2019
import random

complements = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 


def comp(base):
    base.upper()
    return complements.get(base, base)


def revcomp(sequence):
    '''Perform reverse complement.
    '''
    bases = reversed([complements.get(base, base) for base in sequence])
    return ''.join(bases)


def mutate_base(base):
    '''Non complementary mutation of the given base.
    '''
    if base in 'AT':
        return random.choice('GC')
    elif base in 'GC':
        return random.choice('AT')
    else:
        assert False, 'bad base'


def mutate_sequence(sequence, N=1):
    '''Return a copy of the sequence with N mutations at random
    positions.
    '''
    sequence = list(sequence)
    positions = random.sample(range(len(sequence)), N)

    for i in positions:
        sequence[i] = mutate_base(sequence[i])

    return ''.join(sequence)


def mutate_position(sequence, pos):
    '''Return a copy of the sequence with a point mutation
    at the given position.
    '''
    sequence = list(sequence)
    sequence[pos] = mutate_base(sequence[pos])
    return ''.join(sequence)


def get_random_base(bases='ACGT'):
    '''Return a random nucleotide from the standard alphabet.
    '''
    return random.choice(bases)


def _get_random_sequence(length):
    '''Generate a random sequence without care for overlaps.
    '''

    return ''.join((get_random_base() for _ in range(length)))

def get_random_sequence(length, ksize, exclude=None, seen=None):
    '''Generate a random (non-looping) nucleotide sequence.

    To be non-overlapping, the sequence should not include any repeated
    length K-1 k-mers.

    Args:
        exclude (str): If not None, add the k-mers from this sequence to the
        seen set.
        seen (set): An existing set of seen k-mers.

    Returns:
        str: A random non-looping sequence.
    '''

    seen = set() if seen is None else seen.copy()

    if exclude is not None:
        for kmer in kmers(exclude, ksize - 1):
            seen.add(kmer)
    
    # start off the sequence, make sure it doesn't overlap
    # with our exlcusions
    seq = _get_random_sequence(ksize - 1)
    while seq in seen:
        seq = _get_random_sequence(ksize - 1)
    seen.add(seq)
    seq = list(seq)

    # now extend the sequence one base at a time, checking k-1-mer
    # suffix against seen set
    i = 0
    while(len(seq) < length):
        if i > length * 4:
            raise ValueError('K too small for request length.')
        i += 1
        next_base = get_random_base()
        next_kmer = ''.join(seq[-ksize + 2:] + [next_base])
        assert len(next_kmer) == ksize - 1

        if (next_kmer) not in seen:
            seq.append(next_base)
            seen.add(next_kmer)
        else:
            continue

    return ''.join(seq)


def reads(sequence, ksize, L=100, N=100, dbg_cover=False):
    '''Yields simulated reads over the given sequence.
    '''
    positions = list(range(len(sequence) - L))
    if dbg_cover is True:
        for start in range(0, len(sequence), ksize):
            read = sequence[start:start + L]
            if len(read) < ksize:
                read = sequence[-L:]
            yield read
            N -= 1
    if N < 0:
        return
    for i in range(N):
        start = random.choice(positions)
        yield sequence[start:start + L]


def kmers(sequence, K):
    '''Iterate k-mers for the given sequence.
    '''
    for i in range(len(sequence) - K + 1):
        yield sequence[i:i + K]

def left_kmers(kmer):
    for b in 'ACGT':
        yield b + kmer[:-1]


def right_kmers(kmer):
    for b in 'ACGT':
        yield kmer[1:] + b
