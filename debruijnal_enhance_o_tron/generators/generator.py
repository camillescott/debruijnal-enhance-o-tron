#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2019
# File   : generator.py
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 16.05.2019

import random
from typing import Any, Optional, Generator, List, Set

from debruijnal_enhance_o_tron.graph import GraphAdapter

class SequenceGenerator:
    ''' Base driver for random sequence generators.

    '''
    max_tries: int
    graph:     GraphAdapter

    def __init__(self, ksize: int = 13,
                       graph: Optional[GraphAdapter] = None,
                       max_tries: int = 100,
                       rseed: Optional[int] = None,
                       *args,
                       **kwargs) -> None:
        self.max_tries = max_tries
        if graph is None:
            self.graph = GraphAdapter(ksize)
        else:
            if not isinstance(graph, GraphAdapter):
                raise TypeError('graph must be a {0}'.format(GraphAdapter.__name__))
            self.graph = graph

        random.seed(rseed)
    
    def __getattr__(self, attr: Any) -> Any:
        return getattr(self.graph, attr)
    
    def generate(self, length: int, *args, **kwargs) -> str:
        return self.random_lmer(length)
    
    def validate(self, *args, **kwargs):
        raise NotImplementedError()
    
    def random_lexicographic_ordering(self) -> Generator[str, None, None]:
        
        def _incr(base):
            return self.Sigma[self.Sigma.find(base) + 1]
        
        def _next(kmer):
            
            i = self.ksize - 1
            while kmer[i] == self.Sigma[-1] and i >= 0:
                i -= 1
            if i == -1:
                return self.Sigma[0] * self.ksize
            else:
                kmer = list(kmer)
                kmer[i] = _incr(kmer[i])
                
                return ''.join(kmer[:i+1]) + self.Sigma[0] * (self.ksize - i - 1)
        
        def ordering():
            start = self.random_kmer()
            yield start
            kmer = _next(start)
            while kmer != start:
                yield kmer
                kmer = _next(kmer)
        
        return ordering

    def random_base(self) -> str:
        '''Return a random base from the alphabet Sigma
        '''
        return random.choice(self.Sigma)
    
    def mutate_base(self, base: str,
                          exclude: Optional[Set[str]] = None) -> str:
        '''Mutate base, excluding the given options as possible targets.
        '''
        if exclude is None:
            exclude = set()
    
        exclude.add(base)
        if len(exclude) == len(self.graph.Sigma):
            raise ValueError('Tried to exclude all bases.')
        return random.choice(tuple(self.Sigma_set - exclude))
    
    def mutate_position(self, sequence: str,
                              position: int,
                              exclude: Optional[Set[str]] = None) -> str:
        '''Return a copy of sequence mutated at the given position.
        '''
        base = self.mutate_base(sequence[position], exclude=exclude)
        sequence = list(sequence)
        sequence[position] = base
        return ''.join(sequence)
        
    def randomize_sigma(self) -> str:
        '''Returns the underlying alphabet Sigma in randomized order.
        '''
        return random.sample(self.Sigma, len(self.Sigma))
    
    def random_lmer(self, l: int) -> str:
        '''Generate a random sequence of length l
        '''
        return ''.join((self.random_base() for _ in range(l)))
    
    def random_kmer(self) -> str:
        '''Generates a random k-mer of length ksize.
        '''
        return self.random_lmer(self.ksize)
    
    def random_seed(self) -> str:
        '''Generates a random k-mer with no neighbors.
        '''
        seed = self.random_kmer()
        while self.graph.left_degree(seed) or self.graph.right_degree(seed):
            seed = self.random_kmer()
        self.tags.add(seed)
        return seed
    
    def random_unitig_base(self, seed: str,
                                 fw: bool = True) -> str:
        '''Greedily unitig extend a single base in the given direction.
        '''
        if fw:
            neighbors = self.right_neighbors
        else:
            neighbors = self.left_neighbors
        
        extensions = self.randomize_sigma()
        result = None
        for base in extensions:
            if fw:
                kmer = self.suffix(seed) + base
            else:
                kmer = base + self.prefix(seed)
                
            if self.get(kmer):
                continue

            pre_add_decision = list((self.is_decision(neighbor) for neighbor in neighbors(kmer)))
            self.add(kmer)
            post_add_decision = list((self.is_decision(neighbor) for neighbor in neighbors(kmer)))
            if self.is_decision(kmer) or pre_add_decision != post_add_decision:
                self.remove(kmer)
            else:
                result = base
                break
        
        return result
    
    def random_new_base(self, seed: str,
                              fw: bool = True) -> str:
        extensions = self.randomize_sigma()
        result = None
        for base in extensions:
            if fw:
                kmer = self.suffix(seed) + base
            else:
                kmer = base + self.prefix(seed)
            if self.get(kmer):
                continue
            else:
                self.add(kmer)
                result = base
                break

        return result
    
    def random_degree_base(self, seed: str,
                                 degree: int,
                                 fw: bool = True) -> str:
        extensions = self.randomize_sigma()
        result = None
        for base in extensions:
            if fw:
                kmer = self.suffix(seed) + base
            else:
                kmer = base + self.prefix(seed)
            if self.get(kmer):
                continue
            self.add(kmer)
            if self.degree(kmer) != degree:
                self.remove(kmer)
            else:
                result = base
                break
    
        return result
    
    def random_unitig_extend(self, seed: str,
                                   fw: bool = True,
                                   N: int = 1) -> str:
        '''Greedily extend the given seed by N bases such that any added sequence is a unitig. That is,
        added k-mers will not be decision k-mers, nor will they decision-induce any k-mers in the dBG.
        
        May get stuck in a dead-end and raise ValueError.
        '''
        
        # bind some stuff to function scope
        ksize = self.ksize
        random_lmer = self.random_lmer
        graph = self.graph
        get = self.get
        add = self.add
        remove = self.remove
        is_decision = self.is_decision
        
        # First step is to extend the seed out to length K if it isn't already
        _seed = seed

        if len(seed) < ksize:
            seed = _seed + random_lmer(ksize - len(_seed))
            add(seed)
            n_tries = 0
            while is_decision(seed):
                remove(seed)
                seed = _seed + random_lmer(ksize - len(_seed))
                add(seed)
                n_tries += 1
                if n_tries > self.max_tries:
                    raise ValueError('Exceeded max_tries when generating'
                                     'non-decision seed ({0} tries)'.format(n_tries))
       
        self.add(seed)
        while N > 0:
            
            extension = self.random_unitig_base(seed, fw=fw)
            if extension is None:
                raise ValueError('Greedy dead-end, try again')
            else:
                if fw:
                    seed += extension
                else:
                    seed = extension + seed
            
            N = N - 1
            
        return seed
    
    def random_unitig(self, L: int) -> str:
        '''Generate a random unitig of length L >= K.
        '''
        if L < self.ksize:
            raise ValueError('length is < K')
        seed = self.random_seed()
        return self.random_unitig_extend(seed, N=L-self.ksize)
    
    def random_branches(self, seed: str,
                              n_branches: int,
                              n_tail_kmers: int = 1,
                              fw: bool = True) -> Generator[str, None, None]:
        if n_branches > len(self.Sigma) or n_branches < 1:
            raise ValueError('Invalid number of branches.')
        
        yield from (self.random_unitig_extend(seed, fw=fw, N=n_tail_kmers) for _ in range(n_branches))
    
    def get_decision_neighborhood(self, sequence: str) -> List[List[bool]]:
        kmers = list(self.kmers(sequence))
        neighborhood = []
        neighborhood.append(list((self.is_decision(neighbor) for neighbor in self.left_neighbors(kmers[0]))))
        if len(kmers) > 2:
            for kmer in kmers[1:-1]:
                neighborhood.append(list((self.is_decision(neighbor) for neighbor in self.left_neighbors(kmer))))
                neighborhood.append(list((self.is_decision(neighbor) for neighbor in self.right_neighbors(kmer))))
        neighborhood.append(list((self.is_decision(neighbor) for neighbor in self.right_neighbors(kmers[-1]))))
        
        return neighborhood
                
    def has_decision_inducers(self, segment: str) -> bool:
        before = self.get_decision_neighborhood(segment)
        self.add(segment)
        after = self.get_decision_neighborhood(segment)
        self.remove(segment)
        
        return before != after

    def random_unitig_merge(self, left: str,
                                  right: str,
                                  n_kmers: int = 3) -> str:
        '''Given unitigs left and right, merge them together by adding
        n_kmers bases between them such that the merged sequence is still a
        unitig. If this is not possible, raises ValueError.'''
        
        _left = self.suffix(left)
        _right = self.prefix(right)
        
        seeder = self.random_lexicographic_ordering()
        n_tries = 0
        for seed in seeder():
            n_tries += 1
            segment = _left + seed + _right
            print(seed)
            if any(self.get(segment)):
                print('segment had existing k-mers')
                continue
            if not self.has_decision_inducers(segment):
                self.add(segment)
                if list(self.find_decisions(segment)) == []:
                    return left + seed + right
                else:
                    self.remove(segment)
            
            
        raise ValueError('Unable to merge without inducing new decision k-mers, tried', n_tries)
