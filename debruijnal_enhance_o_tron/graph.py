#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2019
# File   : graph.py
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 16.05.2019

from typing import Any, Iterable, Optional, Generator, List, Set, Tuple, TypeVar

T = TypeVar('T', bound='GraphAdapter')

class GraphAdapter:
    '''A basic de Bruijn Graph. End-users can either use as is, or subclass it and provide
    the underscore _query, _insert, and _insert_many methods along with their own store to
    use their own dBG implementations.
    '''

    mask: set
    tags: Set[str]
    ksize: int
    Sigma: str
    Sigma_set: Set[str]
    Sigma_list: List[str]

    def __init__(self, ksize: int,
                       Sigma: str = 'ACGT',
                       store: Optional[Any] = None,
                       *args,
                       **kwargs) -> None:

        self.store = store if store is not None else set()
        self.mask = set()
        self.tags = set()
        self.ksize = ksize
        self.Sigma = Sigma
        self.Sigma_set = set(self.Sigma)
        self.Sigma_list = list(self.Sigma)

    def _query(self, item: str) -> bool:
        '''Override in subclasses'''
        return item in self.store

    def _insert(self, item: str) -> None:
        '''Override in subclasses'''
        self.store.add(item)
    
    def _insert_many(self, items: Iterable[str]) -> None:
        '''Override in subclasses'''
        self.store.update(items)

    def reset(self) -> None:
        '''Reset the dBG to an empty state.
        '''
        self.store.clear()
        self.mask.clear()

    def shallow_clone(self) -> T:
        return GraphAdapter(self.ksize, self.Sigma)

    def get(self, item: str) -> bool:
        '''Return true if the item is in the dBG, false otherwise
        '''
        def query(_item):
            return (_item not in self.mask) and self._query(_item)

        if len(item) < self.ksize:
            raise ValueError(item)
        elif len(item) == self.ksize:
            return query(item)
        else:
            return [query(_item) for _item in self.kmers(item)]


    def add(self, item: str) -> None:
        '''Add the k-mer(s) from the given sequence to the dBG.
        '''
        if len(item) < self.ksize:
            raise ValueError(item)
        elif len(item) == self.ksize:
            self._insert(item)
        else:
            self._insert_many(self.kmers(item))
    
    def tag(self, item: str,
                  dist: int = 3) -> int:
        if len(item) < self.ksize:
            raise ValueError(item)
        count = 0
        n_tagged = 0
        for kmer in self.kmers(item):
            if kmer in self.tags:
                count = 0
            else:
                count += 1
                if count > dist:
                    self.tags.add(kmer)
                    count = 0
                    n_tagged += 1
        return n_tagged
    
    def remove(self, item: str) -> None:
        '''Remove the k-mer(s) from the given sequence from the dBG.

        Seeing as many storage backends don't support removal, we keep
        a mask of removed k-mers and query it as well.
        '''
        if len(item) < self.ksize:
            raise ValueError(item)
        elif len(item) == self.ksize:
            self.mask.add(item)
        else:
            self.mask.add(self.kmers(item))

    def kmers(self, sequence: str) -> Generator[str, None, None]:
        '''Yields the k-mers in the given sequence.
        '''
        for i in range(len(sequence) - self.ksize + 1):
            yield sequence[i:i+self.ksize]
    
    def suffix(self, sequence: str) -> str:
        '''Returns the length k-1 suffix of the sequence.
        '''
        if len(sequence) < self.ksize:
            raise ValueError('Sequence too short')
        return sequence[-self.ksize+1:]
    
    def prefix(self, sequence: str) -> str:
        '''Returns the length k-1 prefix of the sequence.
        '''
        if len(sequence) < self.ksize:
            raise ValueError('Sequence too short')
        return sequence[:self.ksize-1]
    
    def left_neighbors(self, item: str) -> Generator[str, None, None]:
        '''Yields the four left neighbors of the given k-mer.
        '''
        prefix = self.prefix(item)
        for b in self.Sigma:
            yield b + prefix
    
    def right_neighbors(self, item: str) -> Generator[str, None, None]:
        '''Yields the four right neighbors of the given k-mer.
        '''
        suffix = self.suffix(item)
        for b in self.Sigma:
            yield suffix + b

    def left_degree(self, item: str) -> int:
        '''Yields the in-degree (left degree) of the given k-mer.
        '''
        return sum((self.get(neighbor) for neighbor in self.left_neighbors(item)))

    def right_degree(self, item: str) -> int:
        '''Yields the out-degree (right degree) of the given k-mer.
        '''
        return sum((self.get(neighbor) for neighbor in self.right_neighbors(item)))
    
    def degree(self, item: str) -> int:
        return self.left_degree(item), self.right_degree(item)
    
    def is_decision(self, kmer: str) -> bool:
        '''Returns True if the given k-mer is a decision k-mer in the dBG:
        that is, if its in-degree or out-degree is greater than one.
        '''
        return self.left_degree(kmer) > 1 or self.right_degree(kmer) > 1

    def find_decisions(self, sequence: str) -> Generator[Tuple[int, str], None, None]:
        '''Yields the position and k-mer sequence of the decision k-mers in the 
        given sequence.
        '''
        for n, kmer in enumerate(self.kmers(sequence)):
            if self.is_decision(kmer):
                yield n, kmer

    def count_decision_nodes(self, sequence):
        dnodes = {}
        for i, kmer in enumerate(self.kmers(sequence)):
            d = (self.left_degree(kmer), self.right_degree(kmer))
            ld, rd = d
            if ld > 1 or rd > 1:
                dnodes[d] = dnodes.get(d, 0) + 1
                print(i, d)

        return dnodes
