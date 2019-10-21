#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2019
# File   : fork.py
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 16.05.2019

from typing import List, Optional

from debruijnal_enhance_o_tron.generators.generator import SequenceGenerator


class ForkGenerator(SequenceGenerator):
    
    def generate(self, n_branches: int = 2,
                       n_tail_kmers: Optional[int] = None,
                       visualize: bool = False) -> List[str]:
        if n_tail_kmers is None:
            n_tail_kmers = 1
        if n_tail_kmers < 1:
            raise ValueError('n_tail_kmers must be at least 1.')
        
        start = self.random_unitig(self.ksize + n_tail_kmers)
        sequences = [start] * n_branches

        # extend each sequence with a mutant
        for i in range(n_branches):
            sequences[i] = self.random_unitig_extend(sequences[i], N=n_tail_kmers)
        
        if visualize:
            self.visualize(sequences, n_branches, n_tail_kmers)
        
        return sequences
    
    def visualize(self, sequences: List[str],
                        n_branches: int,
                        n_tail_kmers: int) -> None:
        out = '\n' + (' ' * n_tail_kmers) + '|' + (' ' * (self.ksize - 2)) + '|*\n'
        print(out.join(sequences))
