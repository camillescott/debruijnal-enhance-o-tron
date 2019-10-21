#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2019
# File   : bubble.py
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 16.05.2019

from typing import List

from debruijnal_enhance_o_tron.generators.generator import SequenceGenerator

class BubbleGenerator(SequenceGenerator):
    
    def generate(self, n_branches: int = 2,
                       n_tail_kmers: int = 0,
                       visualize: bool = False) -> List[str]:
        if n_branches > len(self.Sigma) or n_branches < 2:
            raise ValueError('Invalid number of branches.')
        
        # Generate a seed k-mer and make n_branches copies
        start = self.random_unitig(self.ksize + n_tail_kmers)
        sequences = [start] * n_branches
        
        # extend each sequence with a mutant
        for i in range(n_branches):
            sequences[i] = self.random_unitig_extend(sequences[i], N=1)
        # extend a new k-mer off the end and append it to each sequence
        end = self.random_unitig_extend(sequences[0], N=self.ksize+n_tail_kmers)[len(sequences[0]):]
        for i in range(len(sequences)):
            self.add(self.suffix(sequences[i]) + end)
            sequences[i] += end
        
        if visualize:
            self.visualize(sequences, n_branches, n_tail_kmers)
        
        return sequences

    def visualize(self, sequences: List[str],
                        n_branches: int,
                        n_tail_kmers: int) -> None:
        out = '\n' + (' ' * n_tail_kmers) + '|' + (' ' * (self.ksize - 2)) + '|*' +\
              '|' + (' ' * (self.ksize - 2)) + '|\n'
        print(out.join(sequences))
