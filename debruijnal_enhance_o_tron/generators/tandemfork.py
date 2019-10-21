#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2019
# File   : tandemfork.py
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 16.05.2019

from typing import List, Optional

from debruijnal_enhance_o_tron.generators.generator import SequenceGenerator


class TandemForkGenerator(SequenceGenerator):
    
    def generate(self, n_forks: int = 2,
                       n_branches: int = 2,
                       n_tail_kmers: Optional[int] = None,
                       visualize: bool = False) -> List[str]:
        if n_tail_kmers is None:
            n_tail_kmers = 1
        if n_tail_kmers < 1:
            raise ValueError('n_tail_kmers must be at least 1.')
        
        if n_forks < 2:
            raise ValueError('Invalid n_forks, must be >= 2 (otherwise they would not be tandem)')
        
        forks = []
        core = self.random_unitig(self.ksize + n_forks + (n_tail_kmers * 2))
        for start_kmer in range(n_tail_kmers, n_tail_kmers + n_forks):
            seed = core[:start_kmer+self.ksize]
            for branch in self.random_branches(seed, n_branches, n_tail_kmers=n_tail_kmers):
                forks.append(branch)
        core = self.random_unitig_extend(core, n_tail_kmers)
        
        return [core] + forks
