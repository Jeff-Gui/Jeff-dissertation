#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  3 19:00:36 2022

@author: jefft
"""

from utils import rand_part, rand_merge_group_layers
import numpy as np
from Bio import AlignIO
#======= Load sample-feature matrix and pre-processing (dimensional reduction?)


#======= Annotate each sample with random partition filters
hash_lookup = []
MIN_SAMPLE_PER_LEAF = 20
time = 0
success = 0
mtx = np.stack([np.random.randint(0,5,(1000,)), 
                np.random.randint(0,5,(1000,)), 
                np.random.randint(0,10,(1000,))], axis=1)

# mtx = np.random.randint(0,5, (1000,3))
for _ in range(10000):
    time += 1
    gp, dic = rand_merge_group_layers(mtx)
    hash_gp = hash(str(gp))
    if hash_gp in hash_lookup:
        continue
    else:
        hash_lookup.append(hash(str(gp)))
        ct = set(np.bincount(gp))
        if 0 in ct:
            ct.remove(0)
        if np.min(list(ct)) >= MIN_SAMPLE_PER_LEAF:
            success += 1
            print('Success %d in %d' %(success, time))
            print(dic)
            print(np.bincount(gp))

#======= Generate final grouping based on annotation, Keep a lookup table


#======= QC: each group should contrain at least n samples


#======= Objective function: ML ROC, Anoism analysis?

