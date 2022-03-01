#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 26 11:09:44 2022

@author: jefft
"""

import pandas as pd
import numpy as np

maf_fp = '/Users/jefft/Desktop/p53_project/datasets/BRCA-TCGA/brca_tcga_pan_can_atlas_2018/data_mutations.txt'
mut_pools = maf['HGVSp']


maf = pd.read_table(maf_fp, comment = '#', sep = '\t', header=0)
    
print(count)