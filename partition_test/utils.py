#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 17 15:05:24 2021

@author: jefft
"""
import numpy as np

def rand_merge_group_layers(mtx):
    """Random merge M layers of groupping on N samples (N, M). Return (N, 1)
    """
    dics = {}
    cur_group = 0
    fn_group = np.array([-1 for _ in range(mtx.shape[0])])
    for i in range(mtx.shape[1]):
        gps = np.unique(mtx[:,i])
        pt = np.random.choice((0,1), len(gps))  # 0: not consider next layer
        for j in range(len(pt)):
            if pt[j] == 0:
                gp_dic = {}
                coll_tp = set()
                idx = np.where((mtx[:,i] == gps[j]) & (fn_group == -1))[0]
                for t in idx:
                    coll_tp.add(tuple(mtx[t,:(i+1)]))
                for tp in coll_tp:
                    gp_dic[tp] = cur_group
                    cur_group += 1
                for t in idx:
                    fn_group[idx] = gp_dic[tuple(mtx[t,:(i+1)])]
                if gp_dic:
                    dics.update(gp_dic)
            
    gp_dic = {}
    coll_tp = set()
    for i in range(len(fn_group)):
        if fn_group[i] == -1:
            coll_tp.add(tuple(mtx[i,:]))
    for tp in coll_tp:
        gp_dic[tp] = cur_group
        cur_group += 1
    for i in range(len(fn_group)):
        if fn_group[i] == -1:
            fn_group[i] = gp_dic[tuple(mtx[i,:])]
    dics.update(gp_dic)
    return fn_group, dics



class rand_part():
    """Random partition object
    """
    
    def __init__(self, lower, upper, resolution, int_bound=False):
        assert lower < upper
        
        self.lower = lower
        self.upper = upper
        self.resolution = resolution
        self.bounds = [self.lower]
        self.partition = []
        
        ptr = lower
        step = (upper - lower) / resolution 
        while ptr < upper:
            ptr = ptr + step
            if int_bound:
                ptr = int(np.round(ptr))
            self.bounds.append(ptr)
        if ptr < upper:
            self.bounds.append(int(upper) if int_bound else upper)
        assert len(self.bounds) > 2   # at least one partition
        self.refresh()
        return
    
    def refresh(self):
        # merge with next segment or not
        pt = np.random.choice((0,1), len(self.bounds) - 1)
        self.partition = [self.bounds[0]]
        for i in range(1, len(self.bounds)-1):
            if pt[i-1] != 1:
                self.partition.append(self.bounds[i])
        self.partition.append(self.bounds[-1])
        return
    
    def get_idx(self, obj):
        """Get partition index of a numerical object
        """
        assert self.lower <= obj <= self.upper
        
        for i in range(len(self.partition)):
            if self.partition[i] >= obj:
                return np.max((0, i-1))
            

    