# -*- coding: utf-8 -*-
"""
Created on Wed May  1 13:58:51 2019

@author: engelen
"""

import numpy as np

def split_list(ls, splits):
    return([ls[i : j] for i, j in zip([0] + 
              splits, splits + [None])])

def find_cumsum_splits(ls, max_sum):
    split_ls = []
    comb_sum = ls[0]
    for i, p in enumerate(ls[1:]):
        if (comb_sum + p) <= max_sum:
            comb_sum += p
        else:
            split_ls.append((i+1,comb_sum))
            comb_sum=p
    
    return list(zip(*split_ls))

def subdivide_time(t, max_perlen):
    #NOTE THIS FUNCTION DOES NOT WORK IF perlen > max_perlen
    perlens = t[1:] - t[:-1]
    splits, comb_perlens = find_cumsum_splits(perlens, max_perlen)
    
    sub_t = split_list(t, list(splits)+[len(perlens)])
    first_el = np.array([sub[0] for sub in sub_t])
    starts, ends = first_el[:-1], first_el[1:]
    
    sub_ends = ends-starts
    sub_t = [sub-start for sub, start in zip(sub_t, starts)]
    
    #TO DO: Needlessly complex, can better use splits with split.prepend(0) and split.append(len(perlen))
    sub_splits = np.cumsum(np.array([0]+[i.shape[0] for i in sub_t]))
    
    return(sub_t, sub_ends, sub_splits)