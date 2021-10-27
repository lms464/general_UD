#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 17 17:25:35 2021

@author: liam
"""

import numpy as np

def all_possible_states():
    import itertools 

    # RUNS, seems fine
    #sorts by anionic, then neutral, then chol
    list_o_list = [np.linspace(0,1,21),np.linspace(0,1,21),np.linspace(0,1,21)]
    all_possibility = list(itertools.product(*list_o_list))
        
    list_100 = []
     
    for li, l in enumerate(all_possibility):
        if np.isclose(np.sum(l),1.0,1e-3,1e-5) == False:
            # print(li,np.sum(l))
            continue
        list_100.append(l)
    list_100 = np.array([*list_100])
    #return list_100[np.argsort(list_100[:,2])][::-1]
    return list_100[np.lexsort((list_100[:,0], list_100[:,1], list_100[:,2]))][::-1]