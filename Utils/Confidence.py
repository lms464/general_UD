#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  6 13:39:55 2021

@author: liam
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import statsmodels.stats.api as sms
import itertools 
from statsmodels.stats.weightstats import ztest

def all_possible_states():
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

def check_states(data_frm, possible_states):
    
    data_frm = np.asarray(data_frm)/np.sum(data_frm)
    holder = []
    for pi, ps in enumerate(possible_states):
        
        rmsd = np.sqrt((data_frm - ps)**2)
        holder.append(rmsd.sum())

    out = np.argmin(holder)
    
    return out


fl = open ("/home/liam/lms464/Active3/borders.txt",'r')
lines = fl.readlines()#.split()
fl.close()
#lines = [int(l) for l in lines]
shell = []
for l in lines[1::3]:
	shell.append([int(l.split()[-3]),int(l.split()[-2]),int(l.split()[-1])])
shell_arr = np.array(shell)

possible_states = all_possible_states()

states = []

for s in shell:
	states.append(check_states(s,possible_states))
states = np.array(states)
n1, bins1 = np.histogram(states, bins=len(possible_states),range=(0,len(possible_states)))
n1 = n1/np.sum(n1)
Z = []
Pat = []
Rat = []
tau = 1
for i in range(0,len(shell),tau):
    print(i)
    states = []
    for s in shell[:i*tau+tau]: 
         states.append(check_states(s,possible_states))
    n, bins = np.histogram(states, bins=len(possible_states),range=(0,len(possible_states)))
    # tmp = sms.DescrStatsW(n).tconfint_mean()
    # N = pd.DataFrame(n/np.sum(n))
    n = n/np.sum(n)
    # N_sort = N.sort_values(0,ascending=False)
    # tmp_N = np.array(N_sort)
    # try:
    #     tmp_ind = [i for i in range(len(N)) if np.sum(tmp_N[:i])>=.945 and np.sum(tmp_N[:i])<=.965][0]
    # except:
    #     tmp_ind = [i for i in range(len(N)) if np.sum(tmp_N[:i])>=.8 and np.sum(tmp_N[:i])<.99][0]
    # inds = list(N_sort.iloc[:tmp_ind].index)
    # other_inds = [i for i in N.index if i not in inds]
    # N = N.T
    # N[other_inds] = 0
    # N = N.T.values
    Pat.append(n)
    ratio= n/n1
    ratio=np.nan_to_num(ratio,0)
    Rat.append(ratio)
    ratio_sum = np.sum(ratio)/len(n1[n1>0])
    Z.append(ratio_sum)
# !    print()
    
    #z,p = ztest(N,n1)
    #Z.append(z)
    #P.append(p)
np.savetxt("Active_cumulative.dat",Z)
plt.plot(np.linspace(0,10,len(Z)),Z,'o')
plt.savefig("Active.pdf")
plt.close()
    # print(tmp)
    #plt.plot(i,tmp[0]/0.00259,"bo--")
