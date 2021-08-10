#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 30 15:00:56 2021

@author: liam
"""


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import itertools 


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

def build_TM(states):
    all_states = all_possible_states()
    TM_m = np.zeros((len(all_states),len(all_states)))
    norm_t = []
    
    if np.ndim(states) == 1:
        for si in range(1, len(states)):
            TM_m[int(states[si]),int(states[si-1])] += 1 
    else:
        for S in states:
            for si in range(1, len(states[S])):
                TM_m[int(states[S][si]),int(states[S][si-1])] += 1 
    TM_sym = 1/2 * (TM_m + TM_m.T)
    norm = TM_m.sum(axis=1)
    TM_norm = np.zeros(np.shape(TM_sym))
    for i,j in enumerate(TM_m):
        TM_norm[i] = j / norm[i]
    
    #TM_norm = np.divide(TM_m, norm)
    TM_norm = np.nan_to_num(TM_norm)
    # TM_test = np.nan_to_num(np.divide(TM_m, norm))
    return TM_norm
    
def solve_pi_eq(P):
    import scipy.sparse.linalg as sla
    evals, evecs = sla.eigs(P.T, k = 1, which='LM')
    evecs = np.real(evecs)
    pi_eg = (evecs/evecs.sum()).real
    pi_eg[pi_eg < 10**(-7)] = 0
    pi_EG = np.array([e[0] for e in pi_eg])
    pi_eg = pi_EG
    return pi_eg, np.linalg.eig(P.T)

# create empty lists for the x and y data
# x = []
# y = []

def UNIX_TEXT(pth2fl):
        fl = open (pth2fl,'r')
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
        
        return shell, n1,bins1, states, possible_states
     
def COMMA_SEPERATED(pth2fl,dt=1):
    possible_states = all_possible_states()
    states = pd.read_csv(pth2fl,index_col=0).T   
    TM_norm = build_TM(states.iloc[0,::dt])
    pi_eq_ref, eigs = solve_pi_eq(TM_norm) 
    
    holder = []
    sim_list = len(states.index)
    for nset in range(0,sim_list,20):
        TM_norm = build_TM(states.iloc[nset:nset+10,::dt])
        pi_eq, eigs = solve_pi_eq(TM_norm)
        holder.append(pi_eq) 
    
    return 0, pi_eq_ref,eigs, holder, possible_states

    
def Generate_Animated_Fig(path="test",file="border.txt"):
    
        # function that draws each frame of the animation
    def animate_histogram_unix(i):
        # pt = randint(1,9) # grab a random integer to be the next y-value in the animation
        # x.append(i)
        # y.append(pt)
        states = []
        for s in shell[i*10:i*10+10]: 
                states.append(check_states(s,possible_states))
        n, bins = np.histogram(states, bins=len(possible_states),range=(0,len(possible_states)))
        
        N = pd.DataFrame(n/np.sum(n))
        N_sort = N.sort_values(0,ascending=False)
        tmp_N = np.array(N_sort)
        try:
            tmp_ind = [i for i in range(len(N)) if np.sum(tmp_N[:i])>=.93 and np.sum(tmp_N[:i])<=.97][0]
        except:
            tmp_ind = [i for i in range(len(N)) if np.sum(tmp_N[:i])>=.8 and np.sum(tmp_N[:i])<.99][0]
        inds = list(N_sort.iloc[:tmp_ind].index)
        other_inds = [i for i in N.index if i not in inds]
        N.T[other_inds] = 0
        
        ax.clear()
        ax.bar(bins1[1:], [i[0] for i in N.values])
        ax.bar(bins1[1:],n1,alpha=.5)
        ax.text(80, 0.2, r'$t(%i)$'%(i*10))
        ax.set_xlim([20,100])
        ax.set_ylim([0,.35])

    def animate_histogram_csv(i):
        # pt = randint(1,9) # grab a random integer to be the next y-value in the animation
        # x.append(i)
        # y.append(pt)

        n, bins = np.histogram(states[i:i+1], bins=len(possible_states),range=(0,len(possible_states)))
        #n = n/n.sum()
        ax.clear()
        ax.bar(bins[1:], states[i])
        ax.bar(bins[1:],n1,alpha=.5)
        ax.set_xlim([100,len(bins)])
        ax.set_ylim([0,.75])
    
    
    pth2fl = str(path+"/"+file)
    shell, n1, bins1, states, possible_states = None, None, None, None, None
    
    if any(file[-4:] == dottype for dottype in [".txt",".log",".dat"]):
        fig, ax = plt.subplots()
        shell, n1, bins1, states, possible_states = UNIX_TEXT(pth2fl)
        ani = FuncAnimation(fig, animate_histogram_unix, frames=125, interval=5, repeat=False)
        #plt.show()

    elif file[-4:] == ".csv":
        fig, ax = plt.subplots()
        shell, n1, bins1, states, possible_states = COMMA_SEPERATED(pth2fl)
        ani = FuncAnimation(fig, animate_histogram_csv, frames=5, interval=10, repeat=False)
        #animate_histogram_csv(0)
        #plt.show()
    ani.save('animation_discrete.gif', writer='imagemagick', fps=1)

    # create the figure and axes objects
    # fig, ax = plt.subplots()
    

        
        # run the animation
    # ani = FuncAnimation(fig, animate_histogram, frames=5, interval=5, repeat=False)

    
#states = Generate_Animated_Fig(path="/home/liam/UDel/resources/test_vor/data",file="SLNEW2.csv")
states = Generate_Animated_Fig()