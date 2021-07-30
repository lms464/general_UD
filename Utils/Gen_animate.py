#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 30 15:00:56 2021

@author: liam
"""


from random import randint
import numpy as np
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

# create empty lists for the x and y data
x = []
y = []

fl = open ("test/border.txt")
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

# create the figure and axes objects
fig, ax = plt.subplots()

# function that draws each frame of the animation
def animate(i):
    # pt = randint(1,9) # grab a random integer to be the next y-value in the animation
    # x.append(i)
    # y.append(pt)
    states = []
    for s in shell[i*250:i*250+250]: 
            states.append(check_states(s,possible_states))
    n, bins = np.histogram(states, bins=len(possible_states),range=(0,len(possible_states)))
    ax.clear()
    ax.bar(bins[1:], n/n.sum())
    ax.set_xlim([20,100])
    ax.set_ylim([0,.5])
    
    # run the animation
ani = FuncAnimation(fig, animate, frames=5, interval=2, repeat=False)

plt.show()
