#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 22 12:57:40 2021

@author: sharplm
"""
from celluloid import Camera
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import itertools 
from scipy.stats import norm
from scipy.stats import f_oneway
from scipy.stats import ttest_ind
####
'''
MAKE the prob a polynomial again!
'''
####

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
 
#### Ternary stuff #####
def get_raw(leaflet_in):
    states = pd.read_csv("%s.txt"%leaflet_in).T   
    return states
def get_hist(states):
    hist,edge = np.histogram(states,bins=len(all_possible_states()),normed=None,range=(0,len(all_possible_states())))
    return hist/hist.sum(),edge

def plot_ticks(start, stop, tick, n, offset=(.0, .0)):
    r = np.linspace(0, 1, n+1)
    x = start[0] * (1 - r) + stop[0] * r
    x = np.vstack((x, x + tick[0]))
    y = start[1] * (1 - r) + stop[1] * r
    y = np.vstack((y, y + tick[1]))
    plt.plot(x, y, 'k', lw=1)
    
    # add tick labels
    for xx, yy, rr in zip(x[1], y[1], r):
        plt.text(xx+offset[0], yy+offset[1], "{:.2}".format(rr))

def run_raw(state):
    states_r = get_raw(state)
    hist ,edge = get_hist(states_r)
    states = np.asarray(all_possible_states())
    return states, hist

def run_pi_eq(state):
    states_r = get_raw(state)#.values[0]
    states = np.asarray(all_possible_states())
    return  states,states_r

#### Runners ####  
def animate_hist():
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
    # mean,std = norm.fit(states[:])
    # x = np.linspace(np.min(states[:]),np.max(states[:]),500)
    # y = norm.pdf(x,mean,std)
    
    # plt.plot(x,y)
    # plt.bar(bins1[:-1],n1)
    # plt.show()
    # plt.close()
    
    fig = plt.figure()
    camera = Camera(fig)
    # plt.plot(bins[:-1],state_dist)
    # plt.show()
    
    for i in range(0,5000,1250):
        states= []
        for s in shell[i:i+99]: 
            states.append(check_states(s,possible_states))
        n, bins = np.histogram(states, bins=len(possible_states),range=(0,len(possible_states)))
        plt.bar(bins[:-1],n/np.sum(n),color='blue')
        plt.bar(bins1[:-1],n1,color='orange',alpha=0.75)
        print(f_oneway(n1,n)[1],ttest_ind(n1, n)[1])
        plt.text(.75,.75, f_oneway(n1,n)[1])
        plt.savefig("hist_%i.png"%i)
        plt.close()
    #     camera.snap()
    # animation = camera.animate()
    # animation.save('test_hist.gif',writer="imagemagick")

def animate_mixture(state, raw=False, pi_eq=False):
    import ternary
    import matplotlib.tri as tri

    n = 4
    tick_size = 0.1
    margin = 0.05
    
    # define corners of triangle    
    left = np.r_[0, 0]
    right = np.r_[1, 0]
    top = np.r_[0.5,  np.sqrt(2.3)*0.576]
    triangle = np.c_[left, right, top, left]
    
    # define vectors for ticks
    bottom_tick = 0.8264*tick_size * (right - top) / n
    right_tick = 0.8264*tick_size * (top - left) / n
    left_tick = 0.8264*tick_size * (left - right) / n
    
    # state = "pi_sl"#"StatesL_I"
    hist = 0
    states = 0
    # if raw == True:
    #     states, hist = run_raw(state)    
    
    # elif pi_eq == True:
    #     states, hist = run_pi_eq(state)    
    
    # elif pi_eq == False and raw == False:
    #     print("Please set raw or pi_eq to True")
    #     return 0
    
    STATES = []
    
    # for si, s in enumerate(states):
    #     STATES.append(np.append(s,hist[si]))
    # STATES = np.asarray(STATES)
    
    
    #Define twin axis
    fig, ax = plt.subplots()
    camera = Camera(fig)
    # Note that the ordering from start to stop is important for the tick labels
    plot_ticks(right, left, bottom_tick, n, offset=(0, -0.04))
    plot_ticks(left, top, left_tick, n, offset=(-0.06, -0.0))
    plot_ticks(top, right, right_tick, n,offset=(0,.01))
    # fig, tax = ternary.figure(scale=100)
    # fig.set_size_inches(5, 4.5)
    # tax.scatter(pd.DataFrame(states)[[0,1,2]].values)
    # tax.gridlines(multiple=20)
    # tax.get_axes().axis('off')
    
    all_states, state_dat = run_pi_eq(state)    

    a = all_states[:,0]
    b = all_states[:,1]
    c = all_states[:,2]
    tau = 1250
    for i in range(0,5000,tau): 
        #elif pi_eq == True:
        hist, edge = get_hist(state_dat.iloc[0,i:i+(tau - 1)])

        
        # # values is stored in the last column
        v = hist
        # t = np.transpose(np.array([[0,0],[1,0],[0,1]]))
        # X,Y = [], []
        # for s in states:
        #     R = t.dot(s)
        #     X.append(R[0]), Y.append(R[1])
        
        # # translate the data to cartesian corrds
        x = 0.5 * ( 2.*b+c ) / ( a+b+c )
        y = 0.5*np.sqrt(3) * c / (a+b+c)
        
        
        # # create a triangulation out of these points
        T = tri.Triangulation(x,y)
        
        # # plot the contour
        print(i)
        plt.tricontourf(x,y,T.triangles,v,cmap='RdBu_r')
        
        
        # create the grid
        corners = np.array([[0, 0], [1, 0], [0.5,  np.sqrt(2.3)*0.576]])
        triangle = tri.Triangulation(corners[:, 0], corners[:, 1])
        
        # creating the grid
        refiner = tri.UniformTriRefiner(triangle)
        trimesh = refiner.refine_triangulation(subdiv=3)
        
        # plt.axis('off')
        
        
        
        #plotting the mesh and caliberate the axis
        plt.triplot(trimesh,'k--')
        #plt.title('Binding energy peratom of Al-Ti-Ni clusters')
        # ax.set_xlabel('Al-Ti',fontsize=12,color='black')
        # ax.set_ylabel('Ti-Ni',fontsize=12,color='black')
        # ax2 = ax.twinx()
        # ax2.set_ylabel('Al-Ni',fontsize=12,color='black')
        # Corners
        fig.text(0.07, 0.05, 'Anionic', fontsize=12, color='black')
        fig.text(0.93, 0.05, 'Chol', fontsize=12, color='black')
        fig.text(0.50, 0.90, 'Neutral', fontsize=12, color='black')
        
        # Connections
        # fig.text(0.47, 0.05, 'Ti-Al', fontsize=12, color='black')  # Note: not sure about
        # fig.text(0.72, 0.50, 'Al-Ni', fontsize=12, color='black')  # the nomenclature;
        # fig.text(0.25, 0.50, 'Ti-Ni', fontsize=12, color='black')  # might be switched
        
        
        #set scale for axis
        # ax.set_xlim(1, 0)
        # ax.set_ylim(0, 1)
        # ax2.set_ylim(1, 0)
        ax.set_axis_off()
        cax = plt.axes([0.75, 0.55, 0.055, 0.3])
        plt.colorbar(cax=cax,format='%.3f')
        plt.savefig("tern_bin%i.png"%i)
        plt.close()
        #animation = camera.animate()
    #animation.save('test_hist.gif',writer="imagemagick")
        # plt.show()
    # plt.savefig("%s_tern.pdf"%state)
    # plt.close()
# animate_mixture("states", False, True)
animate_hist()