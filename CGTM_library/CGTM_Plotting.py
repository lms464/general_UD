#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 18 10:28:15 2021

CGTM Plotting Routines

@author: liam
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import all_possible_states as aps
import choose_path as chp
import matplotlib.colors as mcol

class MidpointNormalize(mcol.Normalize):
	"""
	Normalise the colorbar so that diverging bars work there way either side from a prescribed midpoint value)
	e.g. im=ax1.imshow(array, norm=MidpointNormalize(midpoint=0.,vmin=-100, vmax=100))
	"""
	def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
		self.midpoint = midpoint
		mcol.Normalize.__init__(self, vmin, vmax, clip)

	def __call__(self, value, clip=None):
		# I'm ignoring masked values and all kinds of edge cases to make a
		# simple example...
		x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
		return np.ma.masked_array(np.interp(value, x, y), np.isnan(value))

def plot_CGTM(pi_eq, TM_norm, leaflet_in, kind):

    plt.bar(np.arange(0,len(aps.all_possible_states())),pi_eq)
    plt.ylabel("Prob of State")
    plt.xlabel("State")
    plt.savefig("%sPi_Eq_%s.pdf"%(leaflet_in,kind))
    plt.close()
    
    plt.pcolormesh(TM_norm,cmap="gist_earth_r")
    plt.xlabel("State")
    plt.ylabel("State")
    plt.colorbar(label="Prob of Transition")
    plt.savefig("%s_TMCG_%s.pdf"%(leaflet_in,kind))
    plt.close()
    
    try:
        states = pd.read_csv("/Censere/UDel/resources/test_vor/data/states/%s%s.csv"%(leaflet_in,kind),index_col=0).T   
    
        plt.pcolormesh(states.T,cmap="ocean_r",vmin=0,vmax=220)
        plt.xlabel("Simulation")
        plt.ylabel("Frame")
        plt.colorbar(label="State")
        plt.savefig("%s_State_Change_Sims%s.pdf"%(leaflet_in,kind))
        plt.close()
    except:
        pass
def diff_plot(pi,raw, ax):
    dpi = pi - raw
    #fig, ax = plt.subplots()
    line = ax.bar(np.arange(0,len(aps.all_possible_states())), dpi)
    return line


def plot_sigConverge(sigSU, sigSL,kind):
    sim_list = [2,4,6,8,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,99]
    plt.plot(sim_list,sigSU,"o--",label="Outer Leaflet")
    plt.plot(sim_list,sigSL,"+--",label="Inner Leaflet")
    plt.ylabel(r"$\sigma_{\pi}$")
    plt.xlabel("Simulations Used")
    plt.legend()
    plt.savefig("Convergence_%s.pdf"%kind)
    plt.close()

def Ternary_Heat_Map(leaflet_in,leaflet_in2=None):

    import matplotlib.tri as tri

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
    def get_raw(leaflet_in):
        states = pd.read_csv("%s/%s.csv"%(chp.choose_path(),leaflet_in),index_col=0).T   
        return states       
    def get_hist(states):
        hist,edge = np.histogram(states,bins=len(aps.all_possible_states()),normed=None,range=(0,len(aps.all_possible_states())))
        return hist/hist.sum(),edge
    
    def run_pi_eq(state):
        states_r = get_raw(state).values[0]
        states = np.asarray(aps.all_possible_states())
        return  states,states_r 
    
    def run_raw(state):
        ## Should be removed...
        states_r = get_raw(state)
        hist ,edge = get_hist(states_r)
        states = np.asarray(aps.all_possible_states())
        return states, hist
    
    def run_ternary(state):
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
        states, hist = run_pi_eq(state)    
        
        # elif pi_eq == False and raw == False:
        #     print("Please set raw or pi_eq to True")
        #     return 0
        
        STATES = []
        
        for si, s in enumerate(states):
            STATES.append(np.append(s,hist[si]))
        STATES = np.asarray(STATES)
        
        
        #Define twin axis
        fig, ax = plt.subplots()
        # Note that the ordering from start to stop is important for the tick labels
        plot_ticks(right, left, bottom_tick, n, offset=(0, -0.04))
        plot_ticks(left, top, left_tick, n, offset=(-0.06, -0.0))
        plot_ticks(top, right, right_tick, n,offset=(0,.01))
        # fig, tax = ternary.figure(scale=100)
        # fig.set_size_inches(5, 4.5)
        # tax.scatter(pd.DataFrame(states)[[0,1,2]].values)
        # tax.gridlines(multiple=20)
        # tax.get_axes().axis('off')
        
        
        
        a = states[:,0]
        b = states[:,1]
        c = states[:,2]
        
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
        fig.text(0.15, 0.065, 'Chol', fontsize=12, color='black')
        fig.text(0.92, 0.19, 'Neutral', fontsize=12, color='black')
        fig.text(0.40, 0.89, 'Anionic', fontsize=12, color='black')
        
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
        #plt.show()
        plt.savefig("%s_tern.pdf"%state)
        plt.close()
        
    def run_ternary_diff(state1,state2):
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
        hist1,hist2 = 0,0
        states1,states2 = 0,0
        # if raw == True:
        #     states, hist = run_raw(state)    
        
        # elif pi_eq == True:
        states, hist1 = run_pi_eq(state1)    
        states, hist2 = run_pi_eq(state2)  
        hist = hist1-hist2
        norm2 = MidpointNormalize(midpoint=np.min(hist),vmin=0,vmax=np.max(hist))
        # elif pi_eq == False and raw == False:
        #     print("Please set raw or pi_eq to True")
        #     return 0
        
        STATES = []
        
        for si, s in enumerate(states):
            STATES.append(np.append(s,hist[si]))
        STATES = np.asarray(STATES)
        
        
        #Define twin axis
        fig, ax = plt.subplots()
        # Note that the ordering from start to stop is important for the tick labels
        plot_ticks(right, left, bottom_tick, n, offset=(0, -0.04))
        plot_ticks(left, top, left_tick, n, offset=(-0.06, -0.0))
        plot_ticks(top, right, right_tick, n,offset=(0,.01))
        # fig, tax = ternary.figure(scale=100)
        # fig.set_size_inches(5, 4.5)
        # tax.scatter(pd.DataFrame(states)[[0,1,2]].values)
        # tax.gridlines(multiple=20)
        # tax.get_axes().axis('off')
        
        
        
        a = states[:,0]
        b = states[:,1]
        c = states[:,2]
        
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
        plt.tricontourf(x,y,T.triangles,v,cmap='RdBu_r')
        
        
        # create the grid
        corners = np.array([[0, 0], [1, 0], [0.5,  np.sqrt(2.3)*0.576]])
        triangle = tri.Triangulation(corners[:, 0], corners[:, 1])
        
        # creating the grid
        refiner = tri.UniformTriRefiner(triangle)
        trimesh = refiner.refine_triangulation(subdiv=3)
        
        # plt.axis('off')
        
        
        
        #plotting the mesh and caliberate the axis
        plt.triplot(trimesh,'k--',normalize=norm2)
        #plt.title('Binding energy peratom of Al-Ti-Ni clusters')
        # ax.set_xlabel('Al-Ti',fontsize=12,color='black')
        # ax.set_ylabel('Ti-Ni',fontsize=12,color='black')
        # ax2 = ax.twinx()
        # ax2.set_ylabel('Al-Ni',fontsize=12,color='black')
        # Corners
        fig.text(0.15, 0.065, 'Chol', fontsize=12, color='black')
        fig.text(0.92, 0.19, 'Neutral', fontsize=12, color='black')
        fig.text(0.40, 0.89, 'Anionic', fontsize=12, color='black')
        
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
        # plt.show()
        plt.savefig("%s_tern_diff.pdf"%state1)
        plt.close()    
    
    
    if leaflet_in2 is not None :
        run_ternary_diff(leaflet_in,leaflet_in2)
    else:
        run_ternary(leaflet_in)

        

def Ternary_Scatter(kind, data1, data2=None, rot=None):
    
    import ternary
    
    if rot is not None:
        data1 = np.array([data1[:,1],data1[:,2],data1[:,0]]).T
        
    
    figure, tax = ternary.figure(scale=1)
    #figure.set_size_inches(10, 7.5)
    tax.scatter(data1)
    if data2 is not None:
        tax.scatter(data2)
    tax.boundary(linewidth=2.0)
    tax.gridlines(multiple=1, color="blue")
    tax.ticks(axis='lbr', linewidth=.5, multiple=1)
    tax.clear_matplotlib_ticks()
    tax.get_axes().axis('off')
    # plt.show()
    tax.savefig("Scatter_%s.pdf"%kind)
    tax.close()
    
# Ternary_Heat_Map("pi_eq_SUChainsT")
# Ternary_Heat_Map("pi_eq_SLChainsT")