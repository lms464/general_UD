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

def plot_CGTM(pi_eq, TM_norm, leaflet_in, kind):
    states = pd.read_csv("/Censere/UDel/resources/test_vor/data/%s%s.csv"%(leaflet_in,kind),index_col=0).T   

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
    
    plt.pcolormesh(states.T,cmap="ocean_r",vmin=0,vmax=220)
    plt.xlabel("Simulation")
    plt.ylabel("Frame")
    plt.colorbar(label="State")
    plt.savefig("%s_State_Change_Sims%s.pdf"%(leaflet_in,kind))
    plt.close()

def plot_sigConverge(sigSU, sigSL):
    sim_list = [2,4,6,8,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,99]
    plt.plot(sim_list,sigSU,"o--",label="Outer Leaflet")
    plt.plot(sim_list,sigSL,"+--",label="Inner Leaflet")
    plt.ylabel(r"$\sigma_{\pi}$")
    plt.xlabel("Simulations Used")
    plt.legend()
    plt.savefig("Convergence.pdf")
    plt.close()