#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 11:14:08 2022

@author: liam
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
import all_possible_states as aps

try:
    sys.path.insert(0, '/home/liam/lms464/github/cg_msm_reweighting')
    import traj_analysis as ta

except:
    sys.path.insert(0, '/home/sharplm/github/cg_msm_reweighting')
    import traj_analysis as ta

def running_mean(x, N):
    cumsum = np.cumsum(np.insert(x, 0, 0)) 
    return (cumsum[N:] - cumsum[:-N]) / float(N)

def unmap(remap_dict, data, nstates):
    out = np.zeros(nstates)
    for key in remap_dict:
        out[key] = data[remap_dict[key]]
    return out

def cummulative_difference(pi_ref,pi_eq):
    diff_list = []
    for ref,eq in zip(pi_ref,pi_eq):
        diff_list.append(np.abs(ref - eq))
    return np.sum(diff_list)

def power_law(x, a, b):
    return a*np.power(x, b)

def lag_time_analysis(lagmin,lagmax,lim=200,eig=False):
    import matplotlib.pyplot as plt
    data = pd.read_csv("/home/liam/lms464/CG/data/states/low_short_inactive_tmpNEW3.csv",index_col=0)#low_long_inactive_tmpNEW.csv",index_col=0)#"~/CG/data/states/low_short_inactive_tmp.csv",index_col=0)
    dat_ind = ["%i"%i for i in range(int(len(data.T)/2),len(data.T))]
    dat_ind = ["0"] + dat_ind
    trj,n_states,remp = ta.remap_trajs(data[dat_ind].values)
    hist = np.histogram(trj,bins=n_states,range=(0,n_states))[0]
    hist = hist/hist.sum()
    true_hist = pd.read_csv("/home/liam/lms464/CG/data/pi_eq_low_inactive_long_cgtmp_231.csv",index_col=0).values
    true_hist = true_hist.T[0]

    its = 1
    stps = 1
    out2 = []
    cum_diff_Long = []
    cum_diff_Short = []
    for tau in np.arange(lagmin,lagmax,1):
        dat_num = np.shape(trj[:,::tau])[1]*np.shape(trj[:,::tau])[0]
        print(r"tau = %i, data:%i"%(tau,dat_num))
        if dat_num < lim:
            lagmax = tau
            break
        reweighted_distributions, reweighted_matrices, eign = ta.optimized_resliced_voelz(
            trj[:,::tau],
            its,
            stps,
            n_states,
            return_eig=True,
            return_matrices=True
        )
        if eig==True:
            out2.append(eign)
        reweight_mean = np.mean(reweighted_distributions,axis=0)
        reweight_map = unmap(remp, reweight_mean, 231)
        hist_map = unmap(remp,hist,231)

        cum_diff_Short.append(cummulative_difference(hist_map,reweight_map))
        cum_diff_Long.append(cummulative_difference(true_hist,reweight_map))
        
    if eig==True:
        np.savetxt("eig_OG2.dat",out2)
        diff = np.array(cum_diff_Long)-np.array(cum_diff_Short)
        plt.scatter(np.arange(lagmin,lagmax,1),cum_diff_Short)
        plt.scatter(np.arange(lagmin,lagmax,1),cum_diff_Long)
        plt.scatter(np.arange(lagmin,lagmax,1),diff)
        plt.show()
    else:
        diff = np.array(cum_diff_Long)-np.array(cum_diff_Short)
        plt.scatter(np.arange(lagmin,lagmax,1),cum_diff_Short)
        plt.scatter(np.arange(lagmin,lagmax,1),cum_diff_Long)
        plt.scatter(np.arange(lagmin,lagmax,1),diff)
        plt.show()
        np.savetxt("CGTM_Lag_OG2.dat",np.array([np.arange(lagmin,lagmax,1),cum_diff_Short,cum_diff_Long]).T)
        return np.where(diff == diff[diff>np.min(cum_diff_Short)].min())[0][0]+lagmin
def analyze_eigen():

    def smooth(y, box_pts):
        box = np.ones(box_pts)/box_pts
        y_smooth = np.convolve(y, box, mode='same')
        return y_smooth    

    eigO = np.loadtxt("eig_OG2.dat")
    dt = 1 #ns
    lag_time = np.arange(1,len(eigO)+1)
    tau_O = -1 * lag_time/np.log(eigO)
    count = 0
    test = True
    while count < 11:
        print(count)
        #manipulate to get a usable threshold
        tau_manip = tau_O.copy()
        tau_manip[np.isneginf(tau_manip)] = 0
        tau_means = np.nanmean(tau_manip)
        tau_thresh = tau_means# - tau_means*.333333
        if count == 1:
            tau_thresh = np.nanmean(tau_O)
        elif count > 1:
            tau_thresh = tau_means + tau_means*count/10
        lag_time_manip = lag_time.copy()
        del_ind = np.where(tau_manip==0)[0]
        tau_manip = np.delete(tau_manip,del_ind)
        lag_time_manip = np.delete(lag_time_manip,del_ind)
        del_ind = np.where(tau_manip>=tau_thresh)[0]
        tau_manip = np.delete(tau_manip,del_ind)
        lag_time_manip = np.delete(lag_time_manip,del_ind)
        if len(tau_manip)/len(tau_O)<=.6:
            count = count + 1
        else:
            count = 11
    
    # find the ""flatests sloap""
    m_list = []
    for i in range(1,len(tau_manip)):
        m_list.append(np.log((tau_manip[i]+tau_manip[i-1])/2) )
    mean_m = np.mean(m_list[int(len(m_list)*.25)+3:])
    m_list = np.array(m_list)
    
    out = np.where(np.isclose(m_list,mean_m,rtol=1e-01, atol=1e-03))[0]

    plt.plot((lag_time),np.log(tau_O),"o")
    plt.plot((lag_time_manip),np.log(tau_manip))
    plt.show()
    return [out.min(),out.max()]

def Get_Lag_data(tau,stps):
    data = pd.read_csv("/home/liam/lms464/CG/data/states/low_short_inactive_tmpNEW3.csv",index_col=0)#low_long_inactive_tmpNEW.csv",index_col=0)#"~/CG/data/states/low_short_inactive_tmp.csv",index_col=0)
    dat_ind = ["%i"%i for i in range(int(len(data.T)/2),len(data.T))]
    dat_ind = ["0"] + dat_ind
    trj,n_states,remp = ta.remap_trajs(data[dat_ind].values)
    hist = np.histogram(trj,bins=n_states,range=(0,n_states))[0]
    hist = hist/hist.sum()
    true_hist = pd.read_csv("/home/liam/lms464/CG/data/pi_eq_low_inactive_long_cgtmp_231.csv",index_col=0).values
    true_hist = true_hist.T[0]
    
    its = 300
    # stps = 1
    reweighted_distributions, reweighted_matrices, eign = ta.optimized_resliced_voelz(
        trj[:,::tau],
        its,
        stps,
        n_states,
        return_eig=True,
        return_matrices=True
    )
    reweight_mean = np.mean(reweighted_distributions,axis=0)
    alps = aps.all_possible_states()
    reweight_map = unmap(remp, reweight_mean, 231)
    hist_map = unmap(remp,hist,231)
    print(cummulative_difference(hist_map,reweight_map))
    print(cummulative_difference(true_hist,reweight_map))
    plt.ylim(-0.06,0.06)
    plt.title(r"$\tau$ %i"%tau)
    plt.bar(np.arange(0,231),hist_map-reweight_map)
    plt.bar(np.arange(0,231),true_hist-reweight_map,alpha=.5)
    plt.savefig("Tau_diff_%i.png"%tau)
    plt.close()
    #plt.show()


lag_time_analysis(1,200,lim=400,eig=True)

tau_min,tau_max = analyze_eigen()
lag = lag_time_analysis(tau_min,tau_max)
print(lag)
for t in range(tau_min,tau_max+1):
    Get_Lag_data(t,1)