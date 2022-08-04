#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 28 13:54:02 2022

JD Current Method

@author: liam
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
import all_possible_states as aps

def unmap(remap_dict, data, nstates):
    out = np.zeros(nstates)
    for key in remap_dict:
        out[key] = data[remap_dict[key]]
    return out
        

try:
    sys.path.insert(0, '/home/liam/lms464/github/cg_msm_reweighting')
    import traj_analysis as ta

except:
    sys.path.insert(0, '/home/sharplm/github/cg_msm_reweighting')
    import traj_analysis as ta

def KL(a, b):
    a = np.asarray(a, dtype=np.float)
    b = np.asarray(b, dtype=np.float)

    return np.sum(np.where(a != 0, a * np.log(a / b), 0))

def Find_Window(inpt):
    data = np.array(pd.read_csv("~/CG/data/states/%s.csv"%inpt,index_col=0).values)
    trj,n_states,remp = ta.remap_trajs(data)
    hist = np.histogram(trj,bins=n_states,range=(0,n_states))[0]
    hist = hist/hist.sum()
    
    list_dist = []
    out = []

    # weight = counts / counts.sum()
    # for di, df in enumerate(steps):
    its = 50
    stps = np.arange(20,40)#1,120,1)#,50,75,100
        # try:
    out2 = []
    for stp in stps:        
        reweighted_distributions, reweighted_matrices = ta.optimized_resliced_voelz(
            trj,
            its,
            stp,
            n_states,
            # _initial_weights=weight,
            return_matrices=True
        )
        reweight_mean = np.mean(reweighted_distributions,axis=0)
        out.append(np.sqrt((hist-reweight_mean)**2).sum())
    plt.plot(np.linspace(1,len(out),len(out)),out)
    # plt.savefig("Window_size_out.pdf")
    # plt.close()
    return out
# oot = Find_Window('low_short_inactive_titrate_tmp')


def Find_WindowKL(inpt):
    data = pd.read_csv("~/lms464/CG/data/states/%s.csv"%inpt,index_col=0)
    dat_ind = ["%i"%i for i in range(int(len(data.T)/2),len(data.T))]
    dat_ind = ["0"] + dat_ind
    trj,n_states,remp = ta.remap_trajs(data[dat_ind].values)
    # trj,n_states,remp = ta.remap_trajs(data)
    hist = np.histogram(trj,bins=n_states,range=(0,n_states))[0]
    hist = hist/hist.sum()
    true_hist = pd.read_csv("~/CG/data/pi_eq_low_inactive_long_cgtmp_231.csv",index_col=0).values
    true_hist = true_hist.T[0]
    out = []

    # weight = counts / counts.sum()
    # for di, df in enumerate(steps):
    its = 150
    stps = np.arange(1,50)#1,120,1)#,50,75,100
        # try:
    out2 = []
    for stp in stps:        
        reweighted_distributions, reweighted_matrices = ta.optimized_resliced_voelz(
            trj,
            its,
            stp,
            n_states,
            # _initial_weights=weight,
            return_matrices=True
        )
        reweight_mean = np.mean(reweighted_distributions,axis=0)
        reweight_mean_re = unmap(remp,reweight_mean,231)
        out.append(KL(reweight_mean,hist))#np.sqrt((hist-reweight_mean)**2).sum())
    plt.plot(np.linspace(1,len(out),len(out)),out)
    plt.show()
    # plt.savefig("Window_size_out.pdf")
    # plt.close()
    return out
# oot = Find_WindowKL('low_short_inactive_tmp')

def Find_Iterations():
    lag_out = []
    lag_ids = []
    for lags in np.arange(40,51): 
        data = pd.read_csv("/home/sharplm/CG/data/states/low_short_inactive_tmpNEW.csv",index_col=0,header=0)
        # dat_ind = data.columns[int(len(data.columns)/2):].tolist()  #["%i"%i for i in range(int(len(data.T)/2),len(data.T))]
        # dat_ind = ["0"] + dat_ind
        # trj,n_states,remp = ta.remap_trajs(data.values)#[dat_ind].values)
        # trj,n_states,remp = ta.remap_trajs(data)
        dat_ind = ["%i"%i for i in range(int(len(data.T)/2),len(data.T))]
        dat_ind = ["0"] + dat_ind
        trj,n_states,remp = ta.remap_trajs(data[dat_ind].values)
        hist = np.histogram(trj,bins=n_states,range=(0,n_states))[0]
        hist = hist/hist.sum()
        true_hist = pd.read_csv("/home/sharplm/CG/data/pi_raw_low_inactive_long_cgtmp.csv",index_col=0).values
        true_hist = true_hist.T[0]
        
        list_dist = []
        out = []
    
        # weight = counts / counts.sum()
        # for di, df in enumerate(steps):
        its = np.arange(100,151)#,50,75,100
        stps = np.arange(1,100)#int(2*len(data.columns)/3)-3)
        if len(stps) < 3:
            print("not enought usable windows")
            continue
            # try:
        for it in its:    
            out_iter = []
    
            for stp in stps:#tri in range(3,len(trj[1:])):
                try:
                    reweighted_distributions, reweighted_matrices = ta.optimized_resliced_voelz(
                        trj[:,::lags],
                        it,
                        stp,
                        n_states,
                        # _initial_weights=weight,
                        return_matrices=True
                    )
                    lag_ids.append("tau:%i-iter:%i-window:%i"%(lags,it,stp))
                    reweight_mean = np.mean(reweighted_distributions,axis=0)
                    reweight_mean_re = unmap(remp,reweight_mean,231)
                    out_iter.append(KL(reweight_mean,hist))#np.sqrt((hist-reweight_mean)**2).sum())
                except:
                    continue
            out.append(out_iter)
        pd.DataFrame(out).to_csv("KL_iteration_windows_lag_allframes_maps%i.csv"%lags)
        lag_out.append(out)
        # plt.plot(np.linspace(1,len(out),len(out)),out)
        # plt.show()
    # return lag_out,lag_ids
# Find_Iterations()#'low_short_inactive_tmp')
# pd.DataFrame(oot).to_csv("KL_iteration_windows_lag_map2.csv")
# pool = mpi.Pool()
# resutls = pool.map(Find_Window(), Find_Iterations())
# cul_dif_window_titrate = Find_Window("low_short_inactive_titrate")
# cul_dif_window = Find_Window("low_short_inactive")

# cul_dif_iter = Find_Iterations()

def test_block():
    data = pd.read_csv("/home/liam/lms464/CG/data/states/low_short_inactive_tmpNEW.csv",index_col=0)#"~/CG/data/states/low_short_inactive_tmp.csv",index_col=0)
    dat_ind = ["%i"%i for i in range(int(len(data.T)/2),len(data.T))]
    dat_ind = ["0"] + dat_ind
    trj,n_states,remp = ta.remap_trajs(data[dat_ind].values)
    hist = np.histogram(trj,bins=n_states,range=(0,n_states))[0]
    hist = hist/hist.sum()
    # true_hist = pd.read_csv("/home/liam/lms464/github/general_UD/CGTM_library/pi_raw_low_inactive_long_cgtmp.csv",index_col=0).values
    # true_hist = true_hist.T[0]

    
    # pd.DataFrame(hist).to_csv("~/CG/data/pi_raw_low_inactiv_titratee_short_cgtmp.csv")
    # list_dist = []
    # out = []
    
    # weight = counts / counts.sum()
    # for di, df in enumerate(steps):
    # states,counts = np.unique(trj[:int(1249-200),:], return_counts=True)    
    its = 1
    # its = np.arange(1,120)#1,120,1)#,50,75,100
    stps = 1
        # try:
    out2 = []
    # for it in its: 
    kl_list = []
    # for tau in np.arange(1,201,1):
    reweighted_distributions, reweighted_matrices, eign = ta.optimized_resliced_voelz(
        trj[:,::45],
        its,
        stps,
        n_states,
        return_eig=True,
        # _initial_weights=weight,
        return_matrices=True
    )
    out2.append(eign)
    reweight_mean = np.mean(reweighted_distributions,axis=0)
    kl_list.append(KL(reweight_mean,hist))
        # for rw_i in range(1,len(reweighted_distributions)-1):
        #     if np.allclose(reweighted_distributions[rw_i],reweighted_distributions[rw_i-1])==True:
        #         out2.append(rw_i)
        #         break
            
        #     out.append(np.sqrt((hist-reweight_mean)**2).sum())
        # plt.plot(np.arange(1,120),out)
        # plt.ylabel("Cumulative Difference")
        # plt.xlabel("Iterations")
        # plt.show()
        # reweight_map = unmap(remp, reweight_mean, 231)
        # pd.DataFrame(reweight_mean).to_csv("~/CG/data/pi_eq_low_inactive_titrate_short_cgtmp.csv")
        # print(tau,np.sqrt(np.diff([reweight_mean,hist])**2).sum())
        # plt.bar(np.linspace(0,n_states,n_states), reweight_mean )
        # plt.bar(np.linspace(0,n_states,n_states), hist,alpha=.5 )
        # plt.title("tau=%i"%tau)
        # plt.savefig("tau_%i.pdf"%tau)
        # plt.close()
    # print( np.sqrt((hist-reweight_mean)**2).sum())
    pd.DataFrame(reweight_mean).to_csv("~/lms464/CG/data/pi_eq_low_inactive_short_cgtmpNEW.csv")
    
    # alps = aps.all_possible_states()
    # alps_re_map = []
    # kl = KL(reweight_mean, hist)
    # for key in remp.keys():
    #     alps_re_map.append(alps[key])
    # alps = pd.DataFrame(alps_re_map)
    # alps['3'] = reweight_mean
    # # alps.to_csv("~/CG/data/pi_eq_low_inactive_titrate_short_map.csv")
    # alps['3'] = hist
    # # alps.to_csv("~/CG/data/pi_raw_low_inactive_titrate_short_map.csv")
    # # plt.savefig("Iterations.pdf")
    # # plt.close()
    # return reweighted_matrices
    return out2,kl_list
# test_block()
def test_block_50_per():
    data = pd.read_csv("~/CG/data/states/low_short_inactive_tmp.csv",index_col=0)
    dat_ind = ["%i"%i for i in range(int(len(data.T)/2),len(data.T))]
    dat_ind = ["0"] + dat_ind
    trj,n_states,remp = ta.remap_trajs(data[dat_ind].values)
    hist = np.histogram(trj,bins=n_states,range=(0,n_states))[0]
    hist = hist/hist.sum()
    
    # pd.DataFrame(hist).to_csv("~/CG/data/pi_raw_low_inactiv_titratee_short_cgtmp.csv")
    # list_dist = []
    # out = []
    
    # weight = counts / counts.sum()
    # for di, df in enumerate(steps):
    # states,counts = np.unique(trj[:int(1249-200),:], return_counts=True)    
    its = 300
    # its = np.arange(1,120)#1,120,1)#,50,75,100
    stps = np.arange(1,62,1)
        # try:
    # out2 = []
    # for it in its: 
    for stp in stps:
        reweighted_distributions, reweighted_matrices = ta.optimized_resliced_voelz(
            trj[:,::1],
            its,
            stp,
            n_states,
            return_eig=False,
            # _initial_weights=weight,
            return_matrices=True
        )
        reweight_mean = np.mean(reweighted_distributions,axis=0)
        # for rw_i in range(1,len(reweighted_distributions)-1):
        #     if np.allclose(reweighted_distributions[rw_i],reweighted_distributions[rw_i-1])==True:
        #         out2.append(rw_i)
        #         break
            
        #     out.append(np.sqrt((hist-reweight_mean)**2).sum())
        # plt.plot(np.arange(1,120),out)
        # plt.ylabel("Cumulative Difference")
        # plt.xlabel("Iterations")
        # plt.show()
        reweight_map = unmap(remp, reweight_mean, 231)
        # pd.DataFrame(reweight_mean).to_csv("~/CG/data/pi_eq_low_inactive_titrate_short_cgtmp.csv")
        plt.bar(np.linspace(0,n_states,n_states), reweight_mean )
        plt.bar(np.linspace(0,n_states,n_states), hist, alpha=0.5)
        plt.title("window: %i"%stp)
        plt.show()
        print( np.sqrt((hist-reweight_mean)**2).sum())
        # pd.DataFrame(reweight_mean).to_csv("~/CG/data/pi_eq_low_inactive_long_cgtmp.csv")
        alps = aps.all_possible_states()
        alps_re_map = []
        kl = KL(reweight_mean, hist)
        for key in remp.keys():
            alps_re_map.append(alps[key])
        alps = pd.DataFrame(alps_re_map)
        alps['3'] = reweight_mean
        # alps.to_csv("~/CG/data/pi_eq_low_inactive_titrate_short_map.csv")
        alps['3'] = hist
    # alps.to_csv("~/CG/data/pi_raw_low_inactive_titrate_short_map.csv")
    # plt.savefig("Iterations.pdf")
    # plt.close()
    return reweighted_matrices

# eig = test_block()
# kl = eig[1]
# eign = [e[0] for e in eig[0]]
# lag_analysis = []
# for i,j in enumerate(eign):
#     lag_analysis.append(-i/np.log(j))
# plt.scatter(np.arange(1,len(eign)+1),np.log(lag_analysis))
# plt.savefig("lag_relax.pdf")
# plt.close()
# test_block_50_per()
# fig, ax = plt.subplots(len(matrx),1,figsize=(5,5*len(matrx)))
# for mi, m in enumerate(matrx):
#     ax[mi].pcolor(m)
# plt.savefig("Mat.pdf")
# plt.close()

def weighted_avg(pi_eq):
    all_states = aps.all_possible_states()
    tmp_wsa = [np.sum((pi_eq*all_states[:,0])),np.sum((pi_eq*all_states[:,1])),np.sum((pi_eq*all_states[:,2]))]
    return tmp_wsa

def get_comp_comparison(inpt):
    data = np.array(pd.read_csv("~/CG/data/states/%s.csv"%inpt,index_col=0).values)
    truth = np.array(pd.read_csv("~/CG/data/pi_raw_low_inactive_long_cg.csv",index_col=0).values).T
    truth_max_ind = np.where(truth[0].max()==truth[0])
    max_truth = aps.all_possible_states()[truth_max_ind]
    trj,n_states,remp = ta.remap_trajs(data)
    # hist = np.histogram(truth,bins=231,range=(0,231))[0]
    # hist = hist/hist.sum()
    truth_weighted = weighted_avg(truth)
    list_dist = []
    out = []
    out2 = []
    # weight = counts / counts.sum()
    # for di, df in enumerate(steps):
    its = 50
    stps = 22
        # try:
    out2 = []
    for di in range(2,len(trj)):       
        reweighted_distributions, reweighted_matrices = ta.optimized_resliced_voelz(
            trj[:di,:],
            its,
            stps,
            n_states,
            # _initial_weights=weight,
            return_matrices=True
        )
        
        reweight_mean = np.mean(reweighted_distributions,axis=0)
        tmp = unmap(remp, reweight_mean, 231)
        
        tmp_max_ind = np.where(tmp.max()==tmp)
        tmp_max = aps.all_possible_states()[tmp_max_ind]
        out2.append(np.divide(tmp_max,max_truth))
        
        tmp_weighted_avg = weighted_avg(tmp)
        #out.append(np.sqrt((hist-reweight_mean)**2).sum())
        out.append(np.divide(tmp_weighted_avg,truth_weighted))
    out = np.array(out)   
    out2 = np.array(out2)


def upper_diagonal(cgtm):

    triU = np.triu(cgtm).T
    triL = np.tril(cgtm)
    ds = triU - triL
    ds_abs = np.abs(ds)
    ds_sum = ds_abs.sum()
    # print(ds_sum)
    # tmp = cgtm-ds
    # tmp[tmp==0] = np.nan
    ds[ds == 0] = np.nan
    # cl = plt.pcolor(ds_abs,norm = colors.LogNorm())
    # plt.colorbar(cl)
    # plt.show()
    # plt.pcolor(ds)
    return ds_abs,ds_sum
    # plt.colorbar()
    
    # plt.show()

# import matplotlib.colors as colors
# ds_list = []
# ds_sum_list = []
# reweight = test_block()
# # reweight_50 = test_block_50_per()
# cgtm_me = np.loadtxt("cgtm.dat")
# fig,ax = plt.subplots(5,2,sharex=True,sharey=True,figsize=(5*2,5*5))
# for cgi,cg in enumerate(reweight[0:]):
#     ds,ds_sum = upper_diagonal(cg)
#     ds_sum_list.append(ds_sum)
#     ds[ds==0] = np.nan
#     cg[cg==0] = np.nan
#     if cgi == 1:
#         ds_ax = ax[0,0].pcolor(cg,cmap='RdBu_r')
#         ds_ax1 = ax[0,1].pcolor(ds,cmap='RdBu_r',norm = colors.LogNorm())
#         plt.colorbar(ds_ax1,ax=ax[0,1])
#         plt.colorbar(ds_ax,ax=ax[0,0])
#         ax[0,1].set_title("Iteration 1")
#         ds_list.append(ds)


#     if cgi == 2:
#         ds_ax = ax[1,0].pcolor(cg,cmap='RdBu_r')
#         ds_ax = ax[1,1].pcolor(ds,cmap='RdBu_r',norm = colors.LogNorm())
#         plt.colorbar(ds_ax1,ax=ax[1,1]) 
#         plt.colorbar(ds_ax,ax=ax[1,0])
#         ax[1,1].set_title("Iteration 2")
#         ds_list.append(ds)

#     if cgi == 10:
#         ds_ax = ax[2,0].pcolor(cg,cmap='RdBu_r')
#         ds_ax = ax[2,1].pcolor(ds,cmap='RdBu_r',norm = colors.LogNorm())
#         plt.colorbar(ds_ax1,ax=ax[2,1]) 
#         plt.colorbar(ds_ax,ax=ax[2,0])
#         ax[2,1].set_title("Iteration 10")
#         ds_list.append(ds)
        
#     if cgi == 50:
#         ds_ax = ax[3,0].pcolor(cg,cmap='RdBu_r')
#         ds_ax = ax[3,1].pcolor(ds,cmap='RdBu_r',norm = colors.LogNorm())
#         plt.colorbar(ds_ax1,ax=ax[3,1]) 
#         plt.colorbar(ds_ax,ax=ax[3,0])
#         ax[3,1].set_title("Iteration 50")
#         ds_list.append(ds)

#     if cgi == 149:
#         ds_ax = ax[4,0].pcolor(cg,cmap='RdBu_r')
#         # ds_ax = ax[2,1].pcolor(cg,cmap='RdBu_r')
#         ds_ax = ax[4,1].pcolor(ds,cmap='RdBu_r',norm = colors.LogNorm())
#         plt.colorbar(ds_ax1,ax=ax[4,1])
#         plt.colorbar(ds_ax,ax=ax[4,0])
        # ax[4,1].set_title("Iteration 149")
        # ds_list.append(ds)


# plt.savefig("upper_diag_map.pdf")
# plt.close()
# # plt.show()
# # ds_sum_list[0]= 0#53.2062114947797
# plt.plot(np.linspace(1,len(ds_sum_list),len(ds_sum_list)-1),ds_sum_list[1:])
# plt.savefig("upper_diag.pdf")
# plt.close()
# plt.show()


def KL_Analysis():
    import matplotlib.colors as mc
    
    for i in [40,41,42,43,44,45,46,47,48,49,50]: 
        dat = pd.read_csv('KL_iteration_windows_lag_allframes_maps%i.csv'%i,index_col=0)
        its, wind = np.shape(dat)
        x,y = np.linspace(1,its,its),np.linspace(1,wind,wind)
        X,Y = np.meshgrid(x,y)
        plt.pcolor(X.T,Y.T,dat,edgecolor="face",cmap='RdBu')#,norm=mc.LogNorm(vmin=dat.min().min(), vmax=dat.max().max()))
        # plt.plot(dat)
        plt.xlabel("Iterations")
        plt.ylabel("Windows")
        plt.colorbar()
        plt.show()
        # plt.savefig("LK_Heat%i.pdf"%i)
        # plt.close()
        
        # fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
        # surf = ax.plot_surface(X, Y, dat.values.T)
        # plt.savefig("LK_3Space.pdf")
        # plt.close()
# KL_Analysis() 

# truth = np.array(pd.read_csv("~/CG/data/pi_raw_low_inactive_long_cg.csv",index_col=0).values).T
