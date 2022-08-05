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

def cummulative_difference(pi_ref,pi_eq):
	diff_list = []
	for ref,eq in zip(pi_ref,pi_eq):
		diff_list.append(ref - eq)
	return np.sum(diff_list)

def lag_cumDiff():
    data = pd.read_csv("/home/liam/lms464/CG/data/states/low_short_inactive_tmpNEW.csv",index_col=0)#low_long_inactive_tmpNEW.csv",index_col=0)#"~/CG/data/states/low_short_inactive_tmp.csv",index_col=0)
    dat_ind = ["%i"%i for i in range(int(len(data.T)/2),len(data.T))]
    dat_ind = ["0"] + dat_ind
    trj,n_states,remp = ta.remap_trajs(data[dat_ind].values)
    hist = np.histogram(trj,bins=n_states,range=(0,n_states))[0]
    hist = hist/hist.sum()
    true_hist = pd.read_csv("/home/liam/lms464/CG/data/pi_eq_low_inactive_long_cgtmp_231.csv",index_col=0).values
    true_hist = true_hist.T[0]

    
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
    cum_diff_Long = []
    cum_diff_Short = []
    for tau in np.arange(1,201,1):
	    reweighted_distributions, reweighted_matrices, eign = ta.optimized_resliced_voelz(
	        trj[:,::tau],
	        its,
	        stps,
	        n_states,
	        return_eig=True,
	        # _initial_weights=weight,
	        return_matrices=True
	    )
	    out2.append(eign)
	    reweight_mean = np.mean(reweighted_distributions,axis=0)
	    # kl_list.append(KL(reweight_mean,hist))
	    alps = aps.all_possible_states()
	    reweight_map = unmap(remp, reweight_mean, 231)
	    hist_map = unmap(remp,hist,231)

	    cum_diff_Short.append(cummulative_difference(hist_map,reweight_map))
	    cum_diff_Short.append(cummulative_difference(true_hist,reweight_map))

    # reweigth_out = pd.DataFrame([alps[:,0],alps[:,1],alps[:,2],reweight_map])
    # hist_out = pd.DataFrame([alps[:,0],alps[:,1],alps[:,2],hist_map])
    # kl_list.append(KL(reweight_map,hist_map))
    # hist_out.T.to_csv("~/lms464/CG/data/pi_raw_low_inactive_short_cgtmpNEW.csv")
    # reweigth_out.T.to_csv("~/lms464/CG/data/pi_eq_low_inactive_short_cgtmpNEW.csv")

    return cum_diff_Short,cum_diff_Long
short,long = lag_cumDiff()