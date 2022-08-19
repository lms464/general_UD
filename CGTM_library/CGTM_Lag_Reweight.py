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

def lag_cumDiff(lagmin,lagmax):
    import matplotlib.pyplot as plt
    data = pd.read_csv("/home/sharplm/CG/data/states/low_short_inactive_tmpNEW.csv",index_col=0)#low_long_inactive_tmpNEW.csv",index_col=0)#"~/CG/data/states/low_short_inactive_tmp.csv",index_col=0)
    dat_ind = ["%i"%i for i in range(int(len(data.T)/2),len(data.T))]
    dat_ind = ["0"] + dat_ind
    trj,n_states,remp = ta.remap_trajs(data[dat_ind].values)
    hist = np.histogram(trj,bins=n_states,range=(0,n_states))[0]
    hist = hist/hist.sum()
    true_hist = pd.read_csv("/home/sharplm/CG/data/pi_eq_low_inactive_long_cgtmp_231.csv",index_col=0).values
    true_hist = true_hist.T[0]

    its = 1
    stps = 1
    out2 = []
    cum_diff_Long = []
    cum_diff_Short = []
    for tau in np.arange(lagmin,lagmax,1):
        reweighted_distributions, reweighted_matrices, eign = ta.optimized_resliced_voelz(
            trj[:,::tau],
            its,
            stps,
            n_states,
            return_eig=True,
            return_matrices=True
        )
        out2.append(eign)
        reweight_mean = np.mean(reweighted_distributions,axis=0)
        reweight_map = unmap(remp, reweight_mean, 231)
        hist_map = unmap(remp,hist,231)

        cum_diff_Short.append(cummulative_difference(hist_map,reweight_map))
        cum_diff_Long.append(cummulative_difference(true_hist,reweight_map))
    np.savetxt("CGTM_Lag.dat",np.array([np.arange(lagmin,lagmax,1),cum_diff_Short,cum_diff_Long]).T)
    # plt.plot(np.arange(lagmin,lagmax,1),cum_diff_Short);plt.plot(np.arange(lagmin,lagmax,1),cum_diff_Long)
    #return cum_diff_Short,cum_diff_Long


def window_cumDiff(stpsmin,stpsmax=None):
    import matplotlib.pyplot as plt
    data = pd.read_csv("/home/sharplm/CG/data/states/low_short_inactive_tmpNEW.csv",index_col=0)#low_long_inactive_tmpNEW.csv",index_col=0)#"~/CG/data/states/low_short_inactive_tmp.csv",index_col=0)
    dat_ind = ["%i"%i for i in range(int(len(data.T)/2),len(data.T))]
    dat_ind = ["0"] + dat_ind
    trj,n_states,remp = ta.remap_trajs(data[dat_ind].values)
    hist = np.histogram(trj,bins=n_states,range=(0,n_states))[0]
    hist = hist/hist.sum()
    true_hist = pd.read_csv("/home/sharplm/CG/data/pi_eq_low_inactive_long_cgtmp_231.csv",index_col=0).values
    true_hist = true_hist.T[0]
    tau = 1
    its = 300
    # stps = 1
    out2 = []
    cum_diff_Long = []
    cum_diff_Short = []
    if stpsmax == None or stpsmax > len(trj.T[:,::tau])*2/3:
        stpsmax = int(len(trj.T[:,::int(tau)])*2/3)
    for stp in np.arange(stpsmin,stpsmax,1):
        try:
            reweighted_distributions, reweighted_matrices, eign = ta.optimized_resliced_voelz(
                trj[:,::int(tau)],
                its,
                stp,
                n_states,
                return_eig=True,
                return_matrices=True
            )
        except:
            break
        out2.append(eign)
        reweight_mean = np.mean(reweighted_distributions,axis=0)
        reweight_map = unmap(remp, reweight_mean, 231)
        hist_map = unmap(remp,hist,231)

        cum_diff_Short.append(cummulative_difference(hist_map,reweight_map))
        cum_diff_Long.append(cummulative_difference(true_hist,reweight_map))

    # plt.plot(np.arange(stpsmin,stpsmax,1),cum_diff_Short);
    # plt.plot(np.arange(stpsmin,stpsmax,1),cum_diff_Long)
    np.savetxt("CGTM_Lag.dat",np.array([np.arange(stpsmin,stpsmax,1),cum_diff_Short,cum_diff_Long]).T)
    # return cum_diff_Short,cum_diff_Long

def compair_Lag_window_cumDiff(lagmin, lagmax, stpsmin,stpsmax):
    import matplotlib.pyplot as plt
    data = pd.read_csv("/home/sharplm/CG/data/states/low_short_inactive_tmpNEW.csv",index_col=0)#low_long_inactive_tmpNEW.csv",index_col=0)#"~/CG/data/states/low_short_inactive_tmp.csv",index_col=0)
    dat_ind = ["%i"%i for i in range(int(len(data.T)/2),len(data.T))]
    dat_ind = ["0"] + dat_ind
    trj,n_states,remp = ta.remap_trajs(data[dat_ind].values)
    hist = np.histogram(trj,bins=n_states,range=(0,n_states))[0]
    hist = hist/hist.sum()
    true_hist = pd.read_csv("/home/sharplm/CG/data/pi_eq_low_inactive_long_cgtmp_231.csv",index_col=0).values
    true_hist = true_hist.T[0]
    # tau = 1
    its = 300
    # stps = 1
    # out2 = []
    # cum_diff_Long = []
    # cum_diff_Short = []
    # if stpsmax == None or stpsmax > len(trj.T[:,::tau])*2/3:
    #     stpsmax = int(len(trj.T[:,::int(tau)])*2/3)
    cum_diff_Short = pd.DataFrame(index=np.arange(stpsmin,stpsmax,1),columns=np.arange(lagmin, lagmax, 1))
    cum_diff_Long = cum_diff_Short.copy()
    for tau in np.arange(lagmin, lagmax, 1):
        for stp in np.arange(stpsmin,stpsmax,1):
            print("tau:%i\twind:%i"%(int(tau),int(stp)))
            try:
                reweighted_distributions, reweighted_matrices, eign = ta.optimized_resliced_voelz(
                    trj[:,::int(tau)],
                    its,
                    stp,
                    n_states,
                    return_eig=True,
                    return_matrices=True
                )
            # out2.append(eign)
                reweight_mean = np.mean(reweighted_distributions,axis=0)
                reweight_map = unmap(remp, reweight_mean, 231)
                hist_map = unmap(remp,hist,231)

                cum_diff_Short[tau][stp] = cummulative_difference(hist_map,reweight_map)
                cum_diff_Long[tau][stp] = cummulative_difference(true_hist,reweight_map)
            except:
                cum_diff_Short[tau][stp] = -1
                cum_diff_Long[tau][stp] = -1

    # plt.plot(np.arange(stpsmin,stpsmax,1),cum_diff_Short);
    # plt.plot(np.arange(stpsmin,stpsmax,1),cum_diff_Long)
    cum_diff_Long.to_csv("Long_diffs.csv")
    cum_diff_Short.to_csv("Short_diffs.csv")
    # return cum_diff_Short,cum_diff_Long

def Get_Lag_data(tau):
    data = pd.read_csv("/home/liam/lms464/CG/data/states/low_short_inactive_tmpNEW.csv",index_col=0)#low_long_inactive_tmpNEW.csv",index_col=0)#"~/CG/data/states/low_short_inactive_tmp.csv",index_col=0)
    dat_ind = ["%i"%i for i in range(int(len(data.T)/2),len(data.T))]
    dat_ind = ["0"] + dat_ind
    trj,n_states,remp = ta.remap_trajs(data[dat_ind].values)
    hist = np.histogram(trj,bins=n_states,range=(0,n_states))[0]
    hist = hist/hist.sum()
    true_hist = pd.read_csv("/home/liam/lms464/CG/data/pi_eq_low_inactive_long_cgtmp_231.csv",index_col=0).values
    true_hist = true_hist.T[0]
    
    its = 150
    stps = 1
    reweighted_distributions, reweighted_matrices, eign = ta.optimized_resliced_voelz(
        trj[:,::1],
        its,
        tau,
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

    reweigth_out = pd.DataFrame([alps[:,0],alps[:,1],alps[:,2],reweight_map])
    hist_out = pd.DataFrame([alps[:,0],alps[:,1],alps[:,2],hist_map])
    hist_out.T.to_csv("~/lms464/CG/data/pi_raw_low_inactive_short_cgtmp_iter.csv")
    reweigth_out.T.to_csv("~/lms464/CG/data/pi_eq_low_inactive_short_cgtmp_iter.csv")

# Get_Lag_data(205)
# short,long = 
# lag_cumDiff(1,50)
# window_cumDiff(1,50)
compair_Lag_window_cumDiff(1,50,1,50)