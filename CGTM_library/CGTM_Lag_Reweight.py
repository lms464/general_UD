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

def lag_cumDiff(lagmin,lagmax):
    import matplotlib.pyplot as plt
    data = pd.read_csv("/home/liam/lms464/CG/data/states/low_short_inactive_tmpNEW.csv",index_col=0)#low_long_inactive_tmpNEW.csv",index_col=0)#"~/CG/data/states/low_short_inactive_tmp.csv",index_col=0)
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
        print(r"$\tau$ = %s"%tau)
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
    np.savetxt("eig_OG2.dat",out2)
    np.savetxt("CGTM_Lag_OG2.dat",np.array([np.arange(lagmin,lagmax,1),cum_diff_Short,cum_diff_Long]).T)
    net_diff = np.abs(np.array(cum_diff_Long) - np.array(cum_diff_Short))
    net_diff_avg = running_mean(net_diff,10)
    net_diff[net_diff < .05] = np.nan
    plt.scatter(np.arange(lagmin,lagmax,1),cum_diff_Short)
    plt.scatter(np.arange(lagmin,lagmax,1),cum_diff_Long)    
    plt.scatter(np.arange(lagmin,lagmax,1),net_diff)
    plt.plot(np.arange(lagmin,lagmax-9,1),net_diff_avg,c='orange')
    plt.show()
    #return cum_diff_Short,cum_diff_Long

def Lag_analysis():
    # eigT = np.loadtxt("eig_tit.dat")
    eigO = np.loadtxt("eig_OG2.dat")
    dt = 1 #ns
    lag_time = np.arange(1,len(eigO)+1)
    # tau_T = -1 * lag_time/np.log(eigT)
    tau_O = -1 * lag_time/np.log(eigO)
    
    # plt.plot(lag_time,np.log(tau_T),'*')
    plt.plot(lag_time,np.log(tau_O),"o")
    plt.show()
    # plt.savefig("Lags_relax.pdf")
    # plt.close()
    


def window_cumDiff(stpsmin,stpsmax=None):
    import matplotlib.pyplot as plt
    data = pd.read_csv("/home/liam/lms464/CG/data/states/low_short_inactive_tmpNEW.csv",index_col=0)#low_long_inactive_tmpNEW.csv",index_col=0)#"~/CG/data/states/low_short_inactive_tmp.csv",index_col=0)
    dat_ind = ["%i"%i for i in range(int(len(data.T)/2),len(data.T))]
    dat_ind = ["0"] + dat_ind
    trj,n_states,remp = ta.remap_trajs(data[dat_ind].values)
    hist = np.histogram(trj,bins=n_states,range=(0,n_states))[0]
    hist = hist/hist.sum()
    true_hist = pd.read_csv("/home/liam/lms464//CG/data/pi_eq_low_inactive_long_cgtmp_231.csv",index_col=0).values
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
    net_diff = np.abs(np.array(cum_diff_Long) - np.array(cum_diff_Short))
    net_diff_avg = running_mean(net_diff,10)
    net_diff[net_diff < .05] = np.nan
    plt.scatter(np.arange(stpsmin,stpsmax,1),net_diff)
    plt.plot(np.arange(stpsmin,stpsmax-9,1),net_diff_avg,c='orange')
    plt.show()
    np.savetxt("CGTM_ReW_OG.dat",np.array([np.arange(stpsmin,stpsmax,1),cum_diff_Short,cum_diff_Long]).T)
    # return cum_diff_Short,cum_diff_Long

def compair_Lag_window_cumDiff(lagmin, lagmax, stpsmin,stpsmax):
    import matplotlib.pyplot as plt
    data = pd.read_csv("/home/liam/lms464/CG/data/states/low_short_inactive_tmpNEW.csv",index_col=0)#low_long_inactive_tmpNEW.csv",index_col=0)#"~/CG/data/states/low_short_inactive_tmp.csv",index_col=0)
    dat_ind = ["%i"%i for i in range(int(len(data.T)/2),len(data.T))]
    dat_ind = ["0"] + dat_ind
    trj,n_states,remp = ta.remap_trajs(data[dat_ind].values)
    hist = np.histogram(trj,bins=n_states,range=(0,n_states))[0]
    hist = hist/hist.sum()
    true_hist = pd.read_csv("/home/liam/lms464/CG/data/pi_eq_low_inactive_long_cgtmp_231.csv",index_col=0).values
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
    cum_diff_Long.to_csv("Long_diffs_OG.csv")
    cum_diff_Short.to_csv("Short_diffs_OG.csv")
    # return cum_diff_Short,cum_diff_Long


def plot_lag_iter():
    s = pd.read_csv("Short_diffs.csv",index_col=0,header=0)
    i_min, i_max = 0,60
    fig,ax = plt.subplots(1,1,sharey=True)
    ax.pcolor(s.columns[i_min:i_max],s.index,s[s.columns[i_min:i_max]],vmin=0,vmax=1)
    for i in s.columns[i_min:i_max]:
        ind = np.where(s[i][(s[i]>=0) & (s[i]<2)].max()==s[i])[0][0]+1
        # ax[1].scatter(s[i][ind],ind)
        ax.scatter(i,ind)
    ax.set_xticks(np.arange(-1,i_max-i_min,5))
    ax.set_yticks(np.arange(0,20,5))
    plt.ylim(0,50)

    plt.show()
    
# plot_lag_iter()

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
    plt.bar(np.arange(0,231),hist_map-reweight_map)
    plt.bar(np.arange(0,231),true_hist-reweight_map,alpha=.5)
    plt.show()

    
    reweigth_out = pd.DataFrame([alps[:,0],alps[:,1],alps[:,2],reweight_map])
    hist_out = pd.DataFrame([alps[:,0],alps[:,1],alps[:,2],hist_map])
    hist_out.T.to_csv("~/lms464/CG/data/pi_raw_low_inactive_short_cgtmp_iter.csv")
    reweigth_out.T.to_csv("~/lms464/CG/data/pi_eq_low_inactive_short_cgtmp_iter.csv")
    
    
    '''
    ind_dict_b = dict((k,i) for i,k in enumerate(b))
    ind_dict_a = dict((k,i) for i,k in enumerate(a))
    index = set(ind_dict_a.values()).intersection(ind_dict_b.values())
    inter = pd.DataFrame([alps[:,0],alps[:,1],alps[:,2],a,b])
    intersections are inter[index] --> use this to figure out the "concentration director"
    '''
    
    return hist_map, true_hist

def Get_Wind_data(wnd,tau):
    data = pd.read_csv("/home/liam/lms464/CG/data/states/low_short_inactive_tmpNEW.csv",index_col=0)#low_long_inactive_tmpNEW.csv",index_col=0)#"~/CG/data/states/low_short_inactive_tmp.csv",index_col=0)
    dat_ind = ["%i"%i for i in range(int(len(data.T)/2),len(data.T))]
    dat_ind = ["0"] + dat_ind
    trj,n_states,remp = ta.remap_trajs(data[dat_ind].values)
    hist = np.histogram(trj,bins=n_states,range=(0,n_states))[0]
    hist = hist/hist.sum()
    true_hist = pd.read_csv("/home/liam/lms464/CG/data/pi_eq_low_inactive_long_cgtmp_231.csv",index_col=0).values
    true_hist = true_hist.T[0]
    
    its = 150
    # tau = 6
    reweighted_distributions, reweighted_matrices, eign = ta.optimized_resliced_voelz(
        trj[:,::tau],
        its,
        wnd,
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
    plt.bar(np.arange(0,231),hist_map-reweight_map)
    plt.bar(np.arange(0,231),true_hist-reweight_map,alpha=.5)
    plt.show()

    
    # reweigth_out = pd.DataFrame([alps[:,0],alps[:,1],alps[:,2],reweight_map])
    # hist_out = pd.DataFrame([alps[:,0],alps[:,1],alps[:,2],hist_map])
    # hist_out.T.to_csv("~/lms464/CG/data/pi_raw_low_inactive_short_cgtmp_wind.csv")
    # reweigth_out.T.to_csv("~/lms464/CG/data/pi_eq_low_inactive_short_cgtmp_wind.csv")


# for tau,stps in zip([ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  1,  1,1,  2],[20, 24, 35, 51, 58, 65, 69, 71, 79, 84, 85, 88, 91, 15, 72, 89,96, 68]):
a,b = Get_Lag_data(37,1)
# for i in np.linspace(1,75):
#     print("tau: %i"%int(i))
#     Get_Wind_data(20,int(i))
# Lag_analysis()
# window_cumDiff(1,251)
# lag_cumDiff(1,75)
# Get_Wind_data(60,1)
# # short,long = 
# Get_Lag_data(60,1)
# compair_Lag_window_cumDiff(1,75,1,100)