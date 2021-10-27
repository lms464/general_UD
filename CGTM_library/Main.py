import matplotlib.pyplot as plt
import CGTM_Calculations as cgc
import CGTM_Plotting as cgp
import all_possible_states as aps



def long_state_diff():
    import choose_path as chp
    import pandas as pd
    import numpy as np
    
    i = pd.read_csv("%s/CG/data/pi_raw_inactive_longcg.csv"%chp.choose_path()[1],index_col=0).T
    a = pd.read_csv("%s/CG/data/pi_raw_active_longcg.csv"%chp.choose_path()[1],index_col=0).T    
    act = pd.read_csv("%s/CG/data/pi_raw_active_longcg_time.csv"%chp.choose_path()[1],index_col=0).T
    inact = pd.read_csv("%s/CG/data/pi_raw_inactive_longcg_time.csv"%chp.choose_path()[1],index_col=0).T
    in_data = tmp = pd.DataFrame([i.values[0],inact[0],inact[1],inact[2]]).T
    act_data = tmp = pd.DataFrame([a.values[0],act[0],act[1],act[2]]).T
    fig,ax = plt.subplots(2,4,sharey=True,sharex=True)
    ax[0,0].bar(i.columns,in_data[0])

    ax[0,1].bar(np.arange(0,len(inact)),in_data[0]-in_data[1])
    ax[0,2].bar(np.arange(0,len(inact)),in_data[0]-in_data[2])
    ax[0,3].bar(np.arange(0,len(inact)),in_data[0]-in_data[3])
    #ax[0,3].bar(np.arange(0,len(inact)),inact.mean(axis=1) - ((inact[1]-inact[2])+(inact[0]-inact[1])+(inact[0]-inact[2])))
    ax[1,0].bar(a.columns,a.values[0])
    ax[1,1].bar(np.arange(0,len(act)),act_data[0]-act_data[1])
    ax[1,2].bar(np.arange(0,len(act)),act_data[0]-act_data[2])
    ax[1,3].bar(np.arange(0,len(act)),act_data[0]-act_data[3])
    #ax[1,3].bar(np.arange(0,len(act)),act.mean(axis=1) - ((act[1]-act[2])+(act[0]-act[1])+(act[0]-act[2])))
    plt.tight_layout()
    plt.savefig("States_Long_Iterative_dif.pdf")
    plt.close()
    
# long_state_diff()

def ternary_CG():
    import choose_path as chp
    import pandas as pd
    uot = 0
    oot = 0
    fig,ax = plt.subplots(2,3,figsize=(12,8))
    cgp.Ternary_Heat_Map("CG/data/pi_eq_inactive_shortcg","",ax=ax[0,1],initial="CG/data/init_raw_inactive_shortcg")
    cgp.Ternary_Heat_Map("CG/data/pi_raw_inactive_longcg","",ax=ax[0,0])
    cgp.Ternary_Heat_Map("CG/data/pi_eq_inactive_longcg","","CG/data/pi_raw_inactive_longcg",ax=ax[0,2])
    cgp.Ternary_Heat_Map("CG/data/pi_eq_active_shortcg","",ax=ax[1,1] ,initial="CG/data/init_raw_active_shortcg")
    uot,sm1 = cgp.Ternary_Heat_Map("CG/data/pi_raw_active_shortcg","",ax=ax[1,0], out=uot)
    # cgp.Ternary_Scatter(ax[1,2], "CG/data/pi_eq_active_shortcg")
    oot,sm2 = cgp.Ternary_Heat_Map("CG/data/pi_eq_active_shortcg","","CG/data/pi_raw_active_longcg",ax=ax[1,2],out=oot)
    cax = plt.axes([1, 0.25, 0.025, 0.5])
    plt.colorbar(sm2, cax=cax,format='%.3f')
    cax = plt.axes([0.15,-.025,0.45,0.025])
    plt.colorbar(sm1,cax=cax,format="%.3f",orientation="horizontal")
    plt.tight_layout()
    # plt.show()
    plt.savefig("CG_ternary_short-long.pdf",bbox_inches='tight')
    plt.close()
    # fig,ax = plt.subplots(2,4,figsize=(12,8))
    # i = pd.read_csv("%s/CG/data/pi_raw_inactive_longcg.csv"%chp.choose_path()[1],index_col=0).T
    # a = pd.read_csv("%s/CG/data/pi_raw_active_longcg.csv"%chp.choose_path()[1],index_col=0).T    
    # act = pd.read_csv("%s/CG/data/pi_raw_active_longcg_time.csv"%chp.choose_path()[1],index_col=0).T
    # inact = pd.read_csv("%s/CG/data/pi_raw_inactive_longcg_time.csv"%chp.choose_path()[1],index_col=0).T
    # in_data = tmp = pd.DataFrame([i.values[0],inact[0],inact[1],inact[2]]).T
    # act_data = tmp = pd.DataFrame([a.values[0],act[0],act[1],act[2]]).T
    # uot,sm1 = cgp.Ternary_Heat_Map(in_data[0],fl_name="",ax=ax[0,0],out=uot)
    # oot,sm2 =cgp.Ternary_Heat_Map(in_data[0],"",in_data[1],ax=ax[0,1],)
    # cgp.Ternary_Heat_Map(in_data[0],"",in_data[2],ax=ax[0,2],)
    # cgp.Ternary_Heat_Map(in_data[0],"",in_data[3],ax=ax[0,3],)
    
    # cgp.Ternary_Heat_Map(act_data[0],fl_name="",ax=ax[1,0],out=uot)
    # cgp.Ternary_Heat_Map(act_data[0],"",act_data[1],ax=ax[1,1],)
    # cgp.Ternary_Heat_Map(act_data[0],"",act_data[2],ax=ax[1,2],)
    # cgp.Ternary_Heat_Map(act_data[0],"",act_data[3],ax=ax[1,3],)
    
    # cax = plt.axes([1, 0.25, 0.025, 0.5])
    # plt.colorbar(sm2, cax=cax,format='%.3f')
    # cax = plt.axes([0.15,-.025,0.45,0.025])
    # plt.colorbar(sm1,cax=cax,format="%.3f",orientation="horizontal")
    # plt.tight_layout()
    # plt.savefig("CG_ternary_long-long.pdf",bbox_inches='tight')
    # plt.close()

# ternary_CG() 

def ternary_iterative():
    import choose_path as chp
    import pandas as pd
    import numpy as np
    uot = 0
    oot = 0
    
    long_inact = pd.read_csv("%s/CG/data_dx2/pi_eq_inactive_shortcg.csv"%(chp.choose_path()[1]),index_col=0).T
    long_act = pd.read_csv("%s/CG/data_dx2/pi_eq_active_shortcg.csv"%(chp.choose_path()[1]),index_col=0).T
    
    
    state_act = pd.read_csv("%s/CG/data_dx2/act_short_binned_sim_pi.csv"%(chp.choose_path()[1]),index_col=0).T
    state_inact = pd.read_csv("%s/CG/data_dx2/inact_short_binned_sim_pi.csv"%(chp.choose_path()[1]),index_col=0).T
    
    # dpi_ref = (long_inact - state_inact[0].values).iloc[0,:]
    # import numpy as np
    for si, sj in zip(state_act[-2:-1], state_inact[-2:-1]):
        # d_sig_i = np.sqrt((long_inact.values[0]-state_inact.iloc[:,si])**2).mean()# - np.var()
        # d_sig_a = np.sqrt((long_act.values[0]-state_act.iloc[:,si])**2).mean()
        # # d_sig_a = np.sqrt((long_act.values[0]-state_act.T[si])**2).sum()
        # plt.plot(1+si*.4,d_sig_i,"bo")
        # plt.plot(1+si*.4,d_sig_a,"go")

        # print(si,d_sig_i,d_sig_a)
        fig,ax = plt.subplots(2,3,figsize=(8,8))
        left, width = .0, .95
        bottom, height = .0, .95
        right = left + width
        top = bottom + height
        ax[0,0].text(left, top, '%f'%(si*.4),
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax[0,0].transAxes)
        uot,sm1 = cgp.Ternary_Heat_Map(state_inact[si],fl_name="",ax=ax[1,1],out=uot)
        cgp.Ternary_Heat_Map("CG/data/pi_raw_inactive_shortcg","",ax=ax[1,0])
        oot,sm2 = cgp.Ternary_Heat_Map(long_inact, leaflet_in2=state_inact[si].values,fl_name="",ax=ax[1,2],out=oot)

        # cgp.Ternary_Heat_Map(state_act[sj],fl_name="",ax=ax[0,1])
        # cgp.Ternary_Heat_Map("CG/data/pi_raw_active_shortcg","",ax=ax[0,0])
        # cgp.Ternary_Heat_Map(long_act, leaflet_in2=state_act[sj].values,fl_name="",ax=ax[0,2],out=oot)
        # cax = plt.axes([1, 0.25, 0.025, 0.5])
        # plt.colorbar(sm2, cax=cax,format='%.3f')
        # cax = plt.axes([0.15,-.025,0.45,0.025])
        # plt.colorbar(sm1,cax=cax,format="%.3f",orientation="horizontal")
        # plt.tight_layout()
        # # plt.show()
        # plt.savefig("ternary_short-short_gif%i.png"%si)
        # plt.close()
    plt.savefig("mean_diff_short-short.pdf")
    plt.close()

# ternary_iterative()

def states_CG():
    fig, ax = plt.subplots(2,3,figsize=(8,6),sharey=True,sharex=True)
    cgp.plot_state_dist("~/UDel/CG/data/pi_eq_inactive_shortcg",ax[0,0])
    cgp.plot_state_dist("~/UDel/CG/data/pi_raw_inactive_shortcg",ax[0,1])
    cgp.diff_plot("~/UDel/CG/data/pi_eq_inactive_shortcg","~/UDel/CG/data/pi_raw_inactive_shortcg", ax[0,2])
    cgp.plot_state_dist("~/UDel/CG/data/pi_eq_active_shortcg",ax[1,0])
    cgp.plot_state_dist("~/UDel/CG/data/pi_raw_active_shortcg",ax[1,1])
    cgp.diff_plot("~/UDel/CG/data/pi_eq_active_shortcg","~/UDel/CG/data/pi_raw_active_shortcg", ax[1,2])
    plt.tight_layout()
    plt.savefig("CG_state_dist_short-short.pdf")
    plt.close()
states_CG()


def tmp():
    import pandas as pd
    import matplotlib.pyplot as plt
    import numpy as np
    alps = aps.all_possible_states()
    act_pi = pd.read_csv("~/UDel/CG/data_dx2/pi_eq_active_shortcg.csv",index_col=0)
    inact_pi = pd.read_csv("~/UDel/CG/data_dx2/pi_eq_inactive_shortcg.csv",index_col=0)
    inact_raw = pd.read_csv("~/UDel/CG/data_dx2/pi_raw_inactive_shortcg.csv",index_col=0)
    act_raw = pd.read_csv("~/UDel/CG/data_dx2/pi_raw_active_shortcg.csv",index_col=0)
    plt.bar(np.arange(0,len(act_pi)),act_pi.T.values[0])#-act_raw.T.values[0])
    plt.bar(np.arange(0,len(act_pi)),act_pi.T.values[0]-act_raw.T.values[0])
    plt.savefig("tmp.pdf")
    plt.close()
# tmp()
def state_iterative():
    
    import choose_path as chp
    import pandas as pd
    import numpy as np

    state_act = pd.read_csv("%s/CG/data_dx2/act_short_binned_time_pi.csv"%(chp.choose_path()[1]),index_col=0).T
    state_inact = pd.read_csv("%s/CG/data_dx2/inact_short_binned_time_pi.csv"%(chp.choose_path()[1]),index_col=0).T
    
    
    long_act = pd.read_csv("%s/CG/data_dx2/pi_raw_active_shortcg_time.csv"%(chp.choose_path()[1]),index_col=0).T
    long_inact = pd.read_csv("%s/CG/data_dx2/pi_raw_inactive_shortcg_time.csv"%(chp.choose_path()[1]),index_col=0).T
    
    dsi,dss = 0,0
    
    for si in state_inact:

        fig, ax = plt.subplots(2,3,figsize=(8,6),sharey=True,sharex=True)
        left = .0
        bottom, height = .0, .95
        top = bottom + height
        ax[0,0].text(left, top, '%f'%((si+1)*.4),
            horizontalalignment='left',
            verticalalignment='top',
            transform=ax[0,0].transAxes
        )

        txt = '(0.55:0.15:0.30 DPPC:DOPC:Chol) Martini FF. 50 simulations for 50 ns\nStepsize between Frames: .4 ns. The stepsize between states: 2%'
        # fig.text(.5, -.05, txt, ha='center')
        dsi = np.mean(long_inact.iloc[:,:2+si],axis=1).values - state_inact[si].values
        
        a=ax[0,0].bar(state_inact.index,np.mean(long_inact.iloc[:,:2+si],axis=1))
        b=ax[0,1].bar(state_inact.index, state_inact[si].values)
        c=ax[0,2].bar(state_inact.index, dsi)
        ax[0,0].set_title("Raw From Sim")
        ax[0,0].set_ylabel("Inactive")
        ax[1,0].set_ylabel("Active")
        ax[0,1].set_title("Pi_eq CGTM")
        ax[0,2].set_title("difference")
        dss = np.mean(long_act.iloc[:,:2+si],axis=1).values - state_act[si].values

        a=ax[1,0].bar(state_act.index,np.mean(long_act.iloc[:,:2+si],axis=1))
        b=ax[1,1].bar(state_act.index, state_act[si].values)
        c=ax[1,2].bar(state_act.index, dss)
        del a
        del b
        del c
        dss,dsi = 0,0
        plt.ylim(-.075,.075)
        plt.tight_layout()
        plt.savefig("dx2_state-short-short_frame%s_2.pdf"%si)
        plt.close()
    # return dsi,dss
            
# o = state_iterative()
# del o

def over_lapped_state_network() :
    
    import pandas as pd
    import choose_path as chp
    fig, ax = plt.subplots(1,2,figsize=(36,24))
  

    # long_act = pd.read_csv("%s/CG/data/pi_raw_inactive_longcg.csv"%(chp.choose_path()[1]),index_col=0).T
    # long_inact = pd.read_csv("%s/CG/data/pi_raw_active_longcg.csv"%(chp.choose_path()[1]),index_col=0).T
    
    # state_inact = pd.read_csv("%s/CG/data/pi_eq_inactive_longcg.csv"%(chp.choose_path()[1]),index_col=0).T
    # state_act = pd.read_csv("%s/CG/data/pi_eq_active_longcg.csv"%(chp.choose_path()[1]),index_col=0).T
    
    # d_active = long_act -(long_act - state_act)
    # d_active[d_active <= 0 ] = 0
    
    # d_inactive = long_inact - (long_inact - state_inact)
    # d_inactive[d_active <= 0 ] = 0
    # di = d_inactive.T[d_inactive.T["0"]>0].index
    
    fig, ax = plt.subplots(1,2,figsize=(24,12))

    cgp.network_plot(ax[0],act="inactive")
    cgp.network_plot(ax[1],act="active")

    # cgp.network_plot(ax[0],act="inactive")
    # cgp.network_plot(ax[0],leaflet[0], kind, act)
    # cgp.network_plot(ax[1],leaflet[1], kind, act)


    # cgp.network_plot("inactive", ax[0])
    # cgp.network_plot("active", ax[1])
    plt.tight_layout()
    # plt.show()
    plt.savefig("State_overlap_cutoff0.pdf")
    plt.close()    
    
# over_lapped_state_network()
    
def CGTM():
    fig,ax = plt.subplots(1,2,figsize=(8,4),sharey=True,sharex=True)
    for ii,i in enumerate(["inactive","active"]):
        act = cgc.CGTM_Calculations("",1,"cg",i,"short")
        act_CGTM = act.build_CGTM()[-1]
        cgp.plot_cgtm(act_CGTM,ax[ii])
    plt.tight_layout()
    plt.savefig("CG_CGTM.pdf")
    plt.close()

    
def SI_CG():
    fig, ax = plt.subplots(2,2,figsize=(8,6),sharey='row',sharex="row")
    cgp.plot_state_traj("~/CG/data/states/short_inactive",ax[0,0])
    s1 = cgp.plot_state_traj("~/CG/data/states/short_active",ax[0,1])
    cgp.plot_state_traj("~/CG/data/states/long_inactive",ax[1,0])
    cgp.plot_state_traj("~/CG/data/states/long_active",ax[1,1])
    cax = plt.axes([.96, 0.25, 0.025, 0.5])
    plt.colorbar(s1,cax=cax)
    plt.savefig("CG_SI.pdf",bbox_inches='tight')
    plt.close()

# ternary_CG()
# ternary_iterative()
# states_CG()   
# CGTM()    
# SI_CG()

def Confidence():
    
    # cgc.CGTM_Calculations("",1,"cg","inactive","short").calc_confidence(1)
    # cgc.CGTM_Calculations("",1,"cg","active","short").calc_confidence(1)
    
    
    '''
    Converges after 9 steps, which is the same number of simulations
    suggesting something is wrong with data collection (ie I'm looping
    over the wrong axis)
    '''
    
    path="~/UDel/CG/data"
    fig, ax = plt.subplots(1,2,figsize=(12,6),sharex="row")
    for jj, j in enumerate(["inactive","active"]):
        cgp.plot_convergence(path+"/"+j+"_ratio_pi_sum.csv",ax[jj])
    plt.savefig("confidence.pdf")
    plt.close()

def network(leaflet=None, kind=None, act=None):
    # leaflet and act should be lists with 2 
    # indexes ex ["SU","SL"]
    fig, ax = plt.subplots(1,2,figsize=(24,12))

    if leaflet == None:
        cgp.network_plot(ax[0],leaflet, kind, act[0])
        cgp.network_plot(ax[1],leaflet, kind, act[1])
    elif act == None:
        cgp.network_plot(ax[0],leaflet[0], kind, act)
        cgp.network_plot(ax[1],leaflet[1], kind, act)
    else:
        print("You should not have both leaflet and act")
        print("as 'None' or '<input>'.\n Exiting")
        return 0

    # cgp.network_plot("inactive", ax[0])
    # cgp.network_plot("active", ax[1])
    plt.tight_layout()
    plt.savefig("State_transitions2.pdf")
    plt.close()

def sig_conv(SL,SU,kind):
    cgp.plot_sigConverge(SL,SU,kind)
  
import numpy as np
# test1 = cgc.CGTM_Calculations("",1,"cg","inactive","short")#.sigConverge_time()
# # test2 = cgc.CGTM_Calculations("",1,"cg","inactive","short").sigConverge_simulations()
# # sig_conv(test1,test2,'cg')

test1 = cgc.CGTM_Calculations("SU", 1, "chg",act=None)#.sigConverge_time()
a = test1.build_CGTM()#write_pi_eq()
b = test1.build_raw()#write_pi_raw()

pi_raw = b[0]
pi_eq = a[0]
cgtm = a[-1]
alps = aps.all_possible_states()
empt = []
for i in range(0,len(alps)):
    for j in range(0,len(alps)):
        if cgtm[i,j] == 0 or cgtm[j,i] == 0:
            continue
        if np.isclose(pi_eq[i] *cgtm[i,j],pi_eq[i] * cgtm[j,i])==True:
            empt.append([(i,j),pi_eq[i]*cgtm[i,j],pi_eq[i]*cgtm[j,i]])
        else:
            empt.append([(i,j),0,0])



# pi_rand = pi_raw.copy()
# pi_rand[pi_rand > 0] = 1
pi_rand = np.random.rand(len(pi_eq))*np.random.rand(len(pi_eq))
pi_rand = pi_rand / pi_rand.sum()
pi_rand_tmp = pi_rand.copy()
pi_rand_init = pi_rand.copy()

i = 1
while np.allclose(pi_rand, pi_eq) == False:
    pi_rand_tmp = pi_rand_tmp @ cgtm  
    pi_rand = pi_rand_tmp/pi_rand_tmp.sum()
    i = i + 1
    plt.plot(i,pi_eq[219],'ko')
    plt.plot(i,pi_rand_tmp[219],'bo')
    plt.plot(i,pi_rand[219],'ro')
    print(i)
# plt.ylim([0,1.5])
plt.title("Markov 'simulation'")
plt.legend([r"max($\pi^{eq}$)","Same state Sim (Normalized)","Un-Normalized"])
plt.xlabel(r"Steps to converge to $\pi^{eq}$")
plt.ylabel("Sum of states (should be 1%)")
plt.savefig("Markov_Sim.pdf")
plt.close()

    

# test2 = cgc.CGTM_Calculations("SL", 1, "chg",act=None).sigConverge_time()
# sig_conv(test1,test2,'chg')

# network(["SU","SL"],"sat",None)
# import numpy as np
# test1 = cgc.CGTM_Calculations("",1,"cg","active","short").weighted_avg()
# fig, ax = plt.subplots(1,1)
# oot = 0
# oot,sm2 = cgp.Ternary_Heat_Map("CG/data/pi_eq_active_shortcg","","CG/data/pi_raw_active_shortcg",ax=ax,out=oot)
# plt.colorbar(sm2)