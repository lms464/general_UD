import matplotlib.pyplot as plt
import CGTM_Calculations as cgc
import CGTM_Plotting as cgp
import all_possible_states as aps



def long_state_diff():
    import choose_path as chp
    import pandas as pd
    import numpy as np
    from scipy.stats import entropy as entropy
    
    
    def KL(a, b):
        a = np.asarray(a, dtype=np.float)
        b = np.asarray(b, dtype=np.float)

        return np.sum(np.where(a != 0, a * np.log(a / b), 0))
    
    
    inact = pd.read_csv("%s/CG/data/states/low_long_inactive.csv"%chp.choose_path()[1],index_col=0).T
    inact_long = pd.read_csv("%s/CG/data/pi_raw_low_inactive_long_cg4.csv"%chp.choose_path()[1],index_col=0).values[:,0]
    # a = pd.read_csv("%s/CG/data/pi_raw_active_longcg.csv"%chp.choose_path()[1],index_col=0).T    
    # act = pd.read_csv("%s/CG/data/pi_raw_active_longcg_time.csv"%chp.choose_path()[1],index_col=0).T
    # inact = pd.read_csv("%s/CG/data/pi_raw_inactive_longcg_time.csv"%chp.choose_path()[1],index_col=0).T
    # in_data = tmp = pd.DataFrame([i.values[0],inact[0],inact[1],inact[2]]).T
    # act_data = tmp = pd.DataFrame([a.values[0],act[0],act[1],act[2]]).T
    fig,ax = plt.subplots(3,4,sharey=True,sharex=True)
    hists = [np.histogram(inact[i],bins=231,range=(0,231))[0] for i in inact]
    hists = [hist/np.sum(hist) for hist in hists]
    #print(np.abs(hists[0]-hists[1]).sum(),np.abs(hists[0]-hists[2]).sum(),np.abs(hists[1]-hists[2]).sum())
    #print(np.abs(inact_long-hists[0]).sum(),np.abs(inact_long-hists[1]).sum(),np.abs(inact_long-hists[2]).sum())
    
    print(entropy(hists[0],hists[1]))#KL(hists[0],hists[1]),KL(hists[0],hists[2]),KL(hists[1],hists[2]))
    
    ax[0,0].bar(np.linspace(0,231,len(hists[0])),hists[0])
    ax[0,1].bar(np.linspace(0,231,len(hists[0])),hists[1])
    ax[0,2].bar(np.linspace(0,231,len(hists[0])),hists[2])
    ax[0,3].bar(np.linspace(0,231,len(hists[0])),inact_long)

    ax[1,0].bar(np.linspace(0,231,len(hists[0])),hists[0]-hists[1])
    ax[1,1].bar(np.linspace(0,231,len(hists[0])),hists[0]-hists[2])
    ax[1,2].bar(np.linspace(0,231,len(hists[0])),hists[1]-hists[2])
    
    ax[2,0].bar(np.linspace(0,231,len(hists[0])),inact_long-hists[0])
    ax[2,1].bar(np.linspace(0,231,len(hists[0])),inact_long-hists[1])
    ax[2,2].bar(np.linspace(0,231,len(hists[0])),inact_long-hists[2])

    # ax[0,0].bar(np.arange(0,len(inact)),in_data[0]-in_data[1])
    # ax[0,2].bar(np.arange(0,len(inact)),in_data[0]-in_data[2])
    # ax[0,3].bar(np.arange(0,len(inact)),in_data[0]-in_data[3])
    # #ax[0,3].bar(np.arange(0,len(inact)),inact.mean(axis=1) - ((inact[1]-inact[2])+(inact[0]-inact[1])+(inact[0]-inact[2])))
    # ax[1,0].bar(a.columns,a.values[0])
    # ax[1,1].bar(np.arange(0,len(act)),act_data[0]-act_data[1])
    # ax[1,2].bar(np.arange(0,len(act)),act_data[0]-act_data[2])
    # ax[1,3].bar(np.arange(0,len(act)),act_data[0]-act_data[3])
    #ax[1,3].bar(np.arange(0,len(act)),act.mean(axis=1) - ((act[1]-act[2])+(act[0]-act[1])+(act[0]-act[2])))
    plt.tight_layout()
    # plt.savefig("States_Long_Inactive_dif_all.pdf")
    # plt.close()
    
long_state_diff()

def ternary_CG_leaflet(leaf):
    import choose_path as chp
    import pandas as pd
    uot = 0
    oot = 0
    fig,ax = plt.subplots(2,3,figsize=(12,8))
    uot,sm1 = cgp.Ternary_Heat_Map("CG/data/pi_eq_%s_inactive_short_map"%(leaf),"",ax=ax[0,1],out=uot)#,initial="CG/data/init_raw_%s_inactive_titrate_shortcg"%(leaf))
    cgp.Ternary_Heat_Map("CG/data/pi_raw_%s_inactive_short_map"%(leaf),"",ax=ax[0,0])#,initial="CG/data/init_raw_%s_inactive_longcg"%(leaf))
    oot,sm2 = cgp.Ternary_Heat_Map("CG/data/pi_eq_%s_inactive_short_map"%(leaf),"","CG/data/pi_raw_%s_inactive_short_map"%(leaf),ax=ax[0,2])
    cgp.Ternary_Heat_Map("CG/data/pi_eq_%s_inactive_titrate_short_map"%(leaf),"",ax=ax[1,1])#,initial="CG/data/init_raw_%s_inactive_shortcg"%(leaf))
    uot,sm1 = cgp.Ternary_Heat_Map("CG/data/pi_raw_%s_inactive_titrate_short_map"%(leaf),"",ax=ax[1,0], out=uot)#,initial="CG/data/pi_eq_%s_inactive_titrate_short_cgtmp"%(leaf))
    # cgp.Ternary_Scatter(ax[1,1], "CG/data/pi_slope")
    oot,sm2 = cgp.Ternary_Heat_Map("CG/data/pi_eq_%s_inactive_titrate_short_map"%(leaf),"","CG/data/pi_raw_%s_inactive_titrate_short_map"%(leaf),ax=ax[1,2],out=oot)
    # cgp.Ternary_Heat_Map("CG/data/pi_raw_%s_inactive_titrate_short_cgtmp"%(leaf),"","CG/data/pi_raw_%s_inactive_long_cg"%(leaf),ax=ax[2,0],out=oot)
    # cgp.Ternary_Heat_Map("CG/data/pi_eq_%s_inactive_titrate_short_cgtmp"%(leaf),"","CG/data/pi_raw_%s_inactive_short_cg"%(leaf),ax=ax[2,1],out=oot)

    cax = plt.axes([1, 0.25, 0.025, 0.5])
    plt.colorbar(sm2, cax=cax,format='%.3f')
    cax = plt.axes([0.15,-.025,0.45,0.025])
    plt.colorbar(sm1,cax=cax,format="%.3f",orientation="horizontal")
    plt.tight_layout()
    # plt.show()
    # plt.savefig("CG_ternary_short-long.pdf",bbox_inches='tight')
    # plt.close()
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
    # plt.show()
    # plt.savefig("Extended_Short_%s.pdf"%leaf,bbox_inches='tight')
    # plt.close()

# ternary_CG_leaflet("low") 

def ternary_CG():
    import choose_path as chp
    import pandas as pd
    uot = 0
    oot = 0
    fig,ax = plt.subplots(3,3,figsize=(12,8))
    uot,sm1 = cgp.Ternary_Heat_Map("CG/data/pi_eq_inactive_short_cg_titrate","",ax=ax[0,1],out=uot,initial="CG/data/init_raw_inactive_titrate_shortcg")
    cgp.Ternary_Heat_Map("CG/data/pi_raw_inactive_long_cg","",ax=ax[0,0])#,initial="CG/data/init_raw_inactive_longcg")
    oot,sm2 = cgp.Ternary_Heat_Map("CG/data/pi_eq_inactive_short_cg_titrate","","CG/data/pi_raw_inactive_long_cg",ax=ax[0,2])
    cgp.Ternary_Heat_Map("CG/data/pi_eq_inactive_short_cg_titrate","",ax=ax[1,1],initial="CG/data/init_raw_inactive_titrate_shortcg")
    uot,sm1 = cgp.Ternary_Heat_Map("CG/data/pi_raw_inactive_short_cg","",ax=ax[1,0], out=uot,initial="CG/data/init_raw_inactive_shortcg")
    # cgp.Ternary_Scatter(ax[1,1], "CG/data/pi_slope")
    oot,sm2 = cgp.Ternary_Heat_Map("CG/data/pi_eq_inactive_short_cg_titrate","","CG/data/pi_raw_inactive_short_cg",ax=ax[1,2],out=oot)
    cgp.Ternary_Heat_Map("CG/data/pi_raw_inactive_long_cg","","CG/data/pi_raw_inactive_short_cg",ax=ax[2,0],out=oot)
    cgp.Ternary_Heat_Map("CG/data/pi_eq_inactive_short_cg_titrate","","CG/data/pi_raw_inactive_short_cg_titrate",ax=ax[2,1],out=oot)

    cax = plt.axes([1, 0.25, 0.025, 0.5])
    plt.colorbar(sm2, cax=cax,format='%.3f')
    cax = plt.axes([0.15,-.025,0.45,0.025])
    plt.colorbar(sm1,cax=cax,format="%.3f",orientation="horizontal")
    plt.tight_layout()
    # plt.show()
    # plt.savefig("CG_ternary_short-long.pdf",bbox_inches='tight')
    # plt.close()
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
    # plt.show()
    plt.savefig("Extended_both_leafs.pdf",bbox_inches='tight')
    plt.close()
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
    cgp.plot_state_dist("~//CG/data/pi_eq_inactive_shortcg_titrate",ax[0,0])
    cgp.plot_state_dist("~//CG/data/pi_raw_inactive_shortcg_titrate",ax[0,1])
    cgp.diff_plot("~//CG/data/pi_eq_inactive_shortcg_titrate","~//CG/data/pi_raw_inactive_shortcg_titrate", ax[0,2])
    cgp.plot_state_dist("~//CG/data/pi_raw_inactive_longcg",ax[1,0])
    cgp.plot_state_dist("~//CG/data/pi_eq_inactive_shortcg_titrate",ax[1,1])
    cgp.diff_plot("~//CG/data/pi_raw_inactive_longcg","~//CG/data/pi_raw_inactive_shortcg_titrate", ax[1,2])
    plt.tight_layout()
    plt.savefig("CG_eq_state_dist_short-long-same2.pdf")
    plt.close()
# states_CG()


def tmp():
    import pandas as pd
    import matplotlib.pyplot as plt
    import numpy as np
    
    def make_same_len(sys,alps=aps.all_possible_states()):
        if len(sys.index) != len(alps) or len(sys.columns) != len(alps):
            sys_tmp = pd.DataFrame(columns=np.arange(0,len(alps)))
            sys_tmp[sys.index] = sys.T[sys.index]
            sys_tmp = sys_tmp.T.fillna(0)
            sys = sys_tmp.copy()
        return sys
    
    alps = aps.all_possible_states()
    act_pi = pd.read_csv("~//CG/data/pi_eq_active_shortcg_tmp.csv",index_col=0)
    act_pi = make_same_len(act_pi, alps)
    inact_pi = pd.read_csv("~//CG/data/pi_eq_inactive_shortcg_tmp.csv",index_col=0)
    inact_pi = make_same_len(inact_pi, alps)

    inact_raw = pd.read_csv("~//CG/data/pi_raw_inactive_shortcg.csv",index_col=0)
    act_raw = pd.read_csv("~//CG/data/pi_raw_active_shortcg.csv",index_col=0)
    
    ## make this part a function!
    
    # for sysi, sys in enumerate([act_pi,inact_pi, inact_raw, act_raw]):
    #     if len(sys.index) != len(alps) or len(sys.columns) != len(alps):
    #         sys_tmp = pd.DataFrame(columns=np.arange(0,len(alps)))
    #         sys_tmp[sys.index] = sys.T[sys.index]
    #         sys_tmp = sys_tmp.T.fillna(0)
    #         sys = sys_tmp.copy()
    #         print()
                
    # Cause it'll break otherwise!
    
    
    
    plt.bar(np.arange(0,len(act_pi)),act_pi.T.values[0])
    plt.bar(np.arange(0,len(act_pi)),act_raw.T.values[0],alpha=.5)
    plt.bar(np.arange(0,len(act_pi)),act_pi.T.values[0]-act_raw.T.values[0])

    # plt.savefig("tmp.pdf")
    # plt.close()
# tmp()

def state_iterative():
    
    import choose_path as chp
    import pandas as pd
    import numpy as np

    state_act = pd.read_csv("%s/CG/data/act_short_binned_sim_pi.csv"%(chp.choose_path()[1]),index_col=0).T
    state_inact = pd.read_csv("%s/CG/data/inact_short_binned_sim_pi.csv"%(chp.choose_path()[1]),index_col=0).T
    
    
    long_act = pd.read_csv("%s/CG/data/pi_raw_active_shortcg_sim_seed.csv"%(chp.choose_path()[1]),index_col=0).T
    long_inact = pd.read_csv("%s/CG/data/pi_raw_inactive_shortcg_sim_seed.csv"%(chp.choose_path()[1]),index_col=0).T
    
    # dsi,dss = 0,0
    
    for si in state_inact:

        fig, ax = plt.subplots(2,3,figsize=(8,6),sharey=True,sharex=True)
        left = .0
        bottom, height = .0, .95
        top = bottom + height
        ax[0,0].text(left, top, '%f'%((si+1)*1),
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
        plt.savefig("state-short-short_sim%s_seed.pdf"%si)
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
test1 = cgc.CGTM_Calculations("",1,"cg","inactive","short")#.sigConverge_time()
a = test1.build_CGTM(symitrize=False)#write_pi_eq()
b = test1.build_raw()#write_pi_raw()
pi_raw = b[0]
pi_eq = a[0]
cgtm = a[-1].values

# X = []
# for i in range(0,len(pi_raw)):
#     xi = []
#     for j in range(0,len(pi_raw)):
#         xi.append(pi_eq[i]*cgtm[i,j])
#     X.append(xi)

# X = np.array(X)
# assert np.allclose(X.sum(axis=1),pi_eq,1E-2,1E-10) or np.abs(X.sum(axis=1)-pi_eq).mean()<1E-4, "pi_eq and X_i should the equal"
# assert np.all(X.sum(axis=1) >= 0 ), "X_i is less than 0"
# assert np.any(X[X<0])==False, "Negative values detected"

# tmp = []
# tmp2 = []
# for i in range(0,len(pi_raw)):
#     tmp_i = []
#     tmp_j= []
#     for j in range(0,len(pi_raw)):
#         if cgtm[i,j] == 0: 
#             tmp_i.append(0)
#             tmp_j.append(0)
#             continue
#         heh=X[i,j]/X[i].sum()
#         tmp_i.append(cgtm[i,j] * np.log(heh))
#         tmp_j.append((cgtm[i,j]+cgtm[i,j])/X[i,j] - cgtm[i].sum()/X[i].sum() - cgtm[j].sum()/X[j].sum())
#     tmp2.append(tmp_j)
#     tmp.append(tmp_i)
    
# P_rev = []
# for i in range(0,len(pi_raw)):
#     P_tmp = []
#     for j in range(0,len(pi_raw)):
#         if cgtm[i,j] == 0:
#             P_tmp.append(0)
#             continue
#         P_tmp.append((cgtm[i,j] + cgtm[j,i])*pi_eq[j] / (cgtm[i].sum()*pi_eq[j] + cgtm[j].sum()*pi_eq[i]))
#     P_rev.append(P_tmp)
# P_rev = np.array(P_rev)

# tmp = np.sum(tmp)
# tmp2 = np.array(tmp2)
# print(np.sqrt(np.abs(tmp2@tmp2)).mean())
# tmp2[tmp2==0] = np.nan
# plt.pcolormesh(tmp2,vmax=5,vmin=-5)
# plt.colorbar()

def IDK():

    import numpy as np
    test1 = cgc.CGTM_Calculations("",1,"cg","inactive","short")#.sigConverge_time()
    # # # test2 = cgc.CGTM_Calculations("",1,"cg","inactive","short").sigConverge_simulations()
    # # # sig_conv(test1,test2,'cg')
    
    def make_same_len(sys,alps=aps.all_possible_states()):
        import pandas as pd
        if len(sys.index) != len(alps) or len(sys.columns) != len(alps):
            sys_tmp = pd.DataFrame(columns=np.arange(0,len(alps)))
            sys_tmp[sys.index] = sys.T[sys.index]
            sys_tmp = sys_tmp.T.fillna(0)
            sys = sys_tmp.copy()
        return sys
        
    
    # # test1 = cgc.CGTM_Calculations("SU", 1, "chg",act=None)#.sigConverge_time()
    a = test1.build_CGTM(symitrize=True)#write_pi_eq()
    b = test1.build_raw()#write_pi_raw()
    
    pi_raw = b[0]
    pi_eq = a[0]
    
    # pi_eq = make_same_len(pi_eq)
    cgtm = a[-1]
    hold = []
    alps = aps.all_possible_states()
    DEV = []
    dev_thresh = [0,0.05,0.1,.15,0.2,0.25,1]
    fig,ax = plt.subplots(2,len(dev_thresh),figsize=(5*len(dev_thresh),7.5))
    for ti, dt in enumerate(dev_thresh):
        empt = []
    
        f = open("inactive%.2f.dat"%dt,'w')
        for i in range(0,len(cgtm)):
            for j in range(0,len(cgtm)):
        
                if cgtm.iloc[i,j] == 0 and cgtm.iloc[j,i] == 0:
                    continue        
        
                if cgtm.iloc[i,j] == 1 or  cgtm.iloc[j,i] == 1:
                    if cgtm.iloc[i,j] == 1 and  cgtm.iloc[j,i] == 1:
                        print("sink: ",i,j)
                    else:
                        print("state ",cgtm.index[i],alps[cgtm.index[i]],cgtm.iloc[i,j], "-> states ",cgtm.index[j],alps[cgtm.index[j]])
        
        
                dev = np.sqrt(((pi_eq[i] *cgtm.iloc[i,j])-(pi_eq[j] * cgtm.iloc[j,i]))**2) / ((0.5*pi_eq[i] *cgtm.iloc[i,j]) + (0.5*pi_eq[j] * cgtm.iloc[j,i]))
                DEV.append(dev)
                hold.append([i,j,dev,pi_eq[i]*cgtm.iloc[i,j],pi_eq[j]*cgtm.iloc[j,i]])
                # if np.isclose((pi_eq[i] *cgtm[i,j])/(pi_eq[j] * cgtm[j,i]),1)==True:#np.isclose(pi_eq[i] *cgtm[i,j],pi_eq[i] * cgtm[j,i])==True:
                #     empt.append([(i,j),pi_eq[i]*cgtm[i,j],pi_eq[j]*cgtm[j,i]])
                if dev <= dt:
                    empt.append([(i,j),pi_eq[i]*cgtm.iloc[i,j],pi_eq[j]*cgtm.iloc[j,i]])
                else:
                    empt.append(["x",(i,j),pi_eq[i]*cgtm.iloc[i,j],pi_eq[j]*cgtm.iloc[j,i]])
                    f.write("state %s:%s (%s) <-> (%s) states:%s:%s  \n"%(str(cgtm.index[i]),str(alps[cgtm.index[i]]),str(pi_eq[i]*cgtm.iloc[i,j]),str(pi_eq[j]*cgtm.T.iloc[i,j]),str(cgtm.index[j]),str(alps[cgtm.index[j]])))
    
        f.close()
    
        tmp_len = []
        [tmp_len.append(len(i)) for i in empt]
        h,e = np.histogram(tmp_len,range=(3,4))
    
        if h[0] <= h[-1]:
            print("pi_{eq} @ T = pi_{eq}: ",np.allclose(pi_eq@cgtm,pi_eq))
        ax[1,ti].bar(e[:-1],h/h.sum())
        ax[1,ti].set_xticks([3,4],["ratio reversable","ratio ireversable"])
        
    
    
    
        cgtm_mod = np.zeros_like(cgtm)
        for i in range(0,len(cgtm)):
            for j in range(0,len(cgtm)):
                #if i==j: continue
                dev = np.sqrt(((pi_eq[i] *cgtm.iloc[i,j])-(pi_eq[j] * cgtm.iloc[j,i]))**2) / ((0.5*pi_eq[i] *cgtm.iloc[i,j]) + (0.5*pi_eq[j] * cgtm.iloc[j,i]))
                if dev <= dt:
                    cgtm_mod[i,j] = 5
                elif dev > dt and dev < 1 and dev < 1:
                    cgtm_mod[i,j] = 2.5
                elif dev >= 1:
                    cgtm_mod[i,j] = 0
                else:
                    cgtm_mod[i,j] = np.nan
        ax[0,ti].pcolormesh(cgtm_mod)
        ax[0,ti].set_title("deviation=%0.2f"%dt)
    
    # ax[2].pcolormesh(-np.log(cgtm))
    cgtm[cgtm == 0] = np.nan
    # ax[0,0].pcolormesh(np.log(cgtm))#,cmap="plasma_r")
    ax[0,0].set_title("CGTM")
    ax[0,0].pcolormesh(cgtm)
    plt.tight_layout()
    plt.savefig("cgtm_dev_seed_sym.pdf")
    plt.close()

def state_tracker():
    import pandas as pd
    data = pd.read_csv("/home/liam/lms464/CG/data/states/low_short_inactive_tmp.csv",index_col=0).values
    
    fig, ax = plt.subplots(len(data),1,figsize=(5,2.5*len(data)))
    for i,dat in enumerate(data):
        ax[i].plot(dat)
        ax[i].set_ylim(0,231)
    
    plt.savefig("state_short_change.pdf")
    plt.close()
# state_tracker()
# # pi_rand = pi_raw.copy()
# # pi_rand[pi_rand > 0] = 1
# pi_rand = np.random.rand(len(pi_eq))*np.random.rand(len(pi_eq))
# pi_rand = pi_rand / pi_rand.sum()
# pi_rand_tmp = pi_rand.copy()
# pi_rand_init = pi_rand.copy()
# max_ind = np.where(pi_eq==pi_eq.max())[0][0]
# i = 1
# fig,ax = plt.subplots(1,2)
# while np.allclose(pi_rand, pi_eq) == False:
#     ax[0].plot(i,pi_eq.iloc[max_ind],'ko')
#     ax[0].plot(i,pi_rand[max_ind],'bo')
#     ax[0].plot(i,pi_rand_tmp[max_ind],'ro')
#     pi_rand_tmp = (pi_rand_tmp @ cgtm).values  
#     pi_rand = pi_rand_tmp/pi_rand_tmp.sum()
#     i = i + 1
#     print(i,pi_rand_tmp.sum())
# ax[1].imshow(cgtm)
    
# # plt.ylim([0,1.5])
# ax[0].set_title("Markov 'simulation'")
# # ax[0].legend([r"max($\pi^{eq}$)","Same state Sim (Normalized)","Same state Sim Un-Normalized"])
# ax[0].set_xlabel(r"Steps to converge to $\pi^{eq}$")
# ax[0].set_ylim(-.001,.25)
# ax[0].set_ylabel("Sum of states (should be black line)")
# plt.savefig("Markov_Sim_non-Zero_CGTM_Sym.pdf")
# plt.close()

    

# test2 = cgc.CGTM_Calculations("SL", 1, "chg",act=None).sigConverge_time()
# sig_conv(test1,test2,'chg')

# network(["SU","SL"],"sat",None)
# import numpy as np
# test1 = cgc.CGTM_Calculations("",1,"cg","active","short").weighted_avg()
# fig, ax = plt.subplots(1,1)
# oot = 0
# oot,sm2 = cgp.Ternary_Heat_Map("CG/data/pi_eq_active_shortcg","","CG/data/pi_raw_active_shortcg",ax=ax,out=oot)
# plt.colorbar(sm2)