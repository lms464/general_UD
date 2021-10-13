import matplotlib.pyplot as plt
import CGTM_Calculations as cgc
import CGTM_Plotting as cgp


def ternary_CG():
    uot = 0
    oot = 0
    fig,ax = plt.subplots(2,3,figsize=(12,8))
    cgp.Ternary_Heat_Map("CG/data/pi_eq_inactive_shortcg","",ax=ax[0,1],initial="CG/data/init_raw_inactive_shortcg")
    cgp.Ternary_Heat_Map("CG/data/pi_raw_inactive_longcg","",ax=ax[0,0])
    cgp.Ternary_Heat_Map("CG/data/pi_eq_inactive_shortcg","","CG/data/pi_raw_inactive_longcg",ax=ax[0,2])
    cgp.Ternary_Heat_Map("CG/data/pi_eq_active_longcg","",ax=ax[1,1] ,initial="CG/data/init_raw_inactive_shortcg")
    uot,sm1 = cgp.Ternary_Heat_Map("CG/data/pi_raw_active_shortcg","",ax=ax[1,0], out=uot)
    # cgp.Ternary_Scatter(ax[1,2], "CG/data/pi_eq_active_shortcg")
    oot,sm2 = cgp.Ternary_Heat_Map("CG/data/pi_eq_active_shortcg","","CG/data/pi_raw_active_longcg",ax=ax[1,2],out=oot)
    cax = plt.axes([1, 0.25, 0.025, 0.5])
    plt.colorbar(sm2, cax=cax,format='%.3f')
    cax = plt.axes([0.15,-.025,0.45,0.025])
    plt.colorbar(sm1,cax=cax,format="%.3f",orientation="horizontal")
    plt.tight_layout()
    # plt.show()
    plt.savefig("CG_ternary.pdf",bbox_inches='tight')
    plt.close()

def ternary_iterative():
    import choose_path as chp
    import pandas as pd
    uot = 0
    oot = 0
    
    long_act = pd.read_csv("%s/CG/data/pi_raw_inactive_longcg.csv"%(chp.choose_path()[1]),index_col=0).T
    long_inact = pd.read_csv("%s/CG/data/pi_raw_active_longcg.csv"%(chp.choose_path()[1]),index_col=0).T
    
    
    state_act = pd.read_csv("%s/CG/data/act_short_binned_pi.csv"%(chp.choose_path()[1]),index_col=0).T
    state_inact = pd.read_csv("%s/CG/data/inact_short_binned_pi.csv"%(chp.choose_path()[1]),index_col=0).T
    
    dpi_ref = (long_inact - state_inact[0].values).iloc[0,:]
    import numpy as np
    for si, sj in zip(state_act[-2:-1], state_inact[-2:-1]):
        
        fig,ax = plt.subplots(2,3,figsize=(8,8))
        
        cgp.Ternary_Heat_Map(state_inact[si],fl_name="",ax=ax[1,1],out=uot)
        cgp.Ternary_Heat_Map("CG/data/pi_raw_inactive_longcg","",ax=ax[1,0])
        cgp.Ternary_Heat_Map(long_inact, leaflet_in2=state_inact[si].values,fl_name="",ax=ax[1,2],out=oot)

        cgp.Ternary_Heat_Map(state_act[sj],fl_name="",ax=ax[0,1])
        cgp.Ternary_Heat_Map("CG/data/pi_raw_active_shortcg","",ax=ax[0,0])
        cgp.Ternary_Heat_Map(long_act, leaflet_in2=state_act[sj].values,fl_name="",ax=ax[0,2],out=oot)

        plt.tight_layout()
        plt.savefig("ternary_gif%i.png"%si)
        plt.close()
    
def states_CG():
    fig, ax = plt.subplots(2,3,figsize=(8,6),sharey='col',sharex=True)
    cgp.plot_state_dist("~/UDel/CG/data/pi_eq_inactive_shortcg",ax[0,0])
    cgp.plot_state_dist("~/UDel/CG/data/pi_raw_inactive_shortcg",ax[0,1])
    cgp.diff_plot("~/UDel/CG/data/pi_eq_inactive_shortcg","~/UDel/CG/data/pi_raw_inactive_shortcg", ax[0,2])
    cgp.plot_state_dist("~/UDel/CG/data/pi_eq_active_shortcg",ax[1,0])
    cgp.plot_state_dist("~/UDel/CG/data/pi_raw_active_shortcg",ax[1,1])
    cgp.diff_plot("~/UDel/CG/data/pi_eq_active_shortcg","~/UDel/CG/data/pi_raw_active_shortcg", ax[1,2])
    plt.tight_layout()
    plt.savefig("CG_state_dist.pdf")
    plt.close()

    
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
# # states_CG()   
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
    
    path="~/CG/data"
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
    plt.savefig("State_transitions.pdf")
    plt.close()

def sig_conv(SL,SU,kind):
    cgp.plot_sigConverge(SL,SU,kind)
  


# test1 = cgc.CGTM_Calculations("",1,"cg","inactive","short").sigConverge_simulations()[0]
# test2 = cgc.CGTM_Calculations("",1,"cg","active","short").sigConverge_simulations()[0]
# sig_conv(test1,test2,'cg')

# test1 = cgc.CGTM_Calculations("SU", 1, "chg",act=None).sigConverge_time()
# test2 = cgc.CGTM_Calculations("SL", 1, "chg",act=None).sigConverge_time()
# sig_conv(test1,test2,'chg')

network(["SU","SL"],"sat",None)
# import numpy as np
# test1 = cgc.CGTM_Calculations("",1,"cg","active","short").weighted_avg()
# fig, ax = plt.subplots(1,1)
# oot = 0
# oot,sm2 = cgp.Ternary_Heat_Map("CG/data/pi_eq_active_shortcg","","CG/data/pi_raw_active_shortcg",ax=ax,out=oot)
# plt.colorbar(sm2)