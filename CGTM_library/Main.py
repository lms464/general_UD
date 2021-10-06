import matplotlib.pyplot as plt
import CGTM_Calculations as cgc
import CGTM_Plotting as cgp

def ternary_CG():
    uot = 0
    oot = 0
    fig,ax = plt.subplots(2,3,figsize=(12,8))
    cgp.Ternary_Heat_Map("CG/data/pi_eq_inactive_shortcg","",ax=ax[0,0],initial="CG/data/init_raw_inactive_shortcg")
    cgp.Ternary_Heat_Map("CG/data/pi_raw_inactive_longcg","",ax=ax[0,1])
    cgp.Ternary_Heat_Map("CG/data/pi_eq_inactive_shortcg","","CG/data/pi_raw_inactive_shortcg",ax=ax[0,2])
    cgp.Ternary_Heat_Map("CG/data/pi_eq_active_longcg","",ax=ax[1,0] ,initial="CG/data/init_raw_inactive_shortcg")
    uot,sm1 = cgp.Ternary_Heat_Map("CG/data/pi_raw_active_shortcg","",ax=ax[1,1], out=uot)
    # cgp.Ternary_Scatter(ax[1,2], "CG/data/pi_eq_active_shortcg")
    oot,sm2 = cgp.Ternary_Heat_Map("CG/data/pi_eq_active_shortcg","","CG/data/pi_raw_active_shortcg",ax=ax[1,2],out=oot)
    cax = plt.axes([1, 0.25, 0.025, 0.5])
    plt.colorbar(sm2, cax=cax,format='%.3f')
    cax = plt.axes([0.15,-.025,0.45,0.025])
    plt.colorbar(sm1,cax=cax,format="%.3f",orientation="horizontal")
    plt.tight_layout()
    # plt.show()
    plt.savefig("CG_ternary.pdf",bbox_inches='tight')
    plt.close()

def ternary_scatter():
    #Needs a functin to look at initial distribution, not 
    # the FINAL distribution!
    
    import pandas as pd
    import ternary
    figure, ax = ternary.figure(scale=1)
    figure.set_size_inches(10, 7.5)
    states = pd.read_csv("all_states.csv",index_col=0)
    cgp.Ternary_Scatter(ax,"CG/data/pi_eq_active_shortcg",True)
# ternary_scatter()
    
def states_CG():
    fig, ax = plt.subplots(2,3,figsize=(8,6),sharey='col',sharex=True)
    cgp.plot_state_dist("~/UDel/CG/data/pi_eq_inactive_shortcg",ax[0,0])
    cgp.plot_state_dist("~/UDel/CG/data/pi_raw_inactive_longcg",ax[0,1])
    cgp.diff_plot("~/UDel/CG/data/pi_eq_inactive_shortcg","~/UDel/CG/data/pi_raw_inactive_longcg", ax[0,2])
    cgp.plot_state_dist("~/UDel/CG/data/pi_eq_active_shortcg",ax[1,0])
    cgp.plot_state_dist("~/UDel/CG/data/pi_raw_active_longcg",ax[1,1])
    cgp.diff_plot("~/UDel/CG/data/pi_eq_active_shortcg","~/UDel/CG/data/pi_raw_active_longcg", ax[1,2])
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
ternary_CG()
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
  


# # test1 = cgc.CGTM_Calculations("",1,"cg","inactive","short").sigConverge_time_diff()
# test2 = cgc.CGTM_Calculations("",1,"cg","active","short").sigConverge_time_diff()
# sig_conv(test2,test2,'cg')

# network(["SU","SL"],"charge")
# import numpy as np
# test1 = cgc.CGTM_Calculations("",1,"cg","active","short").weighted_avg()
# fig, ax = plt.subplots(1,1)
# oot = 0
# oot,sm2 = cgp.Ternary_Heat_Map("CG/data/pi_eq_active_shortcg","","CG/data/pi_raw_active_shortcg",ax=ax,out=oot)
# plt.colorbar(sm2)