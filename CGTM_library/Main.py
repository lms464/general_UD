import matplotlib.pyplot as plt
import CGTM_Calculations as cgc
import CGTM_Plotting as cgp

def ternary_CG():
    uot = 0
    oot = 0
    fig,ax = plt.subplots(2,3,figsize=(12,8))
    cgp.Ternary_Heat_Map("CG/data/pi_eq_inactive_shortcg","",ax=ax[0,0])
    cgp.Ternary_Heat_Map("CG/data/pi_raw_inactive_longcg","",ax=ax[0,1])
    cgp.Ternary_Heat_Map("CG/data/pi_eq_inactive_shortcg","","CG/data/pi_raw_inactive_shortcg",ax=ax[0,2])
    cgp.Ternary_Heat_Map("CG/data/pi_eq_active_shortcg","",ax=ax[1,0])
    uot = cgp.Ternary_Heat_Map("CG/data/pi_raw_active_shortcg","",ax=ax[1,1], out=uot)
    # cgp.Ternary_Scatter(ax[1,2], "CG/data/pi_eq_active_shortcg")
    oot = cgp.Ternary_Heat_Map("CG/data/pi_eq_active_shortcg","","CG/data/pi_raw_active_shortcg",ax=ax[1,2],out=oot)
    cax = plt.axes([1, 0.25, 0.025, 0.5])
    plt.colorbar(oot, cax=cax,format='%.3f')
    cax = plt.axes([0.15,-.025,0.45,0.025])
    plt.colorbar(uot,cax=cax,format="%.3f",orientation="horizontal")
    plt.tight_layout()
    #plt.show()
    plt.savefig("CG_ternary.pdf",bbox_inches='tight')
    plt.close()

    
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
    cgp.plot_state_traj("~/UDel/CG/data/states/short_inactive",ax[0,0])
    cgp.plot_state_traj("~/UDel/CG/data/states/short_active",ax[0,1])
    cgp.plot_state_traj("~/UDel/CG/data/states/long_inactive",ax[1,0])
    cgp.plot_state_traj("~/UDel/CG/data/states/long_active",ax[1,1])
    plt.tight_layout()
    plt.savefig("CG_SI.pdf")
    plt.close()
ternary_CG()
# states_CG()   
# CGTM()    
# SI_CG()

def Confidence():
    path="/home/liam/Censere/UDel/CG/data"
    fig, ax = plt.subplots(1,2,figsize=(12,6),sharex="row")
    for jj, j in enumerate(["inactive","active"]):
        cgp.plot_convergence(path+"/"+j+"_ratio_sum.csv",ax[jj])
    
    plt.savefig("confidence.pdf")
    plt.close()
# Confidence()