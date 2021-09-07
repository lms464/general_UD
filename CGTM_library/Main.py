import matplotlib.pyplot as plt
import CGTM_Calculations as cgc
import CGTM_Plotting as cgp

def ternary_CG():
    fig,ax = plt.subplots(2,3,figsize=(8,6))
    cgp.Ternary_Heat_Map("CG/data/pi_eq_inactive_shortcg","",ax=ax[0,0])
    cgp.Ternary_Heat_Map("CG/data/pi_raw_inactive_longcg","",ax=ax[0,1])
    cgp.Ternary_Heat_Map("CG/data/pi_eq_inactive_shortcg","","CG/data/pi_raw_inactive_shortcg",ax=ax[0,2])
    cgp.Ternary_Heat_Map("CG/data/pi_eq_active_shortcg","",ax=ax[1,0])
    cgp.Ternary_Heat_Map("CG/data/pi_raw_active_shortcg","",ax=ax[1,1])
    # cgp.Ternary_Scatter(ax[1,2], "CG/data/pi_eq_active_shortcg")
    cgp.Ternary_Heat_Map("CG/data/pi_eq_active_shortcg","","CG/data/pi_raw_active_shortcg",ax=ax[1,2])
    plt.savefig("CG_ternary.pdf")
    
def states_CG():
    fig, ax = plt.subplots(2,3,figsize=(8,6),sharey='col',sharex=True)
    cgp.plot_state_dist("~/CG/data/pi_eq_inactive_shortcg",ax[0,0])
    cgp.plot_state_dist("~/CG/data/pi_raw_inactive_longcg",ax[0,1])
    cgp.diff_plot("~/CG/data/pi_eq_inactive_shortcg","~/CG/data/pi_raw_inactive_longcg", ax[0,2])
    cgp.plot_state_dist("~/CG/data/pi_eq_active_shortcg",ax[1,0])
    cgp.plot_state_dist("~/CG/data/pi_raw_active_longcg",ax[1,1])
    cgp.diff_plot("~/CG/data/pi_eq_active_shortcg","~/CG/data/pi_raw_active_longcg", ax[1,2])
    plt.tight_layout()
    plt.savefig("CG_state_dist")
    
def SI_CG():
    fig, ax = plt.subplots(2,2,figsize=(8,6),sharey='row',sharex="row")
    cgp.plot_state_traj("~/CG/data/states/short_inactive",ax[0,0])
    cgp.plot_state_traj("~/CG/data/states/short_active",ax[0,1])
    cgp.plot_state_traj("~/CG/data/states/long_inactive",ax[1,0])
    cgp.plot_state_traj("~/CG/data/states/long_active",ax[1,1])
    plt.tight_layout()
    plt.show()

def Confidence():
    for i in ["inactive","active"]:
        con = cgc.CGTM_Calculations("",1,"cg",i,"long")
        con.calc_confidence()
Confidence()
