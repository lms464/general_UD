import matplotlib.pyplot as plt
import CGTM_Calculations as cgc
import CGTM_Plotting as cgp

fig,ax = plt.subplots(2,3,figsize=(8,6))
cgp.Ternary_Heat_Map("CG/data/pi_eq_inactive_shortcg","",ax=ax[0,0])
cgp.Ternary_Heat_Map("CG/data/pi_raw_inactive_longcg","",ax=ax[0,1])
cgp.Ternary_Heat_Map("CG/data/pi_eq_inactive_shortcg","","CG/data/pi_raw_inactive_shortcg",ax=ax[0,2])
cgp.Ternary_Heat_Map("CG/data/pi_eq_active_shortcg","",ax=ax[1,0])
cgp.Ternary_Heat_Map("CG/data/pi_raw_active_shortcg","",ax=ax[1,1])
cgp.Ternary_Scatter(ax[1,2], "CG/data/pi_eq_active_shortcg")
# cgp.Ternary_Heat_Map("CG/data/pi_eq_active_shortcg","","CG/data/pi_raw_active_shortcg",ax=ax[1,2])

plt.show()