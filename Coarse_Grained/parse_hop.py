import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import itertools 
import scipy.sparse.linalg as sla


def all_possible_states():
    # RUNS, seems fine
    #sorts by anionic, then neutral, then chol
    list_o_list = [np.linspace(0,1,21),np.linspace(0,1,21),np.linspace(0,1,21)]
    all_possibility = list(itertools.product(*list_o_list))
        
    list_100 = []
     
    for li, l in enumerate(all_possibility):
        if np.isclose(np.sum(l),1.0,1e-3,1e-5) == False:
            # print(li,np.sum(l))
            continue
        list_100.append(l)
    list_100 = np.array([*list_100])
    #return list_100[np.argsort(list_100[:,2])][::-1]
    return list_100[np.lexsort((list_100[:,0], list_100[:,1], list_100[:,2]))][::-1]

def check_states(data_frm, possible_states):
    
    data_frm = np.asarray(data_frm)/np.sum(data_frm)
    holder = []
    for pi, ps in enumerate(possible_states):
        
        rmsd = np.sqrt((data_frm - ps)**2)
        holder.append(rmsd.sum())

    out = np.argmin(holder)
    
    return out
  
def analysis(path,out):
    fl = open (path+"/borders.txt")
    lines = fl.readlines()#.split()
    fl.close()
    #lines = [int(l) for l in lines]
    shell = []
    for l in lines[1::3]:
    	shell.append([int(l.split()[-3]),int(l.split()[-2]),int(l.split()[-1])])
    shell_arr = np.array(shell)
    
    possible_states = all_possible_states()
    
    states = []
    
    for s in shell:
    	states.append(check_states(s,possible_states))
    states = np.array(states)
    hist, edge = np.histogram(states,range(0,len(possible_states)+1),normed=True)
    plt.bar(edge[:-1],hist)
    plt.show()
    np.savetxt("states.txt",states,fmt="%i")
    pd.DataFrame(hist).to_csv("states.csv")
    dt = 500
    shell_con = (shell_arr.T/ shell_arr.sum(axis = 1)).T
    dppc,dp_var = [],[]
    chol,ch_var = [],[]
    dopc,do_var = [],[]
    for frm in range(0,len(shell_con),dt):
        dppc.append(shell_con[frm:frm+dt,0].mean()), dp_var.append(shell_con[frm:frm+dt,0].var())
        chol.append(shell_con[frm:frm+dt,1].mean()), ch_var.append(shell_con[frm:frm+dt,1].var())
        dopc.append(shell_con[frm:frm+dt,2].mean()), do_var.append(shell_con[frm:frm+dt,2].var())
    
    #plt.plot(range(0,len(shell_con),dt),dppc,'b*-')
    plt.errorbar(np.array(range(0,len(shell_con),dt))/.5, dppc, yerr=dp_var, fmt='bo-',label="DPPC")
    plt.errorbar(np.array(range(0,len(shell_con),dt))/.5,chol,yerr=ch_var,fmt="go-",label="CHOL")
    plt.errorbar(np.array(range(0,len(shell_con),dt))/.5,dopc,yerr=do_var,fmt="ro-",label="DOPC")
    plt.ylim(0,1)
    plt.ylabel("Concentration")
    plt.xlabel("t(ns)")
    plt.legend()
    plt.savefig("%s_Conv.pdf"%out)
    plt.close()
    print("DPPC: %f"%np.mean(dppc[int(len(dppc)/2):]))
    print("DOPC: %f"%np.mean(dopc[int(len(dppc)/2):]))
    print("CHOL: %f"%np.mean(chol[int(len(dppc)/2):]))
    return hist,edge

def build_TM(states):
    all_states = all_possible_states()
    TM_m = np.zeros((len(all_states),len(all_states)))
    # norm_t = []
    
    if np.ndim(states) == 1:
        for si in range(1, len(states)):
            TM_m[int(states[si]),int(states[si-1])] += 1 
    else:
        for S in states:
            for si in range(1, len(states[S])):
                TM_m[int(states[S][si]),int(states[S][si-1])] += 1 
    TM_sym = 1/2 * (TM_m + TM_m.T)
    norm = TM_m.sum(axis=1)
    TM_norm = np.zeros(np.shape(TM_sym))
    for i,j in enumerate(TM_m):
        TM_norm[i] = j / norm[i]
    
    #TM_norm = np.divide(TM_m, norm)
    TM_norm = np.nan_to_num(TM_norm)
    # TM_test = np.nan_to_num(np.divide(TM_m, norm))
    return TM_norm    


def solve_pi_eq(P):
    evals, evecs = sla.eigs(P.T, k = 1, which='LM')
    evecs = np.real(evecs)
    pi_eg = (evecs/evecs.sum()).real
    pi_eg[pi_eg < 10**(-7)] = 0 #TODO wait... why?
    pi_EG = np.array([e[0] for e in pi_eg])
    pi_eg = pi_EG
    return pi_eg, np.linalg.eig(P.T)            

def build_CGTM(path):

    states = pd.read_csv("%s/%s%s.csv"%(path),index_col=0).T   
    TM_norm = build_TM(states)
    pi_eq, eigs = solve_pi_eq(TM_norm)
    pd.DataFrame(pi_eq).to_csv("%s/pi_%s_%s.csv"%(path))
    pd.DataFrame(TM_norm).to_csv("%s/CGTM_%s_%s.csv"%(path))
    #return pi_eq, eigs, TM_norm


def analysis_multi_cgtm(path,init,fint):
    
    full_states = []
    for i in range(init,fint+1):
        fl = open (path+"%i/borders.txt"%i)
        lines = fl.readlines()#.split()
        fl.close()
        #lines = [int(l) for l in lines]
        shell = []
        for l in lines[1::3]:
        	shell.append([int(l.split()[-3]),int(l.split()[-2]),int(l.split()[-1])])
        shell_arr = np.array(shell)
        
        possible_states = all_possible_states()
        
        states = []
        
        for s in shell:
        	states.append(check_states(s,possible_states))
        full_states.append(states)
        
    TM_norm = build_TM(pd.DataFrame(full_states))
    pi_eq, eigs = solve_pi_eq(TM_norm)
    plt.bar(range(0,len(all_possible_states()),len(all_possible_states())),pi_eq)
    plt.show()
    np.savetxt("states_short_raw.txt",states,fmt="%i")
    pd.DataFrame(pi_eq).to_csv("states_short_raw.csv")
    print("DPPC:  %f"%(pi_eq * all_possible_states()[:,0]).sum())
    print("DOPC:  %f"%(pi_eq * all_possible_states()[:,1]).sum())
    print("CHOL:  %f"%(pi_eq * all_possible_states()[:,2]).sum())
    return pi_eq
    

def analysis_multi_raw(path,init,fint):
    
    for i in range(init,fint+1):
        fl = open (path+"%i/borders.txt"%i)
        lines = fl.readlines()#.split()
        fl.close()
        #lines = [int(l) for l in lines]
        shell = []
        for l in lines[1::3]:
        	shell.append([int(l.split()[-3]),int(l.split()[-2]),int(l.split()[-1])])
        shell_arr = np.array(shell)
        
        possible_states = all_possible_states()
        
        states = []
        
        for s in shell:
        	states.append(check_states(s,possible_states))
        states = np.array(states)
    hist, edge = np.histogram(states,range(0,len(possible_states)+1),normed=True)
    plt.bar(edge[:-1],hist)
    plt.show()
    np.savetxt("states_short_raw.txt",states,fmt="%i")
    pd.DataFrame(hist).to_csv("states_short_raw.csv")
    # dt = 500
    # shell_con = (shell_arr.T/ shell_arr.sum(axis = 1)).T
    # dppc,dp_var = [],[]
    # chol,ch_var = [],[]
    # dopc,do_var = [],[]
    # for frm in range(0,len(shell_con),dt):
    #     dppc.append(shell_con[frm:frm+dt,0].mean()), dp_var.append(shell_con[frm:frm+dt,0].var())
    #     chol.append(shell_con[frm:frm+dt,1].mean()), ch_var.append(shell_con[frm:frm+dt,1].var())
    #     dopc.append(shell_con[frm:frm+dt,2].mean()), do_var.append(shell_con[frm:frm+dt,2].var())
    print("DPPC:  %f"%(hist * all_possible_states()[:,0]).sum())
    print("DOPC:  %f"%(hist * all_possible_states()[:,1]).sum())
    print("CHOL:  %f"%(hist * all_possible_states()[:,2]).sum())
    return hist,edge

pi_eq = analysis_multi_cgtm("/home/sharplm/shot_inactive",2,10)
hist,edge=analysis_multi_raw("/home/sharplm/shot_inactive",2,10)
#hist,edge = analysis("/home/sharplm/Inactive3","inactive")

dif = pi_eq - hist
plt.bar(edge[:-1],pi_eq)
plt.savefig("Inactive_pi_eq.pdf")
plt.close()
plt.bar(edge[:-1],hist)
plt.savefig("Inactive_Long.pdf")
plt.close()
plt.bar(edge[:-1],dif)
plt.savefig("Inactive_diff.pdf")
plt.close()