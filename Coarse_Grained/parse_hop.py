import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import itertools 


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

analysis("/home/sharplm/Active3","active")
analysis("/home/sharplm/Inactive3","inactive")