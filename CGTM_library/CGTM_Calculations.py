#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 17 17:43:09 2021

@author: liam
"""
import pandas as pd
import numpy as np
import scipy.sparse.linalg as sla
import all_possible_states as aps
import choose_path as chp


class CGTM_Calculations:
    def __init__(self, leaflet_in, dt, kind, act=None, length=None):
        
        assert dt > 0, "dt MUST be greater than 0"
        
        self.leaflet_in = leaflet_in
        self.dt = dt
        self.kind = None
        self.path = chp.choose_path()
        self.act = act
        self.length = length
        
        if self.act is not None:
            if self.length is None:
                import sys
                print("Please specify simulation length")
                print('length = <long / short>')
                sys.exit()

        if kind == "sat" or kind == "chain":
            self.kind = "ChainsT"
        elif kind == "chg" or kind == "charge":
            self.kind = "New2"
        elif kind == 'cg':
            self.kind = kind
        if self.kind == "cg":
            self.path = self.path[1]
        else:
            self.path = self.path[0]
        
    def __get_leaflet__(self):
        return self.leaflet_in
    
    def __get_kind__(self):
        return self.kind
    
    def __get_length__(self):
        return self.length
    
    def __get_act__(self):
        return self.act
    
    def update_act(self,act):
        self.act = act
    
    def __test__(self):
        
        def build_simplified_CGTM(self):
            '''
            This should really be moved to the test case 
            '''
            states = pd.read_csv("%s/simplified_raw.csv"%self.path,index_col=0).T
            states = states.fillna(0)
            states_use = pd.DataFrame()
            for s in states:
                if states[s][250:].sum() == 0:
                    continue
                else:
                    states_use[s] = states[s]
            TM_norm = self.build_TM(states.iloc[:,::self.dt])
            pi_eq, eigs = self.solve_pi_eq(TM_norm)
            print("DPPC:  %f"%(pi_eq * aps.all_possible_states()[:,0]).sum())
            print("DOPC:  %f"%(pi_eq * aps.all_possible_states()[:,1]).sum())
            print("CHOL:  %f\n\n"%(pi_eq * aps.all_possible_states()[:,2]).sum())
            
            hist,edge = np.histogram(states_use,bins=len(aps.all_possible_states()),normed=None,range=(0,len(aps.all_possible_states())))
            print("DPPC:  %f"%(hist/hist.sum() * aps.all_possible_states()[:,0]).sum())
            print("DOPC:  %f"%(hist/hist.sum() * aps.all_possible_states()[:,1]).sum())
            print("CHOL:  %f"%(hist/hist.sum() * aps.all_possible_states()[:,2]).sum())        
            return pi_eq, eigs, TM_norm,hist  
        
        def get_A(P):
            I = np.eye(len(P))
            row_1 = np.ones((len(P)))
            A = np.vstack([P.T-I,row_1])
            return A
        def get_B(A):
            B = np.zeros(np.shape(A)[0])
            B[-1] = 1.0
            return B
        def LinSolve(A,B):
            return np.linalg.solve(A.T.dot(A),A.T.dot(B)) 
        
        def detailed_balance(states):
            Tij = self.build_TM(states)
            Tji = self.build_TM(states.T)
            pi_ij = self.solve_pi_eq(Tij)[0]
            pi_ji = self.solve_pi_eq(Tji)[0]
            print(np.allclose(pi_ij,pi_ji))
            
        
        states = pd.read_csv("%s/states/%s%s.csv"%(self.path,self.leaflet_in,self.kind),index_col=0).T   
        detailed_balance(states)
        TM_norm = self.build_TM(states.iloc[:,::self.dt])
        A = get_A(TM_norm)
        B = get_B(A)
        pi_lin = LinSolve(A, B)
        pi_lin[pi_lin < 10**(-12)] = 0
        pi_lin = pi_lin / pi_lin.sum()
        pi_eig, eigs = self.solve_pi_eq(TM_norm)
        return pi_lin, pi_eig, np.allclose(pi_lin, pi_eig), TM_norm


    def get_A(self,P):
        I = np.eye(len(P))
        row_1 = np.ones((len(P)))
        A = np.vstack([P.T-I.T,row_1])
        return A
    def get_B(self,A):
        B = np.zeros(np.shape(A)[0])
        B[-1] = 1.0
        return B
    def LinSolve(self,A,B):
        return np.linalg.solve(A.T.dot(A),A.T.dot(B)) 

## To get raw/brute force distributions  
    def build_raw(self,iterate_sims=False,iterate_time=False):
        if self.act is None:
            states = pd.read_csv("%s/states/%s%s.csv"%(self.path,self.leaflet_in,self.kind),index_col=0).T   
        else:
            states = pd.read_csv("%s/CG/data/states/%s_%s.csv"%(self.path,self.length,self.act),index_col=0).T   
        
        if iterate_sims == False and iterate_time == True:
            hist = []
            states = states.T
            for i in states:
                hist_tmp,edge = np.histogram(states[i],bins=len(aps.all_possible_states()),normed=True,range=(0,len(aps.all_possible_states())))
                hist.append(hist_tmp)
            hist = np.array(hist)
        elif iterate_time == False and iterate_sims == True:
            hist = []
            for i in states:
                hist_tmp,edge = np.histogram(states[i],bins=len(aps.all_possible_states()),normed=True,range=(0,len(aps.all_possible_states())))
                hist.append(hist_tmp)
            hist = np.array(hist)
        else:
            hist,edge = np.histogram(states[states.columns[:]],bins=len(aps.all_possible_states()),normed=True,range=(0,len(aps.all_possible_states())))
        # plt.bar(edge[:-1],hist/hist.sum())
        # plt.savefig("%s_state_raw_%s.pdf"%(nm,kind))
        # plt.close()
        return hist,edge
    
    def get_frame_state_composition(self,sel_frm=0):
        
        out = []
        
        states = pd.read_csv("%s/states/%s%s.csv"%(self.path,self.leaflet_in,self.kind),index_col=0).T[sel_frm]
        all_states = aps.all_possible_states()
        for s in states:
            out.append(all_states[s])
        #hist,edge = np.histogram(states,bins=len(aps.all_possible_states()),normed=None,range=(0,len(aps.all_possible_states())))
        # plt.bar(edge[:-1],hist/hist.sum())
        # plt.savefig("%s_state_raw_%s.pdf"%(nm,kind))
        # plt.close()
        return np.array(out)

## CGTM and the likes   
    def build_TM(self,states,overide=False,symitrize=False):
        all_states = aps.all_possible_states()
        TM_m = np.zeros((len(all_states),len(all_states)))
        
        # this isn't realy nessisary is it
        # good intial fail safe though.
        if np.ndim(states) == 1 and overide==False:
            import os
            print("Your states matrix has one row. This does not work")
            os.exit()
        else:
            for S in states:
                for si in range(1, len(states[S])):
                    # starts here int(states[S][si-1])
                    # ends here int(states[S][si])
                    TM_m[int(states[S][si-1]),int(states[S][si])] += 1
        TM_norm = np.zeros(np.shape(TM_m))
        norm = TM_m.sum(axis=1)

        if symitrize == True:
            TM_m = (TM_m + TM_m.T)/2#np.maximum(TM_m, TM_m.transpose())

        for i,j in enumerate(TM_m):
            TM_norm[i] = np.divide(j , norm[i], out=np.zeros_like(j),where=norm[i]!=0)
        # Makes matrix symetrick (sp)

        


        # TM_norm = pd.DataFrame(TM_norm)
        # TM_tmp = TM_norm.loc[(TM_norm.sum(axis=1) != 0), (TM_norm.sum(axis=0) != 0)]
        # TM_index = TM_tmp.index.values
        # TM_cols = TM_tmp.columns.values
        # del TM_tmp
        # state_ind = [c for c in TM_cols]
        # state_ind = [c for c in TM_index]
        # state_ind = np.unique(state_ind)
        # TM_norm = TM_norm.iloc[state_ind,state_ind]
        return TM_norm   


    def solve_pi_eq(self, P):
        #We have to transpose so that Markov transitions correspond to right multiplying by a column vector.  np.linalg.eig finds right eigenvectors.
        evals, evecs = sla.eigs(P.T, k = 1, which='LM')
        evecs = np.real(evecs)
        evecs[np.abs(evecs) < 10**(-12)] = 0
        pi_eg = (evecs/evecs.sum()).real
        #pi_eg[pi_eg < 10**(-12)] = 0 
        pi_EG = np.array([e[0] for e in pi_eg])
        pi_eg = pi_EG
        
        
        # A = self.get_A(P)
        # B = self.get_B(A)
        # pi_lin = self.LinSolve(A, B)
        # pi_lin[pi_lin < 1E-10] = 0
        
        # if np.allclose(pi_lin, pi_eg) == False:
        #     print("Your distribution is wrong. Matrix algebra and eigen algebra breaks.")
        #     return None
        
        return pi_eg, np.linalg.eig(P.T)            
    
    def build_CGTM(self,symitrize=False):
        ''' 
        
        Built a testing funciton, use that please
        
        def get_A(P):
            I = np.eye(len(P))
            row_1 = np.ones((len(P)))
            A = np.vstack([P.T-I,row_1])
            return A
        def get_B(A):
            B = np.zeros(np.shape(A)[0])
            B[-1] = 1.0
            return B
        def LinSolve(A,B):
            return np.linalg.solve(A.T.dot(A),A.T.dot(B)) '''
        # This is a place holder only
        states = 0 
        
        if self.act == None:
            states = pd.read_csv("%s/states/%s%s.csv"%(self.path,self.leaflet_in,self.kind),index_col=0).T   
        else:
            if self.act == "act" or self.act == "Active":
                self.update_act("active")
            elif self.act == "in" or self.act == "inact" or self.act=="Inactive":
                self.update_act("inactive")
            states = pd.read_csv("%s/CG/data/states/%s_%s.csv"%(self.path,self.length,self.act),index_col=0).T
        TM_norm = self.build_TM(states.iloc[:,::self.dt],symitrize=symitrize)
        pi_eq, eigs = self.solve_pi_eq(TM_norm)
        #pi_eq = pd.Series(pi_eq,index=TM_norm.index)
        # A = get_A(TM_norm)
        # B = get_B(A)
        # pi_lin = LinSolve(A,B)
        
        
        return pi_eq, eigs, TM_norm
      
## Calc convergence and sig values
    def develop_lag(self, dt_max,step):
        flin = pd.read_csv("%s/states/%s%s.csv"%(self.path,self.leaflet_in,self.kind),index_col=0).T   
        dts = [1,5,10,20,50,100,200,300,400,500,600]#np.arange(1,dt_max)
        pi = []
        eig = []
        for tau in dts:
            pi_eq, eigs = self.pi_eq_lag(flin, tau)
            pi.append(pi_eq), eig.append(np.real(np.sort(eigs[0])[-2]))
        return eig

    def sig(self, ref, states, sim_list):
        sig_pi = []
        for s,sl in zip(states,sim_list):
            sig_pi.append(np.sum(np.sqrt((s - ref)**2)) / (sl - 1))
        return sig_pi
    
    def sigConverge_simulations(self):
        states = 0 
        
        if self.act == None:
            states = pd.read_csv("%s/states/%s%s.csv"%(self.path,self.leaflet_in,self.kind),index_col=0).T   
        else:
            print("\n########################################################\n\n")
            print("This will run, however there are not enough simulations at")
            print("this time for this to be beneficial")
            print("\n\n\n########################################################")

            if self.act == "act" or self.act == "Active":
                self.update_act("active")
            elif self.act == "in" or self.act == "inact" or self.act=="Inactive":
                self.update_act("inactive")
            states = pd.read_csv("%s/CG/data/states/%s_%s.csv"%(self.path,self.length,self.act),index_col=0) 
        TM_norm = self.build_TM(states.iloc[:,::self.dt])
        pi_eq_ref, eigs = self.solve_pi_eq(TM_norm) 
        pi_sig = []
        sim_list = []
        if self.act is None:
            sim_list = [2,4,6,8,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,99]
        else:
            sim_list = states.columns
        for nset in sim_list[2:]:
            TM_norm = self. build_TM(states.iloc[::self.dt,:int(nset)]) #switch these axis for sims
            pi_eq, eigs = self.solve_pi_eq(TM_norm)
            pi_sig.append(pi_eq)
        sig_ = []#self.sig(pi_eq_ref,pi_sig,sim_list)
        return sig_,pi_sig#[1:]

    def sigConverge_time(self,overide=True):
        states = 0 
        states_ref = pd.read_csv("%s/CG/data/states/%s_%s.csv"%(self.path,"short",self.act),index_col=0).T 
        if self.act == None:
            print("\n########################################################\n\n")
            print("This will run, however these simulations are not enough")
            print("for this to be useful")
            print("\n\n\n########################################################")
            states = pd.read_csv("%s/states/%s%s.csv"%(self.path,self.leaflet_in,self.kind),index_col=0).T   
        else:

            if self.act == "act" or self.act == "Active":
                self.update_act("active")
            elif self.act == "in" or self.act == "inact" or self.act=="Inactive":
                self.update_act("inactive")
            states = pd.read_csv("%s/CG/data/states/%s_%s.csv"%(self.path,self.length,self.act),index_col=0) 
        state_shape = np.shape(states)
        pi_eq_ref = np.histogram(states_ref,bins=len(aps.all_possible_states()),normed=True,range=(0,len(aps.all_possible_states())))[0]
        # pi_eq_ref, eigs = self.solve_pi_eq(TM_norm) 
        pi_sig = []
        sim_list = []
        if self.act is None:
            sim_list = np.arange(1,99*600,10)
        else:
            sim_list = np.arange(1,len(states.T),1)
            
        for nset in sim_list:
            TM_lim = self. build_TM(states.iloc[:,::nset])
            pi_eq, eigs = self.solve_pi_eq(TM_lim)
            pi_sig.append(pi_eq)
            
            # if pi_sig[-1] > np.mean(pi_sig):
            #     return TM_lim, pi_sig, nset
        sig_ = self.sig(pi_eq_ref,pi_sig,sim_list)
        return sig_
    
    def sigConverge_time_diff(self,overide=True):
        states = 0 
        states_ref = 0
        # the assumption is that 
        if self.act == None:
            print("\n########################################################\n\n")
            print("This will run, however these simulations are not enough")
            print("for this to be useful")
            print("\n\n\n########################################################")
            states = pd.read_csv("%s/states/%s%s.csv"%(self.path,self.leaflet_in,self.kind),index_col=0).T   
        else:

            if self.act == "act" or self.act == "Active":
                self.update_act("active")
            elif self.act == "in" or self.act == "inact" or self.act=="Inactive":
                self.update_act("inactive")
            states = pd.read_csv("%s/CG/data/states/%s_%s.csv"%(self.path,'short',self.act),index_col=0) 
            states_ref = pd.read_csv("%s/CG/data/states/%s_%s.csv"%(self.path,'long',self.act),index_col=0) 
            
        #state_shape = np.shape(states_ref)
        TM_norm = self.build_TM(states_ref.iloc[:,::self.dt],overide)
        pi_eq_ref, eigs = self.solve_pi_eq(TM_norm) 
        pi_sig = []
        sim_list = []
        if self.act is None:
            sim_list = np.arange(1,99*600,10)
        else:
            sim_list = np.arange(1,len(states.T),3)
            
        for nset in sim_list:
            TM_lim = self. build_TM(states.iloc[:,::nset])
            pi_eq, eigs = self.solve_pi_eq(TM_lim)
            pi_sig.append(pi_eq)
            
            # if pi_sig[-1] > np.mean(pi_sig):
            #     return TM_lim, pi_sig, nset
        sig_ = self.sig(pi_eq_ref,pi_sig,sim_list)
        return sig_
    
    def calc_confidence(self,tau=1):
        # should only be used for analysis when you have 
        # both a LONG system and a set of short systems
        # TODO fix this to do that!!
        # tau can be played with, but 1 and 2 work best
        if self.act == None:
            states = pd.read_csv("%s/states/%s%s.csv"%(self.path,self.leaflet_in,self.kind),index_col=0).T   
        else:
            if self.act == "act" or self.act == "Active":
                self.update_act("active")
            elif self.act == "in" or self.act == "inact" or self.act=="Inactive":
                self.update_act("inactive")
            states = pd.read_csv("%s/CG/data/states/%s_%s.csv"%(self.path,"long",self.act),index_col=0).T
        possible_states = aps.all_possible_states()
        # states = []
        # for s in shell:
        #     states.append(self.check_states(s,possible_states))
        # states = np.array(states)
        n1, bins1 = np.histogram(states, bins=len(possible_states),range=(0,len(possible_states)))
        # pi1 = self.build_CGTM()[0]
        n1 = n1/np.sum(n1)
        Z = []
        Zp = []

        # tau = 1
        if self.act != None: 
            if self.act == "act" or self.act == "Active":
                self.update_act("active")
            elif self.act == "in" or self.act == "inact" or self.act=="Inactive":
                self.update_act("inactive")
            states = pd.read_csv("%s/CG/data/states/%s_%s.csv"%(self.path,self.length,self.act),index_col=0).T
        iter_max = len(states)
        while len(states)/iter_max < 4:
            iter_max = iter_max - 1
        # bit of a hack, the states.T lets me run through time
        for i in range(0,int(iter_max),tau):
            # ind = states.index.astype(int)[::i*tau+tau]
            if i%10==0:
                print("Frame %i"%i)
            # n, bins = np.histogram(states.iloc[ind,:], bins=len(possible_states),range=(0,len(possible_states)))
            # n = n/np.sum(n)
            TM = self.build_TM(states.iloc[::i*tau+tau,:])
            pi = self.solve_pi_eq(TM)[0]
            # Pat.append(n)

            phi_tau = np.asarray(self.weighted_avg(pi)) / np.asarray(self.weighted_avg(n1))
            print(self.weighted_avg(pi),phi_tau.sum())
            Z.append(phi_tau)
            Zp.append(phi_tau.sum())
            
        try:
            pd.DataFrame(Z).to_csv("%s/CG/data/%s_ratio_sum.csv"%(self.path,self.act))
            (pd.DataFrame(Zp)/3).to_csv("%s/CG/data/%s_ratio_pi_sum.csv"%(self.path,self.act))
        except:
            pd.DataFrame(Z).to_csv("./%s_ratio_sum.csv"%self.act)
            (pd.DataFrame(Zp)/3).to_csv("./%s_ratio_pi_sum.csv"%self.act)
        #return Z,Zp
    def weighted_avg(self, pi_eq=None):
        if pi_eq is None:
            pi_eq = self.build_CGTM()[0]
        all_states = aps.all_possible_states()
        tmp_wsa = [np.sum((pi_eq*all_states[:,0])),np.sum((pi_eq*all_states[:,1])),np.sum((pi_eq*all_states[:,2]))]
        return tmp_wsa

## Routines to write out analyzed data

    def write_initial_states_distribution(self):
        init_states = self.get_initial_states()
        possible_states = aps.all_possible_states()
        init_hist, init_edge = np.histogram(init_states,bins=len(possible_states),range=(0,len(possible_states)))
        init_hist[init_hist>0] = 1 
        init_hist = (possible_states.T * init_hist).T
        print()
        #init_hist = self.weighted_avg(init_hist)
        if self.act==None:
            pd.DataFrame(init_hist).to_csv("%s/init_raw_%s%s.csv"%(self.path,self.leaflet_in,self.kind))
        else:
            pd.DataFrame(init_hist).to_csv("%s/CG/data/init_raw_%s_%s%s.csv"%(self.path,self.act,self.length,self.kind))

    def write_pi_eq(self,symitrize=False):
        pi_eq = self.build_CGTM()[0]
        if self.act==None:
            pd.DataFrame(pi_eq).to_csv("%s/pi_eq_%s%s.csv"%(self.path,self.leaflet_in,self.kind))
        else:
            pd.DataFrame(pi_eq).to_csv("%s/CG/data/pi_eq_%s_%s%s_tmp.csv"%(self.path,self.act,self.length,self.kind))
        
    def write_pi_raw(self,iterate_time=False, iterate_sims=False):
        pi_raw = self.build_raw(iterate_time=iterate_time,iterate_sims=iterate_sims)[0]
        it_val = ""
        if iterate_time == True:
            it_val = "_time"
        if iterate_sims == True:
            it_val = "_sims"
        if self.act == None:
            pd.DataFrame(pi_raw).to_csv("%s/pi_raw_%s%s%s.csv"%(self.path,self.leaflet_in,self.kind,it_val))
        else:
            pd.DataFrame(pi_raw).to_csv("%s/UDel/CG/data/pi_raw_%s_%s%s%s_upd.csv"%(self.path,self.act,self.length,self.kind,it_val))



# import matplotlib.pyplot as plt
# CGTM_Calculations("",1,"cg","inactive","short").write_pi_eq()
# CGTM_Calculations("",1,"cg","active","short").write_pi_eq()

# test1 = CGTM_Calculations("",1,"cg","inactive","short").build_CGTM(symitrize=True)[0]#.write_pi_raw(iterate_time=True)
# test2 = CGTM_Calculations("",1,"cg","inactive","short").build_raw()#(iterate_time=True)
# # test1 = CGTM_Calculations("SU", 1, "sat",act=None).build_CGTM()[0]
# # test2 = CGTM_Calculations("SU", 1, "sat",act=None).build_raw()

# plt.bar(np.arange(0,len(test1.values)),test1.values)
# plt.bar(np.arange(0,len(test1.values)),test2[0],alpha=.5)
# plt.bar(np.arange(0,len(test1.values)),test1.values - test2[0])
# plt.savefig("yarp.pdf")


# CGTM_Calculations("",1,"cg","inactive","short").write_pi_eq()
# CGTM_Calculations("",1,"cg","active","short").write_pi_eq()

# plt.plot(np.linspace(0,50,len(d1)),d1)
# plt.plot(np.linspace(0,50,len(d2)),d2)
# plt.xlabel("time (ns)")
# plt.ylabel(r"$\sigma$")
# plt.savefig("sig_time.pdf")
# plt.close()


# d3 = CGTM_Calculations("",1,"cg","inactive","short").build_raw()[0]

# d2 = CGTM_Calculations("",1,"cg","inactive","short").sigConverge_simulations()[1]
# # # # d1.write_pi_eq()
# d1.write_pi_raw(iterate_time=True)
# # # d2 = CGTM_Calculations("",1,"cg","active","short")
# # # # d2.calc_confidence()
# # # # d2.write_pi_eq()
# # # d2.write_pi_raw(True,False)
# # # d2.write_pi_raw(False,True)

# # # CGTM_Calculations("",1,"cg","active","short").write_pi_raw()
# d3 = CGTM_Calculations("",1,"cg","inactive","long")
# # # d3.write_pi_eq()
# d3.write_pi_raw(iterate_time=True)
# d4 = CGTM_Calculations("",1,"cg","inactive","short")
# d4.calc_confidence()
# d4.write_pi_eq()
# d4.write_pi_raw()
# CGTM_Calculations("",1,"cg","inactive","short").write_pi_raw()

# pd.DataFrame(CGTM_Calculations("",1,"cg","act","short").sigConverge_simulations()[1]).to_csv("act_short_binned_time_pi_upd.csv")
# pd.DataFrame(CGTM_Calculations("",1,"cg","inact","short").sigConverge_simulations()[1]).to_csv("inact_short_binned_time_pi_upd.csv")

# test1.write_pi_eq()
# test1.write_pi_raw()

# # test1.sigConverge()
# import matplotlib.pyplot as plt
# import CGTM_Plotting as cgp
# d1 = CGTM_Calculations("",1,"cg","inactive","long")
# d1.write_pi_raw()

# d1 = CGTM_Calculations("",1,"cg","active","long")
# d1.write_pi_raw()


# d1 = CGTM_Calculations("",1,"cg","inactive","short")
# d1.write_pi_raw()

# d1 = CGTM_Calculations("",1,"cg","active","short")
# d1.write_pi_raw()


# pd.DataFrame(aps.all_possible_states()).to_csv("all_states.csv")