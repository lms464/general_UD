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
            
            print(np.allclose(pi_ij@Tij, pi_ji@Tji))
        
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

    def update_act(self,act):
        self.act = act
        
    def build_raw(self):
        if self.act is None:
            states = pd.read_csv("%s/states/%s%s.csv"%(self.path,self.leaflet_in,self.kind),index_col=0).T   
        else:
            states = pd.read_csv("%s/CG/data/states/%s_%s.csv"%(self.path,self.length,self.act),index_col=0).T   

        hist,edge = np.histogram(states,bins=len(aps.all_possible_states()),normed=None,range=(0,len(aps.all_possible_states())))
        # plt.bar(edge[:-1],hist/hist.sum())
        # plt.savefig("%s_state_raw_%s.pdf"%(nm,kind))
        # plt.close()
        return hist/np.sum(hist),edge

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

    
    def build_TM(self,states,overide=False):
        all_states = aps.all_possible_states()
        TM_m = np.zeros((len(all_states),len(all_states)))
        # norm_t = []
        
        # this isn't realy nessisary is it
        # good intial fail safe though.
        if np.ndim(states) == 1 and overide==False:
            import os
            print("Your states matrix has one row. This does not work")
            os.exit()
            # for si in range(1, len(states)):
            #     TM_m[int(states[si]),int(states[si-1])] += 1 
        else:
            for S in states:
                for si in range(1, len(states[S])):
                    TM_m[int(states[S][si]),int(states[S][si-1])] += 1 
        norm = TM_m.sum(axis=1)
        TM_norm = np.zeros(np.shape(TM_m))
        for i,j in enumerate(TM_m):
            TM_norm[i] = j / norm[i]
        # TM_sym = 1/2 * (TM_norm + TM_norm.T)

        #TM_norm = np.divide(TM_m, norm)
        TM_norm = np.nan_to_num(TM_norm)
        # TM_test = np.nan_to_num(np.divide(TM_m, norm))
        return TM_norm    


    def solve_pi_eq(self, P):
        evals, evecs = sla.eigs(P.T, k = 1, which='LM')
        evecs = np.real(evecs)
        evecs[np.abs(evecs) < 10**(-12)] = 0
        pi_eg = (evecs/evecs.sum()).real
        #pi_eg[pi_eg < 10**(-12)] = 0 
        pi_EG = np.array([e[0] for e in pi_eg])
        pi_eg = pi_EG
        return pi_eg, np.linalg.eig(P.T)            
    
    def build_CGTM(self):
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
            states = pd.read_csv("%s/CG/data/states/%s_%s.csv"%(self.path,self.length,self.act),index_col=0)
        TM_norm = self.build_TM(states.iloc[:,::self.dt])
        pi_eq, eigs = self.solve_pi_eq(TM_norm)
        # A = get_A(TM_norm)
        # B = get_B(A)
        # pi_lin = LinSolve(A,B)
        
        
        return pi_eq, eigs, TM_norm
      

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
            sim_list = states.index
        for nset in sim_list[2:]:
            TM_norm = self. build_TM(states.iloc[:nset,::self.dt])
            pi_eq, eigs = self.solve_pi_eq(TM_norm)
            pi_sig.append(pi_eq)
        sig_ = self.sig(pi_eq_ref,pi_sig,sim_list)
        return sig_,pi_sig#[1:]

    def sigConverge_time(self,overide=True):
        states = 0 
        
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
        TM_norm = self.build_TM(states.iloc[:,::self.dt],overide)
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

    def get_initial_states(self):
        if self.act != None: 
            if self.act == "act" or self.act == "Active":
                self.update_act("active")
            elif self.act == "in" or self.act == "inact" or self.act=="Inactive":
                self.update_act("inactive")
        states = pd.read_csv("%s/CG/data/states/%s_%s.csv"%(self.path,self.length,self.act),index_col=0).T
        return states.iloc[0,:]

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

    def write_pi_eq(self):
        pi_eq = self.build_CGTM()[0]
        if self.act==None:
            pd.DataFrame(pi_eq).to_csv("%s/pi_eq_%s%s.csv"%(self.path,self.leaflet_in,self.kind))
        else:
            pd.DataFrame(pi_eq).to_csv("%s/CG/data/pi_eq_%s_%s%s.csv"%(self.path,self.act,self.length,self.kind))
        
    def write_pi_raw(self):
        pi_raw = self.build_raw()[0]
        if self.act == None:
            pd.DataFrame(pi_raw).to_csv("%s/pi_raw_%s%s.csv"%(self.path,self.leaflet_in,self.kind))
        else:
            pd.DataFrame(pi_raw).to_csv("%s/CG/data/pi_raw_%s_%s%s.csv"%(self.path,self.act,self.length,self.kind))




# d1 = CGTM_Calculations("",1,"cg","active","long")
# d11 = d1.weighted_avg(d1.build_raw()[0])
# d2 = CGTM_Calculations("",1,"cg","active","short")
# d21 = d2.weighted_avg(d2.build_raw()[0])
# d22 = d2.weighted_avg()
# # CGTM_Calculations("",1,"cg","active","short").write_pi_raw()
# d3 = CGTM_Calculations("",1,"cg","inactive","long")
# d31 = d3.weighted_avg(d3.build_raw()[0])
# d4 = CGTM_Calculations("",1,"cg","inactive","short")
# d41 = d4.weighted_avg(d4.build_raw()[0])
# d42 = d4.weighted_avg()

# CGTM_Calculations("",1,"cg","inactive","short").write_pi_raw()

# pd.DataFrame(CGTM_Calculations("",1,"cg","act","short").sigConverge_simulations()[1]).to_csv("act_short_binned_pi.csv")
# pd.DataFrame(CGTM_Calculations("",1,"cg","inact","short").sigConverge_simulations()[1]).to_csv("inact_short_binned_pi.csv")

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