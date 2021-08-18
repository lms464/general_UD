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

class CGTM_Calculations:
    def __init__(self, leaflet_in, dt, kind):
        
        assert dt > 0, "dt MUST be greater than 0"
        
        self.leaflet_in = leaflet_in
        self.dt = dt
        self.kind = None
        
        if kind == "sat" or kind == "chain":
            self.kind = "CH"
        elif kind == "chg" or kind == "charge":
            self.kind = "NEW2"
        
    def __get_leaflet__(self):
        return self.leaflet_in
    
    def __get_kind__(self):
        return self.kind
    
    
    def build_TM(self,states):
        all_states = aps.all_possible_states()
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


    def solve_pi_eq(self, P):
        evals, evecs = sla.eigs(P.T, k = 1, which='LM')
        evecs = np.real(evecs)
        pi_eg = (evecs/evecs.sum()).real
        pi_eg[pi_eg < 10**(-7)] = 0 #TODO wait... why?
        pi_EG = np.array([e[0] for e in pi_eg])
        pi_eg = pi_EG
        return pi_eg, np.linalg.eig(P.T)            
    
    def build_CGTM(self):
    
        states = pd.read_csv("/Censere/UDel/resources/test_vor/data/%s%s.csv"%(self.leaflet_in,self.kind),index_col=0).T   
        TM_norm = self.build_TM(states.iloc[:,::self.dt])
        pi_eq, eigs = self.solve_pi_eq(TM_norm)
        return pi_eq, eigs, TM_norm
        
    def develop_lag(self, dt_max,step):
        flin = pd.read_csv("./data/%s%s.csv"%(self.leaflet_in,self.kind),index_col=0).T   
        dts = [1,5,10,20,50,100,200,300,400,500,600]#np.arange(1,dt_max)
        pi = []
        eig = []
        for tau in dts:
            pi_eq, eigs = self.pi_eq_lag(flin, tau)
            pi.append(pi_eq), eig.append(np.real(np.sort(eigs[0])[-2]))
        return eig

    def sig(self, ref, states,sim_list):
        sig_pi = []
        for s,sl in zip(states,sim_list):
            sig_pi.append(np.sum(np.sqrt((s - ref)**2)) / (sl - 1))
        return sig_pi
    
    def sigConverge(self):
        states = pd.read_csv("/Censere/UDel/resources/test_vor/data/%s%s.csv"%(self.leaflet_in,self.kind),index_col=0).T   
        TM_norm = self.build_TM(states.iloc[0,::self.dt])
        pi_eq_ref, eigs = self.solve_pi_eq(TM_norm) 
        pi_sig = []
        sim_list = [2,4,6,8,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,99]
        for nset in sim_list:
            TM_norm = self. build_TM(states.iloc[:nset,::self.dt])
            pi_eq, eigs = self.solve_pi_eq(TM_norm)
            pi_sig.append(pi_eq)
        sig_ = self.sig(pi_eq_ref,pi_sig,sim_list)
        return sig_
    
    def build_CGTM_series(self):
        states = pd.read_csv("/Censere/UDel/resources/test_vor/data/%s%s.csv"%(self.leaflet_in,self.kind),index_col=0).T   
        pi_eq = []
        eigs = []
        for state in states.T:
            TM_norm = self.build_TM(states.T[state])
            pi_tmp, eigs_tmp = self.solve_pi_eq(TM_norm)
            pi_eq.append(pi_tmp)
            eigs.append(eigs_tmp)
        return np.asarray(pi_eq),eigs

    def series_weighted_avg(self):
        
        wsa = []
        
        pi_eq, eigs = self.build_CGTM_series()
        all_states = aps.all_possible_states()
        for pi in pi_eq:
            tmp_wsa = [np.sum((pi*all_states[:,0])),np.sum((pi*all_states[:,1])),np.sum((pi*all_states[:,2]))]
            wsa.append(tmp_wsa)
        return np.asarray(wsa)

test = CGTM_Calculations("SU",1,"sat")
wsa = test.series_weighted_avg()