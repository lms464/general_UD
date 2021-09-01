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
        states = pd.read_csv("%s/states/%s%s.csv"%(self.path,self.leaflet_in,self.kind),index_col=0).T   
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

    
    def build_TM(self,states):
        all_states = aps.all_possible_states()
        TM_m = np.zeros((len(all_states),len(all_states)))
        # norm_t = []
        
        if np.ndim(states) == 1:
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
            states = pd.read_csv("%s/CG/data/states/shot_%s.csv"%(self.path,self.act),index_col=0)
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

    def sig(self, ref, states,sim_list):
        sig_pi = []
        for s,sl in zip(states,sim_list):
            sig_pi.append(np.sum(np.sqrt((s - ref)**2)) / (sl - 1))
        return sig_pi
    
    def sigConverge(self):
        states = pd.read_csv("%s/states/%s%s.csv"%(self.path,self.leaflet_in,self.kind),index_col=0).T   
        TM_norm = self.build_TM(states.iloc[:,::self.dt])
        pi_eq_ref, eigs = self.solve_pi_eq(TM_norm) 
        pi_sig = []
        sim_list = [2,4,6,8,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,99]
        for nset in sim_list:
            TM_norm = self. build_TM(states.iloc[:nset,::self.dt])
            pi_eq, eigs = self.solve_pi_eq(TM_norm)
            pi_sig.append(pi_eq)
        sig_ = self.sig(pi_eq_ref,pi_sig,sim_list)
        return sig_
    
    # def build_CGTM_series(self):
    # THIS SHOULDN'T BE
    #
    #     states = pd.read_csv("%s/states/%s%s.csv"%(self.path,self.leaflet_in,self.kind),index_col=0).T   
    #     pi_eq = []
    #     eigs = []
    #     for state in states.T:
    #         TM_norm = self.build_TM(states.T[state])
    #         pi_tmp, eigs_tmp = self.solve_pi_eq(TM_norm)
    #         pi_eq.append(pi_tmp)
    #         eigs.append(eigs_tmp)
    #     return np.asarray(pi_eq),eigs
    def calc_confidence(self,tau=1):
        # tau can be played with, but 1 and 2 work best
        
        states = pd.read_csv("%s/states/%s%s.csv"%(self.path,self.leaflet_in,self.kind),index_col=0).T   
        possible_states = aps.all_possible_states()
        # states = []
        # for s in shell:
        #     states.append(self.check_states(s,possible_states))
        # states = np.array(states)
        n1, bins1 = np.histogram(states, bins=len(possible_states),range=(0,len(possible_states)))
        pi1 = self.build_CGTM()[0]
        n1 = n1/np.sum(n1)
        Z = []
        Zp = []
        Pat = []
        Rat = []
        tau = 1
        # bit of a hack, the states.T lets me run through time
        for i in range(0,len(states.T),tau):
            print(i)
            n, bins = np.histogram(states.T[:i*tau+tau], bins=len(possible_states),range=(0,len(possible_states)))
            n = n/np.sum(n)
            TM = self.build_TM(states.T[:i*tau+tau].T)
            pi = self.solve_pi_eq(TM)[0]
            Pat.append(n)
            ratio= n/n1
            ratio=np.nan_to_num(ratio,0)
            Rat.append(ratio)
            ratio_sum = np.sum(ratio)/len(n1[n1>0])
            ratio_pi = pi/pi1
            ratio_pi_sum = (np.sum(ratio_pi)/len(pi1[pi1>0]))
            Z.append(ratio_sum)
            Zp.append(ratio_pi_sum)
        return Z,Zp
    def series_weighted_avg(self):
        
        wsa = []
        pi_eq, eigs = self.build_CGTM_series()
        all_states = aps.all_possible_states()
        for pi in pi_eq:
            tmp_wsa = [np.sum((pi*all_states[:,0])),np.sum((pi*all_states[:,1])),np.sum((pi*all_states[:,2]))]
            wsa.append(tmp_wsa)
        return np.asarray(wsa)

    def write_pi_eq(self):
        pi_eq = self.build_CGTM()[0]
        if self.act==None:
            pd.DataFrame(pi_eq).to_csv("%s/pi_eq_%s%s.csv"%(self.path,self.leaflet_in,self.kind))
        else:
            pd.DataFrame(pi_eq).to_csv("%s/CG/data/pi_eq_%s_%s%s.csv"%(self.path,self.act,self.leaflet_in,self.kind))
        
    def write_pi_raw(self):
        pi_raw = self.build_raw()[0]
        if self.act == None:
            pd.DataFrame(pi_raw).to_csv("%s/pi_raw_%s%s.csv"%(self.path,self.leaflet_in,self.kind))
        else:
            pd.DataFrame(pi_raw).to_csv("%s/CG/data/pi_raw_%s_%s%s.csv"%(self.path,self.act,self.leaflet_in,self.kind))

test1 = CGTM_Calculations("",1,"cg","act")
test1.write_pi_eq()
# test1.sigConverge()
# import CGTM_Plotting as cgp
# cgp.plot_sigConverge( CGTM_Calculations("SU",1,"charge").sigConverge(),CGTM_Calculations("SL",1,"charge").sigConverge(),test1.__get_kind__())
# zed = test1.calc_confidence()
#pi_lin, pi_eig, comp, TM = test1.__test__()
# test1.write_pi_eq()
# test1.write_pi_raw()

# test1 = CGTM_Calculations("SL",1,"charge")
# test1.build_CGTM()
# test1.write_pi_eq()
# test1.write_pi_raw()

# data1 = test1.series_weighted_avg()

# import ternary
# import matplotlib.pyplot as plt


# figure, tax = ternary.figure(scale=1)
# figure.set_size_inches(10, 10)
# tax.scatter(data1)
# # tax.scatter(data2)
# tax.boundary(linewidth=2.0)
# tax.gridlines(multiple=.1, color="blue")
# tax.ticks(axis='lbr', linewidth=.5, multiple=1)
# tax.clear_matplotlib_ticks()
# tax.get_axes().axis('off')
# # tax.savefig("sat.pdf")
# # tax.close()
