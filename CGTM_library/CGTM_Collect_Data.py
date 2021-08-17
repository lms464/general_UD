#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 17 16:24:30 2021

@author: liam
"""

import pandas as pd
import numpy as np
import all_possible_states as all_possible_states

class CGTM_Collect_Data:
    def __init__(self,start_sys, end_sys, ref_frame, counting):
        '''
        

        Parameters
        ----------
        start_sys : int
            starting index 0 - N.
        end_sys : int
            ending index M>0 to N.
        ref_frame : list of int
            what frame did I get the data from?.
        counting : str
            charge, or saturation. needs more work

        Returns
        -------
        None.

        '''
        self.start_sys = start_sys
        self.end_sys = end_sys
        self.ref_frame = ref_frame
        self.counting = counting
        self.check_counting_method()

    def check_counting_method(self):
        if self.counting == "sat" or self.counting == "saturation":
            self.build_ternary_charge_states()
        elif self.counting == "chg" or self.counting == "charge":
            self.build_ternary_saturation_states()
        else:
            print(">>> For now, self.counting must be sat or chg")
            return None
        self.cat_states()
            

    def check_states(self, data_frm, possible_states,fl_comp=None):
        '''

        Parameters
        ----------
        data_frm : data frame or array
            Input data.
        possible_states : List of possible states (seee all possible states)
            DESCRIPTION.
        fl_comp : Redundant, optional
            DESCRIPTION. The default is None.

        Returns
        -------
        out : numpy array
            Array of binned state concentrations.

        '''
        #seems fine
        holder = []
        for pi, ps in enumerate(possible_states):
            
            rmsd = np.sqrt((data_frm - ps)**2)
            holder.append(rmsd.sum())
    
        out = np.argmin(holder)
        try:
            fl_comp.write(str([possible_states[out],data_frm]))
        except:
            pass
        return out

    def build_ternary_charge_states(self):
        
        # Initialize resid's assumption there isn't a lot of movement 
        # due to the amount of time that has passed (30ns)
        # up = np.loadtxt("upp.dat")
        # lo = np.loadtxt("low.dat")
        
        possible_states = all_possible_states()
        
        frames = self.ref_frame
        #resids = np.arange(1,1184+1, 1)
    
        
        # resids = [np.loadtxt('CHL1.dat',dtype='int'),np.loadtxt('neutral.dat',dtype='int'),np.loadtxt('anionic.dat',dtype='int')]
        
        # used this to get the shell frm_ is the frame I get the initial state from
        
        for frm_ in frames:
            for i in range(self.start_sys,self.end_sys):
                for j in range(self.start_sys,self.end_sys):
                    
                    # This looks within a single file
                    
                    try: 
                        up = np.loadtxt("/Censere/UDel/resources/test_vor/data/resids/upp2_%i%i%i.dat"%(i,j,frm_))
                        lo = np.loadtxt("/Censere/UDel/resources/test_vor/data/resids/low2_%i%i%i.dat"%(i,j,frm_))
                        resids = [np.loadtxt('/Censere/UDel/resources/test_vor/data/resids/CHL1_%i%i%i.dat'%(i,j,frm_),dtype='int'),
                              np.loadtxt('/Censere/UDel/resources/test_vor/data/resids/neutral_%i%i%i.dat'%(i,j,frm_),dtype='int'),
                              np.loadtxt('/Censere/UDel/resources/test_vor/data/resids/anionic_%i%i%i.dat'%(i,j,frm_),dtype='int')]
                    except:
                        continue
                    
                    print("Calculatings States for system:  %i%i %i..."%(i,j,frm_))
                    states_u = []
                    states_l = []
                    upp_1 = pd.DataFrame(columns=["CHOL","neutral","anionic"], index=np.arange(0,602,1)).fillna(0)
                    low_1 = pd.DataFrame(columns=["CHOL","neutral","anionic"], index=np.arange(0,602,1)).fillna(0)
    
                    # for lip in ["anionic", "neutral", "CHOL"]:
                    try:
                        dat = open("/Censere/UDel/resources/test_vor/data/vor/shells_%i%i%i.update.log"%(i,j,frm_),'r')
                        d_read = dat.readlines()[::3]
                        if len(d_read) < 300:
                            continue
                        dat.close()
                    except:
                        continue
                    #frm is the frame of the "local" simulation
                    # loop through the file and strip data
                    for frm, dr in enumerate(d_read):
                        # an_check = []
                        res = dr.split()[1::2]
                        shell = dr.split()[0::2]
                        for r,s in zip(res,shell):
                            s = int(s)
                            r = int(r)
                            if s==1:
                                #if lip in ["neutral","pl","sm"]:
                                    # resids[0] ONLY contain cholesterol
                                    # 1 for neutral and 2 for anionic
                                    # lo are resids in the upper leaflet (sim is upside down)
                                    # and up is the lwoer leaflet
                                if r in resids[0]:
                                    if r in lo:
                                        low_1["CHOL"][frm] += 1
                                    if r in up:
                                        upp_1["CHOL"][frm] += 1
                                elif r in resids[1]:
                                    if r in lo:
                                        low_1["neutral"][frm] += 1
                                    if r in up:
                                        upp_1["neutral"][frm] +=1
                                #elif lip == "anionic":
                                elif r in resids[2]:
                                    if r in lo:
                                        low_1["anionic"][frm] +=1
                                        #print("lo",r)
                                        # an_check.append(1)
                                    elif r in up:
                                        upp_1["anionic"][frm] +=1
                                        #print("up",r)
                                        # an_check.append(-1)
                        # print()
                        #upp_1.iloc[frm,:] = upp_1.iloc[frm,:].divide(upp_1.iloc[frm,:].sum()) 
                    # fl_comp = open("check_comp.dat","w")
                    for frm_id in upp_1.index:
                        #upp_1.iloc[frm,:].divide(upp_1.iloc[frm,:].sum())
                        states_u.append(self.check_states(upp_1.iloc[frm_id,:].divide(upp_1.iloc[frm_id,:].sum())
                                                   , possible_states))
                        states_l.append(self.check_states(low_1.iloc[frm_id,:].divide(low_1.iloc[frm_id,:].sum())
                                                   , possible_states))
                    
                    # fl_comp.close()
                    pd.DataFrame(states_u).to_csv("/Censere/UDel/resources/test_vor/data/states/StatesU_%i_%i%iNew2.csv"%(frm_,i,j))
                    pd.DataFrame(states_l).to_csv("/Censere/UDel/resources/test_vor/data/states/StatesL_%i_%i%iNew2.csv"%(frm_,i,j))    
    
    def build_ternary_saturation_states(self):
        
        def chain_sel(r, resids, ref_list, sat_list):
        	tmp = 0
        	test = resids[r==resids.values]
        	for c in test.columns:
        	    if r in test[c].values:
        	        tmp = c
        	return tmp
        
        # Initialize resid's assumption there isn't a lot of movement 
        # due to the amount of time that has passed (30ns)
        # up = np.loadtxt("upp.dat")
        # lo = np.loadtxt("low.dat")
        
        possible_states = all_possible_states()
        
        lipid_list = ["CHL1",
                           "DPPC", "LSM", "NSM", "OAPE", "OAPS", "PAPC", "PAPS",
                           "PDPE", "PLAO", "PLAS", "PLPC", "PLQS", "POPC", "POPS",
                           "POPE", "PSM", "SAPI", "SAPS", "SOPC" ] 
        saturation ={"CHL1":[0,0],"DPPC":[0,0], "LSM":[0,0], "NSM":[0,0], "OAPE":[1,4], "OAPS":[1,4], "PAPC":[0,4], "PAPS":[0,4],
                           "PDPE":[0,6], "PLAO":[1,4], "PLAS":[0,4], "PLPC":[0,0], "PLQS":[0,3], "POPC":[0,1], "POPS":[0,1],
                           "POPE":[0,1], "PSM":[0,0], "SAPI":[0,4], "SAPS":[0,4], "SOPC":[0,1]}
        
        saturation = pd.DataFrame(saturation.values(),index=saturation.keys()).astype(float)
        unsat = saturation[(saturation>0).all(axis=1)]
        sat = saturation[(saturation==0).all(axis=1)]
        frames = self.ref_frame
        
        _count_ = 0
        
        #resids = np.arange(1,1184+1, 1)
    
        
        # resids = [np.loadtxt('CHL1.dat',dtype='int'),np.loadtxt('neutral.dat',dtype='int'),np.loadtxt('anionic.dat',dtype='int')]
        
        # used this to get the shell frm_ is the frame I get the initial state from
        upper_sat = []
        lower_sat = []
        for frm_ in frames:
            for i in range(self.start_sys,self.end_sys):
                for j in range(self.start_sys,self.end_sys):
                    sat_list = []
                    # This looks within a single file
                    
                    try: 
                        up = np.loadtxt("/Censere/UDel/resources/test_vor/data/resids/upp2_%i%i%i.dat"%(i,j,frm_))
                        lo = np.loadtxt("/Censere/UDel/resources/test_vor/data/resids/low2_%i%i%i.dat"%(i,j,frm_))
                        resids = [np.loadtxt('/Censere/UDel/resources/test_vor/data/resids/%s_%i%i%i.dat'%(lip,i,j,frm_),dtype='int') for lip in lipid_list ]
                        _count_ = _count_ + 1
    
                    except:
                        continue
                    resids = pd.DataFrame(resids,index=lipid_list).fillna(-1).T.astype(int)
                    print("Calculatings States for system:  %i%i %i..."%(i,j,frm_))
                    states_u = []
                    states_l = []
                    upp_1 = pd.DataFrame(columns=["CHOL","sat","unsat"], index=np.arange(0,602,1)).fillna(0)
                    low_1 = pd.DataFrame(columns=["CHOL","sat","unsat"], index=np.arange(0,602,1)).fillna(0)
    
                    # for lip in ["anionic", "neutral", "CHOL"]:
                    try:
                        dat = open("/Censere/UDel/resources/test_vor/data/vor/shells_%i%i%i.update.log"%(i,j,frm_),'r')
                        d_read = dat.readlines()[::3]
                        if len(d_read) < 300:
                            continue
                        dat.close()
                    except:
                        continue
                    #frm is the frame of the "local" simulation
                    # loop through the file and strip data
                    for frm, dr in enumerate(d_read):
                        tmp_list = []
                        tmp_low = []
                        tmp_upp = []
    
                        # an_check = []
                        res = dr.split()[1::2]
                        shell = dr.split()[0::2]
                        for r,s in zip(res,shell):
                            s = int(s)
                            r = int(r)
                            if s==1:
                                #if lip in ["neutral","pl","sm"]:
                                    # resids[0] ONLY contain cholesterol
                                    # 1 for neutral and 2 for anionic
                                    # lo are resids in the upper leaflet (sim is upside down)
                                    # and up is the lwoer leaflet
                                if r in lo:
                                    qed=chain_sel(r,resids,sat,saturation)
                                    tmp_low.append(saturation.T[qed].values[0])
                                    tmp_low.append(saturation.T[qed].values[1])
                                if r in up:
                                    qed=chain_sel(r,resids,sat,saturation)
                                    tmp_upp.append(saturation.T[qed].values[0])
                                    tmp_upp.append(saturation.T[qed].values[1])
                                # qed=chain_sel(r,resids,sat,saturation)
                                # tmp_list.append(saturation.T[qed].values[0])
                                if r in resids["CHL1"]:
                                    if r in lo:
                                        low_1["CHOL"][frm] += 1
                                    if r in up:
                                        upp_1["CHOL"][frm] += 1
                                        
                                elif r in resids[sat.index].values:
                                    
                                    if r in lo:
                                        low_1["sat"][frm] += 1
                                    if r in up:
                                        upp_1["sat"][frm] +=1
                                        
                                        
                                elif r in resids[unsat.index].values:
                                    if r in lo:
                                        low_1["unsat"][frm] +=1
    
                                    elif r in up:
                                        upp_1["unsat"][frm] +=1
                        #sat_list.append(np.histogram(tmp_list,range(0,10))[0])
                        upper_sat.append(tmp_upp)
                        lower_sat.append(tmp_low)
                    # sat_list = np.array(sat_list)
    
                    # np.savetxt("/Censere/UDel/resources/test_vor/data/states/saturation_%i_%i%i.txt"%(frm_,i,j),sat_list)
                    # np.savetxt("/Censere/UDel/resources/test_vor/data/states/Upper_sat_%i_%i%i.txt"%(frm_,i,j),upper_sat)
                    # np.savetxt("/Censere/UDel/resources/test_vor/data/states/Lower_sat_%i_%i%i.txt"%(frm_,i,j),lower_sat)
    
                        #upp_1.iloc[frm,:] = upp_1.iloc[frm,:].divide(upp_1.iloc[frm,:].sum()) 
                    for frm_id in upp_1.index:
                        #upp_1.iloc[frm,:].divide(upp_1.iloc[frm,:].sum())
                        states_u.append(self.check_states(upp_1.iloc[frm_id,:].divide(upp_1.iloc[frm_id,:].sum())
                                                   , possible_states))
                        states_l.append(self.check_states(low_1.iloc[frm_id,:].divide(low_1.iloc[frm_id,:].sum())
                                                   , possible_states))
                    
                    pd.DataFrame(states_u).to_csv("/Censere/UDel/resources/test_vor/data/states/StatesU_%i_%i%iChainsT.csv"%(frm_,i,j))
                    pd.DataFrame(states_l).to_csv("/Censere/UDel/resources/test_vor/data/states/StatesL_%i_%i%iChainsT.csv"%(frm_,i,j))
    
    def cat_states(self):
        # import glob    
        
        if self.counting == "sat" or self.counting == "saturation":
            ending = "ChainsT"
        if self.counting == "chg" or self.counting == "charge":
            ending = "New2"
        
        # up_leaf = glob.glob("/home/liam/UDel/resources/test_vor/data/StatesU*%s"%ending)
        # lo_leaf = glob.glob("/home/liam/UDel/resources/test_vor/data/StatesL*%s"%ending)
        
        SU = pd.DataFrame(index=np.arange(0,602))
        SL = pd.DataFrame(index=np.arange(0,602))
        ind = 0
        flInd = []
        for frm in [0,100,105]:
            for i in range(0,8):
                for j in range(0,8):
                    fl = '/home/liam/UDel/resources/test_vor/data/states/States'
                    try:
                        tmp_u = pd.read_csv("%sU_%i_%i%i%s.csv"%(fl,frm,i,j,ending),index_col=0,header=0)
                    except:
                        continue
                    if tmp_u.values[int(len(tmp_u)*3/4):].sum() == 0:
                        print(frm,i,j)
                        print( tmp_u.values[int(len(tmp_u)*3/4):].sum())
                        continue
                    else:
                        tmp_l =  pd.read_csv("%sL_%i_%i%s%s.csv"%(fl,frm,i,j,ending),index_col=0,header=0)
                        SU[ind] = tmp_u
                        SL[ind] = tmp_l
                        ind += 1
                        flInd.append("%i%i_%i"%(i,j,frm))
        SU.columns, SL.columns = flInd, flInd
        SU.to_csv("data/SU%s.csv"%ending,columns=flInd)
        SL.to_csv("data/SL%s.csv"%ending,columns=flInd)
        S = pd.concat([SU,SL],axis=1)
        S.to_csv("data/States%s.csv"%ending)
      
        
test = CGTM_Collect_Data(0,1,[0],'charge')
#test.build_ternary_charge_states()
