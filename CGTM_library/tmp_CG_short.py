import pandas as pd
import numpy as np
import all_possible_states as aps
import choose_path as chp
import logging

class CGTM_Collect_Data:
    def __init__(self,start_sys, end_sys, ref_frame, counting, act=None, length=None,titrate=None,chp_inpt=None):
        '''
        This is really just a complex path chooser
        and state organizer. That's all it does based
        on initial inputs.
        
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

        titrate is reduntant

        Returns
        -------
        None.

        '''


        self.start_sys = start_sys
        self.end_sys = end_sys
        self.ref_frame = ref_frame
        self.counting = counting
        self.path = chp.choose_path(chp_inpt)
        self.act = act
        self.length = length
        # choosing location of directory
        if self.counting == "cg":
            self.path = self.path[1]
        else:
            self.path = self.path[0]
        # if titrate is not None:
        #     self.path = '/home/sharplm/CG/data/shells'
        
        if self.act is not None:
            if self.length is None:
                import sys
                print("Please specify simulation length")
                print('length = <long / short>')
                sys.exit()

    def update_act(self,act):
        self.act = act

    def check_counting_method(self):
        '''Redundant confirm you can delete'''
        if self.counting == "sat" or self.counting == "saturation":
            self.build_ternary_charge_states()
        elif self.counting == "chg" or self.counting == "charge":
            self.build_ternary_saturation_states()
        elif self.counting == "cg":
            print("Implement parse_hop!")
        else:
            print(">>> For now, self.counting must be sat, chg, cg")
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
        # holds rmsd data
        holder = []
        # iterates over possible states and index
        for pi, ps in enumerate(possible_states):
            
            rmsd = np.sqrt((data_frm - ps)**2)
            holder.append(rmsd.sum())
    
        out = np.argmin(holder)
        try:
            # old and can be removed
            fl_comp.write(str([possible_states[out],data_frm]))
        except:
            pass
        return out    

    def build_cg_short_states(self):
        
        if self.length == "long":
            print("Cannot find long simulation data using")
            print("short argument.")
            import sys
            sys.exit()
        
        if self.act == "act" or self.act == "Active":
            self.update_act("active")
        elif self.act == "in" or self.act == "inact" or self.act=="Inactive":
            self.update_act("inactive")
        
        # loads in state space
        possible_states = aps.all_possible_states()

        # holder for states onse determined
        # full_states = []
        
        states_u = []
        states_l = []
        #iterates over each file
        for i in range(self.start_sys,self.end_sys+1):
            print("Running System %i..."%i)
            ##### Load Resid's ##########
            try: 
                # up = np.loadtxt("%s/CG/data/resids/upp2_short_%s_%i.dat"%(self.path,self.act,i))
                # lo = np.loadtxt("%s/CG/data/resids/low2_short_%s_%i.dat"%(self.path,self.act,i))
                
                up_list = []
                lo_list = []
                
                leaflet =  np.loadtxt("%s/CG/data/resids/leaflet_%s%i.txt"%(self.path,self.act,i),dtype=int)
                
                for li in leaflet:
                    lo = []
                    up = []
                    for lr, lf in zip(li[:-1:2],li[1::2]):
                        if lf == 1:
                            lo.append(lr)
                        if lf == 0:
                            up.append(lr)
                    up_list.append(up)
                    lo_list.append(lo)
                    
                up_list = np.array(up_list)
                lo_list = np.array(lo_list)
                
                
                resids = [np.loadtxt('%s/CG/data/resids/DPPC_short_%s_%i.dat'%(self.path,self.act,i),dtype='int'),
                        np.loadtxt('%s/CG/data/resids/CHOL_short_%s_%i.dat'%(self.path,self.act,i),dtype='int'),
                        np.loadtxt('%s/CG/data/resids/DOPC_short_%s_%i.dat'%(self.path,self.act,i),dtype='int')]
            except:
                print("failed at resid load")
            #############################



            fl = open ("%s/CG/data/shells/%s_shell%i.txt"%(self.path,self.act,i))
            lines = fl.readlines()#.split()
            fl.close()

            ###### leaflet splitter #####

            upp_1 = pd.DataFrame(columns=["DPPC","CHOL","DOPC"], index=np.arange(0,int(len(lines)/3),1)).fillna(0) #what will index be?
            low_1 = pd.DataFrame(columns=["DPPC","CHOL","DOPC"], index=np.arange(0,int(len(lines)/3),1)).fillna(0) #double check lipid order
            ##############################

            #### split by leaflet ########

            n_res = []
            for frm, (dr,lo,up) in enumerate(zip(lines[::3],lo_list,up_list)):
                tmp_res = 0
                res = dr.split()[1::2]
                shell = dr.split()[0::2]
                for r,s in zip(res,shell):
                    s = int(s)
                    r = int(r)
                    if s==1:
                        tmp_res = tmp_res + 1
                        #if lip in ["neutral","pl","sm"]:
                            # resids[0] ONLY contain DPPC
                            # 1 for CHOL and 2 for DOPC
                            # lo are resids in the upper leaflet (sim is upside down)
                            # and up is the lwoer leaflet
                        if r in resids[0]:
                            if r in lo:
                                low_1["DPPC"][frm] += 1
                            elif r in up:
                                upp_1["DPPC"][frm] += 1
                        elif r in resids[1]:
                            if r in lo:
                                low_1["CHOL"][frm] += 1
                            elif r in up:
                                upp_1["CHOL"][frm] += 1
                        elif r in resids[2]:
                            if r in lo:
                                low_1["DOPC"][frm] += 1
                            elif r in up:
                                upp_1["DOPC"][frm] += 1
                n_res.append(tmp_res)
            tmp_u = []
            tmp_l = []
            ###############################
            for frm_id in upp_1.index:
                        #upp_1.iloc[frm,:].divide(upp_1.iloc[frm,:].sum())
                # TODO this term here isn't behaving as I thought it should....
                tmp_u.append(self.check_states(upp_1.iloc[frm_id,:].divide(upp_1.iloc[frm_id,:].sum())
                                        , possible_states))
                tmp_l.append(self.check_states(low_1.iloc[frm_id,:].divide(low_1.iloc[frm_id,:].sum())
                                        , possible_states))
            states_u.append(tmp_u)
            states_l.append(tmp_l)
            # #lines = [int(l) for l in lines]
            # shell = []
            # # short, leaflet doesn't mater, gets first shell (every 3 lines form 1)
            # for l in lines[1::3]:
            #     shell.append([int(l.split()[-3]),int(l.split()[-2]),int(l.split()[-1])]) 
                
            #     # Sanitiy check. values should be below 45, allowing up to 55
            #     if np.sum(int(l.split()[-3])+int(l.split()[-2])+int(l.split()[-1])) > 55:
            #         print("Unrealistic shell count. Confirm first shell is being analyzed")
            #         return None
            # # converts to numpy array for math
            # shell_arr = np.array(shell)
            # shell = np.array([shell_arr[:,0] / shell_arr.sum(axis=1),shell_arr[:,1] / shell_arr.sum(axis=1),shell_arr[:,2] / shell_arr.sum(axis=1)])
            # shell = shell.T
            # states = []
            
            # for s in shell:
            #     # determines all states for a file
            # 	states.append(self.check_states(s,possible_states))
            # full_states.append(states)
        pd.DataFrame(states_l).to_csv("%s/CG/data/states/low_%s_%s_tmpNEW.csv"%(self.path, self.length, self.act))
        pd.DataFrame(states_u).to_csv("%s/CG/data/states/upp_%s_%s_tmpNEW.csv"%(self.path, self.length, self.act))
        # hist_low = np.histogram(states_l,bins=231,range=(0,231))[0]
        # hist_low = hist_low/hist_low.sum()
        # pd.DataFrame(hist_low).to_csv("/home/sharplm/CG/data/pi_raw_low_inactive_short_cg.csv")

    def build_cg_short_titrate(self):
        import glob
        logging.basicConfig(filename='extended2.log',  level=logging.DEBUG)
        if self.length == "long":
            print("Cannot find long simulation data using")
            print("short argument.")
            import sys
            sys.exit()
        
        if self.act == "act" or self.act == "Active":
            self.update_act("active")
        elif self.act == "in" or self.act == "inact" or self.act=="Inactive":
            self.update_act("inactive")
        
        # loads in state space
        possible_states = aps.all_possible_states()
        
        # holder for states onse determined
        full_states = []
        #iterates over each file
        states_u, states_l = [],[]
        
        for i in range(1,len(glob.glob("%s/CG/data/shells/titrate*.35*"%self.path))+1):
            for j in [".35",".4",".45",".5",".55",".6",".65",".7",".75"]:#glob.glob("%s/titration1-10/inactive%i/DOPC*"%(self.path,i)):
                print("Running System %s-%s..."%(i,j))
                try: 
                    # up = np.loadtxt("%s/CG/data/resids/upp2_titrate_%i_%s.dat"%(self.path,i,j))
                    # lo = np.loadtxt("%s/CG/data/resids/low2_titrate_%i_%s.dat"%(self.path,i,j))
                    
                    up_list = []
                    lo_list = []
                    
                    leaflet =  np.loadtxt("%s/CG/data/resids/leaflet_%s%i%s_titrate.txt"%(self.path,self.act,i,j),dtype=int)
                    
                    for li in leaflet:
                        lo = []
                        up = []
                        for lr, lf in zip(li[:-1:2],li[1::2]):
                            if lf == 1:
                                lo.append(lr)
                            if lf == 0:
                                up.append(lr)
                        up_list.append(up)
                        lo_list.append(lo)
                        
                    up_list = np.array(up_list)
                    lo_list = np.array(lo_list)
                    
                    resids = [np.loadtxt('%s/CG/data/resids/DPPC_titrate_%i_%s.dat'%(self.path,i,j),dtype='int'),
                            np.loadtxt('%s/CG/data/resids/CHOL_titrate_%i_%s.dat'%(self.path,i,j),dtype='int'),
                            np.loadtxt('%s/CG/data/resids/DOPC_titrate_%i_%s.dat'%(self.path,i,j),dtype='int')]
                except:
                    print("failed at resid load")
                # try:
                upp_1 = pd.DataFrame(columns=["DPPC","CHOL","DOPC"], index=np.arange(0,len(leaflet),1)).fillna(0) #what will index be?
                low_1 = pd.DataFrame(columns=["DPPC","CHOL","DOPC"], index=np.arange(0,len(leaflet),1)).fillna(0) #double check lipid order
                fl = open("%s/CG/data/shells/titrate_inactive_shell%i%s.txt"%(self.path,i,j))
                lines = fl.readlines()#.split()
                fl.close()
                #lines = [int(l) for l in lines]
                shell = []
                
                upp2 = upp_1.copy()
                low2 = low_1.copy()
                
                n_res = []
                for frm, (dr,v_shell) in enumerate(zip(lines[::3],lines[1::3])):
                    v_shell = [int(v) for v in v_shell.split(":")[1].split("\n")[0].split(" ")[-3:]]
                    tmp_res = 0
                    res = dr.split()[1::2]
                    shell = dr.split()[0::2]
                    for r,s in zip(res,shell):
                        s = int(s)
                        r = int(r)
                        if s==1:
                            tmp_res = tmp_res + 1
                            #if lip in ["neutral","pl","sm"]:
                                # resids[0] ONLY contain DPPC
                                # 1 for CHOL and 2 for DOPC
                                # lo are resids in the upper leaflet (sim is upside down)
                                # and up is the lwoer leaflet
                            if r in resids[0]:
                                if r in lo:
                                    low_1["DPPC"][frm] += 1
                                elif r in up:
                                    upp_1["DPPC"][frm] += 1
                                else:
                                    logging.debug("There is no DPPC %i"%r)
                            elif r in resids[1]:
                                if r in lo:
                                    low_1["CHOL"][frm] += 1
                                elif r in up:
                                    upp_1["CHOL"][frm] += 1
                                else:
                                    logging.debug("There is no CHOL %i"%r)

                            elif r in resids[2]:
                                if r in lo:
                                    low_1["DOPC"][frm] += 1
                                elif r in up:
                                    upp_1["DOPC"][frm] += 1
                                else:
                                    logging.debug("There is no DOPC %i"%r)

                                    
                            if (low_1 - low2).sum().sum()>1 or (upp_1 - upp2).sum().sum()>1:                                   
                                print(frm)
                            upp2 = upp_1.copy()
                            low2 = low_1.copy()
                    n_res.append(tmp_res)
                    if (low_1.T[frm]+upp_1.T[frm])["DPPC"] != v_shell[0] or (low_1.T[frm]+upp_1.T[frm])["DOPC"] != v_shell[2]:
                        logging.debug("\nSys: %s%s\nFrame %i\n\n"%(i,j,frm))
                        logging.debug(low_1.T[frm]+upp_1.T[frm].T.values)
                        logging.debug(v_shell)
                    else:
                        logging.debug("all good")

                tmp_u = []
                tmp_l = []
                for frm_id in upp_1.index:
                    print(frm_id)
                    #upp_1.iloc[frm,:].divide(upp_1.iloc[frm,:].sum())
                    # TODO this term here isn't behaving as I thought it should....
                    tmp_u.append(self.check_states(upp_1.iloc[frm_id,:].divide(upp_1.iloc[frm_id,:].sum())
                                            , possible_states))
                    tmp_l.append(self.check_states(low_1.iloc[frm_id,:].divide(low_1.iloc[frm_id,:].sum())
                                            , possible_states))
                states_u.append(tmp_u)
                states_l.append(tmp_l)
                # short, leaflet doesn't mater, gets first shell (every 3 lines form 1)
    #             for l in lines[1::3]:
    #                 shell.append([int(l.split()[-3]),int(l.split()[-2]),int(l.split()[-1])]) 
    #                 # Sanitiy check. values should be below 45, allowing up to 55
    #                 if np.sum(int(l.split()[-3])+int(l.split()[-2])+int(l.split()[-1])) > 55:
    #                     print("Unrealistic shell count. Confirm first shell is being analyzed")
    #                     # return None
    #             # converts to numpy array for math
    #             shell_arr = np.array(shell)
    #             shell = np.array([shell_arr[:,0] / shell_arr.sum(axis=1),shell_arr[:,1] / shell_arr.sum(axis=1),shell_arr[:,2] / shell_arr.sum(axis=1)])
    #             shell = shell.T
    #             states = []
                
    #             for s in shell:
    #                 # determines all states for a file
    #             	states.append(self.check_states(s,possible_states))
    #             full_states.append(states)
                # except:
                #     print("Empty: %s-%s"%(i,j))
       
        ###############################
        pd.DataFrame(states_l).to_csv("%s/CG/data/states/low_%s_%s_titrate_tmpN.csv"%(self.path, self.length, self.act))
        pd.DataFrame(states_u).to_csv("%s/CG/data/states/upp_%s_%s_titrate_tmpN.csv"%(self.path, self.length, self.act))
        # TODO This needs to be moved into the loop, or lipids just keep getting added to the original low_1 variable

        # states_u.append(tmp_u)
        # states_l.append(tmp_l)
    # pd.DataFrame(full_states).to_csv("%s/%s_%s_titrate.csv"%("/home/sharplm/CG/data/states", self.length, self.act))


    def build_cg_long_states(self):
        
        
        if self.act == "act" or self.act == "active" or self.act=="Active":
            self.update_act("Active")
        elif self.act == "in" or self.act == "inact" or self.act=="Inactive":
            self.update_act("Inactive")
        
        # loads in state space
        possible_states = aps.all_possible_states()

        # holder for states onse determined
        # full_states = []
        
        states_u = []
        states_l = []
        #iterates over each file
        for i in range(self.start_sys,self.end_sys+1):
            print("Running System %i..."%i)
            ##### Load Resid's ##########
            try: 
                up = np.loadtxt("%s/CG/data/resids/upp2_long_%s_%i.dat"%(self.path,self.act,i))
                lo = np.loadtxt("%s/CG/data/resids/low2_long_%s_%i.dat"%(self.path,self.act,i))
                
                resids = [np.loadtxt('%s/CG/data/resids/DPPC_long_%s_%i.dat'%(self.path,self.act,i),dtype='int'),
                        np.loadtxt('%s/CG/data/resids/CHOL_long_%s_%i.dat'%(self.path,self.act,i),dtype='int'),
                        np.loadtxt('%s/CG/data/resids/DOPC_long_%s_%i.dat'%(self.path,self.act,i),dtype='int')]
            except:
                print("failed at resid load")
            #############################


            fl = open ("%s/%s%i/borders.txt"%(self.path,self.act,i))
            lines = fl.readlines()#.split()
            fl.close()
            
            ###### leaflet splitter #####

            upp_1 = pd.DataFrame(columns=["DPPC","CHOL","DOPC"], index=np.arange(0,int(len(lines)/3),1)).fillna(0) #what will index be?
            low_1 = pd.DataFrame(columns=["DPPC","CHOL","DOPC"], index=np.arange(0,int(len(lines)/3),1)).fillna(0) #double check lipid order
            ##############################

            n_res = []
            for frm, dr in enumerate(lines[::3]):
                tmp_res = 0
                res = dr.split()[1::2]
                shell = dr.split()[0::2]
                for r,s in zip(res,shell):
                    s = int(s)
                    r = int(r)
                    if s==1:
                        tmp_res = tmp_res + 1
                        #if lip in ["neutral","pl","sm"]:
                            # resids[0] ONLY contain DPPC
                            # 1 for CHOL and 2 for DOPC
                            # lo are resids in the upper leaflet (sim is upside down)
                            # and up is the lwoer leaflet
                        if r in resids[0]:
                            if r in lo:
                                low_1["DPPC"][frm] += 1
                            elif r in up:
                                upp_1["DPPC"][frm] += 1
                        elif r in resids[1]:
                            if r in lo:
                                low_1["CHOL"][frm] += 1
                            elif r in up:
                                upp_1["CHOL"][frm] += 1
                        elif r in resids[2]:
                            if r in lo:
                                low_1["DOPC"][frm] += 1
                            elif r in up:
                                upp_1["DOPC"][frm] += 1
                n_res.append(tmp_res)
            tmp_u = []
            tmp_l = []
            ###############################
            for frm_id in upp_1.index:
                        #upp_1.iloc[frm,:].divide(upp_1.iloc[frm,:].sum())
                # TODO this term here isn't behaving as I thought it should....
                tmp_u.append(self.check_states(upp_1.iloc[frm_id,:].divide(upp_1.iloc[frm_id,:].sum())
                                        , possible_states))
                tmp_l.append(self.check_states(low_1.iloc[frm_id,:].divide(low_1.iloc[frm_id,:].sum())
                                        , possible_states))
            states_u.append(tmp_u)
            states_l.append(tmp_l)
            # #lines = [int(l) for l in lines]
            # shell = []
            # # short, leaflet doesn't mater, gets first shell (every 3 lines form 1)
            # for l in lines[1::3]:
            #     shell.append([int(l.split()[-3]),int(l.split()[-2]),int(l.split()[-1])]) 
                 
            #     # Sanitiy check. values should be below 45, allowing up to 55
            #     if np.sum(int(l.split()[-3])+int(l.split()[-2])+int(l.split()[-1])) > 55:
            #         print("Unrealistic shell count. Confirm first shell is being analyzed")
            #         return None
            # # converts to numpy array for math
            # shell_arr = np.array(shell)
            # shell = np.array([shell_arr[:,0] / shell_arr.sum(axis=1),shell_arr[:,1] / shell_arr.sum(axis=1),shell_arr[:,2] / shell_arr.sum(axis=1)])
            # shell = shell.T
            # states = []
            
            # for s in shell:
            #     # determines all states for a file
            # 	states.append(self.check_states(s,possible_states))
            # full_states.append(states)
        pd.DataFrame(states_l).to_csv("%s/CG/data/states/long_low_%s_%s.csv"%(self.path, self.length, self.act))
        pd.DataFrame(states_u).to_csv("%s/CG/data/states/long_upp_%s_%s.csv"%(self.path, self.length, self.act))


lower = build = CGTM_Collect_Data(1,1,[],"cg","inact", "short",titrate='t').build_cg_short_titrate()

# CGTM_Collect_Data(1,3,[],"cg","in", "long").build_cg_long_states()
