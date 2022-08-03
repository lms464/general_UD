#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 19 12:48:01 2021

@author: sharplm
"""
import os 

def choose_path( inpt = None ):
    
        ## Use inpt for coarse grained systems only
        ## This will not work for the asymetric as of now
    
        if inpt == None:
            censere_asym = "/Censere/UDel/resources/test_vor/data"
            censere_cg = "/home/liam/lms464/"#"/Censere/UDel/"
            censere = [censere_asym, censere_cg]
            deebo_asym = '/home/sharplm/Asym/data/'
            deebo_cg = "/home/sharplm"
            deebo = [deebo_asym,deebo_cg]
            if os.path.isdir(censere_asym):
                return censere
            elif os.path.isdir(deebo_asym) or os.path.isdir(deebo_cg):
                return deebo
            else:
                print("Your are using a computer not yet added to this list!")
                print("Please include new path(s) to check!")
                # sys.exit()
                os._exit(0) 
        else:
            return [inpt]*2
# print(choose_path("sup"))
