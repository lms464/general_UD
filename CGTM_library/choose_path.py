#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 19 12:48:01 2021

@author: sharplm
"""
import os 
# import sys

def choose_path():
    censere_asym = "/Censere/UDel/resources/test_vor/data"
    censere_cg = None #needs path
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
# choose_path()