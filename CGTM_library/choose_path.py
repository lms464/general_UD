#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 19 12:48:01 2021

@author: sharplm
"""
import os 
# import sys

def choose_path():
    censere = "/Censere/UDel/resources/test_vor/data"
    deebo = '/home/sharplm/Asym/data/'
    if os.path.isdir(censere):
        return censere
    elif os.path.isdir(deebo):
        return deebo
    else:
        print("Your are using a computer not yet added to this list!")
        print("Please include new path(s) to check!")
        # sys.exit()
        os._exit(0) 
# choose_path()