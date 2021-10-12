#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 20:11:10 2021

@author: Fu Yin (yinfu@mail.ustc.edu.cn) at USTC

This script:
    1) 



Modify history:
    1) Mar 24 20:11:10 2021    ||    Fu Yin at USTC    ||    The initial release.
    2) ...
"""


import sys
from MCMTpy.sampler.pyfk_MH.sample_MH import sample_MH
from MCMTpy.sampler.grid.sample_grid import sample_grid



#%%##########################################################################
#                         -----------------------
#                          1. sample function
#                         -----------------------
#############################################################################

def sample(filename,method):
    if method=='MH':
        sample_MH(filename)
    elif method=='grid': 
        sample_grid(filename)
    elif method=='Gibbs':                                                       # it not work now
        pass
    else:
        sys.exit('MCMTpy: command error')



























