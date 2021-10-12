#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 20:11:10 2021

@author: Fu Yin (yinfu@mail.ustc.edu.cn) at USTC

This script:
    1) 1D green function syn when method=='pyfk';
    2) 3D green function syn when method=='sem'.  (it not work now)


Modify history:
    1) Mar 24 20:11:10 2021    ||    Fu Yin at USTC    ||    The initial release.
    2) ...
"""


import sys
from MCMTpy.syn.pyfk_syn.syn_pyfk import syn_pyfk



#%%##########################################################################
#                         -----------------------
#                          1. syn function
#                         -----------------------
#############################################################################

def syn(filename,method):
    if method=='pyfk':
        _ = syn_pyfk(filename)
    elif method=='sem':
        pass
    else:
        sys.exit('MCMTpy: command error')



























