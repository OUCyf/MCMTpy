#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 20:11:10 2021

@author: Fu Yin (yinfu@mail.ustc.edu.cn) at USTC

This script:
    1) 1D green function solver when method=='pyfk';
    2) 3D green function solver when method=='sem'.  (it not work now)


Modify history:
    1) Mar 24 20:11:10 2021    ||    Fu Yin at USTC    ||    The initial release.
    2) ...
    
"""


import sys
from MCMTpy.gfs.pyfk_GFs.build_GFs_pyfk import build_GFs_pyfk
from MCMTpy.gfs.sem_GFs.build_GFs_sem import build_GFs_sem




#%%##########################################################################
#                         -----------------------
#                          1. build_GFs function
#                         -----------------------
#############################################################################

def build_GFs(filename,method):
    if method=='pyfk':
        build_GFs_pyfk(filename)
    elif method=='sem':
        build_GFs_sem(filename)
    else:
        sys.exit('MCMTpy: command error')













