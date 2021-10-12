#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 20:11:10 2021

@author: Fu Yin (yinfu@mail.ustc.edu.cn) at USTC

This script:
    1) function build_GFs_pyfk use mpi4py and pyasdf package to cal_gfs.


Modify history:
    1) Mar 24 20:11:10 2021    ||    Fu Yin at USTC    ||    The initial release.
    2) ...
    
"""


import os,sys
from mpi4py import MPI
from MCMTpy.gfs.pyfk_GFs.read_json import read_GFs_json,write_Prepro_para,write_GFs_info,write_Source_Station_info_MPI
from MCMTpy.gfs.pyfk_GFs.cal_gfs import cal_gfs




#%%###########################################

def build_GFs_pyfk(filename):

    filename = os.path.join(filename)                                           # for example: filename = "./S1_built_GFs.json"
    comm     = MPI.COMM_WORLD
    rank     = comm.Get_rank()
    size     = comm.Get_size()


    #%%######################################################################
    #                         -----------------
    #                           1. preprocess
    #                         -----------------
    #########################################################################
    if rank == 0:
        try:
            Prepro_para,GFs_info,Source_Station_info_MPI = read_GFs_json(filename)
            print("Setting: read built_GFs.json successful in path",filename,flush=True)

        except Exception as inst:
            print(inst,flush=True)
            raise ValueError('read_GFs_json error')

        MPI_n = Prepro_para['MPI_n']


        #------------ 1.1. make directory ------------#
        if not os.path.isdir(Prepro_para["DATADIR"]):
            os.mkdir(Prepro_para["DATADIR"])
            print("Setting: mkdir DATADIR successful in path ",Prepro_para["DATADIR"],flush=True)
        if not os.path.isdir(Prepro_para["DATADIR_splits"]):
            os.mkdir(Prepro_para["DATADIR_splits"])
            print("Setting: mkdir DATADIR_splits successful in path ",Prepro_para["DATADIR_splits"],flush=True)


        #------------ 1.2. output parameter info ------------#
        file_path = os.path.join(Prepro_para["DATADIR"],'prepro_para_info.txt') 
        write_Prepro_para(file_path,Prepro_para)
        print("Setting: write_Prepro_para successful in path ",file_path,flush=True)
        
        if Prepro_para["Database_mode"]==True:
            file_path = os.path.join(Prepro_para["DATADIR"],'GFs_info.txt') 
            write_GFs_info(file_path,GFs_info)
            print("Setting: write_GFs_info successful in path ",file_path,flush=True)
        else:
            file_path = os.path.join(Prepro_para["DATADIR"],'write_Source_Station_info_MPI.txt') 
            write_Source_Station_info_MPI(file_path,Source_Station_info_MPI)
            print("Setting: write_Source_Station_info_MPI successful in path ",file_path,flush=True)
        
        print("\n*****************************************************************",flush=True)
        print(f"Now begin to calculate Green's function with {MPI_n} cores...",flush=True)
    
    else:
        Prepro_para,GFs_info,Source_Station_info_MPI = [None for _ in range(3)]





    #%%######################################################################
    #                         ----------------------------
    #                          2. broadcast the variables
    #                         ----------------------------
    #########################################################################
    
    Prepro_para             = comm.bcast(Prepro_para,root=0)
    GFs_info                = comm.bcast(GFs_info,root=0)
    Source_Station_info_MPI = comm.bcast(Source_Station_info_MPI,root=0)





    #%%######################################################################
    #                         ----------------------------------------
    #                          3. MPI cal_gfs loop through each chunk
    #                         ----------------------------------------
    #########################################################################

    try:
        cal_gfs(Prepro_para,GFs_info,Source_Station_info_MPI)
    except Exception as inst:
        print(inst,flush=True)
        raise ValueError('cal_gfs error')





    #%%######################################################################
    #                         ---------------------
    #                          4. MPI comm.barrier
    #                         ---------------------
    #########################################################################

    comm.barrier()
    if rank == 0:
        print("\n\n*****************************************************************",flush=True)
        print("Successful !\nAll of the Green's functions have been computed.\n\n",flush=True)
        sys.exit()












