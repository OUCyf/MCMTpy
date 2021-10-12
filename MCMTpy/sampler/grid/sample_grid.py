#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 15 14:34:20 2021

@author: yf
"""

import numpy as np
import math,os,sys
import random
import shutil
from mpi4py import MPI
from tqdm import tqdm

from MCMTpy.gfs.pyfk_GFs.read_json import read_GFs_json
from MCMTpy.sampler.grid.read_grid_json import read_Inv_json
from MCMTpy.sampler.grid.sampler_module import config_fk,get_GFs_sta_info,get_Sta_rand_data,write_inv_para
from MCMTpy.sampler.grid.sampler_module import get_MISFIT_1








#%%########################################################################
#                           -------------------
#                            1. main-function
#                           -------------------
###########################################################################w

def sample_grid(filename):
    
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()


    #########################################################################
    #                         -----------------
    #                           1.1. preprocess
    #                         -----------------
    #########################################################################

    if rank == 0:
        try:
            # read the inversion superparameter dict, station parameter dict, original data dict, source time function (obspy. Trace)
            Inv_para,Sta_inv_info,Sta_raw_data,Stf,source_info,CPU_splits = read_Inv_json(filename)
            print("Setting: read inversion.json successful in path",filename,flush=True)
        except Exception as inst:
            print(inst,flush=True)
            raise ValueError('read_Inv_json error')

        try:
            # read gfs database info
            Prepro_para,_,Source_Station_info_MPI = read_GFs_json(Inv_para["GFs_json_file"])
            print("Setting: read built_GFs.json successful in path",Inv_para["GFs_json_file"],flush=True)
        except Exception as inst:
            print(inst,flush=True)
            raise ValueError('read_GFs_json error')


        MPI_n = Inv_para['MPI_n']


        #------------ a. make directory ------------#
        if os.path.exists(Inv_para["Output_path"]):
            shutil.rmtree(Inv_para["Output_path"])
        if os.path.exists(Inv_para["Output_path"]) == False:
            os.makedirs(Inv_para["Output_path"])
            print("Setting: mkdir Output_path successful in path ",Inv_para["Output_path"],flush=True)


        #------------ b. output parameter info ------------#
        file_path = os.path.join(Inv_para["Output_path"],'Inv_para_info.txt') 
        file_path_source = os.path.join(Inv_para["Output_path"],'Inv_source_info.txt') 
        write_inv_para(file_path,Inv_para,Sta_inv_info,Sta_raw_data,Stf,file_path_source,source_info,CPU_splits)
        print("Setting: write Inv_para_info.txt successful in path ",Inv_para["Output_path"],flush=True)
        print("\n*****************************************************************",flush=True)
        print(f"Now begin grid-search with {MPI_n} cores ...\n",flush=True)


    else:
        Inv_para,Sta_inv_info,Sta_raw_data,Stf,Prepro_para,Source_Station_info_MPI,source_info,CPU_splits = [None for _ in range(8)]



    #########################################################################
    #                         ------------------------------
    #                          1.2. broadcast the variables
    #                         ------------------------------
    #########################################################################
    Inv_para = comm.bcast(Inv_para,root=0)
    Sta_inv_info = comm.bcast(Sta_inv_info,root=0)
    Sta_raw_data = comm.bcast(Sta_raw_data,root=0)
    Stf = comm.bcast(Stf,root=0)
    Prepro_para = comm.bcast(Prepro_para,root=0)
    Source_Station_info_MPI = comm.bcast(Source_Station_info_MPI,root=0)
    source_info = comm.bcast(source_info,root=0)
    CPU_splits = comm.bcast(CPU_splits,root=0)


    #########################################################################
    #                     -----------------------------------------
    #                     1.3. MPI MCMC_inv loop through each chunk
    #                     -----------------------------------------
    #########################################################################
    try:
        grid_inv(Inv_para,
                 Sta_inv_info,
                 Sta_raw_data,
                 Stf,
                 Prepro_para,
                 Source_Station_info_MPI,
                 source_info,
                 CPU_splits)
    except Exception as inst:
        print(inst,flush=True)
        raise ValueError('grid_inv error')



    #########################################################################
    #                         -----------------------
    #                          1.4. MPI comm.barrier
    #                         -----------------------
    #########################################################################
    comm.barrier()
    if rank == 0:
        print("\n\n*****************************************************************",flush=True)
        print("Successful !\nThe inversion process is completed.\n\n",flush=True)
        sys.exit()







#%%########################################################################
#                           -------------------
#                            2. sub-function
#                           -------------------
###########################################################################

def grid_inv(Inv_para,Sta_inv_info,Sta_raw_data,Stf,Prepro_para,Source_Station_info_MPI,source_info,CPU_splits):

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()


    #########################################################################
    #                         -----------------
    #                           2.1. preprocess
    #                         -----------------
    #########################################################################
    Dt = Inv_para["Dt"]                                                             # sampling rate of the raw data
    dt = Prepro_para["dt"]                                                          # sampling rate of the gfs database

    if Dt!=dt:
        raise ValueError('Error!\nGFs sample rate is not consistent with Raw_data!')

    InvType=Inv_para["InvType"]
    FM0=Inv_para["FM0"]

    source_prem,_,config_prem=config_fk(InvType)
    source_time_function = Stf





    #########################################################################
    #                         -----------------
    #                           2.2. inv
    #                         -----------------
    #########################################################################
    
    index_start = CPU_splits[rank][0]-1
    index_end = CPU_splits[rank][1]
    if rank==0:
        bar = tqdm(range(index_start,index_end,1),position=rank,file=sys.stdout)
        print("\tprint rank 0 as example:",flush=True)
    else:
        bar = range(index_start,index_end,1)

    FM_all=[]
    MISFIT_all=[]
    for jj in bar:
        #%% 2.2.1 inv process
        #------------ a. genrate new FM ------------#
        depth  = source_info[jj][0]
        strike  = source_info[jj][1]
        dip  = source_info[jj][2]
        rake  = source_info[jj][3]
        FM_new=[FM0[0],strike,dip,rake,FM0[4],FM0[5],depth,FM0[7]]                                                           # [m0, str, dip rake, lat, lon, depth, t0] in dc
        GFs_sta_info = get_GFs_sta_info(FM_new[-4:-1],Inv_para["NET_STA"],Prepro_para,Source_Station_info_MPI)


        #------------ % b.2.forward model || give FM and stf ------------#
        source_mechanism = np.array(FM_new[0:4])
        source_prem.update_source_mechanism(source_mechanism)


        #------------ % c.3.forward model syn，and stores all the stations data in a dictionary, and reads the TP Ts ------------#
        Sta_rand_data = get_Sta_rand_data(GFs_sta_info,config_prem,source_time_function)


        #------------ % c.4.MISFIT 1 ------------#
        MISFIT,Mw_mean = get_MISFIT_1(GFs_sta_info,Sta_raw_data,Sta_rand_data,Sta_inv_info,Inv_para,Prepro_para,FM_new[-1])


        #------------ % c.6. MH algorithm ------------#
        FM_new[0] = Mw_mean+FM_new[0]
        FM_all.append(FM_new)
        MISFIT_all.append(MISFIT)





        #------------ c.7 rank 0 display ------------#
        MISFIT_min=min(MISFIT_all)
        index_min=MISFIT_all.index(min(MISFIT_all))
        FM_min=FM_all[index_min]


        if rank == 0:
            bar.write(f"Iter {jj+1}: ")
            bar.write(f'MISFIT:{round(MISFIT,5)} || min MISFIT:{round(min(MISFIT_all),5)}')

            fm1=''
            for rr in range(0,len(FM_new)):
                fm1+=str(round(FM_new[rr],4))+'  '
            bar.write('FM: '+str(fm1))
            
            fm2=''
            for rr in range(0,len(FM_min)):
                fm2+=str(round(FM_min[rr],4))+'  '
            bar.write('FM_min: '+str(fm2))
            bar.write('\n')
            bar.set_description(f"Rank {rank}: now MISFIT={round(MISFIT,5)}") 






    #%% 2.2.3 output
    file_path = os.path.join(Inv_para["Output_path"],'rank'+'_'+str(rank)+'_FM_all')
    fp = open(file_path, 'w')
    for i in  range(0,len(FM_all),1):
        fp.write('{:.8f}'.format(MISFIT_all[i]))
        fp.write('\t')
        for j in  range(0,len(FM_all[0]),1):
            fp.write('{:.8f}'.format(FM_all[i][j]))
            fp.write('\t')
        fp.write('\n')
    fp.close()



    file_path = os.path.join(Inv_para["Output_path"],'rank'+'_'+str(rank)+'_FM_min')
    fp = open(file_path, 'w')
    fp.write('{:.8f}'.format(MISFIT_min))
    fp.write('\t')
    for i in range(0,len(FM_min),1):
        fp.write('{:.8f}'.format(FM_min[i]))
        fp.write('\t')
    fp.close()

