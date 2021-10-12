#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 26 09:15:33 2021

@author: Fu Yin (yinfu@mail.ustc.edu.cn) at USTC

This script:
    1) function cal_gfs.


Modify history:
    1) Mar 25 16:16:07 2021    ||    Fu Yin at USTC    ||    The initial release.
    2) ...
    
"""

import numpy as np
import os
import sys
from obspy import UTCDateTime
from pyfk import SourceModel, SeisModel, Config
from pyfk import calculate_gf
from mpi4py import MPI
from tqdm import tqdm
from MCMTpy.utils.asdf_function import get_StationXML,Add_waveforms_head_info
from MCMTpy.utils.asdf_function import Add_waveforms,Add_stationxml,get_QuakeML,Add_quakeml



############################################
def cal_gfs(Prepro_para,GFs_info=None,Source_Station_info_MPI=None):

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()


    #%%######################################################################
    #                         -------------------
    #                           1. Database mode
    #                         -------------------
    #########################################################################
    if Prepro_para["Database_mode"]==True:
        
        index_start = Prepro_para["CPU_splits"][rank][0]-1
        index_end = Prepro_para["CPU_splits"][rank][1]
        
        # a. Multi-process display || position=rank represents the current number of processes
        # bar = tqdm(range(index_start,index_end,1),position=rank,file=sys.stdout)
        
        # b. rank 0 display
        if rank==0:
            bar = tqdm(range(index_start,index_end,1),position=rank,file=sys.stdout)
            print("\tprint rank 0 as example:",flush=True)
        else:
            bar = range(index_start,index_end,1)
        
        
        for i in bar:
            #--------- 1.1 lat/lon ---------#
            source_name  = GFs_info[i][1]
            source_lat   = 50                                                   # Any given
            source_lon   = 50                                                   # Any given
            source_depth = GFs_info[i][2]
            source_time  = UTCDateTime()
            time_before_first_arrival = Prepro_para["samples_before_first_arrival"]*Prepro_para["dt"]
            
            station_network  = GFs_info[i][3]
            station_name     = GFs_info[i][4]
            station_lat      = 60                                               # Any given
            station_lon      = 60                                               # Any given
            station_depth    = GFs_info[i][6]
            station_distance = GFs_info[i][5]
            
            
            #--------- 1.2 pyfk model ---------#
            source_prem = SourceModel(
                sdep              = source_depth,
                srcType           = Prepro_para["srcType"],
                source_mechanism  = Prepro_para["source_mechanism"])
            
            v_model = SeisModel(
                model             = Prepro_para["Velocity_model_value"],
                flattening        = Prepro_para["flattening"])
            
            config_prem = Config(
                model             = v_model,
                source            = source_prem,
                npt               = Prepro_para["npt"],
                dt                = Prepro_para["dt"],
                receiver_distance = np.array([1]),                              # Any given and change later
                degrees           = Prepro_para["degrees"],
                taper             = Prepro_para["taper"], 
                filter            = Prepro_para["filter"], 
                dk                = Prepro_para["dk"], 
                smth              = Prepro_para["smth"], 
                pmin              = Prepro_para["pmin"], 
                pmax              = Prepro_para["pmax"], 
                kmax              = Prepro_para["kmax"], 
                rdep              = station_depth, 
                updn              = Prepro_para["updn"], 
                samples_before_first_arrival =Prepro_para["samples_before_first_arrival"])
            
            
            #--------- 1.3 calculate gfs ---------#
            if Prepro_para["degrees"]==True:
                config_prem.receiver_distance = np.array([station_distance])*111.19
            else:
                config_prem.receiver_distance = np.array([station_distance])
            
            gf = calculate_gf(config_prem)
            
            tp = gf[0][0].stats.sac['t1']                                       # tp time
            ts = gf[0][0].stats.sac['t2']                                       # ts time
            starttime = source_time + tp - time_before_first_arrival
            
            
            #--------- 1.4 output gfs ---------#
            ff = os.path.join(Prepro_para["DATADIR_splits"],source_name)        # folder name of all Green's functions for one depth
            if not os.path.isdir(ff):os.mkdir(ff)
            file_path = os.path.join( ff,source_name+'_'+station_name+'.h5')    # gfs for a depth earthquake at one station
            
            catalog_create = get_QuakeML(source_name,source_lat, source_lon, source_depth,source_time)
            Add_quakeml(file_path ,catalog_create)
            
            inv = get_StationXML(station_network, station_name,station_lat, 
                               station_lon, station_depth, Prepro_para["srcType"])
            
            Add_waveforms_head_info(gf,station_network,station_name,
                                    Prepro_para["srcType"],starttime)
            
            Add_waveforms(gf, catalog_create,file_path,source_name,tp,ts,station_distance)
            
            Add_stationxml(inv,file_path)
            
            
            #--------- 1.5 display ---------#
            # a. Multi-process
            # bar.set_description(f"Rank {rank}: now get GFs {i}th of [{source_name} {station_name}]") 
            
            # b. rank 0
            if rank == 0:
                bar.write(f"Rank {rank}: iteration number {i}th: ")
                bar.set_description(f"Rank {rank}: now get GFs {i}th of [{source_name} {station_name}]") 





    #%%######################################################################
    #                         ------------------
    #                           2. single mode
    #                         ------------------
    #########################################################################
    else:
        index_start = Prepro_para["CPU_splits"][rank][0]-1
        index_end = Prepro_para["CPU_splits"][rank][1]
        
        # a. Multi-process display || position=rank represents the current number of processes
        # bar = tqdm(range(index_start,index_end,1),position=rank,file=sys.stdout)
        
        # b. rank 0 display
        if rank==0:
            bar = tqdm(range(index_start,index_end,1),position=rank,file=sys.stdout)
            print("\tprint rank 0 as example:",flush=True)
        else:
            bar = range(index_start,index_end,1)
        
        
        for i in bar:
            station_num = len(Source_Station_info_MPI[i])
            source_name = Source_Station_info_MPI[i][0][0]
            source_lat = Source_Station_info_MPI[i][0][1]
            source_lon = Source_Station_info_MPI[i][0][2]
            source_depth = Source_Station_info_MPI[i][0][3]
            source_time = UTCDateTime()
            time_before_first_arrival = Prepro_para["samples_before_first_arrival"]*Prepro_para["dt"]
            
            # Green's function for an earthquake of all stations
            file_path = os.path.join( Prepro_para["DATADIR_splits"],source_name+'.h5')
            catalog_create=get_QuakeML(source_name,source_lat, source_lon, source_depth,source_time)
            Add_quakeml(file_path ,catalog_create)
            
            for j in range(1,station_num,1):
                #--------- 1.1 lat/lon ---------#
                station_network = Source_Station_info_MPI[i][j][0]
                station_name = Source_Station_info_MPI[i][j][1]
                station_lat = Source_Station_info_MPI[i][j][2]
                station_lon = Source_Station_info_MPI[i][j][3]
                station_depth = Source_Station_info_MPI[i][j][4]
                station_distance = Source_Station_info_MPI[i][j][5]
                station_az = Source_Station_info_MPI[i][j][6]
                station_baz = Source_Station_info_MPI[i][j][7]
                
                
                #--------- 1.2 pyfk model ---------#
                source_prem = SourceModel(
                    sdep=source_depth,
                    srcType=Prepro_para["srcType"],
                    source_mechanism=Prepro_para["source_mechanism"])
        
                v_model = SeisModel(
                    model=Prepro_para["Velocity_model_value"],
                    flattening=Prepro_para["flattening"])
        
                config_prem = Config(
                    model=v_model,
                    source=source_prem,
                    npt=Prepro_para["npt"],
                    dt=Prepro_para["dt"],
                    receiver_distance=np.array([1]),                            # Any given and change later
                    degrees=Prepro_para["degrees"],
                    taper=Prepro_para["taper"], 
                    filter=Prepro_para["filter"], 
                    dk=Prepro_para["dk"], 
                    smth=Prepro_para["smth"], 
                    pmin=Prepro_para["pmin"], 
                    pmax=Prepro_para["pmax"], 
                    kmax=Prepro_para["kmax"], 
                    rdep=station_depth, 
                    updn=Prepro_para["updn"], 
                    samples_before_first_arrival=Prepro_para["samples_before_first_arrival"])
                
                
                #--------- 1.3 calculate gfs ---------#
                if Prepro_para["degrees"]==True:
                    config_prem.receiver_distance = np.array([station_distance])*111.19
                else:
                    config_prem.receiver_distance = np.array([station_distance])
                gf = calculate_gf(config_prem)
                tp=gf[0][0].stats.sac['t1']
                ts=gf[0][0].stats.sac['t2']
                starttime = source_time+tp-time_before_first_arrival
                
                
                #--------- 1.4 output gfs ---------#
                inv=get_StationXML(station_network, station_name,station_lat, station_lon, station_depth, Prepro_para["srcType"])
                Add_waveforms_head_info(gf,station_network,station_name,Prepro_para["srcType"],starttime)
                
                Add_waveforms(gf,catalog_create,file_path,source_name,tp,ts,
                              station_distance,station_az,station_baz)
                Add_stationxml(inv,file_path)
            
            
            #--------- 1.5 display ---------#
            # a. Multi-process
            # bar.set_description(f"Rank {rank}: now get GFs {i}th of [{source_name}]") 
            
            # b. rank 0
            if rank == 0:
                bar.write(f"Rank {rank}: iteration number {i}th: ")
                bar.set_description(f"Rank {rank}: now get GFs {i}th of [{source_name}]") 









