#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 17 10:58:27 2021

@author: Fu Yin (yinfu@mail.ustc.edu.cn) at USTC

This script:
    1) function syn_pyfk.


Modify history:
    1) Mar 17 10:58:27 2021    ||    Fu Yin at USTC    ||    The initial release.
    2) ...
    
"""

import numpy as np
import os
import obspy
from obspy import UTCDateTime
from pyfk import generate_source_time_function
import json5

from MCMTpy.gfs.pyfk_GFs.read_json import read_GFs_json
from MCMTpy.utils.asdf_function import get_StationXML,Add_waveforms_head_info
from MCMTpy.utils.asdf_function import Add_waveforms,Add_stationxml,get_QuakeML,Add_quakeml
from MCMTpy.sampler.pyfk_MH.sampler_module import config_fk,get_GFs_sta_info,get_Sta_rand_data
from MCMTpy.utils.rotate import rotate





#%%########################################################################
#                           -------------------
#                            1. main-function
#                           -------------------
###########################################################################

def syn_pyfk(filename):

    syn_json= read_syn_json(filename)

    Sta_data_out = gen_syn(syn_json["GFs_json_file"],syn_json["srcType"], 
                           syn_json["point"], syn_json["source_mechanism"], 
                           syn_json["NET_STA"],
                           
            syn_json["Trapezoid_stf_mode"],syn_json["dura"],
            syn_json["rise"],syn_json["delta"], syn_json["Stf_file"],
            
            syn_json["filter_mode"],syn_json["freqmin"],
            syn_json["freqmax"],
            
            syn_json["Add_noise"], syn_json["noise_level_waveform"], 
            syn_json["noise_level_arrive_time"],
            
            syn_json["Output_mode"], syn_json["Output_path"], 
            syn_json["source_name"], syn_json["ASDF_filename"])


    return Sta_data_out






#%%########################################################################
#                           -------------------
#                            2. sub-function
#                           -------------------
###########################################################################

def read_syn_json(filename):
    with open(filename, 'r',encoding='utf-8') as f1:
        syn_json_r = json5.load(f1)
    syn_json=syn_json_r

    return syn_json



##################################################
def gen_syn(GFs_json_file,srcType, point, source_mechanism, NET_STA,
            Trapezoid_stf_mode=None,dura=None,rise=None,delta=None, Stf_file=None,
            filter_mode=None,freqmin=None,freqmax=None,
            Add_noise=None, noise_level_waveform=None, noise_level_arrive_time=None,
            Output_mode=None, Output_path=None, source_name=None, ASDF_filename=None):

    #------------ 2.1. read gfs information ------------#
    Prepro_para,_,Source_Station_info_MPI = read_GFs_json(GFs_json_file)
    source_prem,v_model,config_prem=config_fk(srcType)


    #------------ 2.2. gets info of all stations for a given single source ------------#
    GFs_sta_info = get_GFs_sta_info(point,NET_STA,Prepro_para,Source_Station_info_MPI)


    #------------ 2.3. give FM and stf ------------#
    if Trapezoid_stf_mode==True:
        source_time_function = generate_source_time_function(dura=dura, rise=rise, delta=delta)
    else:
        source_time_function = obspy.read(Stf_file)[0]

    if srcType=='mt':
        # a):
        #   The seismic moment given by the JSON file is MW magnitude.
        #   pyfk uses dyne-cm units for the synthesis of the moment tensor MT solution, which is converted to M0
        mw = source_mechanism[0]
        source_mechanism[0] = np.power(10., 1.5 * mw + 16.1)
    elif srcType=='dc':
        # b):
        #   when == dc, pyfk uses moment magnitude, which does not need to be changed
        source_mechanism = source_mechanism
    elif srcType=='sf':
        # c):
        #   The seismic moment given by the JSON file is MW magnitude.
        #   pyfk uses dyne units when synthesing the moment tensor sf solution, which is converted to M0
        mw = source_mechanism[0]
        source_mechanism[0] = np.power(10., 1.5 * mw + 16.1)
    elif srcType=='ep':
        # d):
        #   The seismic moment given by the JSON file is MW magnitude.
        #   pyfk uses dyne-cm units when synthesing the moment tensor ep solution, which is converted to M0
        mw = source_mechanism[0]
        source_mechanism[0] = np.power(10., 1.5 * mw + 16.1)

    source_prem.update_source_mechanism(source_mechanism)


    #------------ 2.4.forward model syn，and stores all the stations data in a dictionary, and reads the TP Ts ------------#
    Sta_data = get_Sta_rand_data(GFs_sta_info,config_prem,source_time_function)


    #------------ 2.5. filter ------------#
    Station_num = len(NET_STA)
    if filter_mode==True:
        net_sta_name = list(Sta_data.keys())
        for i in range(0,Station_num,1):
            stream = Sta_data[net_sta_name[i]]['data']
            stream.filter('bandpass',
                          freqmin=freqmin,
                          freqmax=freqmax,
                          corners=4,
                          zerophase=True)


    #------------ 2.6. add noise ------------#
    Sta_data_noise =  Sta_data.copy()
    Station_num = len(NET_STA)
    
    if Add_noise == True:
        net_sta_name = list(Sta_data.keys())
        for i in range(0,Station_num,1):
            tp = Sta_data[net_sta_name[i]]['tp']
            ts = Sta_data[net_sta_name[i]]['ts']
            stream = Sta_data[net_sta_name[i]]['data']
            # a. add noise in waveform
            for j in range(0,3,1):
                max_sac_value = np.max(stream[j].data)
                npts = stream[j].stats.npts
                white_noise = max_sac_value*noise_level_waveform*np.random.standard_normal(size=npts)
                Sta_data_noise[net_sta_name[i]]['data'][j].data = stream[j].data+white_noise
            # b. add noise in travel time
            arrive_time_noise = noise_level_arrive_time*np.random.standard_normal(size=1)
            Sta_data_noise[net_sta_name[i]]['tp'] = tp + arrive_time_noise[0]
            Sta_data_noise[net_sta_name[i]]['ts'] = ts + arrive_time_noise[0]

        Sta_data_out = Sta_data_noise
    else:
        Sta_data_out = Sta_data



    #------------ 2.7. output data ------------#
    if Output_mode == True:
        if not os.path.isdir(Output_path):os.mkdir(Output_path)

        source_lat = point[0]
        source_lon = point[1]
        source_depth = point[2]
        source_time = UTCDateTime()

        Output_file_path = os.path.join(Output_path,ASDF_filename+'.h5')
        catalog_create=get_QuakeML(source_name,source_lat, source_lon, source_depth,source_time)
        Add_quakeml(Output_file_path,catalog_create)

        for i in range(0,Station_num,1):
            net_sta_name = str(NET_STA[i][0]+'.'+NET_STA[i][1])
            station_network = NET_STA[i][0]
            station_name = NET_STA[i][1]
            station_lat = NET_STA[i][2]
            station_lon = NET_STA[i][3]
            station_depth = NET_STA[i][4]
            station_distance = GFs_sta_info[net_sta_name][5]
            station_az = GFs_sta_info[net_sta_name][7]
            station_baz = GFs_sta_info[net_sta_name][8]
            tp = Sta_data_out[net_sta_name]['tp']
            ts = Sta_data_out[net_sta_name]['ts']
            stream = Sta_data_out[net_sta_name]['data']
            time_before_first_arrival = Prepro_para["samples_before_first_arrival"]*Prepro_para["dt"]
            starttime = source_time+tp-time_before_first_arrival

            # rotate ZRT->ENZ
            for chn in range(0,len(stream)):
                stream[chn].stats.back_azimuth = float(station_baz)
            stream = rotate(data_raw=stream, chn=['ZRT','ENZ'])
            
            # write asdf
            inv=get_StationXML(station_network,
                               station_name,
                               station_lat,
                               station_lon,
                               station_depth,
                               srcType='enz' )
            
            Add_waveforms_head_info([stream],
                                    station_network,
                                    station_name,
                                    srcType='enz',
                                    starttime=starttime)
            
            print("now syn "+station_name+':\n')
            print(stream)
            print('\n\n')
            
            Add_waveforms([stream],
                          catalog_create,
                          Output_file_path,
                          source_name,
                          tp,
                          ts,
                          station_distance)
            
            Add_stationxml(inv,Output_file_path)
            
            # write sac
            for j in range(0,3,1):
                CHN = stream[j].stats.channel
                sac_filename = os.path.join(Output_path,
                                            source_name+'.'+station_network+'.'+station_name+'.'+CHN+'.sac')
                stream[j].stats.sac["b"]  = tp-time_before_first_arrival
                stream[j].stats.sac["o"]  = 0.0
                stream[j].stats.sac["t1"] = tp
                stream[j].stats.sac["t2"] = ts
                stream[j].write(sac_filename, format='SAC')


    return Sta_data_out











