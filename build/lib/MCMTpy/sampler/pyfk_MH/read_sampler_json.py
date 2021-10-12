#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  1 11:51:12 2021

@author: Fu Yin (yinfu@mail.ustc.edu.cn) at USTC

This script:
    1) function read_Inv_json.


Modify history:
    1) Apr  1 11:51:12 2021    ||    Fu Yin at USTC    ||    The initial release.
    2) ...
    
"""

import json5
import obspy
import pyasdf
from pyfk import generate_source_time_function




#%%########################################################################
#                           -------------------
#                            1. main-function
#                           -------------------
###########################################################################
def read_Inv_json(filename):
    
    #------------ 1.1. read json parameters ------------#
    with open(filename, 'r',encoding='utf-8') as f1:
        Inversion_json_r = json5.load(f1)
    Inversion_json = Inversion_json_r.copy()


    #------------ 1.2. read the inversion parameter files of each station ——> Sta_inv_info（dict） ------------#
    if Inversion_json["Use_file_mode"] == False:
        write_raw_data_inv_info(Inversion_json)                                 # write raw_data_inv_info.txt based on the inversion parameters
        Sta_inv_info = read_raw_data_inv_info(Inversion_json["Raw_data_inv_info_writed"])
    else:
        Sta_inv_info = read_raw_data_inv_info(Inversion_json["Raw_data_inv_info"]) # read Raw_data_inv_info.txt
    Station_num = len(Sta_inv_info)                                               # number of stations used in inversion


    #------------ 1.3. read raw data and arrival information——> Sta_data（dict） ------------#
    Sta_data = read_raw_data_file(Inversion_json["Raw_data_file"],Inversion_json["Source_tag"])


    #------------ 1.4. read stf——> Stf（obspy.trace) ------------#
    Stf = get_source_time_function(Inversion_json)


    #------------ 1.5. sum all parameter info ——> Inv_para（dict） ------------#
    Inv_para = { 
        # a. path
        "GFs_json_file":Inversion_json["GFs_json_file"],
        "GFs_file": Inversion_json["GFs_file"],
        "Raw_data_file": Inversion_json["Raw_data_file"],
        "Source_tag": Inversion_json["Source_tag"],
        "Raw_p_n0": Inversion_json["Raw_p_n0"],
        "Dt": Inversion_json["Dt"],
        "Output_path": Inversion_json["Output_path"],


        # b. stf
        "Trapezoid_stf_mode":  Inversion_json["Trapezoid_stf_mode"],
        "dura": Inversion_json["dura"],
        "rise": Inversion_json["rise"],      
        "delta": Inversion_json["delta"],
        "Stf_file": Inversion_json["Stf_file"],


        # c. compute times and parallel cores
        "MPI_n": Inversion_json["MPI_n"],
        "Chains_n": Inversion_json["Chains_n"],
        "N": Inversion_json["N"],


        # d. objective function parameter
        "Geometrical_spreading_mode": Inversion_json["Geometrical_spreading_mode"],
        "Distance_0": Inversion_json["Distance_0"],
        "Amplitude_normalization_mode": Inversion_json["Amplitude_normalization_mode"],

        "Use_file_mode": Inversion_json["Use_file_mode"],
        "Raw_data_inv_info": Inversion_json["Use_file_mode"],
        "Raw_data_inv_info_writed": Inversion_json["Use_file_mode"],

        #//--------------- Raw_data_inv_info -----------------
        "NET_STA":Inversion_json["NET_STA"],

        "P_inv_comp":  Inversion_json["P_inv_comp"],
        "S_inv_comp" : Inversion_json["S_inv_comp"],
        "Surf_inv_comp" : Inversion_json["Surf_inv_comp"],

        "P_win":  Inversion_json["P_win"],
        "S_win":  Inversion_json["S_win"],
        "Surf_win":  Inversion_json["Surf_win"],

        "Phase_weight":  Inversion_json["Phase_weight"],

        "Distance_weight":  Inversion_json["Distance_weight"],

        "P_maxlag":  Inversion_json["P_maxlag"],
        "S_maxlag":  Inversion_json["S_maxlag"],
        "Surf_maxlag":  Inversion_json["Surf_maxlag"],

        "P_filter": Inversion_json["P_filter"],
        "S_filter": Inversion_json["S_filter"],
        "Surf_filter":Inversion_json["Surf_filter"],
        #//--------------- Raw_data_inv_info -----------------


        # Inversion parameters
        "InvType": Inversion_json["InvType"],
        "FM0": Inversion_json["FM0"],
        "Fixed_FM": Inversion_json["Fixed_FM"],
        "FM_boundary": Inversion_json["FM_boundary"],
        "Noise_level_waveform": Inversion_json["Noise_level_waveform"],
        "Noise_level_phasetime": Inversion_json["Noise_level_phasetime"],

        "N_k":Inversion_json["N_k"],
        "N_mag":Inversion_json["N_mag"],

        # The standard deviation of the new solution of Mw. (!!******!!)
        "sigma_mw": Inversion_json["sigma_mw"],

        # The standard deviation of the new solution of Moment Tensor.
        "sigma_mt": Inversion_json["sigma_mt"],

        # The standard deviation of the new solution of Strike Dip Rake. (!!******!!)
        "sigma_dc": Inversion_json["sigma_dc"],

        # The standard deviation of the new solution of Single Couple (not test now).
        "sigma_sf": Inversion_json["sigma_sf"],

        # The standard deviation of the new solution of Latitude and Longitude (degree). (!!******!!)
        "sigma_lat_lon": Inversion_json["sigma_lat_lon"],

        # The standard deviation of the new solution of Source Depth (km). (!!******!!)
        "sigma_depth": Inversion_json["sigma_depth"],

        # The standard deviation of the new solution of T0 (s). (!!******!!)
        "sigma_t": Inversion_json["sigma_t"],

        # e. Additional parameters added
        "Station_num":Station_num,

}


    return Inv_para,Sta_inv_info,Sta_data,Stf




#%%########################################################################
#                           -------------------
#                            2. sub-function
#                           -------------------
###########################################################################

def read_raw_data_inv_info(filename):
    '''
    read the inversion parameter files of each station
    '''
    Raw_data_inv_info = []
    one_sta_info=[]                                                             # store a station corresponding to the source
    with open(filename,'r') as f:
        for line in f:
            line = line.strip()
            line = line.strip('\t')
            line = line.strip('\n')
            line = line.split()
            if line == [] and one_sta_info != []:   
                Raw_data_inv_info.append(one_sta_info)
                one_sta_info=[]
            elif len(list(line))>=2:                                            # if the line is not blank, the information is added
                one_sta_info.append( list( line ))

    if one_sta_info != []:
        Raw_data_inv_info.append(one_sta_info)

    Sta_inv_info = {}
    sta_num=len(Raw_data_inv_info)
    for i in range(0,sta_num,1):
        info={}
        for j in range(1,16,1):                                                 #  16 is artificial
            info_list = Raw_data_inv_info[i][j][1:]
            for k in range(0,len(info_list),1):
                if info_list[k]!='no' and info_list[k]!='yes':
                    info_list[k]=float(info_list[k])
            info.update({Raw_data_inv_info[i][j][0]:info_list  }) 

        net_sta_name = Raw_data_inv_info[i][0][0]+'.'+Raw_data_inv_info[i][0][1]
        Sta_inv_info.update({ net_sta_name: info}) 


    return Sta_inv_info




############################################
def read_raw_data_file(Raw_data_file,Source_tag):
    '''
    Read the raw waveform data and add it to the STA_DATA dictionary
    '''
    data_path = Raw_data_file

    ds_raw=pyasdf.ASDFDataSet(data_path,mpi=False,mode='r')                     # make sure that mpi = False
    Station_list = ds_raw.waveforms.list()                                      # a list of all stations in an ASDF file
    Station_num = len(Station_list)

    Sta_data={}                                                                 # defining an empty dictionary
    for i in range(0,Station_num,1):
        info={}
        tags = ds_raw.waveforms[Station_list[i]].get_waveform_tags()            # by default there is only one tag, that is, an ASDF file contains only one event
        if Source_tag in tags:                                                  # if the station contains this event, add it
            raw_sac = ds_raw.waveforms[Station_list[i]][Source_tag]
            tp = float( ds_raw.waveforms[Station_list[i]][Source_tag][0].stats.asdf['labels'][1] )
            ts = float( ds_raw.waveforms[Station_list[i]][Source_tag][0].stats.asdf['labels'][3] )
            info.update({ "tp": tp})  
            info.update({ "ts": ts})
            info.update({ "data": raw_sac})
            Sta_data.update({ Station_list[i]: info})


    return Sta_data




############################################
def get_source_time_function(Inversion_json):
    '''
    read stf
    '''
    Trapezoid_stf_mode= Inversion_json["Trapezoid_stf_mode"]
    if Trapezoid_stf_mode==True:
        dura= Inversion_json["dura"]
        rise=Inversion_json["rise"]
        delta= Inversion_json["delta"]
        Stf=generate_source_time_function(dura=dura, rise=rise, delta=delta)
    
    else:
        Stf_file=Inversion_json["Stf_file"]
        Stf=obspy.read(Stf_file)[0]


    return Stf




############################################
def write_raw_data_inv_info(Inversion_json):
    filename = Inversion_json["Raw_data_inv_info_writed"]
    NET_STA=Inversion_json["NET_STA"]
    Sta_num = len(NET_STA)

    P_inv_comp = Inversion_json["P_inv_comp"]
    S_inv_comp = Inversion_json["S_inv_comp"]
    Surf_inv_comp = Inversion_json["Surf_inv_comp"]

    P_win = Inversion_json["P_win"]
    S_win = Inversion_json["S_win"]
    Surf_win = Inversion_json["Surf_win"]

    Phase_weight = Inversion_json["Phase_weight"]
    Distance_weight= Inversion_json["Distance_weight"]

    P_maxlag = Inversion_json["P_maxlag"]
    S_maxlag = Inversion_json["S_maxlag"]
    Surf_maxlag = Inversion_json["Surf_maxlag"]

    P_filter= Inversion_json["P_filter"]
    S_filter= Inversion_json["S_filter"]
    Surf_filter= Inversion_json["Surf_filter"]

    fp = open(filename, 'w')
    for  i in  range(0,Sta_num,1):
        ######################
        fp.write(NET_STA[i][0])
        fp.write('   ')
        fp.write(NET_STA[i][1])
        fp.write('   ')
        fp.write('\n')

        ######################
        fp.write('Lat_lon_depth')
        fp.write('\t')
        fp.write('{:f}'.format(NET_STA[i][2]))
        fp.write('   ')
        fp.write('{:f}'.format(NET_STA[i][3]))
        fp.write('   ')
        fp.write('{:f}'.format(NET_STA[i][4]))
        fp.write('   ')
        fp.write('\n')

        ######################
        fp.write('P_inv_comp')
        fp.write('\t')
        fp.write('{:s}'.format(P_inv_comp[0]))
        fp.write('\t')
        fp.write('{:s}'.format(P_inv_comp[1]))
        fp.write('\t')
        fp.write('{:s}'.format(P_inv_comp[2]))
        fp.write('\t')
        fp.write('\n')

        ######################
        fp.write('S_inv_comp')
        fp.write('\t')
        fp.write('{:s}'.format(S_inv_comp[0]))
        fp.write('\t')
        fp.write('{:s}'.format(S_inv_comp[1]))
        fp.write('\t')
        fp.write('{:s}'.format(S_inv_comp[2]))
        fp.write('\t')
        fp.write('\n')

        ######################
        fp.write('Surf_inv_comp')
        fp.write('\t')
        fp.write('{:s}'.format(Surf_inv_comp[0]))
        fp.write('\t')
        fp.write('{:s}'.format(Surf_inv_comp[1]))
        fp.write('\t')
        fp.write('{:s}'.format(Surf_inv_comp[2]))
        fp.write('\t')
        fp.write('\n')

        ######################
        fp.write('P_win')
        fp.write('\t')
        fp.write('\t')
        if P_win[0]=='no':
            fp.write('{:s}'.format(P_win[0]))
        else:
            fp.write('{:.2f}'.format(P_win[0]))
            fp.write('\t')
        if P_win[1]=='no':
            fp.write('{:s}'.format(P_win[1]))
        else:
            fp.write('{:.2f}'.format(P_win[1]))
        fp.write('\t')
        fp.write('\n')

        ######################
        fp.write('S_win')
        fp.write('\t')
        fp.write('\t')
        if S_win[0]=='no':
            fp.write('{:s}'.format(S_win[0]))
        else:
            fp.write('{:.2f}'.format(S_win[0]))
        fp.write('\t')
        if S_win[1]=='no':
            fp.write('{:s}'.format(S_win[1]))
        else:
            fp.write('{:.2f}'.format(S_win[1]))
        fp.write('\t')
        fp.write('\n')

        ######################
        fp.write('Surf_win')
        fp.write('\t')
        fp.write('\t')
        if Surf_win[0]=='no':
            fp.write('{:s}'.format(Surf_win[0]))
        else:
            fp.write('{:.2f}'.format(Surf_win[0]))
        fp.write('\t')
        if Surf_win[1]=='no':
            fp.write('{:s}'.format(Surf_win[1]))
        else:
            fp.write('{:.2f}'.format(Surf_win[1]))
        fp.write('\t')
        fp.write('\n')

        ######################
        fp.write('Phase_weight')
        fp.write('\t')
        if Phase_weight[0]=='no':
            fp.write('{:s}'.format(Phase_weight[0]))
        else:
            fp.write('{:.1f}'.format(Phase_weight[0]))
        fp.write('\t')
        if Phase_weight[1]=='no':
            fp.write('{:s}'.format(Phase_weight[1]))
        else:
            fp.write('{:.1f}'.format(Phase_weight[1]))
        fp.write('\t')
        if Phase_weight[2]=='no':
            fp.write('{:s}'.format(Phase_weight[2]))
        else:
            fp.write('{:.1f}'.format(Phase_weight[2]))
        fp.write('\t')
        fp.write('\n')

        ######################
        fp.write('Distance_weight')
        fp.write('\t')
        if Distance_weight[i]=='no':
            fp.write('{:s}'.format(Distance_weight[i]))
        else:
            fp.write('{:.1f}'.format(Distance_weight[i]))
        fp.write('\t')
        fp.write('\n')

        ######################
        fp.write('P_maxlag')
        fp.write('\t')
        fp.write('\t')
        if P_maxlag=='no':
            fp.write('{:s}'.format(P_maxlag))
        else:
            fp.write('{:.4f}'.format(P_maxlag))
        fp.write('\t')
        fp.write('\n')

        ######################
        fp.write('S_maxlag')
        fp.write('\t')
        fp.write('\t')
        if S_maxlag=='no':
            fp.write('{:s}'.format(S_maxlag))
        else:
            fp.write('{:.4f}'.format(S_maxlag))
        fp.write('\t')
        fp.write('\n')

        ######################
        fp.write('Surf_maxlag')
        fp.write('\t')
        if Surf_maxlag=='no':
            fp.write('{:s}'.format(Surf_maxlag))
        else:
            fp.write('{:.4f}'.format(Surf_maxlag))
        fp.write('\t')
        fp.write('\n')

        ######################
        fp.write('P_filter')
        fp.write('\t')
        fp.write('\t')
        if P_filter[0]=='no':
            fp.write('{:s}'.format(P_filter[0]))
        else:
            fp.write('{:.3f}'.format(P_filter[0]))
        fp.write('\t')
        if P_filter[1]=='no':
            fp.write('{:s}'.format(P_filter[1]))
        else:
            fp.write('{:.3f}'.format(P_filter[1]))
        fp.write('\t')
        fp.write('\n')

        ######################
        fp.write('S_filter')
        fp.write('\t')
        fp.write('\t')
        if S_filter[0]=='no':
            fp.write('{:s}'.format(S_filter[0]))
        else:
            fp.write('{:.3f}'.format(S_filter[0]))
        fp.write('\t')
        if S_filter[1]=='no':
            fp.write('{:s}'.format(S_filter[1]))
        else:
            fp.write('{:.3f}'.format(S_filter[1]))
        fp.write('\t')
        fp.write('\n')

        ######################
        fp.write('Surf_filter')
        fp.write('\t')
        if Surf_filter[0]=='no':
            fp.write('{:s}'.format(Surf_filter[0]))
        else:
            fp.write('{:.3f}'.format(Surf_filter[0]))
        fp.write('\t')
        if Surf_filter[1]=='no':
            fp.write('{:s}'.format(Surf_filter[1]))
        else:
            fp.write('{:.3f}'.format(Surf_filter[1]))
        fp.write('\t')
        fp.write('\n')
        fp.write('\n')

    fp.close()














