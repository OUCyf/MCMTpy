#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  1 11:51:12 2021

This script:
    1) ...


Modify history:
    1) Apr  1 11:51:12 2021    ||    Fu Yin at USTC    ||    The initial release.
    2) ...
    
"""


import os
import math
import pyasdf
import numpy as np
from pyfk import SourceModel, SeisModel, Config
from pyfk import calculate_sync
from MCMTpy.utils.distaz import DistAz
from MCMTpy.utils.rotate import rotate





#%% subfunction
############################################
def write_inv_para(file_path,Inv_para,Sta_inv_info,Sta_data,Stf,file_path_source,source_info,CPU_splits):
    fout = open(file_path, 'w')

    # 1.Inv_para
    fout.write('Inv_para:\n')
    for k,v in Inv_para.items():
        fout.write('\t')
        fout.write(str(k)+': '+str(v)+'\n')
    fout.write('\n\n')

    # 2.Sta_inv_info
    fout.write('Sta_inv_info:\n')
    for k1,v1 in Sta_inv_info.items():
        fout.write('\t')
        fout.write(str(k1)+': '+'\n')
        for k2,v2 in Sta_inv_info[k1].items():
            fout.write('\t\t')
            fout.write(str(k2)+': '+str(v2)+'\n')
        fout.write('  \n')
    fout.write('\n\n')

    # 3.Sta_data
    fout.write('Sta_data:\n')
    for k1,v1 in Sta_data.items():
        fout.write('\t')
        fout.write(str(k1)+': '+'\n')
        for k2,v2 in Sta_data[k1].items():
            fout.write('\t\t')
            fout.write(str(k2)+': '+str(v2)+'\n')
        fout.write('  \n')
    fout.write('\n\n')

    # 4.Stf
    fout.write('Stf:\n')
    fout.write('\t')
    fout.write(str(Stf))
    fout.write('\n\n')

    fout.close()


    fout = open(file_path_source, 'w') 
    fout.write('The CPU core mission:\n\n')
    for i in range(0,len(CPU_splits),1):
        fout.write('core_'+str(i)+':  '+str(CPU_splits[i]) )
        fout.write('\t')
        fout.write('\n')
    fout.write('\n**************************************\n')

    fout.write('The source_info format MCMTpy inverted:\n')
    fout.write('number\t[depth]\t[strike]\t[dip]\t[rake]\n\n')
    for i in range(0,len(source_info),1):
        fout.write( str(i)+'\t'+str(source_info[i]) )
        fout.write('\t')
        fout.write('\n')
    fout.close()


############################################
def config_fk(InvType):
    '''
    Initialization pyfkmodel
    '''
    if InvType=='mt' or InvType=='dc':
        srcType='dc'
    else:
        srcType=InvType

    source_prem = SourceModel(
        sdep=1.0111,                                                            # The process of synthesizing the waveform is independent of this parameter
        srcType=srcType)

    v_model = SeisModel(model=np.array( [[0.05, 0.3, 0.6, 1.8, 98, 50],         # The process of synthesizing the waveform is independent of this parameter
                                      [0.05, 0.3, 0.6, 1.8, 98, 50],
                                      [0.05, 0.3, 0.6, 1.8, 98, 50],
                                      [0.05, 0.3, 0.6, 1.8, 98, 50]] ))

    config_prem = Config(
        model=v_model,
        source=source_prem,
        receiver_distance=np.array([1]) )                                       # The process of synthesizing the waveform is independent of this parameter


    return source_prem,v_model,config_prem








############################################
def get_GFs_sta_info(point,NET_STA,prepro_para,Source_Station_info_MPI):
    
    if prepro_para["Database_mode"]==True:
        # 1.gets information about all the stations for a given single earthquake
        source_depth_min=prepro_para["Source_depth_min"]
        source_depth_spacing=prepro_para["Source_depth_spacing"]
        source_distance_radius_min=prepro_para["Station_distance_radius_min"]
        source_distance_spacing=prepro_para["Station_distance_spacing"]
        
        GFs_sta_info={}                                                             # reorganize station information into a dictionary
        for i in range(0,len(NET_STA),1):
            source_lat = point[0]
            source_lon = point[1]
            source_depth = point[2]
            station_lat = NET_STA[i][2]
            station_lon = NET_STA[i][3]
            
            geodist = DistAz(station_lat, station_lon, source_lat, source_lon)
            if prepro_para["degrees"]==True:
                source_distance = geodist.getDelta()                                # get distance
            else:
                source_distance = geodist.getDelta() * 111.19                       # get distance km
            source_az = geodist.getAz()                                             # get azimuth
            source_baz = geodist.getBaz()

            id_dep_a = (source_depth-source_depth_min)//source_depth_spacing
            id_dep_b = (source_depth-source_depth_min)//source_depth_spacing + 1
            depth_a = source_depth_min+id_dep_a*source_depth_spacing
            depth_b = source_depth_min+id_dep_b*source_depth_spacing
            if abs(source_depth-depth_a) <= abs(source_depth-depth_b):
                id_depth_gfs = id_dep_a
                depth_gfs = depth_a
            else:
                id_depth_gfs = id_dep_b
                depth_gfs = depth_b

            id_dis_a = (source_distance-source_distance_radius_min)//source_distance_spacing
            id_dis_b = (source_distance-source_distance_radius_min)//source_distance_spacing + 1
            distance_a = source_distance_radius_min+id_dis_a*source_distance_spacing
            distance_b = source_distance_radius_min+id_dis_b*source_distance_spacing
            if abs(source_distance-distance_a) <= abs(source_distance-distance_b):
                id_distance_gfs = id_dis_a
                distance_gfs = distance_a
            else:
                id_distance_gfs = id_dis_b
                distance_gfs = distance_b

            source_name = prepro_para["Source_name"]+'_'+str(round(id_depth_gfs))
            station_name = prepro_para["Station_name"]+'_'+str(round(id_distance_gfs))
            network_name = prepro_para["Network_name"]
            station_depth_reference = prepro_para["Station_depth_reference"]
            
            ff=os.path.join(prepro_para["DATADIR_splits"],source_name)              # folder name of all Green's functions for a depth
            file_path = os.path.join( ff,source_name+'_'+station_name+'.h5')        # green's function for a deep earthquake at one station
            
            GFs_net_sta_name = NET_STA[i][0]+'.'+NET_STA[i][1]
            one_station_info = [file_path,
                                source_name,depth_gfs,  
                                network_name,station_name,distance_gfs,station_depth_reference,source_az,source_baz ]

            GFs_sta_info.update( {GFs_net_sta_name: one_station_info} ) 


    else:
        Source_num = prepro_para["Source_num"]                                  #  the number of earthquakes in Green's function library
        
        dis = np.zeros(Source_num)                                              # The error between the position and point position of all earthquakes in Green's function library
        for i in range(0,Source_num, 1):
            dis[i] = (point[0]-Source_Station_info_MPI[i][0][1])*111.19**2+(point[1]-Source_Station_info_MPI[i][0][2])*111.19**2\
                +(point[2]-Source_Station_info_MPI[i][0][3])**2
        
        min_tag_index = np.argmin(dis)                                          # the position of the minimum error
        source_name = Source_Station_info_MPI[min_tag_index][0][0]              # minimum error corresponding to source_name
        depth_gfs = Source_Station_info_MPI[min_tag_index][0][3]

        # All station information for a given earthquake (why read all? Because we can't get the specified station name from Source_Station_info_MPI.)
        GFs_sta_info_all={}
        Sta_num_all = len(Source_Station_info_MPI[min_tag_index])-1
        for i in range(0,Sta_num_all,1):
            net_name = Source_Station_info_MPI[min_tag_index][i+1][0]
            sta_name = Source_Station_info_MPI[min_tag_index][i+1][1]
            GFs_net_sta_name = net_name+'.'+sta_name
            
            file_path = os.path.join( prepro_para["DATADIR_splits"],source_name+'.h5')
            network_name = Source_Station_info_MPI[min_tag_index][i+1][0]
            station_name = Source_Station_info_MPI[min_tag_index][i+1][1]
            distance_gfs = Source_Station_info_MPI[min_tag_index][i+1][5]
            station_depth = Source_Station_info_MPI[min_tag_index][i+1][4]
            source_az = Source_Station_info_MPI[min_tag_index][i+1][6]
            source_baz = Source_Station_info_MPI[min_tag_index][i+1][7]
            
            one_station_info = [file_path,
                                source_name,depth_gfs,  
                                network_name,station_name,distance_gfs,station_depth,source_az,source_baz ]
            
            GFs_sta_info_all.update( {GFs_net_sta_name: one_station_info} ) 

        # reorganize the information of the stations used to synthesize the waveform into a dictionary
        GFs_sta_info={}
        for i in range(0,len(NET_STA),1):
            GFs_net_sta_name = NET_STA[i][0]+'.'+NET_STA[i][1]
            one_station_info = GFs_sta_info_all[GFs_net_sta_name]
            GFs_sta_info.update( {GFs_net_sta_name: one_station_info} ) 

    return GFs_sta_info




############################################
def get_Sta_rand_data(GFs_sta_info,config_prem,source_time_function):
    
    Sta_rand_data={}
    net_sta_name = list(GFs_sta_info.keys())
    for i in range(0,len(GFs_sta_info),1):
        info={}
        file_path = GFs_sta_info[net_sta_name[i]][0]
        gfs_path = file_path
        ds_gfs=pyasdf.ASDFDataSet(gfs_path,mpi=False,mode='r')
        Source_tag = GFs_sta_info[net_sta_name[i]][1]                           # find the Source_tag of the Green's function for the first station
        net_sta_tag = GFs_sta_info[net_sta_name[i]][3]+'.'+GFs_sta_info[net_sta_name[i]][4]  # Find the net_sta_tag of Green's function for the first station
        gf = ds_gfs.waveforms[net_sta_tag][Source_tag]
        az = GFs_sta_info[net_sta_name[i]][7]
        baz = GFs_sta_info[net_sta_name[i]][8]
        sync_result = calculate_sync(gf, config_prem, az, source_time_function)  # get syn result data in order of Z R T component
        for chn in range(0,len(sync_result[0])):
            sync_result[0][chn].stats.sac.baz = baz

        tpp = float( ds_gfs.waveforms[net_sta_tag][Source_tag][0].stats.asdf['labels'][1] )
        tss = float( ds_gfs.waveforms[net_sta_tag][Source_tag][0].stats.asdf['labels'][3] )

        info.update({ "tp": tpp})  
        info.update({ "ts": tss})
        info.update({ "data": sync_result[0]})
        Sta_rand_data.update({ net_sta_name[i]: info})


    return Sta_rand_data




############################################
def correlate_maxlag(np1,np2,maxlag,normalized=True):
    
    if len(np1)!=len(np2):
        raise ValueError('Function: correlate_maxlag error of len(np1)!=len(np2)\
One possible reason for the error is that the P_win or S_win or Surf_win \
length of the data is exceeded')
    
    if maxlag < 0:
        maxlag = 0
    if maxlag > len(np1)-1:
        maxlag = len(np1)-1
    
    if normalized==True:
        CC = np.correlate(np1, np2,mode='full')/math.sqrt( np.sum(np1**2) * np.sum(np2**2) )   # -1<= max(CC) <=1
    else:
        CC = np.correlate(np1, np2,mode='full')

    length_CC = CC.shape[0]
    CC_maxlag = CC[ math.floor(length_CC/2)-maxlag:math.floor(length_CC/2)+maxlag+1]
    C_max_maxlag = np.max(CC_maxlag)
    C_max_maxlag_index = np.argmax(CC_maxlag)
    shift = C_max_maxlag_index-maxlag                                           # length_CC_maxlag || The distance the maximum cross - correlation moves when aligning the waveform


    return CC_maxlag,C_max_maxlag,shift




############################################
def zero_padding(np2,shift):
    '''
    It's used in get_MISFIT_1, but not used in plotting.
    '''
    zeros = np.zeros( shape=(abs(shift)) )
    if shift==0:
        np2_add_zeros = np2
    if shift>0:  # np2 || move it back，In front of zero padding
        np2_add_zeros = np.hstack( (zeros,np2[:-shift]) )
    if shift<0:  # np2 || move it front，In back of zero padding
        np2_add_zeros = np.hstack( (np2[-shift:],zeros) )


    return np2_add_zeros







############################################
# T0 is a randomly generated variation of the shock generating time
def get_MISFIT_1(GFs_sta_info,Sta_raw_data,Sta_rand_data,Sta_inv_info,Inv_para,Prepro_para,T0):
    dt = Inv_para["Dt"]
    Distance_0= Inv_para["Distance_0"]
    Raw_p_n0=Inv_para["Raw_p_n0"]
    Rand_p_n0=Prepro_para["samples_before_first_arrival"]

    MISFIT1=0
    M0=[]
    net_sta_name = list(Sta_inv_info.keys())
    for m in range(0,Inv_para["Station_num"],1):
        distance = GFs_sta_info[net_sta_name[m]][5]
        baz = GFs_sta_info[net_sta_name[m]][8]
        Raw_tp = Sta_raw_data[net_sta_name[m]]['tp']
        Raw_ts = Sta_raw_data[net_sta_name[m]]['ts']
        Raw_tsurf = Raw_ts/0.913
        Rand_tp = Sta_rand_data[net_sta_name[m]]['tp']
        Rand_ts = Sta_rand_data[net_sta_name[m]]['ts']
        Rand_tsurf = Rand_ts/0.913
        
        P_win = Sta_inv_info[net_sta_name[m]]['P_win']
        S_win = Sta_inv_info[net_sta_name[m]]['S_win']
        Surf_win = Sta_inv_info[net_sta_name[m]]['Surf_win']
        
        P_inv_comp = Sta_inv_info[net_sta_name[m]]['P_inv_comp']
        S_inv_comp = Sta_inv_info[net_sta_name[m]]['S_inv_comp']
        Surf_inv_comp = Sta_inv_info[net_sta_name[m]]['Surf_inv_comp']
        
        Phase_weight = Sta_inv_info[net_sta_name[m]]['Phase_weight']
        Distance_weight = Sta_inv_info[net_sta_name[m]]['Distance_weight'][0]
        
        P_filter = Sta_inv_info[net_sta_name[m]]['P_filter']
        S_filter = Sta_inv_info[net_sta_name[m]]['S_filter']
        Surf_filter = Sta_inv_info[net_sta_name[m]]['Surf_filter']

        P_maxlag = Sta_inv_info[net_sta_name[m]]['P_maxlag'][0]
        S_maxlag = Sta_inv_info[net_sta_name[m]]['S_maxlag'][0]
        Surf_maxlag = Sta_inv_info[net_sta_name[m]]['Surf_maxlag'][0]
        P_maxlag = round(P_maxlag/dt)
        S_maxlag = round(S_maxlag/dt)
        Surf_maxlag = round(Surf_maxlag/dt)

        Raw_data = Sta_raw_data[net_sta_name[m]]['data']                        # 3 Trace(s) in Stream, order:ENZ
        Rand_data = Sta_rand_data[net_sta_name[m]]['data']                      # order:ZRT
        
        ### rotate 'ENZ->ZRT'
        for chn in range(0,3):
            Raw_data[chn].stats.back_azimuth = float( baz )
        Raw_data = rotate(data_raw=Raw_data, chn=['ENZ','ZRT'])


        #####
        P_win_length = round((P_win[1]-P_win[0])/dt)
        S_win_length = round((S_win[1]-S_win[0])/dt)
        Surf_win_length = round((Surf_win[1]-Surf_win[0])/dt)
        
        raw_P_start = Raw_p_n0 + round(P_win[0]/dt)
        raw_S_start = Raw_p_n0 + round( (Raw_ts+S_win[0]-Raw_tp)/dt )
        raw_Surf_start = Raw_p_n0 + round( (Raw_tsurf+Surf_win[0]-Raw_tp)/dt )
        
        # The estimated seismogenesis time of the original data is used to intercept the time window of Green's function library
        rand_P_start = Rand_p_n0 + round(P_win[0]/dt) + round((Raw_tp-T0-Rand_tp)/dt)  # 用原始数据的估计的发震时间来截取格林函数库的时间窗
        rand_S_start = Rand_p_n0 + round( (Rand_ts+S_win[0]-Rand_tp)/dt ) + round((Raw_ts-T0-Rand_ts)/dt)
        rand_Surf_start = Rand_p_n0 + round( (Rand_tsurf+Surf_win[0]-Rand_tp)/dt ) + round((Raw_tsurf-T0-Rand_tsurf)/dt)
        
        # rand_P_start = Rand_p_n0 + round(P_win[0]/dt)
        # rand_S_start = Rand_p_n0 + round( (Rand_ts+S_win[0]-Rand_tp)/dt )
        # rand_Surf_start = Rand_p_n0 + round( (Rand_tsurf+Surf_win[0]-Rand_tp)/dt )
        
        if (raw_P_start<0) or (raw_S_start<0) or (raw_Surf_start<0) or (rand_P_start<0) or (rand_S_start<0) or (rand_Surf_start<0):
            raise ValueError(f'{net_sta_name[m]} P_win or S_win or Surf_win exceeded the start time of the data with rand_tp: {Rand_tp} raw_tp: {Raw_tp} rand_ts: {Rand_ts} raw_ts: {Raw_ts} rand_tsurf: {Rand_tsurf} raw_tsurf: {Raw_tsurf}')
        
        Raw_P_win_n0 = [raw_P_start, raw_P_start + P_win_length]
        Raw_S_win_n0 = [raw_S_start, raw_S_start + S_win_length]
        Raw_Surf_win_n0 = [raw_Surf_start, raw_Surf_start + Surf_win_length]
        
        Rand_P_win_n0 = [rand_P_start, rand_P_start + P_win_length]
        Rand_S_win_n0 = [rand_S_start, rand_S_start + S_win_length]
        Rand_Surf_win_n0 = [rand_Surf_start, rand_Surf_start + Surf_win_length]

        MISFIT1_one_sta=0
        M0_one_sta=[]
        if Inv_para['Geometrical_spreading_mode']==True:
            #### 1. P wave misfit
            Raw_data_filter = Raw_data.copy()
            Raw_data_filter.filter('bandpass', freqmin=P_filter[0], freqmax=P_filter[1], corners=4, zerophase=True)
            Rand_data_filter = Rand_data.copy()
            Rand_data_filter.filter('bandpass', freqmin=P_filter[0], freqmax=P_filter[1], corners=4, zerophase=True)
            for i in range (0,3,1):                                             # p wave in ZRT comp
                if P_inv_comp[i]=='yes':
                    Raw_data_np = Raw_data_filter[i].data
                    Rand_data_np = Rand_data_filter[i].data
                    
                    Raw_data_np_p = Raw_data_np[ Raw_P_win_n0[0]:Raw_P_win_n0[1] ] 
                    Rand_data_np_p = Rand_data_np[ Rand_P_win_n0[0]:Rand_P_win_n0[1] ] 
                    
                    _, C_max_maxlag, shift=correlate_maxlag(Raw_data_np_p,Rand_data_np_p,P_maxlag,normalized=True)
                    Rand_data_np_p_shift = zero_padding(Rand_data_np_p,shift)
                    # Rand_data_np_p_shift = Rand_data_np[ Rand_P_win_n0[0]-shift:Rand_P_win_n0[1]-shift] 
                    
                    MISFIT1_one_sta += Phase_weight[0] * (distance/Distance_0)**2*math.sqrt(np.mean((Raw_data_np_p-Rand_data_np_p_shift)**2)) # * (1-C_max_maxlag)
                    M0_one_sta.append( np.linalg.norm(Raw_data_np_p) / np.linalg.norm(Rand_data_np_p_shift) )
                
            #### 2. S wave misfit
            Raw_data_filter = Raw_data.copy()
            Raw_data_filter.filter('bandpass', freqmin=S_filter[0], freqmax=S_filter[1], corners=4, zerophase=True)
            Rand_data_filter = Rand_data.copy()
            Rand_data_filter.filter('bandpass', freqmin=S_filter[0], freqmax=S_filter[1], corners=4, zerophase=True)
            for i in range (0,3,1):                                             # S wave in ZRT comp
                if S_inv_comp[i]=='yes':
                    Raw_data_np = Raw_data_filter[i].data
                    Rand_data_np = Rand_data_filter[i].data
                    
                    Raw_data_np_s = Raw_data_np[ Raw_S_win_n0[0]:Raw_S_win_n0[1] ]
                    Rand_data_np_s = Rand_data_np[ Rand_S_win_n0[0]:Rand_S_win_n0[1] ] 
                    
                    _, C_max_maxlag, shift=correlate_maxlag(Raw_data_np_s,Rand_data_np_s,S_maxlag,normalized=True)
                    Rand_data_np_s_shift = zero_padding(Rand_data_np_s,shift)
                    # Rand_data_np_s_shift = Rand_data_np[ Rand_S_win_n0[0]-shift:Rand_S_win_n0[1]-shift ] 
                    
                    MISFIT1_one_sta += Phase_weight[1] * (distance/Distance_0)**2*math.sqrt(np.mean((Raw_data_np_s-Rand_data_np_s_shift)**2)) #  * (1-C_max_maxlag)
                    M0_one_sta.append( np.linalg.norm(Raw_data_np_s) / np.linalg.norm(Rand_data_np_s_shift) )
                
            #### 3. Surf wave misfit
            Raw_data_filter = Raw_data.copy()
            Raw_data_filter.filter('bandpass', freqmin=Surf_filter[0], freqmax=Surf_filter[1], corners=4, zerophase=True)
            Rand_data_filter = Rand_data.copy()
            Rand_data_filter.filter('bandpass', freqmin=Surf_filter[0], freqmax=Surf_filter[1], corners=4, zerophase=True)
            for i in range (0,3,1):                                             # Surf wave in ZRT comp
                if Surf_inv_comp[i]=='yes':
                    Raw_data_np = Raw_data_filter[i].data
                    Rand_data_np = Rand_data_filter[i].data
                    
                    Raw_data_np_surf = Raw_data_np[ Raw_Surf_win_n0[0]:Raw_Surf_win_n0[1] ]  
                    Rand_data_np_surf = Rand_data_np[ Rand_Surf_win_n0[0]:Rand_Surf_win_n0[1] ] 
                    
                    _, C_max_maxlag, shift=correlate_maxlag(Raw_data_np_surf,Rand_data_np_surf,Surf_maxlag,normalized=True)
                    Rand_data_np_surf_shift = zero_padding(Rand_data_np_surf,shift)
                    # Rand_data_np_surf_shift = Rand_data_np[ Rand_Surf_win_n0[0]-shift:Rand_Surf_win_n0[1]-shift ] 
                    
                    MISFIT1_one_sta += Phase_weight[2] * (distance/Distance_0)**1*math.sqrt(np.mean((Raw_data_np_surf-Rand_data_np_surf_shift)**2)) #  * (1-C_max_maxlag)
                    M0_one_sta.append( np.linalg.norm(Raw_data_np_surf) / np.linalg.norm(Rand_data_np_surf_shift) )
        
        else:
            #### 1. P wave misfit
            Raw_data_filter = Raw_data.copy()
            Raw_data_filter.filter('bandpass', freqmin=P_filter[0], freqmax=P_filter[1], corners=4, zerophase=True)
            Rand_data_filter = Rand_data.copy()
            Rand_data_filter.filter('bandpass', freqmin=P_filter[0], freqmax=P_filter[1], corners=4, zerophase=True)
            for i in range (0,3,1):                                             # p wave in ZRT comp
                if P_inv_comp[i]=='yes':
                    Raw_data_np = Raw_data_filter[i].data
                    Rand_data_np = Rand_data_filter[i].data
                    
                    Raw_data_np_p = Raw_data_np[ Raw_P_win_n0[0]:Raw_P_win_n0[1] ] / math.sqrt( np.sum( Raw_data_np[ Raw_P_win_n0[0]:Raw_P_win_n0[1] ]**2, axis=0) )
                    Rand_data_np_p = Rand_data_np[ Rand_P_win_n0[0]:Rand_P_win_n0[1] ] / math.sqrt( np.sum( Rand_data_np[ Rand_P_win_n0[0]:Rand_P_win_n0[1] ]**2, axis=0) )
                    
                    _, C_max_maxlag, shift=correlate_maxlag(Raw_data_np_p,Rand_data_np_p,P_maxlag,normalized=True)
                    np2_zeros = zero_padding(Rand_data_np_p,shift)
                    Rand_data_np_p_shift = np2_zeros / math.sqrt( np.sum(np2_zeros**2, axis=0) )
                    # Rand_data_np_p_shift = Rand_data_np[ Rand_P_win_n0[0]-shift:Rand_P_win_n0[1]-shift ] / math.sqrt( np.sum( Rand_data_np[ Rand_P_win_n0[0]-shift:Rand_P_win_n0[1]-shift ]**2, axis=0) )
                    
                    MISFIT1_one_sta += Phase_weight[0]*math.sqrt(np.mean((Raw_data_np_p-Rand_data_np_p_shift)**2)) # *(1-C_max_maxlag)
                    M0_one_sta.append( np.linalg.norm(Raw_data_np_p) / np.linalg.norm(Rand_data_np_p_shift) )
                
            #### 2. S wave misfit
            Raw_data_filter = Raw_data.copy()
            Raw_data_filter.filter('bandpass', freqmin=S_filter[0], freqmax=S_filter[1], corners=4, zerophase=True)
            Rand_data_filter = Rand_data.copy()
            Rand_data_filter.filter('bandpass', freqmin=S_filter[0], freqmax=S_filter[1], corners=4, zerophase=True)
            for i in range (0,3,1):                                             # S  wave in ZRT comp
                if S_inv_comp[i]=='yes':
                    Raw_data_np = Raw_data_filter[i].data
                    Rand_data_np = Rand_data_filter[i].data
                    
                    Raw_data_np_s = Raw_data_np[ Raw_S_win_n0[0]:Raw_S_win_n0[1] ] / math.sqrt( np.sum( Raw_data_np[ Raw_S_win_n0[0]:Raw_S_win_n0[1] ]**2, axis=0) )
                    Rand_data_np_s = Rand_data_np[ Rand_S_win_n0[0]:Rand_S_win_n0[1] ] / math.sqrt( np.sum( Rand_data_np[ Rand_S_win_n0[0]:Rand_S_win_n0[1] ]**2, axis=0) )
                    
                    _, C_max_maxlag, shift=correlate_maxlag(Raw_data_np_s,Rand_data_np_s,S_maxlag,normalized=True)
                    np2_zeros = zero_padding(Rand_data_np_s,shift)
                    Rand_data_np_s_shift =  np2_zeros / math.sqrt( np.sum( np2_zeros**2, axis=0) )
                    # Rand_data_np_s_shift = Rand_data_np[ Rand_S_win_n0[0]-shift:Rand_S_win_n0[1]-shift ] / math.sqrt( np.sum( Rand_data_np[ Rand_S_win_n0[0]-shift:Rand_S_win_n0[1]-shift ]**2, axis=0) )
                    
                    MISFIT1_one_sta += Phase_weight[1]*math.sqrt(np.mean((Raw_data_np_s-Rand_data_np_s_shift)**2)) # *(1-C_max_maxlag)
                    M0_one_sta.append( np.linalg.norm(Raw_data_np_s) / np.linalg.norm(Rand_data_np_s_shift) )
                
            #### 3. Surf wave misfit
            Raw_data_filter = Raw_data.copy()
            Raw_data_filter.filter('bandpass', freqmin=Surf_filter[0], freqmax=Surf_filter[1], corners=4, zerophase=True)
            Rand_data_filter = Rand_data.copy()
            Rand_data_filter.filter('bandpass', freqmin=Surf_filter[0], freqmax=Surf_filter[1], corners=4, zerophase=True)
            for i in range (0,3,1):                                             # Surf  wave in ZRT comp
                if Surf_inv_comp[i]=='yes':
                    Raw_data_np = Raw_data_filter[i].data
                    Rand_data_np = Rand_data_filter[i].data
                    
                    Raw_data_np_surf = Raw_data_np[ Raw_Surf_win_n0[0]:Raw_Surf_win_n0[1] ] / math.sqrt( np.sum( Raw_data_np[ Raw_Surf_win_n0[0]:Raw_Surf_win_n0[1] ]**2, axis=0) )
                    Rand_data_np_surf = Rand_data_np[ Rand_Surf_win_n0[0]:Rand_Surf_win_n0[1] ] / math.sqrt( np.sum( Rand_data_np[ Rand_Surf_win_n0[0]:Rand_Surf_win_n0[1] ]**2, axis=0) )
                    
                    _, C_max_maxlag, shift=correlate_maxlag(Raw_data_np_surf,Rand_data_np_surf,Surf_maxlag,normalized=True)
                    np2_zeros = zero_padding(Rand_data_np_surf,shift)
                    Rand_data_np_surf_shift =  np2_zeros / math.sqrt( np.sum( np2_zeros**2, axis=0) )
                    # Rand_data_np_surf_shift = Rand_data_np[ Rand_Surf_win_n0[0]-shift:Rand_Surf_win_n0[1]-shift ] / math.sqrt( np.sum( Rand_data_np[ Rand_Surf_win_n0[0]-shift:Rand_Surf_win_n0[1]-shift ]**2, axis=0) )
                    
                    MISFIT1_one_sta += Phase_weight[2]*math.sqrt(np.mean((Raw_data_np_surf-Rand_data_np_surf_shift)**2)) # *(1-C_max_maxlag)
                    M0_one_sta.append( np.linalg.norm(Raw_data_np_surf) / np.linalg.norm(Rand_data_np_surf_shift) )
                
        MISFIT1 += MISFIT1_one_sta*Distance_weight
        M0.append(np.mean(np.array(M0_one_sta)))

    M0_mean = np.mean(np.array(M0))
    Mw_mean = math.log10(M0_mean)/1.5

    return MISFIT1,Mw_mean











