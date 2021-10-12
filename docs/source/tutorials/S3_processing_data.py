#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import glob
import obspy
import shutil
from obspy import Stream
from obspy.taup import TauPyModel
from MCMTpy.utils.asdf_function import Add_waveforms,Add_stationxml,get_QuakeML,Add_quakeml
from MCMTpy.utils.asdf_function import get_StationXML

#----------------------------------------------------#
# 1.read raw data
def read_data(allfiles_path):

    allfiles = sorted( glob.glob(allfiles_path) )
    data_raw = Stream()
    for i in range(0,len(allfiles),1):
        try:
            tr = obspy.read(allfiles[i])
            data_raw += tr
        except Exception:
            print(allfiles[i],': no such file or obspy read error');continue
    data = data_raw.copy()
    
    return data

#----------------------------------------------------#
# 2.get_taup_tp_ts
def get_taup_tp_ts(model,depth,distance,degree=None):
    if degree==False:
        distance = distance/111.19

    time_p = model.get_travel_times(source_depth_in_km=depth,
                                    distance_in_degree=distance,
                                    phase_list=["p", "P"])

    time_s = model.get_travel_times(source_depth_in_km=depth,
                                    distance_in_degree=distance,
                                    phase_list=["s", "S"])

    ray_p = time_p[0].ray_param
    tp = time_p[0].time
    angle_p = time_p[0].incident_angle

    ray_s = time_s[0].ray_param
    ts = time_s[0].time
    angle_s = time_s[0].incident_angle

    return ray_p,tp,angle_p,ray_s,ts,angle_s

#----------------------------------------------------#
# 3.preprocess
def preprocess(data,freqmin,freqmax,samp_freq,p_n0,npts):

    data.detrend(type='demean')
    data.detrend(type='simple')
    data.taper(max_percentage=0.05)

    data.filter('bandpass', freqmin=freqmin, freqmax=freqmax, corners=4, zerophase=True)  # bandpass
    data.resample(samp_freq,window='hanning',no_filter=True, strict_length=False)         # resample

    for i in range(0,len(data),1):
        b = data[i].stats.sac['b']
        t1 = data[i].stats.sac['t1']
        t_start = data[i].stats.starttime + (t1-b) - p_n0/samp_freq
        t_end_npts = t_start+npts/samp_freq
        t_endtime = data[i].stats.endtime
        if t_end_npts > t_endtime:
            t_end = t_endtime
        else:
            t_end = t_end_npts

        data[i].trim(t_start, t_end)

    for i in range(0,round(len(data)/3)):
        t_start_1 = data[3*i].stats.starttime
        t_start_2 = data[3*i+1].stats.starttime
        t_start_3 = data[3*i+2].stats.starttime
        t_end_1 = data[3*i].stats.endtime
        t_end_2 = data[3*i+1].stats.endtime
        t_end_3 = data[3*i+2].stats.endtime
        t_start = max(t_start_1,t_start_2,t_start_3)
        t_end = min(t_end_1,t_end_2,t_end_3)
        data[3*i:3*i+3].trim(t_start, t_end)

    # velocity(nm/s) to velocity(cm/s)
    for i in range(0,len(data),1):
        data[i].data = 1e-6*data[i].data

#----------------------------------------------------#
# 4.write_ASDF
def write_ASDF(data_enz,Output_path,source_name,ASDF_filename):
    source_lat   = data_enz[0].stats.sac.evla
    source_lon   = data_enz[0].stats.sac.evlo
    source_depth = data_enz[0].stats.sac.evdp
    source_time  = data_enz[0].stats.starttime + data_enz[0].stats.sac.b

    Output_file_path = os.path.join(Output_path, ASDF_filename+'.h5')
    catalog_create   = get_QuakeML(source_name, source_lat, source_lon, source_depth,source_time)
    Add_quakeml(Output_file_path,catalog_create)

    Station_num = round(len(data_enz)/3)
    for i in range(0,Station_num,1): 
        net_sta_name     = data_enz[3*i].stats.network+'_'+data[3*i].stats.station
        station_network  = data_enz[3*i].stats.network
        station_name     = data_enz[3*i].stats.station
        station_lat      = data_enz[3*i].stats.sac.stla
        station_lon      = data_enz[3*i].stats.sac.stlo
        station_depth    = 0
        station_distance = data_enz[3*i].stats.sac.dist
        station_az       = data_enz[3*i].stats.sac.az
        station_baz      = data[3*i].stats.sac.baz
        tp               = data_enz[3*i].stats.sac.t1
        ts               = data_enz[3*i].stats.sac.t2
        stream           = data_enz[3*i:3*i+3].copy()
        starttime        = data_enz[3*i].stats.starttime
        time_before_first_arrival = p_n0/samp_freq

        inv=get_StationXML(station_network, 
                           station_name,
                           station_lat, 
                           station_lon, 
                           station_depth, 
                           'enz')

        Add_waveforms([stream], 
                      catalog_create,
                      Output_file_path,
                      source_name,
                      tp,
                      ts,
                      station_distance,
                      station_az,
                      station_baz) 

        Add_stationxml(inv,
                       Output_file_path) 

        for j in range(0,3,1): 
            CHN = stream[j].stats.channel
            sac_filename = os.path.join(Output_path,source_name+'.'+station_network+'.'+station_name+'.'+CHN+'.sac')
            stream[j].write(sac_filename, format='SAC')
#------------------------------------------------------------#




#------------------------------------------------------------#
# 0. main
# The root directory of the project
example_path = '/Users/yf/3.Project/8.MCMTpy/MCMTpy-master/data/example_yunnan'
freqmin       = 0.005                               # pre filtering frequency bandwidth (hz)
freqmax       = 2                                   # note this cannot exceed Nquist frequency (hz)
samp_freq     = 5                                   # targeted sampling rate (hz)
p_n0          = 50                                  # number of sampling points before P wave arrives
npts          = 2048                                # data length
Output_path   = os.path.join(example_path,"YN.202105212148_Inv/YN.202105212148_MCMTpy")
source_name   = "source_enz"                        # asdf‘s source_tag
ASDF_filename = "YN.202105212148"                   # asdf filename





#------------------------------------------------------------#
# we expect no parameters need to be changed below
#------------------------------------------------------------#
# 1. read raw data
allfiles_path = os.path.join(example_path,'YN.202105212148_Inv/YN.202105212148_raw/*.SAC')  
data = read_data(allfiles_path)

# 2. taup ray tracing
model_path = os.path.join(example_path,"v_model/v_model.npz")
model = TauPyModel(model=model_path)
for i in range(0,len(data),1):
    depth = data[i].stats.sac['evdp']
    distance = data[i].stats.sac['dist']
    ray_p,tp,angle_p,ray_s,ts,angle_s = get_taup_tp_ts(model,depth,distance,degree=False)
    data[i].stats.sac["t1"]=tp
    data[i].stats.sac["t2"]=ts

# 3. preprocess
preprocess(data,freqmin,freqmax,samp_freq,p_n0,npts)

# 4. output ASDF and SAC data
shutil.rmtree(Output_path)
os.mkdir(Output_path)
write_ASDF(data,Output_path,source_name,ASDF_filename)