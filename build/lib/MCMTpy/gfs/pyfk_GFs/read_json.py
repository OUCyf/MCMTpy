#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 25 16:16:07 2021

@author: Fu Yin (yinfu@mail.ustc.edu.cn) at USTC

This script:
    1) function read_GFs_json is preprocess of function build_GFs_pyfk.


Modify history:
    1) Mar 25 16:16:07 2021    ||    Fu Yin at USTC    ||    The initial release.
    2) ...
    
"""

import json5
import numpy as np
from MCMTpy.utils.distaz import DistAz


#%%########################################################################
#                           -------------------
#                            1. main-function
#                           -------------------
###########################################################################

def read_GFs_json(filename):
    
    Source_depth_num=None
    Station_distance_num=None
    GFs_num=None
    Source_num=None
    GFs_info=None
    Source_Station_info_MPI=None


    #------------ 1.1. read json parameters ------------#
    with open(filename, 'r',encoding='utf-8') as f1:
        prepro_para_r=json5.load(f1)
    prepro_para = prepro_para_r.copy()                                          # copy()


    #------------ 1.2. read velocity model file ------------#
    Velocity_model = read_Velocity_model(prepro_para["Velocity_model"])
    filter_para = tuple(prepro_para["filter"])                                  # list ——> tuple for pyfk


    #------------ 1.3. compute/read source station info ------------#
    if prepro_para["Database_mode"] == True:
        # a.Databse mode  ——> GFs_num, Source_depth_num, Station_distance_num and GFs_info
        Source_depth_num = round( (prepro_para["Source_depth_max"] 
                                   - prepro_para["Source_depth_min"] )
                                 / prepro_para["Source_depth_spacing"] ) + 1
        
        Station_distance_num = round( (prepro_para["Station_distance_radius_max"] 
                                       - prepro_para["Station_distance_radius_min"] ) 
                                     / prepro_para["Station_distance_spacing"] ) + 1
        
        GFs_num = Source_depth_num*Station_distance_num
        
        GFs_info = []                                                           # storage all gfs info
        now_source = 0
        for i in range(0,Source_depth_num,1):
            for j in range(0,Station_distance_num,1):
                depth = prepro_para["Source_depth_min"] \
                        + i * prepro_para["Source_depth_spacing"]
                distance = prepro_para["Station_distance_radius_min"] \
                            + j * prepro_para["Station_distance_spacing"]
                
                Source_name_now = prepro_para["Source_name"]+'_'+str(i)
                Station_name_now = prepro_para["Station_name"]+'_'+str(j)
                one_gfs_info = [now_source, Source_name_now, depth, prepro_para["Network_name"], 
                                Station_name_now, distance,prepro_para["Station_depth_reference"]]
                GFs_info.append(one_gfs_info)
                now_source+=1 

    else:
        # b.single mode  ——> Source_num, Source_Station_info_MPI (with azimuth info)
        Source_Station_info = read_Source_Station_info(prepro_para["Source_Station_info"])
        Source_num = len(Source_Station_info)
        Source_Station_info_MPI=[]                                              # storage all gfs info
        for i in range(0, Source_num, 1):
            one_source_station_info=[Source_Station_info[i][0]]
            one_station_num = len(Source_Station_info[i])-1
            for j in range(0, one_station_num, 1):
                source_lat = Source_Station_info[i][0][1]
                source_lon = Source_Station_info[i][0][2]
                station_lat = Source_Station_info[i][j+1][2]
                station_lon = Source_Station_info[i][j+1][3]
                geodist = DistAz(station_lat, station_lon, source_lat, source_lon)
                
                dist = geodist.getDelta()                                       # distance
                az = geodist.getAz()                                            # azimuth
                baz = geodist.getBaz()                                          # baz
                
                one_sta_info = Source_Station_info[i][j+1]+[dist,az,baz]
                one_source_station_info.append(one_sta_info)
                
            Source_Station_info_MPI.append(one_source_station_info)


    #------------ 1.4. distribute earthquakes to each core of the CPU ------------#
    if prepro_para["Database_mode"] == True:
        CPU_splits = CPU_split(GFs_num,prepro_para["MPI_n"])
    else:
        CPU_splits = CPU_split(Source_num,prepro_para["MPI_n"])



    #------------ 1.5. sum all parameter info ——> Prepro_para（dict） ------------#
    Prepro_para = {
        # a.1. source station
        "Database_mode":prepro_para["Database_mode"],
        "Source_depth_min":prepro_para["Source_depth_min"],
        "Source_depth_max":prepro_para["Source_depth_max"],
        "Source_depth_spacing":prepro_para["Source_depth_spacing"],
        "Station_distance_radius_min":prepro_para["Station_distance_radius_min"],
        "Station_distance_radius_max":prepro_para["Station_distance_radius_max"],
        "Station_distance_spacing":prepro_para["Station_distance_spacing"],

        "Source_name": prepro_para["Source_name"],
        "Network_name": prepro_para["Network_name"],
        "Station_name": prepro_para["Station_name"],
        "Station_depth_reference":prepro_para["Station_depth_reference"],
        
        "Source_Station_info": prepro_para["Source_Station_info"],

        # a.2. MPI
        # MPI number <= source number
        "MPI_n":prepro_para["MPI_n"],
        
        # a.3. path
        "DATADIR": prepro_para["DATADIR"],
        "DATADIR_splits":prepro_para["DATADIR_splits"] ,

        # a.4. SourceModel,receiver-location
        "source_mechanism":prepro_para["source_mechanism"],
        "srcType":prepro_para["srcType"], 

        # a.5. SeisModel
        "Velocity_model": prepro_para["Velocity_model"],
        "flattening":prepro_para["flattening"],

        # a.6. Config
        "npt": prepro_para["npt"],
        "dt": prepro_para["dt"],
        "degrees": prepro_para["degrees"],
        "taper": prepro_para["taper"], 
        "filter": filter_para,
        "dk": prepro_para["dk"],
        "smth": prepro_para["smth"], 
        "pmin": prepro_para["pmin"], 
        "pmax": prepro_para["pmax"], 
        "kmax": prepro_para["kmax"], 
        "rdep": prepro_para["rdep"], 
        "updn": prepro_para["updn"], 
        "samples_before_first_arrival": prepro_para["samples_before_first_arrival"],

        # a.7. sum station/source parameter
        "Velocity_model_value": Velocity_model,
        "CPU_splits":CPU_splits,

        # Database_mode==true
        "Source_depth_num":Source_depth_num,
        "Station_distance_num":Station_distance_num,
        "GFs_num":GFs_num,

        # Database_mode==false
        "Source_num":Source_num,
    }


    return Prepro_para,GFs_info,Source_Station_info_MPI




#%%########################################################################
#                           -------------------
#                            2. sub-function
#                           -------------------
###########################################################################

def write_Prepro_para(file_path,prepro_para):
    fout = open(file_path, 'w') 
    fout.write('The Prepro_para MCMTpy used:\n\n')
    
    for k,v in prepro_para.items():
        if k=='Velocity_model_value':
            fout.write(str(k)+': '+'\n')
            for i in range(0,len(v),1):
                fout.write('\t')
                fout.write(str(v[i]))
                fout.write('\n')
        else:
            fout.write(str(k)+': '+str(v)+'\n')
             
    fout.close()


############################################
def write_GFs_info(file_path,GFs_info):
    fout = open(file_path, 'w') 
    fout.write('The GFs_info format MCMTpy used:\n')
    fout.write('[number]\t[source_name]\t[source_dep]\t[network_name]\t[station_name]\t[distance]\t[station_depth]\n\n')
    
    for i in range(0,len(GFs_info),1):
        fout.write( str(GFs_info[i]) )
        fout.write('\t')
        fout.write('\n')
        
    fout.close()


############################################
def write_Source_Station_info_MPI(file_path,Source_Station_info_MPI):
    fout = open(file_path, 'w')
    fout.write('The Source_Station_info format MCMTpy used:\n')
    fout.write('\t[source_name]\t[source_lat]\t[source_lon]\t[source_depth]\n')
    fout.write('\t[network_name]\t[station_name]\t[station_lat]\t[source_lon]\t[station_depth]\t[distance]\t[az]\t[baz]\n\n')
    
    for k in range(0,len(Source_Station_info_MPI),1):
        for i in range(0,len(Source_Station_info_MPI[k]),1):
            for j in range(0,len(Source_Station_info_MPI[k][i]),1):
                fout.write(str(Source_Station_info_MPI[k][i][j]))
                fout.write('\t')
            fout.write('\n')
        fout.write('   \n')
        
    fout.close()


############################################
def read_Velocity_model(filename):
    Velocity_model = []
    with open(filename,'r') as f:
        for line in f:
            line = line.strip()                                                 # Remove the leading and trailing space of each line
            line = line.strip('\t')                                             # Remove tabs from each line '\t'
            line = line.strip('\n')                                             # Remove tabs from each newline '\n'
            line = line.split()
            if line == []:
                pass
            else:
                Velocity_model.append( list( line ))
    ### str -->float
    layer_num = len( Velocity_model )
    for i in range(0,layer_num,1):
        for j in range(0,6,1):
            Velocity_model[i][j]=float(Velocity_model[i][j])
    ###
    Velocity_model = np.array(Velocity_model)

    return Velocity_model


############################################
def CPU_split(Source_num,MPI_n):
    '''
    The number of all earthquakes is evenly divided into N calculation cores 
    (MPI_N), and the remaining earthquake numbers are divided into the first 
    several cores
    '''

    cpu_int = Source_num // MPI_n
    cpu_rem = Source_num % MPI_n
    cpu_each_num = np.zeros(shape=(MPI_n,1) )
    
    for i in range(0,MPI_n,1):
        if i < cpu_rem:
            cpu_each_num[i]=cpu_int+1
        else:
            cpu_each_num[i]=cpu_int
    CPU_chunk = np.cumsum(cpu_each_num)
    
    CPU_splits=[]
    for i in range(0,MPI_n,1):
        if i==0:
            index_start = 1
        else:
            index_start = round( CPU_chunk[i-1]+1 )
            
        index_end = round( CPU_chunk[i] )
            
        CPU_splits.append([index_start,index_end])


    return CPU_splits


############################################
def read_Source_Station_info(filename):
    Source_Station_info = []
    one_s_sta_info=[]                                                           # Store a station info corresponding to one source
    with open(filename,'r') as f:
        for line in f:
            line = line.strip()                      
            line = line.strip('\t')               
            line = line.strip('\n')         
            line = line.split()     
            if line == [] and one_s_sta_info != []:   
                Source_Station_info.append(one_s_sta_info)
                one_s_sta_info=[]
            elif len(list(line))==4 or len(list(line))==5:                      # If the line is the source information line or the geophone line, and add the information
                one_s_sta_info.append( list( line ))
                
    if one_s_sta_info != []:
        Source_Station_info.append(one_s_sta_info)
    
    # str --> float
    source_num = len(Source_Station_info)
    for i in range(0,source_num,1):
        station_num = len( Source_Station_info[i] )
        for j in range(0,station_num,1):
            if j==0:
                Source_Station_info[i][0][1]=float(Source_Station_info[i][0][1])
                Source_Station_info[i][0][2]=float(Source_Station_info[i][0][2])
                Source_Station_info[i][0][3]=float(Source_Station_info[i][0][3])
            else:
                Source_Station_info[i][j][2]=float(Source_Station_info[i][j][2])
                Source_Station_info[i][j][3]=float(Source_Station_info[i][j][3])
                Source_Station_info[i][j][4]=float(Source_Station_info[i][j][4])


    return Source_Station_info



