#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 20 21:33:46 2021

@author: Fu Yin (yinfu@mail.ustc.edu.cn) at USTC

This script:
    1) function read_GFs_json is preprocess of function build_GFs_pyfk;


Modify history:
    1) Apr 20 21:33:46 2021    ||    Fu Yin at USTC    ||    The initial release.
    2) ...
    
"""

import obspy
import pyasdf
import os
from obspy import UTCDateTime
from obspy.core.inventory import Inventory, Network, Station, Channel, Site




############################################
def get_StationXML(station_network, station_name,station_lat, station_lon, station_depth, srcType):
    inv = Inventory(networks=[],source="homegrown")

    response = obspy.core.inventory.response.Response()

    net = Network(
        # This is the network code according to the SEED standard.
        code        = station_network,
        stations    = [],
        description = "Synthetic seismic data for example",
        start_date  = [])

    sta = Station(
        code          = station_name,
        latitude      = station_lat,
        longitude     = station_lon,
        elevation     = -station_depth,
        creation_date = UTCDateTime(2012, 4, 4, 14, 18, 37),
        site          = Site(name="First station"))

    if srcType=='dc':
        component = ['grn_0','grn_1','grn_2','grn_3','grn_4','grn_5','grn_6','grn_7','grn_8']
    elif srcType=='ep':
        component = ['grn_a','grn_b','grn_c']
    elif srcType=='sf':
        component = ['grn_0','grn_1','grn_2','grn_3','grn_4','grn_5']
    elif srcType=='syn':                                                        # MCMTpy synthetic seismogram in zrt (it not work now)
        component = ['AZ','BR','CT']
    elif srcType=='enz':                                                        # real data in enz
        component = ['E','N','Z']
    elif srcType=='zne':                                                        # real data in zne
        component = ['Z','N','E']
    else:
        pass

    for ch_code in component:
        cha = Channel(
            code          = ch_code,
            location_code = '',
            latitude      = station_lat,
            longitude     = station_lon,
            elevation     = 0,
            depth         = station_depth,
            azimuth       = 0,
            dip           = 0,
            sample_rate   = 0.01)
        # Now tie it all together.
        cha.response = response
        sta.channels.append(cha)

    net.stations.append(sta)
    inv.networks.append(net)


    return inv




############################################
def Add_waveforms_head_info(gf, station_network, station_name, srcType, starttime):
    '''
    NOTE:
        To write information of StationXML and Waveforms into ASDF files, 
        you must ensure that the station information in the waveform data 
        is consistent with the station information in the waveform data
    '''
    if srcType=='dc':
        component     = ['grn_0','grn_1','grn_2','grn_3','grn_4','grn_5','grn_6','grn_7','grn_8']
        component_num = 9
    elif srcType=='ep':
        component     = ['grn_a','grn_b','grn_c']
        component_num = 3 
    elif srcType=='sf':
        component     = ['grn_0','grn_1','grn_2','grn_3','grn_4','grn_5']
        component_num = 6
    elif srcType=='syn':
        component     = ['AZ','BR','CT']
        component_num = 3
    elif srcType=='enz':
        component=['E','N','Z']
        component_num = 3
    elif srcType=='zne':
        component     = ['Z','N','E']
        component_num = 3
    else:
        pass
    
    for i,ch_code in zip(range(0,component_num,1),component):
        gf[0][i].stats.network   = station_network
        gf[0][i].stats.station   = station_name
        gf[0][i].stats.channel   = ch_code
        gf[0][i].stats.starttime = starttime




############################################
def Add_waveforms(gf, catalog_create, file_path, source_name, tp, ts, station_distance, station_az=None, station_baz=None):
    '''
    NOTE:
        To write information of StationXML and Waveforms into ASDF files, 
        you must ensure that the station information in the waveform data 
        is consistent with the station information in the waveform data
    '''
    if not os.path.isfile(file_path):
        with pyasdf.ASDFDataSet(file_path,mpi=False,compression="gzip-3",mode='w') as ds:
            pass

    with pyasdf.ASDFDataSet(file_path,mpi=False,compression="gzip-3",mode='a') as ds:
        #  tag only supports '^[a-z_0-9]+$'，so all uppercase letters are converted to lowercase
        new_tags = str.lower(source_name)
        new_labels = ["tp",str(tp),
                      "ts",str(ts),
                      "dist",str(station_distance),
                      "az",str(station_az),
                      "baz",str(station_baz)]
        ds.add_waveforms(gf[0],
                         tag=new_tags,
                         labels=new_labels,
                         event_id=catalog_create)





############################################
def Add_stationxml(inv, file_path):
    if not os.path.isfile(file_path):
        with pyasdf.ASDFDataSet(file_path,mpi=False,compression="gzip-3",mode='w') as ds:
            pass
    
    with pyasdf.ASDFDataSet(file_path,mpi=False,compression="gzip-3",mode='a') as ds:
        try:ds.add_stationxml(inv)
        except Exception: pass 





############################################
def get_QuakeML(source_name, source_lat, source_lon, source_depth, source_time):
    catalog_create=obspy.core.event.catalog.Catalog()
    magnitudes=[obspy.core.event.magnitude.Magnitude(resource_id='ResourceIdentifier',
                                                     mag=3,
                                                     magnitude_type='Mw' )]

    origins=[obspy.core.event.origin.Origin(resource_id='ResourceIdentifier',
                                            time=source_time,                   # the current time
                                            longitude=source_lon,
                                            latitude=source_lat,
                                            depth=source_depth,
                                            evaluation_mode='manual')]

    event_descriptions=[obspy.core.event.event.EventDescription(text='Synthetic seismic data')]

    event=obspy.core.event.event.Event(resource_id=source_name,
                                       event_type='not reported',
                                       magnitudes=magnitudes,
                                       origins=origins,
                                       event_descriptions=event_descriptions)

    catalog_create.append(event)


    return catalog_create





############################################
def Add_quakeml(file_path, catalog_create):
    if not os.path.isfile(file_path ):
        with pyasdf.ASDFDataSet(file_path ,mpi=False,compression="gzip-3",mode='w') as ds:
            pass
    
    with pyasdf.ASDFDataSet(file_path ,mpi=False,compression="gzip-3",mode='a') as ds:
        try:ds.add_quakeml(catalog_create)
        except Exception: pass 



