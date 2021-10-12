#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import pyasdf


# The root directory of the project
example_path = '/Users/yf/3.Project/8.MCMTpy/MCMTpy-master/data/example_yunnan'

#------------------------------------------------------------#
# we expect no parameters need to be changed below
#------------------------------------------------------------#
## change the path
h5_path = os.path.join(example_path,'syn/Synthetic/SYN_test.h5')
source_name = 'source_syn'

### read h5 data --> Sta_data (dict)
ds_raw=pyasdf.ASDFDataSet(h5_path,mpi=False,mode='r')
Station_list = ds_raw.waveforms.list()
Station_num = len(Station_list)

Sta_data={}
for i in range(0,Station_num,1):
    info={}
    tags = ds_raw.waveforms[Station_list[i]].get_waveform_tags()
    if source_name in tags:
        raw_sac = ds_raw.waveforms[Station_list[i]][source_name]
        tp = float( ds_raw.waveforms[Station_list[i]][source_name][0].stats.asdf['labels'][1] )
        ts = float( ds_raw.waveforms[Station_list[i]][source_name][0].stats.asdf['labels'][3] )
        info.update({ "tp": tp})  
        info.update({ "ts": ts})
        info.update({ "data": raw_sac})
        Sta_data.update({ Station_list[i]: info})

# plot data
ax = Sta_data['YN.YUM']['data'].plot()
Sta_data['YN.YUM']['data'][0].stats