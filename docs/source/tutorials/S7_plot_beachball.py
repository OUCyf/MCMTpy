#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import glob
import obspy
import numpy as np
import matplotlib.pyplot as plt
from obspy import Stream
from obspy.taup import TauPyModel
from obspy.imaging.mopad_wrapper import beach
from MCMTpy import MomentTensor as MTpy

# Helper function
#----------------------------------------------------#
#%% 1.read raw data
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
#%% 2.taup ray trace
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


#%%------------------------------------------------------------#
# 0. set path
example_path = '/Users/yf/3.Project/8.MCMTpy/MCMTpy-master/data/example_yunnan'


#------------------------------------------------------------#
# we expect no parameters need to be changed below
#------------------------------------------------------------#
FM_path=os.path.join(example_path,"YN.202105212148_Inv/dc_inv/Output_YN.202105212148_dc/rank_0_output/chain_0_FM_accept_all")
allfiles_path = os.path.join(example_path,'YN.202105212148_Inv/YN.202105212148_raw/*.SAC') 

## 1. read FM
N=3 # Three parameters are required to describe the focal mechanism
FM_all = np.loadtxt(FM_path)
FM_raw = FM_all[0:,1:N+1] # Define the number of solutions you want to plot
strike_np = FM_raw[:,0]
dip_np = FM_raw[:,1]
rake_np = FM_raw[:,2]
FM = np.vstack((strike_np, dip_np,  rake_np)).T
FM_mean=np.zeros(shape=(N))
for i in range(0,N,1):
    FM_mean[i]=np.mean(FM[0:,i])

## 2.read raw data
data = read_data(allfiles_path)
data.filter('bandpass', freqmin=0.005, freqmax=0.5, corners=4, zerophase=True)

## 3.ray trace with taup
model_path = os.path.join(example_path,"v_model/v_model.npz") 
model = TauPyModel(model=model_path)    # "iasp91"  "prem"
for i in range(0,len(data),1):
    depth = data[i].stats.sac['evdp']
    distance = data[i].stats.sac['dist']
    ray_p,tp,angle_p,ray_s,ts,angle_s = get_taup_tp_ts(model,depth,distance,degree=False)                                            
    data[i].stats.sac["user1"]=angle_p
    data[i].stats.sac["user2"]=angle_s

## 4.1 plot FM_mean
ax0 = plt.gca()
Length_Ball = 100
beach1 = beach(FM_mean, xy=(50, 50),  linewidth=1,width=Length_Ball-1, alpha=1,\
               facecolor='g',bgcolor='w', edgecolor='k',mopad_basis='NED',nofill=False,zorder=1 )
ax0.add_collection(beach1) 
ax0.set_aspect("equal")

## 4.2 plot FM_all
for i in range(0,FM.shape[0],1):
    beach1 = beach(FM[i,:], xy=(50, 50),  linewidth=1,width=Length_Ball-1, alpha=1,\
                facecolor='b',bgcolor='w', edgecolor='orange',mopad_basis='NED',nofill=True,zorder=1 )
    ax0.add_collection(beach1) 
    ax0.set_aspect("equal")

## 4.3 plot backgroud line
beach1 = beach(FM_mean, xy=(50, 50),  linewidth=1,width=Length_Ball-1, alpha=1,\
                facecolor='w',bgcolor='w', edgecolor='k',mopad_basis='NED',nofill=True,zorder=1 )
ax0.add_collection(beach1) 
ax0.set_aspect("equal")

## 5.plot station and waveform
menthod='schmidt'   # 'schmidt' 'wulff'
for i in range(0,len(data),3):
    AZM = data[i].stats.sac['az'] 
    TKO = data[i].stats.sac['user1']
    net_sta_name = data[i].stats.network+'_'+data[i].stats.station    
    X, Y = MTpy.project_beachball(AZM, TKO, R=Length_Ball/2, menthod=menthod)    
    tt=np.linspace(X, X+10, num=len(data[i].data)) 
    ax0.plot(X, Y, "rv", ms=10,zorder=1) 
    ax0.plot(tt, 5*data[i].data/2000000 + Y,  color='black',lw=0.2,alpha=0.6,zorder=1)
    ax0.text(X, Y,net_sta_name,horizontalalignment='right', verticalalignment='center',\
          fontsize=5, color='black',bbox = dict(facecolor = "r", alpha = 0.0),zorder=1) 

## 6. plot P/T/N axis
MT = MTpy.str_dip_rake2MT(strike=FM_mean[0],dip=FM_mean[1],rake=FM_mean[2])
T_axis, P_axis, N_axis = MTpy.MT2TPN(MT)
T = MTpy.vector2str_dip(T_axis)
P = MTpy.vector2str_dip(P_axis)
N = MTpy.vector2str_dip(N_axis)

Tx, Ty = MTpy.project_beachball(AZM=T.strike, TKO=(90-T.dip), R=Length_Ball/2, menthod=menthod)
ax0.text(Tx,Ty,'T',horizontalalignment='center', verticalalignment='center',\
          fontsize=20, color='k',alpha=0.7,zorder=1) 

Px, Py = MTpy.project_beachball(AZM=P.strike, TKO=(90-P.dip), R=Length_Ball/2, menthod=menthod)
ax0.text(Px,Py,'P',horizontalalignment='center', verticalalignment='center',\
          fontsize=20, color='k',alpha=0.7,zorder=1) 

Nx, Ny = MTpy.project_beachball(AZM=N.strike, TKO=(90-N.dip), R=Length_Ball/2, menthod=menthod)
ax0.text(Nx,Ny,'N',horizontalalignment='center', verticalalignment='center',\
          fontsize=20, color='k',alpha=0.7,zorder=1) 


## 7. save figure
ax0.set_xlim(0,100)
ax0.set_ylim(0,100)
ax0.set_axis_off() 
figurename=os.path.join('./S2_figure/beachball.pdf')
plt.savefig(figurename,dpi=800, format="pdf")