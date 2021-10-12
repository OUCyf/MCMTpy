#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 21 00:20:35 2021

@author: Fu Yin (yinfu@mail.ustc.edu.cn) at USTC

This script:
    1) plot hist


Modify history:
    1) Mar 21 00:20:35 2021   ||    Fu Yin at USTC    ||    The initial release.
    2) Jun 28 11:08:02 2021   ||    Fu Yin at USTC    ||    Add figure format parameter.
"""

import numpy as np
import matplotlib.pyplot as plt
import math
import os
from obspy.imaging.mopad_wrapper import beach
from MCMTpy.gfs.pyfk_GFs.read_json import read_GFs_json
from MCMTpy.sampler.pyfk_MH.read_sampler_json import read_Inv_json
from MCMTpy.sampler.pyfk_MH.sampler_module import config_fk,get_GFs_sta_info,get_Sta_rand_data,correlate_maxlag,zero_padding
from MCMTpy.utils.rotate import rotate




def plot_waveform(Inv_json_file,plot_output_path,NET_STA,FM_best,line_n_sta,max_p_ylim,max_s_ylim,max_surf_ylim,plot_comp,plot_comp_name,fig_format,Geometrical_spreading_mode,Amplitude_normalization_mode):
    
    #%% ###########
    # 1.Initialization parameter
    Station_num = len(NET_STA)
    Inv_para,Sta_inv_info,Sta_raw_data,Stf = read_Inv_json(Inv_json_file)
    Prepro_para,_,Source_Station_info_MPI = read_GFs_json(Inv_para["GFs_json_file"])

    Inv_para["Geometrical_spreading_mode"] = Geometrical_spreading_mode
    Inv_para["Amplitude_normalization_mode"] = Amplitude_normalization_mode
    InvType = Inv_para["InvType"]
    Raw_p_n0 = Inv_para["Raw_p_n0"]
    Rand_p_n0 = Prepro_para["samples_before_first_arrival"]
    dt = Inv_para["Dt"]
    Distance_0= Inv_para["Distance_0"]

    source_prem,_,config_prem=config_fk(InvType)




    #%% ######
    # 2. Generate the most fitting data
    # % c.1.forward model
    GFs_sta_info = get_GFs_sta_info(FM_best[-4:-1],NET_STA,Prepro_para,Source_Station_info_MPI)
    # % c.2.forward model FM
    if InvType=='mt':
        source_mechanism = np.array(FM_best[0:7])
        mw = source_mechanism[0]
        source_mechanism[0] = np.power(10., 1.5 * mw + 16.1)
    elif InvType=='dc':
        mw = FM_best[0]
        source_mechanism = np.array(FM_best[0:4])
    elif InvType=='sf':
        source_mechanism = np.array(FM_best[0:3])
        mw = source_mechanism[0]
        source_mechanism[0] = np.power(10., 1.5 * mw + 16.1)
    source_prem.update_source_mechanism(source_mechanism)
    # % c.3.forward model syn
    Sta_rand_data = get_Sta_rand_data(GFs_sta_info,config_prem,Stf)





    #%%######
    # 4. Paint size Settings
    Comp_num = 0                                                                # Count the number of '1' in plot_comp
    for i in range(0,3,1):
        for j in range(0,3,1):
            if plot_comp[i][j]==1:
                Comp_num += 1

    line_num = math.ceil(Station_num/line_n_sta)                                # How many total rows do I need to draw this
    column_num = line_n_sta*Comp_num                                            # How many total columns do I have to draw

    height_ratios = list(1 for i in range(0,line_num+1))                        # The height ratio of each subgraph max_ylim
    # width_ratios = list(1 for i in range(0,line_num))                         # Adjust according to the length of time
    p_length_all=[];s_length_all=[];surf_length_all=[];
    for ii in range(0,Station_num,1):
        net_sta_name = str(NET_STA[ii][0])+'.'+str(NET_STA[ii][1])
        p_length = Sta_inv_info[net_sta_name]['P_win'][1]-Sta_inv_info[net_sta_name]['P_win'][0]
        s_length = Sta_inv_info[net_sta_name]['S_win'][1]-Sta_inv_info[net_sta_name]['S_win'][0]
        surf_length = Sta_inv_info[net_sta_name]['Surf_win'][1]-Sta_inv_info[net_sta_name]['Surf_win'][0]
        p_length_all.append(p_length)
        s_length_all.append(s_length)
        surf_length_all.append(surf_length)
    
    p_length_max = max(p_length_all)
    s_length_max = max(s_length_all)
    surf_length_max = max(surf_length_all)
    tt_length=[p_length_max,s_length_max,surf_length_max ]
    width_ratios=[]
    for i in range(0,line_n_sta,1):
        for j in range(0,3,1):                                                  # loop P S Surf
            for k in range(0,3,1):                                              # loop Z R T
                if plot_comp[j][k] == 1:
                    width_ratios.append(tt_length[j])






    #%% #############
    # 5. A formal drawing
    # The FigSize for 6 rows and 6 columns is standard (20,20), so scale it up to this size
    fig5, axs = plt.subplots(line_num+1 ,column_num,  dpi=800,figsize=(column_num/6*20, (line_num+1)/6*20),
                              gridspec_kw=dict(height_ratios=height_ratios,width_ratios=width_ratios))
    plt.rc('font',family='Times New Roman')
    for i in range(0,line_num+1):                                               # Turn off all axis information
        for j in range(0,column_num):
            # print(i,j)
            axs[i,j].set_axis_off()
            axs[line_num,j].set_xlabel("Time(s)",fontsize=15)
    
    
    #%% ######
    # 5.1 plot beachball
    number=line_num+1
    ax0=fig5.add_subplot(int(number),1,1)

    if len(source_mechanism[1:])==3:
        MT = source_mechanism[1:]
    else:
        Mxx = source_mechanism[1]
        Mxy = source_mechanism[2]
        Mxz = source_mechanism[3]
        Myy = source_mechanism[4]
        Myz = source_mechanism[5]
        Mzz = source_mechanism[6]
        MT= [Mxx,  Myy,  Mzz,  Mxy,  Mxz,  Myz]

    beach1 = beach(MT, xy=(0.5, 0.5), width=1,facecolor='b',mopad_basis='NED')
    ax0.add_collection(beach1) 
    ax0.set_aspect("equal")
    ax0.set_xlim(0, 1.5*column_num)  
    ax0.set_ylim(0,1)
    ax0.set_axis_off() 
    
    str_FM =''
    for i in range(1,len(source_mechanism)):
        str_FM += str(np.around(float(source_mechanism[i]), decimals=1))+"   "
    str_text="      Lat:   "+str(np.around(FM_best[-4], decimals=3))+"   Lon:   "+str(np.around(FM_best[-3], decimals=3))\
            +"   Dep:   "+str(np.around(FM_best[-2], decimals=1))+"   T0:   "+str(np.around(FM_best[-1], decimals=1))+'\n'\
            +"      Fm:   "+str_FM+"   Mw:   "+str(np.around(mw, decimals=2))
        
    ax0.text(1, 0.7,str_text,horizontalalignment='left', verticalalignment='center',fontsize=25, color='black')


    #%% ######
    # 5.2 Loop through all stations to draw each subgraph
    for m in range(0,Station_num,1):
        j = 0                                                                   # j parameter is used to count how many waveforms each station has to draw
        jj = 0
        net_sta_name = str(NET_STA[m][0])+'.'+str(NET_STA[m][1])
        
        distance = GFs_sta_info[net_sta_name][5]
        baz = GFs_sta_info[net_sta_name][8]
        Raw_tp = Sta_raw_data[net_sta_name]['tp']
        Raw_ts = Sta_raw_data[net_sta_name]['ts']
        Raw_tsurf = Raw_ts/0.913
        Rand_tp = Sta_rand_data[net_sta_name]['tp']
        Rand_ts = Sta_rand_data[net_sta_name]['ts']
        Rand_tsurf = Rand_ts/0.913
        
        P_win = Sta_inv_info[net_sta_name]['P_win']
        S_win = Sta_inv_info[net_sta_name]['S_win']
        Surf_win = Sta_inv_info[net_sta_name]['Surf_win']
        
        P_inv_comp = Sta_inv_info[net_sta_name]['P_inv_comp']
        S_inv_comp = Sta_inv_info[net_sta_name]['S_inv_comp']
        Surf_inv_comp = Sta_inv_info[net_sta_name]['Surf_inv_comp']
        
        Phase_weight = Sta_inv_info[net_sta_name]['Phase_weight']
        Distance_weight = Sta_inv_info[net_sta_name]['Distance_weight'][0]
        
        P_filter = Sta_inv_info[net_sta_name]['P_filter']
        S_filter = Sta_inv_info[net_sta_name]['S_filter']
        Surf_filter = Sta_inv_info[net_sta_name]['Surf_filter']

        P_maxlag = Sta_inv_info[net_sta_name]['P_maxlag'][0]
        S_maxlag = Sta_inv_info[net_sta_name]['S_maxlag'][0]
        Surf_maxlag = Sta_inv_info[net_sta_name]['Surf_maxlag'][0]
        P_maxlag = round(P_maxlag/dt)
        S_maxlag = round(S_maxlag/dt)
        Surf_maxlag = round(Surf_maxlag/dt)

        Raw_data = Sta_raw_data[net_sta_name]['data']                           # 3 Trace(s) in Stream: ENZ
        Rand_data = Sta_rand_data[net_sta_name]['data']
        
        ## rotate 'ENZ->ZRT'
        for cc in range(0,len(Raw_data)):
            Raw_data[cc].stats.back_azimuth = float( baz )
        Raw_data = rotate(data_raw=Raw_data, chn=['ENZ','ZRT'])


        #####
        P_win_length = round((P_win[1]-P_win[0])/dt)
        S_win_length = round((S_win[1]-S_win[0])/dt)
        Surf_win_length = round((Surf_win[1]-Surf_win[0])/dt)
        
        raw_P_start = Raw_p_n0 + round(P_win[0]/dt)
        raw_S_start = Raw_p_n0 + round( (Raw_ts+S_win[0]-Raw_tp)/dt )
        raw_Surf_start = Raw_p_n0 + round( (Raw_tsurf+Surf_win[0]-Raw_tp)/dt )
        
        T0 = FM_best[-1]
        rand_P_start = Rand_p_n0 + round(P_win[0]/dt) + round((Raw_tp-T0-Rand_tp)/dt)
        rand_S_start = Rand_p_n0 + round( (Rand_ts+S_win[0]-Rand_tp)/dt ) + round((Raw_ts-T0-Rand_ts)/dt)
        rand_Surf_start = Rand_p_n0 + round( (Rand_tsurf+Surf_win[0]-Rand_tp)/dt ) + round((Raw_tsurf-T0-Rand_tsurf)/dt)
        # rand_P_start = Rand_p_n0 + round(P_win[0]/dt)
        # rand_S_start = Rand_p_n0 + round( (Rand_ts+S_win[0]-Rand_tp)/dt )
        # rand_Surf_start = Rand_p_n0 + round( (Rand_tsurf+Surf_win[0]-Rand_tp)/dt )

        if (raw_P_start<0) or (raw_S_start<0) or (raw_Surf_start<0) or (rand_P_start<0) or (rand_S_start<0) or (rand_Surf_start<0):
            raise ValueError('P_win or S_win or Surf_win exceeded the start time of the data')
            
        Raw_P_win_n0 = [raw_P_start, raw_P_start + P_win_length]
        Raw_S_win_n0 = [raw_S_start, raw_S_start + S_win_length]
        Raw_Surf_win_n0 = [raw_Surf_start, raw_Surf_start + Surf_win_length]
        
        
        Rand_P_win_n0 = [rand_P_start, rand_P_start + P_win_length]
        Rand_S_win_n0 = [rand_S_start, rand_S_start + S_win_length]
        Rand_Surf_win_n0 = [rand_Surf_start, rand_Surf_start + Surf_win_length]



        # 5.2.1  Geometrical_spreading_mode
        MISFIT1_one_sta=0
        if Inv_para['Geometrical_spreading_mode']==True:
            #### a. P wave
            for i in range (0,3,1):                                             # p ZRT comp
                if plot_comp[0][i] == 1:
                    line_order = m//line_n_sta+1                                # The line in the drawing
                    colume_order = (m%line_n_sta) * Comp_num+j                  # The row in the drawing
                    axs[line_order,colume_order].set_axis_off()
                    j+=1
                    print('Now plotting figure (line_order,colume_order):',line_order,colume_order)
                
            
                    Raw_data_filter = Raw_data.copy()
                    Raw_data_filter.filter('bandpass', freqmin=P_filter[0], freqmax=P_filter[1], corners=4, zerophase=True)
                    Rand_data_filter = Rand_data.copy()
                    Rand_data_filter.filter('bandpass', freqmin=P_filter[0], freqmax=P_filter[1], corners=4, zerophase=True)
                    
                    Raw_data_np = Raw_data_filter[i].data
                    Rand_data_np = Rand_data_filter[i].data
                    Raw_data_np_p = Raw_data_np[ Raw_P_win_n0[0]:Raw_P_win_n0[1] ] 
                    Rand_data_np_p = Rand_data_np[ Rand_P_win_n0[0]:Rand_P_win_n0[1] ] 
                    
                    _, C_max_maxlag, shift=correlate_maxlag(Raw_data_np_p,Rand_data_np_p,P_maxlag,normalized=True)
                    # Rand_data_np_p_shift = zero_padding(Rand_data_np_p,shift)
                    Rand_data_np_p_shift = Rand_data_np[ Rand_P_win_n0[0]-shift:Rand_P_win_n0[1]-shift ] 
                    MISFIT1_one_sta += Phase_weight[0] * (distance/Distance_0)**2*math.sqrt(np.mean((Raw_data_np_p-Rand_data_np_p_shift)**2)) #  * (1-C_max_maxlag)
                    
                    tt = np.linspace(Raw_tp+P_win[0], Raw_tp+P_win[1], num=len(Raw_data_np_p))  # X axis time scale
                    ## plot
                    if line_order==1:
                        str_text=plot_comp_name[0][i]                               # Name comments in subgraphs such as' P_r'
                        axs[line_order, colume_order].set_title(str_text,fontsize=20,fontweight='bold', alpha=0.7)  # Subgraph title
                        # axs[line_order, colume_order].text((tt[-1]+tt[0])/2,1*max_ylim,str_text,horizontalalignment='center', verticalalignment='bottom',
                        #            fontsize=18, color='black')
                        
                    if P_inv_comp[i]=='yes':
                        jj+=1
                        p1,=axs[line_order, colume_order].plot(tt,Raw_data_np_p, "k-",lw=2,alpha=0.6) 
                        p2,=axs[line_order, colume_order].plot(tt,Rand_data_np_p_shift, linestyle='-', color='red', lw=3, alpha=0.6)
                        
                        cc=np.around(C_max_maxlag, decimals=2)
                        # If it is the first submap of each station, mark the station name and epicenter distance on the submap
                        if jj ==1:
                            str_text=net_sta_name+"/"+str(np.around(distance, decimals=2)) 
                            axs[line_order, colume_order].text(tt[0], 0.9*max_p_ylim,str_text,horizontalalignment='left', verticalalignment='top',
                                                               fontsize=15, color='black',bbox = dict(facecolor = "k", alpha = 0.05))
                        cc_shift_text = str(cc)+'/'+str(shift*dt)
                        axs[line_order, colume_order].text(tt[0], -0.8*max_p_ylim,cc_shift_text,horizontalalignment='left', verticalalignment='top',
                                                               fontsize=14, color='black')                      # Mark the correlation coefficient SHIFT
                        
                            
                        axs[line_order, colume_order].set_ylim(-max_p_ylim,max_p_ylim)  # Set the Y-axis range
                            
                        axs[line_order,colume_order].set_axis_on()
                        axs[line_order,colume_order].spines['top'].set_visible(False)
                        axs[line_order,colume_order].spines['right'].set_visible(False)
                        axs[line_order,colume_order].spines['bottom'].set_visible(True)
                        axs[line_order,colume_order].spines['left'].set_visible(False)
                        axs[line_order, colume_order].get_xaxis().set_visible(True)
                        axs[line_order, colume_order].get_yaxis().set_visible(False)
                        # axs[line_order,colume_order].set_xlabel("Time(s)",fontsize=18)
                        
                
                
                
            #### b. S wave
            for i in range (0,3,1):
                if plot_comp[1][i] == 1:
                    line_order = m//line_n_sta+1
                    colume_order = (m%line_n_sta) * Comp_num+j
                    axs[line_order,colume_order].set_axis_off()
                    j+=1
                    print('Now plotting figure (line_order,colume_order):',line_order,colume_order)
                    
                    Raw_data_filter = Raw_data.copy()
                    Raw_data_filter.filter('bandpass', freqmin=S_filter[0], freqmax=S_filter[1], corners=4, zerophase=True)
                    Rand_data_filter = Rand_data.copy()
                    Rand_data_filter.filter('bandpass', freqmin=S_filter[0], freqmax=S_filter[1], corners=4, zerophase=True)
                    
                    Raw_data_np = Raw_data_filter[i].data
                    Rand_data_np = Rand_data_filter[i].data
                    Raw_data_np_s = Raw_data_np[ Raw_S_win_n0[0]:Raw_S_win_n0[1] ]
                    Rand_data_np_s = Rand_data_np[ Rand_S_win_n0[0]:Rand_S_win_n0[1] ] 
                    
                    _, C_max_maxlag, shift=correlate_maxlag(Raw_data_np_s,Rand_data_np_s,S_maxlag,normalized=True)
                    # Rand_data_np_s_shift = zero_padding(Rand_data_np_s,shift)
                    Rand_data_np_s_shift = Rand_data_np[ Rand_S_win_n0[0]-shift:Rand_S_win_n0[1]-shift ] 
                    MISFIT1_one_sta += Phase_weight[1] * (distance/Distance_0)**2*math.sqrt(np.mean((Raw_data_np_s-Rand_data_np_s_shift)**2)) #  * (1-C_max_maxlag)
                    
                    tt = np.linspace(Raw_tp+P_win[0], Raw_tp+P_win[1], num=len(Raw_data_np_p))
                    ## 画图命令
                    if line_order==1:
                        str_text=plot_comp_name[1][i]
                        axs[line_order, colume_order].set_title(str_text,fontsize=20,fontweight='bold', alpha=0.7)
                        # axs[line_order, colume_order].text((tt[-1]+tt[0])/2,1*max_ylim,str_text,horizontalalignment='center', verticalalignment='bottom',
                        #            fontsize=18, color='black')
                        
                    if S_inv_comp[i]=='yes':
                        p1,=axs[line_order, colume_order].plot(tt,Raw_data_np_s, "k-",lw=2,alpha=0.6) 
                        p2,=axs[line_order, colume_order].plot(tt,Rand_data_np_s_shift, linestyle='-', color='red', lw=3, alpha=0.6)
                        
                        cc=np.around(C_max_maxlag, decimals=2)
                        if jj ==1:
                            str_text=net_sta_name+"/"+str(np.around(distance, decimals=2)) 
                            axs[line_order, colume_order].text(tt[0], 0.9*max_s_ylim,str_text,horizontalalignment='left', verticalalignment='top',
                                                               fontsize=15, color='black',bbox = dict(facecolor = "orange", alpha = 0.1))
                        cc_shift_text = str(cc)+'/'+str(shift*dt)
                        axs[line_order, colume_order].text(tt[0], -0.8*max_s_ylim,cc_shift_text,horizontalalignment='left', verticalalignment='top',
                                                               fontsize=14, color='black')
                        
                        axs[line_order, colume_order].set_ylim(-max_s_ylim,max_s_ylim)
                        
                        axs[line_order,colume_order].set_axis_on()
                        axs[line_order,colume_order].spines['top'].set_visible(False)
                        axs[line_order,colume_order].spines['right'].set_visible(False) 
                        axs[line_order,colume_order].spines['bottom'].set_visible(True) 
                        axs[line_order,colume_order].spines['left'].set_visible(False)
                        axs[line_order, colume_order].get_xaxis().set_visible(True)
                        axs[line_order, colume_order].get_yaxis().set_visible(False)
                        # axs[line_order,colume_order].set_xlabel("Time(s)",fontsize=18)


                
                
            #### c. Surf wave
            for i in range (0,3,1):
                if plot_comp[2][i] == 1:
                    line_order = m//line_n_sta+1
                    colume_order = (m%line_n_sta) * Comp_num+j
                    axs[line_order,colume_order].set_axis_off()
                    j+=1
                    print('Now plotting figure (line_order,colume_order):',line_order,colume_order)
                    
                    Raw_data_filter = Raw_data.copy()
                    Raw_data_filter.filter('bandpass', freqmin=Surf_filter[0], freqmax=Surf_filter[1], corners=4, zerophase=True)
                    Rand_data_filter = Rand_data.copy()
                    Rand_data_filter.filter('bandpass', freqmin=Surf_filter[0], freqmax=Surf_filter[1], corners=4, zerophase=True)
                    
                    Raw_data_np = Raw_data_filter[i].data
                    Rand_data_np = Rand_data_filter[i].data
                    Raw_data_np_surf = Raw_data_np[ Raw_Surf_win_n0[0]:Raw_Surf_win_n0[1] ]  
                    Rand_data_np_surf = Rand_data_np[ Rand_Surf_win_n0[0]:Rand_Surf_win_n0[1] ] 
                    
                    _, C_max_maxlag, shift=correlate_maxlag(Raw_data_np_surf,Rand_data_np_surf,Surf_maxlag,normalized=True)
                    # Rand_data_np_surf_shift = zero_padding(Rand_data_np_surf,shift)
                    Rand_data_np_surf_shift = Rand_data_np[ Rand_Surf_win_n0[0]-shift:Rand_Surf_win_n0[1]-shift ] 
                    MISFIT1_one_sta += Phase_weight[2] * (distance/Distance_0)**1*math.sqrt(np.mean((Raw_data_np_surf-Rand_data_np_surf_shift)**2)) #  * (1-C_max_maxlag)
                
                    tt = np.linspace(Raw_tsurf+Surf_win[0], Raw_tsurf+Surf_win[1], num=len(Raw_data_np_surf))
                    
                    if line_order==1:
                        str_text=plot_comp_name[2][i]
                        axs[line_order, colume_order].set_title(str_text,fontsize=20,fontweight='bold', alpha=0.7)
                            # axs[line_order, colume_order].text((tt[-1]+tt[0])/2,1*max_ylim,str_text,horizontalalignment='center', verticalalignment='bottom',
                            #                                    fontsize=18, color='black')
                
                    if Surf_inv_comp[i]=='yes':
                        p1,=axs[line_order, colume_order].plot(tt,Raw_data_np_surf, "k-",lw=2,alpha=0.6) 
                        p2,=axs[line_order, colume_order].plot(tt,Rand_data_np_surf_shift, linestyle='-', color='red', lw=3, alpha=0.6)
                        
                        cc=np.around(C_max_maxlag, decimals=2)
                        if jj ==1:
                            str_text=net_sta_name+"/"+str(np.around(distance, decimals=2)) 
                            axs[line_order, colume_order].text(tt[0], 0.9*max_surf_ylim,str_text,horizontalalignment='left', verticalalignment='top',
                                                               fontsize=15, color='black',bbox = dict(facecolor = "orange", alpha = 0.1))
                        cc_shift_text = str(cc)+'/'+str(shift*dt)
                        axs[line_order, colume_order].text(tt[0], -0.8*max_surf_ylim,cc_shift_text,horizontalalignment='left', verticalalignment='top',
                                                               fontsize=14, color='black')
                        
                        axs[line_order, colume_order].set_ylim(-max_surf_ylim,max_surf_ylim)
                    
                        axs[line_order,colume_order].set_axis_on()
                        axs[line_order,colume_order].spines['top'].set_visible(False)
                        axs[line_order,colume_order].spines['right'].set_visible(False)
                        axs[line_order,colume_order].spines['bottom'].set_visible(True)
                        axs[line_order,colume_order].spines['left'].set_visible(False)
                        axs[line_order, colume_order].get_xaxis().set_visible(True)
                        axs[line_order, colume_order].get_yaxis().set_visible(False)
                        # axs[line_order,colume_order].set_xlabel("Time(s)",fontsize=18)
                    
                    
                
                
                
        # 5.2.2 amplitude normalization
        else:
            #### a. P wave
            for i in range (0,3,1):
                if plot_comp[0][i] == 1:
                    line_order = m//line_n_sta+1
                    colume_order = (m%line_n_sta) * Comp_num+j
                    axs[line_order,colume_order].set_axis_off()
                    j+=1
                    print('Now plotting figure (line_order,colume_order):',line_order,colume_order)
                    
                    Raw_data_filter = Raw_data.copy()
                    Raw_data_filter.filter('bandpass', freqmin=P_filter[0], freqmax=P_filter[1], corners=4, zerophase=True)
                    Rand_data_filter = Rand_data.copy()
                    Rand_data_filter.filter('bandpass', freqmin=P_filter[0], freqmax=P_filter[1], corners=4, zerophase=True)
                    
                    Raw_data_np = Raw_data_filter[i].data
                    Rand_data_np = Rand_data_filter[i].data
                    Raw_data_np_p = Raw_data_np[ Raw_P_win_n0[0]:Raw_P_win_n0[1] ] / math.sqrt( np.sum( Raw_data_np[ Raw_P_win_n0[0]:Raw_P_win_n0[1] ]**2, axis=0) )
                    Rand_data_np_p = Rand_data_np[ Rand_P_win_n0[0]:Rand_P_win_n0[1] ] / math.sqrt( np.sum( Rand_data_np[ Rand_P_win_n0[0]:Rand_P_win_n0[1] ]**2, axis=0) )
                    
                    _, C_max_maxlag, shift=correlate_maxlag(Raw_data_np_p,Rand_data_np_p,P_maxlag,normalized=True)
                    # np2_zeros = zero_padding(Rand_data_np_p,shift)
                    # Rand_data_np_p_shift = np2_zeros / math.sqrt( np.sum(np2_zeros**2, axis=0) )
                    Rand_data_np_p_shift = Rand_data_np[ Rand_P_win_n0[0]-shift:Rand_P_win_n0[1]-shift ] / math.sqrt( np.sum( Rand_data_np[ Rand_P_win_n0[0]-shift:Rand_P_win_n0[1]-shift ]**2, axis=0) )
                    
                    
                    MISFIT1_one_sta += Phase_weight[0]*math.sqrt(np.mean((Raw_data_np_p-Rand_data_np_p_shift)**2)) # *(1-C_max_maxlag)
                    
                    tt = np.linspace(Raw_tp+P_win[0], Raw_tp+P_win[1], num=len(Raw_data_np_p))
                    ## 画图命令
                    if line_order==1:
                        str_text=plot_comp_name[0][i]
                        axs[line_order, colume_order].set_title(str_text,fontsize=20,fontweight='bold', alpha=0.7)
                        # axs[line_order, colume_order].text((tt[-1]+tt[0])/2,1*max_ylim,str_text,horizontalalignment='center', verticalalignment='bottom',
                        #            fontsize=18, color='black')
                        
                    if P_inv_comp[i]=='yes':
                        jj+=1
                        p1,=axs[line_order, colume_order].plot(tt,Raw_data_np_p, "k-",lw=2,alpha=0.6) 
                        p2,=axs[line_order, colume_order].plot(tt,Rand_data_np_p_shift, linestyle='-', color='red', lw=3, alpha=0.6)
                        
                        cc=np.around(C_max_maxlag, decimals=2)
                        if jj ==1:
                            str_text=net_sta_name+"/"+str(np.around(distance, decimals=2)) 
                            axs[line_order, colume_order].text(tt[0], 0.9*max_p_ylim,str_text,horizontalalignment='left', verticalalignment='top',
                                                               fontsize=15, color='black',bbox = dict(facecolor = "orange", alpha = 0.1))
                        cc_shift_text = str(cc)+'/'+str(shift*dt)
                        axs[line_order, colume_order].text(tt[0], -0.8*max_p_ylim,cc_shift_text,horizontalalignment='left', verticalalignment='top',
                                                               fontsize=14, color='black')
                        
                        axs[line_order, colume_order].set_ylim(-max_p_ylim,max_p_ylim)
                        
                        axs[line_order,colume_order].set_axis_on()
                        axs[line_order,colume_order].spines['top'].set_visible(False)
                        axs[line_order,colume_order].spines['right'].set_visible(False)
                        axs[line_order,colume_order].spines['bottom'].set_visible(True)
                        axs[line_order,colume_order].spines['left'].set_visible(False)
                        axs[line_order, colume_order].get_xaxis().set_visible(True)
                        axs[line_order, colume_order].get_yaxis().set_visible(False)
                        # axs[line_order,colume_order].set_xlabel("Time(s)",fontsize=18)

                
            #### b. S wave
            for i in range (0,3,1):
                if plot_comp[1][i] == 1:
                    line_order = m//line_n_sta+1
                    colume_order = (m%line_n_sta) * Comp_num+j
                    axs[line_order,colume_order].set_axis_off()
                    j+=1
                    print('Now plotting figure (line_order,colume_order):',line_order,colume_order)
                    
                    Raw_data_filter = Raw_data.copy()
                    Raw_data_filter.filter('bandpass', freqmin=S_filter[0], freqmax=S_filter[1], corners=4, zerophase=True)
                    Rand_data_filter = Rand_data.copy()
                    Rand_data_filter.filter('bandpass', freqmin=S_filter[0], freqmax=S_filter[1], corners=4, zerophase=True)
                    
                    Raw_data_np = Raw_data_filter[i].data
                    Rand_data_np = Rand_data_filter[i].data
                    Raw_data_np_s = Raw_data_np[ Raw_S_win_n0[0]:Raw_S_win_n0[1] ] / math.sqrt( np.sum( Raw_data_np[ Raw_S_win_n0[0]:Raw_S_win_n0[1] ]**2, axis=0) )
                    Rand_data_np_s = Rand_data_np[ Rand_S_win_n0[0]:Rand_S_win_n0[1] ] / math.sqrt( np.sum( Rand_data_np[ Rand_S_win_n0[0]:Rand_S_win_n0[1] ]**2, axis=0) )
                    
                    _, C_max_maxlag, shift=correlate_maxlag(Raw_data_np_s,Rand_data_np_s,S_maxlag,normalized=True)
                    # np2_zeros = zero_padding(Rand_data_np_s,shift)
                    # Rand_data_np_s_shift =  np2_zeros / math.sqrt( np.sum( np2_zeros**2, axis=0) )
                    Rand_data_np_s_shift = Rand_data_np[ Rand_S_win_n0[0]-shift:Rand_S_win_n0[1]-shift ] / math.sqrt( np.sum( Rand_data_np[ Rand_S_win_n0[0]-shift:Rand_S_win_n0[1]-shift ]**2, axis=0) )
                    
                    MISFIT1_one_sta += Phase_weight[1]*math.sqrt(np.mean((Raw_data_np_s-Rand_data_np_s_shift)**2)) # *(1-C_max_maxlag)
                    # print(shift)
                    tt = np.linspace(Raw_ts+S_win[0], Raw_ts+S_win[1], num=len(Raw_data_np_s))
                    
                    if line_order==1:
                        str_text=plot_comp_name[1][i]
                        axs[line_order, colume_order].set_title(str_text,fontsize=20,fontweight='bold', alpha=0.7)
                        # axs[line_order, colume_order].text((tt[-1]+tt[0])/2,1*max_ylim,str_text,horizontalalignment='center', verticalalignment='bottom',
                        #            fontsize=18, color='black')
                        
                    if S_inv_comp[i]=='yes':
                        jj+=1
                        p1,=axs[line_order, colume_order].plot(tt,Raw_data_np_s, "k-",lw=2,alpha=0.6) 
                        p2,=axs[line_order, colume_order].plot(tt,Rand_data_np_s_shift, linestyle='-', color='red', lw=3, alpha=0.6)
                        
                        cc=np.around(C_max_maxlag, decimals=2)
                        if jj ==1:
                            str_text=net_sta_name+"/"+str(np.around(distance, decimals=2)) 
                            axs[line_order, colume_order].text(tt[0], 0.9*max_s_ylim,str_text,horizontalalignment='left', verticalalignment='top',
                                                               fontsize=15, color='black',bbox = dict(facecolor = "orange", alpha = 0.1))
                        cc_shift_text = str(cc)+'/'+str(shift*dt)
                        axs[line_order, colume_order].text(tt[0], -0.8*max_s_ylim,cc_shift_text,horizontalalignment='left', verticalalignment='top',
                                                               fontsize=14, color='black')
                        
                        axs[line_order, colume_order].set_ylim(-max_s_ylim,max_s_ylim)
                        
                        axs[line_order,colume_order].set_axis_on()
                        axs[line_order,colume_order].spines['top'].set_visible(False)
                        axs[line_order,colume_order].spines['right'].set_visible(False)
                        axs[line_order,colume_order].spines['bottom'].set_visible(True)
                        axs[line_order,colume_order].spines['left'].set_visible(False)
                        axs[line_order, colume_order].get_xaxis().set_visible(True)
                        axs[line_order, colume_order].get_yaxis().set_visible(False)
                        # axs[line_order,colume_order].set_xlabel("Time(s)",fontsize=18)
                        
                
            #### c. Surf wave
            for i in range (0,3,1):
                if plot_comp[2][i] == 1:
                    line_order = m//line_n_sta+1
                    colume_order = (m%line_n_sta) * Comp_num+j
                    axs[line_order,colume_order].set_axis_off()
                    j+=1
                    print('Now plotting figure (line_order,colume_order):',line_order,colume_order)
                    
                    Raw_data_filter = Raw_data.copy()
                    Raw_data_filter.filter('bandpass', freqmin=Surf_filter[0], freqmax=Surf_filter[1], corners=4, zerophase=True)
                    Rand_data_filter = Rand_data.copy()
                    Rand_data_filter.filter('bandpass', freqmin=Surf_filter[0], freqmax=Surf_filter[1], corners=4, zerophase=True)
                    
                    Raw_data_np = Raw_data_filter[i].data
                    Rand_data_np = Rand_data_filter[i].data
                    Raw_data_np_surf = Raw_data_np[ Raw_Surf_win_n0[0]:Raw_Surf_win_n0[1] ] / math.sqrt( np.sum( Raw_data_np[ Raw_Surf_win_n0[0]:Raw_Surf_win_n0[1] ]**2, axis=0) )
                    Rand_data_np_surf = Rand_data_np[ Rand_Surf_win_n0[0]:Rand_Surf_win_n0[1] ] / math.sqrt( np.sum( Rand_data_np[ Rand_Surf_win_n0[0]:Rand_Surf_win_n0[1] ]**2, axis=0) )
                    
                    _, C_max_maxlag, shift=correlate_maxlag(Raw_data_np_surf,Rand_data_np_surf,Surf_maxlag,normalized=True)
                    # np2_zeros = zero_padding(Rand_data_np_surf,shift)
                    # Rand_data_np_surf_shift =  np2_zeros / math.sqrt( np.sum( np2_zeros**2, axis=0) )
                    Rand_data_np_surf_shift = Rand_data_np[ Rand_Surf_win_n0[0]-shift:Rand_Surf_win_n0[1]-shift ] / math.sqrt( np.sum( Rand_data_np[ Rand_Surf_win_n0[0]-shift:Rand_Surf_win_n0[1]-shift ]**2, axis=0) )
                    
                    MISFIT1_one_sta += Phase_weight[2]*math.sqrt(np.mean((Raw_data_np_surf-Rand_data_np_surf_shift)**2)) # *(1-C_max_maxlag)
                    # print(shift) 
                    tt = np.linspace(Raw_tsurf+Surf_win[0], Raw_tsurf+Surf_win[1], num=len(Raw_data_np_surf))
                    
                    if line_order==1:
                        str_text=plot_comp_name[2][i]
                        axs[line_order, colume_order].set_title(str_text,fontsize=20,fontweight='bold', alpha=0.7)
                            # axs[line_order, colume_order].text((tt[-1]+tt[0])/2,1*max_ylim,str_text,horizontalalignment='center', verticalalignment='bottom',
                            #                                    fontsize=18, color='black')
                    if Surf_inv_comp[i]=='yes':
                        jj+=1
                        p1,=axs[line_order, colume_order].plot(tt,Raw_data_np_surf, "k-",lw=2,alpha=0.6) 
                        p2,=axs[line_order, colume_order].plot(tt,Rand_data_np_surf_shift, linestyle='-', color='red', lw=3, alpha=0.6)
                        
                        cc=np.around(C_max_maxlag, decimals=2)
                        if jj ==1:
                            str_text=net_sta_name+"/"+str(np.around(distance, decimals=2)) 
                            axs[line_order, colume_order].text(tt[0], 0.9*max_surf_ylim,str_text,horizontalalignment='left', verticalalignment='top',
                                                               fontsize=15, color='black',bbox = dict(facecolor = "orange", alpha = 0.1))
                        cc_shift_text = str(cc)+'/'+str(shift*dt)
                        axs[line_order, colume_order].text(tt[0], -0.8*max_surf_ylim,cc_shift_text,horizontalalignment='left', verticalalignment='top',
                                                               fontsize=14, color='black')
                        
                        axs[line_order, colume_order].set_ylim(-max_surf_ylim,max_surf_ylim)
                        
                        axs[line_order,colume_order].set_axis_on()
                        axs[line_order,colume_order].spines['top'].set_visible(False)
                        axs[line_order,colume_order].spines['right'].set_visible(False) 
                        axs[line_order,colume_order].spines['bottom'].set_visible(True) 
                        axs[line_order,colume_order].spines['left'].set_visible(False)
                        axs[line_order, colume_order].get_xaxis().set_visible(True)
                        axs[line_order, colume_order].get_yaxis().set_visible(False)
                        # axs[line_order,colume_order].set_xlabel("Time(s)",fontsize=18)



    #%% ######
    # 5.3 Adjust the spacing of subgraphs to add suptitle
    fig5.subplots_adjust(wspace =0.3, hspace =0.5)
    # fig5.suptitle('Synthetic seismograms test',fontsize=30)


    #%% ######
    # 6 save figure
    figurename=os.path.join(plot_output_path,'waveform.'+fig_format)
    fig5.savefig(figurename,dpi=500, format=fig_format)






