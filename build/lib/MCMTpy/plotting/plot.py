#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 21 00:20:35 2021

@author: Fu Yin (yinfu@mail.ustc.edu.cn) at USTC

This script:
    1) plot hist
    2) plot misfit
    3) plot wavaform



Modify history:
    1) Mar 21 00:20:35 2021   ||    Fu Yin at USTC    ||    The initial release.
    2) Jun 28 11:08:02 2021   ||    Fu Yin at USTC    ||    Add figure format parameter.
"""


import numpy as np
import os,sys
import json5
from MCMTpy.plotting.pyfk_plot.plot_hist import plot_hist
from MCMTpy.plotting.pyfk_plot.plot_misfit import plot_misfit
from MCMTpy.plotting.pyfk_plot.plot_waveform import plot_waveform
from MCMTpy.plotting.pyfk_plot.plot_accept_ratio import plot_accept_ratio
from MCMTpy.gfs.pyfk_GFs.read_json import read_GFs_json
from MCMTpy.sampler.pyfk_MH.read_sampler_json import read_Inv_json




def plot(filename,method):
    if method=='pyfk':
        with open(filename, 'r',encoding='utf-8') as f1:
            plot_json_r = json5.load(f1)
        plot_json = plot_json_r.copy()
        
        Inv_para,Sta_inv_info,Sta_raw_data,Stf = read_Inv_json(plot_json["Inv_json_file"])
        Prepro_para,_,Source_Station_info_MPI = read_GFs_json(Inv_para["GFs_json_file"])

        sample_path      = Inv_para["Output_path"]
        plot_output_path = plot_json["plot_output_path"]                        # The path of the output file
        fig_format       = plot_json["fig_format"]                              # history: Jun 28 11:08:02 2021

        if not os.path.isdir(plot_output_path):os.mkdir(plot_output_path)



        #### 1.plot_hist
        if plot_json["plot_hist"]==True:
            print('Plotting hist figure:\n')
            MPI_n          = Inv_para["MPI_n"]
            Chains_n       = Inv_para["Chains_n"]
            num_bins       = plot_json["num_bins"]                                    # plot hist2d: the number of grid
            num_std        = plot_json["num_std"]                                       # range of axes in each subgraph (mean +- several times standard deviation)
            labels_name    = plot_json["labels_name"]                                # labels:['mw','strike/°','dip/°','rake/°','x/km','y/km','z/km','t0/s']
            N_start        = plot_json["N_start"]
            N_start_accept = plot_json["N_start_accept"]
            Fixed_FM       = Inv_para["Fixed_FM"]
            InvType        = Inv_para["InvType"]                                       # 'mt' 'dc' 'sf' 'ep'
            N_k            = Inv_para["N_k"]
            
            for i in range(0,MPI_n,1):
                rank_path = os.path.join(sample_path,'rank_'+str(i)+'_output')
                for j in range(0,Chains_n,1):
                    FM_2_all_path = os.path.join(rank_path,'chain_'+str(j)+'_FM_2_all')
                    FM_2_accept_path = os.path.join(rank_path,'chain_'+str(j)+'_FM_2_accept_all')
                    if i==0 and j==0:
                        FM_2_all        = np.loadtxt(FM_2_all_path)[N_start:,:]
                        FM_2_accept_all = np.loadtxt(FM_2_accept_path)[N_start_accept:,:]
                    else:
                        FM_2_all        = np.vstack( ( FM_2_all, np.loadtxt(FM_2_all_path)[N_start:,:] ) )
                        FM_2_accept_all = np.vstack( ( FM_2_accept_all, np.loadtxt(FM_2_accept_path)[N_start_accept:,:] ) )
                    
            
            plot_hist(plot_output_path,FM_2_all,FM_2_accept_all,num_bins,num_std,labels_name,Fixed_FM,InvType,fig_format)




        #### 2.plot_misfit
        if plot_json["plot_misfit"]==True:
            print('Plotting misfit figure:\n')
            MPI_n_st           = plot_json["MPI_n_mis"]
            Chains_n_st        = plot_json["Chains_n_mis"]
            rank_path          = os.path.join(sample_path,'rank_'+str(MPI_n_st)+'_output')

            MISFIT_1_all_path    = os.path.join(rank_path,'chain_'+str(Chains_n_st)+'_MISFIT_1_all')
            MISFIT_1_accept_path = os.path.join(rank_path,'chain_'+str(Chains_n_st)+'_MISFIT_1_accept_all')
            MISFIT_2_all_path    = os.path.join(rank_path,'chain_'+str(Chains_n_st)+'_MISFIT_2_all')
            MISFIT_2_accept_path = os.path.join(rank_path,'chain_'+str(Chains_n_st)+'_MISFIT_2_accept_all')
            
            MISFIT_1_all         = np.loadtxt(MISFIT_1_all_path)  
            MISFIT_1_accept_all  = np.loadtxt(MISFIT_1_accept_path) 
            MISFIT_2_all         = np.loadtxt(MISFIT_2_all_path)  
            MISFIT_2_accept_all  = np.loadtxt(MISFIT_2_accept_path) 
            
            plot_misfit(plot_output_path,MISFIT_1_all,MISFIT_1_accept_all,MISFIT_2_all,MISFIT_2_accept_all,fig_format)


        #### 3.plot_accept_ratio
        if plot_json["plot_accept_ratio"]==True:
            print('Plotting accept_ratio figure:\n')
            MPI_n_st           = plot_json["MPI_n_acc"]
            Chains_n_st        = plot_json["Chains_n_acc"]
            rank_path          = os.path.join(sample_path,'rank_'+str(MPI_n_st)+'_output')

            accept_ratio_all_path    = os.path.join(rank_path,'chain_'+str(Chains_n_st)+'_accpet_ratio_all')
            accept_ratio_all         = np.loadtxt(accept_ratio_all_path)  

            
            plot_accept_ratio(plot_output_path,accept_ratio_all,fig_format)



        #### 4.plot_waveform
        if plot_json["plot_waveform"]==True:
            print('Plotting waveform figure:\n')
            Inv_json_file  = plot_json["Inv_json_file"]
            NET_STA        = plot_json["NET_STA"]
            FM_best        = plot_json["FM_best"]
            line_n_sta     = plot_json["line_n_sta"]
            max_p_ylim     = plot_json["max_p_ylim"]
            max_s_ylim     = plot_json["max_s_ylim"]
            max_surf_ylim  = plot_json["max_surf_ylim"]
            plot_comp      = plot_json["plot_comp"]
            plot_comp_name = plot_json["plot_comp_name"]
            Geometrical_spreading_mode = plot_json["Geometrical_spreading_mode"]
            Amplitude_normalization_mode = plot_json["Amplitude_normalization_mode"]

            plot_waveform(Inv_json_file,plot_output_path,NET_STA,FM_best,line_n_sta,\
                          max_p_ylim,max_s_ylim,max_surf_ylim,\
                              plot_comp,plot_comp_name,fig_format,Geometrical_spreading_mode,Amplitude_normalization_mode)




            
            
            
    elif method=='sem':
        pass
    else:
        sys.exit('MCMTpy: command error')















