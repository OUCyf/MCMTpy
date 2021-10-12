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

import matplotlib.pyplot as plt
import numpy as np
import os



###################
def plot_misfit(plot_output_path,MISFIT_1_all,MISFIT_1_accept_all,MISFIT_2_all,MISFIT_2_accept_all,fig_format):
    

    Misfit_1_all = MISFIT_1_all
    Misfit_2_all = MISFIT_2_all

    if MISFIT_1_accept_all.any()==True:
        Misfit_accept_1_all = MISFIT_1_accept_all
    else:
        Misfit_accept_1_all = np.array([])

    if MISFIT_1_accept_all.any()==True:
        Misfit_accept_2_all = MISFIT_2_accept_all
    else:
        Misfit_accept_2_all = np.array([])



    #%% 1. plot MISFIT_all curve
    fig, axs = plt.subplots(2, 1, dpi=800,gridspec_kw=dict(height_ratios=[1, 1]),figsize=(12, 6))
    iteration_num_1 = Misfit_1_all.shape[0]
    iteration_num_2 = Misfit_2_all.shape[0]
    iteration_num=Misfit_1_all.shape[0] + Misfit_2_all.shape[0]
    tt=np.arange(0,iteration_num,1)

    #################
    p1,=axs[0].plot(tt[0:iteration_num_1],
                    Misfit_1_all/max(Misfit_1_all),  
                    linestyle='-', 
                    color='k', 
                    alpha=0.6,
                    label="arrive time misift") 
    axs[0].legend(loc='best',fontsize=10, shadow=False)
    axs[0].set_xlim(-1e0,iteration_num+1)
    (1-min(Misfit_1_all)/max(Misfit_1_all))*0.2
    axs[0].set_ylim(min(Misfit_1_all)/max(Misfit_1_all), 1.2)

    #################
    p2,=axs[1].plot(tt[iteration_num_1:iteration_num_1+iteration_num_2],
                    Misfit_2_all/max(Misfit_2_all),  
                    linestyle='-', 
                    color='k', 
                    alpha=0.6,
                    label="waveform misfit") 
    axs[1].legend(loc='best',fontsize=10, shadow=False)
    axs[1].set_xlim(-1e0,iteration_num+1)
    y2_min = min(Misfit_2_all)/max(Misfit_2_all) - 0.2*(1-min(Misfit_2_all)/max(Misfit_2_all))
    axs[1].set_ylim(y2_min,1)

    # ###################
    axs[0].get_xaxis().set_visible(False)
    fig.text(0.5, 0.02, 'Sample number', ha='center',fontsize=20)
    fig.text(0.05, 0.5, 'S(m)', va='center', rotation='vertical',fontsize=20)
    fig.subplots_adjust(wspace =0.1, hspace =0)
    # fig.suptitle('Misfit with iteration',fontsize=20)
    
    ##################
    figurename=os.path.join(plot_output_path,'misfit.'+fig_format)
    fig.savefig(figurename,dpi=800, format=fig_format)




    #%% 2. plot MISFIT_accept_all curve
    fig, axs = plt.subplots(2, 1, dpi=800,gridspec_kw=dict(height_ratios=[1, 1]),figsize=(12, 6))
    iteration_num_1 = Misfit_accept_1_all.shape[0]
    iteration_num_2 = Misfit_accept_2_all.shape[0]
    iteration_num=Misfit_accept_1_all.shape[0] + Misfit_accept_2_all.shape[0]
    tt=np.arange(0,iteration_num,1)

    #################
    p1,=axs[0].plot(tt[0:iteration_num_1],
                    Misfit_accept_1_all/max(Misfit_accept_1_all),  
                    linestyle='-', 
                    color='k', 
                    alpha=0.6,
                    label="arrive time misift") 
    axs[0].legend(loc='best',fontsize=10, shadow=False)
    axs[0].set_xlim(-1e0,iteration_num+1)
    (1-min(Misfit_accept_1_all)/max(Misfit_accept_1_all))*0.2
    axs[0].set_ylim(min(Misfit_accept_1_all)/max(Misfit_accept_1_all), 1.2)

    #################
    p2,=axs[1].plot(tt[iteration_num_1:iteration_num_1+iteration_num_2],
                    Misfit_accept_2_all/max(Misfit_accept_2_all),  
                    linestyle='-', 
                    color='k', 
                    alpha=0.6,
                    label="waveform misfit") 
    axs[1].legend(loc='best',fontsize=10, shadow=False)
    axs[1].set_xlim(-1e0,iteration_num+1)
    y2_min = min(Misfit_accept_2_all)/max(Misfit_accept_2_all) - 0.2*(1-min(Misfit_accept_2_all)/max(Misfit_accept_2_all))
    axs[1].set_ylim(y2_min,1)

    # ###################
    axs[0].get_xaxis().set_visible(False)
    fig.text(0.5, 0.02, 'Sample number', ha='center',fontsize=20)
    fig.text(0.05, 0.5, 'S(m)', va='center', rotation='vertical',fontsize=20)
    fig.subplots_adjust(wspace =0.1, hspace =0)
    # fig.suptitle('Misfit with iteration',fontsize=20)
    
    ##################
    figurename=os.path.join(plot_output_path,'misfit_accept.'+fig_format)
    fig.savefig(figurename,dpi=800, format=fig_format)



