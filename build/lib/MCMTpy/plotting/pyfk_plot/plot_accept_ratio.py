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
def plot_accept_ratio(plot_output_path, accept_ratio_all, fig_format):
    

    fig, axs = plt.subplots(1, 1,figsize=(12, 5))

    xx = np.arange(0, accept_ratio_all.shape[0], 1)
    axs.plot(xx, accept_ratio_all, linestyle='-', color='k', alpha=0.6,label="accept ratio") 

    axs.legend(loc='best',fontsize=12, shadow=False)
    axs.set_xlabel("Sample Number",fontsize=15)
    axs.set_ylabel("Accept Ratio",fontsize=15)
    # axs[0].set_xscale("log")
    # axs[0].set_xlim(0,20000)
    # axs[0].set_ylim(0,1.1)
    # axs[0].get_xaxis().set_visible(False)

 
    
    ##################
    figurename=os.path.join(plot_output_path,'accept_ratio.'+fig_format)
    fig.savefig(figurename,dpi=800, format=fig_format)







