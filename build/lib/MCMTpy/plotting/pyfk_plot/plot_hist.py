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
from math import sqrt




##############################################################################
def plot_hist(plot_output_path,FM_all,FM_accept_all,num_bins,num_std,labels_name,Fixed_FM,InvType,fig_format):

    for ff in range(0,2,1):
        if ff==0:
            FM = FM_all
            figurename=os.path.join(plot_output_path,'hist.'+fig_format)
        else:
            FM = FM_accept_all
            figurename=os.path.join(plot_output_path,'hist_accept.'+fig_format)


        if InvType == 'mt':
            N_all = 11
        elif InvType == 'dc':
            N_all = 8
        elif InvType == 'sf':
            N_all = 8
        elif InvType == 'mt':
            N_all = 5


        N = Fixed_FM.count('variable')                                          # plot: The number of subgraphs of a row

        delete_list=[]
        for i in range(0,N_all,1):
            if Fixed_FM[i]!='variable':
                delete_list.append(i)

        FM = np.delete(FM,delete_list, axis=1)                                  # The delete solution that is a fixed parameter



        FM_mean=np.zeros(shape=(N))
        FM_sigma=np.zeros(shape=(N))
        for i in range(0,N,1):
            FM_mean[i]=np.mean(FM[0:,i])                                        # FM The mean of each parameter
            FM_sigma[i]=np.std(FM[0:,i])                                        # FM The standard deviation of each parameter

        FM_used = FM.copy()                                                     # strike rake cycle skip
        if InvType=='dc':
            if FM_mean[1]>360 and FM_mean[1]<720:
                FM_used[:,1]=FM_used[:,1]-360
                FM_mean[1]=FM_mean[1]-360
            if FM_mean[1]<0 and FM_mean[1]>(-360):
                FM_used[:,1]=FM_used[:,1]+360
                FM_mean[1]=FM_mean[1]+360

            if FM_mean[3]>180 and FM_mean[3]<540:
                FM_used[:,3]=FM_used[:,3]-360
                FM_mean[3]=FM_mean[3]-360
            if FM_mean[3]<(-180) and FM_mean[3]>(-540):
                FM_used[:,3]=FM_used[:,3]+360
                FM_mean[3]=FM_mean[3]+360

        plt.style.use('default')                                                # Use the default plotting style (make sure it's not the ggplot style)
        fig, axs = plt.subplots(N, N, dpi=800,figsize=(16, 16))

        for i in range(0,N,1):
            for j in range(0,N,1):
                x=FM_used[:,i]
                y=FM_used[:,j]
                axs[i,j].get_xaxis().set_visible(False)                         # Do not display the axes
                axs[i,j].get_yaxis().set_visible(False)
                if j==0:
                    axs[i,j].get_yaxis().set_visible(True)
                    axs[i,j].set_ylabel(labels_name[i],fontsize=19)             # The left axis labels
                if i==N-1:
                    axs[i,j].get_xaxis().set_visible(True)
                    axs[i,j].set_xlabel(labels_name[j],fontsize=19)             # The bottom axis labels
                
                if i<j:
                    axs[i,j].set_axis_off()                                     # The upper right part of the image is not displayed
                
                if i>j:
                    h, xx, yy, p=axs[i,j].hist2d( y,x, bins=(50, 50), cmap=plt.cm.gist_earth_r)   # hist2d
                    x_start = FM_mean[j]-num_std*FM_sigma[j]                    # Adjust the range of the coordinate axes
                    x_end = FM_mean[j]+num_std*FM_sigma[j]
                    y_start = FM_mean[i]-num_std*FM_sigma[i]
                    y_end = FM_mean[i]+num_std*FM_sigma[i]
                    axs[i, j].set_xlim(x_start,x_end)
                    axs[i, j].set_ylim(y_start,y_end)
                    axs[i,j].scatter(FM_mean[j],FM_mean[i],s=100,c='tab:red')   # Draw the average value (red dots)
                    # plt.clf()
                    # axs[i,j].imshow(h, origin = "lower", interpolation = "gaussian",cmap=plt.cm.gist_earth_r)
                    cov_yx = np.cov(y.T, x.T)                                   # Covariance matrix
                    cov0_yx = np.around(cov_yx[0,1]/sqrt(cov_yx[0,0]*cov_yx[1,1]), decimals=2)  # The correlation coefficient
                    
                    cov0_yx_text = 'cov0: '+str(cov0_yx)
                    axs[i,j].text(x_start, y_end,cov0_yx_text,horizontalalignment='left', verticalalignment='top',
                                  fontsize=12, color='black',bbox = dict(facecolor = "k", alpha = 0.05))                      # Mark the correlation coefficient SHIFT
                        
                    
                if i==j:
                    axs[i,j].spines['top'].set_visible(False)                   # Coordinate boxes are not displayed
                    axs[i,j].spines['right'].set_visible(False) 
                    axs[i,j].spines['left'].set_visible(False)
                    axs[i,j].spines['bottom'].set_visible(False)
                    axs[i,j].get_yaxis().set_visible(False)
            

                    n, bins, patches = axs[i,j].hist(x, num_bins, density=True,histtype='stepfilled', facecolor='orange',
                                                     alpha=0.5)                 # histgram
        
                    mu=np.around(FM_mean[i], decimals=2)                        # Keep two decimal places
                    sigma =np.around(FM_sigma[i], decimals=2)
                
                    y = ((1 / (np.sqrt(2 * np.pi) * sigma)) *
                         np.exp(-0.5 * (1 / sigma * (bins - mu))**2))
                    axs[i,j].plot(bins, y, '--',color='g')                      # Draw a Gaussian fitting curve
                    str_hist=' $\mu='+str(mu)+'$, $\sigma='+str(sigma)+'$'
                    axs[i,j].set_title(r' $\mu='+str(mu)+'$, $\sigma='+str(sigma)+'$')
                
                    x_start = FM_mean[j]-num_std*FM_sigma[j]                    # Adjust the range of the coordinate axes
                    x_end = FM_mean[j]+num_std*FM_sigma[j]
                    y_start = FM_mean[i]-num_std*FM_sigma[i]
                    y_end = FM_mean[i]+num_std*FM_sigma[i]
                    axs[i,j].set_xlim(x_start,x_end)
                    axs[i,j].spines['bottom'].set_visible(True)                 # Shows the axis of coordinates at the bottom

        fig.subplots_adjust(wspace =0.1, hspace =0.1)
        fig.savefig(figurename,dpi=800, format=fig_format)


##############################################################################





