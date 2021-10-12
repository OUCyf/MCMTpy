#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  7 15:04:31 2021

@author: Fu Yin (yinfu@mail.ustc.edu.cn) at USTC

This script:
    1) function rotate.


Modify history:
    1) Jun  7 15:04:31 2021    ||    Fu Yin at USTC    ||    The initial release.
    2) ...
    
"""

from math import sin,cos,pi




############################################
def rotate(data_raw,chn=['RTZ','ENZ']):
    '''
    Parameters
    ----------
    data : TYPE
        obspy trace.
    chn : TYPE, optional
        The default is ['RTZ','ENZ']. rotate data form chn[0] to chn[1]

    Raises
    ------
    ValueError
        DESCRIPTION.

    Returns
    -------
    data : TYPE
        DESCRIPTION.

    '''
    data = data_raw.copy()
    if ('R' in chn[0]) and ('T' in chn[0]) and ('Z' in chn[0]) and ('E' in chn[1]) and ('N' in chn[1]) and ('Z' in chn[1]):
        r_in_index = chn[0].index('R')
        t_in_index = chn[0].index('T')
        z_in_index = chn[0].index('Z')

        e_out_index = chn[1].index('E')
        n_out_index = chn[1].index('N')
        z_out_index = chn[1].index('Z')

        for i in range(0,len(data),3):
            BAZ = data[i].stats.back_azimuth/180*pi
            R = data[i+r_in_index].data
            T = data[i+t_in_index].data
            Z = data[i+z_in_index].data
            E = -cos(BAZ)*T - sin(BAZ)*R
            N = sin(BAZ)*T - cos(BAZ)*R

            data[i+e_out_index].stats.channel='E'
            data[i+e_out_index].data = E
            data[i+n_out_index].stats.channel='N'
            data[i+n_out_index].data = N
            data[i+z_out_index].stats.channel='Z'
            data[i+z_out_index].data = Z

    elif ('E' in chn[0]) and ('N' in chn[0]) and ('Z' in chn[0]) and ('R' in chn[1]) and ('T' in chn[1]) and ('Z' in chn[1]):
        e_in_index = chn[0].index('E')
        n_in_index = chn[0].index('N')
        z_in_index = chn[0].index('Z')

        r_out_index = chn[1].index('R')
        t_out_index = chn[1].index('T')
        z_out_index = chn[1].index('Z')

        for i in range(0,len(data),3):
            BAZ = data[i].stats.back_azimuth/180*pi
            E = data[i+e_in_index].data
            N = data[i+n_in_index].data
            Z = data[i+z_in_index].data
            R = -cos(BAZ)*N - sin(BAZ)*E
            T = sin(BAZ)*N - cos(BAZ)*E

            data[i+r_out_index].stats.channel='R'
            data[i+r_out_index].data = R
            data[i+t_out_index].stats.channel='T'
            data[i+t_out_index].data = T
            data[i+z_out_index].stats.channel='Z'
            data[i+z_out_index].data = Z

    else:
        raise ValueError('chn format error! abort!')



    return data




