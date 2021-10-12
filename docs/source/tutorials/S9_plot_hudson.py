#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from MCMTpy import MomentTensor as MTpy

# Helper function
#----------------------------------------------------#
#%% 4. plot_Hudson_points
def plot_Hudson_points(FM):
    fig1, ax1 = plt.subplots(1 ,1,  dpi=800)
    plt.rc('font',family='Times New Roman')
    MTpy.Hudson_plot(ax=ax1)
    for i in range(0,len(FM)): 
        MT = MTpy.MTensor(FM[i,:])
        M=MT.mt

        k,T = MTpy.M2kT_space(M)
        U,V = MTpy.kT2UV_space(k,T)
    
        map_vir = cm.get_cmap(name='YlGn')
        colors = map_vir(i/len(FM))
    
        ax1.scatter(U,V, color=colors,marker='o', s=5)


    position=fig1.add_axes([0.85, 0.15, 0.01, 0.5])
    font_colorbar = {'family' : 'Times New Roman',
            'color'  : 'black',
            'weight' : 'normal',
            'size'   : 6,
            }
    sm = cm.ScalarMappable(cmap=map_vir)               
    sm.set_array(np.arange(0,len(FM)+1))    
    cb=plt.colorbar(sm,cax=position,orientation='vertical') 
    cb.set_label('Sample',fontdict=font_colorbar)

    ax1.set_xlim(-4/3-0.1, 4/3+0.3)
    ax1.set_ylim(-1-0.1, 1+0.1)
    
    return fig1


#%%------------------------------------------------------------#
## 0. set path
example_path = '/Users/yf/3.Project/8.MCMTpy/MCMTpy-master/data/example_yunnan'


#------------------------------------------------------------#
# we expect no parameters need to be changed below
#------------------------------------------------------------#
FM_path=os.path.join(example_path,"YN.202105212148_Inv/mt_inv/Output_YN.202105212148_mt/rank_0_output/chain_0_FM_accept_all")
allfiles_path = os.path.join(example_path,'YN.202105212148_Inv/YN.202105212148_raw/*.SAC') 
## 1. read FM
N=6
FM_all = np.loadtxt(FM_path)
FM_raw = FM_all[0:,1:N+1]
Mxx_np = FM_raw[:,0]
Mxy_np = FM_raw[:,1]
Mxz_np = FM_raw[:,2]
Myy_np = FM_raw[:,3]
Myz_np = FM_raw[:,4]
Mzz_np = FM_raw[:,5]
FM = np.vstack((Mxx_np,  Myy_np,  Mzz_np,  Mxy_np,  Mxz_np,  Myz_np)).T

## 2. plot Hudson
fig = plot_Hudson_points(FM)
figurename=os.path.join('./S2_figure/Hudson.pdf')
plt.savefig(figurename,dpi=800, format="pdf")