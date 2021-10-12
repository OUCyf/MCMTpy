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
#%% 1.plot decompose mt
def plot_decompose(FM,FM_DC,FM_CLVD,FM_iso):
    MT = MTpy.MTensor(FM)
    Dec = MTpy.Decompose(MT)
    Dec.decomposition_iso_DC_CLVD()
    
    fig2, ax2 = plt.subplots(1 ,1,  dpi=800)
    Length_Ball = 100

    ###### MT
    beach1 = beach(FM, xy=(50, 50),  linewidth=1,width=Length_Ball-1, alpha=1,\
                   facecolor='g',bgcolor='w', edgecolor='k',mopad_basis='NED',nofill=False,zorder=1 )
    ax2.add_collection(beach1) 
    ax2.set_aspect("equal")

    menthod='schmidt'   # 'schmidt'  # 'wulff'
    T_axis, P_axis, N_axis = MTpy.MT2TPN(MTpy.MTensor(FM))
    T = MTpy.vector2str_dip(T_axis)
    P = MTpy.vector2str_dip(P_axis)
    N = MTpy.vector2str_dip(N_axis)

    Tx, Ty = MTpy.project_beachball(AZM=T.strike, TKO=(90-T.dip), R=Length_Ball/2, menthod=menthod)
    ax2.text(Tx,Ty,'T',horizontalalignment='center', verticalalignment='center',\
              fontsize=10, color='k',alpha=0.7,zorder=1) 

    Px, Py = MTpy.project_beachball(AZM=P.strike, TKO=(90-P.dip), R=Length_Ball/2, menthod=menthod)
    ax2.text(Px,Py,'P',horizontalalignment='center', verticalalignment='center',\
              fontsize=10, color='k',alpha=0.7,zorder=1) 

    Nx, Ny = MTpy.project_beachball(AZM=N.strike, TKO=(90-N.dip), R=Length_Ball/2, menthod=menthod)
    ax2.text(Nx,Ny,'N',horizontalalignment='center', verticalalignment='center',\
              fontsize=10, color='k',alpha=0.7,zorder=1) 

    ###### DC
    beach2 = beach(FM_DC, xy=(200, 50),  linewidth=1,width=Length_Ball-1, alpha=1,\
                   facecolor='g',bgcolor='w', edgecolor='k',mopad_basis='NED',nofill=False,zorder=1 )
    ax2.add_collection(beach2) 
    ax2.set_aspect("equal")

    menthod='schmidt'   # 'schmidt'  # 'wulff'
    T_axis, P_axis, N_axis = MTpy.MT2TPN(MTpy.MTensor(FM_DC))
    T = MTpy.vector2str_dip(T_axis)
    P = MTpy.vector2str_dip(P_axis)
    N = MTpy.vector2str_dip(N_axis)

    Tx, Ty = MTpy.project_beachball(AZM=T.strike, TKO=(90-T.dip), R=Length_Ball/2, menthod=menthod)
    ax2.text(150*1+Tx,Ty,'T',horizontalalignment='center', verticalalignment='center',\
              fontsize=10, color='k',alpha=0.7,zorder=1) 

    Px, Py = MTpy.project_beachball(AZM=P.strike, TKO=(90-P.dip), R=Length_Ball/2, menthod=menthod)
    ax2.text(150*1+Px,Py,'P',horizontalalignment='center', verticalalignment='center',\
              fontsize=10, color='k',alpha=0.7,zorder=1) 

    Nx, Ny = MTpy.project_beachball(AZM=N.strike, TKO=(90-N.dip), R=Length_Ball/2, menthod=menthod)
    ax2.text(150*1+Nx,Ny,'N',horizontalalignment='center', verticalalignment='center',\
              fontsize=10, color='k',alpha=0.7,zorder=1) 

    ####### CLVD
    beach3 = beach(FM_CLVD, xy=(350, 50),  linewidth=1,width=Length_Ball-1, alpha=1,\
                   facecolor='g',bgcolor='w', edgecolor='k',mopad_basis='NED',nofill=False,zorder=1 )
    ax2.add_collection(beach3) 
    ax2.set_aspect("equal")

    menthod='schmidt'   # 'schmidt'  # 'wulff'
    T_axis, P_axis, N_axis = MTpy.MT2TPN(MTpy.MTensor(FM_CLVD))
    T = MTpy.vector2str_dip(T_axis)
    P = MTpy.vector2str_dip(P_axis)
    N = MTpy.vector2str_dip(N_axis)

    Tx, Ty = MTpy.project_beachball(AZM=T.strike, TKO=(90-T.dip), R=Length_Ball/2, menthod=menthod)
    ax2.text(150*2+Tx,Ty,'T',horizontalalignment='center', verticalalignment='center',\
              fontsize=10, color='k',alpha=0.7,zorder=1) 

    Px, Py = MTpy.project_beachball(AZM=P.strike, TKO=(90-P.dip), R=Length_Ball/2, menthod=menthod)
    ax2.text(150*2+Px,Py,'P',horizontalalignment='center', verticalalignment='center',\
              fontsize=10, color='k',alpha=0.7,zorder=1) 

    Nx, Ny = MTpy.project_beachball(AZM=N.strike, TKO=(90-N.dip), R=Length_Ball/2, menthod=menthod)
    ax2.text(150*2+Nx,Ny,'N',horizontalalignment='center', verticalalignment='center',\
              fontsize=10, color='k',alpha=0.7,zorder=1) 

    ####### iso
    beach4 = beach(FM_iso, xy=(500, 50),  linewidth=1,width=Length_Ball-1, alpha=1,\
                   facecolor='g',bgcolor='w', edgecolor='k',mopad_basis='NED',nofill=False,zorder=1 )
    ax2.add_collection(beach4) 
    ax2.set_aspect("equal")

    menthod='schmidt'   # 'schmidt'  # 'wulff'
    T_axis, P_axis, N_axis = MTpy.MT2TPN(MTpy.MTensor(FM_iso))
    T = MTpy.vector2str_dip(T_axis)
    P = MTpy.vector2str_dip(P_axis)
    N = MTpy.vector2str_dip(N_axis)

    Tx, Ty = MTpy.project_beachball(AZM=T.strike, TKO=(90-T.dip), R=Length_Ball/2, menthod=menthod)
    ax2.text(150*3+Tx,Ty,'T',horizontalalignment='center', verticalalignment='center',\
              fontsize=10, color='k',alpha=0.7,zorder=1) 

    Px, Py = MTpy.project_beachball(AZM=P.strike, TKO=(90-P.dip), R=Length_Ball/2, menthod=menthod)
    ax2.text(150*3+Px,Py,'P',horizontalalignment='center', verticalalignment='center',\
              fontsize=10, color='k',alpha=0.7,zorder=1) 

    Nx, Ny = MTpy.project_beachball(AZM=N.strike, TKO=(90-N.dip), R=Length_Ball/2, menthod=menthod)
    ax2.text(150*3+Nx,Ny,'N',horizontalalignment='center', verticalalignment='center',\
              fontsize=10, color='k',alpha=0.7,zorder=1)

    ####### plot '+' and '='
    ax2.text(125,50,'=',horizontalalignment='center', verticalalignment='center',\
              fontsize=20, color='k',alpha=0.7,zorder=1)
    ax2.text(275,50,'+',horizontalalignment='center', verticalalignment='center',\
              fontsize=20, color='k',alpha=0.7,zorder=1)
    ax2.text(425,50,'+',horizontalalignment='center', verticalalignment='center',\
              fontsize=20, color='k',alpha=0.7,zorder=1)
    
    ####### plot percentage
    MT_text = 'MT'
    ax2.text(50,-15,MT_text,horizontalalignment='center', verticalalignment='center',\
              fontsize=8, color='k',alpha=0.7,zorder=1)
    
    iso_text = 'ISO: '+str(round(Dec.M_iso_percentage,2))+'%'
    ax2.text(500,-15,iso_text,horizontalalignment='center', verticalalignment='center',\
              fontsize=8, color='k',alpha=0.7,zorder=1)
    
    DC_text = 'DC: '+str(round(Dec.M_DC_percentage,2))+'%'
    ax2.text(200,-15,DC_text,horizontalalignment='center', verticalalignment='center',\
              fontsize=8, color='k',alpha=0.7,zorder=1)
    
    CLVD_text = 'CLVD: '+str(round(Dec.M_CLVD_percentage,2))+'%'
    ax2.text(350,-15,CLVD_text,horizontalalignment='center', verticalalignment='center',\
              fontsize=8, color='k',alpha=0.7,zorder=1)

    ####### set figure
    plt.title("Decompose")
    ax2.set_xlim(0,550)
    ax2.set_ylim(-30,100)
    ax2.set_axis_off() 

    return fig2

#%%------------------------------------------------------------#
# 1.set FM
# [Mxx_np,  Myy_np,  Mzz_np,  Mxy_np,  Mxz_np,  Myz_np] or [strike,dip,rake]
FM=[150,50,100]


#------------------------------------------------------------#
# we expect no parameters need to be changed below
#------------------------------------------------------------#
MT = MTpy.MTensor(FM)
Dec = MTpy.Decompose(MT)
Dec.decomposition_iso_DC_CLVD()   # Dec.help()
print(Dec.help())
print('\n\n*********************************\n\n')
Dec.print_self()

Mxx = Dec.M_iso[0,0]
Mxy = Dec.M_iso[0,1]
Mxz = Dec.M_iso[0,2]
Myy = Dec.M_iso[1,1]
Myz = Dec.M_iso[1,2]
Mzz = Dec.M_iso[2,2]
FM_iso = np.array((Mxx,  Myy,  Mzz,  Mxy,  Mxz,  Myz))

Mxx = Dec.M_DC[0,0]
Mxy = Dec.M_DC[0,1]
Mxz = Dec.M_DC[0,2]
Myy = Dec.M_DC[1,1]
Myz = Dec.M_DC[1,2]
Mzz = Dec.M_DC[2,2]
FM_DC = np.array((Mxx,  Myy,  Mzz,  Mxy,  Mxz,  Myz))

Mxx = Dec.M_CLVD[0,0]
Mxy = Dec.M_CLVD[0,1]
Mxz = Dec.M_CLVD[0,2]
Myy = Dec.M_CLVD[1,1]
Myz = Dec.M_CLVD[1,2]
Mzz = Dec.M_CLVD[2,2]
FM_CLVD = np.array((Mxx,  Myy,  Mzz,  Mxy,  Mxz,  Myz))

fig = plot_decompose(FM, FM_DC, FM_CLVD, FM_iso)
figurename=os.path.join('./S2_figure/Decompose.pdf')
plt.savefig(figurename,dpi=800, format="pdf")