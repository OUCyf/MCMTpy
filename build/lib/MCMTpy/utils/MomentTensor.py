#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 23 19:53:22 2021
@author: Fu Yin (yinfu@mail.ustc.edu.cn) at USTC

This script refers to some theories and codes of Obspy/MoPaD/Introduction to Seismology (Yongge Wan):
    1) The MoPaD program (https://github.com/geophysics/MoPaD);
    2) str_dip_rake to mt;
    3) The conversion between str_dip_rake and A/N vector;
    4) The conversion between A/N vector and P/T/N vector;
    5) mt to P/T/N vector;
    6) P/T/N vector to P/T/N vector's stirke and dip;
    7) project station to beachball;
    8) Hudson plot
    9) Decompose of mt: isotropic + deviatoric = isotropic + DC + CLVD;
   10) Describe_fault_plane with two str_dip_rake.

Modify history:
    1) May 23 19:53:22 2021    ||    Fu Yin at USTC    ||    The initial release.
    2) ...
    
"""


import numpy as np
import math
from math import pi
from math import sin,cos,tan,atan2,atan,sqrt,acos



#%%######################################################
#               (2) str_dip_rake to mt
#########################################################

def str_dip_rake2MT(strike,dip,rake):
    """
    Input: fault plane' strike dip and rake in degrees.
            strike : [0, 360)
            dip    : [0, 90]
            rake   : [-180, 180)

    Output: a moment tensor object in NED system.
    """
    strike = strike/180*pi
    dip = dip/180*pi
    rake = rake/180*pi
    M0 = 1

    Mxx = -M0*( sin(dip) * cos(rake) * sin(2*strike) + sin(2*dip) * sin(rake) * sin(strike)**2 )
    Myy =  M0*( sin(dip) * cos(rake) * sin(2*strike) - sin(2*dip) * sin(rake) * cos(strike)**2 )
    Mzz =  M0*( sin(2*dip) * sin(rake) )
    Mxy =  M0*( sin(dip) * cos(rake) * cos(2*strike) + 1/2* sin(2*dip) * sin(rake) * sin(2*strike) )
    Mxz = -M0*( cos(dip) * cos(rake) * cos(strike) + cos(2*dip) * sin(rake) * sin(strike) )
    Myz = -M0*( cos(dip) * cos(rake) * sin(strike) - cos(2*dip) * sin(rake) * cos(strike) )

    A = MTensor([Mxx, Myy, Mzz, Mxy, Mxz, Myz])

    return A





#%%###################################################################
#   (3) The conversion between str_dip_rake and A/N vector
######################################################################

# str_dip_rake to A/N vector
def str_dip_rake2AN(strike,dip,rake):
    """
    Input: fault plane' strike dip and rake in degrees.
            strike : [0, 360)
            dip    : [0, 90]
            rake   : [-180, 180)

    Output: slip vector(A) and fault plane's normal vector(N) in NED system.
    """
    strike = strike/180*pi
    dip = dip/180*pi
    rake = rake/180*pi

    A=np.array([  cos(rake)*cos(strike) + sin(rake)*cos(dip)*sin(strike),
                  cos(rake)*sin(strike) - sin(rake)*cos(dip)*cos(strike), 
                 -sin(rake)*sin(dip)] )

    N=np.array([ -sin(strike)*sin(dip),
                  cos(strike)*sin(dip),
                 -cos(dip)] )

    return A,N

# A/N vector to str_dip_rake
def AN2str_dip_rake(A,N):
    """
    Input: slip vector(A) and fault plane's normal vector(N) in NED system.

    Output: fault plane' strike dip and rake.
            strike : [0, 360)
            dip    : [0, 90]
            rake   : [-180, 180)
    """

    if abs(N[2]+1) < 0.00001:                                                   # nz=-1: the fault plane is horizontal
        strike = atan2(A[1],A[0])                                               # The direction of slip is also the strike, because the fault plane is horizontal
        dip = 0.0
    else:
        strike = atan2(-N[0],N[1])
        if abs(N[2]-0) < 0.00001:                                               # nz=-1: the fault plane is vertical
            dip = pi/2
        elif abs(sin(strike)) > abs(cos(strike)):
            dip = atan( (N[0]/sin(strike)) / N[2] )
        else:
            dip = atan( (-N[1]/cos(strike)) / N[2] )

    cos_rake = A[0]*cos(strike) + A[1]*sin(strike)

    if abs(A[2]-0) > 0.0000001:                                                 # az!=0: consider the effect of dip
        if abs(dip-0) > 0.000001:
            rake = atan2(-A[2]/sin(dip),cos_rake)
        else:
            rake = atan2(-100000000.0*A[2],cos_rake)
    else:                                                                       # az=0: don't consider the effect of dip
        if cos_rake > 1:
            cos_rake = 1
        if cos_rake < -1:
            cos_rake = -1
        rake = acos(cos_rake)

    if dip < 0:
        dip = -dip
        strike = strike+pi                                                      # strike need to be in the opposite direction

    if strike >= 2*pi:
        strike = strike-2*pi
    if strike < 0:
        strike = strike+2*pi

    strike = strike*180/pi
    dip = dip*180/pi
    rake = rake*180/pi

    A = str_dip_rake(strike,dip,rake)

    return A




#%%###################################################################
#   (4) The conversion between A/N vector and P/T/N vector
######################################################################

# Calculate the T-axis, P-axis and N-axis according to the slip vector (A) and fault plane direction vector (N)
def AN2TPN(A,N):
    """
    Input: slip vector(A) and fault plane's normal vector(N) in NED system.

    Output: Tension-axis vector(T), Pressure-axis vector(P) and Null-axis 
            vector(Null) in NED system.
    """
    T=sqrt(2)/2*(A+N)
    P=sqrt(2)/2*(A-N)
    Null=np.cross(P,T)

    return T,P,Null

# Calculate the slip vector (A) and fault plane direction vector (N) according to the T-axis and P-axis
def TP2AN(T,P):
    """
    Input: Tension-axis vector(T) and Pressure-axis vector(P) in NED system.

    Output: slip vector(A) and fault plane's normal vector(N) in NED system.
    """
    A=sqrt(2)/2*(T+P)
    N=sqrt(2)/2*(T-P)

    return A,N




#%%#######################################################
#         (5) mt(in NED system) to P/T/N vector
##########################################################

def MT2TPN(MT_raw):
    """
    Input: moment tensor in NED system.

    Output: Tension-axis vector(T), Pressure-axis vector(P) and Null-axis 
            vector(Null) in NED system.
    """

    M = MT_raw.mt
    eigen_val, eigen_vec = np.linalg.eig(M)
    # The TNP axis should be arranged in order of eigenvalues from largest to smallest
    eigen_vec_ord_axis = np.real( np.take(eigen_vec, np.argsort(-eigen_val), 1) )
    
    T_axis = eigen_vec_ord_axis[:, 0] 
    N_axis = eigen_vec_ord_axis[:, 1] 
    P_axis = eigen_vec_ord_axis[:, 2] 

    return T_axis, P_axis, N_axis




#%%##############################################################
#     (6) P/T/N vector to P/T/N vector's stirke and dip
#################################################################

def vector2str_dip(vector):
    """
    Input: a principal axis vector, such as eigenvectors P/T/N of the moment tensor 
            object in NED system.

    Output: a principal axis' strike and dip.
            strike : [0, 360)
            dip    : [0, 90]
    """
    x=vector[0]
    y=vector[1]
    z=vector[2]

    strike = atan2(y,x)*180/pi 
    r = sqrt(x**2+y**2)
    dip = atan2(z,r)*180/pi

    if dip < 0.0:
        dip = -dip
        strike = strike-180
    if strike < 0:
        strike = strike+360
    if strike > 360:
        strike = strike-360

    A = Axis_str_dip(strike,dip)

    return A





#%%##############################################################
#               (7) project station to beachball
#################################################################
# According to the strike(azimuth) and takeoff Angle (or dip, for the PTN 
# axis, you need to make your own Angle transformation pi/2-TKO=dip) 
# projected onto the beachball

def project_beachball(AZM, TKO, R=1, menthod='schmidt'):
    """
    Input in NED system:
        AZM means azimuth (equal to strike) in degrees.
        TKO means takeoff angle ( pi/2-TKO=dip equal to dip) in degrees.
        R means beachball radius that you want to plot.
        note: Takeoff Angle is the angle with the vertical direction, 
            DIP is the angle with the horizontal plane.

    Output: 
        X and Y coordinates in E and N direction respectively, and
        the lower left corner of the circle is the origin.
    """
    AZM = AZM/180*pi
    TKO = TKO/180*pi

    # Schmidt (Lambert, equal-area) default
    if menthod=='schmidt':
        r = math.sqrt(2)*sin(TKO/2)
    # Wulff projection (Stereographic, equal-angle) not recommmended
    elif menthod=='wulff':
        r = tan(TKO/2)
    else:
        raise ValueError('projection error!')

    X = R*r*sin(AZM)+R
    Y = R*r*cos(AZM)+R

    return X,Y





#%%#########################################
#               (9) Hudson plot
############################################
# Hudson, J.A., R.G. Pearce, and R.M.Rogers (1989), "Source type plot for inversion of the moment tensor",\
# J. Geophys. Res., 94, 765?74

def M2kT_space(MT):
    # 1. full-moment
    M = MT
    # M = np.array([Mxx, Mxy, Mxz,    Mxy, Myy, Myz,    Mxz, Myz, Mzz]).reshape(3, 3) 

    # 2.isotropic part
    m_iso = 1./3 * np.trace(M)
    M_iso = np.diag(np.array( [m_iso,m_iso,m_iso] ))

    # 3.deviatoric part
    M_devi = M - M_iso

    # 4.eigenvalues and -vectors of M
    devi_eigen_val, devi_eigen_vec = np.linalg.eig(M_devi)

    # 5.eigenvalues in ascending order:
    devi_eigen_val_ord = np.real( np.take(devi_eigen_val, np.argsort(-devi_eigen_val)) )   # descend order

    if ( abs(m_iso) + max( abs(devi_eigen_val_ord[0]),abs(devi_eigen_val_ord[2]) ) ) == 0 :
        raise TypeError("MomentTensor cannot be translated into [k,T] space.")
    else:
        k = m_iso / ( abs(m_iso) + max(abs(devi_eigen_val_ord[0]), abs(devi_eigen_val_ord[2])) )

    if max(abs(devi_eigen_val_ord[0]), abs(devi_eigen_val_ord[2])) == 0:
        T = 0
    else:
        T = 2*devi_eigen_val_ord[1] / max(abs(devi_eigen_val_ord[0]), abs(devi_eigen_val_ord[2]))

    return k,T


def kT2UV_space(k,T):
    tau = T*(1-abs(k))
    if ( (tau>0) & (k<0) ) | ( (tau<0) & (k>0) ):
        # 2nd and 4th quadrants
        U = tau
        V = k
    elif ( tau < (4*k) ) & ( (tau>=0) & (k>=0) ):
        # First quadrant, Region A
        U = tau/(1-tau/2)
        V = k/(1-tau/2)
    elif ( tau >= (4*k) ) & ( (tau>=0) & (k>=0) ):
        # First quadrant, Region B
        U = tau/(1-2*k)
        V = k/(1-2*k)
    elif ( tau >= (4*k) ) & ( (tau<=0) & (k<=0) ):
        # Third quadrant, Region A
        U = tau/(1+tau/2)
        V = k/(1+tau/2)
    elif ( tau < (4*k) ) & ( (tau<=0) & (k<=0) ):
        # Third quadrant, Region B
        U = tau/(1+2*k)
        V = k/(1+2*k)
    else:
        raise TypeError("def: kT2UV_space(k,T)")

    return U,V


def Hudson_plot(ax, ms=2, marker_ms='o', color_ms='k', alpha_ms=0.5, alpha_text=0.7, fontsize=6):
    ######################
    ## 1. Fill and draw the border
    ax.fill_between(x=[-1,0],y1=[0,0],y2=[0,1],color='k', alpha=0.05)           # fill the second quadrant
    ax.fill_between(x=[0,1],y1=[0,0],y2=[-1,0],color='k', alpha=0.05)           # fill the fourth quadrant
    ax.plot([0, 4/3, 0, -4/3, 0], [1, 1/3, -1, -1/3, 1],
            linestyle='-', color='k', lw=1, alpha=0.6)
    ax.plot([-1, 1], [0, 0], linestyle='-', color='k', lw=1, alpha=0.6)
    ax.plot([0, 0], [-1, 1], linestyle='-', color='k', lw=1, alpha=0.6)

    ######################
    ## 2. Draw the inner dotted line
    U_vector = [];V_vector = []
    for i in np.linspace(-1, 1, num=100):
        k = i
        T = 0.5
        U,V = kT2UV_space(k=k, T=T)
        U_vector.append(U)
        V_vector.append(V)
    ax.plot(U_vector, V_vector, linestyle='--', color='k', lw=1, alpha=0.6)

    U_vector = [];V_vector = []
    for i in np.linspace(-1, 1, num=100):
        k = i
        T = -0.5
        U,V = kT2UV_space(k=k, T=T)
        U_vector.append(U)
        V_vector.append(V)
    ax.plot(U_vector, V_vector, linestyle='--', color='k', lw=1, alpha=0.6)

    U_vector = [];V_vector = []
    for i in np.linspace(-1, 1, num=100):
        k = 0.5
        T = i
        U,V = kT2UV_space(k=k, T=T)
        U_vector.append(U)
        V_vector.append(V)
    ax.plot(U_vector, V_vector, linestyle='--', color='k', lw=1, alpha=0.6)

    U_vector = [];V_vector = []
    for i in np.linspace(-1, 1, num=100):
        k = -0.5
        T = i
        U,V = kT2UV_space(k=k, T=T)
        U_vector.append(U)
        V_vector.append(V)
    ax.plot(U_vector, V_vector, linestyle='--', color='k', lw=1, alpha=0.6)

    ######################
    ## 3. Draw marker points
    # ms=2
    # marker_ms = 'o'
    # color_ms = 'k'
    # alpha_ms = 0.5
    # alpha_text = 0.7
    # fontsize = 7
    U,V = kT2UV_space(k=1, T=1)
    ax.plot(U,V, marker='o', color=color_ms, ms=ms, alpha=alpha_ms) 
    ax.text(U,V,'ISO+ (Explosion)',horizontalalignment='center', verticalalignment='bottom',\
             fontsize=fontsize, color='k',alpha=alpha_text)

    U,V = kT2UV_space(k=-1, T=1)
    ax.plot(U,V, marker='o', color=color_ms, ms=ms, alpha=alpha_ms)  
    ax.text(U,V,'ISO- (Implosion)',horizontalalignment='center', verticalalignment='top',\
             fontsize=fontsize, color='k',alpha=alpha_text)

    U,V = kT2UV_space(k=0, T=1)
    ax.plot(U,V, marker='o', color=color_ms, ms=ms, alpha=alpha_ms) 
    ax.text(U,V,'CLVD (-)',horizontalalignment='left', verticalalignment='top',\
             fontsize=fontsize, color='k',alpha=alpha_text)

    U,V = kT2UV_space(k=-5/9, T=1)
    ax.plot(U,V, marker='o', color=color_ms, ms=ms, alpha=alpha_ms) 
    ax.text(U,V,'Anticrack',horizontalalignment='left', verticalalignment='top',\
              fontsize=fontsize, color='k',alpha=alpha_text)

    U,V = kT2UV_space(k=0, T=-1)
    ax.plot(U,V, marker='o', color=color_ms, ms=ms, alpha=alpha_ms) 
    ax.text(U,V,'CLVD (+)',horizontalalignment='right', verticalalignment='bottom',\
             fontsize=fontsize, color='k',alpha=alpha_text)

    U,V = kT2UV_space(k=5/9, T=-1)
    ax.plot(U,V, marker='o', color=color_ms, ms=ms, alpha=alpha_ms) 
    ax.text(U,V,'Tensile Crack',horizontalalignment='right', verticalalignment='bottom',\
             fontsize=fontsize, color='k',alpha=alpha_text)

    U,V = kT2UV_space(k=0, T=0)
    ax.plot(U,V, marker='o', color=color_ms, ms=ms, alpha=alpha_ms) 
    ax.text(U,V,'DC',horizontalalignment='center', verticalalignment='bottom',\
             fontsize=fontsize, color='k',alpha=alpha_text)

    U,V = kT2UV_space(k=1/3, T=-1)
    ax.plot(U,V, marker='o', color=color_ms, ms=ms, alpha=alpha_ms) 
    ax.text(U,V,'LVD (+)',horizontalalignment='right', verticalalignment='bottom',\
             fontsize=fontsize, color='k',alpha=alpha_text)

    U,V = kT2UV_space(k=-1/3, T=1)
    ax.plot(U,V, marker='o', color=color_ms, ms=ms, alpha=alpha_ms) 
    ax.text(U,V,'LVD (-)',horizontalalignment='left', verticalalignment='top',\
             fontsize=fontsize, color='k',alpha=alpha_text)

    ######################
    ## 4. Set the axes
    ax.set_xlim(-4/3-0.1, 4/3+0.1)
    ax.set_ylim(-1-0.1, 1+0.1)
    ax.set_aspect("equal")
    ax.set_axis_off() 







#%%########################################################### 
# (10) Describe_fault_plane with two str_dip_rake.
############################################################## 

def describe_fault_plane(fm):
    """
    Input: moment tensor object in NED system.

    Output: [strike_1, dip_1, rake_1] and [strike_2, dip_2, rake_2]
    """
    MT = MTensor(fm)

    T_axis, P_axis, N_axis = MT2TPN(MT)
    A,N = TP2AN(T_axis,P_axis)

    a = AN2str_dip_rake(A,N)
    strike_1 = a.strike
    dip_1 = a.dip
    rake_1 = a.rake

    b = AN2str_dip_rake(N,A)
    strike_2 = b.strike
    dip_2 = b.dip
    rake_2 = b.rake

    return np.array([[strike_1,dip_1,rake_1],[strike_2, dip_2, rake_2]])






#%%################################
#          (11) object.
###################################

class MTensor(object):
    """
    Adapted from obspy.A moment tensor in NED system.

    >>> a = MTensor([1, 1, 0, 0, 0, -1])                                          # MTensor(Mxx, Myy, Mzz, Mxy, Mxz, Myz)
    >>> b = MTensor(np.array([[1, 0, 0], [0, 1, -1], [0, -1, 0]]))
    >>> c = MTensor([100,50,30])
    >>> a.mt
    array([[ 1,  0,  0],
           [ 0,  1, -1],
           [ 0, -1,  0]])
    >>> b.yz
    -1
    """
    def __init__(self, a):
        if len(a) == 3 and isinstance(a, list):
            # strike dip rake
            MT = str_dip_rake2MT(a[0],a[1],a[2])
            self.mt = MT.mt
        elif len(a) == 6:
            # six independent components
            self.mt = np.array([[a[0], a[3], a[4]],
                                    [a[3], a[1], a[5]],
                                    [a[4], a[5], a[2]]])
        elif isinstance(a, np.ndarray) and a.shape == (3, 3):
            # full matrix
            self.mt = a
        else:
            raise TypeError("Wrong size of input parameter.")

    @property
    def mt_normalized(self):
        return self.mt / np.linalg.norm(self.mt)

    @property
    def xx(self):
        return self.mt[0][0]

    @property
    def xy(self):
        return self.mt[0][1]

    @property
    def xz(self):
        return self.mt[0][2]

    @property
    def yz(self):
        return self.mt[1][2]

    @property
    def yy(self):
        return self.mt[1][1]

    @property
    def zz(self):
        return self.mt[2][2]


class str_dip_rake(object):
    """
    Describing the faultplanes of the Double Couple

    Strike dip and rake values are in degrees.
        strike : [0, 360)
        dip    : [0, 90]
        rake   : [-180, 180)

    >>> a = str_dip_rake(20, 50, 10)
    >>> a.strike
    20
    >>> a.dip
    50
    >>> a.rake
    10
    """
    def __init__(self, strike=0, dip=0,rake=0):
        self.strike = strike
        self.dip = dip
        self.rake = rake


class Axis_str_dip(object):
    """
    A principal axis' strike and dip.
    Used in P/T/N's axis

    Strike and dip values are in degrees.
        strike : [0, 360)
        dip    : [0, 90]

    >>> a = Axis_str_dip(20, 50)
    >>> a.strike
    20
    >>> a.dip
    50
    """
    def __init__(self, strike=0, dip=0):
        self.strike = strike
        self.dip = dip





#%%###########################################################################
#    (9) Decompose of mt: isotropic + deviatoric = isotropic + DC + CLVD
##############################################################################

class Decompose(object):
    """
    Creates a Decompose object on the basis of a provided MomentTensor object.

    For example:
        m = str_dip_rake2MT(120,50,70)
        AA = Decompose(m)
        AA.decomposition_iso_DC_CLVD()
        AA.M_DC_percentage -> 
    """
    def __init__(self, MT_raw):

        self.M = MT_raw.mt
        self.M_iso = None
        self.M_devi = None
        self.M_DC = None
        self.M_CLVD = None
        self.M0 = None
        self.Mw = None
        self.M_iso_percentage = None
        self.M_DC_percentage = None
        self.M_CLVD_percentage = None
        self.eigen_val = None
        self.eigen_vec = None
        self.F = None

    def decomposition_iso_DC_CLVD(self):
        """
        Input: moment tensor in NED system.

        Output: Tension-axis vector(T), Pressure-axis vector(P) and Null-axis 
                vector(Null) in NED system.

        Decomposition according Aki & Richards and Jost & Herrmann into

        isotropic + deviatoric
        = isotropic + DC + CLVD

        parts of the input moment tensor.

        results are given as attributes, callable via the get_* function:

        DC, CLVD, DC_percentage, seismic_moment, moment_magnitude
        """

        # 1. full-moment
        M = self.M

        # 2.isotropic part
        m_iso = 1./3 * np.trace(M)
        M_iso = np.diag(np.array( [m_iso,m_iso,m_iso] ))
        m0_iso = abs(m_iso)

        # 3.deviatoric part
        M_devi = M - M_iso

        # 4.eigenvalues and -vectors of M
        eigen_val, eigen_vec = np.linalg.eig(M)

        # 5.eigenvalues in ascending order:
        eigen_val_ord = np.real( np.take(eigen_val, np.argsort(abs(eigen_val))) )
        eigen_vec_ord = np.real( np.take(eigen_vec, np.argsort(abs(eigen_val)), 1) )

        # 6.named according to Jost & Herrmann:
        # a1 = eigen_vec_ord[:, 0]
        a2 = eigen_vec_ord[:, 1]
        a3 = eigen_vec_ord[:, 2]
        F = -(eigen_val_ord[0]-m_iso) / (eigen_val_ord[2]-m_iso)

        # 7.decompose
        M_DC = (eigen_val_ord[2]-m_iso) * (1 - 2 * F) * (np.outer(a3, a3) - np.outer(a2, a2))
        M_CLVD = M_devi - M_DC

        # 8.according to Bowers & Hudson:
        M0 = max(abs(eigen_val_ord))                                            # Seismic moment (in Nm)
        Mw = np.log10(M0 * 1.0e7) / 1.5 - 16.1/1.5                              # moment_magnitude unit is Mw

        M_iso_percentage = int(round(m0_iso / M0 * 100, 6))

        M_DC_percentage = int(round((1 - 2 * abs(F)) *
                                    (1 - M_iso_percentage / 100.) * 100, 6))

        M_CLVD_percentage = 100-M_iso_percentage-M_DC_percentage

        self.M_iso = M_iso
        self.M_devi = M_devi
        self.M_DC = M_DC
        self.M_CLVD = M_CLVD
        self.M0 = M0
        self.Mw = Mw
        self.M_iso_percentage = M_iso_percentage
        self.M_DC_percentage = M_DC_percentage
        self.M_CLVD_percentage = M_CLVD_percentage
        self.eigen_val = eigen_val
        self.eigen_vec = eigen_vec
        self.F = F


    def help(self):
        print("Incluing function:\n\
self.M\n\
self.M_iso\n\
self.M_devi\n\
self.M_DC\n\
self.M_CLVD\n\
self.M0\n\
self.Mw\n\
self.M_iso_percentage\n\
self.M_DC_percentage\n\
self.M_CLVD_percentage\n\
self.eigen_val\n\
self.eigen_vec\n\
self.F")


    def print_self(self):
        print("self.M:",self.M,"\n")
        print("self.M_iso:",self.M_iso,"\n")
        print("self.M_devi:",self.M_devi,"\n")
        print("self.M_DC:",self.M_DC,"\n")
        print("self.M_CLVD:",self.M_CLVD,"\n")
        print("self.M0:",self.M0,"\n")
        print("self.Mw:",self.Mw,"\n")
        print("self.M_iso_percentage:",self.M_iso_percentage,"\n")
        print("self.M_DC_percentage:",self.M_DC_percentage,"\n")
        print("self.M_CLVD_percentage:",self.M_CLVD_percentage,"\n")
        print("self.eigen_val:",self.eigen_val,"\n")
        print("self.eigen_vec:",self.eigen_vec,"\n")
        print("self.F:",self.F,"\n")





