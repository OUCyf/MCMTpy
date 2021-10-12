Conversion of Source Parameters
================================

* The MCMTpy provides a series of scripts for source parameter calculation and conversion. The MomentTensor script refers 
  to some theories and codes of `Obspy  <https://github.com/obspy/obspy>`_ `MoPaD  <https://github.com/geophysics/MoPaD>`_ and 
  Introduction to Seismology (Yongge Wan), mainly including:

.. contents::
    :local:
    :depth: 1

* We do not recommend using relative paths because of the possibility of errors. Absolute paths are preferred.






The installation
------------------------------------
* jupyter notebook::

    # import MomentTensor
    from MCMTpy import MomentTensor as MTpy

    # get __doc__
    MTpy.__doc__.strip().split("\n")




str_dip_rake to mt
------------------------------------
* jupyter notebook::
 
    # function: str_dip_rake2MT
    MTpy.str_dip_rake2MT.__doc__.strip().split("\n")

    # input
    strike = 50
    dip = 50
    rake = 100

    A = MTpy.str_dip_rake2MT(strike,dip,rake)
    A.mt




The conversion between str_dip_rake and A/N vector
-----------------------------------------------------
* jupyter notebook::

    #-----------------------------------------------#
    # function: str_dip_rake2AN
    MTpy.str_dip_rake2AN.__doc__.strip().split("\n")

    # str_dip_rake to A/N vector
    # input
    strike = 50
    dip = 50
    rake = 100
    
    A,N = MTpy.str_dip_rake2AN(strike,dip,rake)
    print('A = ', A)
    print('N = ', N)

    #-----------------------------------------------#
    # function AN2str_dip_rake
    MTpy.AN2str_dip_rake.__doc__.strip().split("\n")

    # A/N vector to str_dip_rake
    # input
    A = np.array([0.37330426, -0.53992106, -0.75440651])
    N = np.array([-0.58682409,  0.49240388, -0.64278761])
    
    DD = MTpy.AN2str_dip_rake(A,N)
    print('strike = ', DD.strike)
    print('dip = ', DD.dip)
    print('rake = ', DD.rake)



The conversion between A/N vector and P/T/N vector
-----------------------------------------------------
* jupyter notebook::

    #-----------------------------------------------#
    # function AN2TPN
    MTpy.AN2TPN.__doc__.strip().split("\n")

    # Calculate the T-axis, P-axis and N-axis according to the slip vector (A) and fault plane direction vector (N)
    # input
    A = np.array([0.37330426, -0.53992106, -0.75440651])
    N = np.array([-0.58682409,  0.49240388, -0.64278761])
    
    T,P,Null = MTpy.AN2TPN(A,N)
    print('T = ', T)
    print('P = ', P)
    print('Null = ', Null)

    #-----------------------------------------------#
    # function TP2AN
    MTpy.TP2AN.__doc__.strip().split("\n")

    # Calculate the slip vector (A) and fault plane direction vector (N) according to the T-axis and P-axis
    # input
    T = np.array([-0.15098132, -0.03359972, -0.98796544])
    P = np.array([0.67891327, -0.72996397, -0.07892648])
    
    A,N = MTpy.TP2AN(T,P)
    print('A = ', A)
    print('N = ', N)




mt to P/T/N vector
------------------------------------
* jupyter notebook::

    # function MT2TPN
    MTpy.MT2TPN.__doc__.strip().split("\n")

    # input
    strike = 50
    dip = 50
    rake = 100
    
    A = MTpy.str_dip_rake2MT(strike,dip,rake)
    T, P, Null = MTpy.MT2TPN(A)
    print('T = ', T)
    print('P = ', P)
    print('Null = ', Null)







P/T/N vector to P/T/N vector's stirke and dip
-----------------------------------------------------
* jupyter notebook::

    # function vector2str_dip
    MTpy.vector2str_dip.__doc__.strip().split("\n")

    # input
    A = np.array([0.37330426, -0.53992106, -0.75440651])

    CC = MTpy.vector2str_dip(A)
    print('A.strike = ', CC.strike)
    print('A.dip = ', CC.dip)






Describe_fault_plane with two str_dip_rake
-----------------------------------------------------
* jupyter notebook::

    # function describe_fault_plane
    MTpy.describe_fault_plane.__doc__.strip().split("\n")

    # input
    strike = 50
    dip = 50
    rake = 100
    
    A = MTpy.str_dip_rake2MT(strike,dip,rake)
    CC = MTpy.describe_fault_plane(A.mt)
    print('FM_1 = ', CC[0,:])
    print('FM_2 = ', CC[1,:])

