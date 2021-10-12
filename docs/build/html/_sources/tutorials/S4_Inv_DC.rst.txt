Inversion of  Source
==================================

.. contents::
    :local:
    :depth: 1

* This is the core of the MCMTpy, and all parameters are written into a JSON file named **sample_dc.json**. 
  Some brief descriptions of the parameters are included here following the definination.

* We do not recommend using relative paths because of the possibility of errors. **Absolute paths** are preferred.




Path
------------------------
``GFs_json_file``

* The path to the Json-file used to calculate Green Function Database.
  
``GFs_file``

* The path to the Green Function Database.

``Raw_data_file``

* The path to the Raw data ASDF file.

``Source_tag``

* Raw data ASDF file's source name which setted in the data preprocess function.

``Raw_p_n0``

* The number of points before the first arrival. Be consistent with the GFs' ``samples_before_first_arrival``.

``Dt``

* Sampling rate of the raw data.

``Output_path``

* Inversion result output file path.






Source Time Function
------------------------
``Trapezoid_stf_mode``   

* Whether to use pyfk provided trapezoid shaped source time function. Set to be **.true.** you need to set the following parameters correctly: 

  ``dura``
  ``rise``
  ``delta``  

``dura``    

* Duration time (s).

``rise``  

* Rise time (s).

``delta``  

* The time interval (s).

``Stf_file``    

* Stf_file should be an **SAC** file as the source time function when set ``Trapezoid_stf_mode`` to be **.false.**.
  
.. note::
    To read the data from the ``Stf_file`` file, keep the sampling rate of SAC data consistent with the data and GFs.





MPI and Chains_n
------------------------
``MPI_n``

* CPU number

``Chains_n``

* A total of n Markov Chains are sampled on each core (CPU).

``N``

* Each Markov Chains' sampling number.






Misfit
------------------------
``Geometrical_spreading_mode``

* Whether the geometric diffusion correction is turned on or not. Waveform normalization is used when 
  it set to be **.false.**, and please set the magnitude **M0** not to be inverted in ``Fixed_FM``. 
  Only one of ``Geometrical_spreading_mode`` and ``Amplitude_normalization_mode`` can be **.true.**.

``Distance_0``

* The reference epicentral distance when ``Geometrical_spreading_mode`` set to be **.true.**.

``Amplitude_normalization_mode``

* Waveform normalization mode.Only one of ``Geometrical_spreading_mode`` and 
  ``Amplitude_normalization_mode`` can be **.true.**.

``Use_file_mode``

* When set it to be **.false.**, don't need to provide ``Raw_data_inv_info`` file. Setting following parameters 
  (``Raw_data_inv_info_writed``,
  ``NET_STA``,
  ``P_inv_comp``,
  ``S_inv_comp``,
  ``Surf_inv_comp``,
  ``P_win``,
  ``S_win``,
  ``Surf_win``,
  ``Phase_weight``,
  ``Distance_weight``,
  ``P_maxlag``,
  ``S_maxlag``,
  ``Surf_maxlag``,
  ``P_filter``,
  ``S_filter``,
  ``Surf_filter``)
  correctly, ``MCMTpy`` will generate this file automatically, so you can change individual parameter values 
  and use this file for inversion. ``Raw_data_inv_info_writed`` is the filename that generate this file automatically.
  When set it to be **.true.**, ``Raw_data_inv_info`` is needed which is the path to inversion file.

``Raw_data_inv_info``

* The path to inversion parameters file. It includes the time window length, filtering range, weight 
  and other parameters used for data. When set ``Use_file_mode`` to be **.true.**, it must be provided.

* Raw_data_inv_info.txt format.

.. literalinclude:: S4_Raw_data_inv_info.txt
  :language: none
  :linenos:

``Raw_data_inv_info_writed``

* ``Raw_data_inv_info_writed`` is the filename that generate this file automatically When set ``Use_file_mode`` 
  to be **.true.**.



``NET_STA``

* The information contained in the two-dimensional list: [network_name]	[station_name]	[station_lat]	[source_lon]	[station_depth]

``P_inv_comp``

* P wave Z/R/T component. 'yes': inversion. 'no': no inversion. Select the component you want to invert.

``S_inv_comp``

* S wave Z/R/T component. 'yes': inversion. 'no': no inversion. Select the component you want to invert.

``Surf_inv_comp``


* Surf wave Z/R/T component. 'yes': inversion. 'no': no inversion. Select the component you want to invert.

``P_win``

* P wave P_win=[-1,2] unit/s; The P wave fitting band, and '-1' means 1 second before the P wave, '2' means 2 second after the P wave

``S_win``

* S wave S_win...

``Surf_win``


* Surf wave Surf_win...

``Phase_weight``

* P/S/Surface wave inversion weight

``Distance_weight``


* The weight of each station in ``NET_STA``

``P_maxlag``

* The maximum number of moveable sampling points to fitting waveform of P wave. 
  It is generally half the period of the P phase waveform. unit/s. 

``S_maxlag``

* The maximum number of moveable sampling points to fitting waveform of S wave. 

``Surf_maxlag``


* The maximum number of moveable sampling points to fitting waveform of Surf wave. 

``P_filter``

* P wave filter range, [freqmin, freqmax].

``S_filter``

* S wave filter range, [freqmin, freqmax].

``Surf_filter``

* Surf wave filter range, [freqmin, freqmax].






Inversion parameters
------------------------
``InvType``

* The type of source, mt(moment tensor) dc(double-couple) sf(single-couple) ep(explosion).

``FM0``

* The initial value of the inversion, it's a 1-d list.
* dc:[m0, str, dip rake, lat, lon, depth, t0] || mt:[m0, m_xx, m_xy, m_xz, m_yy, m_yz, m_zz, lat, lon, depth, t0] || ep:[m0]
* The default initial value of t0 is usually 0.

``Fixed_FM``

* The length of the list corresponds to ``FM0``, and 'constant' means fixed to a constant and unchanged in the inversion process. 
  'variable' is represented as a variable and participates in the inversion.

``FM_boundary``



* A 2-d list, FM's boundary.

``Noise_level_waveform``

* Noise level estimation of raw data.

``MISFIT0``



* The initial error of the objective function, usually set to a very large value.

``alpha_max``

* The initial value of **alpha**.

``alpha_min``

* The final value of **alpha**.

``a``

* The alpha function's parameter.

``b``



* The alpha function's parameter.

``sigma_mw``

* The standard deviation of the new solution of **Mw**.

``sigma_mt``

* The standard deviation of the new solution of **Moment Tensor**.

``sigma_dc``

* The standard deviation of the new solution of **Strike Dip Rake**.

``sigma_sf``

* The standard deviation of the new solution of **Single Couple** (not test now).

``sigma_ep``

* The standard deviation of the new solution of **Explosion**  (not test now).

``sigma_lat_lon``

* The standard deviation of the new solution of **Latitude and Longitude** (degree).

``sigma_depth``

* The standard deviation of the new solution of **Source Depth** (km).

``sigma_t``

* The standard deviation of the new solution of **T0** (s).



Alpha
------------------------ 
* The **alpha** function is shown in following script, you can test different parameters 
  ``alpha_max``
  ``alpha_min``
  ``a``
  ``b``
  and plot thefigure.

.. literalinclude:: S4_alpha.py
  :language: python
  :linenos:



Logging file
------------------------ 
* **MCMTpy** will generates 1 logging files when it do **Inversion**, that is **Inv_para_info.txt**. 
  They record all the parameter information of inversion process. You can find them 
  in ``Output_path`` path.




Example Double-Couple Inv
----------------------------  
* The **example_path** need to de changed to your path, and run this notebook to set **sample_dc.json** file.

.. literalinclude:: S4_inv.py
  :language: python
  :linenos:

* Now you can run MCMTpy in bash::
  
  $ MCMTpy  sample MH  -c ./sample_dc_new.json > sample.log



.. note::
    Please follow the above parameter instructions and set all parameters correctly before running the program.
    Otherwise, it is easy to report an error!
