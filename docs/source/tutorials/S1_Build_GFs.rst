Calculate Green Function Database
================================================


.. contents::
    :local:
    :depth: 1

* The first step of **MCMTpy** is to build Green's function database. All parameters are written into a JSON 
  file named **build_GFs.json**. Some brief descriptions of the parameters are included here following the 
  definination. 

* We do not recommend using relative paths because of the possibility of errors. **Absolute paths** are preferred.

* We used ``ASDF`` seismic data format to manage Green Function database, which is convenient to locally build up 
  a database of pre-processed waveforms. ASDF file greatly reduces the number of green functions files, 
  and a single ASDF file can replace thousands of green function files which may threw us into confusion 
  in the management process. We store the Green‘s functions corresponding to **Event_ID** (**QuakeML**) 
  in **“Waveforms”** of ASDF data. In the **“Labels”** of **“Waveforms”** , we store the phase (P and S) travel times 
  which can be computed in advance by a given velocity model.

* The general structure of MCMTpy stores Green‘s function with ASDF file format shown in the image below. 

.. image:: ../../figures/Green_Function_Database.png
    :width: 100%
    :align: center

* You can also find some useful information about parameters in refers to `pyfk <https://github.com/ziyixi/pyfk>`_.





Station and Source
------------------------
``Database_mode``

* Set to be **.true.** when you want to compute the Green's function library for some region. For real 
  examples, the GF databases may require up to several GBs of free disc space. 
  In addition you need to set the following parameters correctly: 

  ``Source_depth_min``,
  ``Source_depth_max``,
  ``Source_depth_spacing``,
  ``Station_distance_radius_min``,
  ``Station_distance_radius_max``,
  ``Station_distance_spacing``,
  ``Source_name``,
  ``Network_name``,
  ``Station_name``,
  ``Station_depth_reference``

* Set to be **.false.** when you want to calculate the Green's function for stations at which the source is given. 
  You must provide an ``.txt`` file named Source_Station_info.txt, which including the longitude latitude and depth 
  information of the source and station. Different sources should be separated by blank lines.

* Source_Station_info.txt format.
  
  [source_name]	[source_lat]	[source_lon]	[source_depth]

  [network_name]	[station_name]	[station_lat]	[source_lon]	[station_depth]

.. literalinclude:: S1_Source_Station_info.txt
  :language: none
  :linenos:


``Source_depth_min``

* The minimum source depth in Green's function databases (km).

``Source_depth_max``

* The maximum source depth in Green's function databases (km).

``Source_depth_spacing``

* Interval of source depths in Green's function databases (km).

``Station_distance_radius_min``

* The minimum distance of the station in Green's function databases (km or degree, depend by ``degrees``).

``Station_distance_radius_max``

* The maximum distance of the station in Green's function databases (km or degree, depend by ``degrees``).

``Station_distance_spacing``

* Horizontal spacing interval in Green's function databases (km or degree, depend by ``degrees``).

``Source_name``

* Source name in Green's function databases, default value is *source*, please do not change.
  Only lowercase letters are supported.

``Network_name``

* Network name in Green's function databases, default value is *NET*, please do not change.

``Station_name``

* Station name in Green's function databases, default value is *STA*, please do not change.

``Station_depth_reference``

* Depth of all virtual stations (km).


.. note::
  When ``Database_mode`` = **.false.**, depths can vary from station to station. But when set to be **.true**, 
  all stations must use the same depth ``Station_depth_reference``.






MPI
------------------------

``MPI_n``

* CPU number. When ``Database_mode`` = **.false.**, the number of call cores should be less than or equal to the number of source.
  But when set to be **.true**, there are no limits.





Path
------------------------

``DATADIR``

* The output path of the Green function databases.

``DATADIR_splits``

* A subfolder under ``DATADIR`` for Green's function. Please keep the path relative to ``DATADIR``.





SourceModel
------------------------

``source_mechanism``

* Focal mechanism, not used. Please set it to **.null.**.

``srcType``

* The type of Green's function, dc(double-couple) sf(single-couple) ep(explosion). Please set it to **.dc.**.





SeisModel
------------------------

``Velocity_model``

* Path to the velocity model file. And the format:
  
  thick(km), vs(km/s), vp(km/s), density(g/cm^3), Qs, Qp

.. literalinclude:: S1_v_model_yunnan.txt
  :language: none
  :linenos:

``flattening``

* Set the flatten status of the velocity model, **.true.** or **.false.**.



Config
------------------------
The following parameters explanation refers to ``pyfk`` <https://github.com/ziyixi/pyfk>.


``npt``

* The sampling points is a multiple of 2.

``dt``

* Sampling interval in seconds

``degrees``

* Use degrees instead of km, defaults to **.true.**. Only when ``Database_mode`` = **.true.**, it can be changed to **.false.**

``taper``

* Taper applies a low-pass cosine filter at fc=(1-taper)*f_Niquest, defaults to 0.3.

``filter``

* Apply a high-pass filter with a cosine transition zone between freq. f1 and f2 in Hz, defaults to (0, 0).

``dk``

* The non-dimensional sampling interval of wavenumber, defaults to 0.3.

``smth``

* Makes the final sampling interval to be dt/smth, defaults to 1.

``pmin``

* The min slownesses in term of 1/vs_at_the_source, defaults to 0.

``pmax``

* The max slownesses in term of 1/vs_at_the_source, defaults to 1.

``kmax``

* The kmax at zero frequency in term of 1/hs, defaults to 15.
  
``rdep``

* Not used now. Please keep the default values of the parameters, defaults to 0 in km.

``updn``

* The "up" for up-going wave only, "down" for down-going wave only, "all" for both "up" and "down", defaults to "all".

``samples_before_first_arrival``

* The number of points before the first arrival, defaults to 50.




Logging file
------------------------ 
* **MCMTpy** will generates 2 logging files when it builds GFs datebase, that is **prepro_para_info.txt** and 
  **write_Source_Station_info_MPI**. They record all the parameter information of GFs datebase. You can find them 
  in ``DATADIR`` path.




Example
------------------------ 
* The **example_path** need to de changed to your path, and run this notebook to set **build_GFs.json** file.

.. literalinclude:: S1_setup_gfs.py
  :language: python
  :linenos:

* Now you can run MCMTpy in bash::
  
  $ mpirun -n 4 MCMTpy build_GFs pyfk -c build_GFs_new.json  > gfs.log



.. note::
    Please follow the above parameter instructions and set all parameters correctly before running the program.
    Otherwise, it is easy to report an error!
