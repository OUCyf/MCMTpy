Synthesize the Test Data
========================

.. contents::
    :local:
    :depth: 1


* The synthesized test data can help us check whether the Green's function database is correctly calculated.

* We do not recommend using relative paths because of the possibility of errors. **Absolute paths** are preferred.



Path
------------------------
``GFs_json_file``

* The path to the Json-file used to calculate Green Function Database.


Source
------------------------
``srcType``

* The type of *Source*, mt(moment tensor) dc(double-couple) sf(single-couple) ep(explosion).

``point`` 

* The latitude longitude and depth of the source.

``source_mechanism`` 

* Focal mechanism. 
* When ``srcType`` = **dc**: ``source_mechanism`` = [Mw, Strike, Dip, Rake]
* When ``srcType`` = **mt**: ``source_mechanism`` = [Mw, Mxx, Mxy, Mxz, Myy, Myz, Mzz]






Station
------------------------  
``NET_STA``    

* The information contained in the two-dimensional list: [network_name]	[station_name]	[station_lat]	[source_lon]	[station_depth]




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



Filter
------------------------  
``filter_mode``  

* Whether to filter.

``freqmin`` 

* The minimum frequency of filtering.

``freqmax``    

* The maximum frequency of filtering.







Add Noise
------------------------
    
``Add_noise``    

* Whether add noise to syn.

``noise_level_waveform`` 

* Noise level added to the waveform.

``noise_level_arrive_time`` 

* Noise level added to the phase time as **tp** and **ts**.






Output
------------------------ 
``Output_mode`` 

* Output data or not.

``Output_path``

* Path of output data

``source_name`` 

* Event name of output data, you can set it any way you want. Only lowercase letters are supported.

``ASDF_filename``    

* The name of the output ASDF file




Example
------------------------ 
* The **example_path** need to de changed to your path, and run this notebook to set **syn.json** file.

.. literalinclude:: S2_syn_waveform.py
  :language: python
  :linenos:

* Now you can run MCMTpy in bash::
  
  $ MCMTpy  syn pyfk  -c ./syn.json



* View syn waveform:

.. literalinclude:: S2_plot_syn.py
  :language: python
  :linenos:

.. note::
    Please follow the above parameter instructions and set all parameters correctly before running the program.
    Otherwise, it is easy to report an error!
