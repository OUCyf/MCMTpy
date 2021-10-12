Prepare Data for MCMTpy
========================

.. contents::
    :local:
    :depth: 1


Format of MCMTpy input data:
-------------------------------

* **ASDF** file format, and the sampling rate ``dt`` is consistent with that of GFs database.

* Must be **ENZ** three-component data and the name of the channel must be 'E'/'N'/'Z' order.

.. note::
    Pyasdf files are automatically sorted by name. So you must make sure that the order of the Trace 
    is still ENZ after the automatic sorting, generally the name of the channel  'BHE'/'BHN'/'BHZ' will be fine.

* **Trace.Stats.sac.t1**/**Trace.Stats.sac.t2** must be included in the trace header file to represent the P/S phase time.

* The trace header file must include the **station** and **network** name.

* The trace header file must include **b/e/o/stla/stlo/stdp/evlo/evla/evdp/mag/dist/az/baz**, time information, 
  event longitude and latitude information, azimuth information.

* The trace header file must include **starttime/endtime**.

* All traces need to retain a uniform number of sampling points before T1, such as **p_n0 = 50**.





Example
------------------------ 
* We have provided some useful scripts for data preprocessing, and the **example_path** need to de changed to your path, 
  and run this notebook.

.. literalinclude:: S3_processing_data.py
  :language: python
  :linenos:


.. note::
    Please follow the above parameter instructions and set all parameters correctly before running the program.
    Otherwise, it is easy to report an error!