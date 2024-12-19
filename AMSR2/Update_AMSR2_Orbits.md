# Update Instruction for AMSR2 Tbs Orbits for ACCESS

This is the first step to performing an update to the AMSR2 Access Data.  Here "Update" means extending the dataset forward in time.

The process here makes the TB orbit files with resampled brightness temperatures geolocated onto various grids.

There are two python driver routines:

~~~
AMSR2_start_orbit_resample_30km_all_regions.py
AMSR2_start_orbit_resample_70km_all_regions.py
~~~

The orbit number in these scripts should be updated, and then run.  

A good place to run them is on ops1p-ifort20-rt-1 on OPS1D, which has a lot of processors (40) and memory (64 GB).  Running with 24 processors is successful.  It might be able to be increased further since it only uses 50% of the memory.

The orbit files are written to the directory structure in

~~~
L:\access\amsr2_tb_orbits\
~~~

They need to be moved to the proper locations:

e.g.
~~~
L:\access\amsr2_tb_orbits\GL_30
L:\access\amsr2_tb_orbits\GL_70
L:\access\amsr2_tb_orbits\NP_30
...etc....
~~~

The scripts to do this are in files like:
~~~
M:\job_access\python\dataset_assembly\util\move_tb_orbit_files_GL_30.py
M:\job_access\python\dataset_assembly\util\move_tb_orbit_files_GL_70.py
M:\job_access\python\dataset_assembly\util\move_tb_orbit_files_NP_30.py
~~~






