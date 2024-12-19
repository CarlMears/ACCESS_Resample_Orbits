# AMSR2 Access Processing
This project performs the AMSR2 resampling onto a fixed grid (either a rectangular lat/lon grid or a 25 km ease2 north pole polar grid)

The starting point is AMSR2 L1A files.
The ending point is gridded, resampled TB files for each channel

Currently, the resampling is hard coded to be 30km - needs to change

## Special execution considerations
This can not be run in a conda prompt, because it finds the wrong netcdf libraries.

## Entry Point
There are is one way to run the code.

executing resample_amsr2_exe_v2.exe, e.g.

~~~
resample_amsr2_exe_v2.exe 829 831 70 0 1 1 2 3 4 5 6 7 8 9 10 11 12

resample_amsr2_exe_v2.exe {start_orbit} {end_orbit} {footprint_size} {do_global} {do_polar} {list_of_channels} 
~~~
arguments are:
* start_orbit:  start_orbit_number
* end_orbit: end_orbit_number
* footprint_size: target footprint diameter in km
* do_global: set to one if global map output is desired, otherwise 0
* do_polar: set to one if north pole output is desired, otherwise 0
* list of channels: list of channels to process (1..12)

code for this entry point is in 
AMSR2\process_amsr2_l2b_orbit_tb_ACCESS.f90

## python driver

The file (and others)
~~~
build\start_orbit_resample_smops_30km_north.py
~~~
contains a python driver which call the .exe file through the OS.  It counts processes with similar names to avoid overloading the processors. I Recommend setting the MAX_PROCESSES parameter to <= 80% of the number of processors available.  Each launch event is delayed by 5s to avoid collisions during the process launch.

Here is the core part of the python driver:
~~~
while len(orbit_list) > 0:
        orbit_to_do = orbit_list.pop()
        while count_processes_with_name(exe_name) >= MAX_PROCESSES:
            if VERBOSE:
                print(f'To many processes: {count_processes_with_name(exe_name)}')
                print(f'Sleeping {SECONDS_TO_SLEEP} seconds')
            time.sleep(SECONDS_TO_SLEEP)
        print(f'{count_processes_with_name(exe_name)} processed running')
        print(f'Starting Process for Orbit: {orbit_to_do}')
        p = start_orbit_resample(sat_name,
                        target_size,
                        orbit_to_do,
                        do_global=False)
        time.sleep(5)
~~~



## building the project

The project is built using meson.
The build is defined by the information in ./meson.build

To build the project, cd to the ./build directory
ninja

Most of the time, any changes in the meson.build file will be automatically implemented

It creates two executables:
- resample_amsr2_exe_v2.exe -> does the resampling
- test_quad_exe.exe -> tests the quadrilateral interpolation


## Wish List




## Calling from python using interface code -- not used and not maintained
## including here in case we want to rejuvenate it

OR

calling directly from python, e.g.

~~~
from interfaces import AMSR2_Lib

def main():
    for iorbit in range(5218,5501):
        channel_list = [i for i in range(5,13)]
        write_swath_tbs = False
        write_swath_geoloc = False
        overwrite_existing_files = False
        #pathname_tb_output = "A/path/set/by/python/A/path/set/by/python/A/path/set/by/python/A/path/set/by/python"
        pathname_tb_output = 'L:/access/amsr2_tb_orbits/'
        AMSR2_lib_interface = AMSR2_Lib()
        AMSR2_lib_interface.process_and_resample_amsr2_l2b_orbit(
            iorbit, channel_list, write_swath_tbs, write_swath_geoloc, overwrite_existing_files, pathname_tb_output
            )
  
if __name__ == "__main__":
    main()
~~~

Python drivers (including ones that use process pools) are in the /python folder

Interface code is in 

AMSR2\amsr2_access_lib.f90

