import os
import subprocess
from pathlib import Path
import time

ACCESS_ROOT = Path("L:/access")

path_to_AMSR2_exe = Path('M:/job_access/resampling/fortran_python/build/resample_amsr2_exe_v2_fast_v2.exe')
AMSR2_exe_name = path_to_AMSR2_exe.name

path_to_SSMI_exe = Path('M:/job_access/resampling/fortran_python/ssmi/build/resample_ssmi_exe.exe')
SSMI_exe_name = path_to_SSMI_exe.name

def get_resampled_file_name(satellite: str,
                            channel,
                            target_size: int,
                            orbit: int,
                            ksat = 13,
                            grid_type: str = 'equirectangular',
                            pole: str = 'north',
                            dataroot: Path = ACCESS_ROOT):
    '''
    Return the filename of the resampled file for a given orbit and channel.
    '''

    orbit_lower, orbit_upper = get_SSMI_orbit_range(orbit)
    if satellite == 'SSMI':
        orbit_dir = dataroot.joinpath(
        f"{satellite}_tb_orbits",f"f{ksat:02d}", f"r{orbit_lower:06d}_{orbit_upper:06d}"
        )
    elif satellite == 'AMSR2':
        orbit_dir = dataroot.joinpath(
        f"{satellite}_tb_orbits", f"r{orbit_lower:05d}_{orbit_upper:05d}"
        )
    else:
        raise ValueError(f'Satellite {satellite} not valid')

    if isinstance(channel,int):
        if ((channel >= 0) and (channel <=99)):
            channel_str = f'ch{channel:02d}'
        else:
            raise ValueError('Channel out of range')
    elif isinstance(channel,str):
        channel_str = channel
    else:
        raise ValueError('Channel type not valid')


    if grid_type == 'equirectangular':
        if satellite == 'SSMI':
            if channel_str == "time":
                filename = orbit_dir / f"r{orbit:06d}.time.nc"
            else:
                filename = (
                    orbit_dir / f"r{orbit:06d}.grid_tb.{channel_str}.{target_size:03d}km.nc"
                )
            return filename
        elif satellite == 'AMSR2':
            if channel_str == "time":
                filename = orbit_dir / f"r{orbit:05d}.polar_grid_time.{pole}.nc"
            else:
                filename = (
                    orbit_dir / f"r{orbit:05d}.grid_tb.{channel_str}.{target_size:03d}km.nc"
                )
            return filename
        else:
            raise ValueError(f'Satellite {satellite} not valid')
    elif grid_type == 'ease2':
        if pole in ['north','south']:
            if satellite == 'SSMI':
                if channel_str == "time":
                    filename = orbit_dir / f"r{orbit:06d}.time.nc"
                else:
                    filename = (
                        orbit_dir / f"r{orbit:06d}.grid_tb.{channel_str}.{target_size:03d}km.nc"
                    )
            elif satellite == 'AMSR2':
                if channel_str == "time":
                    filename = orbit_dir / f"r{orbit:05d}.polar_grid_time.{pole}.nc"
                else:
                    filename = (
                        orbit_dir / f"r{orbit:05d}.polar_grid_tb.{pole}."
                                    f"{channel_str}.{target_size:03d}km.nc"
                    )
            else:
                raise ValueError(f'Satellite {satellite} not valid')
        else:
            raise ValueError(f'Pole {pole} must be north or south')
        return filename
    else:
        raise ValueError(f'Grid type {grid_type} not valid')

def get_SSMI_orbit_range(orbit: int):
    """Return the lower/upper bounds to an orbit.

    >>> get_orbit_range(1)
    (1, 5000)
    >>> get_orbit_range(5000)
    (1, 5000)
    >>> get_orbit_range(5001)
    (5001, 10000)
    """
    BIN_WIDTH = 5000
    j = int((orbit - 1) / BIN_WIDTH)
    orbit_lower = 1 + j * BIN_WIDTH
    orbit_upper = (j + 1) * BIN_WIDTH
    return orbit_lower, orbit_upper

def count_processes_with_name(name: str):
    processes = os.popen('wmic process get description, processid').read()
    return processes.count(name)

def start_orbit_resample(sat_name,
                        target_size,
                        orbit_num,
                        ksat=13,
                        path_to_exe = '',
                        do_global = True,
                        do_polar=True,
                        pole_list = ['north','south'],
                        minimized = True):
    '''
    Start the orbit resample process
    input variables:
        sat_name: 'SSMI' or 'AMSR2'
        target_size: 30 or 70
        orbit_num: orbit number to resample
        ksat: 8-15 for SSMI, not needed for AMSR2
        path_to_exe: path to the resample executable
        do_global: True or False
        do_polar: True or False
        pole_list: list of poles to resample
    returns:
        p: subprocess object
    '''
    do_global_int = 0
    if do_global:
        do_global_int = 1

    do_NP_int = 0
    do_SP_int = 0
    if do_polar:
        if 'north' in pole_list:
            do_NP_int = 1
        if 'south' in pole_list:
            do_SP_int = 1

    window_name = f'{orbit_num} {target_size}km'
    min_str = ''
    if minimized:
        min_str = '/min'
    if sat_name == 'AMSR2':
        path_to_exe = path_to_AMSR2_exe
        if target_size == 30:
            # don't do the 6.9 and 7.2 channels
            command_str = f'start {min_str} "{window_name}" cmd /C {path_to_exe} {orbit_num} {orbit_num} {target_size} {do_global_int} {do_NP_int} {do_SP_int} 5 6 7 8 9 10 11 12'
        elif target_size == 70:
            command_str = f'start {min_str} "{window_name}" cmd /C {path_to_exe} {orbit_num} {orbit_num} {target_size} {do_global_int} {do_NP_int} {do_SP_int} 1 2 3 4 5 6 7 8 9 10 11 12'
        else:
            raise ValueError(f'target_size: {target_size} is not valid')
    elif sat_name == 'SSMI':
        path_to_exe = path_to_SSMI_exe
        num_overlap_scans=350
        if target_size == 70:
            command_str = f'start {min_str} "{window_name}" cmd /C {path_to_exe} {orbit_num} {orbit_num} {ksat} {target_size} {num_overlap_scans} {do_global_int} {do_NP_int} {do_SP_int} 1 2 3 4 5 6'
        else:
            raise ValueError(f'target_size: {target_size} is not valid')
    else:
        raise ValueError(f'Satellite {sat_name} not valid')
    print(command_str)
    
    p = subprocess.Popen(command_str,
                         stdin=subprocess.PIPE,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT,
                         close_fds=False,
                         shell=True)

    return p

def get_file_times(satellite_name,
                    start_orbit,
                    end_orbit,
                    target_size,
                    ksat=13,
                    grid_type='equirectangular',
                    pole=None):

    orbit_list = []
    file_times = []

    for test_orbit in range(start_orbit,end_orbit+1):

        filename = get_resampled_file_name(satellite=satellite_name,
                ksat=ksat,
                channel=5,
                target_size=target_size,
                orbit=test_orbit,
                grid_type = grid_type,
                pole=pole)

        if filename.is_file():
            if test_orbit%100 == 0:
                print(f'{filename} exists')
            file_time = filename.stat().st_mtime
            file_times.append(file_time)
            orbit_list.append(test_orbit)

    return orbit_list,file_times
        
    
if __name__ =='__main__':
    import numpy as np
    import matplotlib.pyplot as plt
    import scipy.stats
    import datetime

    
    sat_name = 'AMSR2'
    ksat=13
    grid_type = 'equirectangular'
    target_size = 70
    pole=None
    
    start_orbit = 56255
    end_orbit = 65700

    print('Getting file modification times')
    orbit_list,file_times = get_file_times(sat_name,
                                       start_orbit,
                                       end_orbit,
                                       target_size,
                                       ksat=ksat,
                                       grid_type=grid_type)

    file_times = np.array(file_times,dtype=np.float64)
    orbit_list = np.array(orbit_list,dtype=np.int32)

    file_times = file_times - file_times[0]

    fig,ax = plt.subplots()
    ax.plot(file_times,orbit_list)
    ax.set_xlabel('Time when file was written (s)')
    ax.set_ylabel('Orbit Number Processed')

    fit = scipy.stats.linregress(file_times[-200:-1],orbit_list[-200:-1])
    fitted_line = fit.intercept + fit.slope*(file_times[-200]+np.arange(0,350000,1000))
    ax.plot((file_times[-200]+np.arange(0,350000,1000)),fitted_line,'r--',label=f'Orbit = {fit.intercept:.1f} + {fit.slope:.1f} * Time')
    
    est_time_left = (end_orbit-orbit_list[-1])/fit.slope
    est_ending_time = file_times[-1]+est_time_left
    ax.plot((0,est_ending_time),(end_orbit,end_orbit),'g--',label='End Orbit')

    print(f'Estimated orbits per hour: {fit.slope*3600:.1f}')
    print(f'Estimated time to finish: {est_time_left/3600:.1f} hours')
    print(f'Estimated ending time/date: {datetime.datetime.today()+datetime.timedelta(seconds=est_time_left)}')
    plt.show()
    print('Done')





    





