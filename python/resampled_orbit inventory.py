import numpy as np
from pathlib import Path
import os
import matplotlib.pyplot as plt

def get_orbit_range(orbit):
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

def check_for_missing_amsr2_tb_orbit_files(pathname_tb_output,
                               channel_list,iorbit):

    orbit_lower, orbit_upper = get_orbit_range(iorbit)
    orbit_dir = Path(pathname_tb_output) / f"r{orbit_lower:05d}_{orbit_upper:05d}"
    
    #check time file
    filename = orbit_dir / f"r{iorbit:05d}.time.nc"

    if not(os.path.exists(filename)):
        return True

def check_for_missing_l1a(pathname_l1a,iorbit):

    orbit_lower, orbit_upper = get_orbit_range(iorbit)
    orbit_dir = Path(pathname_l1a) / f"r{orbit_lower:05d}_{orbit_upper:05d}"
    
    #check l1a file
    filename = orbit_dir / f"r{iorbit:05d}.gz"
    if not(os.path.exists(filename)):
        return True
    else:
        return False

if __name__ == '__main__':

    pathname_tb_output = Path('L:/access/amsr2_tb_orbits/')
    pathname_l1a = Path('J:/AMSR2/L1A/')
    channel_list = range(5,13)

    processed = np.zeros((52000),dtype = np.int32)
    processed.fill(-9)
    
    for iorbit in range(0,52000):
        if iorbit%1000 == 1:
            print(iorbit)
        missing = check_for_missing_amsr2_tb_orbit_files(pathname_tb_output,channel_list,iorbit)
        if missing:
            l1a_missing = check_for_missing_l1a(pathname_l1a,iorbit)
            if l1a_missing:
                 processed[iorbit] = 0
            else:
                processed[iorbit] = -1
                print(f'Missing orbit with L1A exists: {iorbit}')
        else:
            processed[iorbit] = 1

    fig,ax = plt.subplots()
    ax.plot(processed)
    ax.set_ylim(-2.0,2.0)
    num_complete = np.sum(processed == 1)
    num_missing = np.sum(processed == 0)
    num_incomplete = np.sum(processed == -1)

    print(f'Num Complete: {num_complete}')
    print(f'Num To Go:    {num_incomplete}')
    print(f'Num Missing:  {num_missing}')
    plt.show()
    print
