from pathlib import Path
import os

root_path = Path('B:/_access_temp/ssmi_tb_orbits/')

for region_str in ['GL','NP','SP']:
    for target_size in [70]:
        for ksat in [14]:
            for orbit_set in range(0,15):
                start_orbit = orbit_set * 5000 + 1
                end_orbit = start_orbit + 4999
                path =root_path / f'f{ksat:02d}' / f'{region_str}_{target_size:02d}' / f'r{start_orbit:06d}_{end_orbit:06d}'
                print(path)
                os.makedirs(path,exist_ok=True)






