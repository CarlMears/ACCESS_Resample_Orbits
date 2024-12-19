
# file = 'M:/job_access/resampling/fortran_python/build/make_bat_file_ops1d.89.bat'

# with open(file,'w') as f:
#     for i in range(52200,53000,30):
#         cmd = f'start "{i} 30km" cmd /K resample_amsr2_exe_89_4.exe {i} {i+29} 30 1 1 1 13 14\n'
#         f.write(cmd)
#         cmd = 'timeout 10\n'
#         f.write(cmd)


file = 'M:/job_access/resampling/fortran_python/build/make_bat_file_smops.89.bat'

with open(file,'w') as f:
    for i in range(28971,28990,3):
        cmd = f'start "{i} 30km" cmd /K resample_amsr2_exe_89_4.exe {i} {i+2} 30 1 1 1 13 14\n'
        f.write(cmd)
        cmd = 'timeout 10\n'
        f.write(cmd)