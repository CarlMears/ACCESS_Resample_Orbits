import os
from typing import List
from ctypes import WinDLL, Structure, POINTER, pointer, c_int32, c_bool

from FortranTools.python_types import iop_type, py_string

class AMSR2_Lib:
    DLL_SEARCH_PATHS = [
        r"M:\job_access\resampling\fortran_python\build\amsr2.dll"
    ]

    DLL_FUNC_NAMES = [
        "process_and_resample_amsr2_l2b_orbit_PY_HOOK"
    ]

    class process_and_resample_amsr2_l2b_orbit_call(Structure):
        _fields_ = [
            ("iorbit", c_int32), 
            ("footprint_size", c_int32),
            ("channel_list", c_int32 * 20),
            ("write_swath_tbs", c_bool),
            ("write_swath_geoloc", c_bool),
            ("overwrite_existing_files", c_bool),
            ("pathname_tb_output", iop_type)
        ]

    def __init__(self):
        self.dll = None
        for path in AMSR2_Lib.DLL_SEARCH_PATHS:
            try:
                # Ugly hack to force proper DLL loading order/paths
                # Have a programatic approach in the works.
                WinDLL(r"X:\fortran_tools\builds\win\XE_v19.1.0.166\third_party\netcdf4.7.4\bin\netcdf.dll")
                WinDLL(r"X:\fortran_tools\builds\win\XE_v19.1.0.166\third_party\zlib-1.2.11\bin\zlibd.dll")
                WinDLL(r"C:\Program Files (x86)\IntelSWTools\parallel_studio_xe_2020.0.075\compilers_and_libraries_2020\windows\redist\intel64\compiler\libifcoremd.dll")
                WinDLL(r"C:\Program Files (x86)\IntelSWTools\parallel_studio_xe_2020.0.075\compilers_and_libraries_2020\windows\redist\intel64\compiler\libifcoremdd.dll")
                WinDLL(r"C:\Program Files (x86)\Common Files\Intel\Shared Libraries\redist\intel64_win\compiler\libifportmd.dll")
                WinDLL(r"C:\Program Files (x86)\IntelSWTools\parallel_studio_xe_2020.0.075\compilers_and_libraries_2020\windows\redist\intel64\compiler\libmmdd.dll")
                self.dll = WinDLL(path)
                break
            except:
                if os.path.exists(path):
                    from FortranTools.dll_tools import find_problem_dll
                    find_problem_dll(path)
                else:
                    raise
            finally:
                if self.dll is None: 
                    count = 0
                    while self.dll is None:
                        self.dll = WinDLL(path)
                        count += 1
                        if count > 50:
                            raise RuntimeError(
                                f'Did not find logger DLL in {AMSR2_Lib.DLL_SEARCH_PATHS}')

                self._process_amsr2_hook = self.dll.process_and_resample_amsr2_l2b_orbit_PY_HOOK
                self._process_amsr2_hook.argtypes = [POINTER(AMSR2_Lib.process_and_resample_amsr2_l2b_orbit_call)]
                self._process_amsr2_hook.restype = None

    def process_and_resample_amsr2_l2b_orbit(
        self,
        iorbit: int,
        footprint_size: int,
        channel_list: List[int],
        write_swath_tbs: bool,
        write_swath_geoloc: bool,
        overwrite_existing_files: bool,
        pathname_tb_output: str
    ) -> None: 
        call_args = AMSR2_Lib.process_and_resample_amsr2_l2b_orbit_call()
        call_args.iorbit = c_int32(iorbit)
        call_args.footprint_size = c_int32(footprint_size)
        call_args.channel_list = ( c_int32 * 20)(*channel_list)
        call_args.write_swath_tbs = c_bool(write_swath_tbs)
        call_args.write_swath_geoloc = c_bool(write_swath_geoloc)
        call_args.overwrite_existing_files = c_bool(overwrite_existing_files)
        call_args.pathname_tb_output = py_string.serialize(pathname_tb_output)
        self._process_amsr2_hook(pointer(call_args))
