project('RSS ACCESS Project for AMSR2', ['fortran'], version : '1.0.0')

# Setup global arguments.
dialect = ['/free', '/fpp'] #, '/debug:full']
bounds = '/CB'

fortran_args = ['/assume:byterecl']

add_global_arguments(dialect, language : 'fortran')

# Find dependencies.
nc_f_dep = dependency('netcdf-fortran-xe20')
log_dep = dependency('logger-xe20')
legacy_gzip_dep = dependency('legacy-gzip-xe20')

# Add source.
date_utilities_src  = ['../date_utilities/date_utilities.f90',
                       '../date_utilities/time_conversions.f90',
                       '../date_utilities/time_conversions_extra.f90']

gridding_src = ['../gridding/earth_ref_maps.f90',
                '../gridding/map_type_constants.f90']
                       
l2_module_src       = ['l2_tb_module_ACCESS.f90']

land_ice_src       = ['../land_masks/percent_land.f90',
                      '../land_masks/fdice.f90']

math_src            = ['../math/interpolate.f90',
                       '../math/math_routines.f90',
                       '../math/quaternion_routines.f90']

nc_src              =   ['../io_nc.f90',
                         '../CF_metadata.f90']

amsr2_src          =    ['satellite_constants_amsr2.f90',
                         'l2_module_amsr_v8.f90',
                         'amsr2_L2_Tb_io.f90',
                         'RSS_L1A.f90',
                         'filename_routines.f90',
                         'amsr2_define_filenames_v8.f90',
                         'amsr2_data_files_v8.f90',
                         'amsr_geolocation.f90',
                         'amsr_id_module.f90',
                         'io_filenames_v8.f90',
                         'AMSR2_resample_version.f90',
                         ]

system_src            = ['../syspgm/get_cmd_arg.f90',
                         '../syspgm/open_file_routines.f90',
                         '../syspgm/fileexist.f90',
                        ]

nan_src = ['../compatibility/NAN_support_ieee.f90']

amsr2_sat_id_module_src = ['amsr_id_module.f90']
amsr2_constants_src = ['satellite_constants_amsr2.f90']
amsr2_l2_module_src = ['l2_module_amsr_v8.f90']
amsr2_l2_tb_module_ACCESS_src = ['l2_tb_module_ACCESS.f90']

# These are mostly unchanged from franks versions.
counts_to_tb_src      = ['counts_to_tb/avg_cal_counts.f90',
                         'counts_to_tb/ckocntb.f',   
                         'counts_to_tb/convert_ta_to_tb.f90',
                         'counts_to_tb/fd_cold_mirror_rfi_v8.f90',
                         'counts_to_tb/fd_geostationary_latlon.f90',
                         'counts_to_tb/fd_leap_sec.f',
                         'counts_to_tb/fd_monthly_wind.f90',
                         'counts_to_tb/fd_ta_adjusted_v8.f90',
                         'counts_to_tb/fd_tbsss_clm.f90',
                         'counts_to_tb/fd_tcold.f90',
                         'counts_to_tb/fd_thot.f',
                         'counts_to_tb/field_routines.f',
                         'counts_to_tb/find_therm.f',
                         'counts_to_tb/moon_in_cold_view.f',
                         'counts_to_tb/pre_nut_routines.f90',
                         'counts_to_tb/resample_89.f',
                         'counts_to_tb/resample_ta.f90',
                         'counts_to_tb/reverse_agc.f90',
                         'counts_to_tb/rfi_flag.f90',
                         'counts_to_tb/spillover_routines.f90',
                         'counts_to_tb/tatb_amsr2.f90',
                         'counts_to_tb/tatb_routines_ops.f',
                        ]

io_nc_src              =   ['../io_nc.f90', '../CF_metadata.f90']

gridding_src = ['../gridding/equirectangular_maps.f90',
                '../gridding/equirectangular_maps64.f90',
                '../gridding/equirectangular_map_constants.f90',
               ]

polar_gridding_src = ['../gridding/read_NSIDC_polar_maps.f90',
                     '../gridding/earth_ref_maps.f90',
                     '../gridding/ease2.f90',
                     '../gridding/init_polar_grids.f90',
                     '../gridding/map_type_constants.f90',
                ]

land_ice_src       = ['../land_masks/percent_land.f90',
                      '../land_masks/fdice.f90']

math_src            = ['../math/interpolate.f90',
                       '../math/math_routines.f90',
                       '../math/quaternion_routines.f90',
                       '../math/trig_degrees.f90']

sun_moon_src          = ['sun_moon/sun_moon.f90']
                        
amsr2_internal_src          =    [
                         'amsr2_L2_Tb_io.f90',
                         'RSS_L1A.f90',
                         'filename_routines.f90',
                         'amsr2_define_filenames_v8.f90',
                         'amsr2_data_files_v8.f90',
                         'amsr_geolocation.f90',
                         'io_filenames_v8.f90',
                         'AMSR2_resample_version.f90',
                         ]

amsr2_entry_src = ['amsr2_access_entry.f90']

resample_src = ['../resample/resample_ta_CM.f90',
                '../resample/resample_ta_polar_map.f90',
                '../resample/wgs84.f90']


quad_interp_src = ['../quad_interp/QuadrilateralUtilties.f90',
                   '../quad_interp/QuadrilateralInterp.f90',
                   '../quad_interp/QuadrilateralTransform.f90']

amrs2_lib_src = ['amsr2_access_lib.f90']

python_interfaces_src = ['../python/interfaces/f90/access_py_interfaces.f90']

amsr2_exe_src = ['process_amsr2_l2b_orbit_tb_ACCESS.f90']

test_quad_source = ['../quad_interp/test_quadrilateral.f90']
test_map_source = ['../gridding/test_new_maps.f90']

AMSR2_src = [date_utilities_src,
            gridding_src,
            l2_module_src,
            land_ice_src,
            math_src,
            nc_src,
            resample_src,
            amsr2_src,
            system_src,
            nan_src,
            amsr2_sat_id_module_src,
            amsr2_constants_src,
            amsr2_l2_module_src,
            amsr2_l2_tb_module_ACCESS_src,
            counts_to_tb_src,
            io_nc_src,
            polar_gridding_src,
            math_src,
            sun_moon_src,
            amsr2_internal_src,
            resample_src,
            quad_interp_src,
            amsr2_exe_src
            ]

amsr2_exe = executable('resample_amsr2_exe_89_4', 
   AMSR2_src,
   dependencies : [
      nc_f_dep,
      legacy_gzip_dep,
   ]
)

test_quad_exe = executable('test_quad_exe', 
   [test_quad_source,
   quad_interp_src],
)
