project('RSS ACCESS Project for SSMI', ['fortran'], version : '1.0.0')

# Setup global arguments.
dialect = ['/nologo', '/free', '/fpp', '/debug:full']
bounds = '/CB'

fortran_args = ['/assume:byterecl']

add_global_arguments(dialect, language : 'fortran')
add_global_arguments(bounds, language : 'fortran')

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

l2_module_src       = ['amsr2/l2_tb_module_ACCESS.f90']

land_ice_src       = ['../land_masks/percent_land.f90',
                      '../land_masks/fdice.f90']

system_src            = ['../syspgm/get_cmd_arg.f90',
                         '../syspgm/open_file_routines.f90',
                         '../syspgm/fileexist.f90',
                        ]

nan_src = ['../compatibility/NAN_support_ieee.f90']

# amsr2_sat_id_module_src = ['../amsr2/amsr_id_module.f90']
# amsr2_constants_src = ['../amsr2/satellite_constants_amsr2.f90']
# amsr2_l2_module_src = ['../amsr2/l2_module_amsr_v8.f90']
# amsr2_l2_tb_module_ACCESS_src = ['../amsr2/l2_tb_module_ACCESS.f90']

io_nc_src    =   ['../io_nc.f90', 
                  '../CF_metadata.f90',
                  'ssmi_L2_tb_io.f90',]

geolocation_src = ['../geolocation/ecef_conversions.f90',
                   '../geolocation/local_proj.f90',
                   'extra_locs.f90',]

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

ssmi_internal_src   =    [
                         'ssmi_filename_routines.f90',
                         'module_ssmi_l1c.f90',
                         'ssmi_id_module.f90',
                         ] 

resample_src = ['../resample/resample_tb_rect_map.f90',
                '../resample/resample_tb_polar_map.f90',
                '../resample/wgs84.f90']

quad_interp_src = ['../quad_interp/QuadrilateralUtilties.f90',
                   '../quad_interp/QuadrilateralInterp.f90',
                   '../quad_interp/QuadrilateralTransform.f90']

ssmi_exe_src = ['process_ssmi_l2b_orbit_tb_ACCESS.f90']

test_quad_source = ['../quad_interp/test_quadrilateral.f90']

ssmi_src = [ssmi_exe_src,
            ssmi_internal_src,
            system_src,
            nan_src,
            io_nc_src,
            geolocation_src,
            math_src,
            quad_interp_src,
            date_utilities_src,
            resample_src,
            gridding_src,
            polar_gridding_src,
            ]

ssmi_exe = executable('resample_ssmi_exe', 
                       ssmi_src,
   dependencies : [
       nc_f_dep,
    #   log_dep,
    #   legacy_gzip_dep,
    #   py_integrations_dep
   ]
)


test_quad_exe = executable('test_quad_exe', 
                           [test_quad_source,
                            quad_interp_src],
   dependencies : [
#      nc_f_dep,
#      log_dep,
#      legacy_gzip_dep,
#      py_integrations_dep
   ]
)
