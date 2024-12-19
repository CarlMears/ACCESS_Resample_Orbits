! This module contains the MWI (pseudo-) TDR derived type and the
! routines to read and write it via netCDF
module L2_Tb_io
  use, intrinsic :: iso_fortran_env, only: int8, int32, real32, real64
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_finite
  use NAN_support, only: nan_f32, nan_f64
    
  !use simulation_param, only: SimulationParameterData_f, FOVS_PER_SCAN, MAXSCAN, NBAND, NPOLBAND, NEXTRAFOV,NEXTRASCAN
  
  !use satellite_constants, only: nChannel,nBand,nPol !,nFOV,nScan
  use sat_id_module, only: maxcel,maxscan,nChannel
  use io_nc, only: handle_nc_err, minmax_iso8601
  use netcdf

  !use version, only: RESAMPLE_VERSION, RESAMPLE_VERSION_DATE
  implicit none
  private
  public :: write_geoloc_extra_netcdf

contains


  ! -------------------------------------------------------------------------------
  ! Create a file with resampled brightness temperatures in swath geometry
  ! -------------------------------------------------------------------------------
  subroutine write_geoloc_extra_netcdf(filename_geo,orbit_number,overwrite)

    use ssmi_l1c, only: scan_time,cellat,cellon
    !use l2_module, only: scan_time,cellat,cellon
    use l2_module_extra

    character(len=*), intent(in) :: filename_geo
    integer(4),intent(in)        :: orbit_number
    integer(4),intent(in)        :: overwrite

    integer :: ncid, varid
    integer :: i

    character(len=80) :: nc_version_full
    character(len=10) :: nc_version
    character(len=5) :: time_zone
    integer, dimension(8) :: time_vals
    character(len=24) :: timestamp
    character(len=20) :: time_start, time_end
    character(len=:), allocatable :: history
    !integer :: cmdline_len

    integer :: dim_scan, dim_fov, dim_band, dim_pol, dim_chan, dim_xyz, dim_q

    nc_version_full = trim(nf90_inq_libvers())
    ! Only use the first whitespace-delimited word, which is the version number
    nc_version = nc_version_full(1:index(nc_version_full, " "))

    history = 'unknown'
    call date_and_time(zone=time_zone, values=time_vals)
    write(timestamp, '(I4, "-", I2.2, "-", I2.2, " ", I2.2, ":", I2.2, ":", I2.2, A5)') &
         time_vals(1), time_vals(2), time_vals(3), time_vals(5), time_vals(6), time_vals(7), time_zone

    !call minmax_iso8601(scan_time_utc, tdr_data%leap_seconds, time_start, time_end)
         !not sure tis is exactly correct
    call minmax_iso8601(scan_time, time_start, time_end)
    print *,'time_start: ',time_start
    print *,'time_end: ',time_end
    print *,'Writing: ',filename_geo
    if (overwrite == 1) then
        call handle_nc_err(nf90_create(filename_geo, ior(NF90_CLOBBER, NF90_NETCDF4), ncid))
    else
        call handle_nc_err(nf90_create(filename_geo, ior(NF90_NOCLOBBER, NF90_NETCDF4), ncid))
    endif

    ! Define global attributes
    call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "Conventions", "CF-1.8,ACDD-1.3"))
    call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "title", "AMSR2 Geolocation File"))
    call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "institution", "REMSS"))
    call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "history", trim(history)))
    call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "netcdf_version_id", trim(nc_version)))
    !call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "git_revision", trim(RESAMPLE_VERSION)))
    !call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "git_date", trim(RESAMPLE_VERSION_DATE)))
    call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "date_created", timestamp))
    call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "source", "AMSR2"))
    call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "platform", "GCOM-W"))
    call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "sensor", "AMSR2"))
    call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "creator_name", "Remote Sensing Systems"))
    call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "creator_email", "support@remss.com"))
    call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "creator_url", "http://www.remss.com"))
    call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "orbit", orbit_number))
    call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "geospatial_lat_min", minval(cellat)))
    call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "geospatial_lat_max", maxval(cellat)))
    call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "geospatial_lon_min", minval(cellon)))
    call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "geospatial_lon_max", maxval(cellon)))
    call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "time_coverage_start", time_start))
    call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "time_coverage_end", time_end))

    ! Define dimensions
    call handle_nc_err(nf90_def_dim(ncid, "scan", maxscan_extra, dim_scan))
    call handle_nc_err(nf90_def_dim(ncid, "fov",  maxcel_extra, dim_fov))
    call handle_nc_err(nf90_def_dim(ncid, "chan", nChannel, dim_chan))
    
    
    !Define and write coordinate variables
    !use satellite_constants, only: nChannel,nBand,nPol,nFOV,nScan
  
    call handle_nc_err(nf90_def_var(ncid, "scan", NF90_INT, dim_scan, varid))
    call handle_nc_err(nf90_put_var(ncid, varid, [(i, i = 1, maxscan_extra)]))

    call handle_nc_err(nf90_def_var(ncid, "fov", NF90_INT, dim_fov, varid))
    call handle_nc_err(nf90_put_var(ncid, varid, [(i, i = 1, maxcel_extra)]))

    call handle_nc_err(nf90_def_var(ncid, "chan", NF90_INT, dim_chan, varid))
    call handle_nc_err(nf90_put_var(ncid, varid, [(i, i = 1, nChannel)]))
    
    !print *,maxscan_extra,maxcel_extra,nChannel
    ! print *,size(zang_extra)
    ! Define and write main variables

    ! call handle_nc_err(nf90_def_var(ncid, "zang", NF90_DOUBLE, dim_scan, varid, &
    !  deflate_level=2, shuffle=.true.))
    ! call handle_nc_err(nf90_put_var(ncid, varid, zang_extra))

    call handle_nc_err(nf90_def_var(ncid, "scan_time", NF90_DOUBLE, dim_scan, varid, &
          deflate_level=2, shuffle=.true.))
    call handle_nc_err(nf90_def_var_fill(ncid, varid, 0, -999999999.0D0))
    call handle_nc_err(nf90_put_var(ncid, varid, scan_time_extra))

    ! call handle_nc_err(nf90_def_var(ncid, "scan_index", NF90_DOUBLE, dim_scan, varid, &
    !      deflate_level=2, shuffle=.true.))
    ! call handle_nc_err(nf90_put_var(ncid, varid, scan_index_extra))

    ! call handle_nc_err(nf90_def_var(ncid, "scan_rate", NF90_FLOAT, dim_scan, varid, &
    !      deflate_level=2, shuffle=.true.))
    ! call handle_nc_err(nf90_put_var(ncid, varid, omega_extra))


    call handle_nc_err(nf90_def_var(ncid, "latitude", NF90_FLOAT, [dim_fov, dim_scan], varid, &
         deflate_level=2, shuffle=.true.))
    call handle_nc_err(nf90_def_var_fill(ncid, varid, 0, -999.0))
    call handle_nc_err(nf90_put_var(ncid, varid, cellat_extra))

    call handle_nc_err(nf90_def_var(ncid, "longitude", NF90_FLOAT, [dim_fov, dim_scan], varid, &
         deflate_level=2, shuffle=.true.))
    call handle_nc_err(nf90_def_var_fill(ncid, varid, 0, -999.0))
    call handle_nc_err(nf90_put_var(ncid, varid, cellon_extra))

    ! call handle_nc_err(nf90_def_var(ncid, "range", NF90_FLOAT, [dim_fov, dim_scan], varid, &
    ! deflate_level=2, shuffle=.true.))
    ! call handle_nc_err(nf90_put_var(ncid, varid, celrange_extra))

    ! call handle_nc_err(nf90_def_var(ncid, "incidence", NF90_FLOAT, [dim_fov, dim_scan], varid, &
    !      deflate_level=2, shuffle=.true.))
    ! call handle_nc_err(nf90_put_var(ncid, varid, celtht_extra))

    call handle_nc_err(nf90_def_var(ncid, "azimuth", NF90_FLOAT, [dim_fov, dim_scan], varid, &
         deflate_level=2, shuffle=.true.))
    call handle_nc_err(nf90_def_var_fill(ncid, varid,0 , -999.0))
    call handle_nc_err(nf90_put_var(ncid, varid, celphi_extra))

    print *,'writing tb'
    print *,tbr_extra(5,100,:)
    call handle_nc_err(nf90_def_var(ncid, "tb_resampled", NF90_FLOAT, [dim_fov, dim_scan, dim_chan], varid, &
         deflate_level=2, shuffle=.true.))
    call handle_nc_err(nf90_def_var_fill(ncid, varid, 0, -999.0))
    call handle_nc_err(nf90_put_var(ncid, varid, tbr_extra))

    ! Set attributes on datasets
    ! call handle_nc_err(nf90_inq_varid(ncid, "scan_index", varid))
    ! call handle_nc_err(nf90_put_att(ncid, varid, "long_name", "scan index"))
    ! call handle_nc_err(nf90_put_att(ncid, varid, "axis", "Y"))

    call handle_nc_err(nf90_inq_varid(ncid, "latitude", varid))
    call handle_nc_err(nf90_put_att(ncid, varid, "standard_name", "latitude"))
    call handle_nc_err(nf90_put_att(ncid, varid, "units", "degrees_north"))

    call handle_nc_err(nf90_inq_varid(ncid, "longitude", varid))
    call handle_nc_err(nf90_put_att(ncid, varid, "standard_name", "longitude"))
    call handle_nc_err(nf90_put_att(ncid, varid, "units", "degrees_east"))

    call handle_nc_err(nf90_inq_varid(ncid, "azimuth", varid))
    call handle_nc_err(nf90_put_att(ncid, varid, "long_name", "azimuth angle"))
    call handle_nc_err(nf90_put_att(ncid, varid, "standard_name", "sensor_azimuth_angle"))
    call handle_nc_err(nf90_put_att(ncid, varid, "units", "degree"))

    call handle_nc_err(nf90_inq_varid(ncid, "scan_time", varid))
    call handle_nc_err(nf90_put_att(ncid, varid, "long_name", "scan time"))
    call handle_nc_err(nf90_put_att(ncid, varid, "standard_name", "time"))
    call handle_nc_err(nf90_put_att(ncid, varid, "units", "seconds since 2000-01-01 00:00:00 UTC"))

    ! call handle_nc_err(nf90_inq_varid(ncid, "incidence", varid))
    ! call handle_nc_err(nf90_put_att(ncid, varid, "long_name", "Earth incidence angle"))
    ! call handle_nc_err(nf90_put_att(ncid, varid, "standard_name", "sensor_view_angle"))
    ! call handle_nc_err(nf90_put_att(ncid, varid, "units", "degree"))

    ! call handle_nc_err(nf90_inq_varid(ncid, "scan_rate", varid))
    ! call handle_nc_err(nf90_put_att(ncid, varid, "long_name", "antenna revolutions per minute"))
    ! call handle_nc_err(nf90_put_att(ncid, varid, "units", "rpm"))

    ! call handle_nc_err(nf90_inq_varid(ncid, "zang", varid))
    ! call handle_nc_err(nf90_put_att(ncid, varid, "long_name", "spacecraft orbit angle"))
    ! call handle_nc_err(nf90_put_att(ncid, varid, "units", "degree"))

    call handle_nc_err(nf90_close(ncid))
  end subroutine write_geoloc_extra_netcdf

end module L2_Tb_io
