! This module contains the MWI (pseudo-) TDR derived type and the
! routines to read and write it via netCDF
module L2_Tb_io
  use, intrinsic :: iso_fortran_env, only: int8, int32, real32, real64
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_finite
  use NAN_support, only: nan_f32, nan_f64
    
  !use simulation_param, only: SimulationParameterData_f, FOVS_PER_SCAN, MAXSCAN, NBAND, NPOLBAND, NEXTRAFOV,NEXTRASCAN
  use satellite_constants, only: nChannel,nBand,nPol !,nFOV,nScan
  use sat_id_module, only: maxcel,maxscan
  use io_nc, only: handle_nc_err, minmax_iso8601
  use netcdf

  use version, only: RESAMPLE_VERSION, RESAMPLE_VERSION_DATE
  implicit none
  private
  public :: write_tb_native_netcdf, write_geoloc_extra_netcdf

contains

  ! -------------------------------------------------------------------------------
  ! Create a new Tb file
  subroutine write_tb_native_netcdf(filename_tb,orbit_number,overwrite)

    use l2_module

    character(len=*), intent(in) :: filename_tb
    integer(4),intent(in)        :: orbit_number
    integer(4),intent(in)        :: overwrite

     integer :: ncid, varid
     integer :: i,jscan1,jscan2
    
     character(len=80) :: nc_version_full
     character(len=10) :: nc_version
    character(len=5) :: time_zone
    integer, dimension(8) :: time_vals
    character(len=24) :: timestamp
    character(len=20) :: time_start, time_end
    character(len=:), allocatable :: history
    !integer :: cmdline_len

    integer :: dim_scan, dim_fov, dim_band, dim_pol, dim_chan, dim_xyz, dim_q
    integer :: dim_scan_89, dim_fov_89, dim_chan_89
    do i = 1,maxscan
        if (.not. ieee_is_finite(scan_time(i))) then
            cellat(:,i) = nan_f32
            cellon(:,i) = nan_f32
            celtht(:,i) = nan_f32
            celphi(:,i) = nan_f32
            !this is copied from amsr_geolocation.f90
            jscan1=2*i-3  !89b assign to this to make it adjacent and before to 89a
            jscan2=2*i
            if (jscan1 > 0) then
                cellat_89(:,jscan1) = nan_f32
                cellon_89(:,jscan1) = nan_f32
                celtht_89(:,jscan1) = nan_f32
                celphi_89(:,jscan1) = nan_f32
                celrng_89(:,jscan1) = nan_f32
            endif
            cellat_89(:,jscan2) = nan_f32
            cellon_89(:,jscan2) = nan_f32 
            celtht_89(:,jscan2) = nan_f32
            celphi_89(:,jscan2) = nan_f32
            celrng_89(:,jscan2) = nan_f32
            ta(:,i,:)   = nan_f32
            scpos(:,i)  = nan_f64
            scvel(:,i)  = nan_f64
            scan_index(i) = nan_f64
            omega(i) = nan_f64
            zang(i) = nan_f64
            scrpy(:,i)  = nan_f64
        endif
    enddo

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

    if (overwrite == 1) then
        call handle_nc_err(nf90_create(filename_tb, ior(NF90_CLOBBER, NF90_NETCDF4), ncid))
    else
        call handle_nc_err(nf90_create(filename_tb, ior(NF90_NOCLOBBER, NF90_NETCDF4), ncid))
    endif

    ! Define global attributes
    call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "Conventions", "CF-1.8,ACDD-1.3"))
    call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "title", "AMSR2"))
    call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "institution", "REMSS"))
    !call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "history", trim(history)))
    call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "netcdf_version_id", trim(nc_version)))
    call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "git_revision", trim(RESAMPLE_VERSION)))
    call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "git_date", trim(RESAMPLE_VERSION_DATE)))
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
    call handle_nc_err(nf90_def_dim(ncid, "scan", maxscan, dim_scan))
    call handle_nc_err(nf90_def_dim(ncid, "fov",  maxcel, dim_fov))
    call handle_nc_err(nf90_def_dim(ncid, "scan89", maxscan_89, dim_scan_89))
    call handle_nc_err(nf90_def_dim(ncid, "fov89",  maxcel_89, dim_fov_89))
    !call handle_nc_err(nf90_def_dim(ncid, "band", nBand, dim_band))
    !call handle_nc_err(nf90_def_dim(ncid, "pol", nPol, dim_pol))
    call handle_nc_err(nf90_def_dim(ncid, "chan", nChannel, dim_chan))
    call handle_nc_err(nf90_def_dim(ncid, "chan89",2, dim_chan_89))
    
    call handle_nc_err(nf90_def_dim(ncid, "xyz", 3, dim_xyz))
    call handle_nc_err(nf90_def_dim(ncid, "q", 4, dim_q))

    !Define and write coordinate variables
    !use satellite_constants, only: nChannel,nBand,nPol,nFOV,nScan
  
    call handle_nc_err(nf90_def_var(ncid, "scan", NF90_INT, dim_scan, varid))
    call handle_nc_err(nf90_put_var(ncid, varid, [(i, i = 1, maxscan)]))

    call handle_nc_err(nf90_def_var(ncid, "fov", NF90_INT, dim_fov, varid))
    call handle_nc_err(nf90_put_var(ncid, varid, [(i, i = 1, maxcel)]))

    call handle_nc_err(nf90_def_var(ncid, "chan", NF90_INT, dim_chan, varid))
    call handle_nc_err(nf90_put_var(ncid, varid, [(i, i = 1, nChannel)]))

    call handle_nc_err(nf90_def_var(ncid, "chan89", NF90_INT, dim_chan_89, varid))
    call handle_nc_err(nf90_put_var(ncid, varid, [(i, i = 13, 14)]))

    call handle_nc_err(nf90_def_var(ncid, "xyz", NF90_INT, dim_xyz, varid))
    call handle_nc_err(nf90_put_var(ncid, varid, [(i, i = 1, 3)]))
 
    
    print *,maxscan
    print *,size(scpos)
    ! Define and write main variables
     call handle_nc_err(nf90_def_var(ncid, "scpos", NF90_DOUBLE, [dim_scan, dim_xyz], varid, &
          deflate_level=2, shuffle=.true.))
     call handle_nc_err(nf90_put_var(ncid, varid, transpose(scpos)))

     call handle_nc_err(nf90_def_var(ncid, "scvel", NF90_DOUBLE, [dim_scan, dim_xyz], varid, &
          deflate_level=2, shuffle=.true.))
     call handle_nc_err(nf90_put_var(ncid, varid, transpose(scvel)))

    call handle_nc_err(nf90_def_var(ncid, "zang", NF90_DOUBLE, dim_scan, varid, &
     deflate_level=2, shuffle=.true.))
    call handle_nc_err(nf90_put_var(ncid, varid, zang))

    call handle_nc_err(nf90_def_var(ncid, "scan_time", NF90_DOUBLE, dim_scan, varid, &
         deflate_level=2, shuffle=.true.))
    call handle_nc_err(nf90_put_var(ncid, varid, scan_time))

    call handle_nc_err(nf90_def_var(ncid, "scan_index", NF90_DOUBLE, dim_scan, varid, &
         deflate_level=2, shuffle=.true.))
    call handle_nc_err(nf90_put_var(ncid, varid, scan_index))

    call handle_nc_err(nf90_def_var(ncid, "scan_rate", NF90_FLOAT, dim_scan, varid, &
         deflate_level=2, shuffle=.true.))
    call handle_nc_err(nf90_put_var(ncid, varid, omega))

    call handle_nc_err(nf90_def_var(ncid, "scrpy", NF90_DOUBLE, [dim_scan, dim_xyz], varid, &
         deflate_level=2, shuffle=.true.))
    call handle_nc_err(nf90_put_var(ncid, varid, transpose(scrpy)))

    call handle_nc_err(nf90_def_var(ncid, "latitude", NF90_FLOAT, [dim_fov, dim_scan], varid, &
         deflate_level=2, shuffle=.true.))
    call handle_nc_err(nf90_put_var(ncid, varid, cellat))

    call handle_nc_err(nf90_def_var(ncid, "longitude", NF90_FLOAT, [dim_fov, dim_scan], varid, &
         deflate_level=2, shuffle=.true.))
    call handle_nc_err(nf90_put_var(ncid, varid, cellon))

    call handle_nc_err(nf90_def_var(ncid, "incidence", NF90_FLOAT, [dim_fov, dim_scan], varid, &
         deflate_level=2, shuffle=.true.))
    call handle_nc_err(nf90_put_var(ncid, varid, celtht))

    call handle_nc_err(nf90_def_var(ncid, "azimuth", NF90_FLOAT, [dim_fov, dim_scan], varid, &
         deflate_level=2, shuffle=.true.))
    call handle_nc_err(nf90_put_var(ncid, varid, celphi))

    call handle_nc_err(nf90_def_var(ncid, "latitude_89", NF90_FLOAT, [dim_fov_89, dim_scan_89], varid, &
         deflate_level=2, shuffle=.true.))
    call handle_nc_err(nf90_put_var(ncid, varid, cellat_89))

    call handle_nc_err(nf90_def_var(ncid, "longitude_89", NF90_FLOAT, [dim_fov_89, dim_scan_89], varid, &
           deflate_level=2, shuffle=.true.))
    call handle_nc_err(nf90_put_var(ncid, varid, cellon_89))

    call handle_nc_err(nf90_def_var(ncid, "incidence_89", NF90_FLOAT, [dim_fov_89, dim_scan_89], varid, &
           deflate_level=2, shuffle=.true.))
    call handle_nc_err(nf90_put_var(ncid, varid, celtht_89))

    call handle_nc_err(nf90_def_var(ncid, "azimuth_89", NF90_FLOAT, [dim_fov_89, dim_scan_89], varid, &
           deflate_level=2, shuffle=.true.))
    call handle_nc_err(nf90_put_var(ncid, varid, celphi_89))

    call handle_nc_err(nf90_def_var(ncid, "range_89", NF90_FLOAT, [dim_fov_89, dim_scan_89], varid, &
           deflate_level=2, shuffle=.true.))
    call handle_nc_err(nf90_put_var(ncid, varid, celrng_89))

    call handle_nc_err(nf90_def_var(ncid, "native_tb", NF90_FLOAT, [dim_fov, dim_scan, dim_chan], varid, &
         deflate_level=2, shuffle=.true.))
    call handle_nc_err(nf90_put_var(ncid, varid, ta))

    call handle_nc_err(nf90_def_var(ncid, "native_tb_89", NF90_FLOAT, [dim_fov_89, dim_scan_89, dim_chan_89], varid, &
         deflate_level=2, shuffle=.true.))
    call handle_nc_err(nf90_put_var(ncid, varid, ta_89))

    ! Set attributes on datasets
    call handle_nc_err(nf90_inq_varid(ncid, "scan_index", varid))
    call handle_nc_err(nf90_put_att(ncid, varid, "long_name", "scan index"))
    call handle_nc_err(nf90_put_att(ncid, varid, "axis", "Y"))

    ! call handle_nc_err(nf90_inq_varid(ncid, "fov", varid))
    ! call handle_nc_err(nf90_put_att(ncid, varid, "long_name", "Earth-view sample index"))
    ! call handle_nc_err(nf90_put_att(ncid, varid, "axis", "X"))

    ! call handle_nc_err(nf90_inq_varid(ncid, "xyz", varid))
    ! call handle_nc_err(nf90_put_att(ncid, varid, "long_name", "x, y, and z"))

    ! call handle_nc_err(nf90_inq_varid(ncid, "band", varid))
    ! call handle_nc_err(nf90_put_att(ncid, varid, "long_name", "band frequency"))
    ! call handle_nc_err(nf90_put_att(ncid, varid, "units", "GHz"))

    ! call handle_nc_err(nf90_inq_varid(ncid, "scan_time_utc", varid))
    ! call handle_nc_err(nf90_put_att(ncid, varid, "standard_name", "time"))
    ! call handle_nc_err(nf90_put_att(ncid, varid, "units", "seconds since 1958-01-01 00:00:00"))
    ! call handle_nc_err(nf90_put_att(ncid, varid, "comment", "Includes leap seconds"))
    ! call handle_nc_err(nf90_put_att(ncid, varid, "coordinates", "sclat sclon"))

    ! call handle_nc_err(nf90_inq_varid(ncid, "leap_seconds", varid))
    ! call handle_nc_err(nf90_put_att(ncid, varid, "long_name", "accumulated leap seconds"))
    ! call handle_nc_err(nf90_put_att(ncid, varid, "units", "s"))
    ! call handle_nc_err(nf90_put_att(ncid, varid, "comment", "The sum of all leap seconds since the 1958 epoch"))
    ! call handle_nc_err(nf90_put_att(ncid, varid, "coordinates", "sclat sclon"))

    ! call handle_nc_err(nf90_inq_varid(ncid, "delta_ut1utc", varid))
    ! call handle_nc_err(nf90_put_att(ncid, varid, "long_name", "delta UT1-UTC"))
    ! call handle_nc_err(nf90_put_att(ncid, varid, "units", "s"))
    ! call handle_nc_err(nf90_put_att(ncid, varid, "coordinates", "sclat sclon"))

    ! call handle_nc_err(nf90_inq_varid(ncid, "elevation", varid))
    ! call handle_nc_err(nf90_put_att(ncid, varid, "standard_name", "height_above_reference_ellipsoid"))
    ! call handle_nc_err(nf90_put_att(ncid, varid, "units", "m"))

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

    call handle_nc_err(nf90_inq_varid(ncid, "incidence", varid))
    call handle_nc_err(nf90_put_att(ncid, varid, "long_name", "Earth incidence angle"))
    call handle_nc_err(nf90_put_att(ncid, varid, "standard_name", "sensor_view_angle"))
    call handle_nc_err(nf90_put_att(ncid, varid, "units", "degree"))

    call handle_nc_err(nf90_inq_varid(ncid, "scan_rate", varid))
    call handle_nc_err(nf90_put_att(ncid, varid, "long_name", "antenna revolutions per minute"))
    call handle_nc_err(nf90_put_att(ncid, varid, "units", "rpm"))

    call handle_nc_err(nf90_inq_varid(ncid, "scrpy", varid))
    call handle_nc_err(nf90_put_att(ncid, varid, "long_name", "spacecraft roll, pitch, and yaw"))
    call handle_nc_err(nf90_put_att(ncid, varid, "units", "degree"))
    call handle_nc_err(nf90_put_att(ncid, varid, "coordinates", "sclat sclon"))

    call handle_nc_err(nf90_inq_varid(ncid, "scvel", varid))
    call handle_nc_err(nf90_put_att(ncid, varid, "long_name", "spacecraft velocity"))
    call handle_nc_err(nf90_put_att(ncid, varid, "units", "m/s"))
    call handle_nc_err(nf90_put_att(ncid, varid, "comment", "ECI J2000 coordinates"))
    call handle_nc_err(nf90_put_att(ncid, varid, "coordinates", "sclat sclon"))

    ! call handle_nc_err(nf90_inq_varid(ncid, "sun_alpha", varid))
    ! call handle_nc_err(nf90_put_att(ncid, varid, "long_name", "scan angle toward Sun"))
    ! call handle_nc_err(nf90_put_att(ncid, varid, "units", "degree"))

    call handle_nc_err(nf90_inq_varid(ncid, "zang", varid))
    call handle_nc_err(nf90_put_att(ncid, varid, "long_name", "spacecraft orbit angle"))
    call handle_nc_err(nf90_put_att(ncid, varid, "units", "degree"))


    
    ! call handle_nc_err(nf90_inq_varid(ncid, "geoloc_qc", varid))
    ! call handle_nc_err(nf90_put_att(ncid, varid, "long_name", "Geolocation quality control flags"))
    ! call handle_nc_err(nf90_put_att(ncid, varid, "flag_masks", [integer(int8) :: 1, 2, 4, 8]))
    ! call handle_nc_err(nf90_put_att(ncid, varid, "flag_meanings", &
    !      "IsNotValid IsExclusionCondition IsDegradationCondition IsLimitedUtilityValidationCondition"))

    ! call handle_nc_err(nf90_inq_varid(ncid, "validation_condition", varid))
    ! call handle_nc_err(nf90_put_att(ncid, varid, "long_name", "Geolocation validation condition flags"))
    ! call handle_nc_err(nf90_put_att(ncid, varid, "flag_masks", [integer(int8) :: 1, 2]))
    ! call handle_nc_err(nf90_put_att(ncid, varid, "flag_meanings", "IsLand Coastline"))

    ! It seems like these attributes should be written, but doing so
    ! confuses Panoply. Not writing them results in "Geo2D" plottable
    ! variables, which is what is intended.
    ! call write_att_string(fid, "azimuth", "coordinates", "latitude longitude")
    ! call write_att_string(fid, "elevation", "coordinates", "latitude longitude")
    ! call write_att_string(fid, "faraday_rotation", "coordinates", "latitude longitude")
    ! call write_att_string(fid, "faraday_rotation_true", "coordinates", "latitude longitude")
    ! call write_att_string(fid, "geoloc_qc", "coordinates", "latitude longitude")
    ! call write_att_string(fid, "incidence", "coordinates", "latitude longitude")
    ! call write_att_string(fid, "latitude", "coordinates", "latitude longitude")
    ! call write_att_string(fid, "longitude", "coordinates", "latitude longitude")
    ! call write_att_string(fid, "ipp_mag", "coordinates", "latitude longitude")
    ! call write_att_string(fid, "pra", "coordinates", "latitude longitude")
    ! call write_att_string(fid, "sun_azimuth", "coordinates", "latitude longitude")
    ! call write_att_string(fid, "sun_glint", "coordinates", "latitude longitude")
    ! call write_att_string(fid, "sun_zenith", "coordinates", "latitude longitude")
    ! call write_att_string(fid, "tec", "coordinates", "latitude longitude")
    ! call write_att_string(fid, "tec_true", "coordinates", "latitude longitude")
    ! call write_att_string(fid, "land_fraction", "coordinates", "latitude longitude")

    call handle_nc_err(nf90_close(ncid))
  end subroutine write_tb_native_netcdf

    ! -------------------------------------------------------------------------------
  ! Create a new Tb file
  subroutine write_geoloc_extra_netcdf(filename_geo,orbit_number,overwrite)

    use l2_module, only: scan_time,cellat,cellon
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
    call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "git_revision", trim(RESAMPLE_VERSION)))
    call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "git_date", trim(RESAMPLE_VERSION_DATE)))
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
    !call handle_nc_err(nf90_def_dim(ncid, "band", nBand, dim_band))
    !call handle_nc_err(nf90_def_dim(ncid, "pol", nPol, dim_pol))
    call handle_nc_err(nf90_def_dim(ncid, "chan", nChannel, dim_chan))
    
    
    !Define and write coordinate variables
    !use satellite_constants, only: nChannel,nBand,nPol,nFOV,nScan
  
    call handle_nc_err(nf90_def_var(ncid, "scan", NF90_INT, dim_scan, varid))
    call handle_nc_err(nf90_put_var(ncid, varid, [(i, i = 1, maxscan_extra)]))

    call handle_nc_err(nf90_def_var(ncid, "fov", NF90_INT, dim_fov, varid))
    call handle_nc_err(nf90_put_var(ncid, varid, [(i, i = 1, maxcel_extra)]))

    call handle_nc_err(nf90_def_var(ncid, "chan", NF90_INT, dim_chan, varid))
    call handle_nc_err(nf90_put_var(ncid, varid, [(i, i = 1, nChannel)]))

    !call handle_nc_err(nf90_def_var(ncid, "xyz", NF90_INT, dim_xyz, varid))
    !call handle_nc_err(nf90_put_var(ncid, varid, [(i, i = 1, 3)]))
 
    
    print *,maxscan_extra
    print *,size(zang_extra)
    ! Define and write main variables

    call handle_nc_err(nf90_def_var(ncid, "zang", NF90_DOUBLE, dim_scan, varid, &
     deflate_level=2, shuffle=.true.))
    call handle_nc_err(nf90_put_var(ncid, varid, zang_extra))

    call handle_nc_err(nf90_def_var(ncid, "scan_time", NF90_DOUBLE, dim_scan, varid, &
         deflate_level=2, shuffle=.true.))
    call handle_nc_err(nf90_put_var(ncid, varid, scan_time_extra))

    call handle_nc_err(nf90_def_var(ncid, "scan_index", NF90_DOUBLE, dim_scan, varid, &
         deflate_level=2, shuffle=.true.))
    call handle_nc_err(nf90_put_var(ncid, varid, scan_index_extra))

    call handle_nc_err(nf90_def_var(ncid, "scan_rate", NF90_FLOAT, dim_scan, varid, &
         deflate_level=2, shuffle=.true.))
    call handle_nc_err(nf90_put_var(ncid, varid, omega_extra))


    call handle_nc_err(nf90_def_var(ncid, "latitude", NF90_FLOAT, [dim_fov, dim_scan], varid, &
         deflate_level=2, shuffle=.true.))
    call handle_nc_err(nf90_put_var(ncid, varid, cellat_extra))

    call handle_nc_err(nf90_def_var(ncid, "longitude", NF90_FLOAT, [dim_fov, dim_scan], varid, &
         deflate_level=2, shuffle=.true.))
    call handle_nc_err(nf90_put_var(ncid, varid, cellon_extra))

    call handle_nc_err(nf90_def_var(ncid, "range", NF90_FLOAT, [dim_fov, dim_scan], varid, &
    deflate_level=2, shuffle=.true.))
    call handle_nc_err(nf90_put_var(ncid, varid, celrange_extra))

    call handle_nc_err(nf90_def_var(ncid, "incidence", NF90_FLOAT, [dim_fov, dim_scan], varid, &
         deflate_level=2, shuffle=.true.))
    call handle_nc_err(nf90_put_var(ncid, varid, celtht_extra))

    call handle_nc_err(nf90_def_var(ncid, "azimuth", NF90_FLOAT, [dim_fov, dim_scan], varid, &
         deflate_level=2, shuffle=.true.))
    call handle_nc_err(nf90_put_var(ncid, varid, celphi_extra))

    call handle_nc_err(nf90_def_var(ncid, "tb_resampled", NF90_FLOAT, [dim_fov, dim_scan, dim_chan], varid, &
         deflate_level=2, shuffle=.true.))
    call handle_nc_err(nf90_put_var(ncid, varid, ta_extra))


    ! Set attributes on datasets
    call handle_nc_err(nf90_inq_varid(ncid, "scan_index", varid))
    call handle_nc_err(nf90_put_att(ncid, varid, "long_name", "scan index"))
    call handle_nc_err(nf90_put_att(ncid, varid, "axis", "Y"))

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

    call handle_nc_err(nf90_inq_varid(ncid, "incidence", varid))
    call handle_nc_err(nf90_put_att(ncid, varid, "long_name", "Earth incidence angle"))
    call handle_nc_err(nf90_put_att(ncid, varid, "standard_name", "sensor_view_angle"))
    call handle_nc_err(nf90_put_att(ncid, varid, "units", "degree"))

    call handle_nc_err(nf90_inq_varid(ncid, "scan_rate", varid))
    call handle_nc_err(nf90_put_att(ncid, varid, "long_name", "antenna revolutions per minute"))
    call handle_nc_err(nf90_put_att(ncid, varid, "units", "rpm"))

    call handle_nc_err(nf90_inq_varid(ncid, "zang", varid))
    call handle_nc_err(nf90_put_att(ncid, varid, "long_name", "spacecraft orbit angle"))
    call handle_nc_err(nf90_put_att(ncid, varid, "units", "degree"))

    call handle_nc_err(nf90_close(ncid))
  end subroutine write_geoloc_extra_netcdf

end module L2_Tb_io
