    module NSIDC_grids

        private
        public :: read_NSIDC_polargrid

    contains
    
        subroutine read_NSIDC_polargrid(pole,resolution,grid_type,num_x,num_y,lat_array,lon_array,error)

            use, intrinsic :: iso_fortran_env, only: real32, real64, int64, int32
            use netcdf, only: nf90_inq_varid, nf90_get_var, nf90_close, nf90_open, nf90_nowrite !nf90_noerr, 
            !use rssnetcdf, only: nf_return_err_msg_no_stop, nf_handle_err
            use io_nc, only: handle_nc_err
            implicit none

            character(*), intent(in)  :: pole
            character(*), intent(in)  :: resolution
            character(*), intent(in)  :: grid_type
            integer(int32),intent(out) :: num_x
            integer(int32),intent(out) :: num_y

            real(real64),intent(inout),allocatable :: lat_array(:,:)
            real(real64),intent(inout),allocatable :: lon_array(:,:)
            integer(int32),intent(out) :: error

            character(100) :: polar_grid_root
            character(150) :: file
            logical :: match_found
            integer(int32)                      :: status
            integer(int32)                      :: ivid
            integer(int32)                      :: ncid

            match_found = .false.
            !print *,'grid_type = ',grid_type
            !print *,'pole      = ',pole
            !print *,'resol.    = ',resolution

#ifdef _WIN32
    polar_grid_root = 'M:\job_access\polar_grids\'
#else
    polar_grid_root = '/mnt/ops1p-ren/m/job_access/polar_grids/'
#endif
            
            if (trim(grid_type) .eq. 'ease2') then 
                if (trim(pole) .eq. 'north') then
                    if (trim(resolution) .eq. "25km") then
                        file = trim(polar_grid_root) // 'NSIDC0772_LatLon_EASE2_N25km_v1.0.nc'
                        num_x = 720
                        num_y = 720
                        match_found = .true.
                    endif

                    if (trim(resolution) .eq. "12.5km") then
                        file = trim(polar_grid_root) // 'NSIDC0772_LatLon_EASE2_N12.5km_v1.0.nc'
                        num_x = 1440
                        num_y = 1440
                        match_found = .true.
                    endif
                endif
                if (trim(pole) == "south") then
                    if (trim(resolution) .eq. "25km") then
                        file = trim(polar_grid_root) // 'NSIDC0772_LatLon_EASE2_S25km_v1.0.nc'
                        num_x = 720
                        num_y = 720
                        match_found = .true.
                    endif

                    if (trim(resolution) .eq. "12.5km") then
                        file = trim(polar_grid_root) // 'NSIDC0772_LatLon_EASE2_S12.5km_v1.0.nc'
                        num_x = 1440
                        num_y = 1440
                        match_found = .true.
                    endif
                endif
            endif

            if (match_found) then 

                allocate(lat_array(num_x,num_y))
                allocate(lon_array(num_x,num_y))

                status = nf90_open(file, NF90_NOWRITE, ncid)
                status = nf90_inq_varid(ncid,"latitude", ivid)
                status = nf90_get_var(ncid, ivid, lat_array)
                status = nf90_inq_varid(ncid,"longitude", ivid)
                status = nf90_get_var(ncid, ivid, lon_array)
                
                error = 0
            else
                error = -1
            endif

        end subroutine read_NSIDC_polargrid
    end module NSIDC_grids