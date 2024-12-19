module resample_polar

    use, intrinsic :: iso_fortran_env, only: int8, int32, real32, real64, ERROR_UNIT
    use ieee_arithmetic, only: ieee_is_finite

    use map_type_constants, only:CELL_CENTERED_2D_MAP,EDGE_CENTERED_2D_MAP,SMALL_NUMBER
    use earth_ref_maps, only: earth_ref_map_str32,earth_ref_map_str64
    use ease2_calcs, only: ease2
    use polar_grids, only: init_polar_grid32,init_polar_grid64
    use NSIDC_grids, only: read_NSIDC_polargrid
    
    implicit none 

    integer(int32),parameter :: max_select=5000

    private
    public resample_onto_polar_map_interp
    !public write_polar_resampled_time_map
    !public write_polar_resampled_maps
    
contains


    subroutine resample_onto_polar_map_interp(pole,  &
                                              resolution, &
                                              channel, &
                                              do_time, &
                                              gridded_time, &
                                              gridded_tbs, &
                                              min_valid_extra_scan, &
                                              max_valid_extra_scan, &
                                              num_overlap_scans, &
                                              error)

        ! This version of the resampling routine interpolates to the desired grid point using a general 
        ! quadrilateral interpolation method.

        ! It is different from resample_onto_lat_lon_grid_interp because the lats and lons in the map
        ! structure do not need to be on a rectangular grid.

        ! The looping is a little different than the previous versions
        !
        ! Here, the outside loops are over scans/fovs
        ! Then, the center point of the quadrilateral formed by (iscan,ifov),(iscan,ifov+1),(iscan+1,ifov+1),(iscan+1,ifov)
        ! is found.
        ! Then, the closest target point to this point is found and checked to see if it inside the quadrilateral
        ! If it is, then the interpolation to this point is done.
        !
        ! The algorithm also searches nearby fixed grid points to see if they are inside the quadrilateral
        ! The search range in the E-W direction expands near the pole so all possible grid points should be
        ! filled when using our standard 0.25 degree grids.

        use l2_module, only: re  !radius of Earth, not regular expressions (!)
        use l2_module_extra      ! this is where all the data are
        use QuadrilateralUtilties, only: check_convex, &
                                         is_inside_quadrilateral, &
                                         order_vertices_quadrilateral

        use InterpBilinearQuadrilateral, only: init_InterpBilinearQuadrilateral,eval_InterpBilinearQuadrilateral
        use InterpBilinearQuadrilateral, only: quad_interpolator_type

        implicit none

        character(*),intent(in)     :: pole
        character(*),intent(in)     :: resolution
        integer(int32),intent(in)   :: channel
        logical,intent(in)          :: do_time

        type(earth_ref_map_str64),intent(inout) :: gridded_time
        type(earth_ref_map_str32),intent(inout) :: gridded_tbs

        integer(int32),intent(in)   :: min_valid_extra_scan
        integer(int32),intent(in)   :: max_valid_extra_scan
        integer(int32),intent(in)   :: num_overlap_scans

        integer(int32),intent(out) :: error

        !locals
        integer(int32)              :: ilat,ilon,closest_scan,closest_fov
        real(real32)                :: xlat,xlon,closest_distance
        integer(int32)              :: num_grid_points_filled
        integer(int32)              :: max_valid_scan

        integer(int32)              :: jscan,jfov,kscan,kfov,jlat,jlon,klat,klon,klat0,klon0,ilatx,ilonx,iquad,j,k,dk
        integer(int32)              :: ix0,iy0,iy,ix
        integer(int32)              :: ilon_closest,ilat_closest
        real(real32)                :: distance
        logical                     :: bad_loc,is_convex,is_inside
        integer(int32)              :: interp_error
        integer(int32)              :: i
        type(quad_interpolator_type):: qi

        real(real64)                :: tot_lat,tot_lon,first_lat,first_lon
        real(real32)                :: grid_lat,grid_lon,x,y
        real(real32),dimension(0:3) :: quad_lats,quad_lons,quad_ta,quad_x,quad_y
        real(real64),dimension(0:3) :: quad_times
        real(real64),dimension(0:3) :: quad_lats_sorted,quad_lons_sorted,quad_ta_sorted,quad_times_sorted
        real(real64),dimension(0:3) :: quad_x_sorted,quad_y_sorted
        
        real(real64)                :: mean_lat,mean_lon,tb_test,mean_x,mean_y
        integer(int32),dimension(0:3) :: sorted_indices
        
        character(20)               :: grid_type
        logical                     :: no_polar_data
        real(real32),parameter      :: polar_lat_limit = 30.0
        real(real32),parameter      :: fill_value32 = -999.0
        real(real64),parameter      :: fill_value64 = -999.0

        logical,parameter           :: verify_lat_lon = .true.
        real(real64)                :: lat_interp, lon_interp,x_verify, y_verify
        integer(int32)              :: num_locs_verified,num_locs_failed

        ! find number of valid scans
        max_valid_scan = 0
        do jscan = 1,maxscan_extra
            if (.not. ieee_is_finite(scan_time_extra(jscan))) cycle
            max_valid_scan = jscan
        enddo
        if (max_valid_scan == 0) then
            error = 1
            print *,'No Valid Scans'
            return
        endif

        no_polar_data = .true.
        scan_loop: do jscan = min_valid_extra_scan+num_overlap_scans,max_valid_extra_scan-num_overlap_scans
            fov_loop: do jfov = 1,maxcel_extra-1
                if (ieee_is_finite(cellat_extra(jfov,jscan))) then
                    if (pole == 'north') then
                        if (cellat_extra(jfov,jscan) > polar_lat_limit) then
                            no_polar_data = .false.
                            exit scan_loop
                        endif
                    endif
                    if (pole == 'south') then
                        if (cellat_extra(jfov,jscan) < -polar_lat_limit) then
                            no_polar_data = .false.
                            exit scan_loop
                        endif
                    endif
                endif
            enddo fov_loop
        enddo scan_loop

        if (no_polar_data) then
            error = 1
            print *,'No Scans with Polar Data for ',pole
            return
        endif

        grid_type = 'ease2'
        !print *,pole,resolution,grid_type
        call init_polar_grid32(pole,        &
                               resolution,  &
                               grid_type,   &
                               gridded_tbs, &
                               error)

        if (do_time) then
            call init_polar_grid64(pole,        &
                               resolution,  &
                               grid_type,   &
                               gridded_time, &
                               error)
        endif

        call gridded_tbs%set_to_constant(fill_value32,error)
        call gridded_time%set_to_constant(fill_value64,error)
        num_locs_verified = 0
        num_locs_failed = 0
        ! loop over resampled observations
        do jscan = min_valid_extra_scan+num_overlap_scans,max_valid_extra_scan-num_overlap_scans  
            do jfov = 1,maxcel_extra-1

                bad_loc = .false.
                do kscan=jscan,jscan+1
                    do kfov=jfov,jfov+1
                        if (.not. ieee_is_finite(cellat_extra(kfov,kscan))) bad_loc = .true.
                        if (.not. ieee_is_finite(cellon_extra(kfov,kscan))) bad_loc = .true.
                        if (.not. ieee_is_finite(tbr_extra(kfov,kscan,channel))) bad_loc = .true.
                        if (tbr_extra(kfov,kscan,channel) < 20.0) bad_loc = .true.
                    enddo
                enddo

                if (bad_loc) cycle  !either the lat/lon or Ta is bad for one of vertices of the quadrilateral

                !find lat/lon of center of the quadrilateral
                !and store vertices
                tot_lat = 0.0
                tot_lon = 0.0
                iquad = 0
                first_lon = cellon_extra(jfov,jscan)
                do kscan=jscan,jscan+1
                    do kfov=jfov,jfov+1
                        tot_lat = tot_lat + cellat_extra(kfov,kscan)
                        quad_lats(iquad) = cellat_extra(kfov,kscan)

                        if (abs(first_lon - cellon_extra(kfov,kscan)) .lt. 180.0) then
                            tot_lon = tot_lon + cellon_extra(kfov,kscan)
                            quad_lons(iquad) = cellon_extra(kfov,kscan)
                        else if (first_lon - cellon_extra(kfov,kscan) .gt. 180.0) then
                            ! this is the case where first_lon is near 360.0 and other lons are near 0.0
                            tot_lon = tot_lon + 360.0 + cellon_extra(kfov,kscan)
                            quad_lons(iquad) = cellon_extra(kfov,kscan)+360.0
                        else
                            ! this is the case where first lon is near 0.0, and cellon_extra(kfov,kscan) > 360.0
                            tot_lon = tot_lon - 360.0 + cellon_extra(kfov,kscan)
                            quad_lons(iquad) = cellon_extra(kfov,kscan)-360.0
                        endif
                        quad_ta(iquad) = tbr_extra(kfov,kscan,channel)
                        if (do_time) quad_times(iquad) = scan_time_extra(kscan)
                                                

                        ! these splices were for verification purposes
                        ! the results are in
                        !
                        !L:\access\amsr2_tb_orbits\r00001_05000\r00676.grid_interp_lat.05.nc
                        !L:\access\amsr2_tb_orbits\r00001_05000\r00676.grid_interp_lon.05.nc
                        !
                        ! verification case splice latitude
                        !quad_ta(iquad) = cellat_extra(kfov,kscan)
                        ! end verification case splice

                        ! verification case splice latitude
                        !quad_ta(iquad) = cellon_extra(kfov,kscan)
                        ! end verification case splice

                        iquad=iquad+1
                    enddo
                enddo

                if (pole == 'north') then
                    if (maxval(quad_lats) < 30.0) cycle
                endif

                if (pole == 'south') then
                    if (minval(quad_lats) > -30.0) cycle
                endif
                
                ! print *,quad_lats
                ! print *,quad_lons
                

                !this is the center of the quadrilateral
                mean_lat = tot_lat/4.0
                mean_lon = modulo((tot_lon/4.0) + 1080.0,360.0)

                ! print *,mean_lat
                ! print *,mean_lon

                !convert the lat/lons of vertices to x,y in the polar projection
                do iquad = 0,3
                    call ease2(quad_lats(iquad),quad_lons(iquad), quad_x(iquad), quad_y(iquad))
                enddo

                ! print *,quad_x
                ! print *,quad_y


                !sort the quadrilateral points and the associated Ta's
                call order_vertices_quadrilateral(quad_x,quad_y,sorted_indices)
                do iquad=0,3
                    quad_x_sorted(iquad) = quad_x(sorted_indices(iquad))
                    quad_y_sorted(iquad) = quad_y(sorted_indices(iquad))
                    quad_ta_sorted(iquad)   = quad_ta(sorted_indices(iquad))
                    if (do_time) quad_times_sorted(iquad) = quad_times(sorted_indices(iquad))
                    if (verify_lat_lon) then
                        quad_lats_sorted(iquad) = quad_lats(sorted_indices(iquad))
                        quad_lons_sorted(iquad) = quad_lons(sorted_indices(iquad))
                    endif
                enddo

                is_convex = check_convex(quad_x_sorted,quad_y_sorted)
                if (.not. is_convex) then
                    print *,quad_x_sorted
                    print *,quad_y_sorted
                    print *,'Quadrilateral Not Convex - Should NOT happen -- skipping'
                    cycle
                endif

                !find closest point in the target polar grid 
                call ease2(mean_lat,mean_lon,mean_x,mean_y)
                call gridded_tbs%get_grid_index(real(mean_x,real32),real(mean_y,real32),ix0,iy0,error)
                ! print *,mean_lat,mean_lon
                ! print *,gridded_tbs%lats(ix0,iy0),gridded_tbs%lons(ix0,iy0)
                ! print *,ix0,iy0

                ! loop over the area to find target locations inside the quadrilateral
                dk = 1
                ! ! search a larger longitude range close to the poles for 
                ! ! grid points inside the quadrilateral
                ! if (abs(mean_lat) .gt. 84.0) dk = 3
                ! if (abs(mean_lat) .gt. 87.0) dk = 7
                ! if (abs(mean_lat) .gt. 88.0) dk = 15
                do j =-1,1
                    do k = -1,1
                        ix = ix0+j 
                        iy = iy0+k 

                        call gridded_tbs%get_grid_location(ix,iy,x,y,error)

                        if (error == 0) then
                            ! print *, j,k,ix,iy
                            ! print *, quad_x_sorted
                            ! print *, x 
                            ! print *, quad_y_sorted
                            ! print *, y 

                            !check to make sure that x(ix,iy),y(ix,iy) is inside the quadrilateral
                            is_inside = is_inside_quadrilateral(quad_x_sorted,       &
                                                                quad_y_sorted,       &
                                                                real(x,real64),      &
                                                                real(y,real64),.false.) ! the false is for return_reordered

                            if (is_inside) then
                                !do the interpolation
                                call init_InterpBilinearQuadrilateral(     &
                                                quad_x_sorted,   & ! array of 4 x locations
                                                quad_y_sorted,   & ! array of 4 y locations
                                                quad_ta_sorted,     & ! array of 4 z locations
                                                qi,                 & ! interpolator instance
                                                interp_error)         ! error flag    
                                if (interp_error == 0) then 
                                    gridded_tbs%dat(ix,iy) = eval_InterpBilinearQuadrilateral(real(x,real64),real(y,real64),qi,error)
                                endif

                                !do the time interpolation
                                if (do_time) then
                                    call init_InterpBilinearQuadrilateral(     &
                                                    quad_x_sorted,   & ! array of 4 x locations
                                                    quad_y_sorted,   & ! array of 4 y locations
                                                    quad_times_sorted,  & ! array of 4 z locations
                                                    qi,                 & ! interpolator instance
                                                    interp_error)         ! error flag    
                                    if (interp_error == 0) then 
                                        gridded_time%dat(ix,iy) = eval_InterpBilinearQuadrilateral(real(x,real64),real(y,real64),qi,error)
                                    endif
                                endif

                                if (verify_lat_lon) then
                                    lat_interp = -999.0
                                    call init_InterpBilinearQuadrilateral(     &
                                                    quad_x_sorted,   & ! array of 4 x locations
                                                    quad_y_sorted,   & ! array of 4 y locations
                                                    quad_lats_sorted,  & ! array of 4 z locations
                                                    qi,                 & ! interpolator instance
                                                    interp_error)         ! error flag
                                    if (interp_error == 0) then  
                                        lat_interp = eval_InterpBilinearQuadrilateral(real(x,real64),real(y,real64),qi,error)
                                    endif

                                    lon_interp = -999.0
                                    call init_InterpBilinearQuadrilateral(     &
                                                    quad_x_sorted,   & ! array of 4 x locations
                                                    quad_y_sorted,   & ! array of 4 y locations
                                                    quad_lons_sorted,  & ! array of 4 z locations
                                                    qi,                 & ! interpolator instance
                                                    interp_error)         ! error flag
                                    if (interp_error == 0) then  
                                        lon_interp = eval_InterpBilinearQuadrilateral(real(x,real64),real(y,real64),qi,error)
                                    endif
                                    call ease2(lat_interp,lon_interp, x_verify, y_verify)
                                    num_locs_verified = num_locs_verified + 1
                                    !looks for locations errors larger than 0.05 km, or 50 m in x or y
                                    if ((abs(x-x_verify) > 0.05) .or. (abs(y-y_verify) > 0.05)) then
                                        num_locs_failed = num_locs_failed+1
                                        print *,x,x_verify,x-x_verify
                                        print *,y,y_verify,y-y_verify
                                        print *,gridded_tbs%lats(ix,iy),lat_interp
                                        print *,gridded_tbs%lons(ix,iy),lon_interp
                                        print *,num_locs_verified,num_locs_failed
                                        print *
                                    endif
                                endif
                            endif
                        endif
                    enddo
                enddo
            enddo !jfov
        enddo !jscan
        if (verify_lat_lon) then
            print *,'Number Locations Verified: ',num_locs_verified
            print *,'Number locations Failed: ',num_locs_failed
            print *
        endif
    end subroutine resample_onto_polar_map_interp
end module resample_polar