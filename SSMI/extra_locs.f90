module extra_locs

    public :: ssmi_geolocation_extra_lf

contains

    subroutine ssmi_geolocation_extra_lf(min_valid_extra_scan,max_valid_extra_scan)

        use,intrinsic :: iso_fortran_env, only: int32,real32,real64
        use,intrinsic :: ieee_arithmetic, only: ieee_is_finite
        
        use NAN_support, only: nan_f32, nan_f64
        use sat_id_module
        use ssmi_l1c,only: numscan,maxcel,maxscan
        use ssmi_l1c,only: cellat
        use ssmi_l1c,only: cellon
        use ssmi_l1c,only: celphi
        use ssmi_l1c,only: iqual_flag
        use ssmi_l1c,only: scan_time
        
        !use ssmi_l1c,only: numscan,maxcel,maxscan, &
        !                   cellat,cellon,celphi
        use l2_module  !everything!
        use l2_module_extra,only: num_extra_scans,num_extra_cels, &
                                   maxcel_extra,maxscan_extra, &
                                   cellat_extra,cellon_extra, &
                                   celphi_extra,scan_time_extra, &
                                   iscan_flag_extra
        use local_proj,only: lat_lon_to_local_xy,local_xy_to_lat_lon
        use local_proj,only: lat_lon_to_local_xy_real32,local_xy_to_lat_lon_real32
        use interpolate,only: wt_mean_angle

        implicit none

        integer(int32),intent(out) :: min_valid_extra_scan
        integer(int32),intent(out) :: max_valid_extra_scan

        logical :: min_valid_extra_scan_found
        real(real32),dimension(maxcel,maxscan) :: cellat2,cellon2,celphi2
        real(real64),dimension(maxscan) :: scan_time2

        real(real32),dimension(maxcel,2) :: x2,y2,phi2
        real(real64),dimension(2) :: time2

        real(real32),dimension(maxcel)   :: est_error_x2,est_error_y2

        real(real32),dimension(maxcel_extra,2+num_extra_scans) :: x_extra,y_extra,azimuth_extra
        real(real64),dimension(2+num_extra_scans) :: time_extra
        
        integer(int32) :: icel,iscan
        integer(int32) :: icel_extra,iscan_extra
        integer(int32) :: iscan1,icel2,iscan2,iscan_out
        
        real(real32) :: lat0,lon0
        real(real32) :: wt_upper,wt_lower,error_factor

        real(real32),dimension(2) :: azimuth_arr
        real(real32),dimension(2) :: weight_arr

        

        cellat_extra = -999.0
        cellon_extra = -999.0
        celphi_extra = -999.0
        scan_time_extra = -999999999.0
        iscan_flag_extra = 1
        min_valid_extra_scan = 0
        max_valid_extra_scan = 0
        min_valid_extra_scan_found = .false.

        ! transfer geolocation with the native footprints to
        ! local array


        do iscan = 1,numscan/2  !numscan is from module ssmi_l1c
            iscan2 = (iscan-1)*(2)+1
            if (iqual_flag(iscan2) .eq. 0) then
                scan_time2(iscan) = scan_time(iscan2)
            else
                scan_time2(iscan) = nan_f64
            endif
            do icel = 1,maxcel  !numcel is from module ssmi_l1c
                !odd scans and cels are the native footprints
                
                icel2 = (icel-1)*(2)+1

                if ((abs(cellat(icel2,iscan2)) < 0.001) .and. &
                    (abs(cellon(icel2,iscan2)) < 0.001)) then
                    ! print *,'bad location',iscan,icel,iscan2,icel2,cellat(icel2,iscan2),cellon(icel2,iscan2)
                    cellat2(icel,iscan) = nan_f32
                    cellon2(icel,iscan) = nan_f32
                    celphi2(icel,iscan) = nan_f32
                else
                    cellat2(icel,iscan) = cellat(icel2,iscan2)
                    cellon2(icel,iscan) = cellon(icel2,iscan2)
                    celphi2(icel,iscan) = celphi(icel2,iscan2)
                endif 
            enddo
            ! print *,iscan,celphi2(maxcel,iscan)
        enddo

        ! convert to x,y in a local projection
        

        do iscan = 1,numscan/2-1  !numscan is from module ssmi_l1c
            lat0 = cellat2(32,iscan)
            lon0 = cellon2(32,iscan)
            time2(1) = scan_time2(iscan)
            time2(2) = scan_time2(iscan+1)
            do icel = 1,maxcel  !numcel is from module ssmi_l1c
                call lat_lon_to_local_xy(cellat2(icel,iscan),  &
                                         cellon2(icel,iscan),  &
                                         lat0,                 &
                                         lon0,                 &
                                         x2(icel,1),            &
                                         y2(icel,1))
                phi2(icel,1) = celphi2(icel,iscan)
                
                ! print *,lat0,lon0
                ! print *,cellat2(icel,iscan),cellon2(icel,iscan),x2(icel,1),y2(icel,1)
                ! print *,icel,iscan
                ! print *
                call lat_lon_to_local_xy(cellat2(icel,iscan+1),  &
                                         cellon2(icel,iscan+1),  &
                                         lat0,                 &
                                         lon0,                 &
                                         x2(icel,2),            &
                                         y2(icel,2))
                phi2(icel,2) = celphi2(icel,iscan+1)

            enddo

            do icel = 2,maxcel-1 
                est_error_x2(icel) = 0.5*(x2(icel+1,1) + x2(icel-1,1)) - x2(icel,1)
                est_error_y2(icel) = 0.5*(y2(icel+1,1) + y2(icel-1,1)) - y2(icel,1)
            enddo

            est_error_x2(1) = est_error_x2(2)-(est_error_x2(3)-est_error_x2(2))
            est_error_y2(1) = est_error_y2(2)-(est_error_y2(3)-est_error_y2(2))

            est_error_x2(maxcel) = est_error_x2(maxcel-1)-(est_error_x2(maxcel-2)-est_error_x2(maxcel-1))
            est_error_y2(maxcel) = est_error_y2(maxcel-1)-(est_error_y2(maxcel-2)-est_error_y2(maxcel-1))
            
        !     ! transfer the native location the output array in x/y space for the native scans
        !     ! and interpolated the extra fovs in the scans
            do iscan1 =1,2
                iscan2 = (iscan1-1)*(num_extra_scans+1)+1
                time_extra(iscan2) = time2(iscan1)
                do icel_extra = 1,maxcel_extra
                    if (modulo(icel_extra,1+num_extra_cels) == 1) then
                        icel = 1+(icel_extra-1)/(num_extra_cels + 1)
                        x_extra(icel_extra,iscan2) = x2(icel,iscan1)
                        y_extra(icel_extra,iscan2) = y2(icel,iscan1)
                        azimuth_extra(icel_extra,iscan2) = phi2(icel,iscan1)
                        ! print *,icel_extra,iscan2,x_extra(icel_extra,iscan2),y_extra(icel_extra,iscan2),azimuth_extra(icel_extra,iscan2)
                    else
                        icel = 1 + ((icel_extra-1)/(num_extra_cels + 1))
                        wt_upper = (icel_extra - ((icel-1.0)*(num_extra_cels + 1.0) + 1))/ &
                                   (1+num_extra_cels)
                        wt_lower = 1.0-(wt_upper)
                        error_factor = wt_upper*wt_lower
                        ! print *,wt_lower,wt_upper,error_factor
                        x_extra(icel_extra,iscan2) = wt_lower*x2(icel,iscan1) + &
                                                     wt_upper*x2(icel+1,iscan1) + &
                                                     error_factor*est_error_x2(icel)

                        y_extra(icel_extra,iscan2) = wt_lower*y2(icel,iscan1) + &
                                                     wt_upper*y2(icel+1,iscan1) + &
                                                     error_factor*est_error_y2(icel)

                        azimuth_arr(1) = phi2(icel,iscan1)
                        azimuth_arr(2) = phi2(icel+1,iscan1)
                        ! print *,azimuth_arr(1),azimuth_arr(2)
                        weight_arr(1) = wt_lower
                        weight_arr(2) = wt_upper
                        azimuth_extra(icel_extra,iscan2) = wt_mean_angle(azimuth_arr,weight_arr)
                        ! if (icel_extra == 2) then
                        !      print *,azimuth_arr(1),azimuth_arr(2)
                        !      print *,icel_extra,iscan2,azimuth_extra(icel_extra,iscan2)
                        ! endif
                        
                    endif
                enddo
            enddo

        !     ! fill in the extra scans

            do iscan2 = 2,num_extra_scans+1
                wt_upper = (iscan2-1.0)/(num_extra_scans+1)
                wt_lower = 1.0-wt_upper
                ! print *,wt_lower,wt_upper,iscan2
                time_extra(iscan2) = wt_lower*time_extra(1) + &
                                     wt_upper*time_extra(2+num_extra_scans)
                    
                ! print *,time_extra(1),time_extra(2+num_extra_scans)
                ! print *,wt_lower,wt_upper
                ! print *, time_extra(iscan2)
                ! print *
                do icel_extra = 1,maxcel_extra
                    x_extra(icel_extra,iscan2) = wt_lower*x_extra(icel_extra,1) + &
                                                 wt_upper*x_extra(icel_extra,2+num_extra_scans)

                    y_extra(icel_extra,iscan2)       = wt_lower*y_extra(icel_extra,1) + &
                                    wt_upper*y_extra(icel_extra,2+num_extra_scans)
                    
                    azimuth_arr(1) = azimuth_extra(icel_extra,1)
                    azimuth_arr(2) = azimuth_extra(icel_extra,2+num_extra_scans)  
                    weight_arr(1) = wt_lower
                    weight_arr(2) = wt_upper

                    azimuth_extra(icel_extra,iscan2) = wt_mean_angle(azimuth_arr,weight_arr)
                enddo
            enddo

            ! convert back to lat/lon
            do iscan2 = 1,num_extra_scans+1
                iscan_out = (iscan-1)*(num_extra_scans+1) + iscan2
                scan_time_extra(iscan_out) = time_extra(iscan2)
                if (ieee_is_finite(scan_time_extra(iscan_out))) then
                    iscan_flag_extra(iscan_out) = 0
                    if (.not. min_valid_extra_scan_found) then
                        min_valid_extra_scan = iscan_out
                        min_valid_extra_scan_found = .true.
                    endif
                    max_valid_extra_scan = iscan_out
                endif
                ! print *,iscan,iscan2,iscan_out
                do icel_extra = 1,maxcel_extra
                    call local_xy_to_lat_lon(x_extra(icel_extra,iscan2),        &
                                             y_extra(icel_extra,iscan2),        &
                                             lat0,                              &
                                             lon0,                              &
                                             cellat_extra(icel_extra,iscan_out),   &
                                             cellon_extra(icel_extra,iscan_out))
                    ! tranfer the interpolated azimuth angle
                    celphi_extra(icel_extra,iscan_out) = azimuth_extra(icel_extra,iscan2)
                enddo
            enddo
        enddo
    end subroutine ssmi_geolocation_extra_lf
end module extra_locs