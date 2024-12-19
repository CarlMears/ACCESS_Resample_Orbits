module amsr2_access_entry
    use, intrinsic :: iso_fortran_env, only: int8, int32, real32, real64 
    use, intrinsic :: iso_c_binding, only: c_int, c_bool
    implicit none
    private

    !public :: process_and_resample_amsr2_l2b_orbit

contains
    
    subroutine process_and_resample_amsr2_l2b_orbit(iorbit,  &                  ! orbit to process
                                                    channel_list,  &            ! list of AMSR2 channels to process
                                                    footprint_size, &           ! footprint size to resample to
                                                    write_swath_tbs,  &         ! set true to write swath tb file
                                                    write_swath_geoloc, &       ! set true to write swath tb file with extra locations
                                                    overwrite_existing_files, & ! if true, existing files are overwritten
                                                    pathname_tb_output)         ! path to output files
        use RSS_L1A, only:  read_rss_l1a,read_binary_head, &
                            level1a_rss_header, level1a_rss, & ! types
                            L1arss_hd,L1arss                ! module variables

        use sat_id_module, only: amsr2_nfreq => nfreq
        use l2_module        
        use l2_tb_module_ACCESS
        use Geolocation, only: amsr_geolocation,amsr2_geolocation_extra
        !use filename_routines
        use percent_land, only: fd_percent_land
        use findice, only: fdice
        use filename_routines, only: get_l2b_filename_base
        use FileExist, only: file_exist
        use open_file_routines, only: open_binary
        use date_utilities, only: find_month_day,fd_time_2000,fd_date_2000,days_in_month,doy_from_year_month_day
        use netcdf
        use DataFiles, only: readin_data_files
        use L2_Tb_io, only: write_tb_native_netcdf, write_geoloc_extra_netcdf
        use l2_module_extra
        use resample, only: init_resample_weights,ResampleWeights,resample_ta_extra,free_resample_weights
        use resample, only: find_closest_resample_point
        use resample, only: resample_onto_lat_lon_grid,resample_onto_lat_lon_grid_interp
        use resample, only: resample_time_onto_lat_lon_grid,write_resampled_time_map
        use resample, only: write_resampled_maps,write_wt_slice_netcdf
        use NAN_support, only: nan_f32, nan_f64
        use equirectangular_maps, only: map_str,write_map_netcdf
        !use logger, only: log_msg_string, initialize_logger, close_log
        !use logger, only: set_filter, set_screen_filter, set_file_filter
        !use logger, only: log_line, log_debug, log_info, log_error, log_crit
                                                
        implicit none    

        integer(int32),intent(in)  :: iorbit
        integer(int32),dimension(20),intent(in) :: channel_list
        integer(int32),intent(in)  :: footprint_size
        logical(c_bool),intent(in) :: write_swath_tbs
        logical(c_bool),intent(in) :: write_swath_geoloc
        logical(c_bool),intent(in) :: overwrite_existing_files
        character(len=*),intent(in) :: pathname_tb_output

        logical,dimension(amsr2_nfreq) :: do_freq

        character(:), allocatable :: str                                               
        real(real64)     :: start_time
        real(real64)     :: orbit_time                     

        integer(int32)  :: iscan
        integer(int32)  :: icel
        integer(int32)  :: iexist
        integer(int32)  :: ierror1
        integer(int32)  :: ierror2
        integer(int32)  :: ibad
        integer(int32)  :: kbad

        logical  :: lexist
        logical  :: tbs_ok_this_scan


        real(real64)     :: secyr
        real(real64)     :: secdy                     
        integer(int32)  :: lyear
        integer(int32)  :: idayjl
        integer(int32)  :: imon
        integer(int32)  :: idaymo
        real(real32)     :: frcrev
        integer(int32)  :: ifreq
        integer(int32)  :: ich1
        integer(int32)  :: ich2
        integer(int32)  :: ichan

        integer(int32)  :: ix,iy
        integer(int32)  :: NP_error    
        integer(int32)  :: SP_error
        real(real32)     :: xlat0
        real(real32)     :: xlon0
        real(real32),dimension(16) :: tax
        real(real32) ta89_resp(maxcel,maxscan,13:16)
            
        character(256) :: acmd(100)
        integer(int32)     :: numarg

        character(150)  :: test_filename
        character(150)  :: filename_l2b_tb_base
        character(150)  :: filename_l2b_tb
        character(150)  :: filename_l2b_tb_nc
        character(150)  :: filename_l2b_tb_resamp
        character(150)  :: filename_rss_l1a
        character(180)  :: filename_l2b_tb_grid 
        character(180)  :: filename_l2b_tb_dist 


        integer(int32) :: lyearfx
        integer(int32) :: nleap
        integer(int32) :: imonfx
        integer(int32) :: idayfx
        integer(int32) :: fice

        real(real32),dimension(2) :: tb_mea
        real(real32),dimension(2) :: ta_mea
        real(real32),dimension(5) :: perc_land
            
        integer(int32) :: overwrite
        integer(int32) :: write_error
        integer(int32), parameter :: ilu_l2b_tb  =  55

        type(tb_scn) :: tb_scan
        integer(int32)   :: num_scans

        character(10) :: overwrite_command

        type(ResampleWeights) :: amsr2_wts
        type(map_str) :: gridded_Tb_map
        type(map_str) :: gridded_time_map
        logical       :: do_time_interpolation

        real(real32)  :: grid_lon,grid_lat,search_size
        integer(int32) :: scan_closest,fov_closest
        real(real32) :: distance
        integer(int32) :: channel
        common/resampled89/ ta89_resp

        str = nf90_strerror(0)


        ksat=34  
        call define_filenames

        if (overwrite_existing_files) then
            overwrite = 1
        else 
            overwrite = 0
        endif


        !call get_l2b_filename_base(pathname_tb_output,1,34,iorbit, filename_l2b_tb_base,iexist) !1 means final
        !print *,pathname_tb_output
        call get_l2b_filename_base(pathname_tb_output,iorbit, filename_l2b_tb_base)
        filename_l2b_tb = trim(filename_l2b_tb_base)//'.tb'
        filename_l2b_tb_nc = trim(filename_l2b_tb_base)//'.nc'
        filename_l2b_tb_resamp = trim(filename_l2b_tb_base)//'.geo.nc'

        call file_exist(filename_l2b_tb, lexist)
        if ((overwrite .eq. 0) .and. (lexist)) then
            print *,'Tb File Exists, skipping orbit ',iorbit
            return
        endif


        ! read in ancillary data and initilize allocatable arrays
        ice_filename='x:\land\ice4.dat'  !use this ice mask when you do dynamic ice detection,it includes iceberg alley 
        call readin_data_files

        scan_time = nan_f64
        call read_rss_l1a(iorbit,0, iexist,start_time)    !0 means standard read (not quick read) 
        if(iexist.eq.0) then
            print *,'L1A file does not exist, skipping orbit ',iorbit
            return
        else
            print *,'Loaded L1A file for ',iorbit
        endif
            
        call fd_date_2000(start_time, secyr,lyear,idayjl,imon,idaymo,secdy)
        !write(*,'(i3,i6,i5,2i3,f7.2)') ksat,iorbit,lyear,imon,idaymo,secdy/3600.
            
        if(numscan.lt.35) then
            print *,'not enough scans to do scan averaging and resampling'
            stop
        endif
            
        ! orbit seems worth processing.  Open the L2B brightness temperature file

        ! start L2 TB processing
        !print *,'Doing Geolocation'
        call amsr_geolocation  
        !print *,'Doing Reverse AGC'  
        call reverse_agc

        !print *,'Finding Cal Temps'  
        call fd_thot      !must be called first and only once becuase it initializes  iflag_cal
        call fd_tcold(0)  !0 means use climate ta map for cold mirror spillover correction 

        !print *,'Doing Cold Mirror Calcs'  
        call moon_in_cold_view
        call fd_cold_mirror_rfi(iorbit)
        !print *,'Averaging Cal Counts'  
        call avg_cal_counts
        !print *,'Adjusting TAs'  
        call fd_ta_adjusted(6) !do first 6 channels
        call mk_spill_tamap
        call fd_tcold(1)        !1 means use actual ta measurements for cold mirror spillover correction
        call fd_ta_adjusted(16) !do all channels

        !Above tis point, it is the same as regular AMSR2 processing
        !print *,'Convert to Tb'  
        !Not sure if we have to do the next line....
        !call convert_ta_to_tb  !does not do 89 ghz, this is done in call stats

        call convert_ta_to_tb_native_only

        if (write_swath_tbs) then
            ! print *,'Writing Native Tbs to ',  trim(filename_l2b_tb_nc)
            ! call write_tb_native_netcdf(trim(filename_l2b_tb_nc),iorbit,overwrite)
        endif


        !print *,'Calculating Locations of Extra FOVs'  
        call amsr2_geolocation_extra

        !determine which frequencies are needed
        do_freq = .false.
        do ichan = 1,size(channel_list)
            if (channel_list(ichan) .gt. 0) then
                !print *,(channel_list(ichan)+1)/2
                do_freq((channel_list(ichan)+1)/2) = .true.
            endif
        enddo

        !print *,do_freq

        do ifreq = 1,amsr2_nfreq
            if (.not. (do_freq(ifreq))) cycle
            print *,'Resampling to ',footprint_size,' km'
            print *,'Resampling Channel: ',(2*ifreq)-1,' orbit',iorbit
            call init_resample_weights('AMSR2',band_number(ifreq),footprint_size,amsr2_wts)
            !call write_wt_slice_netcdf(amsr2_wts,trim(freq_name(ifreq)),1,243)
           
            call resample_ta_extra(amsr2_wts,(2*ifreq)-1)  ! Vpol
            print *,'Resampling Channel: ',(2*ifreq),' orbit',iorbit
            call resample_ta_extra(amsr2_wts,2*ifreq)      ! Hpol
            call free_resample_weights(amsr2_wts)
        enddo

        if (write_swath_geoloc) then
            print *,'Writing Resampled Orbit'
            call write_geoloc_extra_netcdf(trim(filename_l2b_tb_resamp),iorbit,overwrite)
        endif

        ! just put this back in to make the time data
        !print *,'Resampling time onto lat/lon grid, orbit',iorbit
        !call resample_time_onto_lat_lon_grid(0.25,gridded_time_map,gridded_distance_map)

        do ichan = 1,size(channel_list)
            channel = channel_list(ichan)
            if (channel == 0) cycle
            print *,'Resampling onto lat/lon grid using interpolation, channel: ',channel
            if (ichan .eq. 1) then
                do_time_interpolation = .true.
                print *,'Resampling time onto lat/lon grid using interpolation'
            else
                do_time_interpolation = .false.
            endif

            !interpolate onto the fixed grid
            call resample_onto_lat_lon_grid_interp(0.25,channel,do_time_interpolation,gridded_time_map,gridded_Tb_map)

            write(filename_l2b_tb_grid,161) trim(filename_l2b_tb_base),channel,footprint_size
            161 format(a,'.grid_tb.ch',i2.2,'.',i3.3,'km.nc')
            print *, filename_l2b_tb_grid
            call write_resampled_maps(gridded_Tb_map,filename_l2b_tb_grid)

            if (ichan .eq. 1) then
                write(filename_l2b_tb_grid,261) trim(filename_l2b_tb_base)
                261 format(a,'.time.nc')
                print *,'Writing Resampled Time, orbit',iorbit
                print *, filename_l2b_tb_grid
                call write_resampled_time_map(gridded_time_map,filename_l2b_tb_grid)
            endif
        enddo
    end subroutine
end module amsr2_access_entry
