 
 ! this program processes one amsr2 L1A orbit, and outputs various 
 ! files with brightness temperatures 
 !
 !
 ! example command
 !      
 ! start "17000 70km" cmd /K amsr2_exe_test_logic.exe 17000 17999 70 1 1 1 2 3 4 5 6 7 8 9 10 11 12                                  
 !                                                     strt  end  sz r p Channels
 !
 ! sz = footprint size
 ! r means do rectangular grid
 ! p means do polar grid
 ! channels: 1 = 6.9 GHz V
 !           2 = 6.9 GHz H
 !           3 = 7.2 GHz V, etc

 program process_amsr2_l2b_orbit_tb_ACCESS   

    !use, intrinsic :: iso_c_binding, only: c_bool
    use, intrinsic :: iso_fortran_env, only: int8, int32, real32, real64 
    !use amsr2_access_entry, only: process_and_resample_amsr2_l2b_orbit
    use QuadrilateralUtilties, only: triangle_area,is_inside_triangle

    implicit none 

    character(256)  :: acmd(100)
    integer(int32)  :: numarg
    integer(int32)  :: iorbit,iorbit_start,iorbit_end
    integer(int32)  :: footprint_size
    integer(int32)  :: nchannel
    integer(int32)  :: channel,i
    integer(int32)  :: use_rectangular_grid_int
    integer(int32)  :: use_NP_grid_int, use_SP_grid_int
    
    integer(int32),dimension(20) :: channel_list
    logical(4)         :: write_swath_tbs
    logical(4)         :: write_swath_geoloc
    logical(4)         :: overwrite_existing_files
    logical(4)         :: use_rectangular_grid  
    logical(4)         :: use_NP_grid, use_SP_grid
    integer(int32)  :: num_overlap_scans
    character(150)  :: pathname_tb_output

    ! read command line
    call get_cmd_arg(numarg,acmd)
    read(acmd(1),'(i6)') iorbit_start
    read(acmd(2),'(i6)') iorbit_end
    read(acmd(3),'(i3)') footprint_size
    read(acmd(4),'(i1)') use_rectangular_grid_int
    read(acmd(5),'(i1)') use_NP_grid_int
    read(acmd(6),'(i1)') use_SP_grid_int
    
    nchannel = numarg-6
    channel_list = 0
    if (nchannel .gt. 0) then
        do i = 1,numarg-6
            read(acmd(i+6),'(i6)') channel
            channel_list(i) = channel
        enddo
    endif

    use_rectangular_grid = .false.
    use_NP_grid = .false.
    use_SP_grid = .false.
    if (use_rectangular_grid_int > 0) use_rectangular_grid=.true.
    if (use_NP_grid_int > 0) use_NP_grid=.true.
    if (use_SP_grid_int > 0) use_SP_grid=.true.
    
    print *, iorbit_start
    print *, iorbit_end
    print *,channel_list

    if (use_rectangular_grid) print *,'resampling to rectangular grid'
    if (use_NP_grid) print *,'resampling to north polar grid'
    if (use_SP_grid) print *,'resampling to south polar grid'

    write_swath_tbs = .true.
    write_swath_geoloc = .true.
    overwrite_existing_files = .true.
    num_overlap_scans = 100
    
    pathname_tb_output = 'L:\access\amsr2_tb_orbits_test\  '
    do iorbit = iorbit_start,iorbit_end
        call process_and_resample_amsr2_l2b_orbit(  iorbit,             &
                                                channel_list,       &
                                                footprint_size,     &
                                                write_swath_tbs,    &
                                                write_swath_geoloc, &
                                                overwrite_existing_files, &
                                                use_rectangular_grid,  &    
                                                use_NP_grid,  &
                                                use_SP_grid,  &
                                                num_overlap_scans, &
                                                trim(pathname_tb_output))
    enddo
 end program process_amsr2_l2b_orbit_tb_ACCESS 

 subroutine process_and_resample_amsr2_l2b_orbit(iorbit,  &                  ! orbit to process
                                                 channel_list,  &            ! list of AMSR2 channels to process
                                                 footprint_size,  &
                                                 write_swath_tbs,  &         ! set true to write swath tb file
                                                 write_swath_geoloc, &       ! set true to write swath tb file with extra locations
                                                 overwrite_existing_files, & ! if true, existing files are overwritten
                                                 use_rectangular_grid,  &    ! set true to resample to rectangular grid
                                                 use_NP_grid,        &       ! set true to resample to north polar grid
                                                 use_SP_grid,        &       ! set true to resample to south polar grid
                                                 num_overlap_scans, &        ! number of scans to ignore at beginning and end of orbit
                                                 pathname_tb_output)         ! path to output files
                       

    use, intrinsic :: iso_fortran_env, only: int8, int32, real32, real64
    !use, intrinsic :: iso_c_binding, only: c_bool
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan,ieee_is_finite
  
    use RSS_L1A, only:  read_rss_l1a,read_binary_head, &
                        level1a_rss_header, level1a_rss, & ! types
                        L1arss_hd,L1arss                ! module variables

    use sat_id_module, only: amsr2_nfreq => nfreq,amsr2_channel_name => channel_name
    use l2_module        
    use l2_tb_module_ACCESS
    use Geolocation, only: amsr_geolocation,amsr2_geolocation_extra
    !use filename_routines
    use percent_land, only: fd_percent_land
    use findice, only: fdice
    use filename_routines, only: get_l2b_filename_base
    use FileExist, only: file_exist,file_exists  !file_exists is the funciton form
    use open_file_routines, only: open_binary
    use date_utilities, only: find_month_day,fd_time_2000,fd_date_2000,days_in_month,doy_from_year_month_day
    use netcdf
    use DataFiles, only: readin_data_files
    use L2_Tb_io, only: write_tb_native_netcdf, write_geoloc_extra_netcdf
    use l2_module_extra
    use resample, only: init_resample_weights,ResampleWeights,resample_ta_extra,free_resample_weights
    use resample, only: find_closest_resample_point
    use resample, only: resample_onto_lat_lon_grid_interp
    use resample, only: write_resampled_maps,write_wt_slice_netcdf
    use resample, only: resample_time_onto_lat_lon_grid,write_resampled_time_map
    use resample_polar, only: resample_onto_polar_map_interp
    use NAN_support, only: nan_f32, nan_f64
    use equirectangular_maps, only: map_str,write_map_netcdf
    use equirectangular_maps64, only: map_str64,write_map_netcdf64
    use earth_ref_maps, only: earth_ref_map_str32,earth_ref_map_str64
    use NSIDC_grids, only: read_NSIDC_polargrid
                                          
    implicit none    

    integer(int32),intent(in)  :: iorbit
    integer(int32),dimension(20),intent(in) :: channel_list
    integer(int32),intent(in) :: footprint_size
    logical(4),intent(in) :: write_swath_tbs
    logical(4),intent(in) :: write_swath_geoloc
    logical(4),intent(in) :: overwrite_existing_files
    logical(4),intent(in) :: use_rectangular_grid
    logical(4),intent(in) :: use_NP_grid
    logical(4),intent(in) :: use_SP_grid
    integer(int32) :: num_overlap_scans
    character(len=*) :: pathname_tb_output
    
    logical,dimension(amsr2_nfreq) :: do_freq_requested
    logical,dimension(amsr2_nfreq) :: do_freq_any
    logical,dimension(amsr2_nfreq) :: do_freq_rect
    logical,dimension(amsr2_nfreq) :: do_freq_polar
    
    
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
    type(map_str)   :: gridded_Tb_map
    type(map_str)   :: gridded_distance_map
    type(map_str64) :: gridded_time_map

    type(earth_ref_map_str32) :: polar_gridded_Tb_map
    type(earth_ref_map_str64) :: polar_gridded_time_map

    real(real32)  :: grid_lon,grid_lat,search_size
    integer(int32) :: scan_closest,fov_closest
    real(real32) :: distance
    integer(int32) :: channel
    common/resampled89/ ta89_resp

    logical :: do_time

    real(real64),allocatable :: lat_array(:,:)
    real(real64),allocatable :: lon_array(:,:)
    integer(int32)            :: num_x,num_y,error

    ksat=34  
    call define_filenames
    
    if (overwrite_existing_files) then
        overwrite = 1
    else 
        overwrite = 0
    endif

    call get_l2b_filename_base(pathname_tb_output,iorbit, filename_l2b_tb_base)
    
    filename_l2b_tb_nc = trim(filename_l2b_tb_base)//'.nc'
    filename_l2b_tb_resamp = trim(filename_l2b_tb_base)//'.geo.nc'
    
    ! The following logic figures out what needs to be done depending
    ! on the frequencies and output grids requested, and the existence
    ! of any pre-existing output files.
    !
    ! if do_freq_any(ichan) is true for a channel, then the tbs are required
    ! if do_freq_rect(ichan) is true, then the resampling to the rectangular grid is needed
    ! if do_freq_polar(ichan) is true, then the resampling to the polar grid is needed

    !determine which frequencies are requested
    do_freq_requested = .false.
    do ichan = 1,size(channel_list)
        if (channel_list(ichan) .gt. 0) then
            do_freq_requested((channel_list(ichan)+1)/2) = .true.
        endif
    enddo

    do_freq_any = .false.
    do_freq_rect = .false.
    do_freq_polar = .false.
    !print *,do_freq_requested
    !print *,channel_list
    do ichan = 1,size(channel_list)
        if (channel_list(ichan) .lt. 1) cycle
        if (do_freq_requested((channel_list(ichan)+1)/2)) then
            if (overwrite) then
                do_freq_rect((channel_list(ichan)+1)/2) = .true.
                do_freq_polar((channel_list(ichan)+1)/2) = .true.
                do_freq_any((channel_list(ichan)+1)/2) = .true.
            else
                channel = channel_list(ichan)
                write(filename_l2b_tb_grid,161) trim(filename_l2b_tb_base),channel,footprint_size
                !print *, filename_l2b_tb_grid
                if (.not. file_exists(filename_l2b_tb_grid)) then
                    do_freq_rect((channel_list(ichan)+1)/2) = .true.
                    do_freq_any((channel_list(ichan)+1)/2) = .true.
                endif

                write(filename_l2b_tb_grid,461) trim(filename_l2b_tb_base),channel,footprint_size
                !print *, filename_l2b_tb_grid
                if (.not. file_exists(filename_l2b_tb_grid)) then
                    do_freq_polar((channel_list(ichan)+1)/2) = .true.
                    do_freq_any((channel_list(ichan)+1)/2) = .true.
                endif
            endif
        endif
    enddo

    print *,do_freq_requested
    print *,do_freq_any
    print *,do_freq_rect
    print *,do_freq_polar

    if (.not. any(do_freq_any)) then
        print *,'All gridded files already exist, skipping orbit ',iorbit
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
    endif
        
    call fd_date_2000(start_time, secyr,lyear,idayjl,imon,idaymo,secdy)
    write(*,'(i3,i6,i5,2i3,f7.2)') ksat,iorbit,lyear,imon,idaymo,secdy/3600.
        
    if(numscan.lt.35) then
        print *,'not enough scans to do scan averaging and resampling'
        return
    endif
        
    ! orbit seems worth processing.  Open the L2B brightness temperature file

    ! start L2 TB processing
    print *,'Doing Geolocation'
    call amsr_geolocation  
    print *,'Doing Reverse AGC'  
    call reverse_agc

    print *,'Finding Cal Temps'  
    call fd_thot      !must be called first and only once becuase it initializes  iflag_cal
    call fd_tcold(0)  !0 means use climate ta map for cold mirror spillover correction 

    print *,'Doing Cold Mirror Calcs'  
    call moon_in_cold_view
    call fd_cold_mirror_rfi(iorbit)
    print *,'Averaging Cal Counts'  
    call avg_cal_counts
    print *,'Adjusting TAs'  
    call fd_ta_adjusted(6) !do first 6 channels
    call mk_spill_tamap
    call fd_tcold(1)        !1 means use actual ta measurements for cold mirror spillover correction
    call fd_ta_adjusted(16) !do all channels

    !Above tis point, it is the same as regular AMSR2 processing
    print *,'Convert to Tb'  
    !Not sure if we have to do the next line....
    !call convert_ta_to_tb  !does not do 89 ghz, this is done in call stats

    call convert_ta_to_tb_native_only

    if (write_swath_tbs) then
        print *,'Writing Native Tbs to ',  trim(filename_l2b_tb_nc)
        call write_tb_native_netcdf(trim(filename_l2b_tb_nc),iorbit,overwrite)
    endif

    ! this is the end of calculating the brightness temperatures at the 
    ! native swath locations

    print *,'Calculating Locations of Extra FOVs'  
    call amsr2_geolocation_extra

    print *,'Resampling Tbs to the Extra FOVs'  

    do ifreq = 1,amsr2_nfreq
        if (.not. (do_freq_any(ifreq))) cycle
        call init_resample_weights('AMSR2',band_number(ifreq),footprint_size,amsr2_wts)
        !call write_wt_slice_netcdf(amsr2_wts,trim(freq_name(ifreq)),1,243)

        print *,'Resampling Channel: ',(2*ifreq)-1
        call resample_ta_extra(amsr2_wts,(2*ifreq)-1)  ! Vpol
        print *,'Resampling Channel: ',(2*ifreq)
        call resample_ta_extra(amsr2_wts,2*ifreq)      ! Hpol
        call free_resample_weights(amsr2_wts)
    enddo

    ! write the resampled orbit -- this is done for debugging
    if (write_swath_geoloc) then
        print *,'Writing Resampled Orbit'
        call write_geoloc_extra_netcdf(trim(filename_l2b_tb_resamp),iorbit,overwrite)
    endif
    
    if ((use_rectangular_grid) .and. any(do_freq_rect)) then
        print *,'Resampling time onto rectangular lat/lon grid'
        call resample_time_onto_lat_lon_grid(0.25,gridded_time_map,gridded_distance_map)

        write(filename_l2b_tb_grid,261) trim(filename_l2b_tb_base)
        261 format(a,'.time.nc')
        print *, filename_l2b_tb_grid
        call write_resampled_time_map(gridded_time_map,filename_l2b_tb_grid)

        do ichan = 1,size(channel_list)
            if (channel_list(ichan) .lt. 1) cycle
            if (do_freq_rect((channel_list(ichan)+1)/2)) then
                channel = channel_list(ichan)
                if (channel == 0) cycle
                print *,'Resampling onto lat/lon grid, channel: ',channel
                call resample_onto_lat_lon_grid_interp(0.25,channel,.true.,gridded_time_map,gridded_Tb_map)

                write(filename_l2b_tb_grid,161) trim(filename_l2b_tb_base),channel,footprint_size
                161 format(a,'.grid_tb.ch',i2.2,'.',i3.3,'km.nc')
                print *, filename_l2b_tb_grid

                call write_resampled_maps(gridded_Tb_map,       &
                                        filename_l2b_tb_grid)
            endif
        enddo
    else
        print *,'All rect grid files exist - skipping'
    endif !use rectangular grid

    if ((use_NP_grid) .and. any(do_freq_polar)) then
        do_time = .true.
        do ichan = 1,size(channel_list)
            if (channel_list(ichan) == 0) cycle
            if (do_freq_polar((channel_list(ichan)+1)/2)) then
                channel = channel_list(ichan)
                if (channel == 0) cycle
                print *,'Resampling onto north polar grid, channel: ',channel
                call resample_onto_polar_map_interp('north',  &
                                                    '25km', &
                                                    channel, &
                                                    50, & !num_overlap_scans
                                                    do_time, &
                                                    polar_gridded_time_map, &
                                                    polar_gridded_tb_map, &
                                                    error)
                
                if (do_time) then
                    ! write the time map
                    print *,'Writing the time map'
                    write(filename_l2b_tb_grid,471) trim(filename_l2b_tb_base)
                    471 format(a,'.polar_grid_time.north.nc')
                    print *, trim(filename_l2b_tb_grid)
                    call polar_gridded_time_map%write_map_netcdf(filename_l2b_tb_grid,-999.0D0,error)
                    do_time = .false. !time is done for this orbit, so don't do it again for other channels
                endif

                ! write the Tb map
                print *,'writing the polar-gridded tb map'
                write(filename_l2b_tb_grid,461) trim(filename_l2b_tb_base),channel,footprint_size
                461 format(a,'.polar_grid_tb.north.ch',i2.2,'.',i3.3,'km.nc')
                
                print *, trim(filename_l2b_tb_grid)
                call polar_gridded_tb_map%write_map_netcdf(filename_l2b_tb_grid,-999.0,error)
            endif
        enddo
    else
        print *,'All polar grid files exist - skipping'
    endif

    if ((use_SP_grid) .and. any(do_freq_polar)) then
        do_time = .true.
        do ichan = 1,size(channel_list)
            if (channel_list(ichan) == 0) cycle
            if (do_freq_polar((channel_list(ichan)+1)/2)) then
                channel = channel_list(ichan)
                if (channel == 0) cycle
                print *,'Resampling onto south polar grid, channel: ',channel
                call resample_onto_polar_map_interp('south',  &
                                                    '25km', &
                                                    channel, &
                                                    50, & !num_overlap_scans
                                                    do_time, &
                                                    polar_gridded_time_map, &
                                                    polar_gridded_tb_map, &
                                                    error)
                
                if (do_time) then
                    ! write the time map
                    print *,'Writing the time map'
                    write(filename_l2b_tb_grid,571) trim(filename_l2b_tb_base)
                    571 format(a,'.polar_grid_time.south.nc')
                    print *, trim(filename_l2b_tb_grid)
                    call polar_gridded_time_map%write_map_netcdf(filename_l2b_tb_grid,-999.0D0,error)
                    do_time = .false. !time is done for this orbit, so don't do it again for other channels
                endif

                ! write the Tb map
                print *,'writing the polar-gridded tb map'
                write(filename_l2b_tb_grid,561) trim(filename_l2b_tb_base),channel,footprint_size
                561 format(a,'.polar_grid_tb.south.ch',i2.2,'.',i3.3,'km.nc')
                
                print *, trim(filename_l2b_tb_grid)
                call polar_gridded_tb_map%write_map_netcdf(filename_l2b_tb_grid,-999.0,error)
            endif
        enddo
    else
        print *,'All polar grid files exist - skipping'
    endif

    end subroutine
    