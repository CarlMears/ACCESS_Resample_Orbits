 
 ! this program processes one amsr2 L1A orbit, and outputs various 
 ! files with brightness temperatures 
 !
 !
 ! example command
 !      
 ! start "17000 70km" cmd /K amsr2_exe_test_logic.exe 17000 17999 13  70 1  1  1 1 2 3 4 5 6                           
 !                                                     strt  end ksat sz r np sp Channels
 !
 ! sz = footprint size
 ! r means do rectangular grid
 ! np means do north polar grid
 ! sp means do south polar grid
 ! channels: 1 = 19.3 GHz V
 !           2 = 19.3 GHz H
 !           3 = 22.235 GHz H
 !           4 = 22.235 GHz V
 !			 5 = 37.0 GHz H
 !			 6 = 37.0 GHz V

 program process_ssmi_l1c_orbit_tb_ACCESS   

    use, intrinsic :: iso_fortran_env, only: int8, int32, real32, real64 
    use QuadrilateralUtilties, only: triangle_area,is_inside_triangle

    implicit none 

    character(256)  :: acmd(100)
    integer(int32)  :: numarg
    integer(int32)  :: iorbit,iorbit_start,iorbit_end
	integer(int32)  :: ksat
    integer(int32)  :: footprint_size
    integer(int32)  :: nchannel
    integer(int32)  :: channel,i
    integer(int32)  :: use_rectangular_grid_int
    integer(int32)  :: use_NP_grid_int, use_SP_grid_int
    
    integer(int32),dimension(20) :: channel_list
    logical(4)         :: write_swath_tbs
    logical(4)         :: overwrite_existing_files
    logical(4)         :: use_rectangular_grid  
    logical(4)         :: use_NP_grid, use_SP_grid
    integer(int32)  :: num_overlap_scans
    character(150)  :: pathname_tb_output
    character(150)  :: pathname_ssmi_l1c_input

    ! read command line
    call get_cmd_arg(numarg,acmd)
    read(acmd(1),'(i6)') iorbit_start
    read(acmd(2),'(i6)') iorbit_end
	read(acmd(3),'(i6)') ksat
    read(acmd(4),'(i3)') footprint_size
    read(acmd(5),'(i3)') num_overlap_scans
    read(acmd(6),'(i1)') use_rectangular_grid_int
    read(acmd(7),'(i1)') use_NP_grid_int
    read(acmd(8),'(i1)') use_SP_grid_int
    
    nchannel = numarg-6
    channel_list = 0
    if (nchannel .gt. 0) then
        do i = 1,numarg-8
            read(acmd(i+8),'(i6)') channel
            channel_list(i) = channel
        enddo
    endif

    use_rectangular_grid = .false.
    use_NP_grid = .false.
    use_SP_grid = .false.
    if (use_rectangular_grid_int > 0) use_rectangular_grid=.true.
    if (use_NP_grid_int > 0) use_NP_grid=.true.
    if (use_SP_grid_int > 0) use_SP_grid=.true.
    
    print *, 'start_orbit: ',iorbit_start
    print *, 'end_orbit: ',iorbit_end
    print *, 'ksat: ',ksat
    print *, 'footprint_size: ',footprint_size
    print *, 'num_overlap_scans: ',num_overlap_scans
    print *, 'channel_list'
    print *,channel_list

    if (use_rectangular_grid) print *,'resampling to rectangular grid'
    if (use_NP_grid) print *,'resampling to north polar grid'
    if (use_SP_grid) print *,'resampling to south polar grid'

    write_swath_tbs = .false.
    overwrite_existing_files = .false.
    
    !pathname_tb_output = 'L:\access\ssmi_tb_orbits\  '
    pathname_tb_output = 'B:\_access_temp\ssmi_tb_orbits\'
    write(pathname_ssmi_l1c_input,111) ksat
    !111 format('S:\SSMI\F',i2.2,'\L1C\')
    111 format('\\ocean\S\SSMI\F',i2.2,'\L1C_V08\')
    do iorbit = iorbit_start,iorbit_end
        call process_and_resample_ssmi_l1c_orbit(ksat,              &
                                                iorbit,             &
                                                channel_list,       &
                                                footprint_size,     &
                                                write_swath_tbs,    &
                                                overwrite_existing_files, &
                                                use_rectangular_grid,  &    
                                                use_NP_grid,  &
                                                use_SP_grid,  &
                                                num_overlap_scans, &
                                                trim(pathname_tb_output), &
                                                trim(pathname_ssmi_l1c_input))
    enddo
 end program process_ssmi_l1c_orbit_tb_ACCESS   

 subroutine process_and_resample_ssmi_l1c_orbit( ksat,          &
                                                 iorbit,        &                  ! orbit to process
                                                 channel_list,  &            ! list of AMSR2 channels to process
                                                 footprint_size,  &
                                                 write_swath_tbs,  &         ! set true to write swath tb file
                                                 overwrite_existing_files, & ! if true, existing files are overwritten
                                                 use_rectangular_grid,  &    ! set true to resample to rectangular grid
                                                 use_NP_grid,        &       ! set true to resample to north polar grid
                                                 use_SP_grid,        &       ! set true to resample to south polar grid
                                                 num_overlap_scans, &        ! number of scans to ignore at beginning and end of orbit
                                                 pathname_tb_output, &       ! path to output files
                                                 pathname_ssmi_l1c_input)    ! path to input files

    use, intrinsic :: iso_fortran_env, only: int8, int32, real32, real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan,ieee_is_finite
  
    ! use RSS_SSMI_L1C, only:  read_rss_ssmi_l1c,read_ssmi_l1c_binary_head, &
    !                     level1a_rss_header, level1a_rss, & ! types
    !                     L1arss_hd,L1arss                ! module variables

    use sat_id_module, only: ssmi_nfreq => nfreq,ssmi_channel_name => channel_name
     
    !use l2_tb_module_ACCESS

	use ssmi_filename_routines, only: get_ssmi_l1c_filename,get_ssmi_l2b_filename_base
    !use filename_routines, only: get_l2b_filename_base
    use FileExist, only: file_exist,file_exists  !file_exists is the function form
    use open_file_routines, only: open_binary
    use date_utilities, only: find_month_day,fd_time_2000,fd_date_2000,days_in_month,doy_from_year_month_day
    use netcdf
    use L2_Tb_io, only: write_geoloc_extra_netcdf
    !use l2_module_extra
    use resample, only: init_resample_weights,ResampleWeights
    use resample, only: resample_tb_extra,resample_tb_extra_85
    use resample, only: free_resample_weights
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

    use ssmi_l1c, only: read_l1c_file
    use ssmi_l1c, only: numscan
    use ssmi_l1c, only: iqual_flag,scan_time,orbit,zang,scpos,scvel
    use ssmi_l1c, only: cellat,cellon,celtht,celphi,celsun,celrfi
    use ssmi_l1c, only: ta_nat,ta_rsp,tb_nat,tb_rsp  !(maxchn,maxcel,maxscan)
    use extra_locs, only: ssmi_geolocation_extra_lf                             
    implicit none    

    integer(int32),intent(in) :: ksat
    integer(int32),intent(in)  :: iorbit
    integer(int32),dimension(20),intent(in) :: channel_list
    integer(int32),intent(in) :: footprint_size
    logical(4),intent(in) :: write_swath_tbs
    logical(4),intent(in) :: overwrite_existing_files
    logical(4),intent(in) :: use_rectangular_grid
    logical(4),intent(in) :: use_NP_grid
    logical(4),intent(in) :: use_SP_grid
    integer(int32) :: num_overlap_scans
    character(len=*) :: pathname_tb_output
    character(len=*) :: pathname_ssmi_l1c_input
    
    logical,dimension(ssmi_nfreq) :: do_freq_requested
    logical,dimension(ssmi_nfreq) :: do_freq_any
    logical,dimension(ssmi_nfreq) :: do_freq_rect
    logical,dimension(ssmi_nfreq) :: do_freq_polar
    
    
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

    integer(int32)  :: l1c_read_error
    
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
    !real(real32) ta89_resp(maxcel,maxscan,13:16)
     
    character(256) :: acmd(100)
    integer(int32)     :: numarg
    
    character(150)  :: test_filename

    character(150)  :: filename_l1c_tb
   
    character(150)  :: filename_l2b_tb_base_GL
    character(150)  :: filename_l2b_tb_base_NP
    character(150)  :: filename_l2b_tb_base_SP

    character(150)  :: filename_l2b_tb_GL_nc
    character(150)  :: filename_l2b_tb_NP_nc
    character(150)  :: filename_l2b_tb_SP_nc

    character(150)  :: filename_l2b_tb_GL_resamp
    character(150)  :: filename_l2b_tb_NP_resamp
    character(150)  :: filename_l2b_tb_SP_resamp

    character(150)  :: filename_l2b_tb_GL_grid
    character(150)  :: filename_l2b_tb_NP_grid
    character(150)  :: filename_l2b_tb_SP_grid

    character(150)  :: filename_l2b_tb_GL_dist
    character(150)  :: filename_l2b_tb_NP_dist
    character(150)  :: filename_l2b_tb_SP_dist

    ! character(150)  :: filename_l2b_tb_nc
    ! character(150)  :: filename_l2b_tb_resamp
    ! character(180)  :: filename_l2b_tb_grid 
    ! character(180)  :: filename_l2b_tb_dist 
    

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
    
    !type(tb_scn) :: tb_scan
    integer(int32)   :: num_scans
    integer(int32)   :: min_valid_extra_scan
    integer(int32)   :: max_valid_extra_scan
    
    character(10) :: overwrite_command

    type(ResampleWeights) :: ssmi_wts

    type(map_str)   :: gridded_Tb_map
    type(map_str)   :: gridded_distance_map
    type(map_str64) :: gridded_time_map

    type(earth_ref_map_str32) :: polar_gridded_Tb_map
    type(earth_ref_map_str64) :: polar_gridded_time_map

    real(real32)  :: grid_lon,grid_lat,search_size
    integer(int32) :: scan_closest,fov_closest
    real(real32) :: distance
    integer(int32) :: channel
    ! common/resampled89/ ta89_resp

    logical :: do_time

    real(real64),allocatable :: lat_array(:,:)
    real(real64),allocatable :: lon_array(:,:)
    integer(int32)            :: num_x,num_y,error,ok
  
    ! call define_filenames
    
    if (overwrite_existing_files) then
        overwrite = 1
    else 
        overwrite = 0
    endif


    call get_ssmi_l1c_filename(pathname_ssmi_l1c_input, ksat, iorbit, filename_l1c_tb)
    print *,pathname_ssmi_l1c_input
    print *,filename_l1c_tb
    !call get_ssmi_l2b_filename_base(pathname_tb_output, ksat, iorbit, filename_l2b_tb_base)
    call get_ssmi_l2b_filename_base(pathname_tb_output, &
                                    ksat,               &
                                    iorbit,             &
                                    'GL',               &
                                    footprint_size,                 &
                                    filename_l2b_tb_base_GL)

    call get_ssmi_l2b_filename_base(pathname_tb_output, &
                                ksat,               &
                                iorbit,             &
                                'NP',               &
                                footprint_size,                 &
                                filename_l2b_tb_base_NP)

    call get_ssmi_l2b_filename_base(pathname_tb_output, &
                                ksat,               &
                                iorbit,             &
                                'SP',               &
                                footprint_size,                 &
                                filename_l2b_tb_base_SP)

                      
    filename_l2b_tb_GL_nc = trim(filename_l2b_tb_base_GL)//'.nc'
    filename_l2b_tb_NP_nc = trim(filename_l2b_tb_base_NP)//'.nc'
    filename_l2b_tb_SP_nc = trim(filename_l2b_tb_base_SP)//'.nc'

    filename_l2b_tb_GL_resamp = trim(filename_l2b_tb_base_GL)//'.resamp.nc'
    
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
    print *,do_freq_requested
    print *,channel_list
    do ichan = 1,size(channel_list)
        if (channel_list(ichan) .lt. 1) cycle
        if (do_freq_requested((channel_list(ichan)+1)/2)) then
            if (overwrite) then
                do_freq_rect((channel_list(ichan)+1)/2) = .true.
                do_freq_polar((channel_list(ichan)+1)/2) = .true.
                do_freq_any((channel_list(ichan)+1)/2) = .true.
            else
                channel = channel_list(ichan)

                ! Global Grid
                write(filename_l2b_tb_GL_grid,&
                        "(a,'.grid_tb.ch',i2.2,'.',i3.3,'km.nc')") &
                        trim(filename_l2b_tb_base_GL),channel,footprint_size
                print *, "GL Grid: ",filename_l2b_tb_GL_grid
                if (.not. file_exists(filename_l2b_tb_GL_grid)) then
                    do_freq_rect((channel_list(ichan)+1)/2) = .true.
                    do_freq_any((channel_list(ichan)+1)/2) = .true.
                endif

                ! North Polar Grid
                write(filename_l2b_tb_NP_grid, &
                                "(a,'.polar_grid_tb.north.ch',i2.2,'.',i3.3,'km.nc')") &
                                trim(filename_l2b_tb_base_NP),channel,footprint_size
                print *,"NP Grid: ",filename_l2b_tb_NP_grid
                if (.not. file_exists(filename_l2b_tb_NP_grid)) then
                    do_freq_polar((channel_list(ichan)+1)/2) = .true.
                    do_freq_any((channel_list(ichan)+1)/2) = .true.
                endif

                ! South Polar Grid
                write(filename_l2b_tb_SP_grid, &
                                "(a,'.polar_grid_tb.south.ch',i2.2,'.',i3.3,'km.nc')") &
                                trim(filename_l2b_tb_base_SP),channel,footprint_size
                print *,"SP Grid: ",filename_l2b_tb_SP_grid
                if (.not. file_exists(filename_l2b_tb_SP_grid)) then
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

    print *,'Reading L1C file: ',trim(filename_l1c_tb)
    call read_l1c_file(ksat,             &
                       iorbit,           &
                       filename_l1c_tb,  &
                       start_time,       &
                       l1c_read_error)
        
    call fd_date_2000(start_time, secyr,lyear,idayjl,imon,idaymo,secdy)
    write(*,'(i3,i6,i5,2i3,f7.2)') ksat,iorbit,lyear,imon,idaymo,secdy/3600.
        
    if(numscan.lt.35) then
        print *,'not enough scans to do scan averaging and resampling'
        return
    endif
        
    ! orbit seems worth processing.  Open the L2B brightness temperature file

    print *,'Calculating Locations of Extra FOVs'  
    call ssmi_geolocation_extra_lf(min_valid_extra_scan, &
                                   max_valid_extra_scan)

    print *,'Min and Max Valid Extra Scans: ',min_valid_extra_scan,max_valid_extra_scan
    print *
    print *,'Resampling Tbs to the Extra FOVs'  

    do ifreq = 1,4
        if (.not. (do_freq_any(ifreq))) cycle
        call init_resample_weights('SSMI',ifreq-1,footprint_size,ssmi_wts)
        !call write_wt_slice_netcdf(amsr2_wts,trim(freq_name(ifreq)),1,243)
        if (ifreq <= 3) then
            print *,'Resampling Channel: ',(2*ifreq)-1
            call resample_tb_extra(ssmi_wts,(2*ifreq)-1)  ! Vpol
            print *,'Resampling Channel: ',(2*ifreq)
            call resample_tb_extra(ssmi_wts,2*ifreq)      ! Hpol
        else
            print *,'Resampling Channel: ',(2*ifreq)-1
            call resample_tb_extra_85(ssmi_wts,(2*ifreq)-1)  ! Vpol
            print *,'Resampling Channel: ',(2*ifreq)
            call resample_tb_extra_85(ssmi_wts,2*ifreq)      ! Hpol
        endif
        call free_resample_weights(ssmi_wts)
    enddo

    ! write the resampled orbit -- this is done for debugging
    if (write_swath_tbs) then
        print *,'Writing Resampled Orbit'
        print *,trim(filename_l2b_tb_GL_resamp)
        call write_geoloc_extra_netcdf(trim(filename_l2b_tb_GL_resamp),iorbit,overwrite)
    endif
    
    if ((use_rectangular_grid) .and. any(do_freq_rect)) then
        print *,'Resampling time onto rectangular lat/lon grid'
        call resample_time_onto_lat_lon_grid(0.25, &
                                            gridded_time_map, &
                                            gridded_distance_map, &
                                            min_valid_extra_scan, &
                                            max_valid_extra_scan, &
                                            num_overlap_scans)  

        write(filename_l2b_tb_GL_grid,261) trim(filename_l2b_tb_base_GL)
        261 format(a,'.time.nc')
        print *, filename_l2b_tb_GL_grid
        call write_resampled_time_map(gridded_time_map,filename_l2b_tb_GL_grid)

        do ichan = 1,size(channel_list)
            if (channel_list(ichan) .lt. 1) cycle
            if (do_freq_rect((channel_list(ichan)+1)/2)) then
                channel = channel_list(ichan)
                if (channel == 0) cycle
                print *,'Resampling onto lat/lon grid, channel: ',channel
                call resample_onto_lat_lon_grid_interp( 0.25, &
                                                        channel,  &
                                                        .true.,   &
                                                        gridded_time_map, &
                                                        gridded_Tb_map, &
                                                        min_valid_extra_scan, &
                                                        max_valid_extra_scan, &
                                                        num_overlap_scans)

                write(filename_l2b_tb_GL_grid,161) trim(filename_l2b_tb_base_GL),channel,footprint_size
                161 format(a,'.grid_tb.ch',i2.2,'.',i3.3,'km.nc')
                print *, filename_l2b_tb_GL_grid

                call write_resampled_maps(gridded_Tb_map,       &
                                          filename_l2b_tb_GL_grid)
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
                                                    do_time, &
                                                    polar_gridded_time_map, &
                                                    polar_gridded_tb_map, &
                                                    min_valid_extra_scan, &
                                                    max_valid_extra_scan, &
                                                    num_overlap_scans, &
                                                    error)
                
                if (do_time) then
                    ! write the time map
                    print *,'Writing the time map'
                    write(filename_l2b_tb_NP_grid,471) trim(filename_l2b_tb_base_NP)
                    471 format(a,'.polar_grid_time.north.nc')
                    print *, trim(filename_l2b_tb_NP_grid)
                    call polar_gridded_time_map%write_map_netcdf(filename_l2b_tb_NP_grid,-999.0D0,error)
                    do_time = .false. !time is done for this orbit, so don't do it again for other channels
                endif

                ! write the Tb map
                print *,'writing the polar-gridded tb map'
                write(filename_l2b_tb_NP_grid,461) trim(filename_l2b_tb_base_NP),channel,footprint_size
                461 format(a,'.polar_grid_tb.north.ch',i2.2,'.',i3.3,'km.nc')
                
                print *, trim(filename_l2b_tb_NP_grid)
                call polar_gridded_tb_map%write_map_netcdf(filename_l2b_tb_NP_grid,-999.0,error)
            endif
        enddo
    else
        print *,'All NP polar grid files exist - skipping'
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
                                                    do_time, &
                                                    polar_gridded_time_map, &
                                                    polar_gridded_tb_map, &
                                                    min_valid_extra_scan, &
                                                    max_valid_extra_scan, &
                                                    num_overlap_scans, &
                                                    error)
                
                if (do_time) then
                    ! write the time map
                    print *,'Writing the time map'
                    write(filename_l2b_tb_SP_grid,571) trim(filename_l2b_tb_base_SP)
                    571 format(a,'.polar_grid_time.south.nc')
                    print *, trim(filename_l2b_tb_SP_grid)
                    call polar_gridded_time_map%write_map_netcdf(filename_l2b_tb_SP_grid,-999.0D0,error)
                    do_time = .false. !time is done for this orbit, so don't do it again for other channels
                endif

                ! write the Tb map
                print *,'writing the polar-gridded tb map'
                write(filename_l2b_tb_SP_grid,561) trim(filename_l2b_tb_base_SP),channel,footprint_size
                561 format(a,'.polar_grid_tb.south.ch',i2.2,'.',i3.3,'km.nc')
                
                print *, trim(filename_l2b_tb_SP_grid)
                call polar_gridded_tb_map%write_map_netcdf(filename_l2b_tb_SP_grid,-999.0,error)
            endif
        enddo
    else
        print *,'All SP polar grid files exist - skipping'
    endif

    end subroutine
    