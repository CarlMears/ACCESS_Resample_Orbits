!      3/8/2017 changed 3/26/2017. filename_zang added. dta_zang(0:99,2,nfreq) added


!      3/7/2017 changed 3/8/2017 values for tcold_offset and spillover_coldmirror changed

!      3/2/2017 changed 3/7/2016.  rfi_cold_mirror_map(360,80,2,2) changed to rfi_cold_mirror_map(360,80,6)

!      below are comments before adding suffix v8 to filename
!      7/17/2014 changed 3/2/2017.  tht_offset changed


!      7/10/2014 version changed 7/17/2014.  filename_2be and filename_inv added

!      1/4/2014 version updated 7/10/2014.  array corr_tab_11vlna(0:359,0:359) and filename_lna added

!      11/9/2013 version changed 1/4/2014.  rss_l1a_tobe_processed(maxrev) removed. it was doing nothing

!      11/4/2013 version changed 11/9/2013.  Values for thot_offset changed

!      9/23/2013 version changed 11/4/2013.  delta_apc and chi_apc added


!      9/20/2013 version changed 9/23/2013.  counts_per_deg added
!      tcold_offset is now a function of freq with 0.1 being added to 19-89 GHz

!      9/19/2013 version changed 9/20/2013.  arrays reflat and reflon added.
!      geosync_lon and nsat_rfi and array celrfi removed

!      9/18/2013 version changed 9/19/2013.  filename_crfi, rfi_cold_mirror_map(360,80,2,2) and rfi_cold_count(2,maxscan) added

!      9/10/2013 version changed 9/18/2013.  filename_moon added 
!      moon in cold mirror parameters added and iflag_moon
!      spill_adj and spill_xscan added
!      maxcel_max replaced by maxcel everywhere

!      9/9/2013 version changed 9/10/2013. spillover_coldmirror values fined tuned
!      base on analysis. see 'O:\amsr2\cold_counts\memo2.txt'



!      8/28/2013 version changed 9/9/2013.  array ta_spill and gainest are added and
!      a trial value of 0.0015 is now being used for 11 ghz spillover (was 0)
!      iflag_cal had been mistakenly called integer(1).  it was change to integer(2)
!      need to determine the effect of this mistake on previous processing

!      8/23/2013 version change 8/28/2013.  tamap(360,180,12,6,0:1) and ta_map(360,180,6,0:1) added
!      filename_tamap replaces filename_cband_map; filename_moon removed


!      8/22/2013 version change 8/23/2013
!      geoloc_jaxa(4,maxcel_89,maxscan_89) changed to geoloc_jaxa(4,maxcel_89,maxscan)
!      this was a mistake.  it had no effect but was wasting a lot of memory

!      6/26/2013 version changed 8/22/2013
!      a number of variables and arrays were removed and added
!      see 'O:\amsr2\L2_processing\memo3.txt' for a list of changes


!      May 3 2013 version updated 6/26/2013. tht_offset changed. see 'O:\amsr2\abs_cal\memo6.txt'
!      tht_offset only used by ta_inter_calibration_v7.f


!      april 25 2013 version changed may 3 2013. earthcounts_jaxa added

!      march 13 2013 version changed april 25 2013. irx_temp_cnt added

!      march 6 2013 version changed march 13 2013
!      dimension of xscan_tab is now maxcel_89 rather than maxcel
!      this change has no effect on previous processing because xscan had been set to zero


!      march 6, 2013 version changed on march 8 2013.  following arrays were added
!      geoloc_jaxa, interp_flag,thot_jaxa   

!      feb 11 2013 version changed march 6 2013
!      1. beta_term array no longer need and is removed
!      2. acoef_nl replaces beta_quad and is readin via 'o:\amsr2\routines\readin_data_files.f'
!      3. hot target offset thot_offset(nfreq) is added
!      4. filename_nl added
!      5. nfreq and nch move to amsr_id_module
!      5. therm_sps removed. it was not being used


    module l2_module
        use sat_id_module
        implicit none

        !    earth dimensions match wgs84    
        real(8),    parameter :: re=6378.137d3        !meters
        real(8),    parameter :: rp=6356.752d3        !meters
        real(8),    parameter :: ffac=(rp/re)*(rp/re)
        real(8),    parameter :: geosync_alt   =42164.0e3      !satellite orbits, montenbruck and gill, page 294
        !      integer(4), parameter :: nsat_rfi=8
        !      real(8),    parameter :: geosync_lon(nsat_rfi)=(/13., 19.2, 28.2, 317.0, 352.8, 10.0, 257.2, 260.8/)

        integer(4), parameter :: ilu_l2b  =  4

        integer(4), parameter :: ngrp=3

        !    integer(4), parameter :: maxcel_max=maxcel !if you change this to maxcel=maxcel_89 be careful with the array below
        integer(4), parameter :: ncel_rsp=maxcel  !number of cells for resampling program    
        
        integer(4), parameter :: iminc=-2047  !min for radiometer counts
        integer(4), parameter :: imaxc= 2046  !max for radiometer counts
        
        real(4), parameter :: ta_earth_min(nch)=(/ 55., 55., 55., 55., 55., 55., 55., 55., 55., 55., 55., 55., 55., 55., 55., 55./)
        real(4), parameter :: ta_earth_max(nch)=(/330.,330.,330.,330.,330.,330.,330.,330.,330.,330.,330.,330.,330.,330.,330.,330./)

        real(4), parameter :: tcld_plk(nfreq)    =(/2.73, 2.73, 2.74, 2.75, 2.77, 2.82, 3.24, 3.24/) !tcold for apc see 'o:\emiss_model\planck.txt'

        !      on 9/23/2013 i added 0.1 to tcold_offset for 19-89 ghz because these channels have no spillover adjustment      
        !      real(4), parameter :: tcold_offset(nfreq)=(/0.30, 0.30, 0.30, 0.40, 0.40, 0.40, 0.40, 0.40/) 
        real(4), parameter :: tcold_offset(nfreq)=(/0.00, 0.00, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10/) !inserted 3/8/2017

        !      real(4), parameter :: thot_offset(nfreq)=(/-1.8,-1.8,-1.8,-1.8,-1.8,-1.8,-0.8,-0.8/)  !values before 11/9/2013 change
        real(4), parameter :: thot_offset(nfreq)=(/-1.3,-1.3,-1.3,-1.3,-1.3,-1.3,-1.3,-1.3/)

        !      this is now defined a little dif from amsre.  it is only used by the ta_inter_calibration routine
        !      it is the dif between the average eia for given channel minus the average eia for 37 ghz.  
        !      currently set to 0 for amsr2 except for 89b where the eia offsets was computed to be -0.51 using chi formula
        real(4), parameter :: tht_offset(nfreq)=(/0., 0., 0., 0., 0., 0., 0., -0.51/)    

        !      tht_offset for 6.9 and 10.7 set to amsr-e values.  this is necessary because amsr-e algo tables
        !      assume this value for the tht offset.  7.3 ghz is set to 6.9 value for now
        !      see 'O:\amsr2\abs_cal\memo6.txt'
        !    real(4), parameter :: tht_offset(nfreq)=(/0.1316, 0.1316, 0.0920, 0., 0., 0., 0., -0.65/)    !this was used before 3/2/2017

        !      these limits are currently wide open for amsr2.  need to constrain them more.
        real(4), parameter :: therm_min=280.0
        real(4), parameter :: therm_max=300.0
        real(4), parameter :: therm_var=90.0 ! 90.0 is equivalent to sqrt(var)*sqrt(10./9.)>10
    
        !      see 'O:\amsr2\cold_counts\memo3.txt' for following values
        real(4), parameter :: counts_per_deg(nch)= &
            (/11.10, 10.63, 11.54, 11.07, 10.13, 10.12,  9.93,  9.51, 10.63,  9.81,  9.67,  9.87,  9.97,  9.69, 10.58,  9.66/)

        !    real(4), parameter :: spillover_coldmirror(nch)=(/ 0.0042,0.0035,0.0040,0.0030,0.0005,0.0005,0.,0.,0.,0.,0.,0.,0.,0.,0.,0./)    
        real(4), parameter :: spillover_coldmirror(nch)=(/ 0.0042,0.0035,0.0040,0.0030,0.0000,0.0000,0.,0.,0.,0.,0.,0.,0.,0.,0.,0./)    !inserted 3/8/2017

        real(4), parameter :: spill_xscan(low_cal_obs)= &
            (/0.840,0.840,0.841,0.845,0.851,0.862,0.878,0.901,0.931,0.970,1.018,1.077,1.147,1.231,1.328,1.440/)  ! 'O:\amsr2\cold_counts\memo2.txt'
      
        integer(4), parameter ::  iradius_fieldavg=35  !see 'O:\ssmis\L2_processing\memo9.txt' for discussion of search radius


        character(120) filename_leap
        character(120) filename_l1a_inventory,filename_rss_l1a_inventory,filename_orb,filename_node

        character(120) filename_bad,filename_xsc,filename_apc,filename_ta_resample_wts,filename_tbl
        character(120) filename_geo,filename_cnt
        character(120) filename_ta89_resample_wts
        character(120) filename_tamap
        character(120) filename_nl,filename_zang
        character(120) land_filename,ice_filename,climate_filename,sss_filename,filename_percent_land
        character(120) filename_moon,filename_crfi,filename_lna,filename_2be,filename_inv

        character( 60) pathname_l1a,pathname_l2b
        character(120) pathname_epmaps

        real(8) ut1_1993(maxscan) !really tai_1993. it is a placeholder for heritage code

        !      ====================================================================================================
        !      ========================= tables read in by readin_data_files ======================================
        !      ====================================================================================================
        integer(4) itai_1993_leap(100) !table of tai_1993 when leap seconds were added
        integer(4) nbad,iorbit_bad(10000)
        real(8) eia_convert_37to19,beta(ngrp),azoffset(ngrp),start_time_a,start_time_b,spillover_time,rolloffset,rpy_error
        real(4) xscan_cal(2,nch,high_cal_obs),stddv_avg(2,nch)
        real(8) acoef_nl(5,nch)
        real(4) xscan_tab(nch,maxcel_89)
        real(4) tbmin(2,nfreq),tbmax(2,nfreq)
        character(1) percent_land_tab(4320,2160,5)  !land table
        character(1) amask_ice(45,180,12)           !monthly climatology sea ice map
        integer(1) isal_corr(10,360,180,12)
        real(4)  tamap(360,180,12,6,0:1)
        real(4) amp_moon(nfreq),phase_factor_moon(nfreq),bw_cold(high_cal_obs,nfreq)
        real(4) alpha_cold(high_cal_obs,nfreq),beta_cold(high_cal_obs,nfreq)
        real(4) rfi_cold_mirror_map(360,80,6)
        real(4) delta_apc(2,nfreq),chi_apc(2,nfreq)
        real(4) corr_tab_11vlna(0:359,0:359)
        real(4) dta_zang(0:99,nch)
      
        ! l1a file content
        integer(4) ksat,numscan
        real(8) scan_time(maxscan),scpos(3,maxscan),scvel(3,maxscan),orbit(maxscan)
        real(8) scrpy(3,maxscan)
        real(4) therm(10,maxscan)  !amsre was 8
        integer(4) irx_temp_cnt(nfreq,maxscan)
    
        ! jaxa data need to testing
        integer(2) interp_flag(32,nch,maxscan)
        real(4)    thot_jaxa(nch,maxscan)
        real(4) geoloc_jaxa(4,maxcel_89,maxscan)  !4 denotes lat89a, lat89b, lon89a, lona9b

        integer(4) iscan_flag(maxscan)


        !      geolocation output
        real(4) omega(maxscan),omega_smooth(maxscan),scan_index(maxscan)
        real(8) scpos_eci0(3,maxscan), x_eci0(3,maxscan), y_eci0(3,maxscan), z_eci0(3,maxscan),sunvec0(3,maxscan),moonvec0(3,maxscan)
        real(8) zang(maxscan),alpha_sun(maxscan),beta_sun(maxscan)
        integer(1) kflag_sun(maxscan)

        real(4) sunlat(maxscan),sunlon(maxscan),sundis(maxscan)
        real(4) scloc(3,maxscan)
        real(4) cellat(maxcel,maxscan),cellon(maxcel,maxscan),celtht(maxcel,maxscan)
        real(4) celphi(maxcel,maxscan),celsun(maxcel,maxscan)
        real(4) reflat(maxcel,maxscan),reflon(maxcel,maxscan) !geostationary lat and lon of reflected ray

        real(4) cellat_89(maxcel_89,maxscan_89),cellon_89(maxcel_89,maxscan_89),cellat_cold(maxscan),cellon_cold(maxscan)
        real(4) celtht_89(maxcel_89,maxscan_89),celphi_89(maxcel_89,maxscan_89),celrng_89(maxcel_89,maxscan_89)
        real(4) ta_cold(nch,maxscan),ta_hot(nch,maxscan),ta_spill(nch,maxscan),gainest(nch,maxscan)
        
    !      moon in cold mirror parameters
        integer(1) iflag_moon(nch,high_cal_obs,maxscan)

    !      these arrays relate to cold count filtering and averaging.  leading index of 2 means:  1=cold, 2=hot
        integer(1) lflag_counts(nch,high_cal_obs,maxscan)
        integer(4) num_count(2,nch,maxscan)
        real(4)    avg_count(2,nch,maxscan),rms_count(2,nch,maxscan),spill_adj(nch,maxscan)
        real(4) rfi_cold_count(2,maxscan)  

        integer(2) rx_offset_gain(nch, 2,maxscan)
        integer(2) earthcounts_jaxa(nch,maxcel_89,maxscan)
        integer(2) earthcounts(nch,maxcel_89,maxscan)
        integer(2) csmcounts(nch,high_cal_obs,maxscan),htscounts(nch,high_cal_obs,maxscan)
        
    !    real(4) jaxa_tb(nch,maxcel_89,maxscan)

        integer(1) iflag_l0(maxscan),iflag_rx(nch,maxscan)
        integer(2) iflag_cal(nch,maxscan)

        real(4) ta(maxcel,maxscan,12),tar(maxcel,maxscan,24),ta_89(maxcel_89,maxscan_89,2)

    !      ocean  block
        integer(1) isurcel(5,maxcel,maxscan),ice_flag(3,maxcel,maxscan), ice_flag2(maxcel,maxscan)
        integer(1) ioob_flag(5,maxcel,maxscan), iret_flag(4,maxcel,maxscan), irainadj(4,maxcel,maxscan)

        real(4) sstcl(maxcel,maxscan),wincl(maxcel,maxscan),vapcl(maxcel,maxscan)
        real(4) phir(maxcel,maxscan) !wgcm no longer use, if ncep is available ncep speed is stored in wincl
        real(4) ep(0:10,maxcel,maxscan),sst_reg(maxcel,maxscan)
        
        real(4) tbmd_sav(6,maxcel,maxscan),tbhi_sav(2,maxcel,maxscan)

        integer(1) ice_map(1440,720,2)
        real(4) ta_map(360,180,6,0:1)
    end module l2_module

    ! module l2_module_extra
    !     use sat_id_module
    !     implicit none

    !     !    earth dimensions match wgs84    
    !     !real(8),    parameter :: re=6378.137d3        !meters
    !     !real(8),    parameter :: rp=6356.752d3        !meters
    !     !real(8),    parameter :: ffac=(rp/re)*(rp/re)
    !     !real(8),    parameter :: geosync_alt   =42164.0e3      !satellite orbits, montenbruck and gill, page 294
    !     !!      integer(4), parameter :: nsat_rfi=8
    !     !!      real(8),    parameter :: geosync_lon(nsat_rfi)=(/13., 19.2, 28.2, 317.0, 352.8, 10.0, 257.2, 260.8/)

    !     integer(4), parameter :: num_extra_scans = 1
    !     integer(4), parameter :: num_extra_fovs = 1
    !     integer(4), parameter :: maxscan_extra = maxscan*(1+num_extra_scans)
    !     integer(4), parameter :: maxcel_extra  = maxcel*(1+num_extra_fovs) - num_extra_fovs

    !     integer(4), parameter :: ilu_l2b_extra  =  7

    !     integer(4), parameter :: ncel_rsp = maxcel_extra  !number of cells for resampling program    

    !     integer(4), parameter ::  iradius_fieldavg=35  !see 'O:\ssmis\L2_processing\memo9.txt' for discussion of search radius

    !     character(120) filename_tb_native,filename_tb_resampled
    !     character(120) filename_ta_resample_wts
    !     character( 60) pathname_tb_native,pathname_tb_resampled

    !     real(8) ut1_1993_extra(maxscan_extra) !really tai_1993. it is a placeholder for heritage code
    !     integer(4) iscan_flag_extra(maxscan_extra)

    !     ! geolocation info -- need to keep track of these when geolocating the extra locations
    !     real(4) omega_extra(maxscan_extra),scan_index_extra(maxscan_extra)
    !     real(8) scpos_eci0_extra(3,maxscan_extra), x_eci0_extra(3,maxscan_extra), y_eci0_extra(3,maxscan_extra), z_eci0_extra(3,maxscan_extra)
    !     real(8) zang_extra(maxscan_extra),scan_time_extra(maxscan_extra)

    !     !not sure if we need to reconsider all the sun and moon things....
    !     integer(1) kflag_sun_extra(maxscan_extra)
    !     real(8) sunvec0_extra(3,maxscan_extra),moonvec0_extra(3,maxscan_extra),alpha_sun_extra(maxscan_extra),beta_sun_extra(maxscan_extra)
    !     real(4) sunlat_extra(maxscan_extra),sunlon_extra(maxscan_extra),sundis_extra(maxscan_extra)

    !     real(4) scloc_extra(3,maxscan_extra)
    !     real(4) cellat_extra(maxcel_extra,maxscan_extra),cellon_extra(maxcel_extra,maxscan_extra)
    !     real(4) celrange_extra(maxcel_extra,maxscan_extra)
    !     real(4) celphi_extra(maxcel_extra,maxscan_extra),celtht_extra(maxcel_extra,maxscan_extra)
    !     real(4) celsun_extra(maxcel_extra,maxscan_extra)

    !     !we'll use the rfi flag from the tb file
    !     real(4) reflat_extra(maxcel_extra,maxscan_extra),reflon_extra(maxcel_extra,maxscan_extra) !geostationary lat and lon of reflected ray

    !     !not doing 89 right now....
    !     !real(4) cellat_89(maxcel_89,maxscan_89),cellon_89(maxcel_89,maxscan_89),cellat_cold(maxscan),cellon_cold(maxscan)
    !     !real(4) ta_cold(nch,maxscan),ta_hot(nch,maxscan),ta_spill(nch,maxscan),gainest(nch,maxscan)
        

    !     ! this is where the resampled TBs go
    !     real(4) tbr_extra(maxcel_extra,maxscan_extra,12) !,ta_89(maxcel_89,maxscan_89,2)

    !     !   ocean  block
    !     integer(1) isurcel_extra(5,maxcel_extra,maxscan_extra)
    !     integer(1) ice_flag_extra(3,maxcel_extra,maxscan_extra)
    !     integer(1) ice_flag2_extra(maxcel_extra,maxscan_extra)
    !     integer(1) ioob_flag_extra(5,maxcel_extra,maxscan_extra)

    !     real(4) ta_extra(maxcel_extra,maxscan_extra,12)

    ! end module l2_module_extra
    
 module l2_module_extra
    use sat_id_module
    implicit none

    ! Parameters remain unchanged
    integer(4), parameter :: num_extra_scans = 1
    integer(4), parameter :: num_extra_fovs = 1
    integer(4), parameter :: maxscan_extra = maxscan*(1+num_extra_scans)
    integer(4), parameter :: maxcel_extra  = maxcel*(1+num_extra_fovs) - num_extra_fovs
    integer(4), parameter :: ilu_l2b_extra  =  7
    integer(4), parameter :: ncel_rsp = maxcel_extra
    integer(4), parameter :: iradius_fieldavg = 35

    ! File paths and names remain as regular variables
    character(120) :: filename_tb_native, filename_tb_resampled
    character(120) :: filename_ta_resample_wts
    character(60)  :: pathname_tb_native, pathname_tb_resampled

    ! Convert fixed arrays to allocatable arrays
    real(8), allocatable :: ut1_1993_extra(:)
    integer(4), allocatable :: iscan_flag_extra(:)

    ! Geolocation arrays
    real(4), allocatable :: omega_extra(:), scan_index_extra(:)
    real(8), allocatable :: scpos_eci0_extra(:,:)
    real(8), allocatable :: x_eci0_extra(:,:), y_eci0_extra(:,:), z_eci0_extra(:,:)
    real(8), allocatable :: zang_extra(:), scan_time_extra(:)

    ! Sun and moon related arrays
    integer(1), allocatable :: kflag_sun_extra(:)
    real(8), allocatable :: sunvec0_extra(:,:), moonvec0_extra(:,:)
    real(8), allocatable :: alpha_sun_extra(:), beta_sun_extra(:)
    real(4), allocatable :: sunlat_extra(:), sunlon_extra(:), sundis_extra(:)

    ! Cell location and geometry arrays
    real(4), allocatable :: scloc_extra(:,:)
    real(4), allocatable :: cellat_extra(:,:), cellon_extra(:,:)
    real(4), allocatable :: celrange_extra(:,:)
    real(4), allocatable :: celphi_extra(:,:), celtht_extra(:,:)
    real(4), allocatable :: celsun_extra(:,:)

    ! Reflection arrays
    real(4), allocatable :: reflat_extra(:,:), reflon_extra(:,:)

    ! Temperature arrays
    real(4), allocatable :: tbr_extra(:,:,:)

    ! Ocean block arrays
    integer(1), allocatable :: isurcel_extra(:,:,:)
    integer(1), allocatable :: ice_flag_extra(:,:,:)
    integer(1), allocatable :: ice_flag2_extra(:,:)
    integer(1), allocatable :: ioob_flag_extra(:,:,:)
    real(4), allocatable :: ta_extra(:,:,:)

contains
    subroutine initialize_l2_extra_arrays()
        implicit none
        integer :: alloc_stat

        ! Allocate 1D arrays
        allocate(ut1_1993_extra(maxscan_extra), &
                 iscan_flag_extra(maxscan_extra), &
                 omega_extra(maxscan_extra), &
                 scan_index_extra(maxscan_extra), &
                 zang_extra(maxscan_extra), &
                 scan_time_extra(maxscan_extra), &
                 kflag_sun_extra(maxscan_extra), &
                 alpha_sun_extra(maxscan_extra), &
                 beta_sun_extra(maxscan_extra), &
                 sunlat_extra(maxscan_extra), &
                 sunlon_extra(maxscan_extra), &
                 sundis_extra(maxscan_extra), &
                 stat=alloc_stat)
        if (alloc_stat /= 0) then
            write(*, *) "Error allocating 1D arrays. Status: ", alloc_stat
            error stop
        endif

        ! Allocate 2D arrays
        allocate(scpos_eci0_extra(3,maxscan_extra), &
                 x_eci0_extra(3,maxscan_extra), &
                 y_eci0_extra(3,maxscan_extra), &
                 z_eci0_extra(3,maxscan_extra), &
                 sunvec0_extra(3,maxscan_extra), &
                 moonvec0_extra(3,maxscan_extra), &
                 scloc_extra(3,maxscan_extra), &
                 cellat_extra(maxcel_extra,maxscan_extra), &
                 cellon_extra(maxcel_extra,maxscan_extra), &
                 celrange_extra(maxcel_extra,maxscan_extra), &
                 celphi_extra(maxcel_extra,maxscan_extra), &
                 celtht_extra(maxcel_extra,maxscan_extra), &
                 celsun_extra(maxcel_extra,maxscan_extra), &
                 reflat_extra(maxcel_extra,maxscan_extra), &
                 reflon_extra(maxcel_extra,maxscan_extra), &
                 ice_flag2_extra(maxcel_extra,maxscan_extra), &
                 stat=alloc_stat)
        if (alloc_stat /= 0) then
            write(*, *) "Error allocating 2D arrays. Status: ", alloc_stat
            error stop
        endif

        ! Allocate 3D arrays
        allocate(tbr_extra(maxcel_extra,maxscan_extra,14), &
                 isurcel_extra(5,maxcel_extra,maxscan_extra), &
                 ice_flag_extra(3,maxcel_extra,maxscan_extra), &
                 ioob_flag_extra(5,maxcel_extra,maxscan_extra), &
                 ta_extra(maxcel_extra,maxscan_extra,14), &
                 stat=alloc_stat)
        if (alloc_stat /= 0) then
            write(*, *) "Error allocating 3D arrays. Status: ", alloc_stat
            error stop
        endif

        ! Initialize arrays to zero/default values
        ut1_1993_extra = 0.0d0
        iscan_flag_extra = 0
        omega_extra = 0.0
        scan_index_extra = 0.0
        scpos_eci0_extra = 0.0d0
        x_eci0_extra = 0.0d0
        y_eci0_extra = 0.0d0
        z_eci0_extra = 0.0d0
        zang_extra = 0.0d0
        scan_time_extra = 0.0d0
        kflag_sun_extra = 0
        sunvec0_extra = 0.0d0
        moonvec0_extra = 0.0d0
        alpha_sun_extra = 0.0d0
        beta_sun_extra = 0.0d0
        sunlat_extra = 0.0
        sunlon_extra = 0.0
        sundis_extra = 0.0
        scloc_extra = 0.0
        cellat_extra = 0.0
        cellon_extra = 0.0
        celrange_extra = 0.0
        celphi_extra = 0.0
        celtht_extra = 0.0
        celsun_extra = 0.0
        reflat_extra = 0.0
        reflon_extra = 0.0
        tbr_extra = 0.0
        isurcel_extra = 0
        ice_flag_extra = 0
        ice_flag2_extra = 0
        ioob_flag_extra = 0
        ta_extra = 0.0
    end subroutine initialize_l2_extra_arrays

    subroutine cleanup_l2_extra_arrays()
        implicit none
        
        if (allocated(ut1_1993_extra)) deallocate(ut1_1993_extra)
        if (allocated(iscan_flag_extra)) deallocate(iscan_flag_extra)
        if (allocated(omega_extra)) deallocate(omega_extra)
        if (allocated(scan_index_extra)) deallocate(scan_index_extra)
        if (allocated(scpos_eci0_extra)) deallocate(scpos_eci0_extra)
        if (allocated(x_eci0_extra)) deallocate(x_eci0_extra)
        if (allocated(y_eci0_extra)) deallocate(y_eci0_extra)
        if (allocated(z_eci0_extra)) deallocate(z_eci0_extra)
        if (allocated(zang_extra)) deallocate(zang_extra)
        if (allocated(scan_time_extra)) deallocate(scan_time_extra)
        if (allocated(kflag_sun_extra)) deallocate(kflag_sun_extra)
        if (allocated(sunvec0_extra)) deallocate(sunvec0_extra)
        if (allocated(moonvec0_extra)) deallocate(moonvec0_extra)
        if (allocated(alpha_sun_extra)) deallocate(alpha_sun_extra)
        if (allocated(beta_sun_extra)) deallocate(beta_sun_extra)
        if (allocated(sunlat_extra)) deallocate(sunlat_extra)
        if (allocated(sunlon_extra)) deallocate(sunlon_extra)
        if (allocated(sundis_extra)) deallocate(sundis_extra)
        if (allocated(scloc_extra)) deallocate(scloc_extra)
        if (allocated(cellat_extra)) deallocate(cellat_extra)
        if (allocated(cellon_extra)) deallocate(cellon_extra)
        if (allocated(celrange_extra)) deallocate(celrange_extra)
        if (allocated(celphi_extra)) deallocate(celphi_extra)
        if (allocated(celtht_extra)) deallocate(celtht_extra)
        if (allocated(celsun_extra)) deallocate(celsun_extra)
        if (allocated(reflat_extra)) deallocate(reflat_extra)
        if (allocated(reflon_extra)) deallocate(reflon_extra)
        if (allocated(tbr_extra)) deallocate(tbr_extra)
        if (allocated(isurcel_extra)) deallocate(isurcel_extra)
        if (allocated(ice_flag_extra)) deallocate(ice_flag_extra)
        if (allocated(ice_flag2_extra)) deallocate(ice_flag2_extra)
        if (allocated(ioob_flag_extra)) deallocate(ioob_flag_extra)
        if (allocated(ta_extra)) deallocate(ta_extra)
    end subroutine cleanup_l2_extra_arrays

end module l2_module_extra   

