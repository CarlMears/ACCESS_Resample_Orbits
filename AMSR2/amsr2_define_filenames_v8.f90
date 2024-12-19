!     3/26/2017 changed 8/31/2017.  filename_leap changed to 'x:\time\ut1.dat'


!     3/7/2017 changed 3/26/2017.  filename_xsc updated to v8.  filename_zang added

!     2/28/2017 changed 3/7/2017. filename_crfi updated to v8 version

!     the v8 version is in the beginning process of getting update

!     7/11/2014 version changed 7/17/2014.  filename_2be and filename_inv added

!     9/19/2013 version updated 7/11/2014.  filename_lna added

!     9/17/2013 version changed 9/19/2013.  filename_crfi added

!     9/9/2013 version changed 9/17/2013.  filename_moon added

!     8/28/2013 version changed 9/9/2013 name of filename_cnt changed to name of updated file

!     8/22/2013 versin changed 8/28/2013  filename_tamap replaces filename_cband_map; filename_moon removed

!     march 13 2013 version changed 8/22/2013.
!     filename_ut1 replaced by filename_leap
!     four filenames are hardwired rather than being read from 'x:\geomod_data_files\ocean_input_files.txt'
!     filename_percent_land added


!     march 6 2013 version changed march 13 2013.  xscan filename has changed

!     feb 8 2013 version changed march 6 2013, filename_nl added

!     this routine is now updated to the v7 standard
    subroutine define_filenames
        use l2_module                                                      
        implicit none    
    
        character(2) asat

!   filename_leap='x:\time\tai_leap_table.dat'      
        filename_leap='x:\time\ut1.dat'      
        
        filename_l1a_inventory=      'j:\amsr2\tables\l1a_inventory.dat'
        filename_rss_l1a_inventory = 'j:\amsr2\tables\rss_l1a_inventory.dat'
        filename_node=               'j:\amsr2\tables\node.dat'
        filename_orb=                'j:\amsr2\tables\orbit_times.dat'
        filename_bad=                'j:\amsr2\tables\bad_orbits.txt'
            
        pathname_l1a=        'j:\amsr2\l1a\' 
        pathname_l2b=        'j:\amsr2\l2b_v08\' 
    
        !     these are all v7 files
        filename_ta_resample_wts=  'o:\amsr2\resampling\make_processing_weights.dat'
        filename_ta89_resample_wts='o:\amsr_l2\v05\resampling\make_processing_weights_hi.dat'
        filename_tamap=            'o:\amsr2\cold_counts\spillover_tamaps.dat'
        filename_nl=               'o:\amsr2\abs_cal\acoef_nonlinear_v8.dat'
        filename_zang=             'o:\amsr2\abs_cal\smooth_zang.dat'
    
    !    filename_tef=              'o:\amsr\abs_cal\fd_dteff_corr_new_f32_yescorr.dt2'

        write(asat,9000) ksat
        9000 format(i2.2)
        filename_tbl='o:\algo99\mk_algo\mk_tables_v8\data\fd_tb_range_f'//asat//'.lis'  
        filename_apc='j:\amsr2\tables\apc_coefs_v8_f'//asat//'.txt'
        filename_geo='j:\amsr2\tables\geo_coefs_f'//asat//'.txt'
        filename_cnt='j:\amsr2\tables\average_xscan_std_stats.dat'
        filename_xsc='j:\amsr2\tables\xscan_table_v8_f'//asat//'.dat'             

        land_filename=   'x:\land\land_mask_01.txt'
        ice_filename=    'x:\land\ice3.dat'
        climate_filename='x:\land\mergwin_vap_v2.txt'
        sss_filename=    'x:\geomod_data_files\mk_tb_salinity_climate.txt'  ! o:\salinity\mk_tb_salinity_climate.txt

        filename_percent_land =    'x:\land\land_tables_amsra.dat'
        
        filename_moon = 'o:\amsr2\moon_in_cold_mirror\mk_moon_tables_amsr2.lis'
        filename_crfi = 'o:\amsr2\cold_counts\data\mk_cold_mirror_rfi_maps_v8.dat'
        filename_lna=   'o:\amsr2\abs_cal\11ghzlna\fd_11vlna_corr.dat'
        filename_2be=   'j:\amsr2\tables\l2b_files_tobe_processed.dat'
        filename_inv=   'j:\amsr2\tables\l2b_inventory.dat'
        return
    end subroutine define_filenames
