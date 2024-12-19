!     7/13/2021 Convert to module, changed to openbig_try
!     3/26/2017 changed 8/31/2017.  itai_1993_leap is no longer read in
!     it is read in by 'read_rss_l1a.f'.  
!     see 'O:\amsr2\L2_processing\memo7.txt'


!     v8 routine establish on 3/26/2017 when i put in read dta_zang


!     the following are v7 comments

!     11/4/2013 version updated 7/10/2014.  array corr_tab_11vlna is read in

!     9/19/2013 version changed 11/4/2013.  apc delta,chi is now read in

!     9/17/2013 version changed 9/19/2013.  rfi_cold_mirror_map is now read in

!     9/9/2013 version changed 9/17/2013.  moon parameter are now read in

!     8/28/2013 version changed 9/9/2013.  an updated filename_cnt is read in

!     8/22/2013 version changed 8/28/2013. tamap read in  

!     march 13 2013 version change 8/22/2013
!     itai_1993_leap is now read in
!     tbmin and tbmax read in
!     percent_land_tab, amask_ice and isal_corr now read in


!     march 6 2013 version changed march 13 2013.  xscan table read here rather than in fd_ta_adjusted

!     feb 8 2013 version changed march 6 2013. acoef_nl read in

!     reference xscan_cal to middle cell, ie make middle cell zero
!     use same stddv_cal for all cells    

module DataFiles 

contains
    subroutine readin_data_files 
    use open_file_routines, only: openbig_try                                           
    use l2_module                                                      
    implicit none                                               
 
    real(8) period0
    integer(4) ibad,nfreq_in,ifreq

!     ===============================================================================================================
!     ============================== read in table of tai leap seconds ==============================================   
!     ===============================================================================================================
!     call openbig_try(3,filename_leap,'old') 
!     read(3) itai_1993_leap
!     close(3)
    
!     ===============================================================================================================
!     ======================================= read in  bad orbit data ===============================================   
!     ===============================================================================================================
    open(3,file=filename_bad,status='old')    
    read(3,'(i5)') nbad
    if(nbad.gt.10000) stop 'pgm stopped, nbad.gt.10000'
    do ibad=1,nbad
        read(3,'(i5)') iorbit_bad(ibad)
    enddo
    close(3)
    !if (nbad .gt. 0) then
    !    write(*,*) 'nbad = ',nbad
    !endif
    

!   ===============================================================================================================
!   ======================================= read in geolocatoin data ==============================================   
!   ===============================================================================================================
    open(3,file=filename_geo,status='old')
    read(3,3001) period0,eia_convert_37to19
    read(3,3002) beta
    read(3,3002) azoffset
    read(3,3003) start_time_a,start_time_b,spillover_time
    read(3,3002) rolloffset,rpy_error
    3001 format(f10.2,f10.6)
    3002 format(10f10.3)
    3003 format(10f10.7)
    close(3)

    
!     ===============================================================================================================
!     ============================== read in cold count filtering data ==============================================   
!     ===============================================================================================================
    call openbig_try(3,filename_cnt,'old') 
    read(3) xscan_cal
    read(3) stddv_avg
    close(3)
    
!     ===============================================================================================================
!     ======================================= read in  non-linearity coefs ==========================================   
!     ===============================================================================================================
    call openbig_try(3,filename_nl,'old') 
    read(3) acoef_nl
    close(3)
      
!     ===============================================================================================================
!     ===================================== read in along-scan correction table =====================================   
!     ===============================================================================================================
    call openbig_try(3,filename_xsc,'old')
    read(3) xscan_tab
    close(3)
    
    
!     ===============================================================================================================
!     ============================================= read in tb limits  ==============================================   
!     ===============================================================================================================
    open(3,file=filename_tbl,status='old',action='read')
    read(3,'(i9)') nfreq_in
    if(nfreq_in.gt.nfreq) stop 'nfreq oob in ckconta, pgm stopped'
    read(3,'(16f9.3)') tbmin
    read(3,'(16f9.3)') tbmax
    close(3)

!     ===============================================================================================================
!     =============================== read in land percent map and climate ice mask =================================   
!     ===============================================================================================================
    open(3,file=filename_percent_land,status='old',form='binary')
    read(3) percent_land_tab
    close(3)

    open(3,file=ice_filename,status='old',form='binary')
    read(3) amask_ice
    close(3)


!   ===============================================================================================================
!   ======================================== read sss correction table ============================================   
!   ===============================================================================================================
    open(3,file=sss_filename,status='old')
    read(3,3006) isal_corr
    3006 format(30i4)
    close(3)
      
!     ===============================================================================================================
!     ======================================== read tamaps for spillover correction  ================================   
!     ===============================================================================================================
    call openbig_try(3,filename_tamap,'old')  !contains maps for first 6 channels
    read(3) tamap
    close(3)
      

!     ===============================================================================================================
!     ========================================== read moon parameters ===============================================   
!     ===============================================================================================================
    open(3,file=filename_moon,status='old')
    read(3,3011) amp_moon
    read(3,3012) phase_factor_moon
!     the ifreq do loops are required so that i can use the same format for both amsre and amsrj
    do ifreq=1,nfreq
        read(3,3013) bw_cold(:,ifreq)
    enddo
    do ifreq=1,nfreq
        read(3,3013) alpha_cold(:,ifreq) 
    enddo
    do ifreq=1,nfreq
        read(3,3013)  beta_cold(:,ifreq)
    enddo
    3011 format( 8f8.3)
    3012 format( 8f8.5)
    3013 format(32f8.3)
    close(3)
    
!   ===============================================================================================================
!   =============================== read rfi map for residual cold mirror contamination ===========================   
!   ===============================================================================================================
    call openbig_try(3,filename_crfi,'old') 
    read(3) rfi_cold_mirror_map
    close(3)
    
!   ===============================================================================================================
!   ================================================ read apc =====================================================   
!   ===============================================================================================================
    open(3,file=filename_apc,status='old')
    read(3,3014) delta_apc
    read(3,3014) chi_apc
    3014 format(16f10.6)
    close(3)
      
!   ===============================================================================================================
!   ============================= read correction table for 11V lna probleem =====================================   
!   ===============================================================================================================
    call openbig_try(3,filename_lna,'old')
    read(3) corr_tab_11vlna
    close(3)

!   ===============================================================================================================
!   ============================= read correction table dta versus zang ==========================================   
!   ===============================================================================================================
    call openbig_try(3,filename_zang,'old')
    read(3) dta_zang
    close(3)
    return
    end subroutine readin_data_files

end module DataFiles
        