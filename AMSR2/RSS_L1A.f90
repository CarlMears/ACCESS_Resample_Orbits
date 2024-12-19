!         7/13/2021 Converted to module

!         8/31/2017 changed 11/27/2017.
!         when calling find_therm, i now use the start orbit orbit(1) rather than orbit(iscan)
!         but this is only done for orbit.ge.29392 so that previous processing is not affect when reprocessing is done
!         this change was in response to another jaxa change in cspc.


!         12/4/2015 changed 8/31/2017.  problem with leap seconds.
!         it was not getting updated.  last update was the 2012 183  7  1 update
!         i had to reprocess starting at orbit 16591 (june 30 2015) using this update
!         see 'O:\amsr2\L2_processing\memo7.txt'


!         12/5/2013 updated 12/4/2015
!         gap error message now writes iscan_last and gap tolerance change from 0.5 to 0.5 + 2*abs((iscan-iscan_last)/4000.)

!         9/19/2013 version changed 12/5/2013.  argument orbit(iscan) added of find_therm

!         9/9/2013  version changed 9/19/2013  setup for navy

!         8/28/2013 versions changed 9/9/2013.  routine now checks for time steps to be correct.
!         test run was done for all orbits and so far there is no problem with time steps.

!         8/19/2013 version modified 8/29/2013. nleap is now found by a separate routine: fd_leap_sec
!         jscan1 is now 2*iscan-3 to make 89b adjacent and before to 89a

!         may 3 2013 version changed 8/19/2013
!         ut1 code was removed.  i dont need it because jaxa provide earth-fixed sc location
!         timeofscan and leap_sec not used are were removed
!         ut1_1993 simply set to tai_1993 (see below for explanation)
!         nleap found in a more simple way

!         april 25 2013 version changed on may 3 2013.  array earthcounts_jaxa added.
!         now the earth counts i recover from jaxa l1b files, which have the suffix _rss,
!         are stored in earthcounts and the l1a jaxa values are stored in earthcounts_jaxa.
!         earthcounts_jaxa is not used for anything other than may some diag or qc later

!         mar 8 2013 version change april 25 2013. irx_temp_cnt array added

!         march 1 2013 version changed on mar 8 2013:  major changes to update to v1 formats

!         feb 11 2013 version change mar 1 2013,  Computation of scvel is now done correctly.

!         1/12/2010 version changed on 10/8/2010.  ut1 info now comes from a static file.
!         i did not want v7 processing to interfere with the standard processing.
!         i also wanted to be able to run mnay l2b programs at once so i did
!         not want to use the standard routines that update the table as time goes on.
!         this routine will let you know when you need a fresh copy of the ut1 files    

!         this routine is dif from  'o:\amsr_l2\v05\routines\read_amsr_l0.f'     in the following respects
!         1. start_time is relative to jan 1 2000      rather than 1993 
!         2. variable scan_time (relaive to jan 1 2000) is added so i now have two scan times reference to 2000 and 1993.
!         3. call new get_l1a_filename

module RSS_L1A 

    use sat_id_module
    implicit none

    integer(4), parameter :: reclen =39816 !you can determine this by writing out L1arss and looking at it size
    integer(4), parameter :: npad=reclen-16

    character(reclen) apar(0:maxscan)
    character(reclen) abuf1
    character(reclen) abuf2
     
    !     ==========================================================================
    !     ======================= header structure =================================
    !     ==========================================================================
    type level1a_rss_header
        sequence
        real(8) scan_time
        integer(4) iorbit
        integer(4) numscan
        character(npad) apad
    end type level1a_rss_header
    

    !     ==========================================================================
    !     ======================= rss l1a data structure ===========================
    !     ==========================================================================

    type level1a_rss
        sequence
        
        real(8) scan_time
        real(8) rev
        
        real(4) omega
        real(4) navigation(6)
        real(4) attitude(3)

        integer(4) iscan
        integer(4) spc_temp_count(spc_length)
        integer(4) sps_temp_count(sps_length)
        integer(4) rx_offset_gain(2,nch)
        integer(4) ipad  !this makes the structure a multiple of 8

        real(4) lat89a(maxcel_89)
        real(4) lat89b(maxcel_89)
        real(4) lon89a(maxcel_89)
        real(4) lon89b(maxcel_89)
        real(4) thot(nch)
        
        integer(2) loEarthCounts_rss(maxcel,   l1a_lo_chan)
        integer(2) hiEarthCounts_rss(maxcel_89,l1a_hi_chan)

        integer(2) loEarthCounts(maxcel,   l1a_lo_chan)
        integer(2) hiEarthCounts(maxcel_89,l1a_hi_chan)
        integer(2) cld_lo_counts(l1a_lo_chan, low_cal_obs)
        integer(2) hot_lo_counts(l1a_lo_chan, low_cal_obs)
        integer(2) cld_hi_counts(l1a_hi_chan,high_cal_obs)
        integer(2) hot_hi_counts(l1a_hi_chan,high_cal_obs) 
        integer(2) loTb(maxcel,   l1a_lo_chan)  
        integer(2) hiTb(maxcel_89,l1a_hi_chan)  
        
        integer(2) interp_flag(32,nch)
    
    end type level1a_rss

    type (level1a_rss_header) :: L1arss_hd
    type (level1a_rss) :: L1arss

     
    equivalence(L1arss_hd ,abuf1)
    equivalence(L1arss,    abuf2)

    public level1a_rss_header, level1a_rss ! types
    public L1arss_hd,L1arss                ! module variables

    !these are subroutines
    public read_rss_l1a,convert_earthfix_to_2000,get_sc_axes,calcrollpitchyaw,read_binary_head

contains

    subroutine read_rss_l1a(iorbit,iquick, iexist,start_time)
        !use amsr_rss_l1a_input_module - now in this module                                                                  
        use l2_module    
        use io_filenames, only: get_amsr2_rss_l1a_filename 
        use open_file_routines, only: open_direct       
        use pre_nut_routines, only: get_gm_angle, fd_precession, fd_nutation  
        use math_routines, only: cross_norm,dot_product_unit8,invert_3by3,fixang4,fixang8,minang4,minang8                                                    
        implicit none

        real(8), parameter :: d_time=1.5d0  !arbitary, could have be 1 sec
        real(8), parameter :: scan_period=1.5d0 

        character(100) filename_rss_l1a 
        real(8) xtime_first,start_time,time_last
        real(8) scpos0(3),scvel0(3)
        real(8) ut1_2000,rotangle,days,np(3,3),np_inv(3,3),det 
        real(8) scposx(3),rotanglex,deltim
        
        real(4) gap_tol
        
        integer(4) itype,i,iscan_last
        integer(4) iscan
        integer(4) icel,ich,iload
        integer(4) iorbit,iquick,iexist
        integer(4) itempy
        integer(4) jscan1,jscan2
        integer(4) nleap
        integer(4) irec
        
        call get_amsr2_rss_l1a_filename(iorbit, filename_rss_l1a,iexist) 
        if(iexist.eq.0) return
        !print *,'Reading: ',filename_rss_l1a

        call read_binary_head(iorbit,filename_rss_l1a, numscan,xtime_first)
        !print *,'Num Scans in L1A file = ',numscan
        !========================================================================================================
        !======================= new code to update leap sec inserted 8/31/2017 =================================     
        !========================================================================================================
        call open_direct(3, filename_leap, 'old', 'read', 'denywr', 4, 10)   !10  is 10 seconds of trying
        do irec=1,100
            read(3,rec=irec) itai_1993_leap(irec)
            if(itai_1993_leap(irec).eq.-2147483648) itai_1993_leap(irec)=2147483647
        enddo
        close(3)


        iflag_rx=0
        ta_89=0
        ta_89=0
        iscan_last=-1
        time_last=-1

        iflag_l0 = 10
        do iscan=1,numscan
            jscan1=2*iscan-3  !89b assign to this to make it adjacent and before to 89a
            jscan2=2*iscan    
            abuf2=apar(iscan)
            if(L1arss.iscan.ne.iscan) stop 'pgm stopped, iscan out of sync'
            if(iscan.eq.1 .and. l1arss.scan_time.eq.0) stop 'problem with first scan, pgm stopped'

            iflag_l0(iscan)=0
            if(l1arss.scan_time.eq.0) then  !missing scans
                iflag_l0(iscan)=10      !no science data
                scan_time(iscan)=-1.d30
                ut1_1993(iscan)=0
                orbit(iscan)=0
                omega(iscan)=0
                scpos(:,iscan)=0
                scvel(:,iscan)=0
                scrpy(:,iscan)=0
                therm(:,iscan)=0
                rx_offset_gain(:,:,iscan)=0
                earthcounts(   :,:,iscan)=0
                earthcounts_jaxa(   :,:,iscan)=0
                csmcounts(     :,:,iscan)=0
                htscounts(     :,:,iscan)=0
                ta(:,iscan,:)=0
                cycle
            endif
        
            if(iscan_last.ne.-1) then
                deltim=l1arss.scan_time-time_last
                gap_tol=0.5 + 2.*abs((iscan-iscan_last)/4000.)
                if(abs((iscan-iscan_last)*scan_period-deltim).gt.gap_tol) then
                    write(*,*) iscan_last,iscan,abs((iscan-iscan_last)*scan_period-deltim)
                    stop 'scan gap error, pgm stopped'
                endif
            endif
        
            iscan_last=iscan
            time_last=l1arss.scan_time

            call fd_leap_sec(L1arss.scan_time, nleap)
        
            scan_time(iscan)=L1arss.scan_time - nleap  - 220838400.d0 !convert from tai 93 to utc begin 2000
            if(iscan.eq.1) start_time=scan_time(iscan)
        
            orbit(iscan)=L1arss.rev
            omega(iscan)=L1arss.omega
        
                !         ut1_1993 is actually tai_1993 but for this program it does not matter because
                !         ut1_1993 is just used to convert jaxa fixed earth refernced to j2000 and then it
                !         is used to convert it back to fixed earth reference in amsr geolocation.
                !         i verfied its value does not matter because of the circular usage.
                !         in fact, you could simple use ut1_1993(iscan)=0, and the program works fine.
                !         i keep the name ut1_1993 to reference back to legacy code.
            ut1_1993(iscan)=L1arss.scan_time

            if(maxval(abs(L1arss.navigation(1:3))).eq.0) iflag_l0(iscan)=9      !no geolocation info

            if(iquick.eq.1) cycle

            if(iflag_l0(iscan).ne.0) then
                scpos(:,iscan)=0
                scvel(:,iscan)=0
                scrpy(:,iscan)=0
            else
                scpos(:,iscan)=L1arss.navigation(1:3) 
                scvel(:,iscan)=L1arss.navigation(4:6)  
                scrpy(:,iscan)=L1arss.attitude  
            endif

            !         =====================================================================================================
            !         =============== begin special block for amsr2, which uses earth fixed reference system ==============     
            !         =====================================================================================================
            ut1_2000= ut1_1993(iscan)  - 220838400.d0 !convert from begin 93 to begin 2000
            days=ut1_2000/86400.d0 -0.5d0
            call get_gm_angle(ut1_2000,days, rotangle)
            
            call fd_precession(days, np)  
            call invert_3by3(np, np_inv,det)
            if(abs(det-1).gt.1.e-6) stop 'error in np_inv,pgm stopped'

            ut1_2000= ut1_2000 + d_time !convert from begin 93 to begin 2000
            days=ut1_2000/86400.d0 -0.5d0
            call get_gm_angle(ut1_2000,days, rotanglex)

            scpos0=scpos(:,iscan)
            scvel0=scvel(:,iscan)
            scposx=scpos0 + d_time*scvel0
        
            call convert_earthfix_to_2000(rotangle,np_inv, scpos0) 
            call convert_earthfix_to_2000(rotanglex,np_inv, scposx)
            scvel0=(scposx-scpos0)/d_time

            scpos(:,iscan)=scpos0
            scvel(:,iscan)=scvel0

            !         =====================================================================================================
            !         ================  end special block for amsr2, which uses earth fixed reference system ==============     
            !         =====================================================================================================
                
            !         the following check was done during v0 processing from July 2012 to Jan 2013
            !         no error was detected
            !         i have now commented this out, so that it would not stop the program if future data is oob
            !     do i=1,sps_length
            !     if(L1arss.sps_temp_count(i).eq.65535) cycle  !valid error
            !     if(L1arss.sps_temp_count(i).lt.0 .or. L1arss.sps_temp_count(i).gt.8191) stop 'sps count oob, pgm stopped'
            !     enddo
        
            do itype=1,10
                i=6+itype  !position in spc array
                if(orbit(1).lt.29392) then
                    call find_therm(orbit(iscan),itype,L1arss.spc_temp_count(i), therm(itype,iscan))
                else
                    call find_therm(orbit(1),    itype,L1arss.spc_temp_count(i), therm(itype,iscan))
                endif
            enddo

            irx_temp_cnt(:,iscan)=L1arss.sps_temp_count(33:40)

            do ich=1,nch
                do iload=1,2
                    rx_offset_gain(ich,iload,iscan)=L1arss.rx_offset_gain(iload,ich)
                    if(rx_offset_gain(ich,iload,iscan).le.0 .or. rx_offset_gain(ich,iload,iscan).ge.255) then
                        iflag_rx(ich,iscan) = ibset(iflag_rx(ich,iscan),1)
                    endif
                enddo
            enddo

            do icel=1,maxcel
                earthcounts(     1:12,icel,iscan)=L1arss.loEarthCounts(icel,:)
                earthcounts_jaxa(1:12,icel,iscan)=L1arss.loEarthCounts_rss(icel,:)
            enddo
        
            do icel=1,maxcel_89
                earthcounts(     13:16,icel,iscan) =L1arss.hiEarthCounts(icel,:)
                earthcounts_jaxa(13:16,icel,iscan) =L1arss.hiEarthCounts_rss(icel,:)
            enddo
        
            do icel=1,low_cal_obs
                csmcounts(1:12,icel,iscan)=L1arss.cld_lo_counts(:,icel) 
                htscounts(1:12,icel,iscan)=L1arss.hot_lo_counts(:,icel) 
            enddo
        
            do icel=1,high_cal_obs
                csmcounts(13:16,icel,iscan)=L1arss.cld_hi_counts(:,icel) 
                htscounts(13:16,icel,iscan)=L1arss.hot_hi_counts(:,icel) 
            enddo
        
            do icel=1,maxcel
                do ich=1,12
                    itempy=L1arss.loTb(icel,ich)
                    if(itempy.lt.0) itempy=itempy+65536
                    ta(icel,iscan,ich)=0.01*itempy
                enddo  !ich
            enddo  !icel
        
            do icel=1,maxcel_89
                do ich=1,2
                    itempy=L1arss.hiTb(icel,ich)
                    if(itempy.lt.0) itempy=itempy+65536
                    ta_89(icel,jscan2,ich)=0.01*itempy
                    itempy=L1arss.hiTb(icel,ich+2)
                    if(itempy.lt.0) itempy=itempy+65536
                    if(jscan1.ge.1) ta_89(icel,jscan1,ich)=0.01*itempy
                enddo  !ich
            enddo  !icel
        
            geoloc_jaxa(1,:,iscan)=L1arss.lat89a
            geoloc_jaxa(2,:,iscan)=L1arss.lat89b
            geoloc_jaxa(3,:,iscan)=L1arss.lon89a
            geoloc_jaxa(4,:,iscan)=L1arss.lon89b
        
            thot_jaxa(:,iscan)=L1arss.thot
            interp_flag(:,:,iscan)=L1arss.interp_flag
        enddo      !iscan

        ! =       =======================================================================================================
        ! =       =================================== begin block to find iflag_rx ======================================
        ! =       ============================ note bit 1 could have been set above =====================================
        ! =       =======================================================================================================
     
        do iscan=1,numscan
            do ich=1,16
            
                if(iflag_l0(iscan).eq.10) then
                    iflag_rx(ich,iscan)=ibset(iflag_rx(ich,iscan),2) !no science data
                    cycle
                endif
                
                if(iscan.eq.1 .or. iscan.eq.numscan) then       !exclude end scan because i dont know about changing rx
                    iflag_rx(ich,iscan)=ibset(iflag_rx(ich,iscan),0) 
                    cycle
                endif
                
                if(rx_offset_gain(ich, 1,iscan).ne.rx_offset_gain(ich, 1,iscan-1)) iflag_rx(ich,iscan)=ibset(iflag_rx(ich,iscan),0)
                if(rx_offset_gain(ich, 1,iscan).ne.rx_offset_gain(ich, 1,iscan+1)) iflag_rx(ich,iscan)=ibset(iflag_rx(ich,iscan),0) 
                if(rx_offset_gain(ich, 2,iscan).ne.rx_offset_gain(ich, 2,iscan-1)) iflag_rx(ich,iscan)=ibset(iflag_rx(ich,iscan),0)
                if(rx_offset_gain(ich, 2,iscan).ne.rx_offset_gain(ich, 2,iscan+1)) iflag_rx(ich,iscan)=ibset(iflag_rx(ich,iscan),0) 

            enddo  !ich
        enddo  !iscan
        ! =       =======================================================================================================
        ! =       ===================================  end  block to find iflag_rx ======================================
        ! =       =======================================================================================================
        return
    end subroutine read_rss_l1a

    subroutine convert_earthfix_to_2000(rotangle,np_inv, x)
        implicit none
        real(8) rotangle,np_inv(3,3),x(3)
        real(8) cosr,sinr,y(3)
        
        cosr=cosd(rotangle)
        sinr=sind(rotangle)

        !        xlon=xlon + rotangle     !convert from earth to inertial
        y(1)=cosr*x(1) - sinr*x(2)
        y(2)=cosr*x(2) + sinr*x(1)
        y(3)=x(3)
     
        x=matmul(np_inv,y) !convert eci back to j2000
        return
     end subroutine convert_earthfix_to_2000




    subroutine get_sc_axes(attdq, x,y,z)
        implicit none

        real(8), intent(in)::attdq(4)
        real(8), intent(out):: x(3),y(3),z(3)

        x(1) = (attdq(1)**2+ attdq(4)**2) - (attdq(2)**2 + attdq(3)**2) 
        x(2) = 2.d0 * (attdq(1) * attdq(2) + attdq(3) * attdq(4))
        x(3) = 2.d0 * (attdq(1) * attdq(3) - attdq(2) * attdq(4))
                    
        y(1) = 2.d0 * (attdq(1) * attdq(2) - attdq(3) * attdq(4))
        y(2) = (attdq(2)**2 + attdq(4)**2) - (attdq(1)**2 + attdq(3)**2)
        y(3) = 2.d0 * (attdq(2) * attdq(3) + attdq(1) * attdq(4))
                
        z(1) = 2.d0 * (attdq(1) * attdq(3) + attdq(2) * attdq(4))
        z(2) = 2.d0 * (attdq(2) * attdq(3) - attdq(1) * attdq(4))
        z(3) = (attdq(3)**2 + attdq(4)**2) - (attdq(1)**2 + attdq(2)**2)

        return
    end subroutine get_sc_axes

    subroutine calcrollpitchyaw(scpos,scvel,x,y,z, rpy)
        use math_routines, only: cross_norm,dot_product_unit8,invert_3by3,fixang4,fixang8,minang4,minang8                                                    
        
        implicit none

        !arguments
        real(8), intent(in):: scpos(3), scvel(3), x(3), y(3), z(3)
        real(8), intent(out):: rpy(3)
        
        !local variables
        real(8)::attdmatrix(3,3), x0(3), y0(3), z0(3)
        real(8)::orbtoeci(3,3), rotmat(3,3)

        attdmatrix(1,:) = x 
        attdmatrix(2,:) = y
        attdmatrix(3,:) = z
                    
    !         x0,y0,z0 are s/c axis for 0 roll,pitch,yaw
        z0=-scpos/sqrt(dot_product(scpos,scpos))
        call cross_norm(z0,scvel, y0)
        call cross_norm(y0,z0, x0)
               
        !transformation from orbital to eci frame
        orbtoeci(:,1)=x0
        orbtoeci(:,2)=y0
        orbtoeci(:,3)=z0
               
        rotmat = matmul(attdmatrix, orbtoeci)

    !         following is the 'rpy' or  '123' convention
        rpy(1) = -datan2d(rotmat(3,2),rotmat(3,3))
        rpy(2) = dasind(rotmat(3,1))                
        rpy(3) =  -datan2d(rotmat(2,1),rotmat(1,1))
               
        return
    end subroutine calcrollpitchyaw

    subroutine read_binary_head(iorbit,filename_rss_l1a, numscan,scan_time)
        !use amsr_rss_l1a_input_module                                                      
        implicit none
    
        character(100) filename_rss_l1a   
        real(8) scan_time
        integer(4) iorbit,numscan,isize,ierr,nrec
    
        call fastunzip(filename_rss_l1a, isize,apar,ierr)
        if(ierr.ne.0) stop 'error in fastunzip when reading amsr l0 file'

        nrec=int(isize/reclen)
        if(nrec.le.1 .or. nrec.gt.maxscan+1 .or. reclen*nrec.ne.isize) stop 'pgm stopped, error1 read_binary_header'
        numscan=nrec-1  !subtract header scan

        abuf1=apar(0)
        if(l1arss_hd.iorbit.ne.iorbit .or. l1arss_hd.numscan.ne.numscan) stop 'wrong orbit or numscan, pgm stopped'
        scan_time=l1arss_hd.scan_time
    
        return
    end subroutine read_binary_head
    
end module RSS_L1A 


