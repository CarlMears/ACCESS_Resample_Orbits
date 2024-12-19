!     9/19/2013 version changed on 9/23/2013.  counts_per_deg now is function of channel
!     and is specified in the l2 module.  Also the correction for rfi in the cold counts
!     is now scaled by gainest.


!     9/18/2013 version changed on 9/19/2013.  rfi_cold_count correction implemented

!     i have now impletment reverse_agc which sets oob counts to -32768.  when reverse_agc is done
!     the cal count can fall outside the iminc,imaxc range and this is fine.
!     iflag_moon is now checked
!     spill_adj is now computed

!     aug 28 2013 version changed 9/10/2013.  fixed value for counts_per_deg is no longer used except for a default
!     the actual hot and cold counts for each scan, thot, and tcold are now used to find this.
!     in doing this gainest is a new l2 module array and gainest is also used to convert stddv_avg, 
!     which is now a temperature rms, to a count rms. 
!     ta_spill is no longer found in this routine.  it is found by fd_tcold which must be called first
!    icel .gt. 10 are now excluded for hot counts for channel 9  !see 'O:\amsr2\cold_counts\memo1.txt'
!     old routine is saved with dated added as suffix


!     feb 4 2013 version changed aug 28 2013. i had to change find_tamea_cband to find_ta_spill to
!     accomodate changes to that routine.  new version can handle spillover for all channes, not just 6 and 7 ghz.
!     there was no functional change other than if spillover location is not avaiable ta is set simply to 200
!     but i dont think this has ever happened

!     see 'O:\amsr2\cold_counts\memo1.txt' for a description of this routine and what needs to be done

    subroutine avg_cal_counts
        use l2_module
        use math_routines, only: fd_median
        implicit none                                            

        integer(4), parameter :: iscan_del1=100  !window for cold count filterint
        integer(4), parameter :: iscan_del2=  5  !window for count averaging
        integer(4), parameter :: nval=high_cal_obs*(2*iscan_del1+1)
        
        integer(4) ich,iscan,kscan,icel,kcel,iload,ncel,iradius,istep
        integer(4) i,imid,isum0,isum1

        integer(2)   lflag_tab(0:1,high_cal_obs,maxscan)
    
        integer(4) num,ival(nval),kscansv(nval),icelsv(nval),ivalmid
        
        real(4) xcount,xcount_tab(high_cal_obs,maxscan)
        real(4) ccavg,hcavg
        real(8) xsum1,xsum2,xsum3,xsum4
      
        !     ================================================================================================
        !     ===== make table of bad cold counts: integer(1) lflag_counts(nch,high_cal_obs,maxscan)  ========
        !     ================================================================================================
        lflag_counts=1
          
        do ich=1,nch 
            if(ich.le.12) then
                ncel=low_cal_obs
            else
                ncel=high_cal_obs
            endif
   
            !     =================================================================================================
            !     ====== find gainest, correct cold counts for xscan and spillover, store in xcount_tab  ==========
            !     =================================================================================================
            xcount_tab=1.e30    

            do iscan=1,numscan
                gainest(ich,iscan) = counts_per_deg(ich)    
                if(iflag_l0(iscan).ge.8) cycle ! no science or no geolocation data 
    
                ! find cold count average
                xsum1=0; xsum2=0
                do icel=1,ncel
                    xcount=csmcounts(ich,icel,iscan)
                    if(nint(xcount).eq.-32768) cycle
                    xcount=xcount - xscan_cal(1,ich,icel) !remove cell dependent error
                    xsum1=xsum1 + 1
                    xsum2=xsum2 + xcount
                enddo  !icel
                if(xsum1.eq.ncel) ccavg=xsum2/xsum1
    
                !     find hot count average
                xsum3=0; xsum4=0
                do icel=1,ncel
                    xcount=htscounts(ich,icel,iscan)
                    if(nint(xcount).eq.-32768) cycle
                    xcount=xcount - xscan_cal(2,ich,icel) !remove cell dependent error
                    xsum3=xsum3 + 1
                    xsum4=xsum4 + xcount
                enddo  !icel
                if(xsum3.eq.ncel) hcavg=xsum4/xsum3
    
                if(xsum1.eq.ncel .and. xsum3.eq.ncel) gainest(ich,iscan)=(hcavg-ccavg)/(ta_hot(ich,iscan)-ta_cold(ich,iscan))
    
                !     remove spillover so filter will work better    
                do icel=1,ncel
                    xcount=csmcounts(ich,icel,iscan)
                    if(nint(xcount).eq.-32768) cycle
                    xcount=xcount - xscan_cal(1,ich,icel) !remove cell dependent error
                    xcount_tab(icel,iscan)=xcount - spillover_coldmirror(ich)*gainest(ich,iscan)*ta_spill(ich,iscan)
                enddo  !icel  
            enddo  !iscan


            !     ==============================================================================================
            !     ================== find candidates for bad cold counts:  lflag_tab ===========================
            !     ==============================================================================================
            lflag_tab=0

            do iscan=1,numscan
                if(iflag_l0(iscan).ge.8)     cycle ! no science or no geolocation data 
                if(iflag_rx(ich,iscan).ne.0) cycle !scan next to agc change    

                num=0

                do kscan=iscan-iscan_del1,iscan+iscan_del1
                    if(kscan.lt.1 .or. kscan.gt.numscan)                             cycle
                    if(iflag_l0(kscan).ge.8)                                         cycle
                    if(iflag_rx(ich,kscan).ne.0)                                     cycle    
                    if(rx_offset_gain(ich, 1,kscan).ne.rx_offset_gain(ich, 1,iscan)) cycle
                    if(rx_offset_gain(ich, 2,kscan).ne.rx_offset_gain(ich, 2,iscan)) cycle     
    
                    do icel=1,ncel
                        xcount=xcount_tab(icel,kscan)
                        if(xcount.gt.1.e29) cycle
                        num=num+1
                        ival(   num)=nint(10*xcount)  !time 10 to get more precession out of median fileter
                        kscansv(num)=kscan
                        icelsv( num)=icel
                    enddo  !icel
                enddo  !kscan

                if(num.lt.100) cycle  !not enough obs to do meaning stats, should rarely if ever happen

                imid=(num+1)/4    !lower quarter
                ivalmid=fd_median(imid,num,ival)  !ivalmid is ten time counts

                do i=1,num
                    kscan=kscansv(i)
                    icel=  icelsv(i)
                    lflag_tab(0,icel,kscan)=lflag_tab(0,icel,kscan) + 1 !will not exceed 1+2*iscan_del1 so integer*2 can be used
                    if(0.1*(ival(i)-ivalmid).gt.2.0*stddv_avg(1,ich)*gainest(ich,kscan))  lflag_tab(1,icel,kscan)=lflag_tab(1,icel,kscan) + 1
                enddo  !num

            enddo  !iscan


            !     ==============================================================================================
            !     =============only flag as bad if suspected bad cc is a member of a group of bad c!     =======
            !     ==============================================================================================
            do iscan=1,numscan
                do icel=1,ncel
                    isum0=0; isum1=0
                    do kscan=iscan-2,iscan+2
                        if(kscan.lt.1 .or. kscan.gt.numscan) cycle
                        do kcel=icel-2,icel+2
                            if(kcel.lt.1 .or. kcel.gt.ncel)      cycle
                            isum0=isum0  + lflag_tab(0,kcel,kscan)
                            isum1=isum1  + lflag_tab(1,kcel,kscan)
                        enddo  !kcel
                    enddo  !kscan
                    if(isum1.lt.0.5*isum0) lflag_counts(ich,icel,iscan)=0
                enddo  !icel
            enddo  !iscan
    
    
        enddo  !ich



        !     =================================================================================================
        !     ================= find avg cold and hot counts via scan averaging  ==============================
        !     ==== find arrays: num_count(2,nch,maxscan), avg_count(2,nch,maxscan),rms_count(2,nch,maxscan) ===
        !     =================================================================================================
        num_count=0; avg_count=0; rms_count=0

        do ich=1,nch
      
            if(ich.le.12) then
                ncel=low_cal_obs
            else
                ncel=high_cal_obs
            endif
 
            do iload=1,2
                do iscan=1,numscan
                    if(iflag_l0(iscan).ge.8)     cycle ! no science or no geolocation data 
                    if(iflag_rx(ich,iscan).ne.0) cycle ! scan is at agc step    

                    num=0; isum0=0; xsum1=0; xsum2=0
                    
                    do iradius=0,iscan_del1
                        istep=2*iradius
                        if(istep.eq.0) istep=1

                        do kscan=iscan-iradius,iscan+iradius,istep
                        
                            if(kscan.lt.1 .or. kscan.gt.numscan)                           cycle
                            if(iflag_l0(kscan).ge.8)                                       cycle
                            if(iflag_rx(ich,kscan).ne.0)                                   cycle    
                            if(rx_offset_gain(ich,1,kscan).ne.rx_offset_gain(ich,1,iscan)) cycle
                            if(rx_offset_gain(ich,2,kscan).ne.rx_offset_gain(ich,2,iscan)) cycle     
                    
                            do icel=1,ncel
                                if(iload.eq.1) xcount=csmcounts(ich,icel,kscan)
                                if(iload.eq.2) xcount=htscounts(ich,icel,kscan)
                                if(nint(xcount).eq.-32768) cycle
                                xcount=xcount - xscan_cal(iload,ich,icel)
                                isum0=isum0 + 1
                                xsum1=xsum1 + xcount
                                xsum2=xsum2 + xcount**2
                                if(iload.eq.1 .and. lflag_counts(ich,icel,kscan).eq.1) cycle
                                if(iload.eq.1 .and. ich.ge.7 .and. iflag_moon(ich,icel,kscan).ge.2) cycle  !> 0.15K see 'O:\amsr2\moon_in_cold_mirror\memo3.txt'
                                if(iload.eq.2 .and. ich.eq.9 .and. icel.gt.10) cycle !see 'O:\amsr2\cold_counts\memo1.txt'
                                num=num+1
                                ival(num)=nint(10*xcount)  !time 10 to get more precession out of median fileter
                                icelsv(num)=icel
                            enddo  !icel
                
                        enddo  !kscan
                
                        if(iradius.ge.iscan_del2 .and. num.ge.100) exit !as soon as you get 100 counts for default radius, your done
            
                    enddo  !iradius

                    if(isum0.lt.5) cycle  !too few obs even to get rms_count, should rarely if ever happen

                    rms_count(iload,ich,iscan)=sqrt(abs(xsum2 - xsum1**2/isum0)/(isum0-1))

                    if(num.eq.0) cycle  !no good obs to get avg_count, should rarely if ever happen

                    imid=(num+1)/2
                    ivalmid=fd_median(imid,num,ival)  !ivalmid is ten time counts

                    isum0=0; xsum1=0; xsum2=0
                    do i=1,num
                        if(abs(0.1*(ival(i)-ivalmid)).gt.3.0*stddv_avg(iload,ich)*gainest(ich,iscan)) cycle
                        isum0=isum0 + 1
                        xsum1=xsum1 + 0.1*ival(i)
                        if(iload.eq.1 .and. ich.le.4) xsum2=xsum2 + spill_xscan(icelsv(i))
                    enddo !num
            
                    if(isum0.eq.0) stop 'isum0=0 in avg_cal_counts, pgm stopped'

                    num_count(iload,ich,iscan)=isum0 
                    avg_count(iload,ich,iscan)=xsum1/isum0

                    if(iload.eq.1) then
                    if(ich.le.4) spill_adj(ich,iscan)=xsum2/isum0 ! 'O:\amsr2\cold_counts\memo2.txt'
                    if(ich.ge.3 .and. ich.le.4) &
                        avg_count(iload,ich,iscan)=avg_count(iload,ich,iscan) - rfi_cold_count(1,iscan)*gainest(ich,iscan)/counts_per_deg(ich)    
                                                                        
                    if(ich.ge.5 .and. ich.le.6) &
                        avg_count(iload,ich,iscan)=avg_count(iload,ich,iscan) - rfi_cold_count(2,iscan)*gainest(ich,iscan)/counts_per_deg(ich)
                    endif                                                           

                enddo  !iscan
            enddo  !iload
        enddo  !ich

        return
    end subroutine avg_cal_counts

