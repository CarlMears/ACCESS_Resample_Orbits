!     9/18/2013 version change 9/23/2013 tcold_offset is now a function of freq

!     9/7/2013 version change 9/18/2013 spill_adj added

    subroutine fd_tcold(iopt_cld_spill)
        use, intrinsic :: IEEE_ARITHMETIC, only: ieee_is_finite
        use l2_module
        use date_utilities, only: fd_date_2000
        implicit none

        integer(4) iopt_cld_spill

        integer(4) iscan,ifreq,ich
        real(4) ta_spillx(nch)
        integer(4) iasc

        real(8) secyr,secdy    
        integer(4)     lyear,idayjl,imon,idaymo
    

        !     the first time this routine is call (ie iopt_cld_spill.eq.0) it sets spill_adj=1
        !     then avg_cal_counts will compute its actual value for average of the cold counts 
        !     the extreme values of spill_adj will be between 0.8 and 1.3, but more typically it will
        !     be near 1, depending on the distribution of icel. An even distribution of cells gives spill_adj=1.
        if(iopt_cld_spill.eq.0) spill_adj=1  

        do iscan=1,numscan

            !     ============================================================================================
            !     ================================== see if scan is ok =======================================
            !     ============================================================================================
            if(iflag_l0(iscan).ge.8) then ! no science, spacecraft or geoloc data
                ta_cold(:,iscan)= -999
                cycle
            endif

            if ((.not. ieee_is_finite(cellat_cold(iscan))) .or. &
                (.not. ieee_is_finite(cellon_cold(iscan))) .or. &
                (.not. ieee_is_finite(scan_time(iscan)))) then
                ta_cold(:,iscan) = -999
                cycle
            endif



            !     ============================================================================================
            !     ========================== correct 7 ghz tcold for cold mirror spillover ===================
            !     ============================================================================================
            call fd_date_2000(scan_time(iscan), secyr,lyear,idayjl,imon,idaymo,secdy)
    
            if(iflag_l0(iscan).le.4 .or. iflag_l0(iscan).eq.6) then !location for cold mirror measurement exists
                iasc=0 
                if(orbit(iscan)-int(orbit(iscan)).gt.0.5) iasc=1 
                call find_ta_spill(iopt_cld_spill,iasc,secyr,cellat_cold(iscan),cellon_cold(iscan), ta_spillx) 
                ta_spill(:,iscan)=ta_spillx
            else
                ta_spill(:,iscan)=150  !this is an extremely rare event and doesnt matter
            endif
    

            !     ============================================================================================
            !     ================================ loop through channels =====================================
            !     ============================================================================================
            do ich=1,nch
                ifreq=1 + int((ich-1)/2)
                ! amsre used a constant 3.4K rather than tcld_plank. amsre also used tcold_offset    
                ta_cold(ich,iscan)=tcld_plk(ifreq) + spill_adj(ich,iscan)*spillover_coldmirror(ich)*ta_spill(ich,iscan) + tcold_offset(ifreq)
            enddo  !ich
            
        enddo  !iscan

        return
    end subroutine fd_tcold
