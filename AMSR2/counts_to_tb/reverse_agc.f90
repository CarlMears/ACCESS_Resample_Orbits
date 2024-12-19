!     9/12/2013 version change 9/17/2013.  if agc adjusted count falls outside -32767 to 32737, count is flagged as bad
!     but i dont think this every happens

    subroutine reverse_agc
        use l2_module
        implicit none
    
        integer(4), parameter :: ioffset0(nch)= &
            (/129, 148,  91,  99,  72,  93,  68,  65,  53,  68,  78,  79, 158, 146, 143, 145/)
        integer(4), parameter ::   igain0(nch)= &
            (/ 56,  64,  79,  75,  65,  40,  50,  83, 101,  77,  62,  72,  58,  61,  66,  66/)

        !     see 'O:\amsr2\agc\memo1.txt' for following values.      
        real(8), parameter :: d_offset(nch)=(/20.5,  21.7,  16.1,  17.7,  31.8,  24.7,  66.5,  62.8,  89.1,  62.2,  79.1,  76.7, &
                                            100.2, 107.7, 102.6, 113.3/)

        real(8), parameter :: d_gain(nch)=(/1.0068,1.0062,1.0064,1.0066,1.0107,1.0143,1.0080,1.0070,1.0056,1.0063,1.0080,1.0069, &
                                            1.0092,1.0103,1.0097,1.0092/)

        integer(4) iscan,ich,icel,ncel
        real(4) doffset,dgain,xcount
    
        do ich=1,nch

            do iscan=1,numscan
                if(rx_offset_gain(ich,1,iscan).le.0 .or. rx_offset_gain(ich,1,iscan).ge.255) cycle
                if(rx_offset_gain(ich,2,iscan).le.0 .or. rx_offset_gain(ich,2,iscan).ge.255) cycle
                doffset=d_offset(ich)*( rx_offset_gain(ich,1,iscan)-ioffset0(ich))
                dgain  =  d_gain(ich)**(rx_offset_gain(ich,2,iscan)-  igain0(ich))
            
            
                !     =====================================================================================================================
                !     ==========================================   cal counts ============================================================
                !     =====================================================================================================================
                if(ich.le.12) then
                ncel=low_cal_obs
                else
                ncel=high_cal_obs
                endif

                do icel=1,ncel

                    if(csmcounts(ich,icel,iscan).lt.iminc .or. csmcounts(ich,icel,iscan).gt.imaxc) then
                        csmcounts(ich,icel,iscan)=-32768
                    else
                        xcount=csmcounts(ich,icel,iscan)/dgain + doffset  
                        if(abs(xcount).gt.32767) xcount=-32768
                        csmcounts(ich,icel,iscan)=nint(xcount)  
                    endif
    
                    if(htscounts(ich,icel,iscan).lt.iminc .or. htscounts(ich,icel,iscan).gt.imaxc) then
                        htscounts(ich,icel,iscan)=-32768
                    else
                        xcount=htscounts(ich,icel,iscan)/dgain + doffset 
                        if(abs(xcount).gt.32767) xcount=-32768
                        htscounts(ich,icel,iscan)=nint(xcount) 
                    endif
    
                enddo  !icel
    
                !     =====================================================================================================================
                !     ==========================================  earth counts ============================================================
                !     =====================================================================================================================
                if(ich.le.12) then
                ncel=maxcel
                else
                ncel=maxcel_89
                endif

                do icel=1,ncel
                    if(earthcounts(ich,icel,iscan).lt.iminc .or. earthcounts(ich,icel,iscan).gt.imaxc) then
                        earthcounts(ich,icel,iscan)=-32768
                    else
                        xcount=earthcounts(ich,icel,iscan)/dgain + doffset
                        if(abs(xcount).gt.32767) xcount=-32768
                        earthcounts(ich,icel,iscan)=nint(xcount)
                    endif
                enddo  !icel
            enddo  !iscan
        enddo  !ich    
        return
        end subroutine reverse_agc
      
