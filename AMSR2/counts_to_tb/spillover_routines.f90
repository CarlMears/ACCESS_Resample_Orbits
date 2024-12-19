!     feb 5 2013 version change aug 28 2013.  This program is now setup to do the first 6 channels
!     also each channels is done separately.  previously i required both 6.9 vpol and hpol to be available
!     but now each channels is checked separately.  There was also signficant restructuring that put arrays into the l2_module
!     and reading tables into read_in

!     december 12, 2003, 03:38:21 pm version changed aug 28 2004:     coslat is not allowed to exceed 0.02 when computing nlon. will have neglible effect on previous runs.

!     these routine are the same as it mk_cold_cnt_map.f version 10/23/2003
!     except that i fixed the bug mentioned in memo1 and the upper bound on ta
!     has been increased to 400 to handle rfi.

    subroutine mk_spill_tamap
        use l2_module
        implicit none
        
        integer(4), parameter :: nlat = 180, nlon = 360
        integer(4)  iscan,icel,ilat,ilon,ich

        ta_map=0

        do iscan=1,numscan
            if(iflag_l0(iscan).eq.5 .or. iflag_l0(iscan).ge.7) cycle !earth location not done for this scan or no science data

            do icel=1,maxcel
                ilon=1+int(cellon(icel,iscan))
                ilat=1+int(cellat(icel,iscan)+90.)
                if(ilon.lt.1 .or. ilon.gt.nlon .or. ilat.lt.1 .or. ilat.gt.nlat) stop 'error in mk_spill_tamap, pgm stopped'
                
                do ich=1,6
                    if(iflag_cal(ich,iscan).gt.0) cycle
                    if(ta(icel,iscan,ich).le.55 .or. ta(icel,iscan,ich).ge.400) cycle  !large value of 400 use for rfi
                    ta_map(ilon,ilat,ich,0)=ta_map(ilon,ilat,ich,0) + 1
                    ta_map(ilon,ilat,ich,1)=ta_map(ilon,ilat,ich,1) + ta(icel,iscan,ich)
                enddo  !ich
            enddo !icel
        enddo !iscan

        return
    end subroutine mk_spill_tamap

    subroutine find_ta_spill(iopt,iasc,secyr,xlat,xlon, tax) !iopt=0 means only do climate map
        use l2_module
        implicit none

        integer(4), parameter :: nlat = 180, nlon = 360
        integer(4) iopt,iasc,nstep,ilat,ilon,klon,jlon,jlat,idellat,idellon,nnnlon,isum(nch),ich
        real(4) xlat,xlon,tax(nch)
        real(4) coslat,coslatsq,xterm,dissq
        real(8) secyr
        
        data nstep/3/

        call fd_tamap_climate(iasc,secyr,xlat,xlon, tax)
        if(iopt.eq.0) return
      
        !     for full orbits isum ragnes from 3000 to 7000.  for climate data i set isum= 300
        !     which gives it no more than 10% weight
        isum=300
        tax=tax*isum

        ilon=1+int(xlon)
        ilat=1+int(xlat+90.)
        if(ilon.lt.1.or.ilon.gt.nlon.or.ilat.lt.1.or.ilat.gt.nlat) stop '******* error1 in ilon,ilat'

        coslat=cosd(xlat)
        coslatsq=coslat*coslat

        if(coslat.gt.0.02) then
            nnnlon=nint((nstep+1)/coslat)
        else
            nnnlon=nint((nstep+1)/0.02)
        endif

        if(nnnlon.gt.0.5*nlon) nnnlon=0.5*nlon

        do klon=ilon-nnnlon,ilon+nnnlon
            idellon=klon-ilon
            if(abs(idellon*coslat).gt.nstep) cycle
            jlon=klon
            if(jlon.lt.   1) jlon=jlon+nlon 
            if(jlon.gt.nlon) jlon=jlon-nlon
            xterm=coslatsq*idellon*idellon 
            do jlat=ilat-nstep,ilat+nstep
                if(jlat.lt.1 .or. jlat.gt.nlat) cycle
                idellat=jlat-ilat
                dissq=idellat*idellat + xterm !
                if(dissq.gt.7.305) cycle       ! dis=111*sqrt(dissq), where max dis is 300 km
                do ich=1,6
                    isum(ich)=isum(ich) + nint(ta_map(jlon,jlat,ich,0))
                    tax(ich) =tax(ich)   +      ta_map(jlon,jlat,ich,1)
                enddo
            enddo !jlat

        enddo !klon

        tax(1:6)=tax(1:6)/isum(1:6)
        tax(7:nch)=0  !value does not matter since spillover for these channels is 0                 

        return
    end subroutine find_ta_spill

    subroutine fd_tamap_climate(iasc,secyr,xlat,xlon, tax)                            
        use l2_module
        implicit none

        !     character(*) filename_cband_map
        real(8) secyr                                                        
        real(4) xlat,xlon,tax(nch)                                                                 
        real(4) a1,a2,b1,b2,c1,c2,brief
        integer(4) iasc,i1,i2,j1,j2,k1,k2,ich

        !     check inputs                                                      
                                                                        
        if(secyr.lt.0 .or. secyr.gt.31622400) stop 'error1 in fdtamap'                     
        if(xlon.lt.0 .or. xlon.gt.360)        stop 'error2 in fdtamap'       
        if(abs(xlat).gt.90)                   stop 'error3 in fdtamap'   
        
        

        !     do time,lat,lon interpolation                                    
        !     2629800 is the average num sec in a month, 86400*365.25/12        

        brief=(secyr-1314900)/2629800.d0                                  
        i1=1+brief                                                        
        i2=i1+1                                                           
        a1=i1-brief                                                       
        a2=1-a1                                                          
        if(i1.eq. 0) i1=12
        if(i2.eq.13) i2= 1

        brief=xlat+89.5
        j1=1+brief                                                        
        j2=j1+1                                                           
        b1=j1-brief                                                       
        b2=1-b1
        if(j1.eq.  0) j1=  1
        if(j2.eq.181) j2=180

        brief=xlon-0.5
        k1=1+brief                                                       
        k2=k1+1                                                           
        c1=k1-brief                                                       
        c2=1-c1                                                          
        if(k1.eq.  0) k1=360
        if(k2.eq.361) k2=  1 
   
        do ich=1,6 !6.9,7.3,10.7                                                                  
            tax(ich)=    &          
                a1*b1*(c1*tamap(k1,j1,i1,ich,iasc)+c2*tamap(k2,j1,i1,ich,iasc))+  &               
                a1*b2*(c1*tamap(k1,j2,i1,ich,iasc)+c2*tamap(k2,j2,i1,ich,iasc))+  &               
                a2*b1*(c1*tamap(k1,j1,i2,ich,iasc)+c2*tamap(k2,j1,i2,ich,iasc))+  &                 
                a2*b2*(c1*tamap(k1,j2,i2,ich,iasc)+c2*tamap(k2,j2,i2,ich,iasc))
        enddo
      
        tax(7:nch)=0  !value does not matter since spillover for these channels is 0                 

        return                                                            
      end subroutine fd_tamap_climate
