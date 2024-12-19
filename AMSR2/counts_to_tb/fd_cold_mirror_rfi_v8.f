        subroutine fd_cold_mirror_rfi(iorbit)
        use l2_module   
        use intrinsic
      implicit none
      
      integer(4) iorbit
      integer(4) j1,j2,k1,k2,ihemi,iscan,ifreq,iseg
      real(4) geolat,geolon
      real(4) b1,b2,c1,c2,brief
      real(4) rsq,dcnt
      
      rfi_cold_count=0  
      
	do iscan=1,numscan
	if(iflag_l0(iscan).ge.8) cycle ! no science or no geolocation data 
	
	if(zang(iscan).ge. 30 .and. zang(iscan).le.100) cycle  !geolat always .gt. 17 for these zangs
	if(zang(iscan).ge.180 .and. zang(iscan).le.310) cycle  !geolat always .gt. 17 for these zangs
	
	ihemi=1
	if(zang(iscan).gt.90..and. zang(iscan).le.270.) ihemi=2
	
	do ifreq=2,3
    call fd_geostationary_latlon(iscan,ifreq, geolat,geolon)
      if ((.not. ieee_is_finite(geolat)) .or. &
          (.not. ieee_is_finite(geolon))) then
            rfi_cold_count(1,iscan) = -999
            rfi_cold_count(2,iscan) = -999
            cycle
      endif
      if(geolon.lt.0. .or. geolon.gt.360.)     stop 'geolon ooob in fd_cold_mirror_rfi, pgm stopped'                       
	 
      if(abs(geolat).gt.15) cycle  !table is all zeros outside of +- 15 deg
       
      brief=2*(geolat+19.75)
      j1=1+brief                                                        
      j2=j1+1                                                           
      b1=j1-brief                                                       
      b2=1-b1

      brief=geolon-0.5 
      k1=1+brief                                                       
      k2=k1+1                                                           
      c1=k1-brief                                                       
      c2=1-c1                                                          
      if(k1.eq.  0) k1=360
      if(k2.eq.361) k2=  1 
      
                                                                       
      if(ifreq.eq.2) then
      iseg=ihemi
      rfi_cold_count(1,iscan)=b1*(c1*rfi_cold_mirror_map(k1,j1,iseg) + c2*rfi_cold_mirror_map(k2,j1,iseg)) + &
                             b2*(c1*rfi_cold_mirror_map(k1,j2,iseg) + c2*rfi_cold_mirror_map(k2,j2,iseg)) 
      endif
      
      
      if(ifreq.eq.3 ) then
      dcnt=0
      
      if(ihemi.eq.1) then
      iseg=3              
      if(iorbit.ge.23001) then
      rsq=(2*(geolat+2.75))**2 + (geolon-56.5)**2
      dcnt=1.0*(1 -rsq/50.)
      if(dcnt.lt.0) dcnt=0
      endif
      endif
      
      if(ihemi.eq.2) then
      iseg=4
      if(iorbit.ge.10001) iseg=5  
      if(iorbit.ge.20001) iseg=6  
      if(iorbit.ge.19001 .and. iorbit.le.20000) then
      rsq=(2*(geolat+0.25))**2 + (geolon-56.5)**2
      dcnt=2.0*(1 -rsq/50.)
      if(dcnt.lt.0) dcnt=0
      endif
      endif
      
      rfi_cold_count(2,iscan)=b1*(c1*rfi_cold_mirror_map(k1,j1,iseg) + c2*rfi_cold_mirror_map(k2,j1,iseg)) + &
                              b2*(c1*rfi_cold_mirror_map(k1,j2,iseg) + c2*rfi_cold_mirror_map(k2,j2,iseg))  + dcnt
     
      endif
      
      enddo  !ifreq
      enddo  !iscan
      
      return                                                            
      end   
      
      
