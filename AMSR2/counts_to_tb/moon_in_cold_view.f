!     this routine corrects cold counts for 6.9, 7.3, and 10.7 ghz channels
!     it also makes a table of moon contamination in units of 0.1 K.
    subroutine moon_in_cold_view
    use l2_module
    implicit none
    real(4),    parameter :: moondis0=370000.     !km
    
    integer(4) iscan,ifreq,ich,ichv,ichh,icel,ncel

    real(4) moon_angle,arg,cc,ch,dif_count,ta_moon,cc_moon
    real(8) b(3),cosbeta,sinbeta,bx,by,bz
    real(8) sun_to_moon(3),sc_to_moon(3),sun_angle,moondis_km_sq
     
    do iscan=1,numscan
    if(iflag_l0(iscan).ge.8) cycle ! no science or no geolocation data 

    sun_to_moon=moonvec0(:,iscan) - sunvec0(:,iscan) 
    sun_to_moon=sun_to_moon/sqrt(dot_product(sun_to_moon,sun_to_moon))
    
    sc_to_moon=moonvec0(:,iscan) - 1.d-3*scpos_eci0(:,iscan)  !scpos is in meters, it is convert to km here
    sc_to_moon=sc_to_moon/sqrt(dot_product(sc_to_moon,sc_to_moon))

    moondis_km_sq=dot_product(moonvec0(:,iscan),moonvec0(:,iscan))

    arg=sqrt(dot_product(sun_to_moon-sc_to_moon,sun_to_moon-sc_to_moon))/2.
    if(arg.gt.1) arg=1
    sun_angle =2*asind(arg)

    do ifreq=1,nfreq

    if(ifreq.le.6) then
    ncel=low_cal_obs
    else
    ncel=high_cal_obs
    endif

    do icel=1,ncel
    
      cosbeta=cosd(beta_cold(icel,ifreq))                                                                      
      sinbeta=sind(beta_cold(icel,ifreq))                                                                      
      bx= sinbeta*cosd(alpha_cold(icel,ifreq)) !input alpha defined relative to the -z axes                                       
      by=-sinbeta*sind(alpha_cold(icel,ifreq)) !input alpha defined relative to the -z axes                                       
      bz=cosbeta                 
      b(1)=bx*x_eci0(1,iscan)+by*y_eci0(1,iscan)+bz*z_eci0(1,iscan)                            
      b(2)=bx*x_eci0(2,iscan)+by*y_eci0(2,iscan)+bz*z_eci0(2,iscan)                            
      b(3)=bx*x_eci0(3,iscan)+by*y_eci0(3,iscan)+bz*z_eci0(3,iscan)

    arg=sqrt(dot_product(b-sc_to_moon,b-sc_to_moon))/2.
    if(arg.gt.1) arg=1
    moon_angle=2*asind(arg)

    arg=-(moon_angle/bw_cold(icel,ifreq))**2
    if(arg.lt.-30) arg=-30

    ta_moon= (1-phase_factor_moon(ifreq)*(sun_angle-80.))*(amp_moon(ifreq)/bw_cold(icel,ifreq)**2)*exp(arg)*moondis0**2/moondis_km_sq
     
    ichv=2*ifreq-1
    ichh=2*ifreq

    
      if(ta_moon.le.12) then
      iflag_moon(ichv:ichh,icel,iscan)=nint(10*ta_moon)
      else
      iflag_moon(ichv:ichh,icel,iscan)=120
      endif

    if(ifreq.ge.4 .or. ta_moon.lt.0.05) cycle  !do not do the correction for 19 ghz and above see 'O:\amsr2\moon_in_cold_mirror\memo3.txt'

    do ich=ichv,ichh

    cc=csmcounts(ich,icel,iscan)
    ch=htscounts(ich,icel,iscan)
    
    if(cc.eq.-32768)   cycle
    if(ch.eq.-32768)   cycle
    dif_count=ch-cc
    if(dif_count.lt.100) cycle

    cc_moon=dif_count*ta_moon/(ta_hot(ich,iscan)-(ta_cold(ich,iscan) + ta_moon))
    
    csmcounts(ich,icel,iscan)=nint(csmcounts(ich,icel,iscan)-cc_moon)
    
    enddo  !ich
    enddo  !icel
    enddo  !ifreq
    enddo  !iscan
    
      return
    end
