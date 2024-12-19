      subroutine get_diurnal_sst_v2(sunlat,sunlon,sundis,hour,xlat,xlon,wind, diurnal_sst)
      implicit none

	real(4), parameter :: solar_flux_surface=750 !watts/m**2, typical surface value
	real(4), parameter :: cp=4185.5              !heat capacity joules/(kg C)
	real(4), parameter :: rho=998.               !water density kg/m**2
	real(4), parameter :: d0=0.37                 !subskin depth meters
	real(4), parameter :: coef=1./(rho*cp)
      real(4), parameter :: dt=360.  !time step in seconds
      
	real(4) sunlat,sunlon,sundis,hour,xlat,xlon,wind
      real(4) diurnal_sst
      
      integer(4) jhour
      real(4) sunlonx,flux,term1,term2,zhour,hour_lag,dlag,hour_local,dang,sunlon_sunrise
      real(4) dhour,dhr,d
      real(8) tsum
      
      hour_local=hour + xlon/15.
      hour_lag=2.75
      dlag=cosd(15.*( hour_local-1.))
      if(dlag.gt.0) hour_lag=2.75 + 1.1*dlag**2

      term1=sind(xlat)*sind(sunlat)
      term2=cosd(xlat)*cosd(sunlat)
      
      
c     =============================================================================================
c     ======================================= find hour of sunrise ================================
c     =============================================================================================
      if(abs(term1/term2).gt.1) then
      dhr=-999  !there is no sunrise, either all day or all  night

      else
      dang=acosd(-term1/term2)  !0 to 180
      sunlon_sunrise=xlon + dang !sun longitude decreases with time, thus to go back in time and find longitude at sunrise, you need to add dang
      dhr=(sunlon_sunrise-sunlon)/15.  !hour after sunrise
      if(dhr.ge.24) dhr=dhr-24.
      if(dhr.ge.24) dhr=dhr-24.
      if(dhr.lt. 0) dhr=dhr+24.
	endif
     
c     =============================================================================================
c     ==================================== backward integration ===================================
c     =============================================================================================
      tsum=0
      do jhour=0,120
      zhour=0.1*jhour
      sunlonx=sunlon + 15.*zhour
      flux=term1 + term2*cosd(sunlonx-xlon) 
      
      if(flux.le.0) then
      flux=0
      
      else
      
      if(flux.lt.0.01) then
      flux=0
      else
      flux=flux*exp(-0.2/flux)
      endif
      
      if(dhr.gt.-998) then
      dhour=dhr-zhour
 	if(dhour.lt.0) dhour=dhour+24.
      if(dhour.le.4) flux=(1-exp(-0.25*dhour**2))*flux  !speed up program, dhour>4 represents essentially no change
      endif
      
      endif
     
      tsum=tsum + exp(-zhour/hour_lag)*flux
      enddo  !jhour
      
      d=d0*exp(0.46*wind)
      
	diurnal_sst=(coef/d)*solar_flux_surface*tsum*dt/sundis**2 
	
	return
	end
	