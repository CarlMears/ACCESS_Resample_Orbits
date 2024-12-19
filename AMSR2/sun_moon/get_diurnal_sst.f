      subroutine get_diurnal_sst(sunlat,sundis,xhour,xlat,xlon,wind, diurnal_sst)
      implicit none

      real(4), parameter :: pi=3.141592654
	real(4), parameter :: solar_constant=1360   !watts/m**2

	real(4) sunlat,sundis,xhour,xlon,xlat,wind
      real(4) diurnal_sst, sinH, Q, H, rlmt,rr
      real(4) a0,a1,a2,a3,a4,a5,b1,b2,b3,b4,b5,w	 

      data a0,a1,a2,a3,a4,a5/4.118E-3, -4.132E-3, 0.8746E-3, -0.2460E-3, 0.2762E-3,-0.0609E-3/
      data    b1,b2,b3,b4,b5/          -5.093E-3, 2.5830E-3, -0.5143E-3,-0.3355E-3, 0.2269E-3/
      data w/0.261799388/  !2*pi/24.

c     calculate solar insolation
       
      rr=-tand(xlat)*tand(sunlat);
      if(rr> 1) rr= 1.
      if(rr<-1) rr=-1.

      H = acos(rr);	 
      sinH=sqrt(1 - rr**2) !sqrt(1-cos(h)**2) =sin(h)
      Q = (solar_constant/pi)*(sind(xlat)*sind(sunlat)*H + cosd(xlat)*cosd(sunlat)*sinH)/sundis**2

      rlmt=xhour + (xlon/360.)*24.

      if (Q>154) then
      diurnal_sst=(1.23459589668*(Q-154.)-1.56088414175378e-3*(Q-154.)**2)*exp(-.44*wind)*
     &				(a0 + a1*cos(   rlmt*w) + b1*sin(   rlmt*w) + a2*cos(2.*rlmt*w) + b2*sin(2.*rlmt*w)+
     &				      a3*cos(3.*rlmt*w) + b3*sin(3.*rlmt*w) + a4*cos(4.*rlmt*w) + b4*sin(4.*rlmt*w)+
     &                      a5*cos(5.*rlmt*w) + b5*sin(5.*rlmt*w))
	else
	diurnal_sst=0;
	endif

	RETURN
	END
