!     2/5/2013 version changed on 8/22/2013.  subroutine fd_tbsss_clm removed.
!     this subroutine is being modifed and I wanted it in a separate file.
!     all with this 8/22/2013 update, the name of this module was changed from ocean_routines to field_routines

!     for amsr2 i took out subroutine salinity_adjustment so that it would not be used improperly
!     the freq indexing in subroutine salinity_adjustment was for amsre and i no longer need this routine



!     same as 8/27/2005 version except:
!     1.  module name changed to l2_module
!     2.  iscan_flag used in place of iflag_l0
!     3.  also ckocntb and chkpol have been removed
!     4.  routine salinity_adjustment changed a little to avoid apply sal correction to negative tb

!    12/12/03 version changed on 5/19/2004.  routine tbspace_icealgo was removed and coslat not allowed to go to zero

	subroutine fieldavg(numr,icel,iscan, wfield,vfield)
	use l2_module
	implicit none

	integer(4) numr,icel,jcel,iscan,jscan,iflagw,iflagv,iradius
	real(4) wfield,vfield,xlat0,xlon0,coslat
	real(4)	wsum1,wsum2,vsum1,vsum2,wt

      if(irainadj(1,icel,iscan).eq.0) then
	wfield=ep(4,icel,iscan)
	vfield=ep(5,icel,iscan)
      goto 200
	endif

	xlat0=cellat(icel,iscan)
	xlon0=cellon(icel,iscan)
	coslat=cosd(xlat0)

      wsum1=0; wsum2=0; vsum1=0; vsum2=0;	iflagw=0; iflagv=0

	do iradius=0,numr

	if(iradius.eq.0) then
	call sumep(icel,iscan,xlat0,xlon0,coslat,iflagw,iflagv, wsum1,wsum2,vsum1,vsum2)

	else

	do jcel=icel-iradius,icel+iradius
	if(jcel.lt.1 .or. jcel.gt.maxcel)    cycle

	do  jscan=iscan-iradius,iscan+iradius,2*iradius
	if(jscan.lt.1 .or. jscan.gt.numscan) cycle
	if(iscan_flag(jscan).ne.0)             cycle
      if(isurcel(3,jcel,jscan).ne.0)         cycle !only do open ocean
      if(ice_flag2(jcel,jscan).ne.0)       cycle !skip ice

	call sumep(jcel,jscan,xlat0,xlon0,coslat,iflagw,iflagv, wsum1,wsum2,vsum1,vsum2)
	enddo
	enddo

	do jscan=iscan-iradius+1,iscan+iradius-1
	if(jscan.lt.1 .or. jscan.gt.numscan) cycle
	if(iscan_flag(jscan).ne.0)             cycle

	do jcel=icel-iradius,icel+iradius,2*iradius
	if(jcel.lt.1 .or. jcel.gt.maxcel)    cycle
      if(isurcel(3,jcel,jscan).ne.0)         cycle !only do open ocean
      if(ice_flag2(jcel,jscan).ne.0)       cycle !skip ice

	call sumep(jcel,jscan,xlat0,xlon0,coslat,iflagw,iflagv, wsum1,wsum2,vsum1,vsum2)
	enddo
	enddo
	endif

	if(wsum1.gt.150) iflagw=1
	if(vsum1.gt.150) iflagv=1
	if(iflagw.eq.1 .and. iflagv.eq.1) exit
	enddo

!     add in climate if there is little adjacent data
	if(iflagw.eq.0 .or. iflagv.eq.0) then
	if(iflagw.eq.0) then
	wt=10
      wsum1=wsum1+wt
      wsum2=wsum2+wt*wincl(icel,iscan)
	endif
	if(iflagv.eq.0) then
	wt=10
      vsum1=vsum1+wt
      vsum2=vsum2+wt*vapcl(icel,iscan)
	endif
	endif

      wfield=wsum2/wsum1
      vfield=vsum2/vsum1

  200 continue
      if(wfield.lt.  0) wfield= 0
      if(wfield.gt. 25) wfield=25
      if(vfield.lt.  0) vfield= 0
      if(vfield.gt. 68) vfield=68

	return
	end



	subroutine sumep(jcel,jscan,xlat0,xlon0,coslat,iflagw,iflagv, wsum1,wsum2,vsum1,vsum2)
	use l2_module
	implicit none

	integer(4) jcel,jscan,iflagw,iflagv
	real(4) xlat0,xlon0,coslat,wsum1,wsum2,vsum1,vsum2,diflat,diflon,dissq,wt,wind,vapor,totliq

      if(iret_flag(3,jcel,jscan).ne.0) return 

	diflat=cellat(jcel,jscan)-xlat0
	diflon=cellon(jcel,jscan)-xlon0
	if(diflon.lt.-180) diflon=diflon+360
	if(diflon.gt. 180) diflon=diflon-360
	diflon=coslat*diflon
	dissq=diflat*diflat + diflon*diflon
	if(dissq.gt.3.2465)   return !.gt. 200 km

	if(dissq.gt.0.0001) then
	wt=1./dissq
	else 
	wt=10000
	endif

	wind=  ep(4,jcel,jscan)
	vapor= ep(5,jcel,jscan)
      totliq=ep(6,jcel,jscan)

      if(irainadj(1,jcel,jscan).eq.0 .and. iflagw.eq.0) then
      wsum1=wsum1+wt
      wsum2=wsum2+wt*wind
	endif

      if(totliq.le.1.00 .and. iflagv.eq.0) then
      vsum1=vsum1+wt
      vsum2=vsum2+wt*vapor
	endif

	return
	end


!     i include ice cells (i.e., ice_flag2=1) because ice sometimes causes rain
!     and hence adjrain can also serve as a near to ice flag

	subroutine adjrain(icase,iorbit)
	use l2_module
	implicit none

	integer(4) iscan,icel,ilat,ilon,jlat,jlon,klon,nlon,idellat,idellon,iorbit,ihalf,iscansv,iscan1,icase,irain
	real(4) xlat0,xlon0,coslat,coslatsq,xterm,frcrev,dissq,vapor,cloud,totliq

	integer(1) irain_map(8640,4320)

      irainadj(:,:,1:numscan)=0 

	do 250 ihalf=1,2

	if(ihalf.eq.1) then
	iscan1=1
	iscansv=999999
	else
	if(iscansv.eq.999999) goto 250	 !short orbit containing no second half
	iscan1=iscansv
	endif

	irain_map=-2	 ! set to default value indicating no obs

	do iscan=1,numscan
	if(iscan_flag(iscan).ne.0) cycle

	frcrev=orbit(iscan)-iorbit
	if(ihalf.eq.1 .and. frcrev.gt.0.75) cycle
	if(ihalf.eq.2 .and. frcrev.lt.0.25) cycle

	do icel=1,maxcel

      ilat=1+nint(24.*(cellat(icel,iscan)+89.97917))
      ilon=1+nint(24.*(cellon(icel,iscan)- 0.02083))
      if(ilat.lt.1) ilat=1
      if(ilon.lt.1) ilon=1
      if(ilat.gt.4320) ilat=4320
      if(ilon.gt.8640) ilon=8640

      if(icase.eq.1) then
      if(iret_flag(3,icel,iscan).gt.0) cycle
      totliq= ep(6,icel,iscan)
	if(totliq.le.0.18) then
      cloud=totliq
	else
	cloud=totliq -0.45*(totliq-0.18) 
	endif

	else
      if(iret_flag(4,icel,iscan).gt.1) cycle
      cloud= ep(7,icel,iscan) 
	endif

	if(cloud.le.0.18) then
	irain=-1  !this means no rain
	else
      vapor = ep(5,icel,iscan)
	call fdirain(vapor,cloud, irain)  !0 mean very light rain, 1 means heavier rian
	endif

	if(irain.gt.irain_map(ilon,ilat)) irain_map(ilon,ilat)=irain
	
	enddo  !icel
	enddo  !scan




	do 200 iscan=iscan1,numscan
	if(iscan_flag(iscan).ne.0) cycle

	frcrev=orbit(iscan)-iorbit
	if(ihalf.eq.1 .and. frcrev.gt.0.50) then
	iscansv=iscan
	exit
	endif


      do 100 icel=1,maxcel
	xlat0=cellat(icel,iscan)
	xlon0=cellon(icel,iscan)
	coslat=cosd(xlat0)
	coslatsq=coslat*coslat

	if(coslat.gt.0.02) then
	nlon=nint(10./coslat)
	else
	nlon=nint(10./0.02)
	endif

      ilat=1+nint(24.*(xlat0+89.97917))
      ilon=1+nint(24.*(xlon0- 0.02083))
      if(ilat.lt.1) ilat=1
      if(ilon.lt.1) ilon=1
      if(ilat.gt.4320) ilat=4320
      if(ilon.gt.8640) ilon=8640


	do klon=ilon-nlon,ilon+nlon
	idellon=klon-ilon
	if(abs(idellon*coslat).gt.9) cycle
	jlon=klon
	if(jlon.lt.   1) jlon=jlon+8640
	if(jlon.gt.8640) jlon=jlon-8640	
	xterm=coslatsq*idellon*idellon 

	do jlat=ilat-9,ilat+9
	if(jlat.lt.1 .or. jlat.gt.4320) cycle
	idellat=jlat-ilat
	dissq=idellat*idellat + xterm

	if(dissq.gt.74.799) cycle	  ! dis=111*sqrt(dissq/576.), where max dis is 40 km
 
	if(irain_map(jlon,jlat).lt.0)   cycle

	if(dissq.le.29.218) then !.le. 25 km
	if(irain_map(jlon,jlat).eq.1) irainadj(1,icel,iscan)=1  !.le. 25 km, light rain
	irainadj(2,icel,iscan)=1      !.le. 25 km, any rain
	endif

	if(irain_map(jlon,jlat).eq.1) irainadj(3,icel,iscan)=1  !.le. 40 km, light rain
	irainadj(4,icel,iscan)=1      !.le. 40 km, any rain

	enddo !jlat
	enddo !klon

			
  100 continue
  200 continue
  250 continue  

      return
      end



	subroutine fdirain(vapor,cloud, irain)
	implicit none

	integer(4) irain
	real(4) vapor,cloud,u,wt,cloud_limit

	if(cloud.gt.0.31) then !0.31 corresponds to rr=0.5 km mm/hr
      irain=1
      return
	endif

      u=(vapor-30)/20.
      if(u.ge.0 .and. u.le.1) then
      wt=1-u*u*(3-2*u)
      else
      wt=1
      if(u.gt.0.5) wt=0
      endif

      cloud_limit=0.18+0.13*wt
	if(cloud.gt.cloud_limit) then
	irain=1
	else
	irain=0
	endif

	return
	end


