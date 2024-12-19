!     feb 7 2013 version changed 8/29/2013.  with the new jscan convention, jscan=2*numscan-1 has no data.
!     thus this scan is skipped.  it would have been skipped anyway because the tas are -999.

!     i verified that there is never more that 4 obs in a lat/lon cell except sometimes during rpy manuevers.
!     here i simply limit the number of obs/cell to be 4.

      subroutine resample_89(iorbit)
	use l2_module						 					
	implicit none

	integer(4) iscan,jscan,icel,jcel,ilat,ilon,jlat,jlon,klon,nlon,idellat,idellon,ihorn,iorbit,ihalf,iscansv,iscan1,ich1,ich2
	integer(4) jsumv,jsumh,num
	real(4) xlat0,xlon0,coslat,coslatsq,xterm,frcrev,tax,dissq
	real(8) xsumv,xsumh

	integer(1) isumv(8640,4320),isumh(8640,4320)
	integer(2) sumv(8640,4320),sumh(8640,4320)

	real(4) ta89_resp(maxcel,maxscan,13:16)
	common/resampled89/ ta89_resp


	ta89_resp=0	 ! set ta89_resp to default values

	do 300 ihorn=1,2

	if(ihorn.eq.1) then
	ich1=13; ich2=14
	else
	ich1=15; ich2=16
	endif

	do 250 ihalf=1,2

	if(ihalf.eq.1) then
	iscan1=1
	iscansv=999999
	else
	if(iscansv.eq.999999) goto 250	 !short orbit containing no second half
	iscan1=iscansv
	endif


	isumv=0; isumh=0; sumv=0; sumh=0

	do iscan=1,numscan
	if(iflag_l0(iscan).ne.0) cycle

	frcrev=orbit(iscan)-iorbit
	if(ihalf.eq.1 .and. frcrev.gt.0.75) cycle
	if(ihalf.eq.2 .and. frcrev.lt.0.25) cycle

!     amsr2 uses even scans for 89a, odd scans for 89b.  this is the reverse of amsre
	if(ihorn.eq.2) then 
	jscan=2*iscan-1
	else
	jscan=2*iscan
	endif
	
	if(jscan.eq.2*numscan-1) cycle  !no 89b data for this scan

	do jcel=1,maxcel_89

      ilat=1+nint(24.*(cellat_89(jcel,jscan)+89.97917))
      ilon=1+nint(24.*(cellon_89(jcel,jscan)- 0.02083))
      if(ilat.lt.1) ilat=1
      if(ilon.lt.1) ilon=1
      if(ilat.gt.4320) ilat=4320
      if(ilon.gt.8640) ilon=8640
	
	tax=ta_89(jcel,jscan,1)
	if(tax.ge.55 .and. tax.le.330 .and. iflag_cal(ich1,iscan).eq.0 .and. isumv(ilon,ilat).lt.4) then 
	isumv(ilon,ilat)=isumv(ilon,ilat) + 1
	 sumv(ilon,ilat)= sumv(ilon,ilat) + nint(50*(tax-200))
	endif

	tax=ta_89(jcel,jscan,2)
	if(tax.ge.55 .and. tax.le.330 .and. iflag_cal(ich2,iscan).eq.0 .and. isumh(ilon,ilat).lt.4) then 
	isumh(ilon,ilat)=isumh(ilon,ilat) + 1
	if(isumh(ilon,ilat).eq.5) stop 'sumh oob in resample_89, pgm stopped'
	 sumh(ilon,ilat)= sumh(ilon,ilat) + nint(50*(tax-200))
	endif

	enddo  !jcel
	enddo  !iscan


	do 200 iscan=iscan1,numscan
	if(iflag_l0(iscan).ne.0) cycle

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
	nlon=nint(5./coslat)
	else
	nlon=nint(5./0.02)
	endif

      ilat=1+nint(24.*(xlat0+89.97917))
      ilon=1+nint(24.*(xlon0- 0.02083))
      if(ilat.lt.1) ilat=1
      if(ilon.lt.1) ilon=1
      if(ilat.gt.4320) ilat=4320
      if(ilon.gt.8640) ilon=8640

	jsumv=0; jsumh=0; xsumv=0; xsumh=0

	do klon=ilon-nlon,ilon+nlon
	idellon=klon-ilon
	if(abs(idellon*coslat).gt.4) cycle
	jlon=klon
	if(jlon.lt.   1) jlon=jlon+8640
	if(jlon.gt.8640) jlon=jlon-8640	
	xterm=coslatsq*idellon*idellon 

	do jlat=ilat-4,ilat+4
	if(jlat.lt.1 .or. jlat.gt.4320) cycle
	idellat=jlat-ilat
	dissq=idellat*idellat + xterm
	if(dissq.gt.10.518) cycle	  ! dis=111*sqrt(dissq/576.), where max dis is 15 km

	num=isumv(jlon,jlat)
	if(num.ne.0) then
      jsumv=jsumv + num
	xsumv=xsumv + dble(sumv(jlon,jlat))
	endif

	num=isumh(jlon,jlat)
	if(num.ne.0) then
      jsumh=jsumh + num
	xsumh=xsumh + dble(sumh(jlon,jlat))
	endif

	enddo !klon
	enddo !jlat


	if(jsumv.ne.0) ta89_resp(icel,iscan,ich1)=0.02*xsumv/jsumv	+ 200.
	if(jsumh.ne.0) ta89_resp(icel,iscan,ich2)=0.02*xsumh/jsumh	+ 200.
			
  100 continue
  200 continue
  250 continue  
  300 continue

      return
      end
