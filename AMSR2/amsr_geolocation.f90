!    7/13/2021 Convert to module
!
!    8/29/2013 version updated 9/20/2103
!    lots of code was being used to compute celrfi which is only used for writing to the l2b data file.
!    and ta_inter_calibration does not use it.  And surprisingly as far as I can tell no program that
!    reads the l2b files (like amsr_L2b_stats.f and byte_maps.f) uses it either.  Rather the call to
!    check_min_rfi is done which mostly redoes what is done here.  Also this routine is too low level
!    to be worrying about the changing configuration for geostationary satellites.
!    Also celrfi contains a very limited amount of information.
!    So I decided to make the following changes
!    1.  The computation of celrfi is remove as well as everything need to do the computation
!    2.  Now computing and storing reflat and reflon

!    feb 3 2013 version changed 8/29/2013.
!    jscan1 is now 2*iscan-3 to make 89b adjacent and before to 89a
!    with this change jscan=2*numscan-1 does not have lats and lons.  it also does not have ta (-999).


!    5/15/2010 version changed on 12/17/2010. additional qc check of abs(omega_avg-240.).gt.10 inserted
!    i dont know if this has every happened but this qc is used in sips processing so i thought it best to include here 
!    also rfi flagging is improved by finding the minimum value of rfiang-rfimap by cycling over all 8 sats


!    same as amsr_geolocation data 4/29/2007 except that:
!    1.  it has been updated based on ssmi module
!    2.  functionally, rp is now sligthly different and comes from amsr_module
!    3.  computation of subpoint lat,lon,alt is every so slightly dif.  this version uses the updated method
!    4.  zang and the alpha,beta,kflag for the sun angles are computed
!    5.  the geoloc parameters that were specified by data statements and parameter statements here in are now
!        read in and passed to these routines via l2_module
!    6.  now checks 4 geostationary sats

!    nov 29 2006 version changed on april 29 2007.  amsrj spillover_time change from -0.0880458d0 to -0.0858375d0 to
!    compenstate change to azoffset_l0.  this will have neglibable effect and is just done as a formality.
 
 
!    december 12, 2005, 02:35:16 pm version changed on nov 29, 2006.  old version stored as amsr_geolocation_preamsrj.f.
!    this change only affects the iamsr=2 option.  it has no effect on amsr-e
!    this new version is based on the nov 2006 geoloc analysis for amsr-j done in folder o:\amsr_l2\midori2\geoloc
 
!    9/28/2005 version changed on dec 12 2005.  when rpy is changing too fast for midori-2 i set the
!    quality flag to 6 rather than 8, and i go ahead with the geolocation processing.  midori-2 has yet to be processed
!    so no reprocessing is necessary
 
module Geolocation
    use NAN_support, only: nan_f32, nan_f64
contains
    subroutine amsr_geolocation

        use l2_module  !everything!
        use pre_nut_routines, only: get_gm_angle, fd_precession, fd_nutation
        use math_routines, only:  cross_norm,dot_product_unit8,invert_3by3,fixang4,fixang8,minang4,minang8
        use sun_moon, only: sun_vector,moon_vector
        
        implicit none
    
        integer(4), parameter:: imode=1      ! 1 means only due precession; otherwise do both precession and nutation
        real(8), parameter :: alphax=0  !dummy value for alpha
        real(8), parameter :: betax =0  !dummy value for beta
    
        !real(8) dot_product_unit8
        real(8) high_res_step_time,low_res_step_time
        real(8) x_eci2000(3),y_eci2000(3),z_eci2000(3),x_eci(3),y_eci(3),z_eci(3),alpha,r3,v3
        real(8) scpos_eci2000(3),scpos_eci(3),scvel_eci(3),rpy(3),rotangle,days
        real(8) scpos0(3),scvel0(3),sunvec_eci2000(3),sunvec(3),sundis_km,moonvec_eci2000(3),moonvec(3),moondis_km
        real(8) observation_time, time_since_scan_start
        real(8) pitch_corr
        real(8) x0(3),y0(3),z0(3)
        real(8) precession(3,3),nutation(3,3),np(3,3)
        real(8) b(3),refl(3),xlat,xlon_rot,range,rearth,thtinc,azim,pra,suninc,sunazm,sunglt
        real(8) sc_to_sun(3),dot1,dot2,cossunang
        real(8) geolat,geolon
        real(4) omega_avg,asum
    
        integer(4) iscan,jcel,icel,kscan,iasum,ibad
        integer(4) ierror,jflag
        integer(4) jscan1,jscan2
 
        data high_res_step_time/1.3d-3/
        data low_res_step_time/2.6d-3/
        
        pitch_corr=0  !no pitch correction for amsr2
        
        cellat_89=0  !i put the default here because of the jscan1,jscan2 complication
        cellon_89=0  !i put the default here because of the jscan1,jscan2 complication
        celtht_89=0
        celphi_89=0
        do 200 iscan=1,numscan
            scan_index(iscan) = iscan
            scpos0=scpos(:,iscan)
            scvel0=scvel(:,iscan)
            rpy=   scrpy(:,iscan)
            z0= -scpos0
            x0=  scvel0
            call cross_norm(z0,x0, y0)
 
            jscan1=2*iscan-3  !89b assign to this to make it adjacent and before to 89a
            jscan2=2*iscan    
 
            !    ====================================================================================================================
            !    ==========================================  find average spin rate omega ===========================================
            !    ====================================================================================================================

            ibad=0; iasum=0; asum=0
            do kscan=iscan-5,iscan+5
                if(kscan.lt.1 .or. kscan.gt.numscan) cycle
                if(iflag_l0(kscan).eq.10) cycle     !no science data
                if(omega(kscan).le.0) then         !bad omega
                    ibad=ibad+1
                else
                    iasum=iasum+1
                    asum=asum+omega(kscan)
                endif
            enddo
    
            if(ibad.gt.1) iasum=0 !if there is more than one bad omega near scan, do not process

            if(iasum.ne.0) then
                omega_avg=asum/iasum
            else
                omega_avg=9999.
            endif

            omega_smooth(iscan) = omega_avg

            !    ====================================================================================================================
            !    ========================================END find average spin rate omega ===========================================
            !    ====================================================================================================================

            if(abs(omega_avg-240.).gt.10 .and. iflag_l0(iscan).eq.0) iflag_l0(iscan)=8

            if(iflag_l0(iscan).ne.0) then
                zang(iscan)=0
                beta_sun(iscan)=0
                alpha_sun(iscan)=0
                kflag_sun(iscan)=0
                sunlat(iscan)=0
                sunlon(iscan)=0
                sundis(iscan)=0
                scloc(:,iscan)=0
                cellat(:,iscan)=0
                cellon(:,iscan)=0
                celtht(:,iscan)=0
                celphi(:,iscan)=0
                celsun(:,iscan)=0
                cellat_cold(iscan)=0
                cellon_cold(iscan)=0
                cycle
            endif
 
            if(maxval(abs(rpy)).gt.rpy_error) iflag_l0(iscan)=6
 
            jflag=0
 
            do 100 jcel=1,maxcel_89
                time_since_scan_start = (jcel-1)*high_res_step_time
                observation_time = ut1_1993(iscan) + time_since_scan_start - 220838400.d0 !convert from begin 93 to begin 2000
                days=observation_time/86400.d0 -0.5d0

                scpos_eci2000 = scpos0 + time_since_scan_start*scvel0
        
                call get_sc_axes_from_rpy(scpos_eci2000,y0,rpy(1),rpy(2),rpy(3) ,x_eci2000,y_eci2000,z_eci2000)
                call pitch_sc(pitch_corr, x_eci2000,y_eci2000,z_eci2000)
                call roll_sc(rolloffset, x_eci2000,y_eci2000,z_eci2000)
                call get_gm_angle(observation_time,days, rotangle)
        
                if(imode.eq.1) then
                    call fd_precession(days, np)
                else
                    call fd_precession(days, precession)
                    call fd_nutation(  days, nutation)
                    np=matmul(nutation,precession)
                endif
        
                scpos_eci=matmul(np,scpos_eci2000)
                x_eci    =matmul(np,x_eci2000)
                y_eci    =matmul(np,y_eci2000)
                z_eci    =matmul(np,z_eci2000)
 
                if(jcel.eq.1) then
 
                    scvel_eci=matmul(np,scvel0)
        
                    !    it is used to compute zang in the eci2000 system (windsat still does i think).  but its better to compute it in the eci system
                    r3=scpos_eci(3)/sqrt(dot_product(scpos_eci,scpos_eci))
                    v3=scvel_eci(3)/sqrt(dot_product(scvel_eci,scvel_eci))
                    zang(iscan)=datan2d(r3,v3) + 90.
                    call fixang8(zang(iscan))
        
                    call  sun_vector(days,  sunvec_eci2000, sundis_km)     !returns vector in j2000 system, dis is km
                    sunvec    =matmul(np,sunvec_eci2000)
        
                    call moon_vector(days, moonvec_eci2000,moondis_km)     !returns vector in j2000 system, dis is km
                    moonvec    =matmul(np,moonvec_eci2000)
        
                    sunlat(iscan)=asind(sunvec(3))
                    sunlon(iscan)=atan2d(sunvec(2),sunvec(1)) - rotangle
                    call fixang4( sunlon(iscan))
                    sundis(iscan)=sundis_km/149.619d6
        
                    sc_to_sun=sundis_km*sunvec - 1.d-3*scpos_eci
                    sc_to_sun=sc_to_sun/sqrt(dot_product(sc_to_sun,sc_to_sun))
 
                    cossunang=-dot_product_unit8(z_eci,sc_to_sun)    !minus sign accounts for dif between windsat x,y,z and amsr x,y,z
                    beta_sun(iscan)=acosd(cossunang)
        
                    dot1= dot_product(x_eci,sc_to_sun)
                    dot2=-dot_product(y_eci,sc_to_sun)      !minus sign accounts for dif between windsat x,y,z and amsr x,y,z
                    alpha_sun(iscan)=atan2d(dot2,dot1)
        
                    !     icase=1 means compute satellite subpoint location
                    !     icase=2 means compute just cell locations given alpha and beta
                    !     icase=3 means compute all  cell parameters given alpha and beta
                    !     icase=4 means compute cell locations given boresight
        
                    ! check to see if SC can see the sun by geolocating the sunvec
                    call fd_cell_parameters(4,re,rp,ffac,geosync_alt, &
                                    rotangle,scpos_eci,x_eci,y_eci,z_eci,betax,alphax,sunvec,sc_to_sun, &
                                    refl,xlat,xlon_rot,range,rearth,thtinc,azim,pra, &
                                    suninc,sunazm,sunglt,geolat,geolon,ierror)

                    kflag_sun(iscan)=ierror
    
                    !find the satellite subpoint.
                    call fd_cell_parameters(1,re,rp,ffac,geosync_alt,&
                                        rotangle,scpos_eci,x_eci,y_eci,z_eci,betax,alphax,sunvec,&
                                        b,refl,xlat,xlon_rot,range,rearth,thtinc,azim,pra,&
                                        suninc,sunazm,sunglt,geolat,geolon,ierror)

                    if(ierror.ne.0) jflag=1
                    scloc(1,iscan)=xlat
                    scloc(2,iscan)=xlon_rot
                    scloc(3,iscan)=range
                    scpos_eci0(:,iscan)=scpos_eci
                    x_eci0(:,iscan)=x_eci
                    y_eci0(:,iscan)=y_eci
                    z_eci0(:,iscan)=z_eci
                    sunvec0(:,iscan)= sundis_km*sunvec
                    moonvec0(:,iscan)=moondis_km*moonvec
                endif
 
                icel=1 + int((jcel-1)/2)
    
                if(2*icel-1.eq.jcel) then  !do low frequency

                    alpha=omega_avg*(time_since_scan_start + start_time_a) + azoffset(1) - 180.
    
                    call fd_cell_parameters(3,re,rp,ffac,geosync_alt,&
                                rotangle,scpos_eci,x_eci,y_eci,z_eci,beta(1),alpha,sunvec,&
                                b,refl,xlat,xlon_rot,range,rearth,thtinc,azim,pra, &
                                suninc,sunazm,sunglt,geolat,geolon,ierror)
    
                    if(ierror.ne.0) jflag=1
                    cellat(icel,iscan)=xlat
                    cellon(icel,iscan)=xlon_rot
                    call fixang4( cellon(icel,iscan))  !needed because of real(8) to real(4) conversion
                    celtht(icel,iscan)=thtinc
                    celphi(icel,iscan)=azim
                    celsun(icel,iscan)=sunglt
                    reflat(icel,iscan)=geolat
                    reflon(icel,iscan)=geolon
                    call fixang4( reflon(icel,iscan))  

                endif
 
                alpha=omega_avg*(time_since_scan_start + start_time_a) + azoffset(2) - 180.

                !changed icase from 2 to 3 so that all cell parameters are computed
                call fd_cell_parameters(3,re,rp,ffac,geosync_alt, &
                                rotangle,scpos_eci,x_eci,y_eci,z_eci,beta(2),alpha,sunvec, &
                                b,refl,xlat,xlon_rot,range,rearth,thtinc,azim,pra,suninc,sunazm,sunglt,geolat,geolon,ierror)
                if(ierror.ne.0) jflag=1
                cellat_89(jcel,jscan2)=xlat
                cellon_89(jcel,jscan2)=xlon_rot
                celtht_89(jcel,jscan2)=thtinc
                celphi_89(jcel,jscan2)=azim
                celrng_89(jcel,jscan2)=range
                call fixang4( cellon_89(jcel,jscan2))  !needed because of real(8) to real(4) conversion
                ! print *, beta(2),thtinc
                ! print *, 'jcel, jscan2, cellat_89(jcel,jscan2), cellon_89(jcel,jscan2)', jcel, jscan2, cellat_89(jcel,jscan2), cellon_89(jcel,jscan2)
                ! print *, 'jcel, jscan2, celtht_89(jcel,jscan2), celphi_89(jcel,jscan2)', jcel, jscan2, celtht_89(jcel,jscan2), celphi_89(jcel,jscan2)
                ! print *, '-----'

                alpha=omega_avg*(time_since_scan_start + start_time_b) + azoffset(3) - 180.
                !changed icase from 2 to 3 so that all cell parameters are computed
                call fd_cell_parameters(3,re,rp,ffac,geosync_alt, &
                                rotangle,scpos_eci,x_eci,y_eci,z_eci,beta(3),alpha,sunvec, &
                                b,refl,xlat,xlon_rot,range,rearth,thtinc,azim,pra, &
                                suninc,sunazm,sunglt,geolat,geolon,ierror)
                if(ierror.ne.0) jflag=1
                if(jscan1.ge.1) then
                    cellat_89(jcel,jscan1)=xlat
                    cellon_89(jcel,jscan1)=xlon_rot
                    celtht_89(jcel,jscan1)=thtinc
                    celphi_89(jcel,jscan1)=azim
                    celrng_89(jcel,jscan1)=range
                    call fixang4( cellon_89(jcel,jscan1))  !needed because of real(8) to real(4) conversion
                    ! print *, beta(3),thtinc
                    ! print *, 'jcel, jscan1, cellat_89(jcel,jscan1), cellon_89(jcel,jscan1)', jcel, jscan1, cellat_89(jcel,jscan1), cellon_89(jcel,jscan1)
                    ! print *, 'jcel, jscan1, celtht_89(jcel,jscan1), celphi_89(jcel,jscan1)', jcel, jscan1, celtht_89(jcel,jscan1), celphi_89(jcel,jscan1)
                    ! print *, '-----'        
                endif
            100 continue
 
            if(jflag.eq.1) then
                if(iflag_l0(iscan).eq.6) then
                    iflag_l0(iscan)=7  !rpy error plus no earth intersection error
                else
                    iflag_l0(iscan)=5  !rpy ok but no earth intersection error occured
                endif
            endif
 
            !    ===========================================================================================
            !    ===================== new block for computing location of cold spillover ==================
            !    ===========================================================================================
    
            time_since_scan_start=spillover_time
            observation_time = ut1_1993(iscan) + time_since_scan_start - 220838400.d0 !convert from begin 93 to begin 2000
            days=observation_time/86400.d0 -0.5d0
        
            scpos_eci2000 = scpos0 + time_since_scan_start*scvel0
 
            call get_sc_axes_from_rpy(scpos_eci2000,y0,rpy(1),rpy(2),rpy(3) ,x_eci2000,y_eci2000,z_eci2000)
        
            call pitch_sc(pitch_corr, x_eci2000,y_eci2000,z_eci2000)
            call  roll_sc(rolloffset, x_eci2000,y_eci2000,z_eci2000)
        
            call get_gm_angle(observation_time,days, rotangle)
        
            if(imode.eq.1) then
                call fd_precession(days, np)
            else
                call fd_precession(days, precession)
                call fd_nutation(  days, nutation)
                np=matmul(nutation,precession)
            endif
        
            scpos_eci=matmul(np,scpos_eci2000)
            x_eci    =matmul(np,x_eci2000)
            y_eci    =matmul(np,y_eci2000)
            z_eci    =matmul(np,z_eci2000)
    
            alpha=omega_avg*(time_since_scan_start + start_time_a) + azoffset(1) - 180.
        
            call fd_cell_parameters(2,re,rp,ffac,geosync_alt, &
                                    rotangle,scpos_eci,x_eci,y_eci,z_eci,beta(1),alpha,sunvec, &
                                    b,refl,xlat,xlon_rot,range,rearth,thtinc,azim,pra,suninc,sunazm,sunglt,geolat,geolon,ierror)
            
            if(ierror.ne.0 .and. iflag_l0(iscan).eq.0) iflag_l0(iscan)=5  !strange error: everything ok but cold location in error
        
            cellat_cold(iscan)=xlat
            cellon_cold(iscan)=xlon_rot
            call fixang4(cellon_cold(iscan))  !needed because of real(8) to real(4) conversion
    
        200 continue
    return
    end subroutine amsr_geolocation

!    x, y, z  s/c vectors
!    x vector points in general direction of s/c velocity
!    z vector points is s/c nadir and points towards earth for normal attitude
!    y vector is given by z cross x
 

!    icase=1 means compute satellite subpoint location
!    icase=2 means compute just cell locations given alpha and beta
!    icase=3 means compute all  cell parameters given alpha and beta
!    icase=4 means compute cell locations given boresight

    subroutine fd_cell_parameters(icase,re,rp,ffac,geosync_alt, &
                                  rotangle,scpos,x,y,z,beta,alpha,sunvec, &
                                  b,refl,xlat,xlon_rot,range,rearth,thtinc,azim,pra,suninc,sunazm,sunglt,geolat,geolon,ierror)
                                  
        use math_routines, only:  cross_norm,dot_product_unit8,invert_3by3,fixang4,fixang8,minang4,minang8
        
        implicit none

        real(8),    intent(in)    :: re,rp,ffac,geosync_alt,rotangle,scpos(3),x(3),y(3),z(3),beta,alpha,sunvec(3)
        integer(4), intent(in)    :: icase
        real(8),    intent(out)   :: refl(3),xlat,xlon_rot,range,rearth,thtinc,azim,pra,suninc,sunazm,sunglt,geolat,geolon
        integer(4), intent(out)   :: ierror
        real(8),    intent(inout) :: b(3)
    
        !real(8) dot_product_unit8
    
        real(8) r,ru(3),sinlat,coslat,costht,xlon
        real(8) sinbeta,cosbeta
        real(8) bx,by,bz,cel1,cel2,cel3
        real(8) coslon,sinlon
        real(8) geoid(3),cosazim,cossunglt
        real(8) deltau
        real(8) ha(3),he(3),ve(3),dot1,dot2
        real(8) rcel(3),earth_to_geo(3),dist
 
        r=sqrt(dot_product(scpos,scpos))
        ru=scpos/r
    
        if(icase.eq.1) then     !find pointing vector b(3) to geoid
            sinlat=ru(3)
            coslat=sqrt(ru(1)*ru(1)+ru(2)*ru(2))
            rearth=rp/sqrt(ffac*coslat**2+sinlat**2)  
            deltau=((r-rearth)/re)*sqrt(1 - (1-ffac)*sinlat**2/((ffac*coslat)**2+sinlat**2))
            b(1)=-ru(1)
            b(2)=-ru(2)
            b(3)=-ru(3)*(1+deltau)/(ffac+deltau)
            b=b/sqrt(dot_product(b,b))    
        endif
 
        if(icase.eq.2 .or. icase.eq.3) then     !compute b(3) from alpha and beta
            cosbeta=cosd(beta)
            sinbeta=sind(beta)
            bx= sinbeta*cosd(alpha) !input alpha defined relative to the x axes,  alpha=-90 denotes bx=0, by=1
            by=-sinbeta*sind(alpha) !input alpha defined relative to the x axes , alpha=-90 denotes bx=0, by=1
            bz=cosbeta                
            b=bx*x + by*y + bz*z
        endif
 
 
        !    for icase=4, b(3) is an input
    
        call loccel_range(r,ru(1),ru(2),ru(3),b(1),b(2),b(3),re,rp, cel1,cel2,cel3,rearth,range,ierror)
 
        if(ierror.ne.0) then
            refl=0
            xlat=0
            xlon_rot=0
            range=0
            rearth=0
            thtinc=0
            azim=0
            pra=0
            suninc=0
            sunazm=0
            sunglt=0
            return
        endif

        if(icase.eq.4) return  !case 4 is to see if s/c can see sun, so there is nothing else to computed
 
        coslat=sqrt(cel1*cel1+cel2*cel2)     !cosine of geocentric lat
        xlat=atand(cel3/(ffac*coslat))       !geodetic cell latitude
        xlon=atan2d(cel2,cel1)
        call fixang8( xlon)

        xlon_rot=xlon - rotangle    !convert from inertial to earth
        call fixang8( xlon_rot)

        if(icase.ne.3) return  !nothing else to do for case 1,2 and 4
 
 
        !=================================================================================================================
        !======= compute earth inc. ang. (thtinc) and earth azimuth angle wrt. north (azim) ==============================
        !=================================================================================================================
 
        coslon=cosd(xlon)
        sinlon=sind(xlon)
        coslat=cosd(xlat)
        geoid(1)=coslat*coslon
        geoid(2)=coslat*sinlon
        geoid(3)=sind(xlat)
 
        costht=-dot_product_unit8(b,geoid)
        thtinc=acosd(costht)
 
        if(abs(costht).gt.0.9999999999d0) then    !too close to nadir to get azim
            azim=0
        else
            cosazim=(-sinlon*b(1)+coslon*b(2))/sqrt(1.-costht*costht)
            if(cosazim.lt.-1.) cosazim=-1.
            if(cosazim.gt. 1.) cosazim= 1.
            azim=acosd(cosazim)
            if(b(3)+costht*geoid(3).lt.0.) azim=-azim
            azim=90.-azim     ! convert to relative to clockwise from north
            if(azim.lt.  0.) azim=azim+360.
            if(azim.gt.360.) azim=azim-360.
        endif
 
        !    =================================================================================================================
        !    ======= compute sun glint angle (sunglt), sun inc. ang. (suninc) and sun azimuth angle (sunazm) =================
        !    =================================================================================================================
 
        refl=b + 2*costht*geoid
 
        cossunglt=dot_product_unit8(refl,sunvec)
        sunglt=acosd(cossunglt)
 
        !    compute sun tht and azimuth angles
        costht=dot_product_unit8(sunvec,geoid)
        suninc=acosd(costht)
 
        if(abs(costht).gt.0.9999999999d0) then    !too close to nadir to get sunazm
            sunazm=0
        else
            cosazim=(-sinlon*sunvec(1)+coslon*sunvec(2))/sqrt(1.-costht*costht)
            if(cosazim.lt.-1.) cosazim=-1.
            if(cosazim.gt. 1.) cosazim= 1.
            sunazm=acosd(cosazim)
            if(sunvec(3)-costht*geoid(3).lt.0.) sunazm=-sunazm     !minus sign used here becuase sun is pointing away from surface
            sunazm=90-sunazm     ! convert to relative to clockwise from north
            call fixang8(sunazm)
        endif

        !    =================================================================================================================
        !    ============================= compute polarization rotation angle (pra) =========================================
        !    =================================================================================================================
        call cross_norm(b,z, ha)      !ha=kx(-z)=(-b)x(-z)=bxz
 
        if(thtinc.gt.0.01) then
            call cross_norm(geoid,b, he) !he=kxg=(-b)xg=gxb
        else
            he=ha
        endif
 
        call cross_norm(b,he, ve)     !ve=hexk=bxhe
    
        dot1=dot_product(ve,ha)
        dot2=dot_product(he,ha)
        pra=atan2d(dot1,dot2)
    

        !    ========================================================================================================================
        !    ==================================== 9/20/2013 update: geo lat/lon are computed ========================================
        !    ========================================================================================================================
      
        rcel(1)=cel1
        rcel(2)=cel2
        rcel(3)=cel3
        
        dot1=dot_product(rcel,refl)
        
        dist=-rearth*dot1 + sqrt((rearth*dot1)**2 + (geosync_alt**2 - rearth**2))
        
        earth_to_geo=rearth*rcel + dist*refl
        
        dist=sqrt(earth_to_geo(1)**2 +earth_to_geo(2)**2 +earth_to_geo(3)**2)
        
        if(abs(dist-geosync_alt).gt.10.) stop 'earth_to_geo in error, pgm stopped' !exceeds 10 meter
        
        earth_to_geo=earth_to_geo/dist

        geolat=asind(earth_to_geo(3))
        geolon=atan2d(earth_to_geo(2),earth_to_geo(1)) - rotangle

        return
    end subroutine fd_cell_parameters
 
 
    subroutine loccel_range(r,r1,r2,r3,s1,s2,s3,re,rp, cel1,cel2,cel3,rearth,range,ierror)
        implicit none
    
        real(8), intent(in)  :: r,r1,r2,r3,s1,s2,s3
        real(8), intent(out) :: cel1,cel2,cel3,rearth,range
        real(8) delta,a,cosx,b,c,celmag,arg
        real(8) re,rp
        integer(4) ierror
    
        delta=(re/rp)*(re/rp)-1.
        a=1.+delta*s3*s3
        cosx=-r1*s1-r2*s2-r3*s3
        b=cosx-delta*r3*s3
        c=1.-(re/r)*(re/r)+delta*r3*r3
    
        arg=b*b-a*c
        if(arg.lt.0) then !vector does not intercept earth
            ierror=1
            return
        else
            ierror=0
        endif
 
        range=(b-sqrt(arg))/a
        cel1=r1+range*s1
        cel2=r2+range*s2
        cel3=r3+range*s3
        celmag=sqrt(cel1*cel1+cel2*cel2+cel3*cel3)
        cel1=cel1/celmag
        cel2=cel2/celmag
        cel3=cel3/celmag
        rearth=r*celmag
        range=range*r
        if(range.le.0) ierror=1
        return
    end subroutine loccel_range
 
    subroutine get_sc_axes_from_rpy(scpos,y0,roll,pitch,yaw, x,y,z)

        use math_routines, only:  cross_norm,dot_product_unit8,invert_3by3,fixang4,fixang8,minang4,minang8
        
        implicit none
    
        !arguments
        real(8), intent(in):: scpos(3), y0(3), roll,pitch,yaw
        real(8), intent(out)::  x(3), y(3), z(3)
        
        !local variables
        real(8) x0(3), z0(3)
        real(8) att_mat(3,3),a0(3,3), a(3,3)
        real(8) cos_roll,sin_roll,cos_pitch,sin_pitch,cos_yaw,sin_yaw
    
                
    !    x0,y0,z0 are s/c axis for 0 roll,pitch,yaw
        z0=-scpos/sqrt(dot_product(scpos,scpos))
        call cross_norm(y0,z0, x0)
    
        a0(1,:)=x0
        a0(2,:)=y0
        a0(3,:)=z0
    
        cos_roll= cosd(roll)
        sin_roll= sind(roll)
        cos_pitch=cosd(pitch)
        sin_pitch=sind(pitch)
        cos_yaw=  cosd(yaw)
        sin_yaw=  sind(yaw)
    
    !    following is the 'rpy' or  '123' convention
 
        att_mat(1,1)= cos_pitch*cos_yaw    
        att_mat(1,2)= sin_roll*sin_pitch*cos_yaw + cos_roll*sin_yaw    
        att_mat(1,3)=-cos_roll*sin_pitch*cos_yaw + sin_roll*sin_yaw    
        att_mat(2,1)=-cos_pitch*sin_yaw    
        att_mat(2,2)=-sin_roll*sin_pitch*sin_yaw + cos_roll*cos_yaw    
        att_mat(2,3)= cos_roll*sin_pitch*sin_yaw + sin_roll*cos_yaw    
        att_mat(3,1)= sin_pitch    
        att_mat(3,2)=-sin_roll*cos_pitch    
        att_mat(3,3)= cos_roll*cos_pitch
    
        a=matmul(att_mat,a0)
    
        x=a(1,:)
        y=a(2,:)
        z=a(3,:)
    
        return
    end subroutine get_sc_axes_from_rpy

    subroutine get_att_from_rpy(roll,pitch,yaw,att_mat)

        use math_routines, only:  cross_norm,dot_product_unit8,invert_3by3,fixang4,fixang8,minang4,minang8
        
        implicit none
    
        !arguments
        real(8), intent(in)::  roll,pitch,yaw
        real(8), intent(out)::  att_mat(3,3)
        
        !local variables
        real(8) cos_roll,sin_roll,cos_pitch,sin_pitch,cos_yaw,sin_yaw
    
        cos_roll= cosd(roll)
        sin_roll= sind(roll)
        cos_pitch=cosd(pitch)
        sin_pitch=sind(pitch)
        cos_yaw=  cosd(yaw)
        sin_yaw=  sind(yaw)
    
    !    following is the 'rpy' or  '123' convention
 
        att_mat(1,1)= cos_pitch*cos_yaw    
        att_mat(1,2)= sin_roll*sin_pitch*cos_yaw + cos_roll*sin_yaw    
        att_mat(1,3)=-cos_roll*sin_pitch*cos_yaw + sin_roll*sin_yaw    
        att_mat(2,1)=-cos_pitch*sin_yaw    
        att_mat(2,2)=-sin_roll*sin_pitch*sin_yaw + cos_roll*cos_yaw    
        att_mat(2,3)= cos_roll*sin_pitch*sin_yaw + sin_roll*cos_yaw    
        att_mat(3,1)= sin_pitch    
        att_mat(3,2)=-sin_roll*cos_pitch    
        att_mat(3,3)= cos_roll*cos_pitch
 
        return
    end subroutine get_att_from_rpy

    function slerp_rpy(rpy0,rpy1,t)
        use quaternion_routines, only: slerp, rotation_matrix_to_quaternion, quaternion_to_rotation_matrix,att_q_to_rpy

        implicit none 

        real(8),intent(in) :: rpy0(3)
        real(8),intent(in) :: rpy1(3)
        real(8),intent(in) :: t
        real(8)     :: slerp_rpy(3)
        real(8)     :: att0(3,3),att1(3,3)
        real(8)     :: q0(4),q1(4),q_interp(4)

        call get_att_from_rpy(rpy0(1),rpy0(2),rpy0(2),att0)
        call get_att_from_rpy(rpy1(1),rpy1(2),rpy1(2),att1)

        q0 = rotation_matrix_to_quaternion(att0)
        q1 = rotation_matrix_to_quaternion(att1)

        q_interp = slerp(q0,q1,t)

        call att_q_to_rpy(q_interp,slerp_rpy)

        

    end function slerp_rpy

 
    subroutine pitch_sc(pitch, x,y,z)
 
        implicit none
        
        real(8) x(3), y(3), z(3)
        real(8) att_mat(3,3),a0(3,3), a(3,3)
        real(8) pitch,cos_pitch,sin_pitch
 
        if(pitch.eq.0) return
        
        a0(1,:)=x
        a0(2,:)=y
        a0(3,:)=z
    
        cos_pitch=cosd(pitch)
        sin_pitch=sind(pitch)
    
        att_mat(1,1)= cos_pitch    
        att_mat(1,2)= 0    
        att_mat(1,3)=-sin_pitch    
        att_mat(2,1)= 0    
        att_mat(2,2)= 1    
        att_mat(2,3)= 0    
        att_mat(3,1)= sin_pitch    
        att_mat(3,2)= 0    
        att_mat(3,3)= cos_pitch
    
        a=matmul(att_mat,a0)
    
        x=a(1,:)
        y=a(2,:)
        z=a(3,:)
    
        return
    end subroutine pitch_sc
 
 
    subroutine roll_sc(roll, x,y,z)
 
    implicit none
    
        real(8) x(3), y(3), z(3)
        real(8) att_mat(3,3),a0(3,3), a(3,3)
        real(8) roll,cos_roll,sin_roll
        
        if(roll.eq.0) return
    
        a0(1,:)=x
        a0(2,:)=y
        a0(3,:)=z
    
        cos_roll=cosd(roll)
        sin_roll=sind(roll)
    
        att_mat(1,1)= 1    
        att_mat(1,2)= 0    
        att_mat(1,3)= 0    
        att_mat(2,1)= 0
        att_mat(2,2)= cos_roll    
        att_mat(2,3)= sin_roll    
        att_mat(3,1)= 0
        att_mat(3,2)= -sin_roll    
        att_mat(3,3)= cos_roll
    
        a=matmul(att_mat,a0)
    
        x=a(1,:)
        y=a(2,:)
        z=a(3,:)
    
        return
    end subroutine roll_sc

    subroutine amsr2_geolocation_extra

        use l2_module  !everything!
        use l2_module_extra
        use pre_nut_routines, only: get_gm_angle, fd_precession, fd_nutation
        use math_routines, only:  cross_norm,dot_product_unit8,invert_3by3,fixang4,fixang8,minang4,minang8

        !eventually won't need this:
        use sun_moon, only: sun_vector,moon_vector

        implicit none
    
        integer(4), parameter:: imode=1      ! 1 means only due precession; otherwise do both precession and nutation
        real(8), parameter :: alphax=0  !dummy value for alpha
        real(8), parameter :: betax =0  !dummy value for beta
    
        real(8) high_res_step_time,low_res_step_time
        real(8) x_eci2000(3),y_eci2000(3),z_eci2000(3),x_eci(3),y_eci(3),z_eci(3),alpha,r3,v3
        real(8) scpos_eci2000(3),scpos_eci(3),scvel_eci(3),rpy(3),rotangle,days
        real(8) scpos0(3),scvel0(3),sunvec_eci2000(3),sunvec(3),sundis_km,moonvec_eci2000(3),moonvec(3),moondis_km
        real(8) scpos_before(3),scvel_before(3),scpos_after(3),scvel_after(3),rpy_before(3),rpy_after(3)

        real(8) observation_time, time_since_scan_start
        real(8) pitch_corr
        real(8) x0(3),y0(3),z0(3)
        real(8) precession(3,3),nutation(3,3),np(3,3)
        real(8) b(3),refl(3),xlat,xlon_rot,range,rearth,thtinc,azim,pra,suninc,sunazm,sunglt
        real(8) sc_to_sun(3),dot1,dot2,cossunang
        real(8) geolat,geolon
        real(4) omega_avg,asum
        real(8) scan_period,delta_time_between_subscans,kscan_delta_time
        real(8) t_norm,norm_after,norm_before,avg_norm,avg_spd
        real(8) step_time_extra
    
        integer(4) iscan,jcel,kscan !,iasum,ibad
        integer(4) ierror,jflag
        integer(4) kscan_start,kscan_end
        integer(4) output_scan_number
 
        !everything in the extra geolocation is based on the low res sample time
        data low_res_step_time/2.6d-3/
        data high_res_step_time/1.3d-3/

        cellat_extra = nan_f32
        cellon_extra = nan_f32
        celtht_extra = nan_f32
        celphi_extra = nan_f32
        scan_index_extra = nan_f64
        omega_extra = nan_f64
        scan_time_extra = nan_f64
        zang_extra = nan_f64

        step_time_extra = low_res_step_time/(1+num_extra_fovs)

        pitch_corr=0  !no pitch correction for amsr2
        
        kscan_start = 0
        kscan_end = 0
        if (num_extra_scans > 0) then
            kscan_start = -1*floor(num_extra_scans/2.0)
            kscan_end = floor((1+num_extra_scans)/2.0)
        endif
        !print *,kscan_start,kscan_end
        scan_period = 60.0/omega_avg

        scans: do iscan=2,numscan-1
          ! Use average spin rate from native L2 data
          omega_avg = omega_smooth(iscan)
          !print *,'omega_avg = ',omega_avg
          if(abs(omega_avg-240.).gt.10 .and. iflag_l0(iscan).eq.0) iflag_l0(iscan)=8
          scan_period = 60.0/omega_avg   ! in seconds
          delta_time_between_subscans = scan_period/(1+num_extra_scans)

          kscans: do kscan=kscan_start,kscan_end
            output_scan_number = ((iscan-1)*(1+num_extra_scans))+kscan-kscan_start+1
            scan_index_extra(output_scan_number) = iscan+kscan*(1.0/(1+num_extra_scans))
            scan_time_extra(output_scan_number) = scan_time(iscan)+kscan*delta_time_between_subscans
            !print *,'output_scan_number:',output_scan_number
            !---------------------------------------------------------------------------------
            ! Interpolate the SC position, velocity and attitude to the virtual scan if needed
            !---------------------------------------------------------------------------------
             if (kscan == 0) then
                scpos0=scpos(:,iscan)
                scvel0=scvel(:,iscan)
                rpy=   scrpy(:,iscan)
                kscan_delta_time = 0.0D0
            else 

                kscan_delta_time = kscan*delta_time_between_subscans
                if (kscan < 0) then
                    scpos_before = scpos(:,iscan-1)
                    scvel_before = scvel(:,iscan-1)
                    scpos_after  = scpos(:,iscan)
                    scvel_after  = scvel(:,iscan)
                    rpy_after    = scrpy(:,iscan)
                    rpy_before   = scrpy(:,iscan-1)
                    t_norm       = 1.0 + kscan_delta_time/scan_period !remember - kscan_delta_time is negative
                    norm_before  = norm2(scpos_before)
                    norm_after   = norm2(scpos_after)
                else
                    scpos_before = scpos(:,iscan)
                    scvel_before = scvel(:,iscan)
                    scpos_after  = scpos(:,iscan+1)
                    scvel_after  = scvel(:,iscan+1)
                    rpy_after    = scrpy(:,iscan+1)
                    rpy_before   = scrpy(:,iscan)
                    t_norm       = (kscan_delta_time/scan_period)
                endif
                ! print *, t_norm
                avg_norm  = (norm2(scpos_before)+norm2(scpos_after))/2.0
                scpos0 = (1.0-t_norm)*scpos_before + t_norm*scpos_after
                !correct the distance from earth center
                ! print *, avg_norm,norm2(scpos0)

                scpos0 = scpos0*avg_norm/norm2(scpos0)

                avg_spd = (norm2(scvel_before)+norm2(scvel_after))/2.0
                scvel0 = (1.0-t_norm)*scvel_before + t_norm*scvel_after
                ! correct for spped change due to interpolation
                ! print *, avg_spd,norm2(scvel0)
                scvel0 = scvel0*avg_spd/norm2(scvel0)
                ! print *,scvel_before,scvel_after,scvel0
                !interpolate the rpy
                rpy = slerp_rpy(rpy_before,rpy_after,t_norm)
                ! print *,rpy_before,rpy_after,rpy
                ! print *
            endif
            !---------------------------------------------------------------------------------
            ! END: Interpolate the SC position, velocity and attitude to the virtual scan if needed
            !---------------------------------------------------------------------------------
            
            z0= -scpos0
            x0=  scvel0
            call cross_norm(z0,x0, y0)
 
            if(iflag_l0(iscan).ne.0) then
                zang(iscan)=0
                beta_sun(iscan)=0
                alpha_sun(iscan)=0
                kflag_sun(iscan)=0
                sunlat(iscan)=0
                sunlon(iscan)=0
                sundis(iscan)=0
                scloc(:,iscan)=0
                cellat(:,iscan)=0
                cellon(:,iscan)=0
                celtht(:,iscan)=0
                celphi(:,iscan)=0
                celsun(:,iscan)=0
                cellat_cold(iscan)=0
                cellon_cold(iscan)=0
                cycle
            endif
 
            if(maxval(abs(rpy)).gt.rpy_error) iflag_l0(iscan)=6
 
            jflag=0
 
            fovs : do jcel=1,maxcel_extra
                time_since_scan_start = (jcel-1)*step_time_extra
                ! print *, time_since_scan_start
                observation_time = ut1_1993(iscan) + time_since_scan_start + kscan_delta_time - 220838400.d0 !convert from begin 93 to begin 2000
                days=observation_time/86400.d0 -0.5d0

                scpos_eci2000 = scpos0 + time_since_scan_start*scvel0
        
                call get_sc_axes_from_rpy(scpos_eci2000,y0,rpy(1),rpy(2),rpy(3) ,x_eci2000,y_eci2000,z_eci2000)
                call pitch_sc(pitch_corr, x_eci2000,y_eci2000,z_eci2000)
                call roll_sc(rolloffset, x_eci2000,y_eci2000,z_eci2000)
                call get_gm_angle(observation_time,days, rotangle)

                if(imode.eq.1) then
                    call fd_precession(days, np)
                else
                    call fd_precession(days, precession)
                    call fd_nutation(  days, nutation)
                    np=matmul(nutation,precession)
                endif
        
                scpos_eci=matmul(np,scpos_eci2000)
                x_eci    =matmul(np,x_eci2000)
                y_eci    =matmul(np,y_eci2000)
                z_eci    =matmul(np,z_eci2000)
 
                if(jcel.eq.1) then
                    scvel_eci=matmul(np,scvel0)
                    ! it is used to compute zang in the eci2000 system (windsat still does i think).  but its better to compute it in the eci system
                    r3=scpos_eci(3)/sqrt(dot_product(scpos_eci,scpos_eci))
                    v3=scvel_eci(3)/sqrt(dot_product(scvel_eci,scvel_eci))

                    zang_extra(output_scan_number)=datan2d(r3,v3) + 90.
                    call fixang8(zang_extra(output_scan_number))
        
                    call  sun_vector(days,  sunvec_eci2000, sundis_km)     !returns vector in j2000 system, dis is km
                    sunvec=matmul(np,sunvec_eci2000)
        
                    call moon_vector(days, moonvec_eci2000,moondis_km)     !returns vector in j2000 system, dis is km
                    moonvec    =matmul(np,moonvec_eci2000)
        
                    sunlat_extra(output_scan_number)=asind(sunvec(3))
                    sunlon_extra(output_scan_number)=atan2d(sunvec(2),sunvec(1)) - rotangle
                    call fixang4(sunlon_extra(output_scan_number))
                    sundis_extra(output_scan_number)=sundis_km/149.619d6
        
                    sc_to_sun=sundis_km*sunvec - 1.d-3*scpos_eci
                    sc_to_sun=sc_to_sun/sqrt(dot_product(sc_to_sun,sc_to_sun))
 
                    cossunang=-dot_product_unit8(z_eci,sc_to_sun)    !minus sign accounts for dif between windsat x,y,z and amsr x,y,z
                    beta_sun_extra(output_scan_number) = acosd(cossunang)
        
                    dot1= dot_product(x_eci,sc_to_sun)
                    dot2=-dot_product(y_eci,sc_to_sun)      !minus sign accounts for dif between windsat x,y,z and amsr x,y,z
                    alpha_sun_extra(output_scan_number) = atan2d(dot2,dot1)
        
                    !     icase=1 means compute satellite subpoint location
                    !     icase=2 means compute just cell locations given alpha and beta
                    !     icase=3 means compute all  cell parameters given alpha and beta
                    !     icase=4 means compute cell locations given boresight
        
                    ! check to see if SC can see the sun by geolocating the sunvec
                    call fd_cell_parameters(4,re,rp,ffac,geosync_alt, &
                                    rotangle,scpos_eci,x_eci,y_eci,z_eci,betax,alphax,sunvec,sc_to_sun, &
                                    refl,xlat,xlon_rot,range,rearth,thtinc,azim,pra, &
                                    suninc,sunazm,sunglt,geolat,geolon,ierror)

                    kflag_sun_extra(output_scan_number) = ierror

                    !find the satellite subpoint.
                    call fd_cell_parameters(1,re,rp,ffac,geosync_alt,&
                                        rotangle,scpos_eci,x_eci,y_eci,z_eci,betax,alphax,sunvec,&
                                        b,refl,xlat,xlon_rot,range,rearth,thtinc,azim,pra,&
                                        suninc,sunazm,sunglt,geolat,geolon,ierror)

                    if(ierror.ne.0) jflag=1
                    scloc_extra(1,output_scan_number)=xlat
                    scloc_extra(2,output_scan_number)=xlon_rot
                    scloc_extra(3,output_scan_number)=range
                    
                    scpos_eci0_extra(:,output_scan_number)=scpos_eci
                    x_eci0_extra(:,output_scan_number)=x_eci
                    y_eci0_extra(:,output_scan_number)=y_eci
                    z_eci0_extra(:,output_scan_number)=z_eci
                    sunvec0_extra(:,output_scan_number)= sundis_km*sunvec
                    moonvec0_extra(:,output_scan_number)=moondis_km*moonvec
                endif

                alpha=omega_avg*(time_since_scan_start + start_time_a) + azoffset(1) - 180.
    
                call fd_cell_parameters(3,re,rp,ffac,geosync_alt,&
                            rotangle,scpos_eci,x_eci,y_eci,z_eci,beta(1),alpha,sunvec,&
                            b,refl,xlat,xlon_rot,range,rearth,thtinc,azim,pra, &
                            suninc,sunazm,sunglt,geolat,geolon,ierror)
    
                !print *,xlat,xlon_rot
                if(ierror.ne.0) jflag=1
                cellat_extra(jcel,output_scan_number)=xlat
                cellon_extra(jcel,output_scan_number)=xlon_rot
                celrange_extra(jcel,output_scan_number) = range
                call fixang4(cellon_extra(jcel,output_scan_number))  !needed because of real(8) to real(4) conversion
                celtht_extra(jcel,output_scan_number)=thtinc
                celphi_extra(jcel,output_scan_number)=azim
                celsun_extra(jcel,output_scan_number)=sunglt
                reflat_extra(jcel,output_scan_number)=geolat
                reflon_extra(jcel,output_scan_number)=geolon
                call fixang4(reflon_extra(jcel,output_scan_number))  
            enddo fovs
          enddo kscans
        enddo scans
    return
    end subroutine amsr2_geolocation_extra
end module Geolocation
