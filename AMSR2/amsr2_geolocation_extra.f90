module xGeolocation_extra
    
    private
    public amsr2_geolocation_extra

contains

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
            print *,'output_scan_number:',output_scan_number
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

                    kflag_sun(output_scan_number) = ierror
                    
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
module xGeolocation_extra