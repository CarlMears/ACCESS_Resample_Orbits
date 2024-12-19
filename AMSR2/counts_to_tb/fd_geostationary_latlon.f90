    subroutine fd_geostationary_latlon(iscan,ifreq, geolat,geolon)
        use l2_module
        use pre_nut_routines, only: get_gm_angle
        use math_routines, only: fixang4
        implicit none

        integer(4) iscan,ifreq,icel
        real(4) geolat,geolon

        real(8) b(3),cosbeta,sinbeta,bx,by,bz

        real(8) observation_time,rotangle,days,dotx,dist
        real(8) r(3),earth_to_geo(3)
        real(8) amag
      
        !     use middle cell for specifying geolat,geolon
        icel=8
        if(ifreq.ge.7) icel=16
      
        observation_time = ut1_1993(iscan) + spillover_time - 220838400.d0 !convert from begin 93 to begin 2000
        days=observation_time/86400.d0 -0.5d0
        call get_gm_angle(observation_time,days, rotangle)
    
        r=scpos_eci0(:,iscan)
        amag=sqrt(r(1)**2 + r(2)**2 + r(3)**2)

        cosbeta=cosd(beta_cold(icel,ifreq))                                                                      
        sinbeta=sind(beta_cold(icel,ifreq))                                                                      
        bx= sinbeta*cosd(alpha_cold(icel,ifreq)) !input alpha defined relative to the -z axes                                       
        by=-sinbeta*sind(alpha_cold(icel,ifreq)) !input alpha defined relative to the -z axes                                       
        bz=cosbeta                 
        b(1)=bx*x_eci0(1,iscan)+by*y_eci0(1,iscan)+bz*z_eci0(1,iscan)                            
        b(2)=bx*x_eci0(2,iscan)+by*y_eci0(2,iscan)+bz*z_eci0(2,iscan)                            
        b(3)=bx*x_eci0(3,iscan)+by*y_eci0(3,iscan)+bz*z_eci0(3,iscan)

        dotx=dot_product(r,b)/amag
      
        dist=(-2*amag*dotx  +sqrt((2*amag*dotx)**2 + 4*(geosync_alt**2 - amag**2)))/2.
        
        earth_to_geo=r + dist*b
      
        if(abs(sqrt(dot_product(earth_to_geo,earth_to_geo))-geosync_alt).gt.10.) stop 'earth_to_geo in error, pgm stopped' !exceeds 10 meter
      
        earth_to_geo=earth_to_geo/sqrt(earth_to_geo(1)**2 +earth_to_geo(2)**2 +earth_to_geo(3)**2)

        geolat=asind(earth_to_geo(3))
        geolon=atan2d(earth_to_geo(2),earth_to_geo(1)) - rotangle
        call fixang4( geolon)
    
        return
    end
