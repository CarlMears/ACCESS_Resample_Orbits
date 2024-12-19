!     july 2 2005 version changed 2/11/2014.  the routine was added.
!     also converted to lower case
!     also    ut1sec93 in subroutine get_gm_angle(ut1sec93,days, angle) renamed ut1sec
!     this renaming is to make it clear that ut1sec is just used to find seconds in day
!     but remember the input days must be from ut12000 


!     june 15, 2005, 01:30:17 pm version change on july 2, 2005.  only comments were changed regarding utc vs ut1.
!     no functional changes were made.
!     
module pre_nut_routines

    use math_routines, only: fixang8
    private
    public fd_precession,fd_nutation,get_gm_angle,fd_gast

contains
    subroutine fd_precession(days, precession)
        implicit none
        
        real(8) days,t,t2,t3,chi,zzz,tht,cosc,cosz,cost,sinc,sinz,sint,precession(3,3)

        !     table 3.211.1 with capital t equal 0

        !    t=(utcsec93/86400.d0 + dif_epoch)/dayspercentury 
        t=days/36525.d0
        t2=t*t
        t3=t*t2
        chi=(2306.2181d0*t + 0.30188d0*t2 + 0.017998d0*t3)/3600.d0    
        zzz=(2306.2181d0*t + 1.09468d0*t2 + 0.018203d0*t3)/3600.d0
        tht=(2004.3109d0*t - 0.42665d0*t2 - 0.041833d0*t3)/3600.d0

        cosc=cosd(chi)
        cosz=cosd(zzz)
        cost=cosd(tht)
        sinc=sind(chi)
        sinz=sind(zzz)
        sint=sind(tht)

        !    eq.  3.21-8
        precession(1,1)= cosz*cost*cosc - sinz*sinc
        precession(1,2)=-cosz*cost*sinc - sinz*cosc
        precession(1,3)=-cosz*sint
        precession(2,1)= sinz*cost*cosc + cosz*sinc
        precession(2,2)=-sinz*cost*sinc + cosz*cosc
        precession(2,3)=-sinz*sint
        precession(3,1)= sint*cosc
        precession(3,2)=-sint*sinc
        precession(3,3)= cost

        return
    end subroutine fd_precession


    !note that this routine does one additioal rotation to put the x axis at the mean equinox
    subroutine fd_nutation(days, nutation)
        implicit none
        
        real(8) days,t,t2,t3,e0,e,dpsi,deps,cosp,cose,cos0,sinp,sine,sin0,nutation(3,3)
        real(8) r_angle,cosr,sinr 
        real(8) n11,n12,n13,n21,n22,n23,n31,n32,n33

        t=days/36525.d0
        t2=t*t
        t3=t*t2

        !     eq 3.222-1
        e0=(84381.448d0 - 46.8150d0*t -0.00059d0*t2 + 0.001813d0*t3)/3600.d0  !deg

        !     eq. 3.225-4
        dpsi=-0.0048d0*sind(125.d0-0.05295d0*days) - 0.0004d0*sind(200.9d0 + 1.97129d0*days)   !deg
        deps= 0.0026d0*cosd(125.d0-0.05295d0*days) + 0.0002d0*cosd(200.9d0 + 1.97129d0*days)   !deg

!      call iau_nut00b(2451545.d0,days, dpsi, deps)
!      dpsi=dpsi*180.d0/3.141592654d0
!      deps=deps*180.d0/3.141592654d0

        e=e0 + deps

        cos0=cosd(e0)
        cose=cosd(e)
        cosp=cosd(dpsi)
        sin0=sind(e0)
        sine=sind(e)
        sinp=sind(dpsi)

!    eq.  3.222-4
        n11= cosp 
        n12=-sinp*cos0
        n13=-sinp*sin0
        n21= sinp*cose
        n22= cosp*cose*cos0 + sine*sin0
        n23= cosp*cose*sin0 - sine*cos0
        n31= sinp*sine
        n32= cosp*sine*cos0 - cose*sin0
        n33= cosp*sine*sin0 + cose*cos0

        r_angle=dpsi*cose

        cosr=cosd(r_angle)
        sinr=sind(r_angle)

        nutation(1,1)= cosr*n11 + sinr*n21
        nutation(1,2)= cosr*n12 + sinr*n22
        nutation(1,3)= cosr*n13 + sinr*n23

        nutation(2,1)=-sinr*n11 + cosr*n21
        nutation(2,2)=-sinr*n12 + cosr*n22
        nutation(2,3)=-sinr*n13 + cosr*n23
    
        nutation(3,1)= n31
        nutation(3,2)= n32
        nutation(3,3)= n33
        return
    end subroutine fd_nutation



    !     ut1sec is just used to find seconds in day.  it can be reference to any years
    !     days should be refereced to jan 1 2000 12z
    subroutine get_gm_angle(ut1sec,days, angle)
        implicit none

        real(8), intent(in)::ut1sec,days
        real(8), intent(out)::angle

        real(8), parameter:: dayspercentury = 36525.d0
        real(8), parameter:: c0 = 24110.54841d0, c1 = 8640184.812866d0, c2 = 0.093104d0, c3 = -6.2d-6
        real(8)  t,t2,t3,ut1

        t=days/36525.d0
        t2=t*t
        t3=t*t2

        ut1=ut1sec - 86400.d0*floor(ut1sec/86400.d0)
        angle = c0 + c1*t  + c2*t2 + c3*t3 +  ut1 !seconds
        angle = angle/240.d0       !convert seconds to degrees
        angle=  angle-360.d0*floor(angle/360.d0)
        return
    end subroutine get_gm_angle 

!     this code converts apparent greenwich angle to mean greenwich angle
!     gast = gmst + eqeq    
    subroutine fd_gast(days,gmst, gast)
        implicit none
        
        real(8) days
        real(8) gast, gmst
        real(8) e,l,omega,dphi,eqeq

        e=23.4393d0 - 0.0000004d0*days
        l=280.47d0 + 0.98565d0*days
        omega=125.04d0 - 0.052954d0*days
        dphi=-0.000319d0*sind(omega) -0.000024d0*sind(2*l)
        eqeq=dphi*cosd(e)  !hours
        gast = gmst + 360.d0*eqeq/24.d0    !degrees
        call fixang8( gmst)
        return
    end subroutine fd_gast
end module pre_nut_routines