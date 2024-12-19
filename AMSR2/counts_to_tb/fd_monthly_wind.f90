    subroutine fd_monthly_wind(lyear,idayjl,xhour,xlat,xlon, win)
        
        use open_file_routines, only: openbig_try
        implicit none

        integer(4) lyear,idayjl
        real(4) xhour,xlat,xlon,win

        real(4) wintab(1440,720,12)
        integer(4) istart
        real(4) brief,xmon
        real(4) a1,a2,b1,b2,c1,c2
        integer(4) i1,i2,j1,j2,k1,k2

        data istart/1/

        if(istart.eq.1) then
            istart=0
            call openbig_try(3,'O:\ssmi_level3\monthly_wind.dat','old')
            read(3) wintab
            close(3)
        endif
                                                                       
        !     check inputs                                                      

        if(xhour.lt.0  .or. xhour.gt.24.001)   stop 'error1 in fd_monthly_wind, pgm stopped'
        if(idayjl.lt.1    .or. idayjl.gt. 366) stop 'error2 in fd_monthly_wind, pgm stopped'                                                                       
        if(xlon.lt.0..or.xlon.gt.360.)         stop 'error3 in fd_monthly_wind, pgm stopped'                       
        if(abs(xlat).gt.90.)                   stop 'error4 in fd_monthly_wind, pgm stopped'                                   

        !     do time,lat,lon interpolation,  xmon denotes semi-month, not month
 
        if(4*int(lyear/4).eq.lyear) then
            xmon=12*(idayjl-1+xhour/24.)/366.
        else
            if(idayjl.gt.365) stop 'error5 in fd_monthly_wind, pgm stopped'
            xmon=12*(idayjl-1+xhour/24.)/365.
        endif

        brief=xmon-0.5
        i1=1+brief
        i2=i1+1
        a1=i1-brief
        a2=1-a1
        if(i1.eq. 0) i1=12
        if(i2.eq.13) i2= 1 

        brief=4*(xlat+89.875)
        j1=1+brief                                                        
        j2=j1+1                                                           
        b1=j1-brief                                                       
        b2=1-b1
        if(j1.eq.  0) j1=  1
        if(j2.eq.721) j2=720

        brief=4*(xlon-0.125) 
        k1=1+brief                                                       
        k2=k1+1                                                           
        c1=k1-brief                                                       
        c2=1-c1                                                          
        if(k1.eq.   0) k1=1440
        if(k2.eq.1441) k2=   1 
                                                                       
        win=   &          
            a1*b1*(c1*wintab(k1,j1,i1) + c2*wintab(k2,j1,i1))+ a1*b2*(c1*wintab(k1,j2,i1) + c2*wintab(k2,j2,i1))+      &           
            a2*b1*(c1*wintab(k1,j1,i2) + c2*wintab(k2,j1,i2))+ a2*b2*(c1*wintab(k1,j2,i2) + c2*wintab(k2,j2,i2)) 
                      
        return                                                            
    end subroutine fd_monthly_wind                                                               
