module date_utilities

    !use logger,       only: log_debug, log_info, log_error, log_crit, log_msg_string
!  
    implicit none
    
    private
    
    public find_month_day,fd_time_2000,fd_date_2000,days_in_month,doy_from_year_month_day
    
    integer(4),dimension(12,2),parameter :: start_day_array = reshape((/1,32,60,91,121,152,182,213,244,274,305,335, &           
                                                                1,32,61,92,122,153,183,214,245,275,306,336/), (/12,2/))
                                           
    integer(4),dimension(12,2),parameter :: days_in_month_arr =reshape((/31,28,31,30,31,30,31,31,30,31,30,31,    &
                                                                 31,29,31,30,31,30,31,31,30,31,30,31/), (/12,2/))


contains
    
    subroutine find_month_day(lyear,idayjl, imon,idaymo)
    
        implicit none
        
        integer(4),intent(in)      :: lyear
        integer(4),intent(in)      :: idayjl
        integer(4),intent(out)      :: imon
        integer(4),intent(out)      :: idaymo
        
        integer(4)                 :: ileap
        integer(4)                 :: jmon
        
        

        ileap=1
        if(lyear.eq.4*int(lyear/4)) ileap=2
        do 10 jmon=2,12
            if(start_day_array(jmon,ileap).gt.idayjl) then
            imon=jmon-1
            go to 20
        endif
       10 continue
          imon=12
       20 continue

       idaymo=1+idayjl-start_day_array(imon,ileap)
       return
    end subroutine find_month_day
    
    subroutine fd_time_2000(lyear,imon,idaymo,secdy,time2000,secyr)
      
        implicit none
        
        integer(4),intent(in)   :: lyear
        integer(4),intent(in)   :: imon
        integer(4),intent(in)   :: idaymo
        real(8),intent(in)      :: secdy
        real(8),intent(out)     :: time2000
        real(8),intent(out)     :: secyr
        
        integer(4)  :: ileap
        real(8)     :: time1987
        
        ileap = 1
        if(lyear .eq. (4*int(lyear/4))) then
            ileap = 2    
        endif    
        
        secyr = secdy + 86400.0D0*(start_day_array(imon,ileap) + idaymo - 2)
        time1987 = secyr + 31536000.D0*(lyear-1987)+86400.D0*int((lyear-1985)/4)
        time2000 = time1987 - 410227200.D0
     end subroutine fd_time_2000
        
     subroutine fd_date_2000(time_2000, secyr,lyear,idayjl,imon,idaymo,secdy)
     
        implicit none
        
        ! input
        !     time_2000 = sec from begin of 2000
        
        ! outputs
        !     isecyr    = sec from begin of lyear
        !     lyear     = year, full for digits, i.e., 1987
        !     idayjl    = julian dat of year
        !     imon      = month (1 to 12)                   
        !     idaymo    = day of month (1 to 31)
        !     secdy     = seconds of day (0 to 86399.99999)
        
        real(8),intent(in)      ::  time_2000
      
        real(8),intent(out)     ::  secyr
        integer(4),intent(out)  ::  lyear
        integer(4),intent(out)  ::  idayjl
        integer(4),intent(out)  ::  imon
        integer(4),intent(out)  ::  idaymo
        real(8),intent(out)     ::  secdy 
        
        ! local variables
        real(8) xtime 
        integer(4) idaytot,idaybg,ileap,jmon
                                                    
        xtime=time_2000 + 410227200.d0 !convert to time87 
    
        if(xtime.lt.0 .or. xtime.gt.3534451200.d0) then
            !call log_error('time oob in fd_date_2000')  !3534451200 is jan 1 2099   
            print *,'time oob in fd_date_2000'  !)  !3534451200 is jan 1 2099   
            return
        endif    
  
        idaytot=1 + int(xtime/86400.d0)
        lyear=1987 + int((idaytot-1)/365)                        !kyear may be 1 year too big
        idaybg= 1 + 365*(lyear-1987) + int((lyear-1985)/4)       !begin day of kyear
        if(idaytot.lt.idaybg) lyear=lyear-1
  
        secyr=xtime-31536000.d0*(lyear-1987)-86400.d0*int((lyear-1985)/4)
                                                  
        ileap=1                                                            
        if(lyear.eq.4*int(lyear/4)) ileap=2   
                                    
        idayjl=1+int(secyr/86400.d0)                                           
        do jmon=2,12                                                  
            imon=jmon-1                                                       
            if(start_day_array(jmon,ileap).gt.idayjl) goto 210  
        enddo                        
        imon=12                                                         
  210   continue                                                          
        idaymo=idayjl-start_day_array(imon,ileap)+1                                  
        secdy=secyr-(idayjl-1)*86400.d0                                         
        return                                                            
    end subroutine fd_date_2000                                                            
  
    function days_in_month(year,month,start_day)

        integer(4),intent(in)  :: year
        integer(4),intent(in)  :: month
        integer(4),optional    :: start_day
        
        integer(4)       :: ileap
        integer(4)       :: days_in_month
        
        if ((month .ge. 1) .and. (month .le. 12)) then   
            ileap = 1
            if (modulo(year,4) .eq. 0) ileap = 2
            if (modulo(year,100) .eq. 0) ileap = 1
            if (modulo(year,400) .eq. 0) ileap = 2
            if (present(start_day)) then
                start_day = start_day_array(month,ileap)
            endif    
            days_in_month = days_in_month_arr(month,ileap)
            return
         else
            !call log_error('Error, month out of range')
            print *,'Error, month out of range'
            if (present(start_day)) then
                start_day = 0
            endif
            days_in_month = 0
            return
        endif
    end function days_in_month

    function doy_from_year_month_day(year,month,day)

        integer(4),intent(in)  :: year
        integer(4),intent(in)  :: month
        integer(4),intent(in)  :: day
        
        integer(4)             :: doy_from_year_month_day
        integer(4)             :: ileap
        !integer(4)             :: days_in_month
        integer(4)             :: start_day
        
        if ((month .ge. 1) .and. (month .le. 12)) then   
            ileap = 1
            if (modulo(year,4) .eq. 0) ileap = 2
            if (modulo(year,100) .eq. 0) ileap = 1
            if (modulo(year,400) .eq. 0) ileap = 2
        
            start_day = start_day_array(month,ileap) 
            doy_from_year_month_day = start_day + day -1
            return
         else
            !call log_error('Error, month out of range')
            print *,'Error, month out of range'
            doy_from_year_month_day = 0
            return
        endif
    end function doy_from_year_month_day
end module date_utilities