!     aug 23 2013 changed 9/6/2017.  The following statement was inserted
!        if(itai_1993_leap(ileap).eq.-2147483648)         exit 
!     this makes the routine more robust.  it will work if itai_1993_leap is
!     read in from x:\land\ut1.dat or from x:\land\tai_leap_table.dat (latter is now obsolete)
!     see 'O:\amsr2\L2_processing\memo7.txt'
!     i also now do a check to make sure the table has been read into time_tai_1993

        subroutine fd_leap_sec(time_tai_1993, nleap)
            use l2_module                                                      
            implicit none
    
            real(8) time_tai_1993
            integer(4) nleap
            
            integer(4) ileap
            
            if(itai_1993_leap(1).ne.15638401) stop 'leap table not read in, pgm stopped'

        !     negative itai_1993_leap(ileap) denote negative leap second 
            nleap=0
            do ileap=1,100
                if(itai_1993_leap(ileap).eq.-2147483648)         exit !missing data
                if(iabs(itai_1993_leap(ileap)).gt.time_tai_1993) exit 
                if(itai_1993_leap(ileap).ge.0) then
                    nleap=nleap + 1
                else
                    nleap=nleap - 1
                endif
            enddo
            
            return
        end subroutine fd_leap_sec
