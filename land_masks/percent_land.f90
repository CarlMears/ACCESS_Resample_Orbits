module percent_land

    character(100)    :: filename(2)
    character(1)      :: percent_land_tab(4320,2160,5)
    integer(4)        :: iamsrsv = 0
    
    data filename /'X:\LAND\LAND_TABLES_AMSRE.DAT','X:\LAND\LAND_TABLES_AMSRA.DAT'/
    
    ! changed to use "open_big_try" to reduce problems with opening files
    ! cm 6/24/2019
    
    ! changed "end" to "end subroutine fd_percent_land" to make it
    ! more compiler friendly   cm 9/25/2012
contains

    subroutine fd_percent_land(iamsr,xlat,xlon,percent_land)

        use open_file_routines, only: openbig_try
        
        implicit none
        integer(4),intent(in)  :: iamsr
        real(4),intent(in)     :: xlat
        real(4),intent(in)     :: xlon
        real(4),intent(out)    :: percent_land(5)
        
        integer(4)             :: ilat
        integer(4)             :: ilon

        if(iamsr .ne. iamsrsv) then
        iamsrsv = iamsr
        call openbig_try(3,filename(iamsr),'old')
        read(3) percent_land_tab
        close(3)
        endif

        ilat=1+nint(12.*(xlat+89.95833))
        ilon=1+nint(12.*(xlon- 0.04167))
        if(ilat.lt.1) ilat=1
        if(ilon.lt.1) ilon=1
        if(ilat.gt.2160) ilat=2160
        if(ilon.gt.4320) ilon=4320

        percent_land=0.4*ichar(percent_land_tab(ilon,ilat,:))
        
        return
    end subroutine fd_percent_land
end module percent_land