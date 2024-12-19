module findice

    character(100),private                       :: ice_filename ='x:\land\ice4.dat'
    character(1),dimension(45,180,12),private    :: amask_ice           !monthly climatology sea ice map
    logical(4),private                           :: ice_mask_loaded = .false.

contains

    subroutine fdice(imon,xlat,xlon, ice)                                                      
        implicit none

        integer(4) imon,ice
        real(4) xlat,xlon

        integer(4) ifac(8),ilat,ilon,jlon,ibit,iii

        data ifac/128,64,32,16,8,4,2,1/
        
        if (.not. ice_mask_loaded) then
            open(3,file=ice_filename,status='old',form='binary')
            read(3) amask_ice
            close(3)
            ice_mask_loaded = .true.
        endif

        ilat=1+nint(xlat+89.5)
        ilon=1+nint(xlon- 0.5)
        if(ilat.lt.1) ilat=1
        if(ilon.lt.1) ilon=1
        if(ilat.gt.180) ilat=180
        if(ilon.gt.360) ilon=360

        jlon=1 + int((ilon-1)/8)
        ibit=ilon-8*(jlon-1)
      
        iii=int(ichar(amask_ice(jlon,ilat,imon))/ifac(ibit))
        ice=1
        if(2*int(iii/2).eq.iii) ice=0
    return
    end subroutine fdice
end module findice