$freeform
module landmask

    use NAN_support 
    use maps

    implicit none


    ! declarations for the array that contains the land mask info

    integer,parameter                ::  NUM_LONG_LANDMASK_RAW        = 540
    integer,parameter                ::  NUM_LAT_LANDMASK_RAW        = 2160
    integer,parameter                ::  NUM_LONG_LANDMASK_INT        = 4320
    integer,parameter                ::  NUM_LAT_LANDMASK_INT        = 2160
    integer,parameter                ::  NUM_LONS_PER_BYTE            = 8
    real(4),parameter                ::    NUM_PER_DEGREE_LANDMASK        = 12.0

    integer,parameter    ::    NUM_SURFACES    = 2

    

    byte,dimension(NUM_LONG_LANDMASK_RAW,NUM_LAT_LANDMASK_RAW)        ::    land_mask_byte
    integer(4),dimension(NUM_LONG_LANDMASK_INT,NUM_LAT_LANDMASK_INT)::  land_mask
    integer(4)                                                        ::  land_mask_read = 0
    character*80,parameter                                            ::     LAND_MASK_FILE = 'x:\land\mimask5.dat'
    integer,parameter                                                ::    LU_LAND_MASK = 15


    integer,parameter                ::  SURFACE_UNKNOWN = 0
    integer,parameter                ::    SEA                = 1
    integer,parameter                ::    LAND            = 2
    integer,parameter                ::  BOTH_SURFS        = 3

    ! declarations for land fraction
    
    character(len = 100),parameter                    :: LAND_FRACTION_FILE = '\\sounder\C\Module_Data\Land_Fraction\land_fraction_144x72_map.dat'
    integer(4),parameter                            :: LU_LAND_FRACTION    = 15
    type(map_str)                                    :: land_frac_map
    integer(4)                                        :: land_fraction_read=0

contains

    subroutine read_land_mask(land_mask,error)

        integer(4),dimension(NUM_LONG_LANDMASK_INT,NUM_LAT_LANDMASK_INT),intent(OUT)::  land_mask
        integer(4), intent(OUT)                                                        ::  error    

        integer(4)        :: ilon
        integer(4)        :: ilat
        integer(4)        :: j
        
        ! begin execution

        error = 0

        open(LU_LAND_MASK,file=LAND_MASK_FILE,form = 'BINARY',iostat = error)

        if (error .eq. 0) then

            read(unit = LU_LAND_MASK,iostat = error) land_mask_byte

            if (error .eq. 0) then

                do ilon = 1, NUM_LONG_LANDMASK_RAW
                    do ilat = 1, NUM_LAT_LANDMASK_RAW
                        do j = 0,NUM_LONS_PER_BYTE-1 
                            land_mask(((ilon-1)*NUM_LONS_PER_BYTE)+j+1,ilat) =  ibits(land_mask_byte(ilon,ilat),j,1)
                        enddo
                    enddo
                enddo
                print *,'done reading landmask'
                close(LU_LAND_MASK)
                return
            else
                print *,' Error reading land mask file!'
                print *,' Advise Stopping.......       '
                return
            endif
        else
            print *,' Error opening land mask file!'
            print *,' Advise Stopping.......       '
        endif
    end subroutine read_land_mask


    integer(4) function find_land_sea(lat,long)

        use NAN_support

        implicit none

        real(4),intent(IN)    :: lat  ! latitude of observation
        real(4),intent(IN)    :: long ! longitude of observation

        integer(4)                :: ilat
        integer(4)                :: ilong

        integer(4)                :: error

        if ((abs(lat) .gt. 90.0) .or. (long .lt. 0.0) .or. (long .gt. 360.0)) then
            find_land_sea = -1
            return
        endif

        ! read in land mask if it hasn't been read yet

        if (land_mask_read .eq. 0) then
            write(6,*)'Reading Land mask ',LAND_MASK_FILE
            call read_land_mask(land_mask,error)
            land_mask_read = 1
        endif

        ilat = int(NUM_PER_DEGREE_LANDMASK * (lat + 90.0)) + 1
        ilong = int(NUM_PER_DEGREE_LANDMASK * long) + 1
        if ((ilong .ge. 1) .and. (ilong .le. NUM_LONG_LANDMASK_INT) .and. (ilat .ge. 1) .and. (ilat .le. NUM_LAT_LANDMASK_INT)) then
            if (land_mask(ilong,ilat) .ne. 0) then
                find_land_sea = LAND
            else
                find_land_sea = SEA
            endif
        else
            find_land_sea = SURFACE_UNKNOWN
        endif

    end function find_land_sea

    real(4) function find_land_fraction(long,lat,no_interp)

        use NAN_support
        use maps

        implicit none

        real(4),intent(IN)                :: long  !longitude
        real(4),intent(IN)                :: lat   !latitude
        logical,intent(IN),optional        :: no_interp

        integer(4)                :: error
        logical                    :: interp_flag
        real(4)                    :: lf

        ! read in the land fraction data if it hasn't been read yet

        if (land_fraction_read .eq. 0) then
            call read_land_fraction_map(land_frac_map,error)
            if (error .eq. 0) land_fraction_read = 1
        endif

        if (.not. present(no_interp)) then 
            interp_flag = .true.
        else
            interp_flag = .not. no_interp
        endif

        if (interp_flag) then
            call get_map_value(land_frac_map,long,lat,lf,error)
        else
            call get_map_value_no_interp(land_frac_map,long,lat,lf,error)
        endif

        find_land_fraction = lf

        return

    end function find_land_fraction

    subroutine read_land_fraction_map(lf_map,error)

        type(map_str),intent(OUT)    :: lf_map
        integer(4),intent(INOUT)        :: error


        open(unit = LU_LAND_FRACTION,file = LAND_FRACTION_FILE,form = 'binary')
        call read_map(LU_LAND_FRACTION,lf_map,error)

        return

    end subroutine read_land_fraction_map
     








    
end module landmask