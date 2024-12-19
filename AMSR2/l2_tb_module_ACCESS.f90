! This is a L2 module for the access project
! This one holds the unresampled brightness temperatures
! I am also storing the incidence angle and look angle
! since they should be part of the ACCESS data

! 6.9 and 7.2 GHz data are included

! so there are 8 channels
! 6.9,7.2,11,19,24,37,89a,89b
! and with 2 polarizations
! 16 brightness temperatures

module l2_tb_module_ACCESS

    use satellite_constants, only: nChannel,nBand,nPol,nFOV,nScan
    implicit none
    private
    public tb_scn,tb_grn,    &                  !types needed by other routines
           init_l2_tb_scan,  &                  !initialize a TB scan
           read_l2_tb_scan,read_l2_tb_file, &   !read routine
           write_l2_tb_scan, &                  !write routine
           init_geoloc_granule,free_geoloc_granule
        
    type tb_scn
        integer(4)      :: numcels
        integer(4)      :: scan_index
        integer(4)      :: scan_qual
        real(8)         :: scan_time
        real(8)         :: orbit
        real(4)         :: cellon(nFOV)
        real(4)         :: cellat(nFOV)
        real(4)         :: tb_measured(nBand,nPol,nFOV)
        integer(4)      :: cel_qual(nFOV)
        real(4)         :: perc_land(nFOV)
        integer(4)      :: i_ice(nFOV)
        real(4)         :: incidence_angle(nFOV)
        real(4)         :: look_angle(nFOV)
    endtype tb_scn
    
    type geoloc_scn
        integer(4)      :: nFOV
        real(8)         :: scan_index
        integer(4)      :: scan_qual

        real(8)         :: scan_time
        real(8)         :: orbit
        real(4), allocatable        :: cellon(:)
        real(4), allocatable        :: cellat(:)
        real(4), allocatable        :: perc_land(:)
        integer(4),allocatable      :: i_ice(:)
        real(4),allocatable         :: incidence_angle(:)
        real(4),allocatable         :: look_angle(:)

    endtype geoloc_scn
    
    type tb_grn
        integer(4),dimension(nScan)    :: scan_valid     ! 1 if valid, 0 is not valid.
        integer(4),dimension(nScan)    :: scan_qual        
        type(tb_scn),dimension(nScan)  :: scans
    endtype
    
    type geoloc_grn
        integer(4)                     :: num_scans
        integer(4),allocatable         :: scan_valid(:)
        integer(4),allocatable         :: scan_qual(:)
        type(geoloc_scn),allocatable   :: scans(:)
    endtype

    
    !integer(4),parameter   :: l2_emiss_module_scan_size_bytes = 19456   ! must recalculate

contains

    subroutine init_geoloc_scan(geoloc_scan,numfov)

        type(geoloc_scn),intent(inout) :: geoloc_scan
        integer(4), intent(in)         :: numfov

        geoloc_scan%nFOV = numfov
        geoloc_scan%scan_index = -999999999.9
        geoloc_scan%scan_qual = 0

        geoloc_scan%scan_time = -999999999.9
        geoloc_scan%orbit = -999999999.9

        allocate(geoloc_scan%cellon(numfov))
        allocate(geoloc_scan%cellat(numfov))
        allocate(geoloc_scan%perc_land(numfov))
        allocate(geoloc_scan%i_ice(numfov))
        allocate(geoloc_scan%incidence_angle(numfov))
        allocate(geoloc_scan%look_angle(numfov))

        geoloc_scan%cellat       = -999.9
        geoloc_scan%cellon       = -999.9
        geoloc_scan%perc_land    = -999.9
        geoloc_scan%i_ice        = -999.9
        geoloc_scan%incidence_angle = -999.0
        geoloc_scan%look_angle   = -999.0
    end subroutine init_geoloc_scan

    subroutine free_geoloc_scan(geoloc_scan)

        type(geoloc_scn),intent(inout) :: geoloc_scan
        

        

        deallocate(geoloc_scan%cellon)
        deallocate(geoloc_scan%cellat)
        deallocate(geoloc_scan%perc_land)
        deallocate(geoloc_scan%i_ice)
        deallocate(geoloc_scan%incidence_angle)
        deallocate(geoloc_scan%look_angle)

    end subroutine free_geoloc_scan

    subroutine init_geoloc_granule(gran,numscan,numfov)

        type(geoloc_grn),intent(inout)      :: gran 
        integer(4),intent(in)               :: numscan
        integer(4),intent(in)               :: numfov

        integer(4)                          :: i

        gran%num_scans = numscan
        allocate(gran%scan_valid(numscan))
        allocate(gran%scan_qual(numscan))
        allocate(gran%scans(numscan))

        gran%scan_valid = 0
        gran%scan_qual = 0
        do i = 1,numscan
            call init_geoloc_scan(gran%scans(i),numfov)
        enddo

    end subroutine init_geoloc_granule

    subroutine free_geoloc_granule(gran)
        type(geoloc_grn),intent(inout)      :: gran  
        integer(4)                          :: i

        deallocate(gran%scan_valid)
        deallocate(gran%scan_qual)
        deallocate(gran%scans)
        
        do i = 1,gran%num_scans
            call free_geoloc_scan(gran%scans(i))
        enddo
    end subroutine free_geoloc_granule

    subroutine init_l2_tb_scan(tb_scan)
    
        type(tb_scn),intent(OUT)  :: tb_scan
     
        tb_scan%numcels      = nFOV
        tb_scan%scan_index   = -999
        tb_scan%scan_qual    = 0
        tb_scan%scan_time    = -9999999999.9
        tb_scan%orbit        = -9999999999.
        tb_scan%cellat       = -999.9
        tb_scan%cellon       = -999.9
        tb_scan%incidence_angle = -999.0
        tb_scan%look_angle   = -999.0
        tb_scan%tb_measured  = -999.0
        tb_scan%cel_qual     = 0
                        
        tb_scan%i_ice = 0
        tb_scan%perc_land = -999.0
    
     end subroutine init_l2_tb_scan
     
     subroutine read_l2_tb_scan(ilu_l2b_tb,tb_scan,read_error)
     
        ! this reads the files in L:\sea_ice\amsr2_tb_orbits\ that were written on or just after June 22
     
        integer(4),intent(IN)     :: ilu_l2b_tb
        type(tb_scn),intent(OUT)  :: tb_scan
        integer(4),intent(OUT)    :: read_error

        read(ilu_l2b_tb,err = 200,end=300) tb_scan%scan_index,tb_scan%scan_time,         &
                                           tb_scan%orbit,                                &
                                           tb_scan%cellat,tb_scan%cellon,                &
                                           tb_scan%tb_measured,tb_scan%i_ice,            &
                                           tb_scan%perc_land,tb_scan%incidence_angle,    &
                                           tb_scan%look_angle
        read_error = 0
        return
        
 200    read_error = -1
        return
        
 300    read_error = 1
        return
     end subroutine read_l2_tb_scan
     
     subroutine read_l2_tb_file(ilu_l2b_tb,tb_granule,num_scans_read,error)
        
        integer(4),intent(IN)     :: ilu_l2b_tb
        type(tb_grn),intent(OUT)  :: tb_granule
        integer(4),intent(OUT)    :: num_scans_read
        integer(4),intent(OUT)    :: error
        
        integer(4)    :: iscan
        integer(4)    :: read_error
        
        tb_granule%scan_valid = 0
        tb_granule%scan_qual  = 0
        num_scans_read        = 0
        error = 1  !error is 1 if file is too big
        
        do iscan = 1,nScan
            call init_l2_tb_scan(tb_granule%scans(iscan))
        enddo
        
        do iscan = 1,nScan
            call read_l2_tb_scan(ilu_l2b_tb,tb_granule%scans(iscan),read_error)
            if (read_error .eq. 0) then
                num_scans_read = num_scans_read+1
                tb_granule%scan_valid(iscan) = 1
                ! might add code here to summarize cel_qual into  scan_qual for now do nothing
            else
                call init_l2_tb_scan(tb_granule%scans(iscan)) ! make sure it stays zeroed
                error = -1
                if (read_error .eq. 1) error = 0  !set error to zero if it is the end of the file
                exit
            endif
        enddo
     end subroutine read_l2_tb_file 
          
     subroutine write_l2_tb_scan(ilu_l2b_tb,tb_scan,write_error)
     
        integer,intent(IN)        :: ilu_l2b_tb
        type(tb_scn),intent(IN)   :: tb_scan
        integer(4),intent(OUT)    :: write_error

        write(ilu_l2b_tb,err = 200) tb_scan%scan_index,tb_scan%scan_time,         &
                                    tb_scan%orbit,                                &
                                    tb_scan%cellat,tb_scan%cellon,                &
                                    tb_scan%tb_measured,tb_scan%i_ice,            &
                                    tb_scan%perc_land,tb_scan%incidence_angle,    &
                                    tb_scan%look_angle
        write_error = 0
        return
        
 200    write_error = -1
        return
               
     end subroutine write_l2_tb_scan
 
end module l2_tb_module_ACCESS
