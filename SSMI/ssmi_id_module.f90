    module sat_id_module
        implicit none
        
        !integer(4), parameter :: iamsr=1  !leave in for now until i remove this everywhere
        character(5), parameter :: satname = 'SSMI'
        integer(4), parameter :: maxcel=64
        integer(4), parameter :: maxcel_85=2*maxcel      
        integer(4), parameter :: maxscan=3600
        integer(4), parameter :: maxrev=1000000

        integer(4), parameter :: nfreq=4
        integer(4), parameter :: nch=2*nfreq
        integer(4), parameter :: nChannel=2*nfreq
        integer(4), parameter :: maxchn=8
        
    
        !maps freq onto band number for resampling
        integer(4), parameter ::band_number(nfreq) = (/1,2,3,4/)
        character(3), parameter :: freq_name(nfreq) = (/'19 ','22 ','37 ','85'/)
        character(4), parameter :: channel_name(nch) = (/'19V ','19H ', &
                                                         '22V ','22H ', &
                                                         '37V ','37H ', &
                                                         '85V','85H'/)
    end module sat_id_module