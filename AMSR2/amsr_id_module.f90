    module sat_id_module
        implicit none
        
        !integer(4), parameter :: iamsr=1  !leave in for now until i remove this everywhere
        character(5), parameter :: satname = 'AMSR2'
        integer(4), parameter :: maxcel=243
        integer(4), parameter :: maxcel_89=2*maxcel      
        integer(4), parameter :: low_cal_obs= 16
        integer(4), parameter :: high_cal_obs=2*low_cal_obs
        integer(4), parameter :: maxscan=4400
        integer(4), parameter :: maxscan_89=2*maxscan
        integer(4), parameter :: maxrev=100000
        integer(4), parameter :: max_jaxa_orbits=200000 
        
        integer(4), parameter :: nfreq=8
        integer(4), parameter :: nch=2*nfreq
        
        integer(4), parameter :: l1a_lo_chan=12
        integer(4), parameter :: l1a_hi_chan=4

        integer(4), parameter :: spc_length=34  !amsre=20
        integer(4), parameter :: sps_length=46  !amsre=32
        integer(4), parameter :: obs_supplement=124  !amsr=27
        
        integer(4), parameter :: max_scans_in_granule = 1300
        
        !maps freq onto band number for resampling
        integer(4), parameter ::band_number(nfreq) = (/1,1,2,3,4,5,6,6/)
        character(3), parameter :: freq_name(nfreq) = (/'6  ','7  ','11 ','19 ','24 ','37 ','89a','89b'/)
        character(4), parameter :: channel_name(nch) = (/'6V  ','6H  ','7V  ','7H  ','11V ','11H ','19V ','19H ', &
                                                         '24V ','24H ','37V ','37H ','89aV','89aH','89bV','89bH'/)
    end module sat_id_module

    
