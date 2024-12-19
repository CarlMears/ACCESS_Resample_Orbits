module satellite_constants

    
     ! prototype use statement for this module
	 !use satellite_constants,only: nBand,nChannel,nFOV, nScan,nFOVR, nScanR,      &  
	 !                              nNorthPoleX,nNorthPoleY,nSouthPoleX,nSouthPoleY,  & 
	 !                              ndays_ice_history,                                &
	 !                              satellite_name,                                   &
	 !                              BandNames,BandsForIce,channel_num_offset,         &
	 !                              history_root,ancillary_root
	 
	implicit none

    !everything in this module is public
    integer(4),parameter  ::      nBand = 8
    integer(4),parameter  ::      nBand_tot = 8
    integer(4),parameter  ::      nPol = 2
    integer(4),parameter  ::      nPol_tot = 2
    integer(4),parameter  ::      nChannel = 16
    integer(4),parameter  ::      nFOV  = 243
    integer(4),parameter  ::      FOV_center = 122
    integer(4),parameter  ::      nScan = 1300
    integer(4),parameter  ::      nFOVR  = 243
    integer(4),parameter  ::      nScanR = 1300
    integer(4),parameter  ::      nNorthPoleX = 304
    integer(4),parameter  ::      nNorthPoleY = 448
    integer(4),parameter  ::      nSouthPoleX = 316
    integer(4),parameter  ::      nSouthPoleY = 332
    integer(4),parameter  ::      ndays_ice_history = 30
    integer(4),parameter  ::      maxrev = 100000 	
    
    integer(4),parameter  ::      ksat = 34 
    
    character(10) :: satellite_name = 'AMSR2'
    
    character(4),dimension(nband) :: BandNames = (/'6.9','7.3','11','19','24','37','89a','89b'/)
    integer(4),dimension(5)       :: BandsForIce = (/3,4,5,6,7/)
    integer(4)                    :: channel_num_offset = 0
    character(150)                :: history_root   = 'L:\sea_ice\amsr2_ice_history\'
    character(150)                :: ancillary_root = 'L:\sea_ice\amsr2_ancillary_data\'
    character(150)                :: rtm_root       = 'L:\sea_ice_test\'
 
 end module satellite_constants
    
    
    