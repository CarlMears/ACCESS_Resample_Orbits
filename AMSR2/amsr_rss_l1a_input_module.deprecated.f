c     3/8/2013 version changed 5/3/2013 arrays loEarthCounts_rss and hiEarthCounts_rss are added
c     reclen had to be increase for 30096 to 39816

	module amsr_rss_l1a_input_module
      use sat_id_module
	implicit none
	
	integer(4), parameter :: reclen =39816 !you can determine this by writing out L1arss and looking at it size
	integer(4), parameter :: npad=reclen-16

	character(reclen) apar(0:maxscan)
 	character(reclen) abuf1
	character(reclen) abuf2
     
c     ==========================================================================
c     ======================= header structure =================================
c     ==========================================================================
      type level1a_rss_header
 	sequence
	real(8) scan_time
	integer(4) iorbit
	integer(4) numscan
	character(npad) apad
 	end type level1a_rss_header
 	type (level1a_rss_header) :: L1arss_hd

c     ==========================================================================
c     ======================= rss l1a data structure ===========================
c     ==========================================================================

	type level1a_rss
	sequence
      
	real(8) scan_time
	real(8) rev
	
	real(4) omega
	real(4) navigation(6)
	real(4) attitude(3)

	integer(4) iscan
	integer(4) spc_temp_count(spc_length)
	integer(4) sps_temp_count(sps_length)
	integer(4) rx_offset_gain(2,nch)
	integer(4) ipad  !this makes the structure a multiple of 8

	real(4) lat89a(maxcel_89)
	real(4) lat89b(maxcel_89)
	real(4) lon89a(maxcel_89)
	real(4) lon89b(maxcel_89)
	real(4) thot(nch)
	
	integer(2) loEarthCounts_rss(maxcel,   l1a_lo_chan)
	integer(2) hiEarthCounts_rss(maxcel_89,l1a_hi_chan)

	integer(2) loEarthCounts(maxcel,   l1a_lo_chan)
	integer(2) hiEarthCounts(maxcel_89,l1a_hi_chan)
	integer(2) cld_lo_counts(l1a_lo_chan, low_cal_obs)
	integer(2) hot_lo_counts(l1a_lo_chan, low_cal_obs)
	integer(2) cld_hi_counts(l1a_hi_chan,high_cal_obs)
	integer(2) hot_hi_counts(l1a_hi_chan,high_cal_obs) 
 	integer(2)          loTb(maxcel,   l1a_lo_chan)  
	integer(2)          hiTb(maxcel_89,l1a_hi_chan)  
	
 	integer(2) interp_flag(32,nch)
	
 	end type level1a_rss
 
 	type (level1a_rss) :: L1arss
 	
	equivalence(L1arss_hd ,abuf1)
	equivalence(L1arss,    abuf2)

	end module amsr_rss_l1a_input_module
	
