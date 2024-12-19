	subroutine read_binary_head(iorbit,filename_rss_l1a, numscan,scan_time)
	use amsr_rss_l1a_input_module						 			 		 	   
      implicit none

 	character(100) filename_rss_l1a   
	real(8) scan_time
      integer(4) iorbit,numscan,isize,ierr,nrec

	call fastunzip(filename_rss_l1a, isize,apar,ierr)
	if(ierr.ne.0) stop 'error in fastunzip when reading amsr l0 file'

	nrec=int(isize/reclen)
	if(nrec.le.1 .or. nrec.gt.maxscan+1 .or. reclen*nrec.ne.isize) stop 'pgm stopped, error1 read_binary_header'
	numscan=nrec-1  !subtract header scan

	abuf1=apar(0)
	if(l1arss_hd.iorbit.ne.iorbit .or. l1arss_hd.numscan.ne.numscan) stop 'wrong orbit or numscan, pgm stopped'
	scan_time=l1arss_hd.scan_time

	return
	end
