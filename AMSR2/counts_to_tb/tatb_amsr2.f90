subroutine tatb_amsr2(icase,ifreq, tax,tb)
	implicit none
	
    integer(4), parameter :: nfreq=8
    real(4), parameter :: tcld_plk(nfreq)=(/2.73, 2.73, 2.74, 2.75, 2.77, 2.82, 3.24, 3.24/) !tcold for apc see 'o:\emiss_model\planck.txt'
	
 
	integer(4) icase,ifreq,jfreq,istart,ich,jch
	real(4) tax(2),tb(2)
      real(4) delta(2,nfreq),chi(2,nfreq)
      real(4) qvv(nfreq),qhv(nfreq),qov(nfreq),qhh(nfreq),qvh(nfreq),qoh(nfreq)
      real(4) avv(nfreq),ahv(nfreq),aov(nfreq),ahh(nfreq),avh(nfreq),aoh(nfreq) 
	real(4) xfac,tv,th
 
!     data initialization
 
      data istart/1/
 
!     begin execution
 
      if(istart.eq.1) then
        istart=0
	    open(3,file='j:\amsr2\tables\apc_coefs_v8_f34.txt',status='old')
	    read(3,3001) delta
	    read(3,3001) chi
 3001   format(16f10.6)
        close(3)

          do jfreq=1,nfreq
     
              qvv(jfreq)=(1-delta(1,jfreq))/(1.+chi(1,jfreq))
              qhv(jfreq)=chi(1,jfreq)*qvv(jfreq)
              qov(jfreq)=delta(1,jfreq)*tcld_plk(jfreq)
              qhh(jfreq)=(1-delta(2,jfreq))/(1.+chi(2,jfreq))
              qvh(jfreq)=chi(2,jfreq)*qhh(jfreq)
              qoh(jfreq)=delta(2,jfreq)*tcld_plk(jfreq)
     
	          xfac=qvv(jfreq)*qhh(jfreq)-qhv(jfreq)*qvh(jfreq)
     
	          avv(jfreq)= qhh(jfreq)/xfac
	          ahv(jfreq)=-qhv(jfreq)/xfac
	          aov(jfreq)=(qhv(jfreq)*qoh(jfreq)-qhh(jfreq)*qov(jfreq))/xfac
     
	          ahh(jfreq)= qvv(jfreq)/xfac
	          avh(jfreq)=-qvh(jfreq)/xfac
	          aoh(jfreq)=(qvh(jfreq)*qov(jfreq)-qvv(jfreq)*qoh(jfreq))/xfac
     
	       enddo
      endif
      
 
    ich=1
	jch=2

	if(icase.eq.1)  then !tb to ta
 
	    if(tb(ich).lt.55 .or. tb(ich).gt.330 .or. tb(jch).lt.55 .or. tb(jch).gt.330) then
	        tax(ich)=tb(ich)
	        tax(jch)=tb(jch)
	    else
	        tv=tb(ich)
	        th=tb(jch)
            tax(ich)=qvv(ifreq)*tv + qhv(ifreq)*th + qov(ifreq)
            tax(jch)=qhh(ifreq)*th + qvh(ifreq)*tv + qoh(ifreq)
	    endif
	else  !ta to tb
 
	    if(tax(ich).lt.55 .or. tax(ich).gt.330 .or. tax(jch).lt.55 .or. tax(jch).gt.330) then
	        tb(ich)=tax(ich)
	        tb(jch)=tax(jch)
	    else
	        tv=tax(ich)
	        th=tax(jch)
            tb(ich)=avv(ifreq)*tv + ahv(ifreq)*th + aov(ifreq)
            tb(jch)=ahh(ifreq)*th + avh(ifreq)*tv + aoh(ifreq)
	    endif
	endif
 
	return
	end