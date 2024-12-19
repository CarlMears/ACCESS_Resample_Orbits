 !        aug 23 2013 version changed 11/19/2013.  bug was fixed. see cbug

 !        tb is 6.9 to 37 ghz without 7.3 ghz
      subroutine ckocntb(tb,ifreq1,ifreq2,  ioob)
    use l2_module
    implicit none

 !        amsr info
 !        tht for amsr (55 deg) is 1.6 deg higher than ssmi
 !        1.6*3 = 4.8 k more v-h seperation for amsr under clear skies
 !        thus diftb must be multipled by 1 -4.8/85 = 0.944, where 85 is indicative of max pol dif
      real(4), parameter ::  tscale=0.944

    integer(4) ifreq1,ifreq2,ioob,ifreq,ipol,ich
    real(4) tb(*)

    ioob=0
      
      ich=0
    do ifreq=ifreq1,ifreq2
!bug    if(ifreq.eq.2) cycle  !do not check 7.3 GHz
    do ipol=1,2
    ich=ich+1
    if(ifreq.eq.2) cycle  !do not check 7.3 GHz
      if(tb(ich).lt.tbmin(ipol,ifreq) .or. tb(ich).gt.280.) ioob=1
    enddo
    enddo

 !        rss processing does not do following check, but i may activate it at some point
 !       if(ioob.eq.0) call chkpol(tscale,tb(7),tb(8),tb(11),tb(12), ioob)

    return
    end




    subroutine chkpol(tscale,tb19v,tb19h,tb37v,tb37h, ioob)
    implicit none

      integer(4) minpol19(31),minpol37(31),maxpol19(31),maxpol37(31),ibin,ioob
    real(4) tscale,tb19v,tb19h,tb37v,tb37h,avgtb,diftb

      data minpol19/71,66,61,57,53,50,47,44,41,39,37,34,31,28,25,23,21,19,17,15,13,11, 9, 7, 5, 3, 1, 0,-1,-1,-1/
      data maxpol19/89,89,89,88,87,85,83,81,79,76,73,70,67,64,61,57,53,50,47,43,39,36,33,30,27,24,21,18,15,12, 9/
      data minpol37/61,61,61,61,61,61,61,55,49,44,39,34,29,22,15, 7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1/
      data maxpol37/85,85,85,85,85,85,85,85,85,84,83,80,77,73,69,64,59,55,51,47,43,39,35,32,29,25,21,18,15,11, 7/

      avgtb=  0.5*(tb19v+tb19h)
      diftb=tscale*(tb19v-tb19h) 
      ibin=nint(0.2*(avgtb-125))
      if(diftb.lt.minpol19(ibin) .or. diftb.gt.maxpol19(ibin)) ioob=1

      avgtb=  0.5*(tb37v+tb37h)
      diftb=tscale*(tb37v-tb37h)
      ibin=nint(0.2*(avgtb-125))
      if(diftb.lt.minpol37(ibin) .or. diftb.gt.maxpol37(ibin)) ioob=1

      return
      end
