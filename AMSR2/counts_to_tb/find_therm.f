!     8/27/2017 changed 11/27/2017.
!     another change in the jaxa therm occured and the following statement was inserted
!     if(orbit.gt.29393.5d0)  cspc=-23.
!     the program that call this routine, read_rss_l1a.f, was also changed.
!     before this 11/27/2017 changed, read_rss_l1a.f sent orbit(iscan)


!     12/5/2013 changed 8/27/2017.  after orbit 20810 changes back to original value of -13. 'see O:\amsr2\abs_cal\hot_load_step.pptx'


!     3/8/2013 version changed 12/5/2013. orbit is now an input argument and cspc is assigned
!     a new value for orbits after 5250.  see 'O:\amsr2\abs_cal\memo10.txt'

!     jan 25 2013 version changed on 3/8/2013
!     check for erroneous icount_spc was changed because for v1 processing an unsigned integer(4) is used for spc counts
!     rather than a signed integer(2).  icount_spc is now integer(4)
!     v0 processing show that sps counts were never .gt.8191 except when they were 65535
!     i assume the same is true for the spc counts that are used for find therm.  a value of 8191 would be way ooob
      subroutine find_therm(orbit,itype,icount_spc, therm)
      implicit none
      
      real(8) orbit
      integer(4) itype
      integer(4) icount_spc
      real(4) therm
      
      real(8) w0(10),w1(10),w2(10),r
      real(4) cspc
      
      data w0/-242.10300,-241.64400,-241.68400,-241.65700,-241.67800,-241.46000,-241.49700,-241.72200,-241.61200,-241.59000/
      data w1/  0.116612,  0.116174  ,0.116131,  0.116165,  0.116180,  0.115852,  0.115839,  0.116029,  0.115961,  0.115949/
      data w2/ 0.2328650, 0.2426100, 0.2420060, 0.2422170, 0.2426980, 0.2465540, 0.2453660, 0.2425070, 0.2414040, 0.2438840/
!     data cspc/-13/ ! cspc: -13 (spc-a), 23 (spc-b)  # spc-a is being used as of december 27, 2012

!     commented out 8/27/2017
!     if(orbit.lt.5250) then
!     cspc=-13
!     else
!     cspc=-20.7
!     endif
      
      if(orbit.ge.5250 .and. orbit.le.20810) then
      cspc=-20.7
      else
      cspc=-13
      endif
      
      if(orbit.gt.29393.5d0)  cspc=-23.
      
      if(icount_spc.lt.0 .or. icount_spc.gt.8191) then  !it should never be .lt.0, but 65535 values do occur when data is missing
      therm=0
      return
      endif
      
      r = 0.196193578d0*(icount_spc-cspc) + 1656.203847756d0
      therm=w0(itype) + w1(itype)*r + 1.e-5*w2(itype)*r*r  + 273.15
      
      return
      end
     


