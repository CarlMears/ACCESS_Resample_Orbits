!     1/30/2014 changed 8/28/2017.  this routine now compares ta_hot with thot_jaxa
!     check is only done for orbits after 21000. orbits prior to this have
!     been verified to have no problems. see 'o:\amsr2\abs_cal\memo22.txt'     


!     12/5/2013 version changed 1/30/2014.  a bug was found. see cbug1
!     i verifed that this check never occured, so its not problem
!     also var.gt.87.5 is replaced by var.gt.therm_var, where therm_var is 90
!     so this has no effect since this is making the check even less strict

!     this 12/5/2013 version is the same as the sept 7 2013 version,
!     a temporary change was made and then removed.


        subroutine fd_thot
            use l2_module
            implicit none

            !     t1 is for orbits 1-5249 and 20811-later, t2 is for orbits 5250-20810
            real(4), parameter :: &
            t1(nch)=(/0.113, 0.113, 0.055, 0.055,-0.406,-0.406,-0.810,-0.810,-0.925,-0.925,-1.061,-1.061,-1.205,-1.205,-1.205,-1.205/)
            real(4), parameter :: &
            t2(nch)=(/0.056, 0.056,-0.002,-0.002,-0.463,-0.463,-0.867,-0.867,-0.982,-0.982,-1.118,-1.118,-1.262,-1.262,-1.262,-1.262/)

            integer(4) iscan,ifreq,ich
            real(4) xmean,var
            real(4) dteff
            real(4) thot_therm(8)
            integer(4) ierr_dteff,igo

            iflag_cal=0

            do iscan=1,numscan
            
                !     ============================================================================================
                !     ================================== see if scan is ok =======================================
                !     ============================================================================================
                if(iflag_l0(iscan).ge.8) then ! no science, spacecraft or geoloc data
                    ta_hot( :,iscan)= -999
                    do ich=1,nch
                        iflag_cal(ich,iscan)=ibset(iflag_cal(ich,iscan),15)    ! no science or geolocation data 
                    enddo
                    cycle
                endif
    
                !     ============================================================================================
                !     ======================== check termistors and get static teff ==============================
                !     ============================================================================================
                xmean=sum(therm(:,iscan))/10.
                var=dot_product(therm(:,iscan),therm(:,iscan))/10. - xmean**2
                if(minval(therm(:,iscan)).lt.therm_min .or. maxval(therm(:,iscan)).gt.therm_max .or. var.gt.therm_var) then 
                    !bug1 iflag_cal(ich,iscan)=ibset(iflag_cal(ich,iscan),9)  !therm data suspect
                    do ich=1,nch
                    iflag_cal(ich,iscan)=ibset(iflag_cal(ich,iscan),9)  !therm data suspect
                    enddo

                endif
      
                thot_therm=xmean  !currently all freqs use same thermistor average
    
                !     ============================================================================================
                !     ================================ loop through channels =====================================
                !     ============================================================================================
                do ich=1,nch
                    ifreq=1 + int((ich-1)/2)
                    
                    !     turn this off for now
                    !    if(iflag_l0(iscan).eq.0 .or. iflag_l0(iscan).eq.1)    then
                    !    call fd_dteff(filename_tef,ifreq,alpha_sun(iscan),beta_sun(iscan), dteff,ierr_dteff)
                    !    else
                    dteff=0
                    ierr_dteff=0
                    !    endif

                    if(ierr_dteff.ne.0) then
                        write(*,1001) iscan,ich,scrpy(:,iscan),alpha_sun(iscan),beta_sun(iscan)
                        1001 format('dteff outside map ',2i6,5f10.5)
                        write(*,*) 'enter 1 to continue, enter 2 to stop'
                        read(*,*) igo
                        if(igo.ne.1) stop 'pgm stopped by user'
                    endif

                    ta_hot(ich,iscan)=thot_therm(ifreq) + thot_offset(ifreq) + dteff

                    !     check is only done for orbits after 21000. orbits prior to this have
                    !     been verified to have no problems. see 'memo22.txt'     
                    if(orbit(iscan).lt.21000)         cycle
                    if(btest(iflag_cal(ich,iscan),9)) cycle
                    
                    if(abs(thot_therm(ifreq) + thot_offset(ifreq) - thot_jaxa(ich,iscan) - t1(ich)).gt.0.002) then
                        write(*,'(f10.2,i4,i3,16f8.3)')  orbit(iscan),iscan,ich, thot_therm(ifreq) + thot_offset(ifreq), thot_jaxa(ich,iscan)
                        write(*,*) 'problem with ta_hot in routine fd_hot, enter 0 to continue, enter 1 to stop program'
                        read(*,*) igo
                        if(igo.eq.1)  stop 'program stopped by user'
                        iflag_cal(ich,iscan)=ibset(iflag_cal(ich,iscan),9)  !therm data suspect  
                    endif    
                enddo  !ich
            enddo  !iscan

            return
        end subroutine fd_thot