!     same as 'o:\amsr2\routines\tatb_routines.f' dated 8/22/2013 except delta_apc,chi_apc are read in by 
!     the routine 'o:\amsr2\routines\readin_data_files.f'

!     feb 5 2013 version changed 8/22/2013.  tbmin and tbmax are now readin by 'o:\amsr2\routines\readin_data_files.f'

    subroutine ckocnta(ifreq, tamin,tamax)
        use l2_module                                             
        implicit none

        integer(4) ifreq
        real(4) tamin(2),tamax(2)
        real(4) tb(2)

        if(ifreq.gt.nfreq) stop 'ifreq oob in ckocnta, pgm stopped'
        tb(1:2)=tbmin(1:2,ifreq)
        call tatb(1,ifreq, tamin,tb)   !tb converted to ta for all channels
        tamax =280.
        return
    end


 
    subroutine tatb(icase,ifreq, tax,tb)
        use l2_module
        implicit none
    
        integer(4) icase,ifreq,jfreq,istart,ich,jch
        real(4) tax(2),tb(2)
        real(4) qvv(nfreq),qhv(nfreq),qov(nfreq),qhh(nfreq),qvh(nfreq),qoh(nfreq)
        real(4) avv(nfreq),ahv(nfreq),aov(nfreq),ahh(nfreq),avh(nfreq),aoh(nfreq)
        real(4) xfac,tv,th
    
    !     data initialization
    
        data istart/1/
    
    !     begin execution
 
        if(istart.eq.1) then
            istart=0
            do jfreq=1,nfreq
        
                qvv(jfreq)=(1-delta_apc(1,jfreq))/(1.+chi_apc(1,jfreq))
                qhv(jfreq)=chi_apc(1,jfreq)*qvv(jfreq)
                qov(jfreq)=delta_apc(1,jfreq)*tcld_plk(jfreq)
                qhh(jfreq)=(1-delta_apc(2,jfreq))/(1.+chi_apc(2,jfreq))
                qvh(jfreq)=chi_apc(2,jfreq)*qhh(jfreq)
                qoh(jfreq)=delta_apc(2,jfreq)*tcld_plk(jfreq)
        
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
    
            if(tb(ich).lt.55 .or. tb(ich).gt.330 .or. &
               tb(jch).lt.55 .or. tb(jch).gt.330) then
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
