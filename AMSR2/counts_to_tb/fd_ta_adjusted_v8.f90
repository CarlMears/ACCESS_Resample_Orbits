!      5/23/2017 updated 6/2/2017.  i rederived the dta_zang correction and turned in on for all channels
!      this routine is now back to its 3/36/2017 state


!      3/26/2017 updated 5/23/2017.  the dta_zang correction is now only done for 89 ghz.


!     i upgraded this to v8 on 3/26/2017 when i added dta_zang

!     below are v7 comments except for the first two regarding turning off the lna correction for later orbits

!     2/28/2017 change 3/2/2017
!     since the lna correction is no longer done after orbit 11000, i turned off the write statement:
!     write(*,1001) alpha,beta
!     1001 format('corr_tab_11vlna table value missing in fd_dta11v_lna ',2f8.3)
!     i no longer need to keep track of this


!     7/11/2014 chnage 2/28/2017. i stopped using lna correction after orbit 11000. see 'O:\amsr2\abs_cal\memo17.txt'

!     9/19/2013 version updated 7/11/2014.  correction for 11v lna put in and subroutine fd_dta11v_lna added

!     9/9/2013 version changed 9/19/2013.  iminc,imaxc check replace by .ne.-32768 check

!     8/29/2013 version changed 9/9/2013.  this was a major change.  the computation of thot was removed
!     and put into fd_thot.  The computation of tcold including the spillover part was removed and put into
!     fd_tcold.  In doing this you must take care about the qc flag iflag_cal.  This routine had been intializing
!     it to zero and then setting all the qc bits.  Now fd_thot set it to zero and sets some bits.  so fd_thot
!     now must be called first.  Then call fd_cold.  The previous version of this routine is save with the date added
!     as a suffix


!     6/26/2013 version changed 8/29/2013.  tamea_cband renamed ta_spill and now includes all 16 channels not just ich=1,2
!     iopt_channels.eq.0 now does first 6 channels rather than first 2
!     jscan1 is now 2*iscan-3 to make 89b adjacent and before to 89a

!     march 13 2013 version updated on 6/26/2013
!     bug found:  xscan was not being applied

!     march 6 2013 version changed march 13 2013
!     the ssmi form of xscan table is now used.  this change has no effect on previous processing
!     because xscan had been set to zero
!     also the xscan table is now readin by subroutine 'o:\amsr2\routines\readin_data_files.f'

!     feb 8 2013 version changed march 6 2013
!     1. beta_term array no longer need and is removed
!     2. acoef_nl replaces beta_quad and is readin via 'o:\amsr2\routines\readin_data_files.f'
!     3. hot target offset thot_offset(nfreq) is added


!     iopt_cld_spill=0 uses climate 6 ghz ta map; iopt_cld_spill=1 uses actual ta measurements

!     if(iopt_channels.eq.0) only the two 6v and 6h channels are process.  this is all mk_7ghz_tamap needs
!     otherwise all channels are process

    
    subroutine fd_ta_adjusted(num_channels)
        use l2_module
        implicit none

        integer(4) num_channels

        integer(4) iscan,jscan1,jscan2,ifreq,ipol,ich,icel,jcel,ncel,ibad,icount
        real(4) tcold,thot,dif_count,offset,gain,tax,cc,ch,dta11v_lna
        real(8) cx,cxx(5)
        
        integer(4) i1,i2
        real(4) brief,a1,a2,qzang

        ta_89=-999  !i put the default here because of the jscan1,jscan2 complication

        do iscan=1,numscan

            jscan1=2*iscan-3  !89b assign to this to make it adjacent and before to 89a
            jscan2=2*iscan    

            !     ============================================================================================
            !     ================================== see if scan is ok =======================================
            !     ============================================================================================

            if(iflag_l0(iscan).ge.8) then ! no science, spacecraft or geoloc data
                ta(:,iscan,:)= -999
                cycle
            endif
    
            call fd_dta11v_lna(corr_tab_11vlna, alpha_sun(iscan), beta_sun(iscan),   dta11v_lna)
    

            !     =================================================================================================
            !     =====================================  dta_zang interpolation coefs =============================
            !     =================================================================================================
            qzang=zang(iscan)
            if(qzang.lt.0 .or. qzang.gt.360.) stop 'pgm stopped, zang oob in fd_ta_adjust'
            if(qzang.gt.359.99) qzang=359.99
            brief=qzang/3.6
            i1=brief
            i2=i1+1
            a1=i2-brief
            a2=1-a1
            if(i2.eq.100) i2=0
            if(i1.lt.0 .or. i2.gt.99) stop 'i oob in fd_ta_adjusted, pgm stopped'


            !     ============================================================================================
            !     ================================ loop through channels =====================================
            !     ============================================================================================
            do ich=1,num_channels

                if(ich.le.12) then
                    ncel=maxcel
                else
                    ncel=maxcel_89
                endif
      
                ifreq=1 + int((ich-1)/2)
                ipol=ich -2*int((ich-1)/2)

                tcold=ta_cold(ich,iscan)
                thot= ta_hot( ich,iscan)

                ibad=0


                !     bit 14 is not currently used for amsr2
                if(iflag_rx(ich,iscan).ne.0) then    
                    iflag_cal(ich,iscan)=ibset(iflag_cal(ich,iscan),13)    
                    ibad=3
                    goto 50
                endif

                if(num_count(1,ich,iscan).eq.0 .or. num_count(2,ich,iscan).eq.0) then 
                    iflag_cal(ich,iscan)=ibset(iflag_cal(ich,iscan),12)    !< one or zero cold or hot counts, no offset,gain found
                    ibad=2
                    goto 50
                endif
    
                !     calibration scan averaging designed to always give at least 100 obs when possible
                !     so the follwing is a pretty serious error, but i still go ahead and find ta
                if(num_count(1,ich,iscan).lt.50) iflag_cal(ich,iscan)=ibset(iflag_cal(ich,iscan),8)    !too few cold counts
                if(num_count(2,ich,iscan).lt.50) iflag_cal(ich,iscan)=ibset(iflag_cal(ich,iscan),7)    !too few hot counts
    
                cc=avg_count(1,ich,iscan)
                ch=avg_count(2,ich,iscan)

                dif_count=ch-cc
                if(dif_count.lt.1) then
                    iflag_cal(ich,iscan)=ibset(iflag_cal(ich,iscan),10)    !cold counts = or > hot counts, no offset,gain found
                    ibad=1
                    goto 50
                endif

                if(dif_count.le. 100)                                            iflag_cal(ich,iscan)=ibset(iflag_cal(ich,iscan),6)  !very small dif_count
                if(dif_count.le.1500)                                            iflag_cal(ich,iscan)=ibset(iflag_cal(ich,iscan),5)  !small dif_count

                offset=(ch*tcold-cc*thot)/dif_count
                gain =(thot-tcold)/dif_count

                50    continue

                do icel=1,ncel

                    if(ich.le.12) then  !find jcel for xscan correction
                        jcel=2*icel-1
                    else
                        jcel=icel
                    endif
    
                    if(ibad.ne.0) then
                        if(ibad.eq.1) tax=-996.
                        if(ibad.eq.2) tax=-997.
                        if(ibad.eq.3) tax=-998.
                    else  ! cal data ok, get earth ta
                        icount=earthcounts(ich,icel,iscan)
    
                        if(icount.ne.-32768) then  !earth counts in bounds
    
                            cx=(icount-cc)/dif_count
                            cxx(1)=cx
                            cxx(2)=cx**2
                            cxx(3)=cx**3
                            cxx(4)=cx**4
                            cxx(5)=cx**5
                            tax=offset + gain*icount  - dot_product(cxx,acoef_nl(:,ich)) - 4.*(a1*dta_zang(i1,ich) + a2*dta_zang(i2,ich))*cx*(1-cx)
                            !x    tax=offset + gain*icount  - dot_product(cxx,acoef_nl(:,ich))
                            !x    if(ich.ge.13)  tax=tax - 4.*(a1*dta_zang(i1,ich) + a2*dta_zang(i2,ich))*cx*(1-cx)  !inserted 5/23/2017

                            !     this is correction for 11v lna problem
                            !x    if(ich.eq.5 .and. orbit(iscan).gt.6000 .and. tax.le.thot) then
                            if(ich.eq.5 .and. orbit(iscan).gt.6000 .and. orbit(iscan).lt.11000 .and. tax.le.thot) then  !inserted 2/28/2017
                                tax=tax - dta11v_lna*(tax-tcold)*(thot-tax)/((163.-tcold)*(thot-163.))  !163 is mean value of 11v in vacinity of problem
                            endif

        
                            !bug    if(tax.lt.0) then !good ta
                            if(tax.gt.0) then !good ta
                                tax=(tax - tcld_plk(ifreq)*xscan_tab(ich,jcel))/(1-xscan_tab(ich,jcel)) 
                            else
                                if(tax.lt.-994) tax=-994.  !so it does not conflict with other error codes above
                            endif
                        else
                            tax=-995  !cal data was good but earth counts oob
                        endif
                    endif  !if(ibad.ne.0)
    
                    if(ich.le.12)   ta(icel,iscan,ich) =tax
                    if(ich.eq.13) ta_89(icel,jscan2,1) =tax
                    if(ich.eq.14) ta_89(icel,jscan2,2) =tax
                    if(jscan1.ge.1) then
                        if(ich.eq.15) ta_89(icel,jscan1,1) =tax
                        if(ich.eq.16) ta_89(icel,jscan1,2) =tax
                    endif

                enddo  !icel
            enddo  !ich
        enddo  !iscan
        return
    end subroutine fd_ta_adjusted


    subroutine fd_dta11v_lna(corr_tab_11vlna,alpha_in,beta,dta11v_lna)
        implicit none
        
        real(4) corr_tab_11vlna(0:359,0:359)
        real(8) alpha_in,beta
        real(4) dta11v_lna

        integer(4) i1,i2,j1,j2
        real(4) brief,a1,a2,b1,b2
        real(4) wt
        real(8) alpha
      
        wt=0.5*(7-abs(beta-107.))  !wt taper correction to zero at edge
        if(wt.gt.1) wt=1
        
        if(wt.le.0) then
            dta11v_lna=0
            return
        endif

        alpha=alpha_in
        if(alpha.lt.0) alpha=alpha+360.d0  !amsr geoloc return -180 to 180
        if(alpha.lt.10 .or. alpha.gt.40) then
            dta11v_lna=0
            return
        endif
    
        brief=alpha
        i1=brief
        i2=i1+1
        a1=i2-brief
        a2=1-a1
        if(i1.lt.0 .or. i2.gt.359) then
            print *,alpha,beta,i1,i2
            stop 'i oob in ffd_dta11v_lna, pgm stopped'
        endif
        brief=2*beta
        j1=brief
        j2=j1+1
        b1=j2-brief
        b2=1-b1

        if(j1.lt.0 .or. j2.gt.359) stop 'j oob in fd_dta11v_lna, pgm stopped'

        if(minval(corr_tab_11vlna(i1:i2,j1:j2)).gt.-998) then
            dta11v_lna=a1*(b1*corr_tab_11vlna(i1,j1) + b2*corr_tab_11vlna(i1,j2)) + &
                       a2*(b1*corr_tab_11vlna(i2,j1) + b2*corr_tab_11vlna(i2,j2))
        else
            i1=nint( alpha)
            j1=nint(2*beta)
            if(i1.lt.0 .or. i1.gt.359) stop 'ialpha oob in fd_dta11v_lna, pgm stopped'
            if(j1.lt.0 .or. j1.gt.359) stop ' ibeta oob in fd_dta11v_lna, pgm stopped'
            dta11v_lna=corr_tab_11vlna(i1,j1)
        endif

        if(dta11v_lna.lt.-998) then
            !x    write(*,1001) alpha,beta
            !x 1001 format('corr_tab_11vlna table value missing in fd_dta11v_lna ',2f8.3)
            dta11v_lna=0
            return
        endif
    
        dta11v_lna=wt*dta11v_lna  !tape to zero at edges
    return
    end subroutine fd_dta11v_lna
