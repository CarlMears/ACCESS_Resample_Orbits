!     feb 5, 2013 versin change 8/29/2013,  the indexing for 89b was changed to make it
!     adjacent to 89a. 'O:\amsr2\resampling\memo1.txt'. the 89 remapping routine has not yet 
!     been used for amsr2, so there is no impact.

!     april 25 2007 was changed on feb 29, 2008.  there was bug in resampling 89h.
!     the wrong ipol index was used (1 rather than 2) see comment cbug

!     the 12/10/2005 version was changed on april 25 2007.  the 89 ghz v vs. h renormalization coefs
!     are now different for amsr-e and amsr-j.  this has no effect on amsr-e

!     the 8/26/2005 version was changed on 12/10/2005.  actual values for pv0,pv1,ph0,ph1
!     were inserted. only amsr_sips_l1a_processing_v05 uses these values and it is yet to be run.  
    subroutine resample_ta(icase)
        use l2_module
        implicit none

        integer(4), parameter :: ncells_max=1073  !37*29

        integer(1) ibadv(maxcel,maxscan),ibadh(maxcel,maxscan)
        integer(4) list1(ncells_max),list2(ncells_max)
        integer(4) iscan,icel,jcel,jscan,ii,jj,num,iii,istart,ifreq,ifreqb,ichv,ichh,kchv,kchh,iset,icase 
        real(4) list3(ncells_max)
        real(8) sum_wtv,sum_wth,sum_tav,sum_tah
        real(4) tavx,tahx,xwt
        real(4) weight(maxcel,-18:18,-14:14,12)  !index order is a bit weird, but this is what peter has been using

        data istart/1/

        if(istart.eq.1) then
            istart=0
            open(3,file=filename_ta_resample_wts,status='old',form='binary')
            read(3) weight
            close(3)
        endif

        kchv=-1; iset=0  !iset denotes which set of resampling weights
    
        do 500 ifreqb=1,4
            if(ifreqb.eq.2) cycle  !7 ghz is not a base freq
            do 490 ifreq=1,6
                if(ifreq.lt.ifreqb) cycle
                if(ifreq.eq.4 .and. ifreqb.eq.4) cycle !no resampling of 19 onto 19

                iset=iset+1
                if(iset.gt.12) stop 'iset.gt.12 in resample_ta, pgm stopped'

                ichv=2*ifreq-1
                ichh=ichv+1
                
                kchv=kchv+2
                kchh=kchv+1

                !     ==========================================================================================================================
                !     ========================================== set quality flags for cells ===================================================
                !     ==========================================================================================================================

                do iscan=1,numscan
                    do icel=1,maxcel
                        tavx=ta(icel,iscan,ichv)
                        tahx=ta(icel,iscan,ichh)

                        if(tavx.lt.ta_earth_min(ichv) .or. tavx.gt.ta_earth_max(ichv) .or. iflag_cal(ichv,iscan).ne.0) then
                            ibadv(icel,iscan)=1
                        else
                            ibadv(icel,iscan)=0
                        endif
                        if(tahx.lt.ta_earth_min(ichh) .or. tahx.gt.ta_earth_max(ichh) .or. iflag_cal(ichh,iscan).ne.0) then
                            ibadh(icel,iscan)=1
                        else
                            ibadh(icel,iscan)=0
                        endif
                    enddo
                enddo 


                do icel=1,maxcel

                    num=0
                    do ii=-18,18
                        do jj=-14,14
                            xwt=weight(icel,ii,jj,iset)
                            if(xwt.ne.0.0) then
                                jcel=ii + icel
                                if(jcel.lt.1 .or. jcel.gt.maxcel) stop 'pgm stopped, error1 in resample_ta'
                                num=num+1
                                list1(num)=jcel
                                list2(num)=jj 
                                list3(num)=xwt
                            endif
                        enddo  !ii
                    enddo  !jj
    
                    do iscan=1,numscan

                        if(iscan.lt.8 .or. iscan.gt.numscan-7) then 
                            tar(icel,iscan,kchv)=ta(icel,iscan,ichv)
                            tar(icel,iscan,kchh)=ta(icel,iscan,ichh)
                            cycle
                        endif


                        sum_wtv=0; sum_wth=0;    sum_tav=0; sum_tah=0; iii=0
                        55  continue    
                        iii=iii+1
                        jcel= list1(iii)
                        jscan=list2(iii) + iscan
                        xwt=list3(iii)

                        if(jscan.ge.1 .and. jscan.le.numscan) then 
                            if(ibadv(jcel,jscan).eq.0) then  
                                sum_wtv=sum_wtv + xwt
                                sum_tav=sum_tav + xwt*ta(jcel,jscan,ichv)
                            endif
                            if(ibadh(jcel,jscan).eq.0) then  
                                sum_wth=sum_wth + xwt
                                sum_tah=sum_tah + xwt*ta(jcel,jscan,ichh)
                            endif
                        endif

                        if(iii.ne.num) goto 55     !this method is faster than a do loop

                        if(ibadv(icel,iscan).eq.0 .and. sum_wtv.ge.0.4) then
                            tar(icel,iscan,kchv)=sum_tav/sum_wtv
                        else
                            tar(icel,iscan,kchv)=ta(icel,iscan,ichv)
                        endif

                        if(ibadh(icel,iscan).eq.0 .and. sum_wth.ge.0.4) then
                            tar(icel,iscan,kchh)=sum_tah/sum_wth
                        else
                            tar(icel,iscan,kchh)=ta(icel,iscan,ichh)
                        endif

                    enddo  !iscan
                enddo  !icel

            490 continue
        500 continue

        if(icase.eq.1) then
            call resample_ta_hi_1
            call resample_ta_hi_2(ncells_max,weight)
        endif

        return
    end subroutine resample_ta



    subroutine resample_ta_hi_1
        use l2_module
        implicit none

        integer(4), parameter :: ncells_max=3770  !65*29*2
    !     following values come from 'o:\amsr_l2\v05\resampling\correlate_89a_89b_11_00416_10000.lis'
        real(4), parameter :: pv0=0.130671, pv1=0.993251, ph0=0.472994, ph1=0.992742

        integer(1) ibadv(maxcel_89,maxscan_89),ibadh(maxcel_89,maxscan_89)
        integer(4) list1(ncells_max),list2(ncells_max),list4(ncells_max)
        integer(4) iscan,icel,jcel,jscan,ii,jj,kk,num,iii,istart,ichv,ichh,kchv,kchh,numscan_89,jscan0,jjscan
        real(4) list3(ncells_max)
        real(8) sum_wtv,sum_wth,sum_tav,sum_tah
        real(4) tavx,tahx,xwt
        real(4) weight(maxcel,-32:32,-14:14,2)

        data istart/1/

        if(istart.eq.1) then
            istart=0
            if(size(tar, dim=3).ne.32) stop 'pgm stopped, tar is not set up for 89 remapping' 
            open(3,file=filename_ta89_resample_wts,status='old',form='binary')
            read(3) weight
            close(3)
        endif

        !     ==========================================================================================================================
        !     ========================================== set quality flags for cells ===================================================
        !     ==========================================================================================================================

        kchv=31
        kchh=32
        numscan_89=2*numscan


        do iscan=1,numscan
            do jcel=1,maxcel_89

                ichv=15
                ichh=16

                do jscan=2*iscan-1, 2*iscan       !89b then 89a, this is dif from amsre

                    tavx=ta_89(jcel,jscan,1)
                    tahx=ta_89(jcel,jscan,2)

                    if(tavx.lt.ta_earth_min(ichv) .or. tavx.gt.ta_earth_max(ichv) .or. iflag_cal(ichv,iscan).ne.0) then
                        ibadv(jcel,jscan)=1
                    else
                        ibadv(jcel,jscan)=0
                    endif

                    if(tahx.lt.ta_earth_min(ichh) .or. tahx.gt.ta_earth_max(ichh) .or. iflag_cal(ichh,iscan).ne.0) then
                        ibadh(jcel,jscan)=1
                    else
                        ibadh(jcel,jscan)=0
                    endif
                    ichv=ichv-2
                    ichh=ichh-2

                enddo  !jscan
            enddo  !jcel
        enddo  !iscan


        do icel=1,maxcel

            num=0
            do ii=-32,32
                do jj=-14,14
                    do kk=1,2    !1 is 89a, 2 is 89b
                        xwt=weight(icel,ii,jj,kk)
                        if(xwt.ne.0.0) then
                            jcel=ii + 2*icel - 1
                            if(jcel.lt.1 .or. jcel.gt.maxcel_89) stop 'pgm stopped, error3 in hi weights'
                            num=num+1
                            list1(num)=jcel
                        !    list2(num)=2*jj + kk - 1   !amsre convention of odd=89a, even=89b 
                            list2(num)=2*jj + 2 - kk   !amsr2 convention of even=89a, odd=89b 
                            list3(num)=xwt
                            list4(num)=kk
                        endif
                    enddo  !ii
                enddo  !jj
            enddo  !kk
    
            do iscan=1,numscan

                if(iscan.lt.8 .or. iscan.gt.numscan-7) then 
                    tar(icel,iscan,kchv)=-993
                    tar(icel,iscan,kchh)=-993
                    cycle
                endif

                jscan0=2*iscan-1

                sum_wtv=0; sum_wth=0; sum_tav=0; sum_tah=0; iii=0
                55  continue    
                iii=iii+1
                jcel= list1(iii)
                jscan=list2(iii) + jscan0
                xwt=  list3(iii)

                if(jscan.ge.2 .and. jscan.le.numscan_89) then !jscan=1 is the first 89b scan that is now thrown way
                
                    jjscan=jscan
                    if(list4(iii).eq.2) jjscan=jscan-2  ! 89b have been shift by 2 scans, see 'O:\amsr2\resampling\memo1.txt'
                    
                    if(jjscan.le.0 .or. jjscan.eq.numscan_89-1) stop 'indexing error in resample_ta_hi_1, pgm stopped'
                    
                    if(ibadv(jcel,jjscan).eq.0) then 
                        sum_wtv=sum_wtv + xwt
                        if(list4(iii).eq.1) then  !convert 89a to 89b
                            sum_tav=sum_tav + xwt*(pv1*ta_89(jcel,jjscan,1) + pv0)
                        else
                            sum_tav=sum_tav + xwt*ta_89(jcel,jjscan,1)
                        endif
                    endif

                    if(ibadh(jcel,jjscan).eq.0) then  
                        sum_wth=sum_wth + xwt
                        if(list4(iii).eq.1) then  !convert 89a to 89b
                            sum_tah=sum_tah + xwt*(ph1*ta_89(jcel,jjscan,2) + ph0)
                        else
                            sum_tah=sum_tah + xwt*ta_89(jcel,jjscan,2)
                        endif
                    endif
                endif

                if(iii.ne.num) goto 55     !this method is faster than a do loop

                if(sum_wtv.ge.0.1) then     !what should the limit really be
                tar(icel,iscan,kchv)=sum_tav/sum_wtv
                else
                tar(icel,iscan,kchv)=-993
                endif

                if(sum_wth.ge.0.1) then     !what should the limit really be
                tar(icel,iscan,kchh)=sum_tah/sum_wth
                else
                tar(icel,iscan,kchh)=-993.
                endif

            enddo  !iscan
        enddo  !icel

        return
    end subroutine resample_ta_hi_1

    subroutine resample_ta_hi_2(ncells_max,weight)
        use l2_module
        implicit none

        integer(4) ncells_max

        integer(1) ibadv(maxcel,maxscan),ibadh(maxcel,maxscan)
        integer(4) list1(ncells_max),list2(ncells_max)
        integer(4) iscan,icel,jcel,jscan,ii,jj,num,iii,ifreqb,ichv,ichh,kchv,kchh,iset 
        real(4) list3(ncells_max)
        real(8) sum_wtv,sum_wth,sum_tav,sum_tah
        real(4) tavx,tahx,xwt
        real(4) weight(maxcel,-18:18,-14:14,12)  !index order is a bit weird, but this is what peter has been using

        ichv=31
        ichh=32

        kchv=23
        do ifreqb=1,3
            if(ifreqb.eq.1) iset= 6    !37 onto  7
            if(ifreqb.eq.2) iset=10    !37 onto 11
            if(ifreqb.eq.3) iset=12    !37 onto 19
    
            kchv=kchv+2
            kchh=kchv+1

            !     ==========================================================================================================================
            !     ========================================== set quality flags for cells ===================================================
            !     ==========================================================================================================================

            do iscan=1,numscan
                do icel=1,maxcel
                    tavx=tar(icel,iscan,ichv)
                    tahx=tar(icel,iscan,ichh)

                    if(tavx.lt.ta_earth_min(15) .or. tavx.gt.ta_earth_max(15)) then
                        ibadv(icel,iscan)=1
                    else
                        ibadv(icel,iscan)=0
                    endif
                    if(tahx.lt.ta_earth_min(16) .or. tahx.gt.ta_earth_max(16)) then
                        ibadh(icel,iscan)=1
                    else
                        ibadh(icel,iscan)=0
                    endif
                enddo
            enddo 


            do icel=1,maxcel

                num=0
                do ii=-18,18
                    do jj=-14,14
                        xwt=weight(icel,ii,jj,iset)
                        if(xwt.ne.0.0) then
                            jcel=ii + icel
                            if(jcel.lt.1 .or. jcel.gt.maxcel) stop 'pgm stopped, error1 in resample_ta'
                            num=num+1
                            list1(num)=jcel
                            list2(num)=jj 
                            list3(num)=xwt
                        endif
                    enddo  !ii
                enddo  !jj
    
                do iscan=1,numscan
                    if(iscan.lt.8 .or. iscan.gt.numscan-7) then 
                        tar(icel,iscan,kchv)=tar(icel,iscan,ichv)
                        tar(icel,iscan,kchh)=tar(icel,iscan,ichh)
                        cycle
                    endif

                    sum_wtv=0; sum_wth=0;    sum_tav=0; sum_tah=0; iii=0
                    55  continue    
                        iii=iii+1
                        jcel= list1(iii)
                        jscan=list2(iii) + iscan
                        xwt=list3(iii)

                        if(jscan.ge.1 .and. jscan.le.numscan) then 
                            if(ibadv(jcel,jscan).eq.0) then  
                                sum_wtv=sum_wtv + xwt
                                sum_tav=sum_tav + xwt*tar(jcel,jscan,ichv)
                            endif
                            if(ibadh(jcel,jscan).eq.0) then  
                                sum_wth=sum_wth + xwt
                                sum_tah=sum_tah + xwt*tar(jcel,jscan,ichh)
                            endif
                        endif

                    if(iii.ne.num) goto 55     !this method is faster than a do loop

                    if(ibadv(icel,iscan).eq.0 .and. sum_wtv.ge.0.4) then
                        tar(icel,iscan,kchv)=sum_tav/sum_wtv
                    else
                        tar(icel,iscan,kchv)=tar(icel,iscan,ichv)
                    endif

                    if(ibadh(icel,iscan).eq.0 .and. sum_wth.ge.0.4) then
                        tar(icel,iscan,kchh)=sum_tah/sum_wth
                    else
                        tar(icel,iscan,kchh)=tar(icel,iscan,ichh)
                    endif

                enddo  !iscan
            enddo  !icel

        enddo ! ifreqb
        return
    end subroutine resample_ta_hi_2




