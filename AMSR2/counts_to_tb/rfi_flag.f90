    !     this is the dec 19 2013 version of 'x:\chelle\amsr2\rfi\rfi_flag.f' with the following changes
    !     1.  routine rfi_flag_nodif is removed.  not used
    !     2.  standard istart for initialization is used
    !     3.  check for existence of ep was wrong.  it should be .lt.999 not >-99

    subroutine rfi_flag(iasc,xlat,xlon,sstvl,sstlo,sstcl,sstrg, &
                         winvl,winlo,winmd,winncep,geo_xlat,geo_xlon,irfi_flag)
        !     clg added change for backup check.  in coastal regions there is no sst data, but winds can still be bad.
        !     so need backup check in case primary check isn't availalbe
        !     irfi_flag is flag
        !     iarea doesn't do max/min test, just looks for area

        implicit none
        real sstvl,sstlo,sstcl,sstrg,winvl,winlo,winmd,rdif(11),xlat,xlon
        real geo_xlat,geo_xlon,winncep,rdif2
        integer(4) irfi_flag,iasc,inum_rfi,ic,i
        integer(4),dimension(100)::lyr_start,idyjl_start,lyr_end,idyjl_end,icheck,iasc_rfi,iground,icheck2
        real(8),dimension(100)::xlon_min,xlon_max,xlat_min,xlat_max,rvalue_min,rvalue_max,rvalue_min2,rvalue_max2
        integer(4) imask(1440,720,2),ilon,ilat, istart
      
        data istart/1/
      
        if(istart.eq.1) then
            istart=0
            print*, 'initialize RFI'
            open(unit=3,file='x:\chelle\amsr2\tables\AMSR2_rfi_inputs.txt',form='formatted')
            inum_rfi=0
            do i=1,100
                read(3,'(4I5,4F6.1,I3,2F6.1,I3,2F6.1,2I3)',end=102) lyr_start(i),idyjl_start(i),lyr_end(i),idyjl_end(i), &
                            xlon_min(i),xlon_max(i),xlat_min(i),xlat_max(i),  &
                            icheck(i),rvalue_min(i),rvalue_max(i),  &
                            icheck2(i),rvalue_min2(i),rvalue_max2(i),  &
                            iground(i),iasc_rfi(i)
                inum_rfi=inum_rfi + 1
            enddo

            102   close(3)
            print*, 'initialize RFI AREAS'
            open(unit=3,file='x:\chelle\amsr2\tables\rfi_mask_final.txt',form='formatted')
            read(3,'(2073600I3)') imask
            close(3)
        endif
      
      
        irfi_flag=0
        rdif=0
        
        if(sstvl.lt.999 .and. sstlo.lt.  999) rdif( 1)=sstvl-sstlo
        if(sstvl.lt.999 .and. sstcl.lt.  999) rdif( 2)=sstvl-sstcl
        if(sstvl.lt.999 .and. sstrg.lt.  999) rdif( 3)=sstvl-sstrg
        if(sstlo.lt.999 .and. sstcl.lt.  999) rdif( 4)=sstlo-sstcl
        if(sstlo.lt.999 .and. sstrg.lt.  999) rdif( 5)=sstlo-sstrg
        if(winvl.lt.999 .and. winlo.lt.  999) rdif( 6)=winvl-winlo
        if(winvl.lt.999 .and. winmd.lt.  999) rdif( 7)=winvl-winmd
        if(winvl.lt.999 .and. winncep.lt.999) rdif( 8)=winvl-winncep
        if(winlo.lt.999 .and. winmd.lt.  999) rdif( 9)=winlo-winmd
        if(winlo.lt.999 .and. winncep.lt.999) rdif(10)=winlo-winncep
        if(winmd.lt.999 .and. winncep.lt.999) rdif(11)=winmd-winncep


        ilat=nint((xlat+89.875)/.25)+1
        ilon=nint((xlon-.125)/.25)+1
        if(ilon<   1)ilon=ilon+1440
        if(ilon>1440)ilon=ilon-1440

        !      if(ilon==876 .and. ilat==589 .and. iasc==2)
        !     1  write(*,'(4I4,8F8.2)') ilon,ilat,imask(ilon,ilat,iasc),icheck(1),rdif(7),rdif(9),
        !     2      xlon_min(1),geo_xlon,xlon_max(1),
        !     2      xlat_min(1),geo_xlat,xlat_max(1)
        !      if(ilon==873 .and. ilat==589 .and. iasc==2)
        !     1  write(*,'(4I4,8F8.2)') ilon,ilat,imask(ilon,ilat,iasc),icheck(1),rdif(7),rdif(9),
        !     2      xlon_min(1),geo_xlon,xlon_max(1),
        !     2      xlat_min(1),geo_xlat,xlat_max(1)
      
        if(imask(ilon,ilat,iasc)==0)return  !not in region with RFI

        do ic=1,inum_rfi
            rdif2=rdif(icheck(ic))
            if(iasc/=iasc_rfi(ic))cycle
            if(iasc==1 .and. iground(ic)==0 .and. xlat>0)cycle !asc souther hemi for sat rfi 
            if(iasc==2 .and. iground(ic)==0 .and. xlat<0)cycle  !dsc northern hemisphere for sat rfi
            if(iground(ic)==0) then !satellite rfi so use geo_sync lat/lon
                if(rdif2>=rvalue_min(ic) .and. rdif2<=rvalue_max(ic) .and. &
                geo_xlon>=xlon_min(ic) .and. geo_xlon<=xlon_max(ic) .and.  &
                geo_xlat>=xlat_min(ic) .and. geo_xlat<=xlat_max(ic)) irfi_flag=ic
                ! iasc1 check southern hemisphere higher varibabilith, so don't flag except near s.america
                if(irfi_flag>0 .and. iasc==1 .and. xlat<-40 .and. xlon<275)irfi_flag=0
                if(irfi_flag>0 .and. iasc==1 .and. xlat<-40 .and. xlon>324)irfi_flag=0
            else   !ground based rfi so use geographic lat/lon
                if(rdif2>=rvalue_min(ic) .and. rdif2<=rvalue_max(ic) .and. &
                xlon>=xlon_min(ic) .and. xlon<=xlon_max(ic) .and. &
                xlat>=xlat_min(ic) .and. xlat<=xlat_max(ic)) irfi_flag=i!      
            endif
        enddo   
      
    end subroutine rfi_flag
      
      

