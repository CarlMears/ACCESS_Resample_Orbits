 module filename_routines
 
 contains
    subroutine get_l1a_filename(pathname_l1a,ksat,iorbit, filename_l1a,iexist) 
        implicit none

        character(*) :: pathname_l1a
        character(*) ::filename_l1a
        integer(4) :: ksat,iorbit,iexist
        logical(4)  :: lexist
        character(100) :: adir
        character(13) :: aorbit

        if(ksat.ne.32 .and. ksat.ne.33) stop 'ksat oob in get_l1a_filename'

        if(iorbit.lt.1 .or. iorbit.gt. 99999) stop 'error1 in get_l1a_filename, pgm stopped'

        if(iorbit.ge.    1 .and. iorbit.le. 5000)   adir=trim(pathname_l1a)//'r00001_05000\r'
        if(iorbit.ge. 5001 .and. iorbit.le.10000)   adir=trim(pathname_l1a)//'r05001_10000\r'
        if(iorbit.ge.10001 .and. iorbit.le.15000)   adir=trim(pathname_l1a)//'r10001_15000\r'
        if(iorbit.ge.15001 .and. iorbit.le.20000)   adir=trim(pathname_l1a)//'r15001_20000\r'
        if(iorbit.ge.20001 .and. iorbit.le.25000)   adir=trim(pathname_l1a)//'r20001_25000\r'
        if(iorbit.ge.25001 .and. iorbit.le.30000)   adir=trim(pathname_l1a)//'r25001_30000\r'
        if(iorbit.ge.30001 .and. iorbit.le.35000)   adir=trim(pathname_l1a)//'r30001_35000\r'
        if(iorbit.ge.35001 .and. iorbit.le.40000)   adir=trim(pathname_l1a)//'r35001_40000\r'
        if(iorbit.ge.40001 .and. iorbit.le.45000)   adir=trim(pathname_l1a)//'r40001_45000\r'
        if(iorbit.ge.45001 .and. iorbit.le.50000)   adir=trim(pathname_l1a)//'r45001_50000\r'
        if(iorbit.ge.50001 .and. iorbit.le.55000)   adir=trim(pathname_l1a)//'r50001_55000\r'
        if(iorbit.ge.55001 .and. iorbit.le.60000)   adir=trim(pathname_l1a)//'r55001_60000\r'
        if(iorbit.ge.60001 .and. iorbit.le.65000)   adir=trim(pathname_l1a)//'r60001_65000\r'
        if(iorbit.ge.65001 .and. iorbit.le.70000)   adir=trim(pathname_l1a)//'r65001_70000\r'
        if(iorbit.ge.70001 .and. iorbit.le.75000)   adir=trim(pathname_l1a)//'r70001_75000\r'
        if(iorbit.ge.75001 .and. iorbit.le.80000)   adir=trim(pathname_l1a)//'r75001_80000\r'
        if(iorbit.ge.80001 .and. iorbit.le.85000)   adir=trim(pathname_l1a)//'r80001_85000\r'
        if(iorbit.ge.85001 .and. iorbit.le.90000)   adir=trim(pathname_l1a)//'r85001_90000\r'
        if(iorbit.ge.90001 .and. iorbit.le.95000)   adir=trim(pathname_l1a)//'r90001_95000\r'
        if(iorbit.ge.95001 .and. iorbit.le.99999)   adir=trim(pathname_l1a)//'r95001_99999\r'

        write(aorbit,9001) iorbit
        9001 format(i5.5,'.gz') 
        filename_l1a=trim(adir)//trim(aorbit)

        inquire(file=filename_l1a,exist=lexist)
        if(lexist) then
            iexist=1
        else
            iexist=0
        endif

    return
    end subroutine get_l1a_filename

    ! subroutine get_ssmi_l1c_filename(pathname_l1c,ksat,iorbit, filename_l1c)
    !     implicit none

    !     character(*),intent(in)    :: pathname_l1c
    !     integer(4),intent(in)      :: ksat
    !     integer(4),intent(in)      :: iorbit
    !     character(*),intent(out)   :: filename_l1c


    !     character(100) adir
    !     character(15) afile

    !     if(iorbit.lt.1 .or. iorbit.gt. 115000) then
    !         stop 'orbit out of bounds in get_ssmi_l1c_filename_base'
    !     endif

  
    !     if(iorbit.ge.    1 .and. iorbit.le. 5000)   adir=trim(pathname_l2b)//'r000001_005000\r'
    !     if(iorbit.ge. 5001 .and. iorbit.le.10000)   adir=trim(pathname_l2b)//'r005001_010000\r'
    !     if(iorbit.ge.10001 .and. iorbit.le.15000)   adir=trim(pathname_l2b)//'r010001_015000\r'
    !     if(iorbit.ge.15001 .and. iorbit.le.20000)   adir=trim(pathname_l2b)//'r015001_020000\r'
    !     if(iorbit.ge.20001 .and. iorbit.le.25000)   adir=trim(pathname_l2b)//'r020001_025000\r'
    !     if(iorbit.ge.25001 .and. iorbit.le.30000)   adir=trim(pathname_l2b)//'r025001_030000\r'
    !     if(iorbit.ge.30001 .and. iorbit.le.35000)   adir=trim(pathname_l2b)//'r030001_035000\r'
    !     if(iorbit.ge.35001 .and. iorbit.le.40000)   adir=trim(pathname_l2b)//'r035001_040000\r'
    !     if(iorbit.ge.40001 .and. iorbit.le.45000)   adir=trim(pathname_l2b)//'r040001_045000\r'
    !     if(iorbit.ge.45001 .and. iorbit.le.50000)   adir=trim(pathname_l2b)//'r045001_050000\r'
    !     if(iorbit.ge.50001 .and. iorbit.le.55000)   adir=trim(pathname_l2b)//'r050001_055000\r'
    !     if(iorbit.ge.55001 .and. iorbit.le.60000)   adir=trim(pathname_l2b)//'r055001_060000\r'
    !     if(iorbit.ge.60001 .and. iorbit.le.65000)   adir=trim(pathname_l2b)//'r060001_065000\r'
    !     if(iorbit.ge.65001 .and. iorbit.le.70000)   adir=trim(pathname_l2b)//'r065001_070000\r'
    !     if(iorbit.ge.70001 .and. iorbit.le.75000)   adir=trim(pathname_l2b)//'r070001_075000\r'
    !     if(iorbit.ge.75001 .and. iorbit.le.80000)   adir=trim(pathname_l2b)//'r075001_080000\r'
    !     if(iorbit.ge.80001 .and. iorbit.le.85000)   adir=trim(pathname_l2b)//'r080001_085000\r'
    !     if(iorbit.ge.85001 .and. iorbit.le.90000)   adir=trim(pathname_l2b)//'r085001_090000\r'
    !     if(iorbit.ge.90001 .and. iorbit.le.95000)   adir=trim(pathname_l2b)//'r090001_095000\r'
    !     if(iorbit.ge.95001 .and. iorbit.le.100000)  adir=trim(pathname_l2b)//'r095001_100000\r'
    !     if(iorbit.ge.100001 .and. iorbit.le.105000) adir=trim(pathname_l2b)//'r100001_105000\r'
    !     if(iorbit.ge.105001 .and. iorbit.le.110000) adir=trim(pathname_l2b)//'r150001_110000\r'
    !     if(iorbit.ge.110001 .and. iorbit.le.115000) adir=trim(pathname_l2b)//'r110001_115000\r'

    !     write(afile,9002) ksat,iorbit
    !     9002 format('f',i2.2,'_',i6.6,'.dat')
    !     filename_l1c=trim(adir)//trim(afile)

    !     return
    ! end subroutine get_ssmi_l1c_filename_base


    subroutine get_l2b_filename(pathname_l2b,ifinal,ksat,iorbit, filename_l2b,iexist)
    implicit none

    character(*) pathname_l2b,filename_l2b
    integer(4) ifinal,ksat,iorbit,iexist

    character(100) adir
    character(13) aorbit
      logical(4) lexist

    if(ksat.ne.32 .and. ksat.ne.33.and. ksat.ne.34) stop 'ksat oob in get_l2b_filename'

    if(iorbit.lt.1 .or. iorbit.gt. 99999) stop 'error1 in get_l2b_filename, pgm stopped'

      if(ifinal.eq.1) then
    if(iorbit.ge.    1 .and. iorbit.le. 5000)   adir=trim(pathname_l2b)//'r00001_05000\r'
    if(iorbit.ge. 5001 .and. iorbit.le.10000)   adir=trim(pathname_l2b)//'r05001_10000\r'
    if(iorbit.ge.10001 .and. iorbit.le.15000)   adir=trim(pathname_l2b)//'r10001_15000\r'
    if(iorbit.ge.15001 .and. iorbit.le.20000)   adir=trim(pathname_l2b)//'r15001_20000\r'
    if(iorbit.ge.20001 .and. iorbit.le.25000)   adir=trim(pathname_l2b)//'r20001_25000\r'
    if(iorbit.ge.25001 .and. iorbit.le.30000)   adir=trim(pathname_l2b)//'r25001_30000\r'
    if(iorbit.ge.30001 .and. iorbit.le.35000)   adir=trim(pathname_l2b)//'r30001_35000\r'
    if(iorbit.ge.35001 .and. iorbit.le.40000)   adir=trim(pathname_l2b)//'r35001_40000\r'
    if(iorbit.ge.40001 .and. iorbit.le.45000)   adir=trim(pathname_l2b)//'r40001_45000\r'
    if(iorbit.ge.45001 .and. iorbit.le.50000)   adir=trim(pathname_l2b)//'r45001_50000\r'
    if(iorbit.ge.50001 .and. iorbit.le.55000)   adir=trim(pathname_l2b)//'r50001_55000\r'
    if(iorbit.ge.55001 .and. iorbit.le.60000)   adir=trim(pathname_l2b)//'r55001_60000\r'
    if(iorbit.ge.60001 .and. iorbit.le.65000)   adir=trim(pathname_l2b)//'r60001_65000\r'
    if(iorbit.ge.65001 .and. iorbit.le.70000)   adir=trim(pathname_l2b)//'r65001_70000\r'
    if(iorbit.ge.70001 .and. iorbit.le.75000)   adir=trim(pathname_l2b)//'r70001_75000\r'
    if(iorbit.ge.75001 .and. iorbit.le.80000)   adir=trim(pathname_l2b)//'r75001_80000\r'
    if(iorbit.ge.80001 .and. iorbit.le.85000)   adir=trim(pathname_l2b)//'r80001_85000\r'
    if(iorbit.ge.85001 .and. iorbit.le.90000)   adir=trim(pathname_l2b)//'r85001_90000\r'
    if(iorbit.ge.90001 .and. iorbit.le.95000)   adir=trim(pathname_l2b)//'r90001_95000\r'
    if(iorbit.ge.95001 .and. iorbit.le.99999)   adir=trim(pathname_l2b)//'r95001_99999\r'

    else
    adir=trim(pathname_l2b)//'rt\'
    endif

    write(aorbit,9002) iorbit
 9002 format(i5.5,'.dat')

    filename_l2b=trim(adir)//trim(aorbit)

    inquire(file=filename_l2b,exist=lexist)
    if(lexist) then
    iexist=1
    else
    iexist=0
    endif

    return
    end subroutine get_l2b_filename

    subroutine get_l2b_filename_base(pathname_l2b,iorbit, filename_l2b)
    implicit none

    character(*) pathname_l2b,filename_l2b
    integer(4) iorbit

    character(100) adir
    character(13) aorbit

    if(iorbit.lt.1 .or. iorbit.gt. 99999) then
        filename_l2b=''
        return
    endif

  
    if(iorbit.ge.    1 .and. iorbit.le. 5000)   adir=trim(pathname_l2b)//'r00001_05000\r'
    if(iorbit.ge. 5001 .and. iorbit.le.10000)   adir=trim(pathname_l2b)//'r05001_10000\r'
    if(iorbit.ge.10001 .and. iorbit.le.15000)   adir=trim(pathname_l2b)//'r10001_15000\r'
    if(iorbit.ge.15001 .and. iorbit.le.20000)   adir=trim(pathname_l2b)//'r15001_20000\r'
    if(iorbit.ge.20001 .and. iorbit.le.25000)   adir=trim(pathname_l2b)//'r20001_25000\r'
    if(iorbit.ge.25001 .and. iorbit.le.30000)   adir=trim(pathname_l2b)//'r25001_30000\r'
    if(iorbit.ge.30001 .and. iorbit.le.35000)   adir=trim(pathname_l2b)//'r30001_35000\r'
    if(iorbit.ge.35001 .and. iorbit.le.40000)   adir=trim(pathname_l2b)//'r35001_40000\r'
    if(iorbit.ge.40001 .and. iorbit.le.45000)   adir=trim(pathname_l2b)//'r40001_45000\r'
    if(iorbit.ge.45001 .and. iorbit.le.50000)   adir=trim(pathname_l2b)//'r45001_50000\r'
    if(iorbit.ge.50001 .and. iorbit.le.55000)   adir=trim(pathname_l2b)//'r50001_55000\r'
    if(iorbit.ge.55001 .and. iorbit.le.60000)   adir=trim(pathname_l2b)//'r55001_60000\r'
    if(iorbit.ge.60001 .and. iorbit.le.65000)   adir=trim(pathname_l2b)//'r60001_65000\r'
    if(iorbit.ge.65001 .and. iorbit.le.70000)   adir=trim(pathname_l2b)//'r65001_70000\r'
    if(iorbit.ge.70001 .and. iorbit.le.75000)   adir=trim(pathname_l2b)//'r70001_75000\r'
    if(iorbit.ge.75001 .and. iorbit.le.80000)   adir=trim(pathname_l2b)//'r75001_80000\r'
    if(iorbit.ge.80001 .and. iorbit.le.85000)   adir=trim(pathname_l2b)//'r80001_85000\r'
    if(iorbit.ge.85001 .and. iorbit.le.90000)   adir=trim(pathname_l2b)//'r85001_90000\r'
    if(iorbit.ge.90001 .and. iorbit.le.95000)   adir=trim(pathname_l2b)//'r90001_95000\r'
    if(iorbit.ge.95001 .and. iorbit.le.99999)   adir=trim(pathname_l2b)//'r95001_99999\r'
    
    write(aorbit,9002) iorbit
 9002 format(i5.5)

    filename_l2b=trim(adir)//trim(aorbit)

    return
    end subroutine get_l2b_filename_base
    
    subroutine get_l2b_filename_local(pathname_l2b,ifinal,ksat,iorbit, filename_l2b,iexist)
    
    implicit none

    character(*) pathname_l2b,filename_l2b
    integer(4) ifinal,ksat,iorbit,iexist

    character(100) adir
    character(13) aorbit
    logical(4) lexist

    if(ksat.ne.32 .and. ksat.ne.33.and. ksat.ne.34) stop 'ksat oob in get_l2b_filename'
    if(iorbit.lt.1 .or. iorbit.gt. 99999) stop 'error1 in get_l2b_filename, pgm stopped'
    ifinal = 0
    adir=trim(pathname_l2b)//'r'
    write(aorbit,9002) iorbit
 9002 format(i5.5,'.dat')

    filename_l2b=trim(adir)//trim(aorbit)

    inquire(file=filename_l2b,exist=lexist)
    if(lexist) then
    iexist=1
    else
    iexist=0
    endif

    return
    end subroutine get_l2b_filename_local


subroutine get_l2b_filename_sim(pathname_l2b,month,iorbit,filename_l2b,iexist)
    implicit none

    character(*) pathname_l2b,filename_l2b
    integer(4) month,iorbit,iexist

    character(20) aorbit
    logical(4) lexist


    if(iorbit.lt.1 .or. iorbit.gt. 99999) stop 'error1 in get_l2b_filename, pgm stopped'


    write(aorbit,9002) month,iorbit
    
 9002 format('m',i2.2,'\r',i5.5,'.sdr.h5')

    filename_l2b=trim(pathname_l2b)//trim(aorbit)

    inquire(file=filename_l2b,exist=lexist)
    if(lexist) then
    iexist=1
    else
    iexist=0
    endif

    return
    end subroutine get_l2b_filename_sim
    
 
end module filename_routines
