!   7/13/2021 Converted to module for ACCESS project

!     7/15/2014 version change 7/17/2014.  subroutine get_amsr2_l2b_filename_rt added
  
!     12/19/2013 version changed 7/15/2014. l2b folder is now 7.2 rather than 7.1
!     also the name of this file is changed from filename_routines to io_filenames

!     3/7/2013 version change 12/19/2013  l2b folder name is now V07.1 rather than V07a

module io_filenames
    use FileExist, only: file_exist
    private
    public get_amsr2_rss_l1a_filename,get_amsr2_l2b_filename,get_amsr2_l2b_filename_rt

contains
    subroutine get_amsr2_rss_l1a_filename(iorbit, filename_l1a,iexist) 
        implicit none

        character*(*) filename_l1a
        character(80) adir
        character(8) aorbit
        logical(4) lexist
        integer(4) iorbit,iexist

        if(iorbit.lt.1 .or. iorbit.gt.99999) stop 'orbit oob in get_get_rss_l1a_filename, pgm stopped'

        if(iorbit.ge.    1 .and. iorbit.le. 5000) adir='j:\amsr2\l1a\r00001_05000\r'
        if(iorbit.ge. 5001 .and. iorbit.le.10000) adir='j:\amsr2\l1a\r05001_10000\r'
        if(iorbit.ge.10001 .and. iorbit.le.15000) adir='j:\amsr2\l1a\r10001_15000\r'
        if(iorbit.ge.15001 .and. iorbit.le.20000) adir='j:\amsr2\l1a\r15001_20000\r'
        if(iorbit.ge.20001 .and. iorbit.le.25000) adir='j:\amsr2\l1a\r20001_25000\r'
        if(iorbit.ge.25001 .and. iorbit.le.30000) adir='j:\amsr2\l1a\r25001_30000\r'
        if(iorbit.ge.30001 .and. iorbit.le.35000) adir='j:\amsr2\l1a\r30001_35000\r'
        if(iorbit.ge.35001 .and. iorbit.le.40000) adir='j:\amsr2\l1a\r35001_40000\r'
        if(iorbit.ge.40001 .and. iorbit.le.45000) adir='j:\amsr2\l1a\r40001_45000\r'
        if(iorbit.ge.45001 .and. iorbit.le.50000) adir='j:\amsr2\l1a\r45001_50000\r'
        if(iorbit.ge.50001 .and. iorbit.le.55000) adir='j:\amsr2\l1a\r50001_55000\r'
        if(iorbit.ge.55001 .and. iorbit.le.60000) adir='j:\amsr2\l1a\r55001_60000\r'
        if(iorbit.ge.60001 .and. iorbit.le.65000) adir='j:\amsr2\l1a\r60001_65000\r'
        if(iorbit.ge.65001 .and. iorbit.le.70000) adir='j:\amsr2\l1a\r65001_70000\r'
        if(iorbit.ge.70001 .and. iorbit.le.75000) adir='j:\amsr2\l1a\r70001_75000\r'
        if(iorbit.ge.75001 .and. iorbit.le.80000) adir='j:\amsr2\l1a\r75001_80000\r'
        if(iorbit.ge.80001 .and. iorbit.le.85000) adir='j:\amsr2\l1a\r80001_85000\r'
        if(iorbit.ge.85001 .and. iorbit.le.90000) adir='j:\amsr2\l1a\r85001_90000\r'
        if(iorbit.ge.90001 .and. iorbit.le.95000) adir='j:\amsr2\l1a\r90001_95000\r'
        if(iorbit.ge.95001 .and. iorbit.le.99999) adir='j:\amsr2\l1a\r95001_100000\r'

        write(aorbit,9001) iorbit
        9001 format(i5.5,'.gz') 
        filename_l1a=trim(adir)//aorbit

        call file_exist(filename_l1a, lexist)
        if(lexist) then
        iexist=1
        else
        iexist=0
        endif

        return
    end subroutine get_amsr2_rss_l1a_filename

    subroutine get_amsr2_l2b_filename(iorbit, filename_l2b,iexist) 
        implicit none

        character*(*) filename_l2b
        character(80) adir
        character(9) aorbit
        logical(4) lexist
        integer(4) iorbit,iexist
    
        if(iorbit.lt.1 .or. iorbit.gt.99999) stop 'orbit oob in get_get_rss_l1a_filename, pgm stopped'

        if(iorbit.ge.    1 .and. iorbit.le. 5000) adir='j:\amsr2\l2b_v08\r00001_05000\r'
        if(iorbit.ge. 5001 .and. iorbit.le.10000) adir='j:\amsr2\l2b_v08\r05001_10000\r'
        if(iorbit.ge.10001 .and. iorbit.le.15000) adir='j:\amsr2\l2b_v08\r10001_15000\r'
        if(iorbit.ge.15001 .and. iorbit.le.20000) adir='j:\amsr2\l2b_v08\r15001_20000\r'
        if(iorbit.ge.20001 .and. iorbit.le.25000) adir='j:\amsr2\l2b_v08\r20001_25000\r'
        if(iorbit.ge.25001 .and. iorbit.le.30000) adir='j:\amsr2\l2b_v08\r25001_30000\r'
        if(iorbit.ge.30001 .and. iorbit.le.35000) adir='j:\amsr2\l2b_v08\r30001_35000\r'
        if(iorbit.ge.35001 .and. iorbit.le.40000) adir='j:\amsr2\l2b_v08\r35001_40000\r'
        if(iorbit.ge.40001 .and. iorbit.le.45000) adir='j:\amsr2\l2b_v08\r40001_45000\r'
        if(iorbit.ge.45001 .and. iorbit.le.50000) adir='j:\amsr2\l2b_v08\r45001_50000\r'
        if(iorbit.ge.50001 .and. iorbit.le.55000) adir='j:\amsr2\l2b_v08\r50001_55000\r'
        if(iorbit.ge.55001 .and. iorbit.le.60000) adir='j:\amsr2\l2b_v08\r55001_60000\r'
        if(iorbit.ge.60001 .and. iorbit.le.65000) adir='j:\amsr2\l2b_v08\r60001_65000\r'
        if(iorbit.ge.65001 .and. iorbit.le.70000) adir='j:\amsr2\l2b_v08\r65001_70000\r'
        if(iorbit.ge.70001 .and. iorbit.le.75000) adir='j:\amsr2\l2b_v08\r70001_75000\r'
        if(iorbit.ge.75001 .and. iorbit.le.80000) adir='j:\amsr2\l2b_v08\r75001_80000\r'
        if(iorbit.ge.80001 .and. iorbit.le.85000) adir='j:\amsr2\l2b_v08\r80001_85000\r'
        if(iorbit.ge.85001 .and. iorbit.le.90000) adir='j:\amsr2\l2b_v08\r85001_90000\r'
        if(iorbit.ge.90001 .and. iorbit.le.95000) adir='j:\amsr2\l2b_v08\r90001_95000\r'
        if(iorbit.ge.95001 .and. iorbit.le.99999) adir='j:\amsr2\l2b_v08\r95001_99999\r'

        write(aorbit,9001) iorbit
        9001 format(i5.5,'.dat') 
        filename_l2b=trim(adir)//aorbit

        call file_exist(filename_l2b, lexist)
        if(lexist) then
            iexist=1
        else
             iexist=0
        endif

        return
    end subroutine get_amsr2_l2b_filename


    subroutine get_amsr2_l2b_filename_rt(iorbit, filename_l2b,iexist) 
        implicit none

        character*(*) filename_l2b
        character(80) adir
        character(9) aorbit
        logical(4) lexist
        integer(4) iorbit,iexist
    
        if(iorbit.lt.1 .or. iorbit.gt.99999) stop 'orbit oob in get_get_rss_l1a_filename, pgm stopped'

        adir='j:\amsr2\l2b_v08\rt\r'

        write(aorbit,9001) iorbit
        9001 format(i5.5,'.dat') 
        filename_l2b=trim(adir)//aorbit

        call file_exist(filename_l2b, lexist)
        if(lexist) then
            iexist=1
        else
            iexist=0
        endif

        return
    end subroutine get_amsr2_l2b_filename_rt

end module io_filenames
