 module ssmi_filename_routines
 
 contains
    subroutine get_ssmi_orbit_subdir(iorbit, subdir)
        
        implicit none
        integer(4),intent(in)        :: iorbit
        character(14),intent(out)    :: subdir

        integer(4) :: lower

        if(iorbit.lt.1 .or. iorbit.gt. 995000) then
            stop 'orbit out of bounds in get_ssmi_l1c_filename_base'
        endif

        lower = ((iorbit-1)/5000)*5000+1
        write(subdir,9001) lower,lower+4999
        9001 format('r',i6.6,'_',i6.6)

        ! if(iorbit.ge.    1 .and. iorbit.le. 5000)   subdir ='r000001_005000'
        ! if(iorbit.ge. 5001 .and. iorbit.le.10000)   subdir ='r005001_010000'
        ! if(iorbit.ge.10001 .and. iorbit.le.15000)   subdir ='r010001_015000'
        ! if(iorbit.ge.15001 .and. iorbit.le.20000)   subdir ='r015001_020000'
        ! if(iorbit.ge.20001 .and. iorbit.le.25000)   subdir ='r020001_025000'
        ! if(iorbit.ge.25001 .and. iorbit.le.30000)   subdir ='r025001_030000'
        ! if(iorbit.ge.30001 .and. iorbit.le.35000)   subdir ='r030001_035000'
        ! if(iorbit.ge.35001 .and. iorbit.le.40000)   subdir ='r035001_040000'
        ! if(iorbit.ge.40001 .and. iorbit.le.45000)   subdir ='r040001_045000'
        ! if(iorbit.ge.45001 .and. iorbit.le.50000)   subdir ='r045001_050000'
        ! if(iorbit.ge.50001 .and. iorbit.le.55000)   subdir ='r050001_055000'
        ! if(iorbit.ge.55001 .and. iorbit.le.60000)   subdir ='r055001_060000'
        ! if(iorbit.ge.60001 .and. iorbit.le.65000)   subdir ='r060001_065000'
        ! if(iorbit.ge.65001 .and. iorbit.le.70000)   subdir ='r065001_070000'
        ! if(iorbit.ge.70001 .and. iorbit.le.75000)   subdir ='r070001_075000'
        ! if(iorbit.ge.75001 .and. iorbit.le.80000)   subdir ='r075001_080000'
        ! if(iorbit.ge.80001 .and. iorbit.le.85000)   subdir ='r080001_085000'
        ! if(iorbit.ge.85001 .and. iorbit.le.90000)   subdir ='r085001_090000'
        ! if(iorbit.ge.90001 .and. iorbit.le.95000)   subdir ='r090001_095000'
        ! if(iorbit.ge.95001 .and. iorbit.le.100000)  subdir ='r095001_100000'
        ! if(iorbit.ge.100001 .and. iorbit.le.105000) subdir ='r100001_105000'
        ! if(iorbit.ge.105001 .and. iorbit.le.110000) subdir ='r150001_110000'
        ! if(iorbit.ge.110001 .and. iorbit.le.115000) subdir ='r110001_115000'
        ! if(iorbit.ge.115001 .and. iorbit.le.120000) subdir ='r115001_120000'
        ! if(iorbit.ge.120001 .and. iorbit.le.125000) subdir ='r120001_125000'
        ! if(iorbit.ge.125001 .and. iorbit.le.130000) subdir ='r125001_130000'
        ! if(iorbit.ge.130001 .and. iorbit.le.135000) subdir ='r130001_135000'
        ! if(iorbit.ge.135001 .and. iorbit.le.140000) subdir ='r135001_140000'
        ! if(iorbit.ge.140001 .and. iorbit.le.145000) subdir ='r140001_145000'
        ! if(iorbit.ge.145001 .and. iorbit.le.150000) subdir ='r145001_150000'
        ! if(iorbit.ge.150001 .and. iorbit.le.155000) subdir ='r150001_155000'
        ! if(iorbit.ge.155001 .and. iorbit.le.160000) subdir ='r155001_160000'
        ! if(iorbit.ge.160001 .and. iorbit.le.165000) subdir ='r160001_165000'
        ! if(iorbit.ge.165001 .and. iorbit.le.170000) subdir ='r165001_170000'
        
    end subroutine get_ssmi_orbit_subdir

    subroutine get_ssmi_l1c_filename(pathname_l1c,ksat,iorbit, filename_l1c)
        implicit none

        character(*),intent(in)     :: pathname_l1c
        integer(4),intent(in)       :: ksat
        integer(4),intent(in)       :: iorbit
        character(*),intent(out)    :: filename_l1c

        character(14)               :: subdir
        character(100)              :: adir
        character(100)              :: afile

        call get_ssmi_orbit_subdir(iorbit,subdir)
        write(afile,9002) subdir,ksat,iorbit
        9002 format(a,'\f',i2.2,'_r',i6.6,'.dat')
        filename_l1c=trim(pathname_l1c)//trim(afile)

        return
    end subroutine get_ssmi_l1c_filename

    subroutine get_ssmi_l2b_filename_base(pathname_l2b,ksat,iorbit,region_str,target_size,filename_l2b)
        implicit none

        character(*),intent(in)  :: pathname_l2b
        integer(4),intent(in)    :: ksat
        integer(4),intent(in)    :: iorbit
        character(2),intent(in)  :: region_str ! GL, NP, SP
        integer(4),intent(in)    :: target_size
        character(*),intent(out) :: filename_l2b

        character(14)  :: subdir
        character(100) :: afile

        call get_ssmi_orbit_subdir(iorbit,subdir)
        print *,subdir

        write(afile,9002) ksat,region_str,target_size,subdir,ksat,iorbit
        print *,afile
        9002 format('f',i2.2,'\',a,'_',i2.2,'\',a,'\f',i2.2,'_',i6.6)
        filename_l2b=trim(pathname_l2b)//trim(afile)
        print *,filename_l2b
        print *,trim(pathname_l2b)  
    return
    end subroutine get_ssmi_l2b_filename_base
end module ssmi_filename_routines
