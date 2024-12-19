	module ssmi_l1c

        use FileExist, only: file_exists
        use sat_id_module, only: maxscan,maxcel,maxchn,maxcel_85
        use open_file_routines, only: openbig
        implicit none

        !l1c file content
        integer(4) ksat_l1c,iorbit_l1c,numscan
        
        integer(4) iqual_flag(maxscan)
        
        real(8) scan_time(maxscan),orbit(maxscan),zang(maxscan),scpos(3,maxscan),scvel(3,maxscan)
        real(4) scloc(3,maxscan),therm(7,maxscan)

        real(4) cellat(maxcel_85,maxscan),cellon(maxcel_85,maxscan),celtht(maxcel_85,maxscan)
        real(4) celphi(maxcel_85,maxscan),celsun(maxcel_85,maxscan),celrfi(maxcel_85,maxscan)
        
        real(4) ta_nat(maxchn,maxcel_85,maxscan)
        real(4) ta_rsp(maxchn,maxcel_85,maxscan)
        real(4) tb_nat(maxchn,maxcel_85,maxscan)
        real(4) tb_rsp(maxchn,maxcel_85,maxscan)
	
    contains

        subroutine read_l1c_file(ksat,      &
                                 iorbit,    &
                                 filename_l1c,  &
                                 start_time,    &
                                 error)
            implicit none
            
            integer(4),intent(in)     :: ksat
            integer(4),intent(in)     :: iorbit
            character(120),intent(in) :: filename_l1c
            real(8),intent(out)       :: start_time
            integer(4),intent(out)    :: error

            !debug
            integer(4) :: scan
            !end debug

            print *, 'reading l1c file: ', filename_l1c
            
            if (file_exists(filename_l1c)) then
                call openbig(3,filename_l1c,'old')
                read(3) ksat_l1c,iorbit_l1c,numscan
                read(3) iqual_flag
                read(3) scan_time,orbit,zang,scpos,scvel,scloc,therm
                read(3) cellat,cellon,celtht,celphi,celsun,celrfi
                read(3) ta_nat,ta_rsp,tb_nat,tb_rsp
                close(3)
	
	            if(numscan.le.0 .or. numscan.gt.maxscan) then
                    stop 'error in ingesting numscan, pgm stopped'
                end if
	            if(ksat.ne.ksat_l1c .or. iorbit.ne.iorbit_l1c) then
                    stop 'error in ingesting ksat or iorbit, pgm stopped'
                end if
	            start_time=scan_time(1)

                !!! debug
                do scan = 1, maxscan, 100
                    print *, 'scan: ',scan,'scan_time: ', scan_time(scan)
                enddo
                !!! end debug
                error = 0
                return
            else
                error = 1
                return
            end if

            
	    end subroutine read_l1c_file
    end module ssmi_l1c

    module l2_module
        
        implicit none

        !    earth dimensions match wgs84    
        real(8),    parameter :: re=6378.137d3        !meters
        real(8),    parameter :: rp=6356.752d3        !meters
        real(8),    parameter :: ffac=(rp/re)*(rp/re)
        real(8),    parameter :: geosync_alt   =42164.0e3 

    end module l2_module 

    module l2_module_extra
        use sat_id_module
        use ssmi_l1c, only: maxscan,maxcel,maxchn,maxcel_85

        implicit none

        !    earth dimensions match wgs84    
        !real(8),    parameter :: re=6378.137d3        !meters
        !real(8),    parameter :: rp=6356.752d3        !meters
        !real(8),    parameter :: ffac=(rp/re)*(rp/re)
        !real(8),    parameter :: geosync_alt   =42164.0e3      !satellite orbits, montenbruck and gill, page 294
        !!      integer(4), parameter :: nsat_rfi=8
        !!      real(8),    parameter :: geosync_lon(nsat_rfi)=(/13., 19.2, 28.2, 317.0, 352.8, 10.0, 257.2, 260.8/)

        integer(4), parameter :: num_extra_scans = 3
        integer(4), parameter :: num_extra_cels = 3
        integer(4), parameter :: num_extra_fovs = 3
        integer(4), parameter :: maxscan_extra = maxscan*(1+num_extra_scans)
        integer(4), parameter :: maxcel_extra  = maxcel*(1+num_extra_cels) - num_extra_cels

        integer(4), parameter :: ilu_l2b_extra  =  7

        integer(4), parameter :: ncel_rsp = maxcel_extra  !number of cells for resampling program    

        integer(4), parameter ::  iradius_fieldavg=35  !see 'O:\ssmis\L2_processing\memo9.txt' for discussion of search radius

        character(120) filename_tb_native,filename_tb_resampled
        character(120) filename_ta_resample_wts
        character( 60) pathname_tb_native,pathname_tb_resampled

        real(8) ut1_1993_extra(maxscan_extra) !really tai_1993. it is a placeholder for heritage code
        integer(4) iscan_flag_extra(maxscan_extra)

        ! geolocation info -- need to keep track of these when geolocating the extra locations
        real(4) omega_extra(maxscan_extra),scan_index_extra(maxscan_extra)
        real(8) scpos_eci0_extra(3,maxscan_extra), x_eci0_extra(3,maxscan_extra), y_eci0_extra(3,maxscan_extra), z_eci0_extra(3,maxscan_extra)
        real(8) zang_extra(maxscan_extra),scan_time_extra(maxscan_extra)

        !not sure if we need to reconsider all the sun and moon things....
        integer(1) kflag_sun_extra(maxscan_extra)
        real(8) sunvec0_extra(3,maxscan_extra),moonvec0_extra(3,maxscan_extra),alpha_sun_extra(maxscan_extra),beta_sun_extra(maxscan_extra)
        real(4) sunlat_extra(maxscan_extra),sunlon_extra(maxscan_extra),sundis_extra(maxscan_extra)

        real(4) scloc_extra(3,maxscan_extra)
        real(4) cellat_extra(maxcel_extra,maxscan_extra),cellon_extra(maxcel_extra,maxscan_extra)
        real(4) celrange_extra(maxcel_extra,maxscan_extra)
        real(4) celphi_extra(maxcel_extra,maxscan_extra),celtht_extra(maxcel_extra,maxscan_extra)
        real(4) celsun_extra(maxcel_extra,maxscan_extra)

        !we'll use the rfi flag from the tb file
        real(4) reflat_extra(maxcel_extra,maxscan_extra),reflon_extra(maxcel_extra,maxscan_extra) !geostationary lat and lon of reflected ray

        !not doing 89 right now....
        !real(4) cellat_89(maxcel_89,maxscan_89),cellon_89(maxcel_89,maxscan_89),cellat_cold(maxscan),cellon_cold(maxscan)
        !real(4) ta_cold(nch,maxscan),ta_hot(nch,maxscan),ta_spill(nch,maxscan),gainest(nch,maxscan)
        

        ! this is where the resampled TBs go
        real(4) tbr_extra(maxcel_extra,maxscan_extra,maxchn) !,ta_89(maxcel_89,maxscan_89,2)

        !   ocean  block
        integer(1) isurcel_extra(5,maxcel_extra,maxscan_extra)
        integer(1) ice_flag_extra(3,maxcel_extra,maxscan_extra)
        integer(1) ice_flag2_extra(maxcel_extra,maxscan_extra)
        integer(1) ioob_flag_extra(5,maxcel_extra,maxscan_extra)

        real(4) ta_extra(maxcel_extra,maxscan_extra,maxchn)

    end module l2_module_extra