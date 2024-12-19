!     feb 5, 2013 versin change 8/29/2013,  the indexing for 89b was changed to make it
!     adjacent to 89a. 'O:\amsr2\resampling\memo1.txt'. the 89 remapping routine has not yet 
!     been used for amsr2, so there is no impact.

!     april 25 2007 was changed on feb 29, 2008.  there was bug in resampling 89h.
!     the wrong ipol index was used (1 rather than 2) see comment cbug

!     the 12/10/2005 version was changed on april 25 2007.  the 89 ghz v vs. h renormalization coefs
!     are now different for amsr-e and amsr-j.  this has no effect on amsr-e

!     the 8/26/2005 version was changed on 12/10/2005.  actual values for pv0,pv1,ph0,ph1
!     were inserted. only amsr_sips_l1a_processing_v05 uses these values and it is yet to be run.  
module resample

    use, intrinsic :: iso_fortran_env, only: int8, int32, real32, real64, ERROR_UNIT
    use ieee_arithmetic, only: ieee_is_finite
    use NAN_support, only: nan_f32

    use map_type_constants
    use equirectangular_maps
    use equirectangular_maps64

    !use sat_id_module
    use sat_id_module, only: maxscan,maxcel,maxchn,maxcel_85
    use ssmi_l1c, only: tb_nat,iqual_flag
    use l2_module, only: re,rp,ffac
    use l2_module_extra, only: maxcel_extra, maxscan_extra, cellat_extra, cellon_extra, tbr_extra, scan_time_extra
    use l2_module_extra, only: num_extra_cels, num_extra_scans
    
    implicit none 

    type ResampleWeights
        integer(int32) :: num_input_cells  ! = maxcel 
        integer(int32) :: num_output_cells ! = (1+num_extra_cels)*maxcel - num_extra_cels 
        integer(int32) :: num_input_scans  ! = 29
        integer(int32) :: num_output_scans ! = (1+num_extra_scans)
        real(real64),allocatable :: wts(:,:,:,:)
    end type ResampleWeights


    integer(int32),parameter :: max_select=5000

    private
    public init_resample_weights,free_resample_weights
    public resample_tb_extra,resample_tb_extra_85
    public ResampleWeights
    public find_closest_resample_point
    public resample_onto_lat_lon_grid
    public resample_time_onto_lat_lon_grid
    public resample_onto_lat_lon_grid_interp
    public write_resampled_time_map
    public write_resampled_maps
    public write_wt_slice_netcdf

contains

    subroutine allocate_resample_weights(weights)

        use sat_id_module
        use l2_module
        use l2_module_extra  ! contains info about extra locations

        type(ResampleWeights),intent(inout) :: weights

        print *,'Allocating resample weights'
        print *,'num_input_cells  = ',weights%num_input_cells
        print *,'num_output_cells = ',weights%num_output_cells
        print *,'num_input_scans  = ',weights%num_input_scans
        print *,'num_output_scans = ',weights%num_output_scans

        allocate(weights%wts(weights%num_input_cells,weights%num_input_scans, &
                             weights%num_output_cells,weights%num_output_scans))

    end subroutine allocate_resample_weights

    subroutine deallocate_resample_weights(weights)

        use sat_id_module
        use l2_module
        use l2_module_extra  ! contains info about extra locations

        type(ResampleWeights),intent(inout) :: weights

        deallocate(weights%wts)

    end subroutine deallocate_resample_weights

    subroutine init_resample_weights(sat_name,bnd_number,target_size,weights)

        character(*),intent(in) :: sat_name
        integer(4),intent(in)   :: bnd_number
        integer(4),intent(in)   :: target_size

        integer(4) :: ksat

        type(ResampleWeights),intent(inout)  :: weights

        character(100),parameter :: wt_path = 'L:\access\resampling\'
        character(150)           :: weight_file_name

        weights%num_input_cells  = maxcel 
        weights%num_output_cells = (1+num_extra_cels)*maxcel - num_extra_cels 
        weights%num_input_scans  = 29
        weights%num_output_scans = (1+num_extra_scans)

        if ((sat_name == 'SSMI') .and. (bnd_number == 3)) then
            weights%num_input_cells  = maxcel_85
            weights%num_input_scans  = 57
        endif

        call allocate_resample_weights(weights)

        select case (sat_name)
            case('AMSR2')
                write(weight_file_name,100)trim(wt_path),trim(sat_name),target_size,bnd_number-1  !band number in files names 1 less
                100 format(a,a,'\resample_weights\circular_',i2.2,'km\resample_weights_band_',i2.2,'.dat')
            case('SSMI')
                ksat = 13
                write(weight_file_name,101)trim(wt_path),trim(sat_name),ksat,target_size,bnd_number  !band number in files names 1 less
                101 format(a,a,'\f',i2.2,'\resample_weights\'i2.2,'km\resample_weights_band_',i2.2,'.dat')
            case default
                print *,'Error in init_resample_weights: unknown sat_name = ',sat_name
                error stop
        end select

        print *,bnd_number
        print *,weight_file_name
        open(3,file=weight_file_name,status='old',form='binary')
        read(3) weights%wts
        close(3)
    
    end subroutine init_resample_weights

    subroutine free_resample_weights(weights)
        use sat_id_module
        use l2_module
        use l2_module_extra

        type(ResampleWeights),intent(inout)  :: weights

        character(100),parameter :: wt_path = 'L:\access\resampling\'
        character(150)           :: weight_file_name

        weights%num_input_cells  = 0
        weights%num_output_cells = 0
        weights%num_input_scans  = 0
        weights%num_output_scans = 0

        call deallocate_resample_weights(weights)
    end subroutine free_resample_weights

    subroutine resample_tb_extra(weights,ichan)

        use l2_module 
        use l2_module_extra

        type(ResampleWeights),intent(in)  :: weights
        integer(int32)                    :: ichan

        real(real64),allocatable          :: wt_slice(:,:)
        ! the current slice of the 4D weight array

        real(real32),allocatable         :: tb_nat_odd(:,:,:)
        
        
        integer(int32) :: iscan,kscan,kscan_start,kscan_end,kfov,kchn
        integer(int32) :: scan_number_out
        real(real64)   :: scan_index_out
        integer(int32) :: jfov,jscan
        real(real64)   :: wt_sum,tb_sum
        real(real64)   :: wt_sum2,tb_sum2
        
        integer(int32) :: num_cells_included

        integer(int32) :: jscan_list(max_select)
        integer(int32) :: jfov_list(max_select)
        real(real64)   :: wt_list(max_select),jwt
        integer(int32) :: num_select,ipos
        
        print *,'ichan = ',ichan
        tbr_extra(:,:,ichan) = -999.0

        kscan_start = 0
        kscan_end = 0
        if (num_extra_scans > 0) then
            kscan_start = -1*floor(num_extra_scans/2.0)
            kscan_end = floor((1+num_extra_scans)/2.0)
        endif
        !print *,'kscan_start,kscan_end = ',kscan_start,kscan_end
        allocate(wt_slice(weights%num_input_cells,weights%num_input_scans))
        allocate(tb_nat_odd(maxchn,maxcel,maxscan/2))

        !For 19, 22, and 37 GHz, tb_nat_odd has the actual observations

        do kscan = 1,maxscan/2
            if (iqual_flag((kscan*2)-1) .ne. 0) then
                tb_nat_odd(:,:,kscan) = nan_f32
            else
                do kfov = 1,maxcel
                    tb_nat_odd(:,kfov,kscan) = tb_nat(:,(kfov*2)-1,(kscan*2)-1)
                    do kchn = 1,6
                        if ((tb_nat_odd(ichan,kfov,kscan)) < 20.0) then
                            tb_nat_odd(ichan,kfov,kscan) = nan_f32
                        endif
                    enddo
                enddo
            endif
        enddo
            
        !The out loops are the position in the output scan....
        do kscan = kscan_start,kscan_end
            do kfov = 1,weights%num_output_cells
                wt_slice = weights%wts(:,:,kfov,kscan-kscan_start+1)
                !Make a list of the non-zero elements
                num_select = 0
                jscan_list = 0
                jfov_list = 0
                wt_list = 0.0D0
                do jscan = -14,14
                    do jfov = 1,maxcel
                        if (abs(wt_slice(jfov,jscan+15)) > 0.00000001) then
                            num_select= num_select+1
                            jscan_list(num_select) = jscan
                            jfov_list(num_select) = jfov
                            wt_list(num_select) = wt_slice(jfov,jscan+15)
                        endif
                    enddo
                enddo

                !for now, loop over all the positions
                do iscan = 14,(maxscan/2)-14 ! loop over the entire orbit, except the ends
                    scan_number_out = ((iscan-1)*(1+num_extra_scans))+kscan-kscan_start
                    
                    ! if ((iscan==14) .and. (kfov == 1)) then
                    !     print *,'iscan,kscan,scan_number_out = ',iscan,kscan,scan_number_out
                    !     print *
                    ! endif
                    wt_sum = 0.0D0
                    tb_sum = 0.0D0
                    num_cells_included = 0

                    !This is a loop over the non-zero elements of the wt array
                    do ipos = 1,num_select
                        jscan = jscan_list(ipos)
                        jfov = jfov_list(ipos)
                        jwt = wt_list(ipos)
                        num_cells_included = num_cells_included + 1

                        if (tb_nat_odd(ichan,jfov,iscan+jscan) > 0.0) then
                            wt_sum = wt_sum + jwt
                            tb_sum = tb_sum + jwt*tb_nat_odd(ichan,jfov,iscan+jscan)
                        endif
                        !print *,num_cells_included,jwt,wt_sum,tb_sum
                    enddo
                    ! print *,'scan_number_out,kfov,tb_sum = ',scan_number_out,kfov,tb_sum
                    if (wt_sum > 0.75) then
                        tbr_extra(kfov,scan_number_out,ichan) = tb_sum/wt_sum
                    else
                        tbr_extra(kfov,scan_number_out,ichan) = nan_f32
                    endif
                enddo !iscan
            enddo !kfov
        enddo !kscan
        deallocate(wt_slice)
        deallocate(tb_nat_odd)
    end subroutine resample_tb_extra

    subroutine resample_tb_extra_85(weights,ichan)

        use l2_module 
        use l2_module_extra

        type(ResampleWeights),intent(in)  :: weights
        integer(int32)                    :: ichan

        real(real64),allocatable          :: wt_slice(:,:)
        ! the current slice of the 4D weight array

        integer(int32) :: iscan,kscan,kscan_start,kscan_end,kfov,kchn
        integer(int32) :: scan_number_out
        real(real64)   :: scan_index_out
        integer(int32) :: jfov,jscan
        real(real64)   :: wt_sum,tb_sum
        real(real64)   :: wt_sum2,tb_sum2
        
        integer(int32) :: num_cells_included

        integer(int32) :: jscan_list(max_select)
        integer(int32) :: jfov_list(max_select)
        real(real64)   :: wt_list(max_select),jwt
        integer(int32) :: num_select,ipos
        
        print *,'ichan = ',ichan
        tbr_extra(:,:,ichan) = -999.0

        kscan_start = 0
        kscan_end = 0
        if (num_extra_scans > 0) then
            kscan_start = -1*floor(num_extra_scans/2.0)
            kscan_end = floor((1+num_extra_scans)/2.0)
        endif
        ! print *,'number of extra scans = ',num_extra_scans
        ! print *,'kscan_start,kscan_end = ',kscan_start,kscan_end
        ! print *,'num_input_cells=',weights%num_input_cells
        ! print *,'num_input_scans=',weights%num_input_scans,weights%num_input_scans/2
        allocate(wt_slice(weights%num_input_cells,weights%num_input_scans))
        ! print *,size(weights%wts,dim = 1),size(weights%wts,dim = 2),size(weights%wts,dim = 3),size(weights%wts,dim = 4)
    
        !The outer loops are the position in the output scan so that the weight
        !arrays are accessed/sliced less often
        do kscan = kscan_start,kscan_end
            do kfov = 1,weights%num_output_cells
                wt_slice = weights%wts(:,:,kfov,kscan-kscan_start+1)
                !Make a list of the non-zero elements
                num_select = 0
                jscan_list = 0
                jfov_list = 0
                wt_list = 0.0D0
                
                do jscan = -weights%num_input_scans/2,weights%num_input_scans/2
                    do jfov = 1,weights%num_input_cells
                        !print *,jfov,jscan,jscan+1+weights%num_input_scans/2
                        !print *,wt_slice(jfov,jscan+1+weights%num_input_scans/2)
                        if (abs(wt_slice(jfov,jscan+1+weights%num_input_scans/2)) > 0.00000001) then
                            num_select= num_select+1
                            jscan_list(num_select) = jscan
                            jfov_list(num_select) = jfov
                            wt_list(num_select) = wt_slice(jfov,jscan+1+weights%num_input_scans/2)
                        endif
                    enddo
                enddo
                !loop over all the positions
                !print *,weights%num_input_scans/2,maxscan-weights%num_input_scans/2

                do iscan = weights%num_input_scans/2,maxscan-(weights%num_input_scans/2) ! loop over the entire orbit, except the ends
                    scan_number_out = ((iscan-1)*(1+num_extra_scans)/2)+kscan-kscan_start
                    ! print *,iscan,scan_number_out
                    ! if ((iscan==28) .and. (kfov == 1)) then
                    !      print *,'iscan,kscan,scan_number_out = ',iscan,kscan,scan_number_out
                    !      print *
                    ! endif
                    wt_sum = 0.0D0
                    tb_sum = 0.0D0
                    num_cells_included = 0

                    !This is a loop over the non-zero elements of the wt array
                    do ipos = 1,num_select
                        jscan = jscan_list(ipos)
                        jfov = jfov_list(ipos)
                        jwt = wt_list(ipos)
                        num_cells_included = num_cells_included + 1
                        
                        
                        ! print *,iscan,jscan,jfov,jwt
                        ! print *,ichan,jfov,iscan+jscan
                        ! print *,size(tb_nat,dim = 1),size(tb_nat,dim = 2),size(tb_nat,dim = 3)
                        ! print *,tb_nat(ichan,jfov,iscan+jscan)
                        ! if ((iscan+jscan) .gt. 1800) then
                        !     print *,iscan,jscan 
                        ! endif
                        if (tb_nat(ichan,jfov,iscan+jscan) > 0.0) then
                            tb_sum = tb_sum + jwt*tb_nat(ichan,jfov,iscan+jscan)
                            wt_sum = wt_sum + jwt
                        endif
                        
                    enddo
                    !print *,'scan_number_out,kfov,tb_sum = ',scan_number_out,kfov,tb_sum
                    if (wt_sum > 0.75) then
                        tbr_extra(kfov,scan_number_out,ichan) = tb_sum/wt_sum
                    else
                        tbr_extra(kfov,scan_number_out,ichan) = nan_f32
                        !print *,scan_number_out,kfov,num_cells_included,wt_sum,tb_sum
                    endif
                    !print *
                enddo !iscan
            enddo !kfov
        enddo !kscan
        deallocate(wt_slice)
    end subroutine resample_tb_extra_85

    subroutine find_closest_resample_point(grid_lat, &
                                           grid_lon, &
                                           search_size, &
                                           channel, &
                                           scan_closest, &
                                           fov_closest, &
                                           closest_distance)

        !use l2_module, only: re
        use l2_module_extra
        use wgs84, only: distance_between_lon_lats
        real(real32),intent(in)    :: grid_lat,grid_lon,search_size
        integer(int32),intent(in)  :: channel
        integer(int32),intent(out) :: scan_closest,fov_closest
        real(real32),intent(out)   :: closest_distance

        integer(int32)  :: jscan,jfov,error
        real(real32)    :: distance

        scan_closest = -999
        fov_closest = -999
        closest_distance = 999.0

        do jscan = 1,maxscan_extra
            if (ieee_is_finite(scan_time_extra(jscan))) then
                do jfov = 1,maxcel_extra
                    if (tbr_extra(jfov,jscan,channel) > 1.0) then
                        ! There is a good tb, worth evaluating distance
                        if (abs(cellat_extra(jfov,jscan) - grid_lat) < search_size) then
                            if (abs(cellon_extra(jfov,jscan) - grid_lon) < search_size) then
                                !calculate distance
                                distance = distance_between_lon_lats(grid_lon,grid_lat, &
                                                                    cellon_extra(jfov,jscan),&
                                                                    cellat_extra(jfov,jscan),&
                                                                    error)
                                
                                if (error == 0) then
                                    if (distance < closest_distance) then
                                        closest_distance = distance
                                        scan_closest = jscan
                                        fov_closest = jfov 
                                    endif
                                else
                                    print *,'Error in distance_between_lon_lats'
                                    error stop
                                endif
                            endif
                        endif
                    endif
                enddo !jfov
            endif
        enddo !jscan
    end subroutine find_closest_resample_point

    subroutine resample_onto_lat_lon_grid_interp(resolution,  &
                                                channel,  &
                                                do_time,  &
                                                gridded_time, &
                                                gridded_tbs, &
                                                min_valid_extra_scan, &
                                                max_valid_extra_scan, &
                                                num_overlap_scans)

        ! This version of the resampling routine interpolates to the desired grid point using a general 
        ! quadrilateral interpolation method.
        ! The looping is a little different than the previous versions
        !
        ! Here, the outside loops are over scans/fovs
        ! Then, the center point of the quadrilateral formed by (iscan,ifov),(iscan,ifov+1),(iscan+1,ifov+1),(iscan+1,ifov)
        ! is found.
        ! Then, the closest target point to this point is found and checked to see if it inside the quadrilateral
        ! If it is, then the interpolation to this point is done.
        !
        ! The algorithm also searches nearby fixed grid points to see if they are inside the quadrilateral
        ! The search range in the E-W direction expands near the pole so all possible grid points should be
        ! filled when using our standard 0.25 degree grids.

        !use l2_module, only: re  !radius of Earth, not regular expressions (!)
        use l2_module_extra
        use QuadrilateralUtilties, only: check_convex, &
                                         is_inside_quadrilateral, &
                                         order_vertices_quadrilateral

        use InterpBilinearQuadrilateral, only: init_InterpBilinearQuadrilateral,eval_InterpBilinearQuadrilateral
        use InterpBilinearQuadrilateral, only: quad_interpolator_type

        implicit none

        real(real32),intent(in)     :: resolution
        integer(int32),intent(in)   :: channel
        logical,intent(in)          :: do_time
        type(map_str64),intent(inout) :: gridded_time
        type(map_str),intent(inout) :: gridded_tbs
        integer(int32),intent(in)   :: min_valid_extra_scan
        integer(int32),intent(in)   :: max_valid_extra_scan
        integer(int32),intent(in)   :: num_overlap_scans

        !locals
        integer(int32)              :: ilat,ilon,closest_scan,closest_fov
        real(real32)                :: xlat,xlon,closest_distance
        integer(int32)              :: num_grid_points_filled
        integer(int32)              :: error
        integer(int32)              :: max_valid_scan

        integer(int32)              :: jscan,jfov,kscan,kfov,jlat,jlon,klat,klon,klat0,klon0,ilatx,ilonx,iquad,j,k,dk
        integer(int32)              :: ilon_closest,ilat_closest
        real(real32)                :: distance
        logical                     :: bad_loc,is_convex,is_inside
        integer(int32)              :: interp_error
        integer(int32)              :: i
        type(quad_interpolator_type):: qi

        real(real64)                :: tot_lat,tot_lon,first_lat,first_lon
        real(real32)                :: grid_lat,grid_lon
        real(real32),dimension(0:3) :: quad_lats,quad_lons,quad_ta
        real(real64),dimension(0:3) :: quad_times
        real(real64),dimension(0:3) :: quad_lats_sorted,quad_lons_sorted,quad_ta_sorted,quad_times_sorted
        real(real64)                :: mean_lat,mean_lon,tb_test
        integer(int32),dimension(0:3) :: sorted_indices
        
        call set_up_map(0.0,0.25,1440,EDGE_CENTERED,-90.0,0.25,721,EDGE_CENTERED,SPHERICAL,gridded_tbs,error)
        if (do_time) call set_up_map64(0.0,0.25,1440,EDGE_CENTERED,-90.0,0.25,721,EDGE_CENTERED,SPHERICAL,gridded_time,error)

        ! do jscan = 1,maxscan_extra
        !     if (.not. ieee_is_finite(scan_time_extra(jscan))) cycle
        !     max_valid_scan = jscan
        ! enddo
            
        do jscan = min_valid_extra_scan+num_overlap_scans,max_valid_extra_scan-num_overlap_scans
            do jfov = 1,maxcel_extra-1

                bad_loc = .false.
                do kscan=jscan,jscan+1
                    do kfov=jfov,jfov+1
                        if (.not. ieee_is_finite(cellat_extra(kfov,kscan))) bad_loc = .true.
                        if (.not. ieee_is_finite(cellon_extra(kfov,kscan))) bad_loc = .true.
                        if (.not. ieee_is_finite(tbr_extra(kfov,kscan,channel))) bad_loc = .true.
                        if (tbr_extra(kfov,kscan,channel) < 20.0) bad_loc = .true.
                    enddo
                enddo

                if (bad_loc) then
                    ! print *,kfov,kscan,channel
                    ! print *,cellat_extra(kfov,kscan)
                    ! print *,cellon_extra(kfov,kscan)
                    ! print *,tbr_extra(kfov,kscan,channel)
                    cycle  !either the lat/lon or Ta is bad
                endif

                !find lat/lon of center of the quadrilateral
                !and store vertices
                tot_lat = 0.0
                tot_lon = 0.0
                iquad = 0
                first_lon = cellon_extra(jfov,jscan)
                do kscan=jscan,jscan+1
                    do kfov=jfov,jfov+1
                        tot_lat = tot_lat + cellat_extra(kfov,kscan)
                        quad_lats(iquad) = cellat_extra(kfov,kscan)

                        if (abs(first_lon - cellon_extra(kfov,kscan)) .lt. 180.0) then
                            tot_lon = tot_lon + cellon_extra(kfov,kscan)
                            quad_lons(iquad) = cellon_extra(kfov,kscan)
                        else if (first_lon - cellon_extra(kfov,kscan) .gt. 180.0) then
                            ! this is the case where first_lon is near 360.0 and other lons are near 0.0
                            tot_lon = tot_lon + 360.0 + cellon_extra(kfov,kscan)
                            quad_lons(iquad) = cellon_extra(kfov,kscan)+360.0
                        else
                            ! this is the case where first lon is near 0.0, and cellon_extra(kfov,kscan) > 360.0
                            tot_lon = tot_lon - 360.0 + cellon_extra(kfov,kscan)
                            quad_lons(iquad) = cellon_extra(kfov,kscan)-360.0
                        endif
                        quad_ta(iquad) = tbr_extra(kfov,kscan,channel)
                        if (do_time) quad_times(iquad) = scan_time_extra(kscan)
                                                

                        ! these splices were for verification purposes
                        ! the results are in
                        !
                        !L:\access\amsr2_tb_orbits\r00001_05000\r00676.grid_interp_lat.05.nc
                        !L:\access\amsr2_tb_orbits\r00001_05000\r00676.grid_interp_lon.05.nc
                        !
                        ! verification case splice latitude
                        !quad_ta(iquad) = cellat_extra(kfov,kscan)
                        ! end verification case splice

                        ! verification case splice latitude
                        !quad_ta(iquad) = cellon_extra(kfov,kscan)
                        ! end verification case splice

                        iquad=iquad+1
                    enddo
                enddo

                mean_lat = tot_lat/4.0
                mean_lon = modulo((tot_lon/4.0) + 1080.0,360.0)

                !sort the quadrilateral points and the associated Ta's
                call order_vertices_quadrilateral(quad_lons,quad_lats,sorted_indices)
                do iquad=0,3
                    quad_lons_sorted(iquad) = quad_lons(sorted_indices(iquad))
                    quad_lats_sorted(iquad) = quad_lats(sorted_indices(iquad))
                    quad_ta_sorted(iquad)   = quad_ta(sorted_indices(iquad))
                    if (do_time) quad_times_sorted(iquad) = quad_times(sorted_indices(iquad))
                enddo

                is_convex = check_convex(quad_lons_sorted,quad_lats_sorted)
                if (.not. is_convex) then
                    print *,quad_lons_sorted
                    print *,quad_lats_sorted
                    print *,'Quadrilateral Not Convex - Should NOT happen - skipping'
                    cycle
                endif

                !find closest point in the target grid (for now a regular grid) to
                !the center lat/lon
                ! for a more complicated grid, this will need to be a search
                ! or an inverse transform
                klat0 = 1+nint(4.0*(90.0+mean_lat))
                klon0 = 1+nint(4.0*mean_lon)

                dk = 1
                ! search a larger longitude range close to the poles for 
                ! grid points inside the quadrilateral
                if (abs(mean_lat) .gt. 84.0) dk = 3
                if (abs(mean_lat) .gt. 87.0) dk = 7
                if (abs(mean_lat) .gt. 88.0) dk = 15
                do j =-1,1
                    do k = -dk,dk
                        klat = klat0+j 
                        klon = klon0+k 
                        call wrap_indices(gridded_tbs,klon0+k,klat0+j ,klon,klat,error)
                        call get_grid_location(gridded_tbs,klon,klat,grid_lon,grid_lat,interp_error)
                        
                        if (grid_lon < 10.0) then
                            do iquad=0,3
                                if (quad_lons_sorted(iquad) > 180.0) then
                                    quad_lons_sorted(iquad) = quad_lons_sorted(iquad) - 360.0
                                endif
                            enddo
                        endif
                        if (grid_lon > 350.0) then
                            do iquad=0,3
                                if (quad_lons_sorted(iquad) < 180.0) then
                                    quad_lons_sorted(iquad) = quad_lons_sorted(iquad) + 360.0
                                endif
                            enddo
                        endif
                        
                        !check to make sure that grid_lat,grid_lon is inside the quadrilateral
                        is_inside = is_inside_quadrilateral(quad_lons_sorted,       &
                                                            quad_lats_sorted,       &
                                                            real(grid_lon,real64),               &
                                                            real(grid_lat,real64),.false.) ! the false is for return_reordered

                        if (is_inside) then
                            !do the interpolation
                            call init_InterpBilinearQuadrilateral(     &
                                            quad_lons_sorted,   & ! array of 4 x locations
                                            quad_lats_sorted,   & ! array of 4 y locations
                                            quad_ta_sorted,     & ! array of 4 z locations
                                            qi,                 & ! interpolator instance
                                            interp_error)         ! error flag    
                            if (interp_error == 0) then 
                                gridded_tbs%dat(klon,klat) = eval_InterpBilinearQuadrilateral(real(grid_lon,real64),real(grid_lat,real64),qi,error)
                            endif

                            !do the time interpolation
                            if (do_time) then
                                call init_InterpBilinearQuadrilateral(     &
                                                quad_lons_sorted,   & ! array of 4 x locations
                                                quad_lats_sorted,   & ! array of 4 y locations
                                                quad_times_sorted,  & ! array of 4 z locations
                                                qi,                 & ! interpolator instance
                                                interp_error)         ! error flag    
                                if (interp_error == 0) then 
                                    gridded_time%dat(klon,klat) = eval_InterpBilinearQuadrilateral(real(grid_lon,real64),real(grid_lat,real64),qi,error)
                                endif
                            endif
                        endif
                    enddo
                enddo
            enddo !jfov
        enddo !jscan
    end subroutine resample_onto_lat_lon_grid_interp

    ! subroutine resample_onto_arbitrary_lat_lon_locations_interp(resolution,channel,do_time,gridded_time,gridded_tbs)

    !     ! This version of the resampling routine interpolates to the desired grid point using a general 
    !     ! quadrilateral interpolation method.

    !     ! It is different from resample_onto_lat_lon_grid_interp because the lats and lons in the map
    !     ! structure do not need to be on a rectangular grid.

    !     ! The looping is a little different than the previous versions
    !     !
    !     ! Here, the outside loops are over scans/fovs
    !     ! Then, the center point of the quadrilateral formed by (iscan,ifov),(iscan,ifov+1),(iscan+1,ifov+1),(iscan+1,ifov)
    !     ! is found.
    !     ! Then, the closest target point to this point is found and checked to see if it inside the quadrilateral
    !     ! If it is, then the interpolation to this point is done.
    !     !
    !     ! The algorithm also searches nearby fixed grid points to see if they are inside the quadrilateral
    !     ! The search range in the E-W direction expands near the pole so all possible grid points should be
    !     ! filled when using our standard 0.25 degree grids.

    !     use l2_module, only: re  !radius of Earth, not regular expressions (!)
    !     use l2_module_extra
    !     use QuadrilateralUtilties, only: check_convex, &
    !                                      is_inside_quadrilateral, &
    !                                      order_vertices_quadrilateral

    !     use InterpBilinearQuadrilateral, only: init_InterpBilinearQuadrilateral,eval_InterpBilinearQuadrilateral
    !     use InterpBilinearQuadrilateral, only: quad_interpolator_type

    !     implicit none

    !     real(real32),intent(in)     :: resolution
    !     integer(int32),intent(in)   :: channel
    !     logical,intent(in)          :: do_time
    !     type(two_d_map_str),intent(inout) :: gridded_time
    !     type(two_d_map_str),intent(inout) :: gridded_tbs

    !     !locals
    !     integer(int32)              :: ilat,ilon,closest_scan,closest_fov
    !     real(real32)                :: xlat,xlon,closest_distance
    !     integer(int32)              :: num_grid_points_filled
    !     integer(int32)              :: error
    !     integer(int32)              :: max_valid_scan

    !     integer(int32)              :: jscan,jfov,kscan,kfov,jlat,jlon,klat,klon,klat0,klon0,ilatx,ilonx,iquad,j,k,dk
    !     integer(int32)              :: ilon_closest,ilat_closest
    !     real(real32)                :: distance
    !     logical                     :: bad_loc,is_convex,is_inside
    !     integer(int32)              :: interp_error
    !     integer(int32)              :: i
    !     type(quad_interpolator_type):: qi

    !     real(real64)                :: tot_lat,tot_lon,first_lat,first_lon
    !     real(real32)                :: grid_lat,grid_lon
    !     real(real32),dimension(0:3) :: quad_lats,quad_lons,quad_ta
    !     real(real64),dimension(0:3) :: quad_times
    !     real(real64),dimension(0:3) :: quad_lats_sorted,quad_lons_sorted,quad_ta_sorted,quad_times_sorted
    !     real(real64)                :: mean_lat,mean_lon,tb_test
    !     integer(int32),dimension(0:3) :: sorted_indices
        

    !     call set_up_map_2d_map(0.0,0.25,1440,EDGE_CENTERED,-90.0,0.25,721,EDGE_CENTERED,SPHERICAL,gridded_tbs,error)
    !     if (do_time) call set_up_map_2d_map(0.0,0.25,1440,EDGE_CENTERED,-90.0,0.25,721,EDGE_CENTERED,SPHERICAL,gridded_time,error)

    !     do jscan = 1,maxscan_extra
    !         if (.not. ieee_is_finite(scan_time_extra(jscan))) cycle
    !         max_valid_scan = jscan
    !     enddo
            
    !     do jscan = 1+10,max_valid_scan-10
    !         do jfov = 1,maxcel_extra-1

    !             bad_loc = .false.
    !             do kscan=jscan,jscan+1
    !                 do kfov=jfov,jfov+1
    !                     if (.not. ieee_is_finite(cellat_extra(kfov,kscan))) bad_loc = .true.
    !                     if (.not. ieee_is_finite(cellon_extra(kfov,kscan))) bad_loc = .true.
    !                     if (.not. ieee_is_finite(tbr_extra(kfov,kscan,channel))) bad_loc = .true.
    !                     if (tbr_extra(kfov,kscan,channel) < 20.0) bad_loc = .true.
    !                 enddo
    !             enddo

    !             if (bad_loc) cycle  !either the lat/lon or Ta is bad

    !             !find lat/lon of center of the quadrilateral
    !             !and store vertices
    !             tot_lat = 0.0
    !             tot_lon = 0.0
    !             iquad = 0
    !             first_lon = cellon_extra(jfov,jscan)
    !             do kscan=jscan,jscan+1
    !                 do kfov=jfov,jfov+1
    !                     tot_lat = tot_lat + cellat_extra(kfov,kscan)
    !                     quad_lats(iquad) = cellat_extra(kfov,kscan)

    !                     if (abs(first_lon - cellon_extra(kfov,kscan)) .lt. 180.0) then
    !                         tot_lon = tot_lon + cellon_extra(kfov,kscan)
    !                         quad_lons(iquad) = cellon_extra(kfov,kscan)
    !                     else if (first_lon - cellon_extra(kfov,kscan) .gt. 180.0) then
    !                         ! this is the case where first_lon is near 360.0 and other lons are near 0.0
    !                         tot_lon = tot_lon + 360.0 + cellon_extra(kfov,kscan)
    !                         quad_lons(iquad) = cellon_extra(kfov,kscan)+360.0
    !                     else
    !                         ! this is the case where first lon is near 0.0, and cellon_extra(kfov,kscan) > 360.0
    !                         tot_lon = tot_lon - 360.0 + cellon_extra(kfov,kscan)
    !                         quad_lons(iquad) = cellon_extra(kfov,kscan)-360.0
    !                     endif
    !                     quad_ta(iquad) = tbr_extra(kfov,kscan,channel)
    !                     if (do_time) quad_times(iquad) = scan_time_extra(kscan)
                                                

    !                     ! these splices were for verification purposes
    !                     ! the results are in
    !                     !
    !                     !L:\access\amsr2_tb_orbits\r00001_05000\r00676.grid_interp_lat.05.nc
    !                     !L:\access\amsr2_tb_orbits\r00001_05000\r00676.grid_interp_lon.05.nc
    !                     !
    !                     ! verification case splice latitude
    !                     !quad_ta(iquad) = cellat_extra(kfov,kscan)
    !                     ! end verification case splice

    !                     ! verification case splice latitude
    !                     !quad_ta(iquad) = cellon_extra(kfov,kscan)
    !                     ! end verification case splice

    !                     iquad=iquad+1
    !                 enddo
    !             enddo

    !             mean_lat = tot_lat/4.0
    !             mean_lon = modulo((tot_lon/4.0) + 1080.0,360.0)

    !             !sort the quadrilateral points and the associated Ta's
    !             call order_vertices_quadrilateral(quad_lons,quad_lats,sorted_indices)
    !             do iquad=0,3
    !                 quad_lons_sorted(iquad) = quad_lons(sorted_indices(iquad))
    !                 quad_lats_sorted(iquad) = quad_lats(sorted_indices(iquad))
    !                 quad_ta_sorted(iquad)   = quad_ta(sorted_indices(iquad))
    !                 if (do_time) quad_times_sorted(iquad) = quad_times(sorted_indices(iquad))
    !             enddo

    !             is_convex = check_convex(quad_lons_sorted,quad_lats_sorted)
    !             if (.not. is_convex) then
    !                 print *,quad_lons_sorted
    !                 print *,quad_lats_sorted
    !                 print *,'Quadrilateral Not Convex - Should NOT happen'
    !                 error stop
    !             endif

    !             !find closest point in the target grid (for now a regular grid) to
    !             !the center lat/lon
    !             ! for a more complicated grid, this will need to be a search
    !             ! or an inverse transform
    !             klat0 = 1+nint(4.0*(90.0+mean_lat))
    !             klon0 = 1+nint(4.0*mean_lon)

    !             dk = 1
    !             ! search a larger longitude range close to the poles for 
    !             ! grid points inside the quadrilateral
    !             if (abs(mean_lat) .gt. 84.0) dk = 3
    !             if (abs(mean_lat) .gt. 87.0) dk = 7
    !             if (abs(mean_lat) .gt. 88.0) dk = 15
    !             do j =-1,1
    !                 do k = -dk,dk
    !                     klat = klat0+j 
    !                     klon = klon0+k 
    !                     !call wrap_indices_2d_map(gridded_tbs,klon0+k,klat0+j ,klon,klat,error)
    !                     call get_grid_location_2d_map(gridded_tbs,klon,klat,grid_lon,grid_lat,interp_error)
                        
    !                     if (grid_lon < 10.0) then
    !                         do iquad=0,3
    !                             if (quad_lons_sorted(iquad) > 180.0) then
    !                                 quad_lons_sorted(iquad) = quad_lons_sorted(iquad) - 360.0
    !                             endif
    !                         enddo
    !                     endif
    !                     if (grid_lon > 350.0) then
    !                         do iquad=0,3
    !                             if (quad_lons_sorted(iquad) < 180.0) then
    !                                 quad_lons_sorted(iquad) = quad_lons_sorted(iquad) + 360.0
    !                             endif
    !                         enddo
    !                     endif
                        
    !                     !check to make sure that grid_lat,grid_lon is inside the quadrilateral
    !                     is_inside = is_inside_quadrilateral(quad_lons_sorted,       &
    !                                                         quad_lats_sorted,       &
    !                                                         real(grid_lon,real64),               &
    !                                                         real(grid_lat,real64),.false.) ! the false is for return_reordered

    !                     if (is_inside) then
    !                         !do the interpolation
    !                         call init_InterpBilinearQuadrilateral(     &
    !                                         quad_lons_sorted,   & ! array of 4 x locations
    !                                         quad_lats_sorted,   & ! array of 4 y locations
    !                                         quad_ta_sorted,     & ! array of 4 z locations
    !                                         qi,                 & ! interpolator instance
    !                                         interp_error)         ! error flag    
    !                         if (interp_error == 0) then 
    !                             gridded_tbs%dat(klon,klat) = eval_InterpBilinearQuadrilateral(real(grid_lon,real64),real(grid_lat,real64),qi,error)
    !                         endif

    !                         !do the time interpolation
    !                         if (do_time) then
    !                             call init_InterpBilinearQuadrilateral(     &
    !                                             quad_lons_sorted,   & ! array of 4 x locations
    !                                             quad_lats_sorted,   & ! array of 4 y locations
    !                                             quad_times_sorted,  & ! array of 4 z locations
    !                                             qi,                 & ! interpolator instance
    !                                             interp_error)         ! error flag    
    !                             if (interp_error == 0) then 
    !                                 gridded_time%dat(klon,klat) = eval_InterpBilinearQuadrilateral(real(grid_lon,real64),real(grid_lat,real64),qi,error)
    !                             endif
    !                         endif
    !                     endif
    !                 enddo
    !             enddo
    !         enddo !jfov
    !     enddo !jscan
    ! end subroutine resample_onto_arbitrary_lat_lon_locations_interp

    subroutine resample_onto_lat_lon_grid(resolution,channel,gridded_time,gridded_tbs,gridded_distance)

        !This version, which finds the closest scan,fov each fixed grid point is deprecated

        !use l2_module, only: re
        use l2_module_extra
        use wgs84, only: distance_between_lon_lats

        real(real32),intent(in)     :: resolution
        integer(int32),intent(in)   :: channel
        type(map_str),intent(in)    :: gridded_time
        type(map_str),intent(inout) :: gridded_tbs
        type(map_str),intent(inout) :: gridded_distance

        !locals
        integer(int32)              :: ilat,ilon,closest_scan,closest_fov
        real(real32)                :: xlat,xlon,closest_distance
        integer(int32)              :: num_grid_points_filled
        integer(int32)              :: error
        integer(int32)              :: max_valid_scan

        integer(int32)              :: jscan,jfov,jlat,jlon,klat,klon,ilatx,ilonx
        integer(int32)              :: ilon_closest,ilat_closest
        real(real32)                :: distance
        

        call set_up_map(0.0,0.25,1440,EDGE_CENTERED,-90.0,0.25,721,EDGE_CENTERED,SPHERICAL,gridded_tbs,error)
        call set_up_map(0.0,0.25,1440,EDGE_CENTERED,-90.0,0.25,721,EDGE_CENTERED,SPHERICAL,gridded_distance,error)
        gridded_distance%dat = 999.9 !set the entire grid to a large number

        do jscan = 1,maxscan_extra
            if (.not. ieee_is_finite(scan_time_extra(jscan))) cycle
            max_valid_scan = jscan
        enddo
            
        do jscan = 1+50,max_valid_scan-350
            do jfov = 1,maxcel_extra
                if (.not. ieee_is_finite(cellat_extra(jfov,jscan))) cycle
                if (.not. ieee_is_finite(cellon_extra(jfov,jscan))) cycle
                

                klat = 1 + nint(4.0*(90.0+cellat_extra(jfov,jscan)))
                klon = 1 + nint(4.0*cellon_extra(jfov,jscan))

                ilon_closest = -999
                ilat_closest = -999
                do ilatx = klat-1,klat+1
                    do ilonx = klon-1,klon+1
                        call wrap_indices(gridded_tbs,ilonx,ilatx,ilon,ilat,error)
                        !find the distance between the scan location (possibly synthetic) and the location of
                        !the nearby fixed grid points
                        distance = distance_between_lon_lats(cellon_extra(jfov,jscan), &
                                                            cellat_extra(jfov,jscan), &
                                                            gridded_tbs%lons(ilon,ilat), &
                                                            gridded_tbs%lats(ilon,ilat), &
                                                            error)

                        if ( scan_time_extra(jscan) > (gridded_time%dat(ilon,ilat)+600.0)) cycle

                        if (distance < gridded_distance%dat(ilon,ilat)) then
                            if ((ieee_is_finite(tbr_extra(jfov,jscan,channel))) .and. &
                                (tbr_extra(jfov,jscan,channel) > 0.001)) then
                                    ilon_closest = ilon
                                    ilat_closest = ilat

                                ! replace values at this location
                                gridded_tbs%dat(ilon,ilat) = tbr_extra(jfov,  &
                                                                     jscan, &
                                                                     channel)
                                gridded_distance%dat(ilon,ilat) = distance
                                
                            endif
                        endif
                    enddo !ilon 
                enddo !ilat
            enddo !jfov
        enddo !jscan
    end subroutine resample_onto_lat_lon_grid

    subroutine write_resampled_maps_w_distance(gridded_tbs,&
                                    gridded_distance,  &
                                    filename_l2b_tb_grid,  &
                                    filename_l2b_tb_dist)

        

        type(map_str),intent(inout) :: gridded_tbs
        type(map_str),intent(inout) :: gridded_distance
        character(*),intent(in)     :: filename_l2b_tb_grid,filename_l2b_tb_dist

        integer(int32)              :: error
        character(150)              :: title

        title = 'Gridded Tbs'
        !call write_map_netcdf(gridded_tbs,filename_l2b_tb_grid,error,trim(title))
        title = 'Distance to regridded swath location'
        !call write_map_netcdf(gridded_distance,filename_l2b_tb_dist,error,trim(title))
    end subroutine write_resampled_maps_w_distance

    subroutine write_resampled_maps(gridded_tbs,&
                                    filename_l2b_tb_grid)

        use equirectangular_maps, only:map_str,write_map_netcdf


        type(map_str),intent(inout) :: gridded_tbs
        character(*),intent(in)     :: filename_l2b_tb_grid

        integer(int32)              :: error
        character(150)              :: title

        title = 'Gridded Tbs'
        call write_map_netcdf(gridded_tbs,filename_l2b_tb_grid,error,trim(title))
    end subroutine write_resampled_maps

    subroutine write_wt_slice_netcdf(weights,freqname,scan,fov)
        use sat_id_module
        use io_nc, only: handle_nc_err, minmax_iso8601
        use netcdf

        type(ResampleWeights),intent(in)  :: weights
        character(*),intent(in)           :: freqname
        integer(int32),intent(in)         :: scan
        integer(int32),intent(in)         :: fov

        character(200)                    :: nc_filename,path
        integer(int32)                    :: ncid,varid,dim_scan,dim_fov,i

        path = 'L:\access\resampling\AMSR2\resample_weights\circular_30km\netcdf\'
        write(nc_filename,111) trim(path),freqname,scan,fov
        111 format(a,'resample_weight_slice_scan_',a,'_',i3.3,'_fov_',i3.3,'.nc')

        call handle_nc_err(nf90_create(nc_filename, ior(NF90_CLOBBER, NF90_NETCDF4), ncid))
        ! Define global attributes
        call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "institution", "REMSS"))
        call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "creator_name", "Remote Sensing Systems"))
        call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "creator_email", "support@remss.com"))
        call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "creator_url", "http://www.remss.com"))
        
        call handle_nc_err(nf90_def_dim(ncid, "scan", weights%num_input_scans, dim_scan))
        call handle_nc_err(nf90_def_dim(ncid, "fov", weights%num_input_cells, dim_fov))

        call handle_nc_err(nf90_def_var(ncid, "scan", NF90_INT,dim_scan, varid))
        call handle_nc_err(nf90_put_var(ncid, varid, [(i, i = -14, 14)]))

        call handle_nc_err(nf90_def_var(ncid, "fov", NF90_INT, dim_fov, varid))
        call handle_nc_err(nf90_put_var(ncid, varid, [(i, i = 1,weights%num_input_cells)]))

        call handle_nc_err(nf90_def_var(ncid, "wts", NF90_DOUBLE, [dim_fov,dim_scan], varid, &
                            deflate_level=2, shuffle=.true.))
        call handle_nc_err(nf90_put_var(ncid, varid, weights%wts(:,:,fov,scan)))
        call handle_nc_err(nf90_close(ncid))
    end subroutine write_wt_slice_netcdf

    subroutine resample_time_onto_lat_lon_grid(resolution, &
                                               gridded_time, &
                                               gridded_distance, &
                                               min_valid_extra_scan, &
                                               max_valid_extra_scan, &
                                               num_overlap_scans)

        !use l2_module, only: re
        use l2_module_extra
        use wgs84, only: distance_between_lon_lats
        use equirectangular_maps64

        real(real32),intent(in)     :: resolution
        type(map_str64),intent(inout) :: gridded_time
        type(map_str),intent(inout) :: gridded_distance
        integer(int32),intent(in)   :: min_valid_extra_scan
        integer(int32),intent(in)   :: max_valid_extra_scan
        integer(int32),intent(in)   :: num_overlap_scans
        

        !locals
        integer(int32)              :: ilat,ilon,closest_scan,closest_fov
        real(real32)                :: xlat,xlon,closest_distance
        integer(int32)              :: num_grid_points_filled
        integer(int32)              :: error

        integer(int32)              :: jscan,jfov,jlat,jlon,klat,klon,ilatx,ilonx,dklon
        real(real32)                :: distance
        integer(int32)              :: max_valid_scan

        type(map_str)               :: filled_during_1st_half

        call set_up_map64(0.0,0.25,1440,EDGE_CENTERED,-90.0,0.25,721,EDGE_CENTERED,SPHERICAL,gridded_time,error)
        call set_up_map(0.0,0.25,1440,EDGE_CENTERED,-90.0,0.25,721,EDGE_CENTERED,SPHERICAL,gridded_distance,error)
        call set_up_map(0.0,0.25,1440,EDGE_CENTERED,-90.0,0.25,721,EDGE_CENTERED,SPHERICAL,filled_during_1st_half,error)
        
        gridded_distance%dat = 999.9
        filled_during_1st_half%dat = 0.0

        do jscan = 1,maxscan_extra
            if (.not. ieee_is_finite(scan_time_extra(jscan))) cycle
            max_valid_scan = jscan
        enddo
        
        do jscan = 1+num_overlap_scans,max_valid_extra_scan-num_overlap_scans
            do jfov = 1,maxcel_extra
                if (.not. ieee_is_finite(cellat_extra(jfov,jscan))) cycle
                if (.not. ieee_is_finite(cellon_extra(jfov,jscan))) cycle
                
                klat = 1 + nint(4.0*(90.0+cellat_extra(jfov,jscan)))
                klon = 1 + nint(4.0*cellon_extra(jfov,jscan))
                dklon = 1
                if (abs(cellat_extra(jfov,jscan)) .gt. 87.) dklon=5
                do ilatx = klat-1,klat+1
                    do ilonx = klon-dklon,klon+dklon
                        call wrap_indices64(gridded_time,ilonx,ilatx,ilon,ilat,error)
                        distance = distance_between_lon_lats(cellon_extra(jfov,jscan), &
                                                            cellat_extra(jfov,jscan), &
                                                            gridded_time%lons(ilon,ilat), &
                                                            gridded_time%lats(ilon,ilat), &
                                                            error)
                        if (distance < gridded_distance%dat(ilon,ilat)) then
                            if (jscan .lt. maxscan_extra/2) then
                                if (ieee_is_finite(scan_time_extra(jscan))) then
                                    ! replace values at this location
                                    gridded_time%dat(ilon,ilat) = scan_time_extra(jscan)
                                    gridded_distance%dat(ilon,ilat) = distance
                                    filled_during_1st_half%dat(ilon,ilat) = 1.0  
                                endif
                            else
                                if ((ieee_is_finite(scan_time_extra(jscan))) .and. &
                                   (filled_during_1st_half%dat(ilon,ilat) .lt. 0.1)) then
                                    ! replace values at this location only if not filled in first half!
                                        gridded_time%dat(ilon,ilat) = scan_time_extra(jscan)
                                        gridded_distance%dat(ilon,ilat) = distance
                                endif
                            endif
                        endif
                    enddo !ilon 
                enddo !ilat
            enddo !jfov
        enddo !jscan
    end subroutine resample_time_onto_lat_lon_grid

    subroutine write_resampled_time_map(gridded_time,&
                                        filename_l2b_time_grid)

        use equirectangular_maps64, only:map_str64,write_map_netcdf64

        type(map_str64),intent(inout) :: gridded_time
        character(*),intent(in)     :: filename_l2b_time_grid

        integer(int32)              :: error
        character(150)              :: title

        title = 'Gridded Time'
        call write_map_netcdf64(gridded_time,filename_l2b_time_grid,error,trim(title))
        !write_map_netcdf(mp,nc_filename,error,cf_var_metadata,cf_global_metadata)
      
    end subroutine write_resampled_time_map

end module resample


