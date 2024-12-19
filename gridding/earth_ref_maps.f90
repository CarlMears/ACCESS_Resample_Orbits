! This is a module that contains generalized 2-D maps on the Earths surface.
! The x and y coordinates are regularly spaced, but
! the latitude and longitude are arbitrary.
! the purpose is to hold things like regaular lat lon maps AND 
! polar-gridded maps.


module earth_ref_maps
    use, intrinsic :: iso_fortran_env, only: int8, int32, real32, real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_finite
    use NAN_support, only: nan_f32, nan_f64
    use QuadrilateralUtilties, only: is_inside_quadrilateral
    use io_nc, only: handle_nc_err, minmax_iso8601
    use netcdf
    !use CF_metadata, only: CF_1_8_attr,init_CF_1_8_attr,put_CF_1_8_attr
    
    implicit none

    integer(int32),parameter                :: CELL_CENTERED_2D_MAP = 1
    integer(int32),parameter                :: EDGE_CENTERED_2D_MAP = 2
    real(real32),parameter                  :: SMALL_NUMBER32  = 1.0e-3
    real(real64),parameter                  :: SMALL_NUMBER64  = 1.0d-3


    type map_axis
        real(real32)    :: min
        real(real32)    :: max
        real(real32)    :: delta
        integer(int32)  :: num_pts
        integer(int32)  :: axis_type
    end type map_axis

    type earth_ref_map_str
        logical         :: initialized = .false.
        logical         :: allocated = .false.
        type(map_axis)  :: x_axis
        type(map_axis)  :: y_axis
        real(real32),dimension(:,:),allocatable :: lons          ! lons
        real(real32),dimension(:,:),allocatable :: lats          ! lats
    contains
        procedure :: init => init_2d_map_base
        procedure :: wrap_indices => wrap_indices_2d_map_base
        procedure :: get_grid_location => get_grid_location_base
        procedure :: get_grid_index => get_grid_index_base
        procedure :: check_x_y => check_x_y_2d_map 
    end type earth_ref_map_str

    type, extends(earth_ref_map_str) :: earth_ref_map_str32
        real(real32),allocatable :: dat(:,:)
    contains
        procedure :: init => init_2d_map32
        procedure :: define_data => define_data32
        procedure :: define_lats_lons => define_lats_lons32
        procedure :: zero_data => zero_2d_map32
        procedure :: set_to_constant => set_2d_map32_to_constant
        procedure :: get_interp_value => get_2d_map_value32
        procedure :: get_closest_value => get_2d_map_value_no_interp32
        procedure :: get_lat_lon => map_lat_lon32
        procedure :: is_inside_quad => is_inside_quad_map32
        procedure :: deinit => deinit_2d_map32
        procedure :: write_map_netcdf => write_map32_netcdf
    end type earth_ref_map_str32

    type, extends(earth_ref_map_str) :: earth_ref_map_str64
        real(real64),allocatable :: dat(:,:)
    contains
        procedure :: init => init_2d_map64
        procedure :: define_data => define_data64
        procedure :: define_lats_lons => define_lats_lons64
        procedure :: zero_data => zero_2d_map64
        procedure :: set_to_constant => set_2d_map64_to_constant
        procedure :: get_interp_value => get_2d_map_value64
        procedure :: get_closest_value => get_2d_map_value_no_interp64
        procedure :: get_lat_lon => map_lat_lon64
        procedure :: is_inside_quad => is_inside_quad_map64
        procedure :: deinit => deinit_2d_map64
        procedure :: write_map_netcdf => write_map64_netcdf
    end type earth_ref_map_str64

    !public earth_ref_map_str  !I don't think the base object needs to be exposed.
    public earth_ref_map_str32
    public earth_ref_map_str64

    ! explanation of the grid type flag:
    !
    !  There are basically two types of evenly spaced rectangular gridded map --
    !
    !        1. Those for which the value represents the value in a cell - I'll call this cell_centered
    !              for these maps, the lower left value represents the value in the cell defined
    !             by long = (0.0:delta_x), lat = (-90,-90+delta_y), and is associated with at point at
    !              (delta_x/2.0, -90+delta_x/2.0)
    !        2. Those for which the value represents the value on the vertex between cells - I'll call this edge_centered
    !              for these maps, the lower left value represents the value at the point 
    !              by long = 0.0, lat = (-90).  A peculiar feature of this sort of grid is
    !              that the row at -90.0 must be identical, since it all refers to the same point!
    !        For now, all routines in the module assume that maps obey spherical boundary conditions
    !        For now, all maps are not redundant at the vertical edges
    !        
    !        This means that for both types of grids, 360.0/(delta_x) longitude points are needed, and the only different is
    !        whether the values are offset by 0.5*delta_x
    !
    !        For cell_centered grids, 180/(delta_y) latitude points are needed, and for 
    !        Edge centered grids 180/(delta_y) + 1 latitude points are needed
    !
    !        note that mixed grid types are allowed, with the x and y coordinates being different types
    contains

        subroutine check_map_parameters(min_x,max_x,num_x,type_x,min_y,max_y,num_y,type_y,parameter_error)
            real(real32),intent(IN)       :: min_x
            real(real32),intent(IN)       :: max_x
            integer(int32),intent(IN)     :: num_x
            integer(int32),intent(IN)     :: type_x

            real(real32),intent(IN)       :: min_y
            real(real32),intent(IN)       :: max_y
            integer(int32),intent(IN)     :: num_y 
            integer(int32),intent(IN)     :: type_y

            integer(int32),intent(OUT)    :: parameter_error

            parameter_error = 0

            if (.not. ieee_is_finite(min_x)) parameter_error = -1
            if (.not. ieee_is_finite(max_x)) parameter_error = -1
            if (num_x .le. 0) parameter_error=-1
            if ((type_x .ne. CELL_CENTERED_2D_MAP) .and. &
                (type_x .ne. EDGE_CENTERED_2D_MAP)) parameter_error = -1

            if (.not. ieee_is_finite(min_y)) parameter_error = -1
            if (.not. ieee_is_finite(max_y)) parameter_error = -1
            if (num_y .le. 0) parameter_error=-1
            if ((type_y .ne. CELL_CENTERED_2D_MAP) .and. &
                (type_y .ne. EDGE_CENTERED_2D_MAP)) parameter_error = -1
            
        end subroutine check_map_parameters

        subroutine init_2d_map_base(mp,min_x,max_x,num_x,type_x,min_y,max_y,num_y,type_y,error)

            use NAN_support

            implicit none

            class(earth_ref_map_str),intent(INOUT)    :: mp

            real(real32),intent(IN)       :: min_x
            real(real32),intent(IN)       :: max_x
            integer(int32),intent(IN)     :: num_x
            integer(int32),intent(IN)     :: type_x

            real(real32),intent(IN)       :: min_y
            real(real32),intent(IN)       :: max_y
            integer(int32),intent(IN)     :: num_y 
            integer(int32),intent(IN)     :: type_y
            integer(int32),intent(OUT)    :: error
            
            integer(int32)                :: parameter_error

            call check_map_parameters(min_x,max_x,num_x,type_x, &
                                      min_y,max_y,num_y,type_y,parameter_error)

            if (parameter_error .ne. 0) then
                error = parameter_error
                return
            endif

            mp%x_axis%min = min_x
            mp%x_axis%max = max_x
            mp%x_axis%num_pts = num_x
            mp%x_axis%axis_type = type_x

            mp%y_axis%min = min_y
            mp%y_axis%max = max_y
            mp%y_axis%num_pts = num_y
            mp%y_axis%axis_type = type_y

            if (mp%x_axis%axis_type .eq. CELL_CENTERED_2D_MAP) then
                mp%x_axis%delta = (mp%x_axis%max - mp%x_axis%min)/(mp%x_axis%num_pts-1)
            endif

            if (mp%x_axis%axis_type .eq. EDGE_CENTERED_2D_MAP) then
                mp%x_axis%delta = (mp%x_axis%max - mp%x_axis%min)/(mp%x_axis%num_pts-1)
            endif

            if (mp%y_axis%axis_type .eq. CELL_CENTERED_2D_MAP) then
                mp%y_axis%delta = (mp%y_axis%max - mp%y_axis%min)/(mp%y_axis%num_pts-1)
            endif

            if (mp%y_axis%axis_type .eq. EDGE_CENTERED_2D_MAP) then
                mp%y_axis%delta = (mp%y_axis%max - mp%y_axis%min)/(mp%y_axis%num_pts-1)
            endif

            allocate(mp%lats(mp%x_axis%num_pts,mp%y_axis%num_pts),STAT = error)
            allocate(mp%lons(mp%x_axis%num_pts,mp%y_axis%num_pts),STAT = error)
        
            mp%initialized = .true.
    
            return

        end subroutine init_2d_map_base

        subroutine wrap_indices_2d_map_base(mp,ix,iy,ix_wrapped,iy_wrapped,error)

            implicit none

            class(earth_ref_map_str),intent(IN)  :: mp
            integer(int32),intent(IN)            :: ix
            integer(int32),intent(IN)            :: iy
            integer(int32),intent(OUT)           :: ix_wrapped
            integer(int32),intent(OUT)           :: iy_wrapped
            integer(int32),intent(OUT)           :: error


            !for the base map object, no wrapping is done
            !anything outside the bounds causes error tp be set to -1
            if ((ix .lt. 1) .or. (ix .gt. mp%x_axis%num_pts) .or. &
               (iy .lt. 1) .or. (iy .gt. mp%y_axis%num_pts)) then
                ix_wrapped = 0
                iy_wrapped = 0     
                error = -1
            else
                ix_wrapped = ix
                iy_wrapped = iy
                error = 0
            endif
            return

        end subroutine wrap_indices_2d_map_base
        
        subroutine init_2d_map32(mp,min_x,max_x,num_x,type_x,min_y,max_y,num_y,type_y,error)

            use NAN_support

            implicit none

            class(earth_ref_map_str32),intent(INOUT)    :: mp

            real(real32),intent(IN)       :: min_x
            real(real32),intent(IN)       :: max_x
            integer(int32),intent(IN)     :: num_x
            integer(int32),intent(IN)     :: type_x

            real(real32),intent(IN)       :: min_y
            real(real32),intent(IN)       :: max_y
            integer(int32),intent(IN)     :: num_y 
            integer(int32),intent(IN)     :: type_y
            integer(int32),intent(OUT)    :: error
            
            call init_2d_map_base(mp,min_x,max_x,num_x,type_x,min_y,max_y,num_y,type_y,error)

            allocate(mp%dat(mp%x_axis%num_pts,mp%y_axis%num_pts),STAT = error)
            mp%allocated = .true.
    
            return

        end subroutine init_2d_map32

        subroutine init_2d_map64(mp,min_x,max_x,num_x,type_x,min_y,max_y,num_y,type_y,error)

            use NAN_support

            implicit none

            class(earth_ref_map_str64),intent(INOUT)    :: mp

            real(real32),intent(IN)       :: min_x
            real(real32),intent(IN)       :: max_x
            integer(int32),intent(IN)     :: num_x
            integer(int32),intent(IN)     :: type_x

            real(real32),intent(IN)       :: min_y
            real(real32),intent(IN)       :: max_y
            integer(int32),intent(IN)     :: num_y 
            integer(int32),intent(IN)     :: type_y
            integer(int32),intent(OUT)    :: error
            
            call init_2d_map_base(mp,min_x,max_x,num_x,type_x,min_y,max_y,num_y,type_y,error)

            allocate(mp%dat(mp%x_axis%num_pts,mp%y_axis%num_pts),STAT = error)
            mp%allocated = .true.
    
            return

        end subroutine init_2d_map64

        
        ! subroutine define_lats_base(mp,lats,lons,error)
        !     class(earth_ref_map_str),intent(INOUT) :: mp
        !     real(real32),dimension(:,:),intent(IN) :: lats
        !     real(real32),dimension(:,:),intent(IN) :: lons
        !     integer(int32),intent(OUT)    :: error

        !     error = 0

        !     mp%lats(:,:) = lats
        !     mp%lons(:,:) = lons

        !     return

        ! end subroutine define_lats_base

        subroutine define_lats_lons32(mp,lats,lons,error)
            class(earth_ref_map_str32),intent(INOUT) :: mp
            real(real32),dimension(:,:),intent(IN) :: lats
            real(real32),dimension(:,:),intent(IN) :: lons
            integer(int32),intent(OUT)    :: error

            error = 0

            mp%lats(:,:) = lats
            mp%lons(:,:) = lons

            return

        end subroutine define_lats_lons32

        subroutine define_lats_lons64(mp,lats,lons,error)
            class(earth_ref_map_str64),intent(INOUT) :: mp
            real(real32),dimension(:,:),intent(IN) :: lats
            real(real32),dimension(:,:),intent(IN) :: lons
            integer(int32),intent(OUT)    :: error

            error = 0

            mp%lats(:,:) = lats
            mp%lons(:,:) = lons

            return

        end subroutine define_lats_lons64

        subroutine define_data32(mp,data_vals,error)
            class(earth_ref_map_str32),intent(INOUT)    :: mp
            real(real32),dimension(:,:),intent(IN)      :: data_vals
            integer(int32),intent(OUT)    :: error

            error = 0

            mp%dat(:,:) = data_vals

            return

        end subroutine define_data32

        subroutine define_data64(mp,data_vals,error)
            class(earth_ref_map_str64),intent(INOUT)    :: mp
            real(real64),dimension(:,:),intent(IN) :: data_vals
            integer(int32),intent(OUT)    :: error

            error = 0

            mp%dat(:,:) = data_vals

            return

        end subroutine define_data64

        subroutine zero_2d_map32(mp,error)

            implicit none

            class(earth_ref_map_str32),intent(INOUT):: mp
            integer(int32),intent(INOUT)            :: error

            mp%dat = 0.0
            error = 0

        end subroutine zero_2d_map32

        subroutine zero_2d_map64(mp,error)

            implicit none

            class(earth_ref_map_str64),intent(INOUT):: mp
            integer(int32),intent(INOUT)            :: error

            mp%dat = 0.0
            error = 0

        end subroutine zero_2d_map64

        subroutine set_2d_map32_to_constant(mp,val,error)

            implicit none

            class(earth_ref_map_str32),intent(INOUT):: mp
            real(real32),intent(IN)                 :: val
            integer(int32),intent(INOUT)            :: error

            mp%dat = val
            error = 0

        end subroutine set_2d_map32_to_constant

        subroutine set_2d_map64_to_constant(mp,val,error)

            implicit none

            class(earth_ref_map_str64),intent(INOUT):: mp
            real(real64),intent(IN)                 :: val
            integer(int32),intent(INOUT)            :: error

            mp%dat = val
            error = 0

        end subroutine set_2d_map64_to_constant

        subroutine deinit_2d_map_base(mp,error)
            
            implicit none

            class(earth_ref_map_str),intent(INOUT)        :: mp
            integer(int32),intent(OUT)         :: error

            mp%x_axis%min    = nan_f32
            mp%x_axis%max    = nan_f32
            mp%x_axis%delta  = nan_f32
            mp%x_axis%num_pts    = 0
            mp%x_axis%axis_type   = 0
            mp%y_axis%min    = nan_f32
            mp%y_axis%max    = nan_f32
            mp%y_axis%delta  = nan_f32
            mp%y_axis%num_pts    = 0
            mp%y_axis%axis_type   = 0

            if (allocated(mp%lats)) then
                deallocate(mp%lats)
            endif 

            if (allocated(mp%lons)) then
                deallocate(mp%lons)
            endif 

            mp%initialized = .false.
            error = 0

            return

        end subroutine deinit_2d_map_base

        subroutine deinit_2d_map32(mp,error)
            
            implicit none

            class(earth_ref_map_str32),intent(INOUT)        :: mp
            integer(int32),intent(OUT)         :: error

            call deinit_2d_map_base(mp,error)
            if (allocated(mp%dat)) then
                deallocate(mp%dat)
            endif

            mp%allocated = .false.

            return
        end subroutine deinit_2d_map32

        subroutine deinit_2d_map64(mp,error)
            
            implicit none

            class(earth_ref_map_str64),intent(INOUT)        :: mp
            integer(int32),intent(OUT)         :: error

            call deinit_2d_map_base(mp,error)

            if (allocated(mp%dat)) then
                deallocate(mp%dat)
            endif
            mp%allocated = .false.

        end subroutine deinit_2d_map64

        subroutine check_x_y_2d_map(mp,x,y,error)
            class(earth_ref_map_str),intent(IN)    :: mp      ! map structure to be interpolated
            real(real32),intent(IN)                :: x       ! x, in x,y space
            real(real32),intent(IN)                :: y       ! y, in x,y space
            integer(int32),intent(OUT)             :: error

            error = 0

            if ((.not. ieee_is_finite(x)) .or. &
                (.not. ieee_is_finite(y))) then
                error = -1
                return
            endif

            if ((x < mp%x_axis%min) .or. &
                (x > mp%x_axis%max) .or. &
                (y < mp%y_axis%min) .or. &
                (y > mp%y_axis%max)) then
                error = -1
                return
            endif

            if (mp%x_axis%axis_type .eq. CELL_CENTERED_2D_MAP) then
                if (((x - mp%x_axis%min) .lt. mp%x_axis%delta/2.0) .or. &
                    ((mp%x_axis%max - x) .lt. mp%x_axis%delta/2.0)) then
                    error = 1 ! warning -- too close to edge x
                endif
            endif

            if (mp%y_axis%axis_type .eq. CELL_CENTERED_2D_MAP) then
                if (((y - mp%y_axis%min) .lt. mp%y_axis%delta/2.0) .or. &
                    ((mp%y_axis%max - y) .lt. mp%y_axis%delta/2.0)) then
                    error = error + 2 ! warning -- too close to edge for y
                endif
            endif

        end subroutine check_x_y_2d_map

        subroutine get_grid_location_base(mp,ix,iy,x,y,error)
            class(earth_ref_map_str),intent(IN)   :: mp
            integer(int32),intent(IN)             :: ix
            integer(int32),intent(IN)             :: iy
            real(real32),intent(OUT)              :: x
            real(real32),intent(OUT)              :: y
            integer(int32),intent(OUT)            :: error

            if ((ix < 1) .or. (ix > mp%x_axis%num_pts) .or. &
                (iy < 1) .or. (iy > mp%y_axis%num_pts)) then
                x = nan_f32
                y = nan_f32
                error = 1
                return
            endif

            x = mp%x_axis%min + (ix-1)*mp%x_axis%delta
            y = mp%y_axis%min + (iy-1)*mp%y_axis%delta
            error = 0
        end subroutine

        subroutine get_grid_index_base(mp,x,y,ix,iy,error)
            class(earth_ref_map_str),intent(IN)   :: mp
            real(real32),intent(IN)              :: x
            real(real32),intent(IN)              :: y
            integer(int32),intent(OUT)             :: ix
            integer(int32),intent(OUT)             :: iy
            integer(int32),intent(OUT)            :: error

            !returns the index of the cell that contains the point x,y

            ix = 1+floor((x-mp%x_axis%min)/mp%x_axis%delta + 0.5)
            iy = 1+floor((y-mp%y_axis%min)/mp%y_axis%delta + 0.5)

            if ((ix < 1) .or. (ix > mp%x_axis%num_pts) .or. &
                (iy < 1) .or. (iy > mp%y_axis%num_pts)) then
                error = 1
                return
            endif

            error = 0
        end subroutine get_grid_index_base




        subroutine get_2d_map_value32(mp,x,y,value,error)

            ! this subroutine returns the interpolated value at a point x,y
            ! note that the interpolation is somewhat different depending on 
            ! what type of grid is used for the map structure

            implicit none

            class(earth_ref_map_str32),intent(IN)  :: mp       ! map structure to be interpolated
            real(real32),intent(IN)                :: x       ! x, in x,y space
            real(real32),intent(IN)                :: y       ! y, in x,y space
            real(real32),intent(OUT)               :: value    ! the interpolated value
            integer(int32),intent(OUT)             :: error

            integer(int32)                         :: ix_left,ix_right        ! map indices of points to the right (left) of desired point
            integer(int32)                         :: iy_lower,iy_upper       ! map indices of points above (below) the desired point
            real(real32)                           :: dx,dy                   ! distance of the desired point to the
                                                                            ! right (above) the lower left grid point
            real(real32)                           :: wt_left,wt_right
            real(real32)                           :: wt_lower,wt_upper
            real(real32)                           :: tot,tot_wt,wt
            real(real32)                           :: v

            error = 0

            ! make sure x and y are finite and within bounds

            call check_x_y_2d_map(mp,x,y,error)
            if (error .ne. 0) then
                value = nan_f32
                return
            endif

            ! find x indices and distance
            select case (mp%x_axis%axis_type)
                case (CELL_CENTERED_2D_MAP) 
                    ix_left        = floor((x-mp%x_axis%min)/mp%x_axis%delta - 0.5) + 1
                    ix_right    = ix_left + 1
                    dx            = x - mp%x_axis%min - ((ix_left - 0.5) * mp%x_axis%delta)
                case (EDGE_CENTERED_2D_MAP) 
                    ix_left  = floor((x-mp%x_axis%min)/mp%x_axis%delta) + 1
                    ix_right = ix_left + 1
                    dx       = x - mp%x_axis%min - ((ix_left - 1) * mp%x_axis%delta)         
                case default
                    error = 1
                    value = nan_f32
                    return
            end select

            ! find y indices and distance

            select case(mp%y_axis%axis_type)
                case (CELL_CENTERED_2D_MAP)
                    iy_lower    = floor((y-mp%y_axis%min)/mp%y_axis%delta - 0.5) + 1
                    iy_upper    = iy_lower + 1
                    dy            = y - mp%y_axis%min - ((iy_lower - 0.5) * mp%y_axis%delta)
                case (EDGE_CENTERED_2D_MAP)
                    iy_lower    = floor((y-mp%y_axis%min)/mp%y_axis%delta) + 1
                    iy_upper    = iy_lower + 1
                    dy            = y - mp%y_axis%min - ((iy_lower - 1) * mp%y_axis%delta)
                case default
                    error = 1
                    return
            end select

            ! now perform the interpolation

            wt_right = dx/mp%x_axis%delta
            wt_left  = 1.0 - wt_right

            wt_upper = dy/mp%y_axis%delta
            wt_lower = 1.0 - wt_upper

            tot = 0.0 
            tot_wt = 0.0

            v = map_value32(mp,ix_left,iy_lower,error) 
            if (ieee_is_finite(v)) then
                wt = wt_left * wt_lower
                tot = tot +( v * wt)
                tot_wt = tot_wt+wt
            endif 

            v = map_value32(mp,ix_left,iy_upper,error) 
            if (ieee_is_finite(v)) then
                wt = wt_left * wt_upper
                tot = tot + (v * wt)
                tot_wt = tot_wt+wt
            endif

            v = map_value32(mp,ix_right,iy_upper,error) 
            if (ieee_is_finite(v)) then
                wt = wt_right*wt_upper
                tot = tot + (v * wt)
                tot_wt = tot_wt + wt
            endif
            
            v = map_value32(mp,ix_right,iy_lower,error) 
            if (ieee_is_finite(v)) then
                wt = wt_right*wt_lower 
                tot = tot + (v * wt)
                tot_wt = tot_wt+wt
            endif 

            value = tot/tot_wt

            return

        end subroutine get_2d_map_value32

        subroutine get_2d_map_value64(mp,x,y,value,error)

            ! this subroutine returns the interpolated value at a point x,y
            ! note that the interpolation is somewhat different depending on 
            ! what type of grid is used for the map structure

            implicit none

            class(earth_ref_map_str64),intent(IN)  :: mp       ! map structure to be interpolated
            real(real32),intent(IN)                :: x       ! x, in x,y space
            real(real32),intent(IN)                :: y       ! y, in x,y space
            real(real64),intent(OUT)               :: value    ! the interpolated value
            integer(int32),intent(OUT)             :: error

            integer(int32)                         :: ix_left,ix_right        ! map indices of points to the right (left) of desired point
            integer(int32)                         :: iy_lower,iy_upper       ! map indices of points above (below) the desired point
            real(real32)                           :: dx,dy                   ! distance of the desired point to the
                                                                            ! right (above) the lower left grid point
            real(real64)                           :: wt_left,wt_right
            real(real64)                           :: wt_lower,wt_upper
            real(real64)                           :: tot,tot_wt,wt
            real(real64)                           :: v

            error = 0

            ! make sure x and y are finite and within bounds

            call check_x_y_2d_map(mp,x,y,error)
            if (error .ne. 0) then
                value = nan_f64
                return
            endif

            ! find x indices and distance
            select case (mp%x_axis%axis_type)
                case (CELL_CENTERED_2D_MAP) 
                    ix_left        = floor((x-mp%x_axis%min)/mp%x_axis%delta - 0.5) + 1
                    ix_right    = ix_left + 1
                    dx            = x - mp%x_axis%min - ((ix_left - 0.5) * mp%x_axis%delta)
                case (EDGE_CENTERED_2D_MAP) 
                    ix_left  = floor((x-mp%x_axis%min)/mp%x_axis%delta) + 1
                    ix_right = ix_left + 1
                    dx       = x - mp%x_axis%min - ((ix_left - 1) * mp%x_axis%delta)         
                case default
                    error = 1
                    value = nan_f64
                    return
            end select

            ! find y indices and distance

            select case(mp%y_axis%axis_type)
                case (CELL_CENTERED_2D_MAP)
                    iy_lower    = floor((y-mp%y_axis%min)/mp%y_axis%delta - 0.5) + 1
                    iy_upper    = iy_lower + 1
                    dy            = y - mp%y_axis%min - ((iy_lower - 0.5) * mp%y_axis%delta)
                case (EDGE_CENTERED_2D_MAP)
                    iy_lower    = floor((y-mp%y_axis%min)/mp%y_axis%delta) + 1
                    iy_upper    = iy_lower + 1
                    dy            = y - mp%y_axis%min - ((iy_lower - 1) * mp%y_axis%delta)
                case default
                    error = 1
                    return
            end select

            ! now perform the interpolation

            wt_right = dx/mp%x_axis%delta
            wt_left  = 1.0 - wt_right

            wt_upper = dy/mp%y_axis%delta
            wt_lower = 1.0 - wt_upper

            tot = 0.0 
            tot_wt = 0.0

            v = map_value64(mp,ix_left,iy_lower,error) 
            if (ieee_is_finite(v)) then
                wt = wt_left * wt_lower
                tot = tot +( v * wt)
                tot_wt = tot_wt+wt
            endif 

            v = map_value64(mp,ix_left,iy_upper,error) 
            if (ieee_is_finite(v)) then
                wt = wt_left * wt_upper
                tot = tot + (v * wt)
                tot_wt = tot_wt+wt
            endif

            v = map_value64(mp,ix_right,iy_upper,error) 
            if (ieee_is_finite(v)) then
                wt = wt_right*wt_upper
                tot = tot + (v * wt)
                tot_wt = tot_wt + wt
            endif
            
            v = map_value64(mp,ix_right,iy_lower,error) 
            if (ieee_is_finite(v)) then
                wt = wt_right*wt_lower 
                tot = tot + (v * wt)
                tot_wt = tot_wt+wt
            endif 

            value = tot/tot_wt

            return

        end subroutine get_2d_map_value64

        subroutine get_2d_map_value_no_interp32(mp,x,y,value,error)

            ! this subroutine returns the non interpolated value closest to a point x,y

            implicit none

            class(earth_ref_map_str32),intent(IN)    :: mp       ! map structure to be interpolated
            real(real32),intent(IN)                :: x        ! x, in x,y space
            real(real32),intent(IN)                :: y        ! y, in x,y space
            real(real32),intent(OUT)               :: value    ! the interpolated value
            integer(int32),intent(OUT)             :: error

            integer(int32)                         :: ix                    ! map index of point closest to the desired point
            integer(int32)                         :: iy                    ! map index of point closest to the desired point

            error = 0

            ! make sure x and y are finite and within bounds

            call check_x_y_2d_map(mp,x,y,error)
            if (error .ne. 0) then
                value = nan_f32
                return
            endif

            ! find closest x index
            select case (mp%x_axis%axis_type)
                case (CELL_CENTERED_2D_MAP) 
                    ix        = floor((x-mp%x_axis%min)/mp%x_axis%delta) + 1
                case (EDGE_CENTERED_2D_MAP) 
                    ix  = floor((x-mp%x_axis%min)/mp%x_axis%delta + 0.5) + 1
                case default
                    error = -1
                    value = nan_f32
                    return
            end select

            ! find closest y index

            select case(mp%y_axis%axis_type)

                case (CELL_CENTERED_2D_MAP)
                    iy    = floor((y-mp%y_axis%min)/mp%y_axis%delta) + 1
                case (EDGE_CENTERED_2D_MAP)
                    iy    = floor((y-mp%y_axis%min)/mp%y_axis%delta + 0.5) + 1
                case default
                    error = -1
                    value = nan_f32
                    return
            end select

            ! now get the data
            value = map_value32(mp,ix,iy,error)
            return

        end subroutine get_2d_map_value_no_interp32

        subroutine get_2d_map_value_no_interp64(mp,x,y,value,error)

            ! this subroutine returns the non interpolated value closest to a point x,y

            implicit none

            class(earth_ref_map_str64),intent(IN)    :: mp       ! map structure to be interpolated
            real(real32),intent(IN)                :: x        ! x, in x,y space
            real(real32),intent(IN)                :: y        ! y, in x,y space
            real(real64),intent(OUT)               :: value    ! the interpolated value
            integer(int32),intent(OUT)             :: error

            integer(int32)                         :: ix                    ! map index of point closest to the desired point
            integer(int32)                         :: iy                    ! map index of point closest to the desired point

            error = 0

            ! make sure x and y are finite and within bounds

            call check_x_y_2d_map(mp,x,y,error)
            if (error .ne. 0) then
                value = nan_f64
                return
            endif

            ! find closest x index
            select case (mp%x_axis%axis_type)
                case (CELL_CENTERED_2D_MAP) 
                    ix        = floor((x-mp%x_axis%min)/mp%x_axis%delta) + 1
                case (EDGE_CENTERED_2D_MAP) 
                    ix  = floor((x-mp%x_axis%min)/mp%x_axis%delta + 0.5) + 1
                case default
                    error = 1
                    value = nan_f64
                    return
            end select

            ! find closest y index

            select case(mp%y_axis%axis_type)

                case (CELL_CENTERED_2D_MAP)
                    iy    = floor((y-mp%y_axis%min)/mp%y_axis%delta) + 1
                case (EDGE_CENTERED_2D_MAP)
                    iy    = floor((y-mp%y_axis%min)/mp%y_axis%delta + 0.5) + 1
                case default
                    error = 1
                    value = nan_f64
                    return
            end select

            ! now get the data
            value = map_value64(mp,ix,iy,error)
            return

        end subroutine get_2d_map_value_no_interp64

        logical function is_inside_quad_map_base(mp,ix,iy,lat,lon,lats_quad,lons_quad)
            class(earth_ref_map_str),intent(IN)     :: mp
            integer(int32),intent(IN)    :: ix
            integer(int32),intent(IN)    :: iy
            real(real32),intent(IN)      :: lat
            real(real32),intent(IN)      :: lon
            real(real32),dimension(4),intent(OUT)    :: lons_quad
            real(real32),dimension(4),intent(OUT)    :: lats_quad
            integer(int32)               :: error


            call map_lat_lon_base(mp,ix,iy,lats_quad(1),lons_quad(1),error)
            if (error .ne. 0) then
                is_inside_quad_map_base = .false.
                return
            endif

            call map_lat_lon_base(mp,ix+1,iy,lats_quad(2),lons_quad(2),error)
            if (error .ne. 0) then
                is_inside_quad_map_base = .false.
                return
            endif

            call map_lat_lon_base(mp,ix+1,iy+1,lats_quad(3),lons_quad(3),error)
            if (error .ne. 0) then
                is_inside_quad_map_base = .false.
                return
            endif

            call map_lat_lon_base(mp,ix,iy+1,lats_quad(4),lons_quad(4),error)
            if (error .ne. 0) then
                is_inside_quad_map_base = .false.
                return
            endif

            is_inside_quad_map_base = is_inside_quadrilateral(lons_quad,lats_quad,lon,lat,.false.)

            if (.not. is_inside_quad_map_base) then
                lats_quad = nan_f32
                lons_quad = nan_f32
            endif

        end function is_inside_quad_map_base

        logical function is_inside_quad_map32(mp,ix,iy,lat,lon,lats_quad,lons_quad,data_quad)
            class(earth_ref_map_str32),intent(IN)     :: mp
            integer(int32),intent(IN)    :: ix
            integer(int32),intent(IN)    :: iy
            real(real32),intent(IN)      :: lat
            real(real32),intent(IN)      :: lon
            real(real32),dimension(4),intent(OUT)    :: lons_quad
            real(real32),dimension(4),intent(OUT)    :: lats_quad
            real(real32),dimension(4),intent(OUT),optional    :: data_quad

            integer(int32)               :: error,error_tot

            is_inside_quad_map32 = is_inside_quad_map_base(mp,ix,iy,lat,lon,lats_quad,lons_quad)
            
            if ((is_inside_quad_map32) .and. (present(data_quad))) then

                error_tot = 0

                data_quad(1) = map_value32(mp,ix,iy,error)
                error_tot = error_tot + error
                data_quad(2) = map_value32(mp,ix+1,iy,error)
                error_tot = error_tot + error
                data_quad(3) = map_value32(mp,ix+1,iy+1,error)
                error_tot = error_tot + error
                data_quad(4) = map_value32(mp,ix,iy+1,error)
                error_tot = error_tot + error

                if (error_tot .ne. 0) then
                    print *,"warning - error_tot not zero in is_inside_quad_map32"
                endif
            endif 

        end function is_inside_quad_map32

        logical function is_inside_quad_map64(mp,ix,iy,lat,lon,lats_quad,lons_quad,data_quad)
            class(earth_ref_map_str64),intent(IN)     :: mp
            integer(int32),intent(IN)    :: ix
            integer(int32),intent(IN)    :: iy
            real(real32),intent(IN)      :: lat
            real(real32),intent(In)      :: lon
            real(real32),dimension(4),intent(OUT)    :: lons_quad
            real(real32),dimension(4),intent(OUT)    :: lats_quad
            real(real64),dimension(4),intent(OUT),optional    :: data_quad

            integer(int32)               :: error,error_tot

           

            is_inside_quad_map64 = is_inside_quad_map_base(mp,ix,iy,lat,lon,lats_quad,lons_quad)
            if ((is_inside_quad_map64) .and. (present(data_quad))) then

                error_tot = 0

                data_quad(1) = map_value64(mp,ix,iy,error)
                error_tot = error_tot + error
                data_quad(2) = map_value64(mp,ix+1,iy,error)
                error_tot = error_tot + error
                data_quad(3) = map_value64(mp,ix+1,iy+1,error)
                error_tot = error_tot + error
                data_quad(4) = map_value64(mp,ix,iy+1,error)
                error_tot = error_tot + error

                if (error_tot .ne. 0) then
                    print *,"warning - error_tot not zero in is_inside_quad_map32"
                endif
            endif  
            
        end function is_inside_quad_map64


        ! real(real32) function map_value(mp,ix,iy,error)

        !     ! this subroutine returns the value of the map at xi,yi even
        !     ! if xi and yi are out of range returns NAN
        !     ! This will be extended when boundary conditions are included

        !     implicit none

        !     class(earth_ref_map_str),intent(IN)     :: mp
        !     integer(int32),intent(IN)    :: ix
        !     integer(int32),intent(IN)    :: iy
        !     integer(int32),intent(OUT)   :: error

        !     if ((ix .ge. 1) .and. &
        !         (ix .le. mp%num_pts_x) .and. &
        !         (iy .ge. 1) .and. &
        !         (iy .le. mp%num_pts_y)) then
        !         map_value = mp%dat(ix,iy)
        !     else
        !         map_value = nan_f32
        !     endif

        !     error = 0
        !     return
        ! end function map_value


        real(real32) function map_value32(mp,ix,iy,error)

            ! this subroutine returns the value of the map at xi,yi even
            ! if xi and yi are out of range by applying the appropriate
            ! boundary conditions

            implicit none

            class(earth_ref_map_str32),intent(IN)     :: mp
            integer(int32),intent(IN)    :: ix
            integer(int32),intent(IN)    :: iy
            integer(int32),intent(OUT)   :: error

            integer(int32)               :: ix_wrapped
            integer(int32)               :: iy_wrapped

            error = 0
            call wrap_indices_2d_map_base(mp,ix,iy,ix_wrapped,iy_wrapped,error) 

            if (error .eq. 0) then
                map_value32 = mp%dat(ix_wrapped,iy_wrapped)
            else
                map_value32 = nan_f32
            endif

        end function map_value32

        real(real64) function map_value64(mp,ix,iy,error)

            ! this subroutine returns the value of the map at xi,yi even
            ! if xi and yi are out of range by applying the appropriate
            ! boundary conditions

            implicit none

            class(earth_ref_map_str64),intent(IN)     :: mp
            integer(int32),intent(IN)    :: ix
            integer(int32),intent(IN)    :: iy
            integer(int32),intent(OUT)   :: error

            integer(int32)               :: ix_wrapped
            integer(int32)               :: iy_wrapped

            error = 0
            call wrap_indices_2d_map_base(mp,ix,iy,ix_wrapped,iy_wrapped,error) 

            if (error .eq. 0) then
                map_value64 = mp%dat(ix_wrapped,iy_wrapped)
            else
                map_value64 = nan_f64
            endif

        end function map_value64

        subroutine map_lat_lon_base(mp,ix,iy,lat,lon,error)
            class(earth_ref_map_str),intent(IN)   :: mp
            integer(int32),intent(IN)    :: ix
            integer(int32),intent(IN)    :: iy
            real(real32),intent(OUT)     :: lat
            real(real32),intent(OUT)     :: lon
            integer(int32),intent(OUT)   :: error

            integer(int32)   :: ix_wrapped
            integer(int32)   :: iy_wrapped


            error = 0
            call wrap_indices_2d_map_base(mp,ix,iy,ix_wrapped,iy_wrapped,error)

            if (error .eq. 0) then
                lat = mp%lats(ix_wrapped,iy_wrapped)
                print *,'lat: ',ix,iy,lat
                if (lat .le. -998) then
                    lat = nan_f32
                    error = 1
                endif
                lon = mp%lons(ix_wrapped,iy_wrapped)
                print *,'lon: ',ix,iy,lon
                if (lon .le. -998) then
                    lon = nan_f32
                    error = 1
                endif
            else
                lat = nan_f32
                lon = nan_f32
            endif

            return

        end subroutine map_lat_lon_base

        subroutine map_lat_lon32(mp,ix,iy,lat,lon,error)
            class(earth_ref_map_str32),intent(IN)   :: mp
            integer(int32),intent(IN)    :: ix
            integer(int32),intent(IN)    :: iy
            real(real32),intent(OUT)     :: lat
            real(real32),intent(OUT)     :: lon
            integer(int32),intent(OUT)   :: error

            call map_lat_lon_base(mp,ix,iy,lat,lon,error)
        end subroutine map_lat_lon32

        subroutine map_lat_lon64(mp,ix,iy,lat,lon,error)
            class(earth_ref_map_str64),intent(IN)   :: mp
            integer(int32),intent(IN)    :: ix
            integer(int32),intent(IN)    :: iy
            real(real32),intent(OUT)     :: lat
            real(real32),intent(OUT)     :: lon
            integer(int32),intent(OUT)   :: error

            call map_lat_lon_base(mp,ix,iy,lat,lon,error)
        end subroutine map_lat_lon64

        subroutine write_map32_netcdf(mp,nc_filename,fill_value,error) !,cf_var_metadata,cf_global_metadata)

            implicit none

            class(earth_ref_map_str32),intent(IN) :: mp
            character(*),intent(IN)              :: nc_filename
            real(real32),intent(IN)              :: fill_value
            integer(int32),intent(OUT)           :: error
            !type(CF_1_8_attr),intent(IN),optional :: cf_var_metadata
            !type(CF_1_8_attr),intent(IN),optional :: cf_global_metadata

            integer(int32)  :: ncid,varid,dim_lon,dim_lat,dim_x,dim_y,ix,iy,loc_error
            character(len = :),allocatable :: data_name

            real(real32),allocatable,dimension(:)  :: x_locations
            real(real32),allocatable,dimension(:)  :: y_locations
            real(real32)                            :: x_loc
            real(real32)                            :: y_loc

            allocate(x_locations(mp%x_axis%num_pts))
            allocate(y_locations(mp%y_axis%num_pts))

            ! do ix = 1,mp%x_axis%num_pts
            !     call mp%get_grid_location(ix,1,x_loc,y_loc,loc_error)
            !     x_locations(ix) = x_loc 
            ! enddo

            ! do iy = 1,mp%y_axis%num_pts
            !     call mp%get_grid_location(1,iy,x_loc,y_loc,loc_error)
            !     y = y_loc 
            ! enddo

            do ix = 1,mp%x_axis%num_pts
                call mp%get_grid_location(ix,1,x_loc,y_loc,loc_error)
                x_locations(ix) = x_loc 
            enddo

            do iy = 1,mp%y_axis%num_pts
                call mp%get_grid_location(1,iy,x_loc,y_loc,loc_error)
                y_locations(iy) = y_loc 
            enddo

            error = 0
        
            call handle_nc_err(nf90_create(nc_filename, ior(NF90_CLOBBER, NF90_NETCDF4), ncid))
        
            ! Define global attributes
            call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "Conventions", "CF-1.8,ACDD-1.3"))
            call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "institution", "REMSS"))
            call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "creator_name", "Remote Sensing Systems"))
            call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "creator_email", "support@remss.com"))
            call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "creator_url", "http://www.remss.com"))
            call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "geospatial_y_min", mp%y_axis%min))
            call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "geospatial_y_max", mp%y_axis%max))
            call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "geospatial_x_min", mp%x_axis%min))
            call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "geospatial_x_max", mp%x_axis%max))

            

            !add other global attributes, if passed
            ! if (present(cf_global_metadata)) then
            !     call put_CF_1_8_attr(cf_global_metadata,ncid)
            ! endif

            call handle_nc_err(nf90_def_dim(ncid, "x", mp%x_axis%num_pts, dim_x))
            call handle_nc_err(nf90_def_dim(ncid, "y", mp%y_axis%num_pts, dim_y))

            call handle_nc_err(nf90_def_var(ncid, "x", NF90_FLOAT,dim_x, varid))
            call handle_nc_err(nf90_put_var(ncid, varid, x_locations(:)))

            call handle_nc_err(nf90_def_var(ncid, "y", NF90_FLOAT, dim_y, varid))
            call handle_nc_err(nf90_put_var(ncid, varid, y_locations(:)))
        
            data_name = 'data'
            ! if (present(cf_var_metadata)) then
            !     if (len(cf_var_metadata%standard_name) > 0) then
            !         data_name = cf_var_metadata%standard_name
            !     endif
            ! endif

            call handle_nc_err(nf90_def_var(ncid, "latitude", NF90_FLOAT, [dim_x,dim_y], varid, &
                                deflate_level=2, shuffle=.true.))
            call handle_nc_err(nf90_put_var(ncid, varid, mp%lats))

            call handle_nc_err(nf90_def_var(ncid, "longitude", NF90_FLOAT, [dim_x,dim_y], varid, &
                                deflate_level=2, shuffle=.true.))
            call handle_nc_err(nf90_put_var(ncid, varid, mp%lons))

            call handle_nc_err(nf90_def_var(ncid, "Data", NF90_FLOAT, [dim_x,dim_y], varid, &
                                deflate_level=2, shuffle=.true.))
            call handle_nc_err(nf90_def_var_fill(ncid,varid,0,fill_value))
            call handle_nc_err(nf90_put_var(ncid, varid, mp%dat))

            !print *,ncid,varid
            ! if (present(cf_var_metadata)) then
            !     call put_CF_1_8_attr(cf_var_metadata,ncid, varid)
            ! endif

        call handle_nc_err(nf90_inq_varid(ncid, "latitude", varid))
        call handle_nc_err(nf90_put_att(ncid, varid, "standard_name", "latitude"))
        call handle_nc_err(nf90_put_att(ncid, varid, "units", "degrees_north"))

        call handle_nc_err(nf90_inq_varid(ncid, "longitude", varid))
        call handle_nc_err(nf90_put_att(ncid, varid, "standard_name", "longitude"))
        call handle_nc_err(nf90_put_att(ncid, varid, "units", "degrees_east"))

        call handle_nc_err(nf90_inq_varid(ncid, "Data", varid))
        
        call handle_nc_err(nf90_close(ncid))

    end subroutine write_map32_netcdf

    subroutine write_map64_netcdf(mp,nc_filename,fill_value,error) !,cf_var_metadata,cf_global_metadata)

            implicit none

            class(earth_ref_map_str64),intent(IN) :: mp
            character(*),intent(IN)               :: nc_filename
            real(real64),intent(IN)               :: fill_value
            integer(int32),intent(OUT)            :: error

            !type(CF_1_8_attr),intent(IN),optional :: cf_var_metadata
            !type(CF_1_8_attr),intent(IN),optional :: cf_global_metadata

            integer(int32)  :: ncid,varid,dim_lon,dim_lat,dim_x,dim_y,ix,iy,loc_error
            character(len = :),allocatable :: data_name

            real(real32),allocatable,dimension(:)   :: x_locations
            real(real32),allocatable,dimension(:)   :: y_locations
            real(real32)                            :: x_loc
            real(real32)                            :: y_loc
            
            allocate(x_locations(mp%x_axis%num_pts))
            allocate(y_locations(mp%y_axis%num_pts))

            do ix = 1,mp%x_axis%num_pts
                call mp%get_grid_location(ix,1,x_loc,y_loc,loc_error)
                x_locations(ix) = x_loc 
            enddo

            do iy = 1,mp%y_axis%num_pts
                call mp%get_grid_location(1,iy,x_loc,y_loc,loc_error)
                y_locations(iy) = y_loc 
            enddo

            error = 0
        
            call handle_nc_err(nf90_create(nc_filename, ior(NF90_CLOBBER, NF90_NETCDF4), ncid))
        
            ! Define global attributes
            call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "Conventions", "CF-1.8,ACDD-1.3"))
            call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "institution", "REMSS"))
            call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "creator_name", "Remote Sensing Systems"))
            call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "creator_email", "support@remss.com"))
            call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "creator_url", "http://www.remss.com"))
            call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "geospatial_y_min", mp%y_axis%min))
            call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "geospatial_y_max", mp%y_axis%max))
            call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "geospatial_x_min", mp%x_axis%min))
            call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "geospatial_x_max", mp%x_axis%max))

            !add other global attributes, if passed
            ! if (present(cf_global_metadata)) then
            !     call put_CF_1_8_attr(cf_global_metadata,ncid)
            ! endif

            call handle_nc_err(nf90_def_dim(ncid, "x", mp%x_axis%num_pts, dim_x))
            call handle_nc_err(nf90_def_dim(ncid, "y", mp%y_axis%num_pts, dim_y))

            call handle_nc_err(nf90_def_var(ncid, "x", NF90_FLOAT,dim_x, varid))
            call handle_nc_err(nf90_put_var(ncid, varid, x_locations(:)))

            call handle_nc_err(nf90_def_var(ncid, "y", NF90_FLOAT, dim_y, varid))
            call handle_nc_err(nf90_put_var(ncid, varid, y_locations(:)))
        
            data_name = 'data'
            ! if (present(cf_var_metadata)) then
            !     if (len(cf_var_metadata%standard_name) > 0) then
            !         data_name = cf_var_metadata%standard_name
            !     endif
            ! endif

            call handle_nc_err(nf90_def_var(ncid, "latitude", NF90_FLOAT, [dim_x,dim_y], varid, &
                                deflate_level=2, shuffle=.true.))
            call handle_nc_err(nf90_put_var(ncid, varid, mp%lats))

            call handle_nc_err(nf90_def_var(ncid, "longitude", NF90_FLOAT, [dim_x,dim_y], varid, &
                                deflate_level=2, shuffle=.true.))
            call handle_nc_err(nf90_put_var(ncid, varid, mp%lons))

            call handle_nc_err(nf90_def_var(ncid, "Data", NF90_DOUBLE, [dim_x,dim_y], varid, &
                                deflate_level=2, shuffle=.true.))
            call handle_nc_err(nf90_def_var_fill(ncid,varid,0,fill_value))
            
            call handle_nc_err(nf90_put_var(ncid, varid, mp%dat))

            
            ! if (present(cf_var_metadata)) then
            !     call put_CF_1_8_attr(cf_var_metadata,ncid, varid)
            ! endif

        call handle_nc_err(nf90_inq_varid(ncid, "latitude", varid))
        call handle_nc_err(nf90_put_att(ncid, varid, "standard_name", "latitude"))
        call handle_nc_err(nf90_put_att(ncid, varid, "units", "degrees_north"))

        call handle_nc_err(nf90_inq_varid(ncid, "longitude", varid))
        call handle_nc_err(nf90_put_att(ncid, varid, "standard_name", "longitude"))
        call handle_nc_err(nf90_put_att(ncid, varid, "units", "degrees_east"))

        call handle_nc_err(nf90_inq_varid(ncid, "Data", varid))
        call handle_nc_err(nf90_close(ncid))

    end subroutine write_map64_netcdf

    
end module earth_ref_maps

