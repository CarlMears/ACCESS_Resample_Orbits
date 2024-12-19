module polar_grids
    use, intrinsic :: iso_fortran_env, only: int8, int32, real32, real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_finite
    use netcdf
    use io_nc, only: handle_nc_err
    use map_type_constants, only: CELL_CENTERED_2D_MAP
    use NSIDC_grids, only: read_NSIDC_polargrid
    use earth_ref_maps, only: earth_ref_map_str32,earth_ref_map_str64
    implicit none

    private
    public init_polar_grid32
    public init_polar_grid64

contains
! The following initializes the Earth Fixed grid using the NSIDC ease2 grids
! These are cell-centered grids
    subroutine init_polar_grid32(pole,        &
                               resolution,  &
                               grid_type,   &
                               polar_map32, &
                               error)

        character(*), intent(in)  :: pole
        character(*), intent(in)  :: resolution
        character(*), intent(in)  :: grid_type
        type(earth_ref_map_str32),intent(inout) :: polar_map32
        integer(int32),intent(OUT)  :: error

        integer(int32) :: num_x
        integer(int32) :: num_y

        real(real64),allocatable :: lat_array(:,:)
        real(real64),allocatable :: lon_array(:,:)

        real(real32),allocatable :: lat_array32(:,:)
        real(real32),allocatable :: lon_array32(:,:)


        real(real32) :: min_x,max_x,min_y,max_y

        call read_NSIDC_polargrid(pole,resolution,grid_type,num_x,num_y,lat_array,lon_array,error)

        allocate(lat_array32(num_x,num_y))
        allocate(lon_array32(num_x,num_y))

        lat_array32 = real(lat_array,real32)
        lon_array32 = real(lon_array,real32)
        

        select case(resolution) 
            case('25km')
                min_x = -8987.5
                min_y = -min_x
                max_x = -min_x
                max_y = -min_y
            case('12.5km')
                min_x = -8993.75
                min_y = -min_x
                max_x = -min_x
                max_y = -min_y
            case default
                error stop
            end select

        call polar_map32%init(min_x,max_x,num_x,CELL_CENTERED_2D_MAP, &
                              min_y,max_y,num_y,CELL_CENTERED_2D_MAP,error)
            
        ! print *,size(lat_array)
        ! print *,size(lon_array)
        ! print *,num_x,num_y
        ! print *,lat_array(1,1)
        ! print *,real(lat_array(1,1),real32)

        ! print *,size(polar_map32%lats)
        ! print *,size(polar_map32%lons)
        

        call polar_map32%define_lats_lons(lat_array32,lon_array32,error)

    end subroutine init_polar_grid32

    subroutine init_polar_grid64(pole,        &
                               resolution,  &
                               grid_type,   &
                               polar_map64, &
                               error)

        character(*), intent(in)  :: pole
        character(*), intent(in)  :: resolution
        character(*), intent(in)  :: grid_type
        type(earth_ref_map_str64),intent(inout) :: polar_map64
        integer(int32),intent(OUT)  :: error

        integer(int32) :: num_x
        integer(int32) :: num_y

        real(real64),allocatable :: lat_array(:,:)
        real(real64),allocatable :: lon_array(:,:)

        real(real32),allocatable :: lat_array32(:,:)
        real(real32),allocatable :: lon_array32(:,:)

        real(real32) :: min_x,max_x,min_y,max_y

        call read_NSIDC_polargrid(pole,resolution,grid_type,num_x,num_y,lat_array,lon_array,error)

        allocate(lat_array32(num_x,num_y))
        allocate(lon_array32(num_x,num_y))

        lat_array32 = real(lat_array,real32)
        lon_array32 = real(lon_array,real32)

        select case(resolution) 
            case('25km')
                min_x = -8987.5
                min_y = -min_x
                max_x = -min_x
                max_y = -min_y
            case('12.5km')
                min_x = -8993.75
                min_y = -min_x
                max_x = -min_x
                max_y = -min_y
            case default
                error stop
            end select

        call polar_map64%init(min_x,max_x,num_x,CELL_CENTERED_2D_MAP, &
                              min_y,max_y,num_y,CELL_CENTERED_2D_MAP,error)
            
        call polar_map64%define_lats_lons(lat_array32,lon_array32,error)

    end subroutine init_polar_grid64

    
end module polar_grids

