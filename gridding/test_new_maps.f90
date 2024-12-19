program test_new_maps

use, intrinsic :: iso_fortran_env, only: int8, int32, real32, real64
use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_finite
use NAN_support, only: nan_f32, nan_f64
use map_type_constants, only: CELL_CENTERED_2D_MAP,EDGE_CENTERED_2D_MAP,SMALL_NUMBER

use earth_ref_maps, only: earth_ref_map_str32
use earth_ref_maps, only: earth_ref_map_str64

use InterpBilinearQuadrilateral, only: init_InterpBilinearQuadrilateral
use InterpBilinearQuadrilateral, only: eval_InterpBilinearQuadrilateral
use InterpBilinearQuadrilateral, only: quad_interpolator_type

use polar_grids, only: init_polar_grid32
use ease2_calcs, only: ease2

implicit none

type(earth_ref_map_str32) :: mp32
type(earth_ref_map_str64) :: mp64
type(earth_ref_map_str32) :: polar_map32

type(quad_interpolator_type) :: qi

real(real32) :: min_x,max_x,min_y,max_y
integer(int32) :: num_x,type_x,num_y,type_y,error

real(real32),allocatable :: lat_array(:,:)
real(real32),allocatable :: lon_array(:,:)
real(real32),allocatable :: test_data32(:,:)
!real(real64),allocatable :: test_data64(:,:)

integer(int32) :: i,j,ix,iy
real(real32) :: value32,lat,lon,x_km,y_km
logical :: inside
real(real32),dimension(4)    :: lons_quad
real(real32),dimension(4)    :: lats_quad
real(real32),dimension(4)    :: data_quad

character(:),allocatable   :: pole
character(:),allocatable   :: resolution
character(:),allocatable   :: grid_type



min_x = 0.0
max_x = 20.0
min_y = 0.0
max_y = 15.0

num_x = 21.0
num_y = 16.0

type_x = EDGE_CENTERED_2D_MAP
type_y = EDGE_CENTERED_2D_MAP

call mp32%init(min_x,max_x,num_x,type_x,min_y,max_y,num_y,type_y,error)
call mp64%init(min_x,max_x,num_x,type_x,min_y,max_y,num_y,type_y,error)

!call set_up_map(mp,min_x,max_x,num_x,type_x,min_y,max_y,num_y,type_y,error)

allocate(lat_array(num_x,num_y))
allocate(lon_array(num_x,num_y))
allocate(test_data32(num_x,num_y))

do i = 1,num_x
    do j = 1,num_y
        test_data32(i,j) = 1.0*i
        lat_array(i,j) = 1.0*j
        lon_array(i,j) = 1.0*i
    enddo
enddo

call mp32%define_lats_lons(lat_array,lon_array,error)
call mp32%define_data(test_data32,error)

print *, mp32%dat(5,5)

!this should be the same..
call mp32%get_interp_value(4.0,4.0,value32,error)
print *,value32
call mp32%get_interp_value(4.3,4.3,value32,error)
print *,value32
call mp32%get_closest_value(4.3,4.3,value32,error)
print *,value32

call mp32%get_lat_lon(4,6,lat,lon,error)
print *,lat,lon
lat = 15.2
lon = 4.2
ix=4
iy=15
inside = mp32%is_inside_quad(ix,iy,lat,lon,lats_quad,lons_quad)
if (inside) then
    print '(f5.2,", ",f5.2," is inside the quadrilateral at ",i3,", ",i3)',lat,lon,ix,iy
else
    print '(f5.2,", ",f5.2," is NOT inside the quadrilateral at ",i3,", ",i3)',lat,lon,ix,iy
endif
print *,lats_quad
print *,lons_quad
lat = 6.8
lon = 4.2
ix=4
iy=6
inside = mp32%is_inside_quad(ix,iy,lat,lon,lats_quad,lons_quad,data_quad)
if (inside) then
    print '(f5.2,", ",f5.2," is inside the quadrilateral at ",i3,", ",i3)',lat,lon,ix,iy
else
    print '(f5.2,", ",f5.2," is NOT inside the quadrilateral at ",i3,", ",i3)',lat,lon,ix,iy
endif
print *,lats_quad
print *,lons_quad
print *,data_quad


call init_InterpBilinearQuadrilateral(     &
                 lons_quad,   & ! array of 4 x locations
                 lats_quad,   & ! array of 4 y locations
                 data_quad,   & ! array of 4 z locations
                 qi,          & ! the interpolator instance
                 error)         ! error flag  
                 
value32 = eval_InterpBilinearQuadrilateral(lon,lat,qi,error)
print *,'Interpolated Value ',value32

call mp32%deinit(error)

print *,error

pole = 'north'
resolution = '25km'
grid_type = 'ease2'

deallocate(lat_array)
deallocate(lon_array)
call init_polar_grid32( &
                        pole,        &
                        resolution,  &
                        grid_type,   &
                        polar_map32,   &
                        error)

do iy = 360,362
    do ix = 360,362
        call polar_map32%get_grid_location(ix,iy,x_km,y_km,error)
        print *,ix,iy,x_km,y_km,error
        print *,ix,iy,polar_map32%lats(ix,721-iy),polar_map32%lons(ix,721-iy)
        call ease2(polar_map32%lats(ix,721-iy),polar_map32%lons(ix,721-iy),x_km,y_km)
        print *,ix,iy,x_km,y_km
        print *
    enddo
enddo



lat = 89.8417272
lon = -45.0
print *,'testing EASE2 transform and index-finding'
call ease2(lat,lon,x_km,y_km)

call polar_map32%get_grid_index(x_km,y_km,ix,iy,error)

print *,'Error: ',error
print *,lat,lon,x_km,y_km,ix,iy

lat = 89.8417272
lon = 45.0
call ease2(lat,lon,x_km,y_km)

call polar_map32%get_grid_index(x_km,y_km,ix,iy,error)

print *,'Error: ',error
print *,lat,lon,x_km,y_km,ix,iy

lat = 89.8417272
lon = 135.0
call ease2(lat,lon,x_km,y_km)

call polar_map32%get_grid_index(x_km,y_km,ix,iy,error)

print *,'Error: ',error
print *,lat,lon,x_km,y_km,ix,iy

lat = 89.8417272
lon = -135.0
call ease2(lat,lon,x_km,y_km)

call polar_map32%get_grid_index(x_km,y_km,ix,iy,error)

print *,'Error: ',error
print *,lat,lon,x_km,y_km,ix,iy

do ix = 359,362
    do iy = 359,362
        inside = polar_map32%is_inside_quad(ix,iy,lat,lon,lats_quad,lons_quad,data_quad)
        if (inside) then
            print '(f7.2,", ",f7.2," is inside the quadrilateral at ",i3,", ",i3)',lat,lon,ix,iy
        else
            print '(f7.2,", ",f7.2," is NOT inside the quadrilateral at ",i3,", ",i3)',lat,lon,ix,iy
        endif
    enddo
enddo


end program
