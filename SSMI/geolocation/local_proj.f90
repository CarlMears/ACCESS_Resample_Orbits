module local_proj

    use, intrinsic :: iso_fortran_env, only: int8, int32, real32, real64
    use trig_degrees, only: DEG2RAD_F32,DEG2RAD_F64, RAD2DEG_F32, RAD2DEG_F64
    
    !WGS84 ellipsoid constants
    real(real64),parameter :: a = 6378137.0
    real(real64),parameter :: b = 6356752.314245
    real(real64),parameter :: f = (a - b) / a
    real(real64),parameter :: e_sq = f * (2.0_real64 - f)

    interface lat_lon_to_local_xy
        module procedure lat_lon_to_local_xy_real64, lat_lon_to_local_xy_real32
    end interface lat_lon_to_local_xy

    interface local_xy_to_lat_lon
        module procedure local_xy_to_lat_lon_real64, local_xy_to_lat_lon_real32
    end interface local_xy_to_lat_lon

contains

    subroutine lat_lon_to_local_xy_real64(lat,lon,lat0,lon0,x,y)

        implicit none

        real(real64),intent(in) :: lat, lon, lat0, lon0
        real(real64),intent(out) :: x, y

        real(real64) :: lat_rad, lon_rad, lat0_rad, lon0_rad
        real(real64) :: delta_lon, delta_lat
        real(real64) :: rho_EW, rho_NS
      
        ! Convert coordinates to radians
        lat_rad = lat * DEG2RAD_F64
        lon_rad = lon * DEG2RAD_F64
        lat0_rad = lat0 * DEG2RAD_F64
        lon0_rad = lon0 * DEG2RAD_F64

        ! Calculate difference in longitude
        delta_lon = lon_rad - lon0_rad
        delta_lat = lat_rad - lat0_rad
     
        ! Calculate auxiliary values
    
        rho_EW = a / sqrt(1 - e_sq * sin(lat0_rad) ** 2)
        rho_NS = a * (1 - e_sq) / ((1 - e_sq * sin(lat0_rad) ** 2) ** (3./2.))
    
        ! Calculate x-coordinate
        x = rho_EW * cos(lat0_rad) * sin(delta_lon)

        ! Calculate y-coordinate
        y = rho_NS * sin(delta_lat)

    end subroutine lat_lon_to_local_xy_real64

    subroutine lat_lon_to_local_xy_real32(lat,lon,lat0,lon0,x,y)

        implicit none

        real(real32),intent(in) :: lat, lon, lat0, lon0
        real(real32),intent(out) :: x, y

        real(real64) :: lat64, lon64, lat064, lon064
        real(real64) :: x64, y64

        lat64 = real(lat,real64)
        lon64 = real(lon,real64)
        lat064 = real(lat0,real64)
        lon064 = real(lon0,real64)

        call lat_lon_to_local_xy_real64(lat64,lon64,lat064,lon064,x64,y64)

        x = real(x64,real32)
        y = real(y64,real32)

    end subroutine lat_lon_to_local_xy_real32

    subroutine local_xy_to_lat_lon_real64(x,y,lat0,lon0,lat,lon)

        real(real64),intent(in) :: x, y, lat0, lon0
        real(real64),intent(out) :: lat, lon

        real(real64) :: lat_rad, lon_rad, lat0_rad, lon0_rad
        real(real64) :: rho_EW, rho_NS

        ! Convert reference latitude to radians
        lat0_rad = lat0*DEG2RAD_F64
        lon0_rad = lon0*DEG2RAD_F64

        ! Calculate auxiliary values
    
        rho_EW = a / sqrt(1 - e_sq * sin(lat0_rad) ** 2)
        rho_NS = a * (1 - e_sq) / ((1 - e_sq * sin(lat0_rad) ** 2) ** (3/2))

        ! Calculate latitude
        lat_rad = lat0_rad + asin((y / rho_NS))
        lat = lat_rad*RAD2DEG_F64

        ! Calculate longitude
        lon_rad = lon0_rad + asin(x / ((rho_EW)*cos(lat0_rad)))
        lon = lon_rad*RAD2DEG_F64
        lon = modulo(lon + 1080.0, 360.0)

    end subroutine local_xy_to_lat_lon_real64

    subroutine local_xy_to_lat_lon_real32(x,y,lat0,lon0,lat,lon)

        real(real32),intent(in) :: x, y, lat0, lon0
        real(real32),intent(out) :: lat, lon

        real(real64) :: x64, y64, lat064, lon064
        real(real64) :: lat64, lon64

        x64 = real(x,real64)
        y64 = real(y,real64)
        lat064 = real(lat0,real64)
        lon064 = real(lon0,real64)

        call local_xy_to_lat_lon_real64(x64,y64,lat064,lon064,lat64,lon64)

        lat = real(lat64,real32)
        lon = real(lon64,real32)

    end subroutine local_xy_to_lat_lon_real32
end module local_proj