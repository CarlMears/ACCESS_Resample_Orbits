module ecef_conversions

    ! This module contains subroutines for converting between
    ! latitude, longitude and ECEF XYZ coordinates at the Earth's
    ! surface.  The subroutines are based on the equations given
    ! in the Wikipedia article on Geodetic System.
    !
    ! note the for the ellipsoid, the conversion from ECEF x,y,z
    ! to latitude, longitude is not exact.  The iterative method
    ! used here is based on the algorithm given in the Wikipedia.  Is
    ! should converge within a few iterations.
    
    use, intrinsic :: iso_fortran_env, only: int8, int32, real32, real64
    use trig_degrees, only: dsind, dcosd, dasind, dacosd, datand, datan2d

    implicit none

    !WGS84 ellipsoid constants
    real(real64),parameter :: a = 6378137.0
    real(real64),parameter :: b = 6356752.314245
    real(real64),parameter :: f = (a - b) / a
    real(real64),parameter :: e_sq = f * (2.0_real64 - f)

    interface lat_lon_to_ecef_xyz
        module procedure lat_lon_to_ecef_xyz_64, lat_lon_to_ecef_xyz_32
    end interface lat_lon_to_ecef_xyz

    interface XYZtoLatLon
        module procedure XYZtoLatLon_64, XYZtoLatLon_32
    end interface XYZtoLatLon

contains

    subroutine lat_lon_to_ecef_xyz_64(latitude, longitude, x, y, z)

        implicit none

        real(real64), intent(in) :: latitude, longitude
        real(real64), intent(out) :: x, y, z

        real(real64)           :: n, cos_lat, sin_lat, cos_lon, sin_lon

        ! Calculate the geocentric radius
        cos_lat = cosd(latitude)
        sin_lat = sind(latitude)
        cos_lon = cosd(longitude)
        sin_lon = sind(longitude)

        n = a / sqrt(1 - e_sq * sin_lat * sin_lat)
        x = n * cos_lat * cos_lon
        y = n * cos_lat * sin_lon
        z = n * (1 - e_sq) * sin_lat

    end subroutine lat_lon_to_ecef_xyz_64

    subroutine lat_lon_to_ecef_xyz_32(latitude, longitude, x, y, z)

        implicit none

        real(real32), intent(in) :: latitude, longitude
        real(real32), intent(out) :: x, y, z

        real(real64) :: latitude64, longitude64
        real(real64) :: x64, y64, z64

        !WGS84 ellipsoid constants
        real(real64),parameter :: a = 6378137.0
        real(real64),parameter :: b = 6356752.314245
        real(real64),parameter :: f = (a - b) / a

        real(real64)           :: e_sq, n, cos_lat, sin_lat, cos_lon, sin_lon

        latitude64 = real(latitude, real64)
        longitude64 = real(longitude, real64)

        call lat_lon_to_ecef_xyz_64(latitude64, longitude64, x64, y64, z64)

        x = real(x64, real32)
        y = real(y64, real32)
        z = real(z64, real32)

    end subroutine lat_lon_to_ecef_xyz_32

    subroutine XYZtoLatLon_64(x, y, z, latitude, longitude)

            implicit none
            real(real64), intent(in) :: x, y, z
            real(real64), intent(out) :: latitude, longitude

            ! WGS84 ellipsoid constants
            real(real64) :: p, theta, lat, lon, prev_lat, N
            real(real64) :: threshold

            integer(int32),parameter :: max_iterations = 10
            integer(int32) :: i

            logical :: warn

            ! Set the WGS84 ellipsoid constants
            ! a = 6378137.0_real64    ! semi-major axis
            ! b = 6356752.314245_real64   ! semi-minor axis
            ! f = (a - b) / a    ! flattening
            ! e_sq = 2.0_real64 * f - f**2   ! squared eccentricity

            !longitude is exactly known
            longitude = atan2(y, x)

            ! find the latitude using the iterative method
            
            p = sqrt(x**2 + y**2)
            latitude = atan2(z, p)

            ! Iteratively refine latitude
            threshold = 1.0e-12_real64
            warn = .true.
            do i = 1, max_iterations
                N = a / sqrt(1.0_real64 - e_sq * sin(latitude)**2)
                prev_lat = latitude
                latitude = atan2(z + e_sq * N * sin(latitude), p)
                if (abs(latitude - prev_lat) < threshold) then
                    warn = .false.
                    exit
                endif
            end do

            ! Convert the latitude and longitude to degrees
            latitude = latitude * 180.0_real64 / acos(0.0_real64)
            longitude = longitude * 180.0_real64 / acos(0.0_real64)
            if (warn) then
                write(*,*) 'Warning: XYZtoLatLon did not converge'
            endif
    end subroutine XYZtoLatLon_64

    subroutine XYZtoLatLon_32(x, y, z, latitude, longitude)

            implicit none
            real(real32), intent(in) :: x, y, z
            real(real32), intent(out) :: latitude, longitude

            real(real64) :: x64, y64, z64
            real(real64) :: latitude64, longitude64

            x64 = real(x, real64)
            y64 = real(y, real64)
            z64 = real(z, real64)

            call XYZtoLatLon_64(x64, y64, z64, latitude64, longitude64)

            latitude = real(latitude64, real32)
            longitude = real(longitude64, real32)

    end subroutine XYZtoLatLon_32
end module ecef_conversions


    

