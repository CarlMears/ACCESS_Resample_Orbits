  ! Forward-project lat/lon into EASE-2 
  ! lat/lon are given in degrees. Returns x/y in units of km from the
  ! projection origin, which is the south pole. x and y range from +- 9000 km.

module ease2_calcs
    use,intrinsic :: iso_fortran_env, only: real32,real64,int32
    use trig_degrees, only: dcosd,dsind

    interface ease2
        module procedure ease2_32
        module procedure ease2_64
    end interface

    private
    public ease2

contains
    subroutine ease2_32(lat,lon, x, y)
        real(real32), intent(in) :: lat, lon
        real(real32), intent(out) :: x, y
        
        ! from WGS84
        real(real64), parameter :: a = 6378137.d0 ! Equatorial radius in meters
        real(real64), parameter :: q_theta_90= 0.1995531087D+01  !q_theta for 90 deg
        real(real64), parameter :: ecc = 0.0818191908426

        integer(int32) isign
        real(real64) rho, radicand
        real(real64) q_theta, sin_lat, e_sin_lat
        
        if(lat.ge.0) then
            isign= -1
        else 
            isign= 1
            ! print *,'Need to check sign for the SH'
        endif
        
        sin_lat = dsind(real(lat,real64))
        e_sin_lat = ecc * sin_lat

        q_theta = (1 - ecc**2) * &
                  (sin_lat / (1 - e_sin_lat**2) - &
                   1 / (2 * ecc) * &
                   log((1 - e_sin_lat)/(1 + e_sin_lat)))

        radicand = q_theta_90 + isign*q_theta
        if (abs(radicand) < 1e-8) then
            rho = 0.
        else
            rho = a * sqrt(radicand)
        endif
        x = real(        rho * dsind(real(lon,real64)) * 1.d-3,real32)
        y = real(isign * rho * dcosd(real(lon,real64)) * 1.d-3,real32)

    end subroutine ease2_32

    subroutine ease2_64(lat,lon, x, y)
        real(real64), intent(in) :: lat, lon
        real(real64), intent(out) :: x, y
        
        ! from WGS84
        real(real64), parameter :: a = 6378137.d0 ! Equatorial radius in meters
        real(real64), parameter :: q_theta_90= 0.1995531087D+01  !q_theta for 90 deg
        real(real64), parameter :: ecc = 0.0818191908426

        integer(int32) isign
        real(real64) rho, radicand
        real(real64) q_theta, sin_lat, e_sin_lat
        
        if(lat.ge.0) then
            isign = -1
        else 
            isign = 1
        ! print *,'Need to check sign for the SH'
        endif
        
        sin_lat = dsind(lat)
        e_sin_lat = ecc * sin_lat

        q_theta = (1 - ecc**2) * &
                  (sin_lat / (1 - e_sin_lat**2) - &
                   1 / (2 * ecc) * &
                   log((1 - e_sin_lat)/(1 + e_sin_lat)))

        radicand = q_theta_90 + isign*q_theta
        if (abs(radicand) < 1e-8) then
            rho = 0.
        else
            rho = a * sqrt(radicand)
        endif
        x =         rho * dsind(lon) * 1.d-3
        y = isign * rho * dcosd(lon) * 1.d-3

    end subroutine ease2_64
end module ease2_calcs

