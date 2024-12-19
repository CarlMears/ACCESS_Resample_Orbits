! Simple netCDF helper routines
module io_nc
  use, intrinsic :: iso_fortran_env, only: int8, int32, real32, real64, ERROR_UNIT
  use netcdf
  use time_conversions, only: time_2000_to_ymd_secs, ccsds_utc_to_time_2000
  use time_conversions_extra, only: acc_leap_secs_2000, ymd_parts_to_2000
  implicit none
  private
  public :: handle_nc_err, minmax_iso8601, d2000_iso8601, iso8601_d2000

contains

  ! -------------------------------------------------------------------------------
  ! Find min/max time strings
  subroutine minmax_iso8601(time_2000, time_start, time_end)
    real(real64), dimension(:), intent(in) :: time_2000
    character(len=20), intent(out) :: time_start, time_end
    time_start = d2000_iso8601(minval(time_2000))
    time_end = d2000_iso8601(maxval(time_2000))
  end subroutine minmax_iso8601

  ! -------------------------------------------------------------------------------
  ! Convert time value to an ISO 8601 string
  !
  ! The time value is the number of seconds since
  ! 2000-01-01T00:00:00Z. Note that leap seconds are not included in
  ! the input so must be accounted for in here.
  pure function d2000_iso8601(time_val)
    real(real64), intent(in) :: time_val
    character(len=20) :: d2000_iso8601

    integer(int32) :: year, month, day, hour, min, sec
    real(real64) :: sec_of_year, sec_of_day
    integer :: leap_seconds

    ! Adjust for leap seconds by removing any that were added between
    ! 2000 and the current time
    leap_seconds = acc_leap_secs_2000(time_val)

    ! Obtain year/month/day and the number of seconds into the day
    call time_2000_to_ymd_secs(time_val - real(leap_seconds, real64), &
         year, month, day, sec_of_year, sec_of_day)

    ! Obtain hour/minute/second
    hour = floor(sec_of_day / 3600)
    sec_of_day = sec_of_day - hour * 3600
    min = floor(sec_of_day / 60)
    sec_of_day = sec_of_day - min * 60
    sec = floor(sec_of_day)

    write(d2000_iso8601, '(I4, "-", I2.2, "-", I2.2, "T", I2.2, ":", I2.2, ":", I2.2, "Z")') &
         year, month, day, hour, min, sec

  end function d2000_iso8601

  ! -------------------------------------------------------------------------------
  ! Convert an ISO 8601 string to a time value
  !
  ! The time value is the number of seconds since
  ! 2000-01-01T00:00:00Z. Leap seconds are accounted for.
  pure function iso8601_d2000(time_str)
    character(len=20), intent(in) :: time_str
    real(real64) :: iso8601_d2000

    integer(int32) :: year, month, day, hour, min, sec

    ! This assumes the string is in this format:
    ! YYYY-MM-DDTHH:MM:SSZ
    read(time_str(1:4), '(i4)') year
    read(time_str(6:7), '(i2)') month
    read(time_str(9:10), '(i2)') day
    read(time_str(12:13), '(i2)') hour
    read(time_str(15:16), '(i2)') min
    read(time_str(18:19), '(i2)') sec

    iso8601_d2000 = ymd_parts_to_2000(year, month, day, hour, min, sec)
  end function iso8601_d2000

  ! -------------------------------------------------------------------------------
  ! This is a very simple netCDF error handler: if there's an error,
  ! then display it and exit with a nonzero code.
  subroutine handle_nc_err(status)
    integer, intent(in) :: status
    if (status /= NF90_NOERR) then
       write (ERROR_UNIT, *) "ERROR: ", nf90_strerror(status)
       error stop "NetCDF error"
    end if
  end subroutine handle_nc_err
end module io_nc
