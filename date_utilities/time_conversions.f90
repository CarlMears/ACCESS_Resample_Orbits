! ---------------------------------------------------------------------------
!
! Ball Aerospace & Technologies Corp.
! ---------------------------------------------------------------------------
! WSF-M PROJECT: MWSDPS - MicroWave Sensor Data Processing Software
! ---------------------------------------------------------------------------
! This is an unpublished work, 2020 Ball Aerospace & Technologies
! Corp. and Remote Sensing Systems Inc. All rights are reserved,
! subject to the U.S. Government license described herein. This work
! was developed in whole or in part with U.S. Government sponsorship
! and was delivered to the Government pursuant to DFARS
! 252.227-13(b)(1) with unlimited rights, including the rights to
! reproduce, prepare derivative works, distribute copies to the
! public, and perform publicly and display publicly, by or on behalf
! of the Government.
! ---------------------------------------------------------------------------
! Government Purpose Rights Notice Contract Number: FA8810-18-C-0002
! Contractor Name: Ball Aerospace & Technologies Corp. Contractor
! Address: 1600 Commerce Street, Boulder, Colorado 80301
! ---------------------------------------------------------------------------
! The Government's rights to use, modify, reproduce, release, perform,
! display, or disclose these technical data are restricted by
! paragraph (b)(2) of the Rights in Technical Data Non-commercial
! Items clause contained in the above identified contract. Any
! reproduction of technical data or portions thereof marked with this
! legend must also reproduce the markings. DISTRIBUTION STATEMENT D:
! Distribution authorized to the Department of Defense and U.S. DoD
! contractors only. Reason: Critical Technology, Export Controlled.
! Date of determination: 15 Nov 2017. Other requests shall be referred
! to SMC/RSK, El Segundo, CA, 90245-2808. WARNING - This document
! contains technical data whose export is restricted by the Arms
! Export Control Act (Title 22, U.S.C., Sec 2751, et seq.) or the
! Export Administration Act of 1979 (Title 50, U.S.C., App. 2401 et
! seq.), as amended. Violations of these export laws are subject to
! severe criminal penalties. Disseminate in accordance with provisions
! of DoD Directive 5230.25.
!
! ---------------------------------------------------------------------------
!
module time_conversions
  use, intrinsic :: iso_fortran_env, only: int32, real32, real64
  use, intrinsic :: iso_c_binding, only: c_double, c_int
  implicit none
  private
  public :: date_parts_to_2000, time_2000_to_frac_year, &
       ccsds_utc_to_j2000_tt, ccsds_utc_to_j2000_ut1, ccsds_utc_to_time_2000, &
       time_2000_to_yo_secs, yo_to_ymd, time_2000_to_ymd_secs, is_leap, &
       SEC_PER_DAY, TAI_TO_TT

  ! The offset between TAI and TT in seconds
  real(real64), parameter :: TAI_TO_TT = 32.184

  ! Number of SI seconds per day
  real(real64), parameter :: SEC_PER_DAY = 86400

  ! This is the difference in seconds between 1958-01-01 00:00:00
  ! UTC/TAI and 2000-01-01 12:00:00 TAI *without* leap seconds included.
  ! It is 42 years and 12 hours. Of the 42 years, 10 of
  ! them are leap years.
  real(real64), parameter :: CCSDS_TO_J2000 = 15340.5 * SEC_PER_DAY

contains

  ! Convert time from seconds since 2000-01-01 00:00:00 to year and
  ! ordinal date, and also the number of seconds into the year and the
  ! number of seconds into the day.
  !
  ! The ordinal date, or day of year, ranges from 1 to 366. The
  ! seconds of year ranges from 0 to 366 * 86400. The seconds of the
  ! day ranges from 0 to 86400.
  !
  ! There is no support for a negative input time (i.e., before the
  ! year 2000).
  pure subroutine time_2000_to_yo_secs(time_2000, year, day_of_year, seconds_of_year, seconds_of_day)
    real(real64), intent(in) :: time_2000
    integer, intent(out) :: year, day_of_year
    real(real64), intent(out) :: seconds_of_year, seconds_of_day

    integer, parameter :: DAYS_PER_YEAR = 365
    integer, parameter :: DAYS_PER_LEAP_YEAR = 366

    integer :: days_to_year_start, days_this_year

    ! No support for inputs before 2000, so exit early if that's the case
    if (time_2000 < 0.) then
       year = 0
       day_of_year = 0
       seconds_of_year = 0.
       seconds_of_day = 0.
       return
    end if

    ! Find the year and day of the year starting at 2000, which is a
    ! leap year
    year = 2000
    days_this_year = DAYS_PER_LEAP_YEAR
    day_of_year = 1 + floor(time_2000 / SEC_PER_DAY)
    days_to_year_start = 0
    do while (day_of_year > days_this_year)
       year = year + 1
       day_of_year = day_of_year - days_this_year
       days_to_year_start = days_to_year_start + days_this_year

       ! This is for the next loop iteration
       if (is_leap(year)) then
          days_this_year = DAYS_PER_LEAP_YEAR
       else
          days_this_year = DAYS_PER_YEAR
       end if
    end do

    ! Then fractional seconds of the year and of the day
    seconds_of_year = (time_2000 - days_to_year_start * SEC_PER_DAY)
    seconds_of_day = (time_2000 - (days_to_year_start + (day_of_year - 1)) * SEC_PER_DAY)
  end subroutine time_2000_to_yo_secs

  ! Convert time from seconds since 2000-01-01 00:00:00 to year,
  ! month, and day. Also, the number of seconds into the year and the
  ! number of seconds into the day.
  !
  ! The seconds of year ranges from 0 to 366 * 86400. The seconds of
  ! the day ranges from 0 to 86400.
  !
  ! There is no support for a negative input time (i.e., before the
  ! year 2000).
  pure subroutine time_2000_to_ymd_secs(time_2000, year, month, day, seconds_of_year, seconds_of_day)
    real(real64), intent(in) :: time_2000
    integer, intent(out) :: year, month, day
    real(real64), intent(out) :: seconds_of_year, seconds_of_day

    integer :: day_of_year

    call time_2000_to_yo_secs(time_2000, year, day_of_year, seconds_of_year, seconds_of_day)
    call yo_to_ymd(year, day_of_year, month, day)
  end subroutine time_2000_to_ymd_secs

  ! Convert an ordinal date to month and day
  pure subroutine yo_to_ymd(year, day_of_year, month, day_of_month)
    integer, intent(in) :: year, day_of_year
    integer, intent(out) :: month, day_of_month

    ! The first column is for non-leap years, the second column for leap years
    integer, dimension(12, 2), parameter :: DAYS_PER_MONTH = &
         reshape([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31, &
         31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31], [12, 2])

    integer :: ileap, imonth

    if (is_leap(year)) then
       ileap = 2
    else
       ileap = 1
    end if

    day_of_month = day_of_year
    do imonth = 1, 12
       if (day_of_month <= DAYS_PER_MONTH(imonth, ileap)) then
          month = imonth
          exit
       else
          day_of_month = day_of_month - DAYS_PER_MONTH(imonth, ileap)
       end if
    end do
  end subroutine yo_to_ymd

  ! This is a C interface to time_2000_to_ymd_secs for testing purposes
  subroutine time_2000_to_ymd_secs_c(time_2000, year, month, day, seconds_of_year, seconds_of_day) bind(c)
    real(c_double), intent(in) :: time_2000
    integer(c_int), intent(out) :: year, month, day
    real(c_double), intent(out) :: seconds_of_year, seconds_of_day

    call time_2000_to_ymd_secs(time_2000, year, month, day, seconds_of_year, seconds_of_day)
  end subroutine time_2000_to_ymd_secs_c

  ! Compute the number of seconds since 2000-01-01 00:00:00 Z.
  !
  ! This assumes the year is 2000 or later!
  pure function date_parts_to_2000(year, day_of_year, seconds_of_day)
    integer(int32), intent(in) :: year, day_of_year
    real(real64), intent(in) :: seconds_of_day
    real(real64) :: date_parts_to_2000

    integer :: current_year, days_per_year

    date_parts_to_2000 = 0

    ! First handle "fast-forwarding" the number of complete years
    do current_year = 2000, year-1
       if (is_leap(current_year)) then
          days_per_year = 366
       else
          days_per_year = 365
       end if
       date_parts_to_2000 = date_parts_to_2000 + days_per_year * SEC_PER_DAY
    end do

    ! Now add in the number of complete days in this year
    date_parts_to_2000 = date_parts_to_2000 + (day_of_year - 1) * SEC_PER_DAY

    ! Finally, the number of complete seconds in this day
    date_parts_to_2000 = date_parts_to_2000 + seconds_of_day
  end function date_parts_to_2000

  ! Convert time from CCSDS epoch using UTC to J2000 epoch using TT
  !
  ! utc_time: the time in UTC, as represented by the number of seconds
  ! since the CCSDS epoch (1958-01-01 00:00:00 Z). This is not a true
  ! duration since leap seconds are not included (i.e., every day is
  ! exactly 86400 seconds).
  !
  ! leap_seconds: the accumulated total number of leap seconds from
  ! the CCSDS epoch to the date. (In principle this is an integer
  ! value but for MWSDPS compatibility it's a real32.)
  !
  ! The CCSDS epoch is 1958-01-01 00:00:00. The J2000 epoch is
  ! 2000-01-01 12:00:00.
  !
  ! The relationship between UTC and TT is this:
  !
  ! TAI = UTC + DeltaAT
  !
  ! TT = TAI + 32.184 s
  !
  ! Where DeltaAT is the accumulated number of leap seconds and TAI
  ! time is a constant offset from TT time. Note that DeltaAT includes
  ! a 10-second offset, which is the difference between UTC and TAI at
  ! 1972 when UTC was defined.
  pure function ccsds_utc_to_j2000_tt(utc_time, leap_seconds)
    real(real64), intent(in) :: utc_time
    real(real32), intent(in) :: leap_seconds
    real(real64) :: ccsds_utc_to_j2000_tt

    real(real64) :: tai_1958, tai_j2000

    ! Add leap seconds to get TAI time, or the duration in seconds
    ! since the 1958 epoch
    tai_1958 = utc_time + real(leap_seconds, real64)
    ! Shift epoch to J2000
    tai_j2000 = tai_1958 - CCSDS_TO_J2000
    ! Convert from TAI to TT
    ccsds_utc_to_j2000_tt = tai_j2000 + TAI_TO_TT
  end function ccsds_utc_to_j2000_tt

  ! Convert time from CCSDS epoch using UTC to J2000 epoch using UT1
  !
  ! utc_time: the time in UTC, as represented by the number of seconds
  ! since the CCSDS epoch (1958-01-01 00:00:00 Z). This is not a true
  ! duration since leap seconds are not included (i.e., every day is
  ! exactly 86400 seconds).
  !
  ! delta_ut1utc: the current difference between UT1 and UTC. Should
  ! always be within +- 1 second.
  !
  ! The CCSDS epoch is 1958-01-01 00:00:00. The J2000 epoch is
  ! 2000-01-01 12:00:00.
  !
  ! The relationship between UTC and UT1 is this:
  !
  ! UT1 = UTC + Delta(UT1-UTC)
  !
  pure function ccsds_utc_to_j2000_ut1(utc_time, delta_ut1utc)
    real(real64), intent(in) :: utc_time
    real(real32), intent(in) :: delta_ut1utc
    real(real64) :: ccsds_utc_to_j2000_ut1

    real(real64) :: utc_j2000

    ! Shift epoch to J2000
    utc_j2000 = utc_time - CCSDS_TO_J2000
    ! Convert to UT1
    ccsds_utc_to_j2000_ut1 = utc_j2000 + real(delta_ut1utc, real64)
  end function ccsds_utc_to_j2000_ut1

  ! Convert time from CCSDS epoch using UTC to seconds since
  ! 2000-01-01 00:00:00 UTC. This is a TAI time since leap seconds are
  ! not included after 2000.
  !
  ! utc_time: the time in UTC, as represented by the number of seconds
  ! since the CCSDS epoch (1958-01-01 00:00:00 Z). This is not a true
  ! duration since leap seconds are not included (i.e., every day is
  ! exactly 86400 seconds).
  !
  ! leap_seconds: the accumulated total number of leap seconds from
  ! the CCSDS epoch to the date. (In principle this is an integer
  ! value but for MWSDPS compatibility it's a real32.)
  !
  ! The CCSDS epoch is 1958-01-01 00:00:00 TAI/UTC. The new epoch is
  ! 2000-01-01 00:00:00 UTC.
  !
  ! The relationship between UTC and TAI is this:
  !
  ! TAI = UTC + DeltaAT
  !
  ! Where DeltaAT is the accumulated number of leap seconds, including
  ! a 10-second offset, which was the difference between the two when
  ! UTC started in 1972.
  pure elemental function ccsds_utc_to_time_2000(utc_time, leap_seconds)
    real(real64), intent(in) :: utc_time
    real(real32), intent(in) :: leap_seconds
    real(real64) :: ccsds_utc_to_time_2000

    real(real64) :: tai_1958, tai_j2000

    ! Add leap seconds to get TAI time, or the duration in seconds
    ! since the 1958 epoch
    tai_1958 = utc_time + real(leap_seconds, real64)
    ! Shift epoch to J2000
    tai_j2000 = tai_1958 - CCSDS_TO_J2000
    ! Add in the 12 hour offset to get time since midnight
    ccsds_utc_to_time_2000 = tai_j2000 + 0.5 * SEC_PER_DAY
  end function ccsds_utc_to_time_2000

  ! Convert time from seconds since 2000-01-01 00:00:00 Z to
  ! fractional years (e.g., 2013.45423)
  pure function time_2000_to_frac_year(time)
    ! Seconds since 2000-01-01 00:00:00Z
    real(real64), intent(in) :: time
    real(real64) :: time_2000_to_frac_year

    real(real64) :: sec_of_year, sec_of_day, day_per_year
    integer(int32) :: year, day_of_year

    call time_2000_to_yo_secs(time, year, day_of_year, sec_of_year, sec_of_day)

    if (is_leap(year)) then
       day_per_year = real(366, real64)
    else
       day_per_year = real(365, real64)
    endif
    time_2000_to_frac_year=real(year + (day_of_year-1)/day_per_year + sec_of_day / SEC_PER_DAY, real64)
  end function time_2000_to_frac_year

  ! Is this a leap year?
  pure function is_leap(year)
    integer, intent(in) :: year
    logical :: is_leap

    ! https://en.wikipedia.org/w/index.php?title=Leap_year&oldid=852830613#Algorithm
    if (modulo(year, 4) /= 0) then
       is_leap = .false.
    else if (modulo(year, 100) /= 0) then
       is_leap = .true.
    else if (modulo(year, 400) /= 0) then
       is_leap = .false.
    else
       is_leap = .true.
    end if
  end function is_leap

end module time_conversions
