
module time_conversions_extra
  use, intrinsic :: iso_fortran_env, only: int32, real32, real64, IOSTAT_END
  use time_conversions, only: SEC_PER_DAY, TAI_TO_TT, is_leap, date_parts_to_2000
  implicit none
  private
  public :: time_2000_to_mwsdps, time_2000_to_ut1, time_2000_to_tt, ymd_parts_to_2000, acc_leap_secs_2000
  public :: DeltaTDatabase, read_deltat_database, free_deltat_database

  ! On 2000-01-01, there were 22 accumulated leap seconds in UTC since
  ! it was defined in 1972, plus TAI was 10 seconds ahead of UTC when
  ! UTC started. So this value is TAI-UTC at 2000-01-01.
  real(real32), parameter :: LEAPS_AT_2000 = 22 + 10


  ! https://web.archive.org/web/20190911154623/http://maia.usno.navy.mil/ser7/deltat.data
  type DeltaTDatabase
     ! The date (YYYY-MM-DD 00:00:00) converted to seconds since
     ! 2000-01-01 00:00:00 Z (without including leap seconds)
     real(real64), dimension(:), allocatable :: time_2000
     ! The Delta T value from the file, in units of seconds. This is a
     ! continuous function of time that will be linearly interpolated.
     real(real32), dimension(:), allocatable :: delta_t
     ! Total number of accumulated leap seconds. This is needed to
     ! account for the discontinuities in DUT.
     real(real32), dimension(:), allocatable :: leaps
  end type DeltaTDatabase

contains

  ! Read the delta T file
  !
  ! The filename is read from the MWI_DELTAT environment variable
  ! (probably set to "/usr/local/share/mwi/deltat.data"), falling back
  ! to "deltat.data" in the current working directory otherwise.
  subroutine read_deltat_database(db)
    type(DeltaTDatabase), intent(inout) :: db

    integer :: lu, i, stat
    integer :: yyyy, mm, dd
    real(real32) :: delta_t
    integer :: total_lines, env_stat
    character(len=*), parameter :: DT_FORMAT = '(1x, i4.4, 1x, i2.1, 1x, i2.1, 1x, f8.4)'
    character(len=255) :: fname

    ! Get the Delta T filename from the environment or use a fallback
    ! value
    call get_environment_variable("MWI_DELTAT", fname, status=env_stat)
    if (env_stat == 1) then
       ! Variable doesn't exist, so use the fallback value
       fname = "deltat.data"
    else if (env_stat == 0) then
       ! Successfully read the environment variable
       continue
    else
       ! No env var support or the string is too short
       error stop "Didn't read MWI_DELTAT environment variable correctly"
    end if

    open(newunit=lu, file=fname, status='old', action='read', form='formatted')

    ! Do an initial scan of the file to count how many lines it
    ! contains after the J2000 epoch
    total_lines = 0
    do
       read(lu, '(1x, i4.4)', iostat=stat) yyyy
       if (stat == IOSTAT_END) then
          exit
       else if (stat /= 0) then
          error stop 1
       else if (yyyy >= 2000) then
          total_lines = total_lines + 1
       end if
    end do

    ! Rewind the file and allocate space for the data
    rewind(lu)
    allocate(db%time_2000(total_lines), db%delta_t(total_lines), db%leaps(total_lines))

    ! Scan again, this time storing the post-2000 lines
    i = 0
    deltat_lines: do
       read(lu, DT_FORMAT, iostat=stat) yyyy, mm, dd, delta_t
       if (stat == IOSTAT_END) then
          exit deltat_lines
       else if (stat /= 0) then
          error stop
       end if
       ! Skip anything before 2000
       if (yyyy < 2000) cycle

       ! Convert times from YYYY-MM-DD to number of seconds since
       ! 2000-01-01 00:00:00 Z. And by pre-scanning the file we've already
       ! determined that i never exceeds total_lines.
       i = i + 1
       db%time_2000(i) = ymd_parts_to_2000(yyyy, mm, dd, 0, 0, 0)
       db%delta_t(i) = delta_t

       db%leaps(i) = LEAPS_AT_2000 + real(acc_leap_secs_2000(db%time_2000(i)), real32)

       ! Since the inputs are in UTC, the leap seconds need to be accounted for
       db%time_2000(i) = db%time_2000(i) - (db%leaps(i) - LEAPS_AT_2000)
    end do deltat_lines
    close(lu)

    write (*, '(1x, i0, A)') total_lines, " post-2000 lines read from: " // fname
  end subroutine read_deltat_database

  subroutine free_deltat_database(db)
    type(DeltaTDatabase), intent(inout) :: db

    deallocate(db%time_2000, db%delta_t, db%leaps)
  end subroutine free_deltat_database

  ! Look up the DUT, or (UT1 - UTC) value from the database
  !
  ! The input time is the number of seconds since 2000-01-01 00:00:00
  ! Z, not including any leap seconds.
  !
  ! The output is the corresponding value of DUT.
  !
  ! Internally the values of Delta T are linearly interpolated. If the
  ! input time is earlier than or later than the internal time bounds,
  ! the nearest value is used. (In other words, the input time is
  ! clamped to the extrema.)
  !
  ! DUT cannot be directly used for interpolation due to the
  ! discontinuities introduced by leap seconds. So Delta T is
  ! interpolated and then the leap seconds used to convert to DUT.
  pure function query_deltat(db, time_2000)
    type(DeltaTDatabase), intent(in) :: db
    real(real64), intent(in) :: time_2000
    real(real32) :: query_deltat

    real(real32) :: delta_t
    integer :: i, N
    real(real64) :: t, w

    ! To linearly interpolate the Delta T data, the corresponding bin
    ! indices need to be determined such that:
    !
    ! db%time_2000(i) <= time_2000 <= db%time_2000(i+1)
    !
    ! The time values aren't quite regularly spaced but they are
    ! monotonically increasing. So a binary search could be used, but
    ! for simplicity I'll use a linear-time search for now.
    !
    ! The input time is clipped so that the index i is always between
    ! 1 and N-1.
    N = size(db%time_2000)
    t = min(max(time_2000, db%time_2000(1)), db%time_2000(N))
    do i = 1, N-1
       if (t < db%time_2000(i+1)) exit
    end do

    ! If the do-loop hasn't exited early, ensure that i is set properly
    if (i > N-1) i = N-1

    ! For the values db%time_2000(i) and db%time_2000(i+1), determine
    ! the linear interpolation weight w, which varies between 0 and 1.
    ! Then apply the interpolation.
    w = (t - db%time_2000(i)) / (db%time_2000(i+1) - db%time_2000(i))
    delta_t = lerp(db%delta_t(i), db%delta_t(i+1), real(w, real32))

    ! With the interpolated Delta T, now convert to DUT, or UT1 - UTC.
    ! The relation is:
    !
    ! Delta T = 32.184 s + (TAI - UTC) - (UT1 - UTC)
    !
    ! Where (TAI - UTC) is the accumulated number of leap seconds.
    ! Note that the leap seconds are *not* interpolated since they
    ! represent discontinuities. Normally the leap seconds on both
    ! sides (index i and index i+1) are the same, but when a leap
    ! second is introduced, the leap seconds from the left side (i.e.,
    ! bin index i) are used.
    query_deltat = -(delta_t - real(TAI_TO_TT, real32) - db%leaps(i))
  end function query_deltat

  ! Compute the number of seconds since 2000-01-01 00:00:00 Z.
  !
  ! This assumes the year is 2000 or later and that the inputs are
  ! valid (e.g., hour is between 0 and 23).
  pure function ymd_parts_to_2000(year, month, day, hour, minute, second)
    integer(int32), intent(in) :: year, month, day, hour, minute, second
    real(real64) :: ymd_parts_to_2000

    real(real64) :: seconds_of_day
    real(real64), parameter :: SEC_PER_MIN = 60., MIN_PER_HOUR = 60.
    integer(int32) :: day_of_year, leap_seconds

    ! The accumulated days of the year as a function of month (1 to 12)
    integer, dimension(12), parameter :: ACC_DAYS = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334]
    integer, dimension(12), parameter :: ACC_DAYS_LEAP = [0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335]

    ! Convert to ordinal day and seconds within the day
    seconds_of_day = second + SEC_PER_MIN * (minute + MIN_PER_HOUR * hour)
    if (is_leap(year)) then
       day_of_year = ACC_DAYS_LEAP(month) + day
    else
       day_of_year = ACC_DAYS(month) + day
    end if

    ymd_parts_to_2000 = date_parts_to_2000(year, day_of_year, seconds_of_day)

    ! We now have the number of seconds since 2000-01-01 but without
    ! including the effect of leap seconds. We cannot reuse the
    ! "LEAP_TIMES" array or the "acc_leap_secs_2000" function since
    ! the inputs require that the leap seconds are included.
    !
    ! But since we have the UTC input time, we can use a different
    ! approach to include leap seconds.
    !
    ! While the logic below isn't complicated, this technique only
    ! includes leap seconds *after* the fact. In other words, if an
    ! input time falls on a leap second (the "second" input is 60)
    ! then it's not counted as a leap second. But that is accounted
    ! for above since the "seconds_of_day" variable would be 86400
    ! instead of 86399.
    leap_seconds = 0

    ! Leap second of 2005-12-31 23:59:60Z
    if (year > 2005) leap_seconds = leap_seconds + 1
    ! Leap second of 2008-12-31 23:59:60Z
    if (year > 2008) leap_seconds = leap_seconds + 1
    ! Leap second of 2012-06-30 23:59:60Z
    if (year > 2012 .or. (year == 2012 .and. month > 6)) then
       leap_seconds = leap_seconds + 1
    end if
    ! Leap second of 2015-06-30 23:59:60Z
    if (year > 2015 .or. (year == 2015 .and. month > 6)) then
       leap_seconds = leap_seconds + 1
    end if
    ! Leap second of 2016-12-31 23:59:60Z
    if (year > 2016) leap_seconds = leap_seconds + 1

    ! Add in the missing leap seconds
    ymd_parts_to_2000 = ymd_parts_to_2000 + real(leap_seconds, real64)
  end function ymd_parts_to_2000

  ! The accumulated number of leap seconds since 2000-01-01
  !
  ! The input time is in seconds since 2000-01-01 00:00:00 Z, not
  ! including any leap seconds.
  !
  ! The output is the number of leap seconds that have happened
  ! between 2000 and the input time.
  pure function acc_leap_secs_2000(time_2000)
    real(real64), intent(in) :: time_2000
    integer :: acc_leap_secs_2000

    ! These are the times that leap seconds were inserted since 2000.
    ! This is indexed as the input time (continuous seconds since
    ! 2000-01-01 00:00:00 UTC). The values are computed using Python:
    !
    ! (datetime(leap_second) - datetime(2000, 1, 1)).total_seconds()
    !
    ! But note that the accumulated number of seconds needs to be
    ! included with each term. And the datetime() class can't accept
    ! 60 as a second so an extra second is added afterward.
    !
    ! (Technically leap seconds may be removed but that's never
    ! happened and so this approach doesn't account for that.)
    real(real64), dimension(5), parameter :: LEAP_TIMES = [ &
         189388800., &  ! 2005-12-31 23:59:60 + 0 seconds
         284083201., &  ! 2008-12-31 23:59:60 + 1 second
         394416002., &  ! 2012-06-30 23:59:60 + 2 seconds
         489024003., &  ! 2015-06-30 23:59:60 + 3 seconds
         536544004.  &  ! 2016-12-31 23:59:60 + 4 seconds
         ]

    acc_leap_secs_2000 = count(time_2000 >= LEAP_TIMES)
  end function acc_leap_secs_2000

  ! Convert time from the 2000 epoch to MWSDPS values.
  !
  ! The input time is the number of SI seconds since 2000-01-01
  ! 00:00:00 UTC (no leap seconds included).
  !
  ! The output values are shifted to the CCSDS epoch of 1958-01-01.
  !
  ! time_ccsds is the number of seconds since the 1958 epoch to obtain
  ! the time in UTC of the measurement. It includes leap seconds.
  !
  ! leap_seconds is the accumulated number of leap seconds from the
  ! epoch to the input time. (It technically only takes on integer
  ! values, but for compatibility with MWSDPS it's defined as a real32.)
  !
  ! dut1utc is the current Delta(UT1-UTC) term. It should be within
  ! +-1 second since UTC is shifted by leap seconds to stay within 1
  ! second of UT1.
  !
  ! An assumption is that we're not dealing with time inputs before
  ! 2000. In other words, the time_2000 input is nonnegative.
  pure subroutine time_2000_to_mwsdps(dut_db, time_2000, time_ccsds, leap_seconds, dut1utc)
    type(DeltaTDatabase), intent(in) :: dut_db
    real(real64), intent(in) :: time_2000
    real(real64), intent(out) :: time_ccsds
    real(real32), intent(out) :: leap_seconds, dut1utc

    ! This is the difference in seconds between 1958-01-01 00:00:00
    ! and 2000-01-01 00:00:00 *without* including leap seconds. It is
    ! 42 years exactly. Of the 42 years, 10 of them are leap years.
    real(real64), parameter :: T2000_TO_CCSDS = 15340.0 * SEC_PER_DAY

    ! Accumulate leap seconds
    leap_seconds = LEAPS_AT_2000 + acc_leap_secs_2000(time_2000)

    ! DUT, or Delta(UT1 - UTC)
    dut1utc = query_deltat(dut_db, time_2000)

    ! Include leap seconds and shift the epoch
    time_ccsds = time_2000 - real(leap_seconds, real64) + T2000_TO_CCSDS
  end subroutine time_2000_to_mwsdps

  ! Convert time from the 2000 epoch to UT1 at the J2000 epoch
  !
  ! The input time is the number of SI seconds since 2000-01-01
  ! 00:00:00 UTC (no leap seconds included).
  !
  ! The output is the number of UT1 seconds since the J2000 epoch
  ! (2000-01-01 12:00:00 UT1).
  pure function time_2000_to_ut1(time_2000)
    real(real64), intent(in) :: time_2000
    real(real64) :: time_2000_to_ut1

    real(real64) :: time_j2000

    ! Subtract the 12-hour phase to get from midnight to noon
    time_j2000 = time_2000 - 0.5 * SEC_PER_DAY

    ! There are two more things to account for: the leap seconds
    ! present in UTC and the difference between UT1 and UTC. The input
    ! time does not include leap seconds.
    !
    ! As of 2017-01-01, there have been 5 leap seconds added since
    ! 2000. UT1 and UTC are kept within a second of each other but
    ! vary as a function of time. So, properly, this function should
    ! include both those terms.
    !
    ! However for the TDR simulation, the UT1 time is only used to
    ! compute the GMST. It's only a small error and for the time being
    ! is left unaccounted for.
    time_2000_to_ut1 = time_j2000
  end function time_2000_to_ut1

  ! Convert time from the 2000 epoch to TT at the J2000 epoch
  !
  ! The input time is the number of SI seconds since 2000-01-01
  ! 00:00:00 UTC (no leap seconds included).
  !
  ! The output is the number of days since the J2000 epoch (2000-01-01
  ! 12:00:00 TT).
  pure function time_2000_to_tt(time_2000)
    real(real64), intent(in) :: time_2000
    real(real64) :: time_2000_to_tt

    real(real64) :: time_j2000
    ! On 2000-01-01, TAI was ahead of UTC by this much (this includes
    ! the 22 accumulated leap seconds up to 2000)
    real(real64), parameter :: DELTA_AT_2000 = 32.

    ! The relationship between TAI, TT, and UTC is this:
    !
    ! TAI = UTC + DeltaAT
    !
    ! TT = TAI + 32.184 s
    !
    ! At the epoch (both 2000 and J2000), DeltaAT was 32 seconds. So
    ! the time in TT from J2000 is computed by adding the two terms,
    ! plus subtracting the 12-hour phase to shift from midnight to
    ! noon.
    time_j2000 = time_2000 + DELTA_AT_2000 + TAI_TO_TT - 0.5 * SEC_PER_DAY

    ! Convert to days
    time_2000_to_tt = time_j2000 / SEC_PER_DAY
  end function time_2000_to_tt

  ! Linear interpolation
  !
  ! t: ranges between 0 and 1
  ! v0: value when t == 0
  ! v1: value when t == 1
  pure function lerp(v0, v1, t)
    real(real32), intent(in) :: v0, v1, t
    real(real32) :: lerp
    lerp = (1 - t) * v0 + t * v1
  end function lerp

end module time_conversions_extra
