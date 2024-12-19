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
! These are additional helper routines that are used by the
! geolocation algorithm but not the resampling algorithms.
module geolocation_extra_routines
  use, intrinsic :: iso_fortran_env, only: real32, real64, int32, int8
  use dem, only: DigitalElevationMap
  use geolocation_core_routines, only: EarthFigure, GeolocationScanlineData, &
       loccel_range, get_sc_axes, point_along_boresight, find_earth_intersection, &
       point_to_nadir
  use math_routines, only: fixang8, dot_product_unit8, cross_norm
  use mwsdps_geolocation_types, only: GeolocationSensorData_f, GeolocationParameterData_f
  use pre_nut_routines, only: get_gm_angle
  use quaternion_routines, only: &
       rotation_matrix_to_quaternion, quaternion_to_rotation_matrix, slerp, normalize_time
  use trig_degrees, only: dsind, dcosd, dasind, dacosd, datand, datan2d
  implicit none
  private
  public :: find_subsat_point, intersects_earth, find_polrot_angle, &
       convert_sc_to_ecef, find_cold_view_outputs
#ifndef MWSDPS
  public :: att_q_to_rpy
#endif

contains

  ! Find the location of the cold sky spillover.
  !
  ! This basically repeats the geolocation approach but for a single
  ! "FOV" with a per-band constant scan angle. One important
  ! difference is that for the normal geolocation implementation, for
  ! each FOV, a time difference is computed and then the corresponding
  ! scan angle. Here it's the other way around: the scan angle is
  ! known but not the time of that scan angle.
  pure subroutine find_cold_view_outputs(scanline, parameters, ellipsoid, input_data, iscan, iband, lat, lon)
    type(GeolocationScanlineData), intent(in) :: scanline
    type(GeolocationParameterData_f), intent(in) :: parameters
    type(EarthFigure), intent(in) :: ellipsoid
    type(GeolocationSensorData_f), intent(in) :: input_data
    integer, intent(in) :: iscan, iband
    real(real64), intent(out) :: lat, lon

    real(real64) :: time_at_fov, scan_angle_at_fov, delta_scan_angle
    real(real64) :: alt, gm_angle, normalized_scan_time, time_since_scan_start
    real(real64), dimension(3) :: scpos_at_fov, boresight, surf_normal
    real(real64), dimension(4) :: att_at_fov
    real(real64), dimension(3) :: sc_x, sc_y, sc_z
    logical :: found

    real(real64), parameter :: RPM_TO_DEG_PER_SEC = 360. / 60.

    ! This is the scan angle corresponding to the cold sky mirror
    scan_angle_at_fov = parameters%ColdSkyScanAngle(iband)

    ! Find the time of the cold sky mirror scan angle based off
    ! the scan rate and the known scan angle and time of the
    ! first sample of the scan
    delta_scan_angle = (scan_angle_at_fov - input_data%FirstSampleScanAngle(iscan, iband))
    time_since_scan_start = delta_scan_angle / (input_data%ScanRate(iscan) * RPM_TO_DEG_PER_SEC) &
         + input_data%FirstSampleDelay(iscan, iband)

    time_at_fov = scanline%time_j2000_ut1 + time_since_scan_start
    call get_gm_angle(time_at_fov, gm_angle)

    ! Propagate the satellite position assuming constant
    ! velocity through the scan
    scpos_at_fov(:) = scanline%scpos(:) + time_since_scan_start * scanline%scvel(:)

    ! Find the SACS-S axes at the mean epoch of date (nominally
    ! oriented pointing toward geodetic nadir), interpolating the
    ! attitude along the scan
    normalized_scan_time = normalize_time(time_at_fov, scanline%time_this_scan, scanline%time_next_scan)
    att_at_fov = slerp(scanline%att_this_scan, scanline%att_next_scan, normalized_scan_time)
    call get_sc_axes(att_at_fov, scanline%precession, sc_x, sc_y, sc_z)

    call point_along_boresight(sc_x, sc_y, sc_z, &
         scan_angle_at_fov, real(parameters%OffNadirAngle(iband), real64), boresight)
    call find_earth_intersection(ellipsoid, scpos_at_fov, boresight, gm_angle, &
         lat, lon, alt, surf_normal, found)
  end subroutine find_cold_view_outputs

#ifndef MWSDPS
  ! Convert the attitude quaternion to satellite roll/pitch/yaw.
  !
  ! Attitude is a quaternion, roll/pitch/yaw are in degrees.
  !
  ! The quaternion convention is xyzw, where the scalar part is after
  ! the vector part.
  pure subroutine att_q_to_rpy(attitude, rpy)
    real(real64), dimension(4), intent(in) :: attitude
    real(real64), dimension(3), intent(out) :: rpy

    real(real64) :: rot_mat_11, rot_mat_21, rot_mat_31, rot_mat_32, rot_mat_33

    ! Rather than building the entire 9-element rotation matrix, only
    ! 5 values from it are needed
    associate(a => attitude(4), b => attitude(1), c => attitude(2), d => attitude(3))
      rot_mat_11 = a**2 + b**2 - c**2 - d**2
      rot_mat_21 = 2*b*c + 2*a*d
      rot_mat_31 = 2*b*d - 2*a*c
      rot_mat_32 = 2*c*d + 2*a*b
      rot_mat_33 = a**2 - b**2 - c**2 + d**2
    end associate

    ! Extract roll, pitch, yaw. Following is the 'rpy' or '123'
    ! convention.
    rpy(1) = -datan2d(rot_mat_32, rot_mat_33)
    rpy(2) = dasind(rot_mat_31)
    rpy(3) = -datan2d(rot_mat_21, rot_mat_11)
  end subroutine att_q_to_rpy
#endif

  ! Checks if the boresight vector intersects the Earth or not
  pure function intersects_earth(ellipsoid, scpos, boresight)
    type(EarthFigure), intent(in) :: ellipsoid
    real(real64), dimension(3), intent(in) :: scpos, boresight
    logical :: intersects_earth

    real(real64), dimension(3) :: ru, cell
    real(real64) :: r, rearth, range

    r=norm2(scpos)
    ru=scpos/r

    call loccel_range(ellipsoid, r,ru,boresight, cell,rearth,range, intersects_earth)
  end function intersects_earth

  ! Find the subsatellite point
  !
  ! ellipsoid: Earth reference ellipsoid parameters
  ! scpos: satellite position in ECI coordinates (meters)
  ! rotangle: Greenwich mean sidereal angle in degrees
  !
  ! lat/lon: the geodetic latitude/longitude
  ! range: the altitude of the satellite in meters
  ! found is true if the intersection exists and false if there isn't an intersection
  pure subroutine find_subsat_point(ellipsoid, scpos, rotangle, lat, lon, range, found)
    type(EarthFigure), intent(in) :: ellipsoid
    real(real64), dimension(3), intent(in) :: scpos
    real(real64), intent(in) :: rotangle
    real(real64), intent(out) :: lat, lon, range
    logical, intent(out) :: found

    real(real64), dimension(3) :: r_u, cell, b
    real(real64) :: r, r_earth, coslat, lon_eci
    real(real64) :: RE, RP, FFAC

    RE = ellipsoid%EarthEquatorialRadius
    RP = ellipsoid%EarthPolarRadius
    FFAC = (RP / RE)**2

    r=norm2(scpos)
    r_u=scpos/r

    call point_to_nadir(ellipsoid, scpos, b)

    call loccel_range(ellipsoid, r,r_u,b, cell,r_earth,range,found)

    if (found) then
       ! cosine of geocentric lat
       coslat=hypot(cell(1), cell(2))
       ! geodetic cell latitude
       lat=datand(cell(3)/(FFAC*coslat))
       lon_eci=datan2d(cell(2), cell(1))
       call fixang8(lon_eci)

       ! convert from inertial longitude to earth-fixed
       lon=lon_eci - rotangle
       call fixang8(lon)
    end if
  end subroutine find_subsat_point

  ! Find the polarization rotation angle in degrees.
  !
  ! boresight: (normalized) boresight vector in ECI frame
  ! nadir: (normalized) vector pointing from spacecraft to nadir
  ! surface_normal: (normalized) Earth surface normal vector in ECI frame
  !
  ! near_nadir: if the incidence angle is very small (< 0.01 degrees
  ! perhaps), the h-pol orientation at the Earth is difficult to
  ! define in this case
  !
  ! Note that the nadir vector is the -Z axis of the SACS-S frame.
  !
  ! Refer to the "Polarization Rotation Angle" section of the MWI
  ! geolocation ATBD.
  pure subroutine find_polrot_angle(boresight, nadir, surface_normal, near_nadir, polarization_rotation)
    real(real64), dimension(3), intent(in) :: boresight, nadir, surface_normal
    logical, intent(in) :: near_nadir
    real(real64), intent(out) :: polarization_rotation

    ! H-pol and V-pol directions at the antenna or at the Earth
    real(real64), dimension(3) :: ha, he, ve

    ! Some important things to note for this routine:
    !
    ! The propagation direction vector (k) is the reverse of the
    ! boresight direction vector (b), or, b = -k. While several of the
    ! polarization directions are defined in terms of k, instead b is
    ! used since that is an input to the routine. k doesn't need to be
    ! explicitly formed here.
    !
    ! The cross product is anticommutative, so -a x b = b x a. In
    ! other words, when one vector is negated the order of the cross
    ! product inputs must be reversed. If both vectors are negated the
    ! order "reverses twice", undoing the effect.

    ! At the antenna, the h-pol direction is defined as the
    ! propagation direction cross the zenith direction. In this case,
    ! both vectors are negated.
    call cross_norm(boresight, nadir, ha)

    if (.not. near_nadir) then
       ! At the Earth surface, the h-pol direction is defined as the
       ! propagation direction cross the surface normal. One vector is
       ! negated so the cross product order must be reversed.
       call cross_norm(surface_normal, boresight, he)
    else
       ! At nadir, the h-pol direction is not defined, so approximate
       ! this case
       he = ha
    end if

    ! V-pol direction relative to the Earth: h-pol cross the
    ! propagation direction. One vector is negated so the cross
    ! product order must be reversed.
    call cross_norm(boresight, he, ve)

    ! The relationship between the four polarization vectors (only
    ! three are needed below) is a rotation using the polarization
    ! rotation angle (pra):
    !
    ! va = ve cos pra - he sin pra
    ! ha = ve sin pra + he cos pra
    !
    ! Due to the orthonormal bases, the dot products for the matching
    ! pairs (antenna and Earth) are:
    !
    ! ve . ve = he . he = va . va = ha . ha = 1
    ! ve . he = va . ha = 0
    !
    ! Taking the dot product between three of the polarization vectors
    ! isolates the angle:
    !
    ! ve . ha = sin pra
    ! he . ha = cos pra
    !
    ! Thus, tan pra = (ve . ha) / (he . ha)
    polarization_rotation = datan2d(dot_product(ve, ha), dot_product(he, ha))
  end subroutine find_polrot_angle


  ! Convert the spacecraft position and orientation from ECI
  ! coordinates at the J2000 epoch to the mean epoch of date, ECEF
  ! coordinates
  pure subroutine convert_sc_to_ecef(scpos_j2000, scvel_j2000, scatt_j2000, &
       precession, gmst, &
       scpos_ecef, scvel_ecef, scatt_ecef)
    real(real32), dimension(3), intent(in) :: scpos_j2000, scvel_j2000
    real(real32), dimension(4), intent(in) :: scatt_j2000
    real(real64), dimension(3, 3), intent(in) :: precession
    real(real64), intent(in) :: gmst
    real(real32), dimension(3), intent(out) :: scpos_ecef, scvel_ecef
    real(real32), dimension(4), intent(out) :: scatt_ecef

    real(real64), dimension(3) :: scpos_eod, scvel_eod
    real(real64), dimension(3, 3) :: gmst_rotation
    real(real64), dimension(3, 3) :: att_mat_j2000, att_mat_eod

    ! Apply precession to obtain mean epoch of date, still in ECI
    scpos_eod = matmul(precession, real(scpos_j2000, real64))
    scvel_eod = matmul(precession, real(scvel_j2000, real64))

    ! Then apply GMST to obtain ECEF (namely, rotation about Z-axis with an
    ! angle of -GMST)
    gmst_rotation(1, :) = [real(real64) :: dcosd(-gmst), -dsind(-gmst), 0]
    gmst_rotation(2, :) = [real(real64) :: dsind(-gmst), dcosd(-gmst), 0]
    gmst_rotation(3, :) = [0, 0, 1]
    scpos_ecef(:) = real(matmul(gmst_rotation, scpos_eod), real32)
    scvel_ecef(:) = real(matmul(gmst_rotation, scvel_eod), real32)

    ! Apply precession to attitude quaternion after converting it to a
    ! rotation matrix
    att_mat_j2000 = quaternion_to_rotation_matrix(real(scatt_j2000, real64))
    att_mat_eod = matmul(precession, att_mat_j2000)

    ! Convert attitude matrix back to a quaternion
    scatt_ecef(:) = real(rotation_matrix_to_quaternion(att_mat_eod), real32)
  end subroutine convert_sc_to_ecef

end module geolocation_extra_routines
