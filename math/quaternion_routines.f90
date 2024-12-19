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
! Small support module for quaternion operations and conversions
!
! The quaternion convention used is "xyzw" where the scalar part (w)
! is after the vector part (xyz).
module quaternion_routines
  use, intrinsic :: iso_fortran_env, only: real64
  implicit none
  private
  public :: slerp, q_normalize, normalize_time, att_q_to_rpy
  public :: rotation_matrix_to_quaternion, quaternion_to_rotation_matrix

contains

  ! Convert the quaternion to a rotation matrix.
  !
  ! The quaternion convention is xyzw, where the scalar part is after
  ! the vector part.
  pure function quaternion_to_rotation_matrix(q) result(mat)
    real(real64), dimension(4), intent(in) :: q
    real(real64), dimension(3,3) :: mat

    associate(a => q(4), b => q(1), c => q(2), d => q(3))
      mat(1,1) = a**2 + b**2 - c**2 - d**2
      mat(2,1) = 2*b*c + 2*a*d
      mat(3,1) = 2*b*d - 2*a*c
      mat(1,2) = 2*b*c - 2*a*d
      mat(2,2) = a**2 - b**2 + c**2 - d**2
      mat(3,2) = 2*c*d + 2*a*b
      mat(1,3) = 2*b*d + 2*a*c
      mat(2,3) = 2*c*d - 2*a*b
      mat(3,3) = a**2 - b**2 - c**2 + d**2
    end associate
  end function quaternion_to_rotation_matrix

  ! Convert a rotation matrix into a quaternion.
  !
  ! This uses the convention where the scalar part is after the vector
  ! part (xyzw).
  !
  ! References:
  !
  ! https://en.wikipedia.org/wiki/Rotation_matrix#Quaternion
  ! http://www.tu-berlin.de/fileadmin/fg169/miscellaneous/Quaternions.pdf
  ! https://d3cw3dd2w32x2b.cloudfront.net/wp-content/uploads/2015/01/matrix-to-quat.pdf
  pure function rotation_matrix_to_quaternion(rot_mat) result(q)
    real(real64), dimension(3, 3), intent(in) :: rot_mat
    real(real64), dimension(4) :: q

    real(real64) :: t

    ! The general approach is that one quaternion element (x, y, z, or
    ! w) is directly computed using the diagonal of the rotation
    ! matrix. The remaining quaternion elements are built up using
    ! off-diagonal entries and all elements are normalized using that
    ! initial element.
    !
    ! If the initial element used is too small, then the normalization
    ! step leads to instabilities or even a division-by-zero. So we
    ! have to pick the initial quaternion element that has the largest
    ! value.
    !
    ! Thus, there are four cases to select from. This uses the
    ! approach in the "matrix-to-quat.pdf" URL above which avoids some
    ! repeated computation that a more straight-forward approach would
    ! take. However, note that its convention uses the transpose of
    ! the rotation matrix. For clarity I'll reuse the variable
    ! notation and then undo the transpose later by taking the
    ! conjugate of the quaternion.
    associate( &
         m00 => rot_mat(1, 1), m11 => rot_mat(2, 2), m22 => rot_mat(3, 3), &
         m01 => rot_mat(1, 2), m02 => rot_mat(1, 3), &
         m10 => rot_mat(2, 1), m12 => rot_mat(2, 3), &
         m20 => rot_mat(3, 1), m21 => rot_mat(3, 2))

      if (m22 < 0) then
         if (m00 > m11) then
            ! The q_x component is the largest
            t = 1 + m00 - m11 - m22
            q(:) = [t, m01 + m10, m20 + m02, m12 - m21]
         else
            ! The q_y component is the largest
            t = 1 - m00 + m11 - m22
            q(:) = [m01 + m10, t, m12 + m21, m20 - m02]
         end if
      else
         if (m00 < -m11) then
            ! The q_z component is the largest
            t = 1 - m00 - m11 + m22
            q(:) = [m20 + m02, m12 + m21, t, m01 - m10]
         else
            ! The q_w component is the largest
            t = 1 + m00 + m11 + m22
            q(:) = [m12 - m21, m20 - m02, m01 - m10, t]
         end if
      end if
    end associate
    ! This is the final normalization step but I add in the
    ! conjugation operation to account for the transposed matrix
    ! convention
    q(:) = [-q(1), -q(2), -q(3), q(4)] * 0.5 / sqrt(t)
  end function rotation_matrix_to_quaternion


  ! Spherical linear interpolation
  !
  ! This interpolates between the two quaternions q0 and q1 using the
  ! parameter t, which varies from 0 to 1.
  !
  ! Adapted from: https://en.wikipedia.org/w/index.php?title=Slerp&oldid=941403106#Source_code

  pure function slerp(q0, q1, t)
    real(real64), dimension(4), intent(in) :: q0, q1
    real(real64), intent(in) :: t
    real(real64), dimension(4) :: slerp

    ! A local copy of q1 is taken here since it is possibly negated. A
    ! copy is required since this function doesn't mutate its inputs.
    ! (This differs from the C++ implementation this is based on,
    ! which takes the arguments by value.)
    real(real64), dimension(4) :: q1_local

    real(real64) :: dot, theta, theta_0, s0, s1
    real(real64) :: sin_theta, sin_theta_0
    real(real64), parameter :: DOT_THRESHOLD = 0.9995

    dot = dot_product(q0, q1)

    ! If the dot product is negative, then the interpolation won't
    ! take the shortest path. So reverse one of the quaternions.
    if (dot < 0) then
       dot = -dot
       q1_local(:) = -q1(:)
    else
       q1_local(:) = q1(:)
    end if

    if (dot > DOT_THRESHOLD) then
       ! If the angle between the quaternions is very small (i.e., the
       ! dot product close to 1) then do a linear interpolation and
       ! then normalize
       slerp = q_normalize(q0 + t * (q1_local - q0))
       return
    end if

    ! At this point the argument to acos() is between 0 and DOT_THRESHOLD
    theta_0 = acos(dot)

    ! This is the interpolated angle between the two quaternions
    theta = theta_0 * t

    sin_theta = sin(theta)
    sin_theta_0 = sin(theta_0)

    s0 = cos(theta) - dot * sin_theta / sin_theta_0
    s1 = sin_theta / sin_theta_0

    slerp(:) = s0 * q0(:) + s1 * q1_local(:)
  end function slerp


  ! Normalize a quaternion
  pure function q_normalize(q0)
    real(real64), dimension(4), intent(in) :: q0
    real(real64), dimension(4) :: q_normalize

    q_normalize(:) = q0(:) / norm2(q0)
  end function q_normalize

  ! For an input time which is between time_start and time_stop,
  ! normalize it so it varies between 0 (if time == time_start) and 1
  ! (if time == time_stop).
  !
  ! If the input time is outside the boundaries, then the output may
  ! be outside 0 and 1. In other words, no checking is performed.
  pure function normalize_time(time, time_start, time_stop)
    real(real64), intent(in) :: time, time_start, time_stop
    real(real64) :: normalize_time

    real(real64) :: dt

    ! To avoid a division by zero, avoid time intervals too small (or
    ! zero)
    dt = time_stop - time_start
    if (abs(dt) < 1e-6) then
       normalize_time = 0
    else
       normalize_time = (time - time_start) / dt
    end if
  end function normalize_time

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

end module quaternion_routines
