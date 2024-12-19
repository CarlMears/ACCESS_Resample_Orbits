$freeform
module interpolate
	use, INTRINSIC :: IEEE_ARITHMETIC
    !use SeaIceCtypes, only: c_int, 4, 8, 4, c_int8_t
    use NAN_support, only: nan_f32
    
	character(len=20),dimension(0:3),parameter		:: INTERP_ERRORS = (/'No Error','Error Given Bounds','Error in X or Y','Error in Logic'/)

contains

	real(4) function interpolate_2D(x,y,mat,x_min,x_max,y_min,y_max,error,x_wrap,y_wrap)

		

		implicit none

		real(4),intent(in)					:: x
		real(4),intent(in)					:: y
		real(4),intent(in),dimension(0:,0:)	:: mat
		real(4),intent(in)					:: x_min
		real(4),intent(in)					:: x_max
		real(4),intent(in)					:: y_min
		real(4),intent(in)					:: y_max
		integer(4),intent(out),optional		:: error
		logical,intent(in),optional			:: x_wrap
		logical,intent(in),optional			:: y_wrap 

		logical								:: xwrap
		logical								:: ywrap
		logical								:: set_error

		integer(4)							:: num_x,num_y
		real(4)								:: dx,dy
		real(4)								:: x_scl,y_scl
		integer(4)							:: x1,x2,y1,y2
		real(4)								:: wt_x_1,wt_x_2,wt_y_1,wt_y_2
		

		! check for presence of optional parameters, set to default if
		! not present

		if (present(x_wrap)) then
			xwrap = x_wrap
		else
			xwrap = .false.
		endif

		if (present(y_wrap)) then
			ywrap = y_wrap
		else
			ywrap = .false.
		endif

		if (present(error)) then
			set_error = .true.
		else
			set_error = .false.
		endif

		! check to make sure x_min lt x_max, etc

		if (x_max .le. x_min) then
			if (set_error) then
				error = 1
			endif
			interpolate_2D = nan_f32
			return
		endif

		if (y_max .le. y_min) then
			if (set_error) then
				error = 1
			endif
			interpolate_2D = nan_f32
			return
		endif

		num_x = ubound(mat,1) + 1
		num_y = ubound(mat,2) + 1

		if ((num_x .le. 1) .or. (num_y .le. 1)) then

			if (set_error) then
				error = 1
			endif

			interpolate_2D = nan_f32
			return
		endif

		dx = (x_max - x_min)/float(num_x - 1)
		dy = (y_max - y_min)/float(num_y - 1)

		! check to make sure x and y are in bounds

		if (xwrap) then
			if ((x .lt. x_min) .or. (x .gt. x_max + dx)) then ! x out of bounds for wrapped case

				if (set_error) then
					error = 2
				endif

				interpolate_2D = nan_f32
				return
			endif
		else
			if ((x .lt. x_min) .or. (x .gt. x_max)) then ! x out of bounds for unwrapped case

				if (set_error) then
					error = 2
				endif

				interpolate_2D = nan_f32
				return
			endif
		endif

 		if (ywrap) then
			if ((y .lt. y_min) .or. (y .gt. y_max + dy)) then ! x out of bounds for wrapped case

				if (set_error) then
					error = 2
				endif

				interpolate_2D = nan_f32
				return
			endif
		else
			if ((y .lt. y_min) .or. (y .gt. y_max)) then ! x out of bounds for unwrapped case

				if (set_error) then
					error = 2
				endif

				interpolate_2D = nan_f32
				return
			endif
		endif

		! Everything seems OK -- do the interpolation

		x_scl = (x - x_min)/dx
		y_scl = (y - y_min)/dy

		x1 = floor(x_scl)
		x2 = x1 + 1
		y1 = floor(y_scl)
		y2 = y1 + 1

		wt_x_2 = (x_scl  - x1)/dx
		wt_x_1 = 1.0 - wt_x_2
		wt_y_2 = (y_scl  - y1)/dx
		wt_y_1 = 1.0 - wt_y_2

		if ((xwrap) .and. (x2 .eq. num_x)) x2 = 0
		if ((ywrap) .and. (y2 .eq. num_y)) y2 = 0

		if ((x1 .ge. 0) .and. (x1 .le. (num_x -1)) .and. (x2 .ge. 0) .and. (x2 .le. (num_x-1)) .and. &
			(y1 .ge. 0) .and. (y1 .le. (num_y -1)) .and. (y2 .ge. 0) .and. (y2 .le. (num_y-1))) then ! everything OK

			interpolate_2D = wt_x_1*wt_y_1*mat(x1,y1) + &
							 wt_x_1*wt_y_2*mat(x1,y2) + &
							 wt_x_2*wt_y_1*mat(x2,y1) + &
							 wt_x_2*wt_y_2*mat(x2,y2)
			if (set_error) then
				error = 0
			endif
		else	
			interpolate_2D = nan_f32
			if (set_error) then
				error = 3
			endif
		endif
	end function interpolate_2D

	real(4) function lin_interpolate_r4_r8(x,x1,x2,y1,y2,error)

		real(4),intent(in)			:: x
		real(4),intent(in)			:: x1
		real(4),intent(in)			:: x2
		real(8),intent(in)			:: y1
		real(8),intent(in)			:: y2
		integer(4),intent(inout)	:: error

		real(8)						:: wt1,wt2
		
		! do some bounds checking
		
		if (x .lt. x1) then 
			error =1
		endif

		if (x2 .lt. x) then 
			error =1
		endif

		if (x2 .le. x1) then 
			error = -1
			lin_interpolate_r4_r8 = nan_f32
			return
		endif

		wt1 = (dble(x2) - dble(x))/(dble(x2) - dble(x1))
		wt2 = 1.0 - wt1

		lin_interpolate_r4_r8 = wt1*y1 + wt2*y2

		return

	end function lin_interpolate_r4_r8


	real(4) function lin_interpolate_r8_r4(x,x1,x2,y1,y2,error)

		real(8),intent(in)			:: x
		real(8),intent(in)			:: x1
		real(8),intent(in)			:: x2
		real(4),intent(in)			:: y1
		real(4),intent(in)			:: y2
		integer(4),intent(inout)	:: error

		real(4)						:: wt1,wt2
		
		! do some bounds checking
		
		if (x .lt. x1) then 
			error =1
		endif

		if (x2 .lt. x) then 
			error =1
		endif

		if (x2 .le. x1) then 
			error = -1
			lin_interpolate_r8_r4 = nan_f32
			return
		endif

		wt1 = (x2 - x)/(x2 - x1)
		wt2 = 1.0 - wt1

		lin_interpolate_r8_r4 = wt1*y1 + wt2*y2

		return

	end function lin_interpolate_r8_r4


	real(4) function lin_interp_1D(x0,x,y,error,extrapolate)

		! this function performs linear interpolation on
		! 1-D arrays of x and y values

		real(4),intent(in)					:: x0
		real(4),intent(in),dimension(:)		:: x
		real(4),intent(in),dimension(:)		:: y
		logical(4),intent(in),optional      :: extrapolate
		integer(4),intent(inout)			:: error

		integer(4)							:: nx,ny,j
		
		logical(4)                          :: extrap
		
		extrap = .false.
		if (present(extrapolate))extrap = extrapolate

		nx = size(x)
		ny = size(y)

		if ((nx <= 1) .or. (nx /= ny)) then
			error = -1
			lin_interp_1D = nan_f32
			return
		endif

		call locate(x,x0,j)

		if ((j == 0) .or. (j >= nx)) then
		    if (extrap .eq. .false.) then
			    error = 1
			    lin_interp_1D = nan_f32
			    return
			 else
			    if (j == 0) then
			        lin_interp_1D = lin_interpolate_no_bounds_check(x0,x(1),x(2),y(1),y(2),error)
			        return
			    endif
			    if (j == nx)then
			        lin_interp_1D = lin_interpolate_no_bounds_check(x0,x(nx-1),x(nx),y(nx-1),y(nx),error)
			        return
			    endif
			 endif
		endif

        if (x(j+1) .gt. x(j)) then
			lin_interp_1D = lin_interpolate(x0,x(j),x(j+1),y(j),y(j+1),error)
		else
			lin_interp_1D = lin_interpolate(x0,x(j+1),x(j),y(j+1),y(j),error)
		endif
		return

	end function lin_interp_1D

	subroutine locate(xx,x,j)

	  ! adapted from numerical recipes

	  real(4),intent(in),dimension(:)	:: xx
	  real(4),intent(in)				:: x
	  integer(4),intent(out)			:: j

	  integer(4)						:: n
	  integer(4)						:: jl
	  integer(4)						:: jm
	  integer(4)						:: ju

	  n = size(xx)
	  jl = 0
	  ju = n+1

      do while(ju-jl > 1)
        jm=(ju+jl)/2
        if((xx(n) >= xx(1)) .eqv. (x >= xx(jm))) then
          jl=jm
        else
          ju=jm
        endif
      enddo

      if(x == xx(1))then
        j=1
      else if(x == xx(n))then
        j=n-1
      else
        j=jl
      endif
      return
    end subroutine locate

	real(4) function lin_interpolate(x,x1,x2,y1,y2,error)

		real(4),intent(in)			:: x
		real(4),intent(in)			:: x1
		real(4),intent(in)			:: x2
		real(4),intent(in)			:: y1
		real(4),intent(in)			:: y2
		integer(4),intent(inout)	:: error

		real(4)						:: wt1,wt2

		error = 0
		if (.not. ieee_is_finite(x))  error = -1
		if (.not. ieee_is_finite(x1)) error = -1
		if (.not. ieee_is_finite(x2)) error = -1
		if (.not. ieee_is_finite(y1)) error = -1
		if (.not. ieee_is_finite(y2)) error = -1
		lin_interpolate = nan_f32
		if (error .ne. 0)  then
			lin_interpolate = nan_f32
			RETURN
		endif


		
		! do some bounds checking
		
		if (x .lt. x1) then 
			error =1
		endif

		if (x2 .lt. x) then 
			error =1
		endif

		if (x2 .le. x1) then 
			error =-1
			lin_interpolate = nan_f32
			return
		endif

		wt1 = (x2 - x)/(x2 - x1)
		wt2 = 1.0 - wt1

		lin_interpolate = wt1*y1 + wt2*y2

		return

	end function lin_interpolate
	
	real(4) function lin_interpolate_no_bounds_check(x,x1,x2,y1,y2,error)

		real(4),intent(in)			:: x
		real(4),intent(in)			:: x1
		real(4),intent(in)			:: x2
		real(4),intent(in)			:: y1
		real(4),intent(in)			:: y2
		integer(4),intent(inout)	:: error

		real(4)						:: wt1,wt2

		error = 0
		if (.not. ieee_is_finite(x))  error = -1
		if (.not. ieee_is_finite(x1)) error = -1
		if (.not. ieee_is_finite(x2)) error = -1
		if (.not. ieee_is_finite(y1)) error = -1
		if (.not. ieee_is_finite(y2)) error = -1
		lin_interpolate_no_bounds_check = nan_f32
		if (error .ne. 0)  then
			lin_interpolate_no_bounds_check = nan_f32
			RETURN
		endif


		
		! do some bounds checking

		wt1 = (x2 - x)/(x2 - x1)
		wt2 = 1.0 - wt1

		lin_interpolate_no_bounds_check = wt1*y1 + wt2*y2

		return

	end function lin_interpolate_no_bounds_check

	function lin_interpolate_directions(x,x1,x2,dir1,dir2,error)

		real(4)			:: x
		real(4)			:: x1
		real(4)			:: x2
		real(4)			:: dir1
		real(4)			:: dir2
		integer(4)		:: error

		real(4)			:: temp_dir1
		real(4)			:: temp_dir2
		real(4)			:: lin_interpolate_directions

		error = 0
		if (.not. ieee_is_finite(x)) error = -1
		if (.not. ieee_is_finite(x1)) error = -1
		if (.not. ieee_is_finite(x2)) error = -1
		if (.not. ieee_is_finite(dir1)) error = -1
		if (.not. ieee_is_finite(dir2)) error = -1
		if (error .ne. 0) return


		temp_dir1 = modulo(dir1,360.0)
		temp_dir2 = modulo(dir2,360.0)

		if (abs(temp_dir1 - temp_dir2) .gt. 180.0) then 
			if ((temp_dir1 .lt. 180) .and. (temp_dir2 .gt. 180.0)) temp_dir2 = temp_dir2 - 360.0
			if ((temp_dir2 .ge. 180) .and. (temp_dir2 .lt. 180.0)) temp_dir2 = temp_dir2 + 360.0
		endif

		lin_interpolate_directions = modulo(lin_interpolate(x,x1,x2,temp_dir1,temp_dir2,error),360.0)

		return
	end function lin_interpolate_directions


	real(4) function exp_interpolate(x,x1,x2,y1,y2,error)

		! interpolates assuming the two endpoints are connected by an exponential
		! function i.e. y = y1*exp(alpha*(x-x1))
		!
		! Error is set to 1 if extrapolation occurs
		! Error is set to -1 if x's out of order, or either y is negative

		real(4),intent(in)			:: x
		real(4),intent(in)			:: x1
		real(4),intent(in)			:: x2
		real(4),intent(in)			:: y1
		real(4),intent(in)			:: y2
		integer(4),intent(inout)	:: error

		real(4)						:: alpha
		
		! do some bounds checking
		
		if (x .lt. x1) then 
			error =1
		endif

		if (x2 .lt. x) then 
			error =1
		endif

		if ((x2 .le. x1) .or. (y1 .le. 0) .or. (y2 .le. 0)) then 
			error =-1
			exp_interpolate = nan_f32
			return
		endif
		alpha = log(y2/y1)/(x2 - x1) 
		exp_interpolate = y1*exp(alpha*(x-x1))

		return

	end function exp_interpolate

	function wt_mean_angle(deg, wts)

		use trig_degrees, only: cosd, sind, atan2d
		real(4) :: wt_mean_angle
		real(4), intent(in) :: deg(:)
		real(4), intent(in) :: wts(:)

		real(4) :: tot_wt_cos,tot_wt_sin,tot_wt

		tot_wt = 0.0
		tot_wt_cos = 0.0
		tot_wt_sin = 0.0

		do i = 1, size(deg)
			tot_wt_cos = tot_wt_cos + wts(i)*cosd(deg(i))
			tot_wt_sin = tot_wt_sin + wts(i)*sind(deg(i))
			tot_wt = tot_wt + wts(i)
		enddo

		wt_mean_angle = atan2d(tot_wt_sin,tot_wt_cos)
		wt_mean_angle = modulo(wt_mean_angle+1080.0,360.0)
	end function wt_mean_angle

end module interpolate









