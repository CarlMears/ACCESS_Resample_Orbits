$freeform
	module NAN_support
        USE, INTRINSIC :: IEEE_ARITHMETIC
        
		implicit none
		
		   
		real(4),parameter  :: NAN = O'017760000000'
        real(8),parameter  :: NAN8 = O'0777610000000000000000'
        
		interface finite
			module procedure finite4
			module procedure finite8
		end interface

	contains

		logical function finite4(x)

		! This function determines if x is a finite number
		! Returns 1 if X is finite, and 0 if x is infinite, or a non-zero denormal

			real(4),intent(in)			::	x
		    finite4 = IEEE_IS_FINITE(x)
			

			return

		end function finite4



		logical function finite8(x)

		! This function determines if x is a finite number
		! Returns 1 if X is finite, and 0 if x is infinite, or a non-zero denormal

			real(8),intent(in)			::	x

			finite8 = IEEE_IS_FINITE(x)
			return

		end function finite8
		
		subroutine init_nan
		
		end subroutine init_nan
	end module NAN_support 


