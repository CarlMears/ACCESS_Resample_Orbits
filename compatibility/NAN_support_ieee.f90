
module NAN_support
    use, intrinsic :: iso_fortran_env, only: real32, real64, int64, int32
    
    implicit none
    
    private

    real(real64), parameter :: nan_f64 =  transfer(-2251799813685248_int64, 1._real64)
    real(real32), parameter :: nan_f32 =  transfer(-4194304_int32, 1._real32)

    public :: nan_f64
    public :: nan_f32

end module NAN_support

	
	