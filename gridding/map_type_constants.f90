module map_type_constants

    use, intrinsic :: iso_fortran_env, only: int8, int32, real32, real64

    integer(int32),parameter                :: CELL_CENTERED_2D_MAP = 1
    integer(int32),parameter                :: EDGE_CENTERED_2D_MAP = 2
    real(real32),parameter                  :: SMALL_NUMBER  = 1.0e-3

end module map_type_constants